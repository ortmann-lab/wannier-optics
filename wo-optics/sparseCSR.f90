!**********************************************************************************!!
!! This module provides routines for setting up and operating on a sparse 
!! Hamiltonian distributed across MPI processes. The Hamiltonian is stored 
!! in Compressed Sparse Row (CSR) format, enabling efficient parallel M*V products.
!!
!! Main author: Frank Ortmann
!! Notes: Requires MPI.
!**********************************************************************************
    module SparseCSR
    use MPI 
    implicit none

    PRIVATE

    integer, save :: p_id, comm_world, istat(MPI_STATUS_SIZE)  ! Local variables related to MPI
    
    integer, save :: nsites, nbar                         ! Total number of sites and number of neighbor couplings for each site
    integer, dimension(:,:), pointer, save :: neiList    
    double complex, dimension(:), pointer, save :: enList 
    double complex, dimension(:,:), pointer, save :: traList  
    double complex, dimension(:), pointer, save :: PsiIn_p, PsiIn_t_p, PsiIn_0_p  ! Arrays for states 
    double complex, dimension(:), pointer, save :: Psi_p, Psi_t_p, Psi_0_p        ! wave function

    integer, dimension(:), pointer, save :: sites_p              ! Site indices handled by process
    integer, save :: sites_p_len, int_sites_len, bound_sites_len ! Numbers of all sites handled, internal sites, and boundary sites belonging to each process
    integer, dimension(:,:), allocatable, save :: worI           ! working array for neiList with corresponding processes
    double complex, dimension(:,:), allocatable, save :: worC    ! working array for traList, enList, PsiIn_p, and PsiIn_t_p etc.
    integer, dimension(:), allocatable, save :: sites_perm_p, old_order ! 
    integer, dimension(:), allocatable, save :: new_order        ! 
    
    double complex, dimension(:), allocatable, save :: HH_int    ! actual Hamiltonian in fully compressed form (internal)
    integer, dimension(:), allocatable, save :: Hj_int           ! CSR index array indicating column
    integer, dimension(:), allocatable, save :: Hi_int           ! CSR index array indicating start position of new row in HH
    integer, save :: HH_int_len
    double complex, dimension(:), allocatable, save :: HH_ext    !actual Hamiltonian in fully compressed form (external matrix elements)
    integer, dimension(:), allocatable, save :: Hj_ext           ! CSR index array indicating column
    integer, dimension(:), allocatable, save :: Hi_ext           ! CSR index array indicating start position of new row in HH
    integer, save :: HH_ext_len
    integer, save :: nnp                                         ! Number of neighboring processes for communication
    
    integer, dimension(:), allocatable, save :: send_pid         ! process id of communication partner
    integer, dimension(:), allocatable, save :: send_ind         ! index array for packed send_sites_p and send_loc_p arrays.
    integer, dimension(:), allocatable, save :: recv_sites_p, recv_i_p  ! site indices received and process id
    
    double complex, dimension(:), allocatable, save :: extVec    !for communication
    double complex, dimension(:), allocatable, save :: extVec1
    double complex, dimension(:), allocatable, save :: extVec2
    integer, dimension(:), allocatable, save :: send_loc_p        
    double complex, dimension(:), allocatable, save :: send_val_p 
    double complex, dimension(:), allocatable, save :: send_val1_p
    double complex, dimension(:), allocatable, save :: send_val2_p
    integer, dimension(:), allocatable, save :: send_rqst, recv_rqst, send_rqst2, recv_rqst2
    double complex, dimension(:), allocatable, save :: extRes   

    double complex, dimension(:), pointer, save :: locRes, psiN, psiNm1 ! Working arrays for time development
    integer, save :: ierr  ! Error variable: if ierr .ne. 0, an error occurred
    
    ! Routines and functions for use from outside this module
    PUBLIC setHamiltPsi
    PUBLIC setupSparseCSRModule, ImplemCSR_IB, allocateArrays
    PUBLIC multHpsi, time_evolution_Psi
    PUBLIC getSites_p, getPsi_p, getPsi_t_p
    PUBLIC deallocateArrays, encapsHamiltonian
    PUBLIC getPsi_0_p

    CONTAINS

    subroutine setupSparseCSRModule(sites_p_a,p_id_a, comm_world_a, nsites_a, nbar_a)
        !***************************************************************************
        !> Initializes module by setting up some module-scope variables.
        !!
        !! Called by all processes
        !! Module-scope variables: sites_p, p_id, comm_world, nsites, nbar
        !***************************************************************************
        implicit none
        integer, intent(IN) :: p_id_a, comm_world_a, nsites_a, nbar_a
        integer, target, intent(IN) :: sites_p_a(sites_p_len)

        p_id = p_id_a; comm_world = comm_world_a
        nsites = nsites_a; nbar = nbar_a
        sites_p => sites_p_a 
    end subroutine setupSparseCSRModule

    function getSites_p()
        !***************************************************************************
        !> Returns number of sites handled by current MPI-process
        !!
        !! Called by all processes
        !! Module-scope variable: sites_p_len
        !****************************************************************************
        implicit none
        integer :: getSites_p

        getSites_p = sites_p_len
    end function getSites_p

    subroutine getPsi_p(Psi_p_a, psiLSize)
        !***************************************************************************
        !> Returns Psi_p and size of pointed vector
        !!
        !! Called by all processes
        !! Module-scope variables: Psi_p, sites_p_len
        !****************************************************************************
        implicit none
        double complex, dimension(:), pointer :: Psi_p_a
	    integer, intent(OUT) :: psiLSize

        Psi_p_a => Psi_p
	    psiLSize = sites_p_len
    end subroutine getPsi_p

    subroutine getPsi_t_p(Psi_t_p_a, psitLSize)
        !***************************************************************************
        !> Returns Psi_t_p and size of pointed vector
        !!
        !! Called by all processes
        !! Module-scope variables: Psi_t_p, sites_p_len
        !****************************************************************************
        implicit none
        double complex, dimension(:), pointer :: Psi_t_p_a
        integer, intent(OUT) :: psitLSize

        Psi_t_p_a => Psi_t_p
        psitLSize = sites_p_len
    end subroutine getPsi_t_p
    
    subroutine getPsi_0_p(Psi_0_p_a, psi0LSize)
        !***************************************************************************
        !> Returns Psi_0_p and size of pointed vector
        !!
        !! Called by all processes
        !! Module-scope variables: Psi_0_p, sites_p_len
        !****************************************************************************
        implicit none
        double complex, dimension(:), pointer :: Psi_0_p_a
        integer, intent(OUT) :: psi0LSize

        Psi_0_p_a => Psi_0_p
        psi0LSize = sites_p_len
    end subroutine getPsi_0_p

    subroutine setHamiltPsi(nsites_a, nbar_a, sites_p_len_a,neiList_a, traList_a, enList_a, PsiIn_p_a, PsiIn_t_p_a,&
                                & PsiIn_0_p_a)
        !***************************************************************************
        !> Sets Hamiltonian data and states in the module
        !! The Hamiltonian is created in a CSR format inside the hamiltonian module.
        !! This subroutine moves it into the module such that we can perform matrix-vector operations
        !!
        !! Called by all processes
        !! Module-scope pointers: neiList, traList, enList
        !! Module-scope arrays and variables: PsiIn_p, PsiIn_t_p, PsiIn_0_p, sites_p_len
        !****************************************************************************
        implicit none
        integer, intent(IN) :: nsites_a, nbar_a, sites_p_len_a
        integer, target, intent(IN) :: neiList_a(sites_p_len_a,nbar_a)
        double complex, target, intent(IN) :: traList_a(sites_p_len_a,nbar_a), enList_a(sites_p_len_a)
        double complex, target, intent(IN) ::  PsiIn_p_a(sites_p_len_a), PsiIn_t_p_a(sites_p_len_a)
        double complex, target, intent(IN) :: PsiIn_0_p_a(sites_p_len_a)

        neiList => neiList_a; traList => traList_a; enList => enList_a 
        PsiIn_p => PsiIn_p_a; PsiIn_t_p => PsiIn_t_p_a; PsiIn_0_p => PsiIn_0_p_a 
        sites_p_len = sites_p_len_a  
    end subroutine setHamiltPsi

    subroutine ImplemCSR_IB(ierr_a)
        !***************************************************************************
        !> Sorts CSR data in internal and boundary sites, including:
        !! creating appropriate data structures at each process
        !! sending/receiving boundary parts of Hamiltonian between processes
        !! 
        !! Called by all processes
        !! Module-scope variables: comm_world, ierr
        !****************************************************************************
        implicit none
        integer, intent(OUT) :: ierr_a

        call MPI_BARRIER(comm_world, ierr_a)
        if (p_id .eq. 0) write(*,*) '  Step 1 of 5 in CSR implementation (reorder local sites)'
        call classifySites()
        call reorderArrays()
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR in classifySites or reorderArrays'            
            ierr_a = ierr
            return
        end if

        call MPI_BARRIER(comm_world, ierr)
        if (p_id .eq. 0) write(*,*) '  Step 2 of 5 CSR implementation (create local CSR part)'
        call CSR_local()
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR in CSR_local'            
            ierr_a = ierr
            return
        end if

        call MPI_BARRIER(comm_world, ierr)
        if (p_id .eq. 0) write(*,*) '  Step 3 of 5 CSR implementation (exchange boundaries)'
        call exchangeBoundaries()
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR in exchangeBoundaries'    
            ierr_a = ierr
            return
        end if

        call MPI_BARRIER(comm_world, ierr)
        if (p_id .eq. 0) write(*,*) '  Step 4 of 5 CSR implementation  (create external CSR part)'
        call CSR_external()
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR in CSR_external'            
            ierr_a = ierr
            return
        end if

        call MPI_BARRIER(comm_world, ierr)
        if (p_id .eq. 0) write(*,*) '  Step 5 of 5 CSR (allocate MatVecCommArrays)'
        call allocateVecCommArrays()
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR in allocateVecCommArrays'
            ierr_a = ierr
            return
        end if

        call MPI_BARRIER(comm_world, ierr)
        if (ierr .ne. 0) STOP "Error in MPI_BARRIER"
    end subroutine ImplemCSR_IB


    subroutine classifySites()
        !***************************************************************************
        !> Identifies internal and boundary sites and length of internal part of matrix
        !! 
        !! Called by all processes
        !! Module-scope variables: worI,Hsites_p_len,HH_int_len,int_sites_len,bound_sites_len
        !****************************************************************************
        implicit none
        integer :: il, jl, rr, rrs, rre
        logical :: flagI

        allocate(sites_perm_p(sites_p_len), STAT=ierr)
        allocate(new_order(sites_p_len), STAT=ierr)
        allocate(old_order(sites_p_len), STAT=ierr)
        if (ierr .ne. 0) STOP "Error in allocation in classifySites"

        ! Initialize
        rrs = 1
        rre = sites_p_len
        int_sites_len = 0
        bound_sites_len = 0
        HH_int_len = sites_p_len

        do il = 1, sites_p_len
            flagI= .TRUE.
            do jl = 1, nbar
                if (worI((nbar+jl), il) .ne. p_id+1) then
                    worI(jl, il) = -1*worI(jl, il)
                    flagI= .FALSE.
                else
                    HH_int_len = HH_int_len + 1
                end if
            end do
            if (flagI) then
                int_sites_len = int_sites_len + 1
                rr = rrs
                rrs = rrs + 1
            else
                bound_sites_len = bound_sites_len + 1
                rr = rre
                rre = rre - 1
            end if
            sites_perm_p(rr) = sites_p(il)
            new_order(il) = rr
            old_order(rr) = il
        end do
    end subroutine classifySites
    
    subroutine reorderArrays()
        !***************************************************************************
        !> brings arrays in  order (new_order)
        !! 
        !! Called by all processes
        !****************************************************************************
        implicit none
        integer :: il

        ! Allocate temporary arrays
        allocate(Psi_p(sites_p_len), STAT=ierr)
        allocate(Psi_t_p(sites_p_len), STAT=ierr)
        allocate(Psi_0_p(sites_p_len), STAT=ierr)
        if (ierr .ne. 0) STOP "Error in allocation in reorderArrays"

        ! Perform permutation
        do il = 1, sites_p_len
            Psi_p(new_order(il)) = worC((nbar+2), il)
            Psi_t_p(new_order(il)) = worC((nbar+3), il)
            Psi_0_p(new_order(il)) = worC((nbar+4), il)
        end do
    end subroutine reorderArrays
        
    subroutine CSR_local()
        USE Geometry
        !***************************************************************************
        !> Creates CSR data structure for internal part of the Hamiltonian HH
        !! 
        !! Called by all processes
        !! Module-scope variables: sites_p_len, nbar, rec_s,  worI, old_order, sites_perm_p, sites_p_len
        !!                         HH_int,  Hj_int, Hi_int
        !! Initializes: HH_int, Hj_int, Hi_int
        !****************************************************************************
        implicit none

        integer :: globCol, locCol
        double complex :: tmp
        integer :: i, j, ii, jj, rr, rs
        logical :: flag
        INTEGER :: locpos
        
        allocate(HH_int(HH_int_len), STAT=ierr)
        allocate(Hj_int(HH_int_len), STAT=ierr)
        allocate(Hi_int(sites_p_len+1), STAT=ierr)
        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR: cannot allocate CSR_loc arrays' 
            return
        end if
        rr = 0   !running index
        HH_int = -1
        Hj_int = sites_p_len + 1
        Hi_int = -1
        do i = 1, sites_p_len
            rr = rr + 1
            Hi_int(i) = rr
            rs = rr
            Hj_int(rr) = i
            HH_int(rr) = worC(nbar+1, old_order(i)) 
            do j = 1, nbar
                globCol = worI(j, old_order(i)) 
                if (globCol > 0) then ! only internals here
                    tmp = worC(j, old_order(i)) 
                    locpos = getLocalId(globCol)
                    locCol = new_order(locpos)
                    
                    if(locCol < 0) then
                        write(*,*)  "ERROR: cannot find local index of internal column" 
                        ierr = 1
                        return
                    end if
                    ii = rs
                    flag = .TRUE.
                    do while(flag)
                        if (Hj_int(ii) > locCol) then
                            if (ii < rr+1) then
                                do jj = rr, ii, -1
                                    Hj_int(jj+1) = Hj_int(jj)
                                    HH_int(jj+1) = HH_int(jj)
                                end do
                            end if
                            Hj_int(ii) = locCol
                            HH_int(ii) = tmp !enter value
                            rr = rr + 1
                            flag = .FALSE.
                        else
                            ii = ii + 1
                        end if
                    end do
                end if
            end do
        end do
        Hi_int(sites_p_len+1) = rr+1 !for completeness
        deallocate(new_order)
    end subroutine CSR_local

    subroutine CSR_external()
        !***************************************************************************
        !> Creates CSR data structure for external part
        !! Called by all processes
        !! Module-scope variables: recv_sites_p, int_sites_len, sites_p_len,nbar, worI,
        !!              recv_i_p, nnp, HH_ext, Hj_ext, Hi_ext, worC, old_order, p_id
        !!
        !! Allocates/Initializes: HH_ext, Hj_ext, Hi_ext
        !! Deallocates: recv_sites_p
        !****************************************************************************
        implicit none
        integer :: extGlobId, oldi
        integer :: i, j, ii, jj, rr
        logical :: flag1

        allocate(HH_ext(HH_ext_len))
        allocate(Hj_ext(HH_ext_len))
        allocate(Hi_ext(bound_sites_len+1))

        Hi_ext(1) = 1
        rr = 1
        jj = 1
        do i = int_sites_len+1, sites_p_len
            oldi=old_order(i)
            do j = 1, nbar
                extGlobId = worI(j, oldi)
                if (extGlobId < 0) then
                    flag1 = .TRUE.
                    ii = 1
                    extGlobId = -1 * extGlobId
                    do while (flag1)
                        if (recv_sites_p(ii) .eq. extGlobId) then
                            flag1 = .FALSE. 
                        else
                            ii = ii + 1
                            if(ii > recv_i_p(nnp+1)) then
                                write(*,*) p_id, "ERROR: cannot find local index (CSR_external)" 
                                return
                            end if
                        end if
                    end do
                    HH_ext(rr) = worC(j, oldi)
                    Hj_ext(rr) = ii
                    rr = rr + 1
                end if
            end do
            Hi_ext(jj+1) = rr
            jj = jj + 1
        end do
        deallocate(recv_sites_p)
    end subroutine CSR_external
    
    subroutine exchangeBoundaries()
        !***************************************************************************
        !> Prepares and exchanges boundary site indices for parallel M*V operations.
        !! Called by all MPI processes.
        !!
        !! Uses module-scope variables:
        !!   sites_p_len, int_sites_len, bound_sites_len, nbar, sites_p, worI,
        !!   old_order, comm_world, p_id, ierr,
        !!   HH_ext_len, send_pid, send_ind, send_loc_p, nnp,
        !!   recv_sites_p, recv_i_p
        !!
        !! Initializes/allocates:
        !!   HH_ext_len, send_pid, send_ind, send_loc_p, nnp,
        !!   recv_sites_p, recv_i_p
        !***************************************************************************

        implicit none

        integer :: oldi, globId, locId, mpipsize
        integer :: i, j, ii, jj
        integer :: rr, rrs, rrlen !set of running indices
        integer :: msg, mm
        logical :: flag0, rowToadd

        integer, dimension(:), allocatable :: send_sites_p, send_len, chckTab                                                             
        integer, dimension(:), allocatable :: recv_len 
        integer, dimension(:), allocatable :: priv_s_rqst, priv_r_rqst

        call MPI_Comm_size(comm_world, mpipsize, ierr)
        if(ierr .ne. 0) stop 'ERROR: cannot get number of processes from MPI'

        msg = 123

        allocate(send_pid(mpipsize), STAT=ierr)
        allocate(send_ind(mpipsize+1), STAT=ierr)
        allocate(chckTab(nbar), STAT=ierr)
        
        ! Allocate maximum possible size: each boundary site might need to be sent to every neighbor,
        ! as determined by worI; actual usage may be less
        allocate(send_sites_p(nbar*bound_sites_len), STAT=ierr)
        allocate(send_loc_p(nbar*bound_sites_len), STAT=ierr) 

        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR: cannot allocate send arraysB1'
            return
        end if
        send_pid = -1
        send_ind = 0    !index array

        nnp = 0
        rr = 1
        HH_ext_len = 0

        do i = int_sites_len+1, sites_p_len !loop over boundary sites
            oldi=old_order(i)
            globId = sites_p(oldi) !global index of boundary site
            locId = i !local index of boundary site
            do mm=1,nbar
               chckTab(mm) = -1
            enddo
            do j = 1, nbar
                if(worI(j, oldi) < 0) then !external site
                    HH_ext_len = HH_ext_len + 1
                    flag0 = .TRUE.
                    ii = 1
                    do while (flag0)
                        if (send_pid(ii) .eq. -1) then !Add the process in slot ii as a new neighbor process to send to
                            flag0 = .FALSE.
                            nnp = nnp + 1
                            send_pid(ii) = worI(nbar+j, oldi)
                            chckTab(j) = worI(nbar+j, oldi)
                            send_ind(nnp) = rr
                            send_ind(nnp+1) = rr + 1
                            send_sites_p(rr) = globId
                            send_loc_p(rr) = locId
                            rr = rr + 1
                        else if (send_pid(ii) .eq. worI(nbar+j, oldi)) then !This site belongs to a neighbor process that is already being scheduled as target
                            flag0 = .FALSE.
                            rowToadd = .TRUE.
                            jj = 1
                            do while (rowToadd .AND. jj < j)
                                if(send_pid(ii) .eq. chckTab(jj)) then
                                    rowToadd = .FALSE.
                                else
                                    jj = jj + 1
                                end if
                            end do
                            if(rowToadd) then
                                chckTab(j) = worI(j+nbar, oldi)
                                if (send_ind(ii+1) .eq. rr) then
                                    send_sites_p(rr) = globId
                                    send_loc_p(rr) = locId
                                    send_ind(ii+1) = rr + 1
                                    rr = rr + 1
                                else
                                    do jj = rr, send_ind(ii+1)+1, -1
                                        send_sites_p(jj) = send_sites_p(jj-1)
                                        send_loc_p(jj) = send_loc_p(jj-1)
                                    end do
                                    send_sites_p(send_ind(ii+1)) = globId
                                    send_loc_p(send_ind(ii+1)) = locId
                                    rr = rr + 1
                                    do jj = ii+1, nnp+1
                                        send_ind(jj) = send_ind(jj) + 1
                                    end do
                                end if
                            end if
                        else
                            ii = ii + 1
                        end if
                    end do
                end if
            end do
        end do

        !exchange with neighbouring processes first the lenght of data array
        !send and receive messages in a deadlock-free mode
        allocate(priv_s_rqst(nnp), STAT=ierr)
        allocate(priv_r_rqst(nnp), STAT=ierr)

        allocate(send_len(nnp), STAT=ierr)
        allocate(recv_len(nnp), STAT=ierr)

        if(ierr .ne. 0) then
            write(*,*) p_id, 'ERROR: cannot allocate send arraysB2'
            return
        end if
        rrlen = 0

        do i = 1, nnp
            send_len(i) = send_ind(i+1) - send_ind(i)

            call MPI_ISEND(send_len(i), 1, MPI_INTEGER, send_pid(i)-1, &
            &  msg,MPI_COMM_WORLD, priv_s_rqst(i), ierr)

            call MPI_IRECV(recv_len(i), 1, MPI_INTEGER, send_pid(i)-1, &
            &  msg,MPI_COMM_WORLD, priv_r_rqst(i), ierr)

        end do
        do i = 1, nnp
            call MPI_WAIT(priv_s_rqst(i), istat, ierr)
            call MPI_WAIT(priv_r_rqst(i), istat, ierr)
            rrlen = rrlen + recv_len(i)   !size of receiving data
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0) STOP "Error in MPI_BARRIER CSR exchangeB3"

        !now exchange neighbouring indices of my boundary sites
        allocate(recv_sites_p(rrlen))
        allocate(recv_i_p(nnp+1))
        rr = 1
        recv_i_p(1) = 1

        do i = 1, nnp   !packed array of sites
            call MPI_ISEND(send_sites_p(send_ind(i) : send_ind(i+1)-1),&
                & send_len(i), MPI_INTEGER, send_pid(i)-1, msg, &
                & MPI_COMM_WORLD, priv_s_rqst(i), ierr)
            if(ierr .ne. 0) then
                write(*,*) p_id, 'ERROR: cannot send site indices in CSR exchangeB4'
                call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
            endif
            rrs = rr
            rr = rr+recv_len(i)

            call MPI_IRECV(recv_sites_p(rrs), recv_len(i), &     ! dont use recv_sites_p(rrs:rr)
                & MPI_INTEGER, send_pid(i)-1, msg, &
                & MPI_COMM_WORLD, priv_r_rqst(i), ierr)
            if(ierr .ne. 0) then
                write(*,*) p_id, 'ERROR: cannot recv site indices in CSR exchangeB5'
                call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
            endif
            recv_i_p(i+1) = rr
        end do
        do i = 1, nnp
            call MPI_WAIT(priv_s_rqst(i), istat, ierr)
            call MPI_WAIT(priv_r_rqst(i), istat, ierr)
        end do

        deallocate(send_len)
        deallocate(recv_len)
        deallocate(send_sites_p)
        deallocate(priv_s_rqst)
        deallocate(priv_r_rqst)
    end subroutine exchangeBoundaries

    subroutine encapsHamiltonian()
        !***************************************************************************
        !> encapsulate Hamiltonian on the process
        !!
        !! Called by all processes
        !! Module-scope variables: nbar, neiList, traList, enList, PsiIn_p, PsiIn_t_p, 
        !!              nP, nP_i, nP_sites, sites_p, sites_p_len, worI, worC
        !! Initializes: worI, worC
        !****************************************************************************
        implicit none
        integer :: ii, jj

        allocate(worI((nbar*2), sites_p_len), STAT=ierr)
        if(ierr .ne. 0) then
            write(*,*) p_id, ' ERROR: cannot allocate MPI receive arrays'                               
            return
        endif      
 
        !This should write to worI contiguously in memory 
        do jj = 1, sites_p_len
            do ii = 1, nbar
               worI(ii, jj) = neiList(jj, ii)
               worI((nbar+ii), jj) = neiList(jj, ii + nbar)                  
            end do
        end do

        allocate(worC((nbar+4), sites_p_len), STAT=ierr)
        if(ierr .ne. 0) then
            write(*,*) p_id, ' ERROR: cannot allocate MPI receive arrays'                               
            return
        endif

        do jj = 1, sites_p_len
            do ii = 1, nbar
                worC(ii, jj) = DCMPLX(traList(jj, ii))
            end do
        end do

        do jj = 1, sites_p_len
            worC((nbar+1), jj) = enList(jj)
            worC((nbar+2), jj) = PsiIn_p(jj)
            worC((nbar+3), jj) = PsiIn_t_p(jj)
            worC((nbar+4), jj) = PsiIn_0_p(jj)
        end do

    end subroutine encapsHamiltonian

    subroutine VecComm(vector)
        !***************************************************************************
        !> Communicates the boundary elements of a local vector to neighboring
        !! processes, used for each M*V operation
        !! Called by all processes
        !! Module-scope variables: nnp, send_ind, send_loc_p, send_val_p, send_pid,
        !!              send_rqst, recv_rqst, recv_i_p, extVec, ierr, MPI_COMM_WORLD 
        !****************************************************************************
        implicit none
        double complex, dimension(:), pointer :: vector
        integer :: sendL, msg
        integer :: ii, rrs, rre

        msg = 111
        !create buffer of data to be sent
        do ii = 1, send_ind(nnp+1)-1
            send_val_p(ii) = vector(send_loc_p(ii))
        end do

        do ii = 1, nnp
            sendL = send_ind(ii+1) - send_ind(ii)
            call MPI_ISEND(send_val_p(send_ind(ii)), sendL, MPI_DOUBLE_COMPLEX, &
                           send_pid(ii)-1, msg, MPI_COMM_WORLD, send_rqst(ii), ierr)

            rrs = recv_i_p(ii)
            rre = recv_i_p(ii+1) - 1
            call MPI_IRECV(extVec(rrs), rre-rrs+1, MPI_DOUBLE_COMPLEX, &
                           send_pid(ii)-1, msg, MPI_COMM_WORLD, recv_rqst(ii), ierr)
        end do
    end subroutine VecComm

    subroutine VecComm_wait()
        !***************************************************************************
        !> Ensures that communication started at VecComm has finished
        !! Called by all processes
        !! Uses module-scope: nnp, send_rqst, recv_rqst, istat, ierr.
        !****************************************************************************
        implicit none
        integer :: ii

        do ii = 1, nnp
            call MPI_WAIT(send_rqst(ii), istat, ierr)
            call MPI_WAIT(recv_rqst(ii), istat, ierr)
        end do
    end subroutine VecComm_wait

    subroutine multHpsi(psiN_p, hPsiN_p, DotSum)
        !***************************************************************************
        !> M*V product: applies Hamiltonian to vector psiN_p, resulting in hPsiN_p. 
        !! Computes additionally the local dot product DotSum.
        !! Module-scope variables: sites_p_len, int_sites_len, bound_sites_len,
        !!                 HH_int, Hj_int, Hi_int, HH_ext, Hj_ext, Hi_ext, extVec, extRes.
        !! Called by all processes
        !****************************************************************************
        implicit none      
        double complex, intent(OUT) :: DotSum    
        double complex, dimension(:), pointer :: psiN_p, hPsiN_p    

        !locals
        double complex :: tmp
        integer :: i, j
        integer :: rr  !running index

        call VecComm(psiN_p) !Boundary parts are exchanged
        !For efficiency: local M*V multiplication performed during communication

        DotSum = (0.d0, 0.d0)
        do i = 1, sites_p_len
            tmp = (0.d0, 0.d0)
            do j = Hi_int(i), Hi_int(i+1)-1
                tmp = tmp + HH_int(j) * psiN_p(Hj_int(j))
            end do
            hPsiN_p(i) = tmp
            if(i .le. int_sites_len) then
            !at this point only M*V prod of internal sites is done
            !and partial dot product can be computed
                DotSum = DotSum + CONJG(psiN_p(i)) * hPsiN_p(i)
            end if
        end do

        call VecComm_wait() !Waits for completion of communication started with VecComm

        !calculate external M*V product
        do i = 1, bound_sites_len
            tmp = (0.d0, 0.d0)
            do j = Hi_ext(i), Hi_ext(i+1)-1
                tmp = tmp + HH_ext(j) * extVec(Hj_ext(j))
            end do
            extRes(i) = tmp
        end do

        !add external M*V product to the boundary sites
        rr = 1
        do i = int_sites_len+1, sites_p_len
            hPsiN_p(i) = hPsiN_p(i) + extRes(rr)
            rr = rr + 1
            !adding partial dot product for boundary sites
            DotSum = DotSum + CONJG(psiN_p(i)) * hPsiN_p(i)
        end do
    end subroutine multHpsi

    subroutine time_evolution_Psi(Npol, ac, bc, cCoeff)
        !***************************************************************************
        !> Time evolution step: Chebyshev propagation using the Hamiltonian and
        !! current vector Psi_t_p.
        !! Called by all processes
        !!
        !! Uses module-scope variables:
        !!   Psi_t_p, psiN, psiNm1, locRes, extRes, extVec,
        !!   sites_p_len, int_sites_len, bound_sites_len,
        !!   HH_int, Hj_int, Hi_int, HH_ext, Hj_ext, Hi_ext,
        !!   recv_i_p, send_pid, send_loc_p, send_reqs, recv_reqs, nnp
        !***************************************************************************

        implicit none
        integer, intent(IN) :: Npol
        double precision, intent(IN) :: ac, bc
        double complex:: cCoeff(Npol)
        integer :: i, j, ii, rr
        double complex :: tmp
        double complex, dimension(:), pointer :: swapp
        
        !separating first and subsequent steps ii=1
        do i = 1, sites_p_len
            psiNm1(i) = Psi_t_p(i)*SQRT(2.d0)
        end do
        call VecComm(psiNm1) !Boundary parts are exchanged
        !For efficiency: local M*V multiplication performed during communication
        do i = 1, sites_p_len
            tmp = (0.d0, 0.d0)
            do j = Hi_int(i), Hi_int(i+1)-1
                if(Hj_int(j) .eq. i) then !remember diagonal element
                    tmp = tmp + DREAL(HH_int(j)-ac) * psiNm1(Hj_int(j))
                else
                    tmp = tmp + HH_int(j) * psiNm1(Hj_int(j))
                end if
            end do
            locRes(i) = tmp
            
            if(i .le. int_sites_len) then
                !at this point only M*V prod of internal sites is done
                !and can be used
                psiN(i) = locRes(i) / (2.d0 * bc)
            end if
            Psi_t_p(i) = cCoeff(1) * Psi_t_p(i)
        end do

        call VecComm_wait() !Waiting for completion of communication started with VecComm

        !calculate external M*V product
        do i = 1, bound_sites_len
            tmp = (0.d0, 0.d0)
            do j = Hi_ext(i), Hi_ext(i+1)-1
                tmp = tmp + HH_ext(j) * extVec(Hj_ext(j))
            end do
            extRes(i) = tmp
        end do

        !add external M*V product to the boundary sites
        rr = 1
        do i = int_sites_len+1, sites_p_len
            locRes(i) = locRes(i) + extRes(rr)
            rr = rr + 1
            !now result for boundary sites can be used
            psiN(i) = locRes(i) / (2.d0 * bc)
        end do

        do ii = 2, Npol
            call VecComm(psiN) !Boundary parts are exchanged
            !For efficiency: local M*V multiplication performed during communication
            do i = 1, sites_p_len
                tmp = (0.d0, 0.d0)
                do j = Hi_int(i), Hi_int(i+1)-1
                    if(Hj_int(j) .eq. i) then !remember diagonal element
                        tmp = tmp + DREAL(HH_int(j)-ac) * psiN(Hj_int(j))
                    else
                        tmp = tmp + HH_int(j) * psiN(Hj_int(j))
                    end if
                end do
                locRes(i) = tmp
                if(i .le. int_sites_len) then
                !only M*V result for internal atoms is completed, so use it
                    psiNm1(i) = (locRes(i) - bc*psiNm1(i)) / bc
                end if
                Psi_t_p(i) = Psi_t_p(i) + cCoeff(ii) * psiN(i)
            end do
           
            call VecComm_wait() !Waiting for completion of communication started at VecComm

        !calculate external M*V product
            do i = 1, bound_sites_len
                tmp = (0.d0, 0.d0)
                do j = Hi_ext(i), Hi_ext(i+1)-1
                    tmp = tmp + HH_ext(j) * extVec(Hj_ext(j))
                end do
                extRes(i) = tmp
            end do

            !add external M*V product to the boundary sites
            rr = 1
            do i = int_sites_len+1, sites_p_len
                locRes(i) = locRes(i) + extRes(rr)
                rr = rr + 1
                !now result for boundary sites can be used
                psiNm1(i) = (locRes(i) - bc*psiNm1(i)) / bc
            end do
            swapp => psiNm1
            psiNm1 => psiN
            psiN   => swapp
        end do
    end subroutine time_evolution_Psi
    
    subroutine allocateArrays(ierr_a)
        !***************************************************************************
        !> Allocates arrays used in time_evolution_Psi
        !!
        !! Called by all processes
        !! Module-scope variables: sites_p_len, locRes, psiN, psiNm1
        !****************************************************************************
        implicit none
        integer, intent(OUT) :: ierr_a
	
        allocate(locRes(sites_p_len), STAT=ierr)
        allocate(psiN(sites_p_len), STAT=ierr)
        allocate(psiNm1(sites_p_len), STAT=ierr)
        if (ierr .ne. 0) then
            write(*,*) p_id, ' ERROR: cannot allocate arrays'
            ierr_a = ierr
            return
        end if
    end subroutine allocateArrays

    subroutine allocateVecCommArrays()
        !***************************************************************************
        !> Allocates arrays for storing data used in communication
        !! Called by all processes
        !****************************************************************************
        implicit none
	    integer :: err
        allocate(send_val_p(send_ind(nnp+1)-1), STAT=err)
        allocate(send_val1_p(send_ind(nnp+1)-1), STAT=err)
        allocate(send_val2_p(send_ind(nnp+1)-1), STAT=err)
        allocate(send_rqst(nnp), STAT=err)
        allocate(recv_rqst(nnp), STAT=err)
        allocate(send_rqst2(nnp), STAT=err)
        allocate(recv_rqst2(nnp), STAT=err)
        allocate(extVec(recv_i_p(nnp+1)-1), STAT=err)
        allocate(extVec1(recv_i_p(nnp+1)-1), STAT=err)
        allocate(extVec2(recv_i_p(nnp+1)-1), STAT=err)
        allocate(extRes(bound_sites_len), STAT=err)
        if(err .ne. 0) then
            write(*,*) p_id, ' ERROR: cannot allocate vec comm arrays' 
            return
        end if
    end subroutine allocateVecCommArrays

	subroutine deallocateArrays()
        !****************************************************************************
        !> Deallocates all arrays
        !!
        !! Called by all processes
        !****************************************************************************
	    implicit none
	    integer :: err
        deallocate(HH_int,STAT=err)
	    deallocate(Hj_int,STAT=err)
	    deallocate(Hi_int,STAT=err)
        deallocate(HH_ext,STAT=err)
	    deallocate(Hj_ext,STAT=err)
	    deallocate(Hi_ext,STAT=err)

	    deallocate(send_ind,STAT=err)
	    deallocate(send_pid,STAT=err)
	    deallocate(send_loc_p,STAT=err)
        deallocate(send_val_p,STAT=err)
        deallocate(send_val1_p,STAT=err) 
        deallocate(send_val2_p,STAT=err)
        deallocate(send_rqst,STAT=err) 
        deallocate(send_rqst2,STAT=err)
	    deallocate(recv_i_p,STAT=err)
	    deallocate(recv_rqst2,STAT=err)
        deallocate(recv_rqst,STAT=err) 

	    deallocate(worI,STAT=err)
	    deallocate(worC,STAT=err)

        deallocate(extRes,STAT=err)
        deallocate(extVec,STAT=err) 
        deallocate(extVec1,STAT=err)
        deallocate(extVec2,STAT=err) 
		if (err.ne.0) then
			write(*,*) p_id,'ERROR: cannot deallocate arrays'
			STOP
		endif

	end subroutine deallocateArrays

end module SparseCSR