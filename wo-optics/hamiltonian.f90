
! This module creates the exciton Hamiltonian in a sparse matrix format
! for each local domain.
!
! Author: Konrad Merkel


module Hamiltonian
    USE MPI
    USE PARFILE_mod
    USE Geometry

    implicit none
    PRIVATE
        
    PUBLIC create_exciton_hamiltonian, estimate_num_neighbors
    contains

    SUBROUTINE estimate_num_neighbors(p_id, NNEIGH)
        !
        ! Estimates the number of neighboring sites (number of off-diagonal matrix entries in the CSR matrix) of the Hamiltonian.
        ! This number has to be known in advance to be able to allocate optimal amount of memory.
        ! The number is estimated by going through the construction of a minimal test Hamiltonian (without materializing it) and
        ! counting the number of sites that would be needed.
        !

        IMPLICIT NONE
        INTEGER, intent(in) ::  p_id
        INTEGER, intent(out) :: NNEIGH

        INTEGER ia, i, numTI
        character(len=1024) :: filename
        integer max_shift_i(3), max_shift_o(3), test_supercell(3)
        integer maxS(3)

        ! Variables for the calculation of neighbors
        INTEGER ALLOC_ERR, elementId, ierr
        integer NumExcitTI, NumExcitTIused
        DOUBLE PRECISION gamma
        DOUBLE COMPLEX, allocatable  :: TIsingle(:,:,:,:,:,:,:), TIexciton(:)
        INTEGER, allocatable :: NNMAP(:,:), NNEIGH_per_site(:,:,:,:,:)
        INTEGER, allocatable :: NNMAPexciton(:,:)

        INTEGER :: c1,c2,v1,v2, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2
        logical :: useTIexciton
        logical, allocatable ::  isUsedTIexciton(:)

        !------------------------ End variable declaration --------------------------

        NumExcitTIused=0


        ! get max shift and allocate TI arrays
        if(p_id.eq.0) THEN
            write(*,*) ""
            write(*,*) "Open files to estimate the number of neighbors..."
        endif
        filename="TINFILE_v"
        call get_max_shift_from_tinfile(p_id, filename, max_shift_o)
        filename="TINFILE_c"
        call get_max_shift_from_tinfile(p_id, filename, max_shift_i)

        ! get maximal shift for all single particle TIs
        do i=1,3
            max_shift_i(i) = max(max_shift_i(i), max_shift_o(i))
        enddo

        allocate(TIsingle(get_N_con(), get_N_val(), get_N_con(), get_N_val(), -max_shift_i(1):max_shift_i(1), &
            & -max_shift_i(2):max_shift_i(2),-max_shift_i(3):max_shift_i(3)))

        ! read all single particle TIs and the nearest neighbor map
        call read_tinfile(p_id, numTI, NNMAP, TIsingle, gamma)

        ! read (non-diagonal) Coulomb and local field effects matrix elements similar as TI
        call get_max_shift_interaction_elements(p_id, maxS, NumExcitTI, useTIexciton)


        ! allocate data structure for excitonic transfer integrals
        if (useTIexciton) then
            allocate(TIexciton(NumExcitTI))
            allocate(NNMAPexciton(NumExcitTI, 10), STAT = ALLOC_ERR )
            IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION NNMAPexciton !!!"
            TIexciton = (0.d0,0.d0)
            NNMAPexciton = 0
            NumExcitTIused = 0

            ! read all excitonic transfer integrals
            call read_localfieldeffects(p_id, gamma, TIexciton, NumExcitTI, NumExcitTIused, NNMAPexciton)
            call read_coulombMatrix(p_id, gamma, TIexciton, NumExcitTI, NumExcitTIused, NNMAPexciton)

            allocate(isUsedTIexciton(NumExcitTIused))
            isUsedTIexciton = .false.

        endif

        do i=1,3
            test_supercell(i) = max(max_shift_i(i), maxS(i))
            test_supercell(i) = ABS(test_supercell(i))
        enddo

        ! if (p_id .eq. 0) then
        !     write(*,*) "Test-Hamiltonian has size: ", test_supercell
        ! endif

        allocate(NNEIGH_per_site(-test_supercell(1):test_supercell(1), &
                                & -test_supercell(2):test_supercell(2), &
                                & -test_supercell(3):test_supercell(3), &
                                &  get_N_con(), get_N_val()))

        NNEIGH_per_site = 0

        if(p_id.eq.0) THEN
            write(*,*) "Estimate number of neighbors (off-diagonals of Hamiltonian)..."
        endif

        ! in the test-Hamiltonian all coordinates are shifted (starting with negative cell)
        ! we cannot use the usual getExcitonIndex functions because they use a different
        ! supercell.
        do Sx1=-test_supercell(1),test_supercell(1)
            if(p_id.eq.0) THEN
                write(*,*) "Progress (%): ", (Sx1+test_supercell(1)) / (test_supercell(1)*2. + 1.)*100.
            endif
            do Sy1=-test_supercell(2),test_supercell(2)
                do Sz1=-test_supercell(3),test_supercell(3)
                    do c1=1,get_N_con()
                        do v1=1,get_N_val()

                            ! estimate max neighbors for single particle part

                            NNEIGH_per_site(Sx1, Sy1, Sz1, c1, v1)=0
                            do ia=1,numTI  ! go over all neighbors (in NNMAP)

                                if((c1.eq.NNMAP(ia,1)).and.(v1.eq.NNMAP(ia,2))) then  ! searching for connections from site c1 to v1

                                    ! get site index for neighbor
                                    c2 = NNMAP(ia,3)
                                    v2 = NNMAP(ia,4)
                                    Sx2 = Sx1 + NNMAP(ia,5)  ! for single particle TIs we have S2-S1 in Hamiltonian
                                    Sy2 = Sy1 + NNMAP(ia,6)
                                    Sz2 = Sz1 + NNMAP(ia,7)

                                    NNEIGH_per_site(Sx1, Sy1, Sz1, c1, v1) = NNEIGH_per_site(Sx1, Sy1, Sz1, c1, v1) +1


                                    ! excitonic TIs that would connect the same sites should be marked as used to prevent double count
                                    if (useTIexciton) then
                                        ! get id in exciton neighbor array
                                        call getNNMAPexcitonId(NNMAPexciton,NumExcitTI, NumExcitTIused, &
                                            & c1,v1,Sx1,Sy1,Sz1, c2,v2,Sx2,Sy2,Sz2, elementId)

                                        if (elementId.gt.0) then
                                            isUsedTIexciton(elementId) = .true.  ! save as used (to prevent double usage)
                                        endif
                                    endif

                                endif  !match c1,v1 index
                            enddo !closes loop over neighbors


                            ! estimate max neighbors for excitonic TIs that do not coincide with single particle TIs

                            if (useTIexciton) then
                                do ia=1,NumExcitTIused  ! go over all neighbors (in NNMAP)

                                    if (isUsedTIexciton(ia)) cycle  ! already set by single particle contributions

                                    if((c1.eq.NNMAPexciton(ia,1)).and.(v1.eq.NNMAPexciton(ia,2)).and.&
                                        &(Sx1.eq.NNMAPexciton(ia,3)).and.(Sy1.eq.NNMAPexciton(ia,4)).and.&
                                        &(Sz1.eq.NNMAPexciton(ia,5))) then  ! only connections from site c1,v1,S1

                                        ! get site index for neighbor
                                        c2 = NNMAPexciton(ia,6)
                                        v2 = NNMAPexciton(ia,7)
                                        Sx2 = NNMAPexciton(ia,8)
                                        Sy2 = NNMAPexciton(ia,9)
                                        Sz2 = NNMAPexciton(ia,10)

                                        if ((c1.eq.c2).and.(v1.eq.v2).and.(Sx1.eq.Sx2).and.(Sy1.eq.Sy2).and.(Sz1.eq.Sz2)) cycle  ! diagonal elements

                                        NNEIGH_per_site(Sx1, Sy1, Sz1, c1, v1) = NNEIGH_per_site(Sx1, Sy1, Sz1, c1, v1) +1

                                    endif  !match c1,v1, S1 index
                                enddo !closes loop over neighbors
                            endif

                        enddo ! v1
                    enddo ! c1
                enddo ! Sz1
            enddo ! Sy1
        enddo ! Sx1


        NNEIGH = MAXVAL(NNEIGH_per_site)  ! output

       
        deallocate(NNMAP)
        deallocate(TIsingle)
        if (allocated(TIexciton)) deallocate(TIexciton)
        if (allocated(NNMAPexciton)) deallocate(NNMAPexciton)
        if (allocated(isUsedTIexciton)) deallocate(isUsedTIexciton)


        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"

        return

    end subroutine

    SUBROUTINE create_exciton_hamiltonian(p_id, sites_p_len, sites_p,neiList,enList,traList, &
        & gamma, config)
        ! This routine sets up the exciton hamiltonian, where electrons and holes can have different
        ! monomentum
        ! This is not the optical hamiltonian

        IMPLICIT NONE
        INTEGER p_id,sites_p_len, sites_p(sites_p_len),ierr,cadd
        type(PARFILE) :: config

        ! Variables which define the potential
        INTEGER ia
        INTEGER globId1, globId2
        INTEGER Sx1, Sy1, Sz1, Sx2, Sy2, Sz2

        integer i, ii,jj, numTI
        character(len=1024) :: filename
        integer max_shift_i(3), max_shift_o(3)
        integer maxS(3)

        ! Variables for the calculation of neighbors
        INTEGER neiList(sites_p_len,2*config%NNEIGH),nbar, nbar_max, ALLOC_ERR, elementId
        integer NprimSites_i, NprimSites_o, NumExcitTI, NumExcitTIused
        DOUBLE PRECISION gamma
        DOUBLE COMPLEX, allocatable  :: TIsingle(:,:,:,:,:,:,:), TIexciton(:)
        INTEGER, allocatable :: NNMAP(:,:)
        INTEGER, allocatable :: NNMAPexciton(:,:)
        DOUBLE COMPLEX traList(sites_p_len,config%NNEIGH), enList(sites_p_len)
        DOUBLE PRECISION, ALLOCATABLE :: onsite_energies_i(:), onsite_energies_o(:)  ! values from the onsite energy file

        INTEGER :: c1,c2,v1,v2, Sx_shifted1, Sy_shifted1, Sz_shifted1, Sx_shifted2, Sy_shifted2, Sz_shifted2
        logical :: useTIexciton
        logical, allocatable ::  isUsedTIexciton(:), isUsedTIexciton_global(:)

        !------------------------ End variable declaration --------------------------

        traList=(0.d0,0.d0)
        neiList=0.d0
        enList=(0.d0,0.d0)
        NumExcitTIused=0

        ! read geometry from POSFILE
        !call generatecoordinates(p_id, ABC, bas_o, bas_i) !all for output

        ! Onsite energy for valence sites
        call read_onsite_energies(p_id, onsite_energies_i, onsite_energies_o, NprimSites_i, NprimSites_o)

        ! get max shift and allocate TI arrays
        if(p_id.eq.0) THEN
            write(*,*) ""
            write(*,*) "Read and extract single particle TIs from TINFILES..."
        endif
        filename="TINFILE_v"
        call get_max_shift_from_tinfile(p_id, filename, max_shift_o)
        filename="TINFILE_c"
        call get_max_shift_from_tinfile(p_id, filename, max_shift_i)

        ! get maximal shift for all single particle TIs
        do i=1,3
            max_shift_i(i) = max(max_shift_i(i), max_shift_o(i))
        enddo

        allocate(TIsingle(get_N_con(), get_N_val(), get_N_con(), get_N_val(), -max_shift_i(1):max_shift_i(1), &
            & -max_shift_i(2):max_shift_i(2),-max_shift_i(3):max_shift_i(3)))

        ! read all single particle TIs and the nearest neighbor map
        call read_tinfile(p_id, numTI, NNMAP, TIsingle, gamma)

        if(p_id.eq.0) THEN
            write(*,*) ""
            write(*,*) "Read and extract excitonic TIs from LOCALFIELDEFFECTS and COULOMB file..."

            ! check duplicates in COULOMB and LOCALFIELDEFFECTS file
            if (config%CHECK_DUPLICATES) THEN
            write(*,*) "Check for duplicates in COULOMB file"            
            filename = "COULOMB"
            call check_duplicates_COULOMB(p_id, filename)

            write(*,*) "Check for duplicates in LOCALFIELDEFFECTS file" 
            filename = "LOCALFIELDEFFECTS"
            call check_duplicates_COULOMB(p_id, filename)
            endif
        endif

        ! read (non-diagonal) Coulomb and local field effects matrix elements similar as TI
        call get_max_shift_interaction_elements(p_id, maxS, NumExcitTI, useTIexciton)


        ! allocate data structure for excitonic transfer integrals
        if (useTIexciton) then
            allocate(TIexciton(NumExcitTI))
            allocate(NNMAPexciton(NumExcitTI, 10), STAT = ALLOC_ERR )
            IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION NNMAPexciton !!!"
            TIexciton = (0.d0,0.d0)
            NNMAPexciton = 0
            NumExcitTIused = 0

            ! read all excitonic transfer integrals
            if (p_id.eq.0) then
            write(*,*) ""
            write(*,*) "Read local field effects"
            endif
            call read_localfieldeffects(p_id, gamma, TIexciton, NumExcitTI, NumExcitTIused, NNMAPexciton)

            if (p_id.eq.0) then
            write(*,*) "Read Coulomb matrix elements"
            endif
            call read_coulombMatrix(p_id, gamma, TIexciton, NumExcitTI, NumExcitTIused, NNMAPexciton)

            ! make excitonic TI hermitian (this is only for testing purposes!!! Do not use it if you dont need it!)
            if (config%FORCE_HERMITIAN) then
            call symmetrize_TIexciton(p_id, TIexciton, NNMAPexciton, NumExcitTI, NumExcitTIused)
            if (p_id.eq.0) write(*,*) "Skip check of hermiticity for exciton TIs."
            else
            ! check if TIexciton is hermitian
            if (p_id.eq.0) then
                write(*,*) "Number of excitonic TIs (combinded LFE and Coulomb):", NumExcitTIused
                write(*,*) "Check hermiticity of exciton TIs..."
                call check_hermicity_TIexciton(TIexciton, NNMAPexciton, NumExcitTI)
            endif
            endif

            allocate(isUsedTIexciton(NumExcitTIused))
            isUsedTIexciton = .false.
        endif

        if (p_id.eq.0) then
            write(*,*) ""
            write(*,*) "Set up optical Exciton Hamiltonian (only vertical excitations!)..."
        endif

        ! Loop to create traList (non-diagonal part of the hamiltonian, supercell) of m only for the corresponding process
        nbar_max = 0
        do cadd=1, sites_p_len  ! go over all sites in the local domain

            ! determine site index for cadd
            globId1 = sites_p(cadd)  ! get global site index from local site index
            call getExcitonIndex(globId1, Sx1,Sy1,Sz1,v1,c1)


            !---------- set onsite energies in the hamiltonian ----------------------------------
            enList(cadd) = enList(cadd) + onsite_energies_i(c1) / gamma
            enList(cadd) = enList(cadd) - onsite_energies_o(v1) / gamma ! minus sign because of exciton hamiltonian


            !---------- set (single particle) transfer integrals in the hamiltonian --------------  
            nbar=0
            do ia=1,numTI  ! go over all neighbors (in NNMAP)

            ! write(*,*) p_id, cadd, ia

            if((c1.eq.NNMAP(ia,1)).and.(v1.eq.NNMAP(ia,2))) then  ! only connections from site c1,v1

                ! get site index for neighbor
                c2 = NNMAP(ia,3)
                v2 = NNMAP(ia,4)
                Sx2 = Sx1 + NNMAP(ia,5)  ! for single particle TIs we have S2-S1 in Hamiltonian
                Sy2 = Sy1 + NNMAP(ia,6)
                Sz2 = Sz1 + NNMAP(ia,7)

                call getGlobIdFromExcitonIndex(Sx2,Sy2,Sz2,v2,c2, globId2)

                if (globId2.gt.0) then
                    nbar=nbar+1
                    if(nbar.gt.config%NNEIGH) then
                        write(*,*) 'Problem with number of neighbors in Hamiltonian:',globId1,ia,nbar,config%NNEIGH
                        STOP
                    endif

                    neiList(cadd,nbar) = globId2   ! set neighbor entry (like in original)
                    neiList(cadd, nbar + config%NNEIGH) = getProcessId(globId2)

                    ! set single particle transfer integral
                    traList(cadd,nbar)=TIsingle(c1,v1,c2,v2, Sx2-Sx1, Sy2-Sy1, Sz2-Sz1)


                    !---------- set (excitonic) transfer integrals in the hamiltonian -------------- 
                    if (useTIexciton) then
                        ! check if also excitonic TI for the very same connection exists
                        call getExcitonIndex_shifted(globId1, Sx_shifted1, Sy_shifted1, Sz_shifted1, v1, c1)  ! it is very important that we use the shifted indexes here!!
                        call getExcitonIndex_shifted(globId2, Sx_shifted2, Sy_shifted2, Sz_shifted2, v2, c2)

                        ! check if S1 and S2 indexes are in the correct range:
                        if ((abs(Sx_shifted1).le.maxS(1)).and.(abs(Sy_shifted1).le.maxS(2)).and.(abs(Sz_shifted1).le.maxS(3))&
                        &.and.(abs(Sx_shifted2).le.maxS(1)).and.(abs(Sy_shifted2).le.maxS(2)).and.&
                        &(abs(Sz_shifted2).le.maxS(3))) THEN

                        ! get id in exciton neighbor array
                        call getNNMAPexcitonId(NNMAPexciton,NumExcitTI, NumExcitTIused, &
                            &c1,v1,Sx_shifted1,Sy_shifted1,Sz_shifted1, &
                            &c2,v2,Sx_shifted2,Sy_shifted2,Sz_shifted2, elementId)


                        if (elementId.gt.0) then
                            traList(cadd,nbar)= traList(cadd,nbar) + TIexciton(elementId) ! set TI
                            isUsedTIexciton(elementId) = .true.  ! save as used (to prevent double usage)

                        endif

                        endif
                    endif

                endif
            endif  !match c1,v1 index
            enddo !closes loop over neighbors


            !---------- set (excitonic) transfer integrals that do not coincide with single particle TIs --------------
            if (useTIexciton) then
            call getExcitonIndex_shifted(globId1, Sx1, Sy1, Sz1, v1, c1)  ! changed S1 to S1_shifted
            do ia=1,NumExcitTIused  ! go over all neighbors (in NNMAP)

                if (isUsedTIexciton(ia)) cycle  ! already set by single particle contributions

                if((c1.eq.NNMAPexciton(ia,1)).and.(v1.eq.NNMAPexciton(ia,2)).and.&
                    &(Sx1.eq.NNMAPexciton(ia,3)).and.(Sy1.eq.NNMAPexciton(ia,4)).and.&
                    &(Sz1.eq.NNMAPexciton(ia,5))) then  ! only connections from site c1,v1,S1

                    ! get site index for neighbor
                    c2 = NNMAPexciton(ia,6)
                    v2 = NNMAPexciton(ia,7)
                    Sx2 = NNMAPexciton(ia,8)
                    Sy2 = NNMAPexciton(ia,9)
                    Sz2 = NNMAPexciton(ia,10)

                    if ((c1.eq.c2).and.(v1.eq.v2).and.(Sx1.eq.Sx2).and.(Sy1.eq.Sy2).and.(Sz1.eq.Sz2)) cycle  ! diagonal elements are set via on-site energies

                    call getGlobIdFromExcitonIndex_shifted(Sx2,Sy2,Sz2,v2,c2, globId2)

                    if (globId2.gt.0) then
                        nbar=nbar+1
                        if(nbar.gt.config%NNEIGH) then
                        write(*,*) 'Problem with number of neighbors in Hamiltonian:',globId1,ia,nbar,config%NNEIGH
                        STOP
                        endif

                        neiList(cadd,nbar) = globId2   ! set neighbor entry (like in original)
                        neiList(cadd, nbar + config%NNEIGH) = getProcessId(globId2)

                        ! set single particle transfer integral
                        traList(cadd,nbar)=TIexciton(ia)
                        isUsedTIexciton(ia) = .true.

                        ! write(*,*) c1,v1,Sx1, Sy1, Sz1, c2,v2,Sx2,Sy2,Sz2, TIexciton(ia), p_id

                    endif
                endif  !match c1,v1 index
            enddo !closes loop over neighbors
            endif
            
            if (nbar.gt.nbar_max) then  ! monitor how many neighbors are actually used
            nbar_max = nbar
            ! if (p_id.eq.0) write(*,*) "new nbar_max = ", nbar_max, nneigh
            endif

        enddo ! loop over all sites in local domain cadd

        ! synchronize all processes
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"

        if (useTIexciton.and.(NumExcitTIused.gt.0)) then
            allocate(isUsedTIexciton_global(NumExcitTIused))
            isUsedTIexciton_global = .false.

            call MPI_Reduce(isUsedTIexciton, isUsedTIexciton_global, NumExcitTIused, &
                &MPI_LOGICAL, MPI_LXOR, 0, MPI_COMM_WORLD, ierr)
                

            if (p_id.eq.0) then
            jj = 0
            do ii=1,NumExcitTIused
                if (.not.isUsedTIexciton_global(ii)) THEN  ! check if some TIs are not used
                    jj = jj+1
                endif
            enddo
            if (jj.gt.0) THEN  ! at least one TI was not used or double used by different processs
                write(*,*) ""
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) jj , "excitonic TIs are not used or double used!!!"
                write(*,*) ""
                do ii=1,NumExcitTIused
                    if (.not.isUsedTIexciton_global(ii)) THEN
                        write(*,*) NNMAPexciton(ii,:), TIexciton(ii) / 2. * gamma
                    endif
                enddo
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) ""
            else
                write(*,*) "All excitonic TIs are used :D"
            endif
            endif
        endif

        call MPI_Reduce(nbar_max, jj, 1, &
            &MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

        if (p_id.eq.0)  then
            write(*,*) "Maximal number of neighbors used = ", jj, "from possible nneigh=", config%NNEIGH

            if (jj.ne.config%NNEIGH) THEN
            write(*,*) ""
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) " You are not using the optimal number of neighbors!"
                write(*,*) ""
                write(*,*) " Please set     NNEIGH = ", jj
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) ""
            endif

        endif

        deallocate(onsite_energies_i)
        deallocate(onsite_energies_o)
        deallocate(NNMAP)
        deallocate(TIsingle)
        if (allocated(TIexciton)) deallocate(TIexciton)
        if (allocated(NNMAPexciton)) deallocate(NNMAPexciton)
        if (allocated(isUsedTIexciton)) deallocate(isUsedTIexciton)
        if (allocated(isUsedTIexciton_global)) deallocate(isUsedTIexciton_global)

        ! density-density interaction via onsite energies
        call coulomb_density_density(sites_p_len, p_id, sites_p, enList, gamma)

        ! add local field effects for S=S'
        call localFieldEffects_onsite(sites_p_len, p_id, sites_p, enList, gamma)

        ! call plotDiagonalHamiltonian(p_id, sites_p_len, sites_p, enList)

        ! add some on-site (local) anderson disorder
        call create_local_disorder(p_id, sites_p_len, enList, gamma, config)


        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"

        ! For debugging of the onsite energies.
        ! call getGlobIdFromExcitonIndex(76,1,1,1,1,globId1)
        ! call getExcitonIndex_shifted(globId1, Sx1,Sy1,Sz1,m1,m2)
        ! if (p_id.eq.0) then
        !    write(*,*) "DEBUG: ", globId1, Sx1,Sy1,Sz1,m1,m2
        ! endif

        ! pId = getProcessId(globId1)
        ! if ((p_id+1).eq.pId) THEN
        !    locId = getLocalId(globId1)
        !    write(*,*) "DEBUG: globId, locId, enList", globId1, locId, enList(locId)
        ! endif

        ! call getGlobIdFromExcitonIndex(76,151,151,1,1,globId1)
        ! call getExcitonIndex_shifted(globId1, Sx1,Sy1,Sz1,m1,m2)
        ! if (p_id.eq.0) then
        !    write(*,*) "DEBUG: ", globId1, Sx1,Sy1,Sz1,m1,m2
        ! endif

        ! pId = getProcessId(globId1)
        ! if ((p_id+1).eq.pId) THEN
        !    locId = getLocalId(globId1)
        !    write(*,*) "DEBUG: globId, locId, enList", globId1, locId, enList(locId)
        ! endif

        ! if(p_id.eq.0) then
        !    write(*,*) "connections of globId=", sites_p(100)
        !    call getExcitonIndex(sites_p(100), Sx1,Sy1,Sz1,m1,m2)
        !    write(*,*) sites_p(100), getProcessId(sites_p(100)), ": ", Sx1,Sy1,Sz1,m1,m2
        !    do ia=1,nneigh
        !       call getExcitonIndex(neiList(100,ia), Sx1,Sy1,Sz1,m1,m2)
        !       write(*,*) traList(100,ia), neiList(100,ia), neiList(100,ia+nneigh), ": ", Sx1,Sy1,Sz1,m1,m2
        !    enddo
        ! endif

        
        if(p_id.eq.0) then
            write(*,*) ""
            write(*,*) 'Hamiltonian finished.'
            write(*,*) "================================================================"
            write(*,*) ""
        endif

        ! call save_onsiteEnergies(sites_p_len, sites_p, p_id, enList, gamma)
        ! call save_alltransferintegrals(sites_p_len,sites_p,neiList, p_id, traList, gamma)

        return

    END SUBROUTINE create_exciton_hamiltonian


    SUBROUTINE create_local_disorder(p_id,sites_p_len, enList, gamma, config)
        ! Reads DISFILE and creates the corresponding on-site (local) disorder
        ! Input : p_id             = id of process
        !         sites_p_len          = number of site in local domain
        !         enList                  = onsite energies for local domain
        !         gamma               = energy scale
        !         config              = parameters from PARFILE
        implicit none
        type(PARFILE) :: config
        integer p_id, sites_p_len !, sites_p(sites_p_len)
        DOUBLE COMPLEX enList(sites_p_len)
        double precision disorder_strength, gamma

        intent(in) p_id,sites_p_len, gamma, config
        intent(inout) enList


        if(TRIM(config%DISORDER_MODEL).eq.'none') then
            if (p_id.eq.0) write(*,*) "Disorder model = none --> calculation without disorder!"
            return
        else if(TRIM(config%DISORDER_MODEL).eq.'ANDERSON') then
            disorder_strength = config%DISORDER_STRENGTH/gamma
            if (abs(disorder_strength) .gt. 1e-10) then
            call generatepot_anderson(sites_p_len, p_id, enList, disorder_strength)
            else
            if (p_id.eq.0) write(*,*) "Disorder strength = 0 --> calculation without disorder!"
            return
            endif
        else
            stop "Disorder model is not known!"
        endif

    END SUBROUTINE


    SUBROUTINE read_tinfile(p_id, numTI, NNMAP, TIsingle, gamma)
        ! Reads transfer integrals from TINFILE
        ! Input :   p_id  = id of process
        ! Ouput :   NAT1     = number of entries in the TINFILE
        !           NNMPA    = connections as stored in TINFILE (first 5 columns)
        !           TIouter  = values of outer transfer integrals
        !           TIinner  = values of inner transfer integrals
        !           gamma    = internal energy unit (usually largest absolute TI)
        USE Geometry
        implicit none
        integer p_id, NAT1outer, NAT1inner, i,j,k, l,ll,ii,jj,kk, numTI
        integer max_shift(3), NBATinner, NBATouter
        CHARACTER dummy


        DOUBLE PRECISION val1,val2, gamma
        INTEGER, allocatable :: NNMAPouter(:,:), NNMAPinner(:,:), NNMAP(:,:)
        DOUBLE COMPLEX, allocatable :: TIouter(:), TIinner(:)
        integer c1,c2,v1,v2,Sx,Sy,Sz
        logical isContained

        DOUBLE COMPLEX, allocatable :: TIsingle(:, :, :, :, :, :, :)


        intent(in) p_id
        intent (inout) TIsingle
        intent(out) NNMAP, gamma, numTI

        ! get dimensions of the array
        NBATinner = size(TIsingle,1)
        NBATouter = size(TIsingle,2)
        max_shift(1) = ubound(TIsingle,5)
        max_shift(2) = ubound(TIsingle,6)
        max_shift(3) = ubound(TIsingle,7)


        if (p_id.eq.0) then
            write(*,*) ""
            write(*,*) "Read valence (outer) TINFILE ..."
        endif

        open(13,FILE='TINFILE_v',ACTION='READ',ERR=140)
        read(13,*,ERR=141) NAT1outer
        read(13,*,ERR=142) dummy
        if(p_id.eq.0) write(*,*) '--> Number of valence TI = ', NAT1outer

        if (mod(NAT1outer,2).ne.0) THEN
            stop "Number of transfer integrals has to be even!&
            & Please check if the complex conjugate is present for each TI"
        endif

        ! read nearest neigbor map (NNMAP) from TINFILE and determine gamma (internal energy scale)
        allocate(NNMAPouter(NAT1outer,5))
        allocate(TIouter(NAT1outer))
        TIouter = (0.d0,0.d0)
        NNMAPouter = 0
        gamma=0.d0
        do i=1,NAT1outer
            read(13,*,ERR=144) l,ll,ii,jj,kk,val1, val2

            if(l.ge.1.and.l.le.NBATouter.and.ll.ge.1.and.ll.le.NBATouter &
            & .and.ii.ge.-max_shift(1).and.ii.le.max_shift(1) &
            & .and.jj.ge.-max_shift(2).and.jj.le.max_shift(2) &
            & .and.kk.ge.-max_shift(3).and.kk.le.max_shift(3))   then

            NNMAPouter(i,1) = l
            NNMAPouter(i,2) = ll
            NNMAPouter(i,3) = ii
            NNMAPouter(i,4) = jj
            NNMAPouter(i,5) = kk
            TIouter(i) = CMPLX(val1, val2)

            gamma=MAX(gamma, abs(TIouter(i)))
            else
            write(*,*) 'Error reading TINFILE_v in line',2+i
            STOP
            endif
        enddo
        close(13)


        if (p_id.eq.0) then
            write(*,*) "Read conduction (inner) TINFILE ..."
        endif
        open(13,FILE='TINFILE_c',ACTION='READ',ERR=140)
        read(13,*,ERR=141) NAT1inner
        read(13,*,ERR=142) dummy
        if(p_id.eq.0) write(*,*) '--> Number of conduciton TI = ', NAT1inner

        if (mod(NAT1inner,2).ne.0) THEN
            stop "Number of transfer integrals has to be even!&
            & Please check if the complex conjugate is present for each TI"
        endif

        ! read nearest neigbor map (NNMAP) from TINFILE and determine gamma (internal energy scale)
        allocate(NNMAPinner(NAT1inner, 5))
        allocate(TIinner(NAT1inner))
        TIinner = (0.d0,0.d0)
        NNMAPinner = 0
        do i=1,NAT1inner
            read(13,*,ERR=144) l,ll,ii,jj,kk,val1, val2
            if(l.ge.1.and.l.le.NBATinner.and.ll.ge.1.and.ll.le.NBATinner &
            & .and.ii.ge.-max_shift(1).and.ii.le.max_shift(1) &
            & .and.jj.ge.-max_shift(2).and.jj.le.max_shift(2) &
            & .and.kk.ge.-max_shift(3).and.kk.le.max_shift(3)) then

            NNMAPinner(i,1) = l
            NNMAPinner(i,2) = ll
            NNMAPinner(i,3) = ii
            NNMAPinner(i,4) = jj
            NNMAPinner(i,5) = kk
            TIinner(i) = CMPLX(val1, val2)

            gamma=MAX(gamma, abs(TIinner(i)))
            else
            write(*,*) l,ll,ii,jj,kk, TIinner(i)
            write(*,*) 'Error reading TINFILE_c in line',2+i
            STOP
            endif
        enddo
        close(13)

        !gamma = 1.0 ! for debugging

        if (abs(gamma).lt.1e-10) then
            write(*,*) "Set gamma to one (it was zero before)."
            gamma = 1.0
        endif


        if(p_id.eq.0) write(*,*) 'Internal energy scale:  gamma = ',gamma


        ! create TIsingel data structure
        ! if(p_id.eq.0) write(*,*) 'Create data structrue (inner)'
        TIsingle=(0.d0,0.d0)
        numTI = NAT1inner*NBATouter + NAT1outer*NBATinner
        if(p_id.eq.0) write(*,*) '--> Estimated number of single-particle TI = ', numTI
        allocate(NNMAP(numTI, 7))

        ! inner (conduction) TI
        j=1
        do i=1,NAT1inner
            do v1=1,NBATouter

            v2 = v1  ! delta(v1-v2) in hamiltonian
            c1 = NNMAPinner(i,1)
            c2 = NNMAPinner(i,2)
            Sx = NNMAPinner(i,3)
            Sy = NNMAPinner(i,4)
            Sz = NNMAPinner(i,5)

            TIsingle(c1,v1,c2,v2,Sx,Sy,Sz) = TIinner(i) / gamma

            ! set NNMAP
            NNMAP(j,1) = c1
            NNMAP(j,2) = v1
            NNMAP(j,3) = c2
            NNMAP(j,4) = v2
            NNMAP(j,5) = Sx
            NNMAP(j,6) = Sy
            NNMAP(j,7) = Sz

            j=j+1

            enddo
        enddo

        ! outer (valence) TI
        ! if(p_id.eq.0) write(*,*) 'Create data structrue (outer)'
        do i=1,NAT1outer
            do c1=1,NBATinner

            c2 = c1  ! delta(c1-c2) in hamiltonian
            v2 = NNMAPouter(i,1)  ! v1 and v2 are exchanged, see Eq. (12) in J. Phys. Mater. 7 (2024) 015001
            v1 = NNMAPouter(i,2)
            Sx = NNMAPouter(i,3)
            Sy = NNMAPouter(i,4)
            Sz = NNMAPouter(i,5)

            TIsingle(c1,v1,c2,v2,Sx,Sy,Sz) = TIsingle(c1,v1,c2,v2,Sx,Sy,Sz) - TIouter(i) / gamma
            

            ! set NNMAP if not already contained
            isContained = .false.
            do k=1,j
                if ((NNMAP(k,1).eq.c1).and.(NNMAP(k,2).eq.v1).and.(NNMAP(k,3).eq.c2).and.(NNMAP(k,4).eq.v2).and.&
                    &(NNMAP(k,5).eq.Sx).and.(NNMAP(k,6).eq.Sy).and.(NNMAP(k,7).eq.Sz)) THEN
                        isContained = .true.
                        exit
                    endif
            enddo

            if (.not.isContained) then
                NNMAP(j,1) = c1
                NNMAP(j,2) = v1
                NNMAP(j,3) = c2
                NNMAP(j,4) = v2
                NNMAP(j,5) = Sx
                NNMAP(j,6) = Sy
                NNMAP(j,7) = Sz

                j=j+1
            endif

            enddo
        enddo

        if(p_id.eq.0) write(*,*) '--> Actual number of single-particle TI    = ', j-1

        ! deallocate(NNMAP(j:numTI,:))
        numTI = j-1
        deallocate(NNMAPouter)
        deallocate(NNMAPinner)
        deallocate(TIinner)
        deallocate(TIouter)

        !debugging
        ! if (p_id.eq.0) then
        !    do i=1,numTI

        !       c1 = NNMAP(i,1)
        !       v1 = NNMAP(i,2)
        !       c2 = NNMAP(i,3)
        !       v2 = NNMAP(i,4)
        !       Sx = NNMAP(i,5)
        !       Sy = NNMAP(i,6)
        !       Sz = NNMAP(i,7)

        !       write(*,*) i, ":", c1,v1,c2,v2,Sx,Sy,Sz, TIsingle(c1,v1,c2,v2,Sx,Sy,Sz)*gamma

        !    enddo
        ! endif



        ! check if hermitian conjugate of TI is also present for each TI
        if (p_id.eq.0) then
            write(*,*) ""
            write(*,*) "Check hermiticity for TIsingle..."
            call check_hermicity_TIsingle(TIsingle)
        endif



        GOTO 149

140   write(*,*) 'Could not read TINFILE 0'
        STOP
141   write(*,*) 'Could not read TINFILE 1'
        STOP
142   write(*,*) 'Could not read TINFILE 2'
        STOP
144   write(*,*) 'Could not read TINFILE 3'
        STOP
149   CONTINUE


        ! TODO check if complex conjugate is present!!! --> Hermiticity

        return

    END SUBROUTINE read_tinfile

    subroutine check_hermicity_TIsingle(TIsingle)
        ! Checks if transferintegrals that are given in the TINFILE produce
        ! a hermitian Hamiltonian.
        ! It does not check the actual Hamiltonian !!!
        USE Geometry
        
        DOUBLE COMPLEX, intent(in), allocatable :: TIsingle(:, :, :, :, :, :, :)
        integer max_shift(3), NBATinner, NBATouter
        double complex val1, val2
        integer ii,jj,kk, v1,v2,c1,c2

        ! get dimensions of the array
        NBATinner = size(TIsingle,1)
        NBATouter = size(TIsingle,2)
        max_shift(1) = ubound(TIsingle,5)
        max_shift(2) = ubound(TIsingle,6)
        max_shift(3) = ubound(TIsingle,7)

        ! RESHAPE(TIsingle,(/ NBATinner, NBATouter,NBATinner, NBATouter,-max_shift(1):max_shift(1),&
        !    & -max_shift(2):max_shift(2),-max_shift(3):max_shift(3) /) )


        do c1=1,NBATinner
            do v1=1,NBATouter
            do c2=1,NBATinner
                do v2=1,NBATouter
                    do ii=0,max_shift(1)
                        do jj=0,max_shift(2)
                        do kk=0,max_shift(3)
                        
                            val1 = TIsingle(c1,v1,c2,v2,ii,jj,kk)
                            val2 = CONJG(TIsingle(c2,v2,c1,v1,-ii,-jj,-kk))

                            if ((abs(real(val1)-real(val2)).gt. 1e-10) .or. &
                                &(abs(imag(val1)-imag(val2)).gt. 1e-10)) THEN
                                write(*,*) "TIsingle do not lead to hermitian Hamiltonian: "
                                write(*,*) c1,v1,c2,v2,ii,jj,kk, val1
                                write(*,*) c2,v2,c1,v1,-ii,-jj,-kk, val2

                                stop
                                
                            endif

                        enddo
                        enddo
                    enddo
                enddo
            enddo
            enddo
        enddo
    

        write(*,*) "--> single particle TIs are hermitian :D"


    end subroutine



    SUBROUTINE get_max_shift_from_tinfile(p_id, filename, max_shift)
        ! Reads TINFILE and determines the larges absolute shift for each dimension.
        ! This is necessary to allocate large enough TIinner/ TIouter arrays.
        ! Input :   p_id   = id of process
        !           filename  = path to the TINFILE
        ! Ouput :   max_shift = 3-dim array of maximum shifts 
        implicit none
        integer max_shift(3)
        integer p_id, NAT1, i, l,ll,ii,jj,kk
        CHARACTER dummy
        character(len=1024) :: filename
        DOUBLE PRECISION val1,val2

        intent(in) p_id, filename
        intent(out) max_shift

        max_shift = 0

        !if (p_id == 0) write(*,*) 'Determine maximum shift in each direction for file', filename
        open(13,FILE=filename,ACTION='READ',ERR=340)
        read(13,*,ERR=341) NAT1
        read(13,*,ERR=342) dummy
        ! if(p_id.eq.0) write(*,*) 'Number of TI = ', NAT1

        ! 
        do i=1,NAT1
            read(13,*,ERR=344) l,ll,ii,jj,kk,val1, val2
            
            max_shift(1) = max(max_shift(1), abs(ii))
            max_shift(2) = max(max_shift(2), abs(jj))
            max_shift(3) = max(max_shift(3), abs(kk))

        enddo
        close(13)

        if(p_id.eq.0) write(*,*) 'Maximum absolute shifts = ', max_shift

        GOTO 349

340   write(*,*) 'Could not read TINFILE 0'
        STOP
341   write(*,*) 'Could not read TINFILE 1'
        STOP
342   write(*,*) 'Could not read TINFILE 2'
        STOP
344   write(*,*) 'Could not read TINFILE 3'
        STOP
349   CONTINUE

        return

    END SUBROUTINE get_max_shift_from_tinfile


    SUBROUTINE read_onsite_energies(p_id, onsite_energies_i, onsite_energies_o, NprimSites_i, NprimSites_o)
        USE Geometry
        implicit none
        DOUBLE PRECISION, ALLOCATABLE :: onsite_energies_i(:), onsite_energies_o(:)  ! values from the onsite energy file
        logical LFI_ONS
        integer NprimSites_i, NprimSites_o, p_id, i

        intent(in) p_id
        intent(out) onsite_energies_i, onsite_energies_o

        if (p_id.eq.0) then
            write(*,*) ""
            write(*,*) "Read onsite energies ..."
        endif

        inquire(FILE='./ONSITE_ENERGY_c',exist=LFI_ONS)  ! TODO change to inner and outer coordinates
        if(.not.LFI_ONS) then
            write(*,*) './ONSITE_ENERGY_c does not exist. Stop.'
            stop
        endif
        inquire(FILE='./ONSITE_ENERGY_v',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) './ONSITE_ENERGY_v does not exist. Stop.'
            stop
        endif

        open(123,FILE='./ONSITE_ENERGY_c',ACTION='READ')
        open(124,FILE='./ONSITE_ENERGY_v',ACTION='READ')
        
        read(123,*) NprimSites_i
        read(124,*) NprimSites_o

        ! if (NprimSites_i.ne.NprimSites_o) then
        !    write(*,*) "Number of onsite energies for inner and outer sites needs to be the same."
        !    stop
        ! endif
        if (NprimSites_o.ne.get_N_val()) then
            write(*,*) "Different number of sites in ONSITE_ENERGY_v file and POSFILE"
            stop
        endif
        if (NprimSites_i.ne.get_N_con()) then
            write(*,*) "Different number of sites in ONSITE_ENERGY_c file and POSFILE"
            stop
        endif

        ALLOCATE(onsite_energies_i(1:NprimSites_i))
        ALLOCATE(onsite_energies_o(1:NprimSites_o))
        do i=1,NprimSites_i
            read(123,*) onsite_energies_i(i)
        enddo
        do i=1,NprimSites_o
            read(124,*) onsite_energies_o(i)
        enddo

        close(123)
        close(124)

        ! print onsite energies
        if (p_id.eq.0) then
            write(*,*) 'Onsite energies for valence (outer) sites:'
            do i=1,NprimSites_o
            write(*,*) i, onsite_energies_o(i)
            enddo
            write(*,*) 'Onsite energies for conduction (inner) sites:'
            do i=1,NprimSites_i
            write(*,*) i, onsite_energies_i(i)
            enddo
        endif

        return
    
    END SUBROUTINE read_onsite_energies

! #ifdef DEBUG
!       subroutine hamilt_hermitianq(neiList,traList,sites_p_len,sites_p,p_id)
!          USE Geometry
!          IMPLICIT NONE
!          INTEGER i,m,n,j,jj,sites_p_len,p_id
!          INTEGER neiList(sites_p_len,2*nneigh),sites_p(sites_p_len)
!          DOUBLE COMPLEX traList(sites_p_len,nneigh)
!          LOGICAL HERMQ

!          HERMQ=.TRUE.
    
!          do i=1,nneigh
!             do m=1,sites_p_len

!                if(neiList(m,i).ne.0) then
!                   !if(neiList(m,nneigh).ne.neiList(m,2*nneigh/3))then
!                   if(neiList(m,i+nneigh)-1.eq.p_id) then 
!                      n=neiList(m,i)
            
!                      do j=1,nneigh
!                         if(neiList(n,j).eq.m) jj=j
!                      enddo

!                      if(DIMAG(traList(m,i))-DIMAG(conjg(traList(n,jj))).ge.1.0d-10.and.dble(traList(m,i))-dble(conjg(traList(n,jj))).ge.1.0d-10) then
!                         write(*,*) m,n,i,traList(m,i),conjg(traList(n,jj))
!                         HERMQ=.FALSE.
!                      endif
!                   endif
!                endif

!             enddo
!          enddo
        
!          write(*,*) 'Hamiltonian HermitianQ=',HERMQ
        
!       END SUBROUTINE hamilt_hermitianq
! #endif


#ifdef DEBUG
    subroutine plotDiagonalHamiltonian(p_id, sites_p_len, sites_p, enList)
        ! Prints the diagonal part of the Hamiltonian to a file for each
        ! process. This is for debugging purposes only.
        USE Geometry
        IMPLICIT NONE
        INTEGER i, p_id, sites_p_len, sites_p(sites_p_len)
        INTEGER Sx,Sy,Sz,lo,li, globId
        DOUBLE COMPLEX enList(sites_p_len)
        character(len=1024) :: filename
        intent(in) p_id, sites_p_len, enList, sites_p

        write (filename, "(A20,I2)") "output/debug_enList.dat_", p_id

        print *, trim(filename)

        open(1, file =filename, status = 'REPLACE')
        write(1,*) "# outerId     innerId     Sx    Sy    Sz    value"
        do i=1,sites_p_len
            globId = sites_p(i)
            call getExcitonIndex_shifted(globId, Sx, Sy, Sz, lo, li)
            write(1,*) lo, li, Sx, Sy, Sz, REAL(enList(i)), IMAG(enList(i))
        end do  
        
        close(1)


    end subroutine plotDiagonalHamiltonian
#endif


    subroutine get_max_shift_interaction_elements(p_id, maxS, NumExcitTI, useTIexciton)
        implicit none
        integer p_id
        integer maxS(3), minS(3)
        INTEGER NumMatEl1, NumMatEl2, stat, NumExcitTI
        character (len=500) :: line
        DOUBLE PRECISION WW
        DOUBLE PRECISION max_dist_coulomb, val1, val2
        INTEGER :: c1,c2,v1,v2
        INTEGER :: RD(3), Rc(3), Rv(3), S1(3), S2(3), ii, jj, kk
        LOGICAL LFI_ONS, useTIexciton

        intent(in) p_id
        intent(out) maxS, useTIexciton, NumExcitTI

        maxS = 0
        minS = 0
        NumExcitTI = 0
        useTIexciton = .false.


        !----------------- READ LOCALFIELDEFFECTS file  -------------------------------
        inquire(FILE='LOCALFIELDEFFECTS',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'LOCALFIELDEFFECTS does not exist. Stop.'
            stop
        endif
        open(114,FILE='LOCALFIELDEFFECTS',ACTION='READ')
        read(114,*) NumMatEl1, max_dist_coulomb, WW

        if (abs(WW) > 1e-8) THEN

            read_loop: do
            read (114, '(A)', iostat=stat)  line
    
            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
    
            ! store data in arrays
            read (line, *, iostat=stat) c1, c2, v1, v2, (RD(ii), ii=1,3), (S1(jj), jj=1,3), (S2(kk), kk=1,3), val1,val2

            ! ignore diagonal elements
            if ((c1.eq.c2).and.(v1.eq.v2).and.(S1(1).eq.S2(1)).and.(S1(2).eq.S2(2)).and.(S1(3).eq.S2(3))) cycle

            NumExcitTI = NumExcitTI+1
            do ii=1,3
                minS(ii) = min(minS(ii), S1(ii), S2(ii))
                maxS(ii) = max(maxS(ii), S1(ii), S2(ii))
            enddo
    
            end do read_loop
        else
            if (p_id.eq.0) write(*,*) "skip local field effects in get_max_shift_interaction_elements"
        endif
        close(114)

        if (p_id.eq.0) THEN
            write(*,*) "minS = ", minS
            write(*,*) "maxS = ", maxS
        endif

        ! do i=1,3
        !    if (.not.(minS(i).eq.-maxS(i))) stop "local field effects are not hermitian!"
        ! enddo

        useTIexciton = NumExcitTI.gt.0

        !-------------------------- READ COULOMB file  --------------------------------
        inquire(FILE='COULOMB',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'COULOMB does not exist. Stop.'
            stop
        endif
        open(114,FILE='COULOMB',ACTION='READ')
        read(114,*) NumMatEl2, max_dist_coulomb, WW

        if (abs(WW) < 1e-8) THEN
            if (p_id.eq.0) write(*,*) "WW = 0 --> without (screened) Coulomb interaction"
            close(114)
            return
        endif

        read_loop2: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop2             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop2          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop2  ! skip comment lines

            ! store data in arrays
            read (line, *, iostat=stat) c1, c2, v1, v2, (RD(ii), ii=1,3), (Rc(jj), jj=1,3), (Rv(kk), kk=1,3), val1,val2

            do ii=1,3
            S1(ii) = -RD(ii)
            S2(ii) = Rc(ii)-RD(ii)-Rv(ii)
            enddo

            ! ignore diagonal elements
            if ((c1.eq.c2).and.(v1.eq.v2).and.(S1(1).eq.S2(1)).and.(S1(2).eq.S2(2)).and.(S1(3).eq.S2(3))) cycle

            NumExcitTI = NumExcitTI+1
            do ii=1,3
            minS(ii) = min(minS(ii), S1(ii), S2(ii))
            maxS(ii) = max(maxS(ii), S1(ii), S2(ii))
            enddo

        end do read_loop2
        close(114)

        if (p_id.eq.0) THEN
            write(*,*) "minS = ", minS
            write(*,*) "maxS = ", maxS
        endif

        ! do i=1,3
        !    if (.not.(minS(i).eq.-maxS(i))) stop "local field effects are not hermitian!"
        ! enddo

        useTIexciton = NumExcitTI.gt.0

    end subroutine get_max_shift_interaction_elements

    subroutine read_localfieldeffects(p_id, gamma, TIexciton, NumMatEl, NumMatElused, NNMAPexciton)
        ! Reads Local field effects file and stores non-diagonal elements
        ! This subroutine overwrites existing entries in NNMAP and TIexciton. It has to be executed before
        ! read_coulombMatrix.
        implicit none
        integer p_id, NNMAPexciton(NumMatEl,10)
        INTEGER :: c1,c2,v1,v2,Sx1,Sx2,Sy1,Sy2,Sz1,Sz2, dummy1, dummy2, dummy3
        double complex :: TIexciton(NumMatEl)
        INTEGER i, NumMatEl, NumMatElused, NumMatElfile, stat
        character (len=500) :: line
        DOUBLE PRECISION WW
        DOUBLE PRECISION max_dist_coulomb, val1, val2, gamma
        
        LOGICAL LFI_ONS

        intent(in) p_id, gamma, NumMatEl
        intent(inout) TIexciton, NNMAPexciton
        intent(out) NumMatElused


        ! READ LOCALFIELDEFFECTS file
        inquire(FILE='LOCALFIELDEFFECTS',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'LOCALFIELDEFFECTS does not exist. Stop.'
            stop
        endif
        open(114,FILE='LOCALFIELDEFFECTS',ACTION='READ')
        read(114,*) NumMatElfile, max_dist_coulomb, WW


        if (abs(WW).lt.1e-10) then
            IF ( p_id .eq. 0 ) write(*,*) "WW_FE = 0 --> calculation without local field effects"
            close(114)
            return
        endif

        if(p_id.eq.0) then
            write(*,*) 'Generate local field effects: WW_FE = ', WW
            if (abs(WW-1.0).gt.1e-8) then 
                write(*,*) '[WARNING] Local field effects should us the bare interaction WW_FE==1 without any screening!!!'
            endif
        endif

        i = 1
        read_loop: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
            if (i.gt.NumMatElfile) STOP "More elements in LOCALFIELDEFFECTS file than expected. "&
                &"Please check the first line of the file!"

            ! store data in arrays
            read (line, *, iostat=stat) c1, c2, v1, v2, dummy1, dummy2, dummy3, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2, val1,val2

            if ( stat /= 0 ) then
                write(*,*) "line ", i
                STOP "Error while reading local field effect matrix element from line."
            endif

            ! ignore diagonal elementsmaxS
            if ((c1.eq.c2).and.(v1.eq.v2).and.(Sx1.eq.Sx2).and.(Sy1.eq.Sy2).and.(Sz1.eq.Sz2)) cycle

            if (i.gt.NumMatEl) stop "Number of excitionic TIs (LFE) is larger than allocation!"

            ! save nearest neighbor list (may overwrite existing values)
            NNMAPexciton(i,1) = c1
            NNMAPexciton(i,2) = v1
            NNMAPexciton(i,3) = Sx1
            NNMAPexciton(i,4) = Sy1
            NNMAPexciton(i,5) = Sz1
            NNMAPexciton(i,6) = c2
            NNMAPexciton(i,7) = v2
            NNMAPexciton(i,8) = Sx2
            NNMAPexciton(i,9) = Sy2
            NNMAPexciton(i,10) = Sz2

            TIexciton(i) = 2* cmplx(val1, val2) * WW / gamma  ! factor 2 because LFE have this factor in the hamiltonian
            i=i+1

        end do read_loop
        close(114)

        NumMatElused = i-1
        if (p_id.eq.0) write(*,*) "Found ", NumMatElused, "non-diagonal matrix elements. Done reading LOCALFIELDEFFECTS file."

    end subroutine read_localfieldeffects


    subroutine read_coulombMatrix(p_id, gamma, TIexciton, NumNNMAP, NumNNMAPused, NNMAPexciton)
        ! Reads COULOMB file and stores non-diagonal elements
        ! This routine appends existing elements in the NNMAP and TIexciton array. If a connection
        ! already exists then the new value is added without creating doublicate connections in NNMAP.
        implicit none
        integer p_id, NNMAPexciton(NumNNMAP,10), NumNNMAPused, elementId
        INTEGER :: c1,c2,v1,v2, i
        double complex :: TIexciton(NumNNMAP)
        INTEGER NumMatEl, stat, ii, jj, kk, RD(3), Rc(3), Rv(3), S1(3), S2(3), NumNNMAP
        character (len=500) :: line
        DOUBLE PRECISION WW
        DOUBLE PRECISION max_dist_coulomb, val1, val2, gamma
        LOGICAL LFI_ONS

        intent(in) p_id, gamma, NumNNMAP
        intent(inout) TIexciton, NumNNMAPused, NNMAPexciton


        ! READ COULOMB file
        inquire(FILE='COULOMB',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'COULOMB does not exist. Stop.'
            stop
        endif
        open(114,FILE='COULOMB',ACTION='READ')
        read(114,*) NumMatEl, max_dist_coulomb, WW

        if (abs(WW).lt.1e-10) then
            IF ( p_id .eq. 0 ) write(*,*) "WW = 0 --> calculation without (screened) Coulomb interaction"
            close(114)
            return
        endif

        if(p_id.eq.0) write(*,*) 'Generate (screened) Coulomb interaction: WW = ', WW

        IF ( p_id .eq. 0 ) write(*,*) "matrix elements with further distance than ", max_dist_coulomb, &
            & "will be ignored."

        i = 1
        read_loop: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
            if (i.gt.NumMatEl) STOP "More elements in COULOMB file than expected. "&
                &"Please check the first line of the file!"

            ! store data in arrays
            read (line, *, iostat=stat) c1, c2, v1, v2, (RD(ii), ii=1,3), (Rc(jj), jj=1,3), (Rv(kk), kk=1,3), val1,val2

            if ( stat /= 0 ) then
                write(*,*) "line ", i
                STOP "Error while reading (screened) Coulomb interaction matrix element from line."
            endif
            i=i+1  ! only count lines (not used matrix elements like for LFE)

            do ii=1,3
            S1(ii) = -RD(ii)
            S2(ii) = Rc(ii)-RD(ii)-Rv(ii)
            enddo

            ! ignore diagonal elements
            if ((c1.eq.c2).and.(v1.eq.v2).and.(S1(1).eq.S2(1)).and.(S1(2).eq.S2(2)).and.(S1(3).eq.S2(3))) cycle

            ! save nearest neighbor list (also contains onsite terms!)
            call appendNNMAPexciton(NNMAPexciton, NumNNMAP, NumNNMAPused, c1,v1,S1(1),S1(2),S1(3), &
            &c2,v2,S2(1),S2(2),S2(3), elementId)

            ! write(*,*) c1,v1,S1(1),S1(2),S1(3), c2,v2,S2(1),S2(2),S2(3), cmplx(val1, val2)

            TIexciton(elementId) = TIexciton(elementId) - cmplx(val1, val2)* WW / gamma  ! minus because of hamiltonian

        end do read_loop
        close(114)

        if ((i-1).ne.NumMatEl) STOP "Different number of elements in COULOMB file than expected. &
            &Please check the first line of the file!"

        if (p_id.eq.0) write(*,*) "Found ", i-1, "matrix elements. Done reading COULOMB file."
    end subroutine read_coulombMatrix


    subroutine getNNMAPexcitonId(NNMAPexciton,NumMatEl, NumMatElused, c1,v1,Sx1,Sy1,Sz1, c2,v2,Sx2,Sy2,Sz2, elementId)
        ! Returns the index of the connection in the nearest neighbor array
        ! If connection is not contained it returns -1
        implicit none
        INTEGER :: NNMAPexciton(NumMatEl,10), elementId
        INTEGER :: ii, NumMatEl, c1,c2,v1,v2, Sx1,Sy1,Sz1, Sx2,Sy2,Sz2, NumMatElused

        intent(in) c1,c2,v1,v2, Sx1,Sy1,Sz1, Sx2,Sy2,Sz2
        intent(in) NNMAPexciton, NumMatElused
        intent(out) elementId

        elementId = -1  ! returns -1 if element is not contained

        do ii=1,NumMatElused  ! only go through elements that are set

            if ((c1.eq.NNMAPexciton(ii,1)).and.(v1.eq.NNMAPexciton(ii,2)).and.&
            &(Sx1.eq.NNMAPexciton(ii,3)).and.(Sy1.eq.NNMAPexciton(ii,4))&
            &.and.(Sz1.eq.NNMAPexciton(ii,5)).and.&
            &(c2.eq.NNMAPexciton(ii,6)).and.(v2.eq.NNMAPexciton(ii,7)).and.&
            &(Sx2.eq.NNMAPexciton(ii,8)).and.(Sy2.eq.NNMAPexciton(ii,9))&
            &.and.(Sz2.eq.NNMAPexciton(ii,10))) THEN

            elementId = ii  ! found element in the array
            return

            endif
        enddo


    end subroutine getNNMAPexcitonId


    subroutine appendNNMAPexciton(NNMAPexciton, NumMatEl, NumMatElused, c1,v1,Sx1,Sy1,Sz1, c2,v2,Sx2,Sy2,Sz2, elementId)
        ! Appends the nearest neighbor list if connection is not already contained in the the list.
        ! Returns element Id of the new or already contained connection
        implicit none
        INTEGER :: NNMAPexciton(NumMatEl,10), elementId
        INTEGER :: NumMatEl, c1,c2,v1,v2, Sx1,Sy1,Sz1, Sx2,Sy2,Sz2, NumMatElused

        intent(in) c1,c2,v1,v2, Sx1,Sy1,Sz1, Sx2,Sy2,Sz2, NumMatEl
        intent(inout) NNMAPexciton, NumMatElused
        intent(out) elementId

        ! check if element is already contained and what its elementId is
        call getNNMAPexcitonId(NNMAPexciton, NumMatEl, NumMatElused, c1,v1,Sx1,Sy1,Sz1, c2,v2,Sx2,Sy2,Sz2, elementId)

        if (elementId.gt.0) return ! already contained

        ! check if enough memory for the new entry is allocated
        if (NumMatEl.lt.(NumMatElused+1)) stop "Allocated NNMAPexciton is too small."

        ! create new entry
        NumMatElused = NumMatElused +1
        elementId = NumMatElused
        NNMAPexciton(NumMatElused,1) = c1
        NNMAPexciton(NumMatElused,2) = v1
        NNMAPexciton(NumMatElused,3) = Sx1
        NNMAPexciton(NumMatElused,4) = Sy1
        NNMAPexciton(NumMatElused,5) = Sz1
        NNMAPexciton(NumMatElused,6) = c2
        NNMAPexciton(NumMatElused,7) = v2
        NNMAPexciton(NumMatElused,8) = Sx2
        NNMAPexciton(NumMatElused,9) = Sy2
        NNMAPexciton(NumMatElused,10) = Sz2


    end subroutine appendNNMAPexciton

    subroutine symmetrize_TIexciton(p_id, TIexciton, NNMAPexciton, NumMatEl, NumMatEl_used)
        ! make the excitonic TIs hermitian.
        USE Geometry
        INTEGER :: NNMAPexciton(NumMatEl,10), NumMatEl, p_id, NumMatEl_used
        INTEGER :: c1,c2,v1,v2,Sx1,Sx2,Sy1,Sy2,Sz1,Sz2, elementId, ii
        !INTEGER :: c1_,c2_,v1_,v2_,Sx1_,Sx2_,Sy1_,Sy2_,Sz1_,Sz2_
        double complex :: TIexciton(NumMatEl)
        double complex val1, val2, sym_val

        intent(in) p_id, NumMatEl
        intent(inout) TIexciton, NNMAPexciton, NumMatEl_used

        if (p_id.eq.0) write(*,*) "Symmetrize excitonic TIs"

        do ii=1,NumMatEl  ! only check TIs that are in the NNMAPexciton list
            c1 = NNMAPexciton(ii,1)
            v1 = NNMAPexciton(ii,2)
            Sx1 = NNMAPexciton(ii,3)
            Sy1 = NNMAPexciton(ii,4)
            Sz1 = NNMAPexciton(ii,5)
            c2 = NNMAPexciton(ii,6)
            v2 = NNMAPexciton(ii,7)
            Sx2 = NNMAPexciton(ii,8)
            Sy2 = NNMAPexciton(ii,9)
            Sz2 = NNMAPexciton(ii,10)

            val1 = TIexciton(ii)

            call getNNMAPexcitonId(NNMAPexciton, NumMatEl, NumMatEl, c2,v2,Sx2,Sy2,Sz2, &
            &c1,v1,Sx1,Sy1,Sz1, elementId)
            if (elementId.le.0) then

            call appendNNMAPexciton(NNMAPexciton,NumMatEl,NumMatEl_used,c2,v2,Sx2,Sy2,Sz2, &
                &c1,v1,Sx1,Sy1,Sz1, elementId)
            write(*,*) "Transpose connection is not contained in NNMAPexciton!"
            write(*,*) c1,v1,Sx1,Sy1,Sz1
            write(*,*) c2,v2,Sx2,Sy2,Sz2
            TIexciton(elementId) = conjg(val1)
            else
            val2 = CONJG(TIexciton(elementId))
            sym_val = cmplx((real(val1)+real(val2))/2 , (imag(val1)+imag(val2))/2)
            TIexciton(ii) = sym_val
            TIexciton(elementId) = conjg(sym_val)
            endif
            
        enddo
    end subroutine symmetrize_TIexciton


    subroutine check_hermicity_TIexciton(TIexciton, NNMAPexciton, NumMatEl)
        ! Checks if transferintegrals that are given in the TINFILE produce
        ! a hermitian Hamiltonian.
        ! It does not check the actual Hamiltonian !!!
        USE Geometry
        INTEGER :: NNMAPexciton(NumMatEl,10), NumMatEl
        INTEGER :: c1,c2,v1,v2,Sx1,Sx2,Sy1,Sy2,Sz1,Sz2, elementId, ii
        !INTEGER :: c1_,c2_,v1_,v2_,Sx1_,Sx2_,Sy1_,Sy2_,Sz1_,Sz2_, ii, jj
        double complex :: TIexciton(NumMatEl)
        double complex val1, val2

        intent(in) TIexciton, NNMAPexciton, NumMatEl

        do ii=1,NumMatEl  ! only check TIs that are in the NNMAPexciton list
            c1 = NNMAPexciton(ii,1)
            v1 = NNMAPexciton(ii,2)
            Sx1 = NNMAPexciton(ii,3)
            Sy1 = NNMAPexciton(ii,4)
            Sz1 = NNMAPexciton(ii,5)
            c2 = NNMAPexciton(ii,6)
            v2 = NNMAPexciton(ii,7)
            Sx2 = NNMAPexciton(ii,8)
            Sy2 = NNMAPexciton(ii,9)
            Sz2 = NNMAPexciton(ii,10)

            val1 = TIexciton(ii)

            call getNNMAPexcitonId(NNMAPexciton, NumMatEl, NumMatEl, c2,v2,Sx2,Sy2,Sz2, &
            &c1,v1,Sx1,Sy1,Sz1, elementId)
            if (elementId.le.0) then
            write(*,*) "Transpose connection is not contained in NNMAPexciton!"
            write(*,*) c1,v1,Sx1,Sy1,Sz1
            write(*,*) c2,v2,Sx2,Sy2,Sz2
            stop
            endif

            val2 = CONJG(TIexciton(elementId))
            if ((abs(real(val1)-real(val2)).gt. 1e-10) .or. &
            &(abs(imag(val1)-imag(val2)).gt. 1e-10)) THEN
            write(*,*) "TIexciton is not hermitian: "
            write(*,*) c1,v1,Sx1,Sy1,Sz1
            write(*,*) c2,v2,Sx2,Sy2,Sz2
            write(*,*) val1, "!=",val2
            stop
            endif
        enddo

        write(*,*) "--> exciton TIs are hermitian :D"
    end subroutine check_hermicity_TIexciton

#ifdef DEBUG
    subroutine save_onsiteEnergies(sites_p_len,sites_p, p_id, enList, gamma)
        ! writes the on-site energies in a output file (mainly for testing purposes)
        USE Geometry
        implicit none
        INTEGER,INTENT(IN)  :: sites_p_len,sites_p(sites_p_len) ! sites in my domain neiList[#Node, #neighbor + rank_of_neighbor]
        INTEGER, INTENT(IN)   :: p_id
        DOUBLE PRECISION, INTENT(IN) :: gamma
        DOUBLE COMPLEX, INTENT(IN):: enList(sites_p_len)

        ! internal variables
        INTEGER :: cadd
        character(len=1024) :: filename
    
        write (filename, "(A7,I2,A9,I2,A4)") "test/enList",p_id,".dat"
        open(42, file=filename)
        
        ! save onsite energies
        write(*,*) "Save on-site energies at ", filename
        do cadd=1,sites_p_len
            write(42,*) sites_p(cadd), DREAL(enList(cadd)*gamma), DIMAG(enList(cadd)*gamma)
        enddo
        
        close(42)
        
    end subroutine save_onsiteEnergies
#endif
        
! #ifdef DEBUG
!       subroutine save_alltransferintegrals(sites_p_len,sites_p,neiList, p_id, traList, gamma)
!          ! writes the transfer integrals in a output file (mainly for testing purposes)
!          USE Geometry
!          implicit none
!          INTEGER,INTENT(IN)  :: sites_p_len,sites_p(sites_p_len),neiList(sites_p_len,2*nneigh)  ! sites in my domain neiList[#Node, #neighbor + rank_of_neighbor]
!          INTEGER, INTENT(IN)   :: p_id
!          DOUBLE PRECISION, INTENT(IN) :: gamma
!          DOUBLE COMPLEX, INTENT(IN):: traList(sites_p_len,nneigh)

!          ! internal variables
!          INTEGER :: cadd, neigh
!          character(len=1024) :: filename
        
!          write (filename, "(A7,I2,A9,I2,A4)") "test/traList",p_id,".dat"
!          open(43, file=filename)
!          ! save transfer integrals energies
!          write(*,*) "Save TI's energies at ", filename
!          do cadd=1,sites_p_len
!                do neigh=1, nneigh
!                   write(43,*) sites_p(cadd), neiList(cadd, neigh), DREAL(traList(cadd, neigh)*gamma) !, DIMAG(traList(cadd, neigh)*gamma)
!                enddo
!          enddo
!          close(43)
        
!       end subroutine save_alltransferintegrals
! #endif

    subroutine check_duplicates_COULOMB(p_id, filename)
        ! checks if the COULOMB or LFE file contains duplicates
        IMPLICIT NONE
        INTEGER j,  p_id,i, ALLOC_ERR
        INTEGER stat, NumMatEl
        INTEGER, allocatable, DIMENSION(:,:) :: CoulombMap
        character (len=1024) :: line, filename
        DOUBLE PRECISION WW, max_dist_coulomb, val1, val2
        LOGICAL LFI_ONS

        intent(in) filename, p_id


        ! READ file
        inquire(FILE=trim(filename),exist=LFI_ONS) ! TODO
        if(.not.LFI_ONS) then
            write(*,*) trim(filename)//' does not exist. Stop.'
            stop
        endif
        open(114,FILE=trim(filename),ACTION='READ')
        read(114,*) NumMatEl, max_dist_coulomb, WW

        if (abs(WW).lt.1e-10) then
            IF ( p_id .eq. 0 ) write(*,*) "WW = 0 --> do not check for duplicates"
            close(114)
            return
        endif


        allocate(CoulombMap(13,NumMatEl), STAT = ALLOC_ERR )
        IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION CoulombMap !!!"

        i = 1
        read_loop: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
            if (i.gt.NumMatEl) STOP "More elements in file than expected. Please check the first line of "//trim(filename)

            ! store data in arrays
            read (line, *, iostat=stat) CoulombMap(1,i),CoulombMap(2,i),CoulombMap(3,i), CoulombMap(4,i),CoulombMap(5,i),&
                & CoulombMap(6,i), CoulombMap(7,i), CoulombMap(8,i), CoulombMap(9,i), CoulombMap(10,i), &
                & CoulombMap(11,i), CoulombMap(12,i), CoulombMap(13,i), val1,val2

            if ( stat /= 0 ) then
                write(*,*) "line ", i
                STOP "Error while reading matrix element from line."
            endif

            i=i+1

        end do read_loop
        close(114)

        if ((i-1).ne.NumMatEl) STOP "Different number of elements in file "//trim(filename)//&
            &" than expected. Please check the first line of the file!"


        ! check if matrix elements are unique
        do i=1, NumMatEl
            do j=1, NumMatEl

                if (i.eq.j) cycle

                if ( (CoulombMap(1,i).eq.CoulombMap(1,j)).and.(CoulombMap(2,i).eq.CoulombMap(2,j)).and.&
                    &(CoulombMap(3,i).eq.CoulombMap(3,j)).and.(CoulombMap(4,i).eq.CoulombMap(4,j)).and.&
                    &(CoulombMap(5,i).eq.CoulombMap(5,j)).and.(CoulombMap(6,i).eq.CoulombMap(6,j)).and.&
                    &(CoulombMap(7,i).eq.CoulombMap(7,j)).and.(CoulombMap(8,i).eq.CoulombMap(8,j)).and.&
                    &(CoulombMap(9,i).eq.CoulombMap(9,j)).and.(CoulombMap(10,i).eq.CoulombMap(10,j)).and.&
                    &(CoulombMap(11,i).eq.CoulombMap(11,j)).and.(CoulombMap(12,i).eq.CoulombMap(12,j)).and.&
                    &(CoulombMap(13,i).eq.CoulombMap(13,j)) ) then

                    write(*,*) i, j
                    STOP "Matrix elements are not unique! Please check for duplicates!"
                endif

            enddo
        enddo

        deallocate(CoulombMap)
    end subroutine check_duplicates_COULOMB


end module Hamiltonian
