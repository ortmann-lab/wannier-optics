!===============================================================================
!! Module: Algorithms
!!
!!   Implements numerical algorithms for Lanczos recursion, Chebyshev time 
!!   evolution, density of states (DOS) calculation, and related routines.
!!   Operates on a distributed memory parallel system using MPI and SparseCSR.
!!
!! Main author: Frank Ortmann
!! Notes: Requires MPI,SparseCSR
!===============================================================================
    module Algorithms
    use MPI
    use SparseCSR
    implicit none

    PRIVATE
    
    integer, save :: p_id               !rank of process
    integer, save :: sites_p_len        !number of sites handled by process
    double precision, save :: ac, bc    !interval parameters
    integer, save :: Npol               !number of polynomial coeficients
    integer, save :: nRecurs            !number of Lanczos recursion steps                       
    double complex, dimension(:), allocatable, save :: cn_priv                  !array containing polynomial coeficients
    double complex, dimension(:), pointer, save :: psiN_p, psiNp1_p, psiNm1_p   !local vectors, used at Lanczos recursion
    double precision, dimension(:), allocatable, save :: aLan, bLan                 !results of Lanczos recursion
                                                                                !within time evolution

    PUBLIC setupAlgorithmsModule, lanczos, iniPolynomCoeff
    PUBLIC evoIter, timeEvolution
    PUBLIC estimate_interval_spectrum, write_dos, calc_chebyshev_coeff, dequal

    CONTAINS

    
    subroutine setupAlgorithmsModule(p_id_a, nRecurs_a, ierr)
        !***************************************************************************
        !> Initializes the algorithm module by setting up required parameters and
        !! allocating memory for parallel evolution arrays.
        !!
        !! @param[in]  p_id_a     Process identifier (MPI-rank).
        !! @param[in]  nRecurs_a  Number of recursion steps for the algorithm.
        !! @param[out] ierr       Error status code (0 on success, non-zero on failure).
        !***************************************************************************
   !     use SparseCSR, only: getSites_p
        implicit none
        integer, intent(IN) :: p_id_a, nRecurs_a
        integer, intent(OUT) :: ierr

        p_id = p_id_a; nRecurs = nRecurs_a;
        sites_p_len = getSites_p()

        allocate(psiN_p(sites_p_len), STAT=ierr)
        allocate(psiNp1_p(sites_p_len), STAT=ierr)
        allocate(psiNm1_p(sites_p_len), STAT=ierr)
        allocate(aLan(nRecurs), STAT=ierr)
        allocate(bLan(nRecurs), STAT=ierr)
        if (ierr /= 0) then
            print *, p_id, 'ERROR: cannot allocate parallel evolution arrays'
            return
        end if
	    deallocate(aLan)
        deallocate(bLan)
    end subroutine setupAlgorithmsModule
    
    subroutine iniPolynomCoeff(ac_ext, bc_ext, cn_ext, Npol_ext) 
        !***************************************************************************
        !> Initializes the array containing Npol polynomial coefficients.
        !! Called by all processes.
        !!
        !! Global arguments: cn_priv, Npol
        !! Initializes: cn_priv, Npol
        !!
        !! @param[in]  ac_ext      Interval parameter a.
        !! @param[in]  bc_ext      Interval parameter b.
        !! @param[in]  cn_ext      array of polynomial coefficients.
        !! @param[in]  Npol_ext    Number of polynomial coefficients.
        !***************************************************************************
        implicit none
        integer, intent(IN) :: Npol_ext
        double precision, intent(IN) :: ac_ext, bc_ext
        double complex :: cn_ext(Npol_ext)

        ac = ac_ext
        bc = bc_ext
        Npol = Npol_ext

        if (Npol .ne. size(cn_ext)) then
            stop 'ERROR: Npol does not match size of cn_ext!'
        end if

        if (allocated(cn_priv)) deallocate(cn_priv)
        allocate(cn_priv(size(cn_ext)))
        cn_priv = cn_ext  
    end subroutine iniPolynomCoeff

    subroutine computePsiNp1(psiN_p_a, psiNm1_p_a, psiNp1_p_a, an, bnm1, InnerProd_p)
        !***************************************************************************
        !> Computes vector used in next Lanczos iteration (after normalization):
        !! |psiNp1> = H|psiN> - an|psiN> - bnm1|psiNm1>,
        !! and local value of dot product for this vector.
        !! As an input psiNp1_p represents H|psiN>.
        !! Called by all processes.
        !***************************************************************************
        implicit none
        double complex, dimension(:), pointer :: psiN_p_a, psiNm1_p_a, psiNp1_p_a
        double precision, intent(IN) :: an, bnm1
        double complex, intent(OUT) :: InnerProd_p !local result of dot product of result vector
        integer :: i

        InnerProd_p = (0.d0, 0.d0)
        do i = 1, sites_p_len
            psiNp1_p_a(i) = psiNp1_p_a(i) - an*psiN_p_a(i) - bnm1*psiNm1_p_a(i)
            InnerProd_p = InnerProd_p + CONJG(psiNp1_p_a(i))*psiNp1_p_a(i)
        end do
    end subroutine computePsiNp1


    subroutine lanczos(Psi_p, aLan, bLan)
        !***************************************************************************
        !> Performs nRecurs Lanczos iterations with Hamiltonian and state
        !! vector psiN_p, assumed as a local, permuted vector, handled by process,
        !! and fills a and b arrays with coefficients of the tridiagonal matrix.
        !! Called by all processes.
        !!
        !! Global arguments: sites_p_len, psiNp1_p, psiNm1_p
        !***************************************************************************
        implicit none
        double complex, dimension(:), pointer :: Psi_p
        double precision, dimension(nRecurs), intent(INOUT) :: aLan, bLan !arrays with coefficients of tridiagonal matrix,
                                                                    !allocated at master process
        double complex, dimension(:), pointer :: swapp
        double complex :: InnerProd_loc, InnerProd_glob
        double precision :: an, bn, normFact
        integer :: i

        do i = 1, sites_p_len
            psiN_p(i) = Psi_p(i)
            psiNm1_p(i) = (0.d0, 0.d0)
        end do
        bn = 0.d0

        do i = 1, nRecurs
	  
            call multHpsi(psiN_p, psiNp1_p, InnerProd_loc)
            call combine_scalar_product(InnerProd_loc, InnerProd_glob)

            an = DBLE(AIMAG((0.d0,1.d0)*InnerProd_glob))
            call computePsiNp1(psiN_p, psiNm1_p, psiNp1_p, an, bn, InnerProd_loc)
            call combine_scalar_product(InnerProd_loc, InnerProd_glob)
            bn = DSQRT(ABS(InnerProd_glob))
	    
            normFact = 1.d0 / bn
            psiNp1_p(1:sites_p_len) = normFact*psiNp1_p(1:sites_p_len)
            aLan(i) = an
            bLan(i) = bn

            !prepare next iteration step (use of pointers allows to avoid copying)
            swapp => psiN_p
            psiN_p => psiNm1_p
            psiNm1_p => swapp
            swapp => psiN_p
            psiN_p => psiNp1_p
            psiNp1_p => swapp
        end do
    end subroutine lanczos


    subroutine timeEvolution(T, nT, outPathStr, ierr, nsites, nTwrite)
        !***************************************************************************
        !> Performs time evolution of the system.
        !! The function computes the overlap of the evolved state with the initial state
        !! and the norm of the evolved state at each time step.
        !! The results are written to files.
        !! Called by all processes.
        !!
        !! Global arguments: sites_p_len, Npol, ac, bc, c
        !***************************************************************************
        implicit none
        integer, intent(IN) :: nsites, nTwrite
        double precision, intent(IN) :: T
        integer, intent(IN) :: nT
        character(Len=*), intent(IN) :: outPathStr
        integer, intent(OUT) :: ierr
        double complex, dimension(:), pointer :: Psi_t_p, Psi_p, Psi_0_p
        double precision, dimension(:), pointer :: aLan, bLan !arrays for storing results
                                                          !of Lanczos recursion for
                                                          !xupsiNorm, at each evolution step
        double complex, dimension(:), allocatable::Psi_t_Psi_0_p
        double precision :: ti, normPsi_t,normPsi,T_i, normPsi_0,normPsi_t_Psi_0
        double complex :: InnerProd_psi_t,InnerProd_psi,InnerProd_psi_0,InnerProd_psi_t_psi_0
        integer :: it
        
        allocate(aLan(nRecurs), STAT=ierr)
        allocate(bLan(nRecurs), STAT=ierr)
        if (ierr /= 0) then
            print *, p_id, 'ERROR: cannot allocate internal timeEvolution arrays'
            return
        end if

        if (p_id == 0) then
	        open(24,file=TRIM(outPathStr)//'T000KNormPsi')
            open(25,file=TRIM(outPathStr)//'T000KOverlap')

            ! head of the array
            write(*,*)
            write(*,*) "# -----------------------------------------------------------------------------------------------------"
            write(*,*) "# Step       | Time                   | Overlap < Psi(t) | Psi(0) >                  ",&
                &"| Norm < Psi(t) | Psi(t) >      <-- TIME"
            write(*,*) "# -----------------------------------------------------------------------------------------------------"

            ! head of array in data files
            write(24,*) "# Time                      | normPsi                | normPsi_t               ",&
                &"| normPsi_0               | normPsi_t_Psi_0"
            write(25,*) "# Time                    | Overlap < Psi(t) | Psi(0) > (real, imag)"
        end if

        !----- Lanczos recursion + CFE of initial vectors -----------------------------------------
        call noevoIter(Psi_t_p, InnerProd_psi_t,Psi_p,InnerProd_psi, Psi_0_p, InnerProd_psi_0,Psi_t_Psi_0_p,&
            InnerProd_psi_t_psi_0)

        normPsi_t = ABS(InnerProd_psi_t)
        normPsi   = ABS(InnerProd_psi)
        normPsi_0 = ABS(InnerProd_psi_0)
        normPsi_t_Psi_0 = ABS(InnerProd_psi_t_psi_0)

        if (p_id == 0) then
            ti = 0.0
            write(*,*) 0, ti, DREAL(InnerProd_psi_t_psi_0), DIMAG(InnerProd_psi_t_psi_0), normPsi_t, " <-- TIME"
            write(25,*) ti, DREAL(InnerProd_psi_t_psi_0), DIMAG(InnerProd_psi_t_psi_0)
        end if


        !-----now perform time evolution--------------------------------------------------------------      
        T_i = T
        do it = 1, nT
            ti=T_i + ti
            call evoIter(Psi_t_p, InnerProd_psi_t,Psi_p,InnerProd_psi, Psi_0_p, InnerProd_psi_0, &
                        & Psi_t_Psi_0_p,InnerProd_psi_t_psi_0)
            normPsi_t = ABS(InnerProd_psi_t)
            normPsi   = ABS(InnerProd_psi)
            normPsi_0 = ABS(InnerProd_psi_0)
            normPsi_t_Psi_0 = ABS(InnerProd_psi_t_psi_0)

            if (p_id == 0) then
                ! write relevant information in std output (as array)
                write(*,*) it, ti, DREAL(InnerProd_psi_t_psi_0), DIMAG(InnerProd_psi_t_psi_0), normPsi_t, " <-- TIME"
                write(24,*) ti, normPsi, normPsi_t, normPsi_0, normPsi_t_Psi_0
                write(25,*) ti, DREAL(InnerProd_psi_t_psi_0), DIMAG(InnerProd_psi_t_psi_0)
                call flush(24)  ! force to write in file (makes sure that values are saved even if program is canceled later)
                call flush(25)
            end if

        end do

        if (p_id == 0) then
	        close(24)
            close(25)
        end if
	 
	    deallocate(psiN_p)
        deallocate(psiNp1_p)
        deallocate(psiNm1_p)
        deallocate(aLan)
        deallocate(bLan)
	
    end subroutine timeEvolution


    subroutine evoIter(Psi_t_p, InnerProd_Psi_t,Psi_p,InnerProd_Psi,Psi_0_p, InnerProd_Psi_0,&
                & Psi_t_Psi_0_p, InnerProd_Psi_t_Psi_0)
        !***************************************************************************
        !> Routine performs one time evolution iteration of the state using the
        !! Chebyshev propagation implemented in SparseCSR. 
        !!
        !! Computes the evolved state Psi_t_p at the current time step and the 
        !! scalar products:
        !!    <Psi_t|Psi_t>, <Psi|Psi>, <Psi_0|Psi_0>, <Psi_t|Psi_0>
        !!
        !! Also returns pointers to the current evolved state, initial state, 
        !! and reference state managed in SparseCSR.
        !!        
        !! Called by all processes.
        !!
        !! Global arguments: sites_p_len, Npol, ac, bc, c
        !***************************************************************************
        implicit none
        double complex, dimension(:), pointer :: Psi_t_p,Psi_p, Psi_0_p
        double complex, intent(OUT) :: InnerProd_Psi_t,InnerProd_Psi,InnerProd_Psi_0,InnerProd_Psi_t_Psi_0
        integer :: i, tmp1
        double complex :: InnerProd_Psi_t_p, InnerProd_Psi_p, InnerProd_Psi_0_p, InnerProd_Psi_t_Psi_0_p
        double complex, dimension(sites_p_len) :: Psi_t_Psi_0_p
   	
        call time_evolution_Psi(Npol, ac, bc, cn_priv)
        call getPsi_t_p(Psi_t_p, tmp1)       ! get multiplied vector
        call getPsi_p(Psi_p, tmp1)	
        call getPsi_0_p(Psi_0_p, tmp1)         

        InnerProd_Psi_t_p = (0.d0, 0.d0)
	    InnerProd_Psi_p = (0.d0, 0.d0)
        InnerProd_Psi_0_p= (0.0d0,0.0d0)
        InnerProd_Psi_t_Psi_0_p= (0.0d0,0.0d0)

        do i = 1, sites_p_len
            InnerProd_Psi_t_p = InnerProd_Psi_t_p + Psi_t_p(i)*CONJG(Psi_t_p(i))
            InnerProd_Psi_p = InnerProd_Psi_p + Psi_p(i)*CONJG(Psi_p(i))
            InnerProd_Psi_0_p= InnerProd_Psi_0_p + Psi_0_p(i)*CONJG(Psi_0_p(i))
            InnerProd_Psi_t_Psi_0_p= InnerProd_Psi_t_Psi_0_p + Psi_0_p(i)*CONJG(Psi_t_p(i))
            
            IF(i.gt.tmp1) then
                STOP "Error in evoIter"
            ENDIF
        end do
       
        call combine_scalar_product(InnerProd_Psi_t_p, InnerProd_Psi_t)
	    call combine_scalar_product(InnerProd_Psi_p,InnerProd_Psi)
        call combine_scalar_product(InnerProd_Psi_0_p,InnerProd_Psi_0)
        call combine_scalar_product(InnerProd_Psi_t_Psi_0_p,InnerProd_Psi_t_Psi_0)	
    end subroutine evoIter


    subroutine noevoIter(Psi_t_p, InnerProd_Psi_t,Psi_p,InnerProd_Psi,Psi_0_p, InnerProd_Psi_0,&
                     & Psi_t_Psi_0_p, InnerProd_Psi_t_Psi_0)
        !***************************************************************************
        !> Performs no time evolution but computes scalar products.
        !!
        !! Called by all processes.
        !***************************************************************************
        implicit none

        double complex, dimension(:), pointer :: Psi_t_p,Psi_p, Psi_0_p
        double complex, intent(OUT) :: InnerProd_Psi_t,InnerProd_Psi,InnerProd_Psi_0,InnerProd_Psi_t_Psi_0
        integer :: i, tmp1
        double complex :: InnerProd_Psi_t_p, InnerProd_Psi_p, InnerProd_Psi_0_p, InnerProd_Psi_t_Psi_0_p
        double complex, dimension(sites_p_len) :: Psi_t_Psi_0_p

        call getPsi_p(Psi_p, tmp1)	
        call getPsi_0_p(Psi_0_p, tmp1)    
        call getPsi_t_p(Psi_t_p, tmp1)       ! get multiplied vector

        InnerProd_Psi_t_p = (0.d0, 0.d0)
	    InnerProd_Psi_p = (0.d0, 0.d0)
        InnerProd_Psi_0_p= (0.0d0,0.0d0)
        InnerProd_Psi_t_Psi_0_p= (0.0d0,0.0d0)
        
        do i = 1, sites_p_len
            InnerProd_Psi_t_p = InnerProd_Psi_t_p + Psi_t_p(i)*CONJG(Psi_t_p(i))
            InnerProd_Psi_p = InnerProd_Psi_p + Psi_p(i)*CONJG(Psi_p(i))
            InnerProd_Psi_0_p= InnerProd_Psi_0_p + Psi_0_p(i)*CONJG(Psi_0_p(i))
            InnerProd_Psi_t_Psi_0_p= InnerProd_Psi_t_Psi_0_p + Psi_0_p(i)*CONJG(Psi_t_p(i))

            IF(i.gt.tmp1) then
                STOP "Error in noevoIter"
            ENDIF
        end do
       
        call combine_scalar_product(InnerProd_Psi_t_p, InnerProd_Psi_t)
        call combine_scalar_product(InnerProd_Psi_p,InnerProd_Psi)
        call combine_scalar_product(InnerProd_Psi_0_p,InnerProd_Psi_0)
        call combine_scalar_product(InnerProd_Psi_t_Psi_0_p,InnerProd_Psi_t_Psi_0)	
    end subroutine noevoIter
    

    subroutine combine_scalar_product(InnerProd_loc, InnerProd_glob)
        !***************************************************************************
        !> Takes local results of inner products and gathers to obtain global,
        !! inner product as a sum of local values
        !!
        !! Called by all processrs
        !****************************************************************************
        implicit none
        double complex, intent(IN) :: InnerProd_loc
        double complex, intent(OUT) :: InnerProd_glob
        integer :: ierr

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"

        call MPI_ALLREDUCE(InnerProd_loc, InnerProd_glob, 1, MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
    end subroutine combine_scalar_product


    subroutine estimate_interval_spectrum(aLan, bLan, spCen, spBan, nrecurs_a)
        !***************************************************************************
        !> Estimates the lower and upper energy bounds of the spectrum of a tridiagonal matrix
        !! Energy bounds are stored in spCen and spBan.
        !! Called by all processes.
        !!
        !! @param[in]  aLan       Array of diagonal elements.
        !! @param[in]  bLan       Array of off-diagonal elements.
        !! @param[in]  nrecurs_a Number of recursion steps.
        !! @param[out] spCen   Center of the spectrum.
        !! @param[out] spBan   Half-width of the spectrum.
        !***************************************************************************
        implicit none
        integer, intent(IN) :: nrecurs_a
        double precision, intent(IN) :: aLan(nrecurs_a), bLan(nrecurs_a)
        double precision, intent(OUT) :: spCen, spBan
        double precision, dimension(:), allocatable :: am, bm
        double precision, dimension(:,:), allocatable :: z
        double precision :: emin,emax
        integer :: i, lrecurs, info
        double precision, allocatable :: work(:)

        lrecurs=Min(5000,nrecurs_a)
        allocate(am(lrecurs))
        allocate(bm(lrecurs))
        allocate(z(lrecurs, lrecurs))
        allocate(work(max(1,2*lrecurs-2)))

        ! Initialization before exact diagonalization
        z = 0.0d0
        do i = 1, lrecurs
            z(i, i) = 1.0d0
        enddo
        am(:)=aLan(:lrecurs)
        bm(:)=bLan(:lrecurs)

        ! Diagonalization of the tridiagonal matrix using LAPACK routine DSTEQR        
        call dsteqr('N', lrecurs, am, bm, z, lrecurs, work, info)
        if (info /= 0) then
            print *, 'ERROR: LAPACK routine dsteqr failed with info = ', info
            stop
        end if
        ! eigenvalues are returned in am in ascending order
        
        emin=am(1)
        emax=am(lrecurs)

        spCen=(emax+emin)/2.d0
        spBan=1.1d0*(emax-emin)/4.d0
        ! this defines an interval [spCen-2*spBan,spCen+2*spBan] containing the entire spectrum (with 10% margin)
        
        deallocate(work, am, bm, z)
    end subroutine estimate_interval_spectrum


    subroutine write_dos(aLan, bLan, nrecurs_a, fileout, norm, eps, nEN, dosEmin, dosEmax, gamma)
        !***************************************************************************
        !> Writes the density of states (DOS) to a file.
        !! The DOS is calculated using the recursion coefficients a and b.
        !! The output is normalized to the specified norm.
        !! Called by all processes.
        !!
        !! @param[in]  aLan         Array of recursion coefficients.
        !! @param[in]  bLan         Array of recursion coefficients.
        !! @param[in]  nrecurs_a Number of recursion steps.
        !! @param[in]  fileout   Output file name.
        !! @param[in]  norm      Normalization factor.
        !! @param[in]  eps       Small imaginary part added to the energy.
        !! @param[in]  nEN       Number of energy points.
        !! @param[in]  dosEmin   Minimum energy for DOS calculation.
        !! @param[in]  dosEmax   Maximum energy for DOS calculation.
        !! @param[in]  gamma     Broadening parameter.
        !***************************************************************************
        implicit none
        integer, intent(IN) :: nrecurs_a, nEN
        double precision, intent(IN) :: aLan(nrecurs_a), bLan(nrecurs_a),  norm, eps, dosEmin, dosEmax, gamma
        character(LEN=*), intent(IN) :: fileout
        double precision DOS, pi, ainf, binf, E, dE
        double complex g, z
        INTEGER i, j

        pi=4.d0*atan(1.d0)

        ! limits of recursion coefficients
        ainf=aLan(nrecurs_a)
        binf=bLan(nrecurs_a)

        dE=(dosEmax-dosEmin)/dble(nEN-1)
        open(21, FILE=TRIM(fileout))
        write(21,*) "# Density of states (DOS)"
        write(21,*) "# The ouput is normalized to ", norm
        write(21,*) "# Energy (eV),   DOS (1/eV)"
    
        do i = 1, nEN
            E=dosEmin+(i-1)*dE
            z=E*(1.d0,0.d0)+eps*(0.d0,1.d0)
            g=(z-ainf-sqrt((z-ainf)**2-4.d0*binf**2))/(2.d0*binf**2)
            if(imag(g).ge.0.d0) then
                g=(z-ainf+sqrt((z-ainf)**2-4.d0*binf**2))/(2.d0*binf**2)
            endif

            ! calculation loop for continued fraction
            do j=nrecurs_a-1,1,-1
                g=1.d0/(z-aLan(j)-bLan(j)**2*g)
            enddo
            DOS=-1.d0/pi*dimag(g)
            write(21,*) E*gamma ,DOS/gamma*norm
        enddo
        close(21)

        return
    end subroutine write_dos


    subroutine calc_chebyshev_coeff(T, ac, bc, c,Npol_a)
        !***************************************************************************
        !> Calculates Chebyshev coefficients for the time evolution of the system.
        !! The coefficients are stored in the array c.
        !! Called by all processes.
        !!
        !****************************************************************************
        implicit none
        integer,intent(in) ::  Npol_a
        integer :: i, j, info
        double precision, intent(IN) :: T, ac, bc
        double complex, intent(INOUT) :: c(Npol_a)
        double precision, dimension(:), allocatable :: d, e, work
        double precision, dimension(:,:), allocatable :: z
        double complex sc

        allocate(d(Npol_a))
        allocate(e(Npol_a))
        allocate(z(Npol_a, Npol_a))
        allocate(work(max(1,2*Npol_a-2)))

        ! Initialization before exact diagonalization
        z = 0.0d0
        do i = 1, Npol_a
            z(i, i) = 1.0d0
        enddo

        ! Initialization of the Chebyshev tridiagonal matrix:
        do i=1,Npol_a
            d(i)=ac
            e(i)=bc
        enddo
        e(1)=sqrt(2.d0)*bc

        ! diagonalization of the Hamiltonian. d(j) eigenvalue, z(i,j) vector
        call dsteqr('I', Npol_a, d, e, z, Npol_a, work, info)
        if (info /= 0) then
            print *, 'ERROR: LAPACK routine dsteqr failed with info = ', info
            stop
        end if

        ! calculation of cn(T)
        do i=1,Npol_a
            sc=(0.d0,0.d0)
            do j=1,Npol_a
                sc=sc+z(1,j)*z(i,j)*exp(-(0.d0,1.d0)*d(j)*T)
            enddo
            c(i)=sc
        enddo

        deallocate(d, e, z, work)
        return
    end subroutine calc_chebyshev_coeff


    function dequal(a, b)
        !***************************************************************************
        !> Compares two double precision numbers for equality.
        !****************************************************************************
        implicit none
        DOUBLE PRECISION, intent(in)  :: a,b
        logical :: dequal

        dequal = dabs(a-b).lt.1e-16
    end function dequal

end module Algorithms
