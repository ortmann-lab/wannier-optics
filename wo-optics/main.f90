
! WannierOptics.x
!
! This programm solves the exction Hamiltonian in a basis of maximally
! localized Wannier functions using a time a time domain approach.
!
! Authors: Konrad Merkel, Frank Ortmann


program WannierOptics
    use MPI
    use SparseCSR
    use Algorithms
    use Geometry
    use test_geometry_mod
    use Hamiltonian
    use PARFILE_mod
    
    implicit none
    integer             :: psize,p_id
    type(PARFILE)       :: config
    character(LEN=200)  :: outPathStr, densPathStr, arg
    logical             :: file_exists

    ! for timing    
    double precision    :: slaveTime, maxTime, timeSum
    double precision    :: progTiming(9), opTiming(2,3)

    ! Decomposition of Sample
    integer, dimension(:), allocatable, target :: sites_p
    integer :: sites_p_len

    ! basic arrays
    integer, dimension(:,:), allocatable, target            :: neiList      ! neighbours list
    double complex, dimension (:,:), allocatable, target    :: traList      ! neighbours energies
    double complex, dimension(:), allocatable, target       :: enList       ! vector of self energies
    double complex, dimension(:), pointer                   :: PsiIn_p, PsiIn_t_p, PsiIn_0_p ! input wavepackets
    double precision, dimension(:), allocatable             :: aLan, bLan   ! results of Lanczos recursion
    double precision                        :: ac, bc   ! spectral parameters
    double complex, dimension(:), pointer   :: cn       ! list of Chebyshev polynomial coefficients
    double complex, dimension(:), pointer   :: Psi_p    ! local part of localized state vector used in
                                                        ! time evolution loop and allocated at SparseCSR
    integer             :: ierr = 0 ! error handling
    integer             :: i, tmp1, NPOL_converged
    integer             :: globId, Sx,Sy,Sz,m,l, n1,n2,n3  ! geometry and indexes
    integer*4           :: nsites
    double complex      :: dipole_norm, dipole_norm_proj
    double precision    :: pi,gamma, VOLUME ! preparing output parameters
    REAL                :: randf
    
    pi=4.d0*atan(1.d0)
    !------------------------ End variable declaration --------------------------

    ! Initialize MPI once
    call MPI_init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, p_id, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
    if(ierr .ne. 0) stop 'ERROR: cannot initialize MPI'

    if (p_id .eq. 0) then

        ! banner
        write(*,*) ""
        write(*,*) " __      __                                               "
        write(*,*) "/\ \  __/\ \                           __                 "
        write(*,*) "\ \ \/\ \ \ \     __      ___     ___ /\_\     __   _ __  "
        write(*,*) " \ \ \ \ \ \ \  /'__`\  /' _ `\ /' _ `\/\ \  /'__`\/\`'__\"
        write(*,*) "  \ \ \_/ \_\ \/\ \L\.\_/\ \/\ \/\ \/\ \ \ \/\  __/\ \ \/ "
        write(*,*) "   \ `\___x___/\ \__/.\_\ \_\ \_\ \_\ \_\ \_\ \____\\ \_\ "
        write(*,*) "    '\/__//__/  \/__/\/_/\/_/\/_/\/_/\/_/\/_/\/____/ \/_/ "
        write(*,*) "                                                          "
        write(*,*) "                                                          "
        write(*,*) "         _____           __                               "
        write(*,*) "        /\  __`\        /\ \__  __                        "
        write(*,*) "        \ \ \/\ \  _____\ \ ,_\/\_\    ___    ____        "
        write(*,*) "         \ \ \ \ \/\ '__`\ \ \/\/\ \  /'___\ /',__\       "
        write(*,*) "          \ \ \_\ \ \ \L\ \ \ \_\ \ \/\ \__//\__, `\      "
        write(*,*) "           \ \_____\ \ ,__/\ \__\\ \_\ \____\/\____/      "
        write(*,*) "            \/_____/\ \ \/  \/__/ \/_/\/____/\/___/       "
        write(*,*) "                     \ \_\                                "
        write(*,*) "                      \/_/                                "
        write(*,*) ""
        write(*,*) "Citation:  Merkel, K. & Ortmann, F. Journal of Physics: Materials 7, 015001 (2024)";
        write(*,*) "           https://dx.doi.org/10.1088/2515-7639/ad06cd";
        write(*,*) ""
        write(*,*) ""
        write(*,*) ""

        write(*,*) "Calculation is performed with ", psize, "processes" 
        progTiming(1) = MPI_Wtime()
    end if


    ! get command line arguments and help
    if (command_argument_count().gt.0) then
        if (p_id.eq.0) then
            write(*,*) ""
            call get_command_argument(1, arg)
            if (arg.eq."-g") then
                write(*,*) "Generate PARFILE"
                config = init_default_parfile()
                call config%save("PARFILE", .False.)
                write(*,*) "Done."
            else
                write(*,*) "This program calculates optical proerties of semiconductors bases on"
                write(*,*) "maximally localized Wannier funcitons."
                write(*,*) "To generate the main configuration file (PARFILE) please run"
                write(*,*) ""
                write(*,*) "  ./WannierOptic.x -g"
                write(*,*) ""
                write(*,*) "You can than edit the PARFILE and adapt to your desires..."
                write(*,*) ""
                write(*,*) "To run actual calculations use:"
                write(*,*) ""
                write(*,*) "  ./WannierOptic.x"
                write(*,*) ""
                write(*,*) ""
                write(*,*) "Have fun :D"
            endif
        endif
        call MPI_Finalize(ierr)
        stop
    end if


    !----------------------------------------------------------------------
    ! Read PARFILE and prepare output directory                         ---
    !----------------------------------------------------------------------

    ! read parfile
    config = parse_parfile("PARFILE")
    if (p_id.eq.0)  call config%print_summary()

    call set_geometry_debugging(config%DEBUGGING)


    ! create output directory
    if (p_id.eq.0) CALL system('mkdir -p ./output')
    outPathStr = './output/'

    ! save copy of PARFILE to make it easier to use
    ! new PARFILE in the transition periode
    call config%save("newPARFILE", .TRUE.)

    ! set seed for each process
    CALL init_random_seed(p_id)


    !-------------------------------------------------------------------------------
    ! Construct supercell / Domain decomposition / Allocation of sparse matrices ---
    !-------------------------------------------------------------------------------

    ! read geometry from POSFILE
    call parse_posfile(p_id, "POSFILE")
    if (p_id .eq. 0) call print_geometry()

    ! domain decomposition
    call getNumSiteInLocalDomain(p_id, sites_p_len)
    allocate(sites_p(sites_p_len), STAT=ierr)

    ! fill sites_p array = mapping of local and gobal indexes for local domain
    call getMappingLocalToGlobal(p_id, sites_p,sites_p_len)


    ! synchronize all processes
    call flush(6)  ! flush standard output for better appearance
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"

    ! optionally perform several tests of the coordinate transforms and process mapping
    if (config%DEBUGGING) then
        call test_indexes(p_id, sites_p, sites_p_len)
        call testProcessMapping(p_id, sites_p, sites_p_len,psize)
        call test_getLocalId(p_id, sites_p, sites_p_len)

        ! synchronize all processes
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr .ne. 0) STOP "Error in MPI_BARRIER"
    endif

    ! set number of neighbors per site
    if (config%NNEIGH .le. 0) THEN
        ! determine number of transfer integrals
        call estimate_num_neighbors(p_id, tmp1)
        config%NNEIGH = tmp1
        if (p_id.eq.0) write(*,*) "Estimate number of neighbors per site to NNEIGH=", config%NNEIGH
    endif

    ! allocate memory (for each process)
    allocate(neiList(sites_p_len,config%NNEIGH * 2), STAT=ierr)
    allocate(traList(sites_p_len,config%NNEIGH), STAT=ierr)
    allocate(enList(sites_p_len), STAT=ierr)
    if(ierr .ne. 0) stop ' ERROR: failure allocating arrays for storing the Hamiltonian'

    allocate(PsiIn_p(sites_p_len), STAT=ierr)
    allocate(PsiIn_t_p(sites_p_len), STAT=ierr)
    allocate(PsiIn_0_p(sites_p_len), STAT=ierr)
    if(ierr .ne. 0) stop ' ERROR: failure allocating arrays for storing the wave functions'



    !----------------------------------------------------------------------
    ! Construction of random-phase state for DOS                        ---
    !----------------------------------------------------------------------
    nsites = get_total_num_sites()
    if (config%DOS_PROJECTION .eq. "none") THEN

        if(p_id.eq.0) then
            write(*,*) 
            write(*,*) 'Generate random phase state (for all sites)...'
        endif

        ! create random phase state
        do i=1,sites_p_len
            call random_number(randf)
            PsiIn_p(i)=exp(2.d0*pi*(0.d0,1.d0)*randf)/dsqrt(1.d0*nsites)
        enddo

    else
        if(p_id.eq.0) then
            write(*,*) 
            write(*,*) 'Generate projected random phase state ...'
            write(*,*) 'Projection onto ', config%DOS_PROJECTION
        endif
        call get_supercell_dimensions(n1,n2,n3)
        do i=1,sites_p_len
    
            globId = sites_p(i)
            call getExcitonIndex_shifted(globId, Sx,Sy,Sz,m,l)
    
            ! we generate random numbers at every step. Even if we dont need them! This makes
            ! sure that we can compare with the unprojected DOS without worries about convergence!
            call random_number(randf)

            if (config%DOS_PROJECTION.eq.'x') then
                if (Sx.eq.0) then
                    PsiIn_p(i)=exp(2.d0*pi*(0.d0,1.d0)*randf) /dsqrt(1.d0*nsites) * dsqrt(1.d0*n1)  ! norm is important for lanczos!
                else
                    PsiIn_p(i)=cmplx(0.0d0,0.0d0)
                endif
            else if (config%DOS_PROJECTION.eq.'y') then
                if (Sy.eq.0) then
                    PsiIn_p(i)=exp(2.d0*pi*(0.d0,1.d0)*randf) /dsqrt(1.d0*nsites) * dsqrt(1.d0*n2)  ! norm is important for lanczos!
                else
                    PsiIn_p(i)=cmplx(0.0d0,0.0d0)
                endif
            else if (config%DOS_PROJECTION.eq.'z') then
                if (Sz.eq.0) then
                    PsiIn_p(i)=exp(2.d0*pi*(0.d0,1.d0)*randf) /dsqrt(1.d0*nsites) * dsqrt(1.d0*n3)  ! norm is important for lanczos!
                else
                    PsiIn_p(i)=cmplx(0.0d0,0.0d0)
                endif
            endif
        enddo
    
    endif
    


    ! for debugging purposes
    ! call read_randomphasestate(sites_p_len, p_id, PsiIn_p)

    ! output first elements of random phase state for every mpi process
    if (config%DEBUGGING) then
        write(*,*) "CPU",p_id+1, ": PsiIn_p=", (PsiIn_p(i), i=1,2), '...'
        !write(*,*) "CPU",p_id+1, ": Norm PsiIn_p=", (abs(PsiIn_p(i)*dsqrt(1.d0*nsites)), i=1,5), '...'
    endif


    !-----------------------------------------------------------------------------
    ! Initialize optical transition dipoles as inital states (xupsi and yupsi) ---
    !-----------------------------------------------------------------------------
    
    if(p_id.eq.0) then
        write(*,*) ""
        write(*,*) "Set up optical transition dipoles..."
    endif

    call setOpticalDipole(p_id, "TRANSITION", PsiIn_0_p, sites_p_len, sites_p, config, dipole_norm)

    inquire(FILE="TRANSITION_proj",exist=file_exists)
    if(file_exists) then
        if (p_id.eq.0) write(*,*) "Set projected dipole moments! Please note that you are calculating a projected spectrum!"
        call setOpticalDipole(p_id, "TRANSITION_proj", PsiIn_t_p, sites_p_len, sites_p, config, dipole_norm_proj)
    else
        if (p_id.eq.0) write(*,*) "Calculate actual autocorrelation function without any projection."
        do i=1,sites_p_len
            PsiIn_t_p(i) = PsiIn_0_p(i)
        enddo
    endif
    if (config%DEBUGGING) then
        write(*,*) "CPU",p_id+1, ": PsiIn_t_p=", (PsiIn_t_p(i), i=1,2), '...'
        !write(*,*) "CPU",p_id+1, ": Norm PsiIn_p=", (abs(PsiIn_p(i)*dsqrt(1.d0*nsites)), i=1,5), '...'
    endif

    progTiming(2) =  MPI_Wtime()

    ! create parallelized Hamiltonian
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (p_id .eq. 0) then
        write(*,*) ""
        write(*,*) "Setup and state preparation finished."
        write(*,*) "================================================================"
        write(*,*) ""
        write(*,*) "Create Hamiltonian ..."
    endif



    !----------------------------------------------------------------------
    ! Construction of Hamiltonian                                       ---
    !----------------------------------------------------------------------

    call create_exciton_hamiltonian(p_id, sites_p_len, sites_p,neiList,enList,traList, gamma, config)

    ! for debugging:
    ! call save_onsiteEnergies(sites_p_len, sites_p, p_id, enList, gamma)
    ! call save_alltransferintegrals(sites_p_len,sites_p,neiList, p_id, traList, gamma)
    ! call save_randomphasestate(sites_p_len,sites_p, p_id, PsiIn_p)

    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! if (ierr .ne. 0) STOP "Error in MPI_BARRIER"

    ! CALL readOpticalDipole(p_id) ! only for test purposes!

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"
    progTiming(3) =  MPI_Wtime()


    !set Hamiltonian and states in SparseCSR module
    call setHamiltPsi(get_total_num_sites(), config%NNEIGH, sites_p_len, neiList, traList, &
    &enList, PsiIn_p, PsiIn_t_p, PsiIn_0_p)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"


    !set local variables at SparseCSR
    if (p_id .eq. 0) write(*,*) 'Setting ParallelMatVec module...'
    call setupSparseCSRModule(sites_p,p_id, MPI_COMM_WORLD, get_total_num_sites(), config%NNEIGH) 




    !----------------------------------------------------------------------
    ! Create OUTFILE                                                    ---
    !----------------------------------------------------------------------
    if (p_id .eq. 0) then

        VOLUME = get_volume_unitcell()
        call get_supercell_dimensions(n1,n2,n3)

        OPEN(111,FILE=TRIM(outPathStr)//'OUTFILE',Status='REPLACE')
        write(111,*) "[output]"
        write(111,*) "Num CPU       = ", psize
        write(111,*) "nT            = ", config%nT
        write(111,*) "NRECURS       = ", config%NRECURS
        write(111,*) "NPOL          = ", config%NPOL
        write(111,*) "N_ENERGIES    = ", config%N_ENERGIES
        write(111,*) "TWRITE        = ", config%TWRITE
        write(111,*) "T             = ", config%T
        write(111,*) "EPS           = ", config%EPS
        write(111,*) "DOSEMIN       = ", config%DOSEMIN
        write(111,*) "DOSEMAX       = ", config%DOSEMAX
        write(111,*) "saveWavPacks  = ", .False. !config%saveWavPacks
        write(111,*) "NNEIGH        = ", config%NNEIGH
        ! write(111,*) "SEED_FILE     = ", TRIM(config%SEED_FILE)
        write(111,*) "DIRECTION     = ", config%DIRECTION
        write(111,*) "cell dim1     = ", n1
        write(111,*) "cell dim2     = ", n2
        write(111,*) "cell dim3     = ", n3
        write(111,*) "gamma         = ", gamma
        write(111,*) "Volume        = ", VOLUME
        write(111,*) "Volume / site = ", VOLUME/(get_N_con()*get_N_val())
        write(111,*) "dipole_norm   = ", DREAL(dipole_norm)
        close(111)

        write(*,*) ""
        write(*,*) 'Performing preprocessing...'
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"




    !----------------------------------------------------------------------
    ! Preprocessing (encapsulate and internally reorder Hamiltonian)    ---
    !----------------------------------------------------------------------

    !prepocessing step of SparseCSR
    call encapsHamiltonian()

    ! free memory
    deallocate(neiList)
    deallocate(traList)
    deallocate(enList)
    deallocate(PsiIn_p)
    deallocate(PsiIn_t_p)
    deallocate(PsiIn_0_p)

    call ImplemCSR_IB(ierr) 
    if (ierr .ne. 0) stop 'ERROR in preprocess'

    if (p_id .eq. 0) then
        write(*,*) ""
        write(*,*) "Preprocessing finished."
        write(*,*) "================================================================"
        write(*,*) ""
    endif




    !----------------------------------------------------------------------
    ! Lanczos recursion to calculate DOS                                ---
    !----------------------------------------------------------------------

    allocate(aLan(config%NRECURS))
    allocate(bLan(config%NRECURS))

    !initializes internal data of Algorithms_mod, needed
    !at Lanczos recursion and each evolution step
    call setupAlgorithmsModule(p_id, config%NRECURS, ierr) 
    if (ierr .ne. 0) stop 'ERROR in setupAlgorithmsModule'

    call getPsi_p(Psi_p, tmp1)    !get local initial state vector
                                        !managed in SparseCSR
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"

    if (p_id .eq. 0) write(*,*) 'Starting DOS calculation...'
    if (p_id .eq. 0) progTiming(4) = MPI_Wtime()
    slaveTime = MPI_Wtime()

    call lanczos(Psi_p, aLan, bLan)  !Lanczos recursion with Hamiltonian
    slaveTime = MPI_Wtime() - slaveTime
    call parallelTiming(slaveTime, maxTime, timeSum)
    opTiming(1, 1) = maxTime
    opTiming(2, 1) = timeSum            !and init state vector

    if (p_id .eq. 0) then
        write(*,*) "DOS-Lanczos finished :D"
        write(*,*) ""
        open(UNIT=38,FILE=TRIM(outPathStr)//'T000Kanbn')
        do i=1,config%NRECURS
            write(38,*) aLan(i),bLan(i)
        enddo
        close(38)
    end if
    
    !----------------------------------------------------------------------
    ! One determines now the interval [ac-2bc,ac+2bc] containing        ---
    ! the whole spectrum (add 10% for security)                         ---
    !----------------------------------------------------------------------
    
    if (p_id .eq. 0) progTiming(5) = MPI_Wtime()
    call estimate_interval_spectrum(aLan, bLan, ac, bc, config%NRECURS)
    if (p_id.eq.0) then
        !write(*,*) 'R',p_id,'ac bc',ac,bc
       densPathStr = TRIM(outPathStr)//'T000Kdensity.dat'
        ! write density of states
        call write_dos(aLan, bLan, config%NRECURS, densPathStr, 1.d0, config%EPS, config%N_ENERGIES, &
            &config%DOSEMIN/gamma, config%DOSEMAX/gamma, gamma)
        progTiming(6) = MPI_Wtime()
    end if



    !----------------------------------------------------------------------
    ! Time evolution                                                    ---
    !----------------------------------------------------------------------

    ! prepare time evolution
    allocate(cn(config%NPOL))

    call calc_chebyshev_coeff(config%T, ac, bc, cn,config%NPOL)
    if (p_id .eq. 0) then
        open(UNIT=654,FILE=TRIM(outPathStr)//'T000Kcn')
        do i=1,config%NPOL
            write(654,*) cn(i)
        enddo
        close(654)
    endif

    ! check if number of polynomials is converged
    if (p_id.eq.0) then
        NPOL_converged = -1
        do i=1,config%NPOL
            if ((ABS(REAL(cn(i))) .lt. 1e-15) .and. (ABS(IMAG(cn(i))) .lt. 1e-15)) THEN
                NPOL_converged = i
                exit
            endif
        enddo
        if (NPOL_converged.gt.0) THEN
            write(*,*) "Number of Chebyshev polynomials is converged for ", NPOL_converged, "polynomials."
            if (NPOL_converged+1.lt.config%NPOL) THEN

                write(*,*) ""
                write(*,*) ""
                write(*,*) "---------------------------- NOTICE --------------------------"
                write(*,*) ""
                write(*,*) "You are not using the optimal number of Chebyshev polynomials!"
                write(*,*) ""
                write(*,*) " Please set     NPOL = ", NPOL_converged
                write(*,*) ""
                write(*,*) "-------------------------------------------------------------"
                write(*,*) ""
                write(*,*) ""

            endif
        else
            write(*,*) ""
            write(*,*) ""
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(*,*) ""
            write(*,*) " You are not converged with the number of Chebyshev polynomials!"
            write(*,*) " Please increase the value of NPOL in the PARFILE."
            write(*,*) ""
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(*,*) ""
            write(*,*) ""
        endif
    endif

    ! sets the spectral bounds and array cn in module Algorithms
    call iniPolynomCoeff(ac, bc, cn, config%NPOL)
    deallocate(cn)
    
    call allocateArrays(ierr)    !allocate first set of arrays used during time evolution at SparseCSR module
    if (ierr .ne. 0) stop 'ERROR in allocateArrays'


    if (p_id .eq. 0) write(*,*) 'Entering time evolution...'
    slaveTime = MPI_Wtime()
    call timeEvolution(config%T, config%nT, outPathStr, ierr, get_total_num_sites(), config%TWRITE)    
    if (ierr .ne. 0) stop 'ERROR in timeEvolution'
    
    slaveTime = MPI_Wtime() - slaveTime
    call parallelTiming(slaveTime, maxTime, timeSum)
    opTiming(1, 2) = maxTime
    opTiming(2, 2) = timeSum

    if (p_id .eq. 0) then
        progTiming(7) = MPI_Wtime()
        write(*,*) "Time evolution finished :D"
        progTiming(8) = MPI_Wtime()
    end if


    deallocate(aLan)
    deallocate(bLan)

    !------deallocates dynamic fields in MODULE parallelMatVec   --------------------
    call deallocateArrays()
    ! deallocate(sites_p)
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0) STOP "Error in MPI_BARRIER"



    !----------------------------------------------------------------------
    !   Calculations finished --> create summary and outfile            ---
    !----------------------------------------------------------------------

    progTiming(8) = MPI_Wtime()
    call MPI_Finalize(ierr)
    if (ierr .ne. 0) stop 'ERROR during MPI_Finalize'

    if (p_id .eq. 0) then
        print *, ""
        print *, "All calculations finished."
        print *, "================================================================"
        print *, ""
        print *, "Timing info (in seconds)"
        print *, "whole execution               : ", progTiming(8)-progTiming(1)
        print *, "creating Hamiltonian          : ", progTiming(2)-progTiming(1)
        print *, "PSI/XUPSI reading or creating : ", progTiming(3)-progTiming(2)
        print *, "all preprocessing up to DOS   : ", progTiming(4)-progTiming(3)
        print *, "Lanczos DOS                   : ", progTiming(5)-progTiming(4)
        print *, "Interval + Density            : ", progTiming(6)-progTiming(5)
        print *, "Time evolution                : ", progTiming(7)-progTiming(6)
        print *, ""
        print *, "one Lanczos recursion max     : ", opTiming(1,1), ' , ', 'sum :', opTiming(2,1)
        print *, "timeEvolution ", opTiming(1,2), ' , ', 'sum :', opTiming(2,2)
    endif
	
end program WannierOptics



!----------------------------------------------------------------------
SUBROUTINE init_random_seed(p_id)
    !
    !   Reads seed file and sets the seed for every process
    !   such that every process has a unique series of pseudorandom
    !   numbers, which are not correlated and non-overlaping.
    !
    use seed_data
    IMPLICIT NONE
    INTEGER, allocatable :: seed(:)
    DOUBLE PRECISION     :: random1, random2
    INTEGER :: i,j, n, p_id
    intent(in) p_id

    !how large is the necessary size of the seed  
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    if (p_id.eq.0) then
        write(*,*) 'Initialize random numbers generator with seed0 (hardcoded values)'
        write(*,*) 'Seed size: ', n
        write(*,*) ''
    endif

    if (n.gt.SIZE(seed0)) then
        STOP "Cannot initialize random numbers generator because we don't have enough seed0 values!"
    endif

    ! set same seed for all processes
    do i=1,n
        seed(i)=seed0(i)
    enddo
    CALL RANDOM_SEED(PUT = seed)


    ! make sure every process has a different sequence of random numbers
    do j=1,2*n*(p_id+1)      ! every process ends up with a different number
        call random_number(random1)
    enddo
    
    ! fill seed array
    do j=1,n
        call random_number(random1)
        call random_number(random2)  ! prefactor
        if (random2.ge.0.5) then
            seed(j)  = int(random1*1000000000)
        else
            seed(j)  = -int(random1*1000000000)
        endif
    end do
    call random_seed(put=seed)  ! set seed
    call random_seed(get=seed)

    ! if (p_id.le.10) write(*,*) 'Process: ', p_id, 'Seed: ', seed

END SUBROUTINE init_random_seed



subroutine read_randomphasestate(sites_p_len, p_id, PsiIn_p)
    ! reads the random phase state (mainly for testing purposes)
    USE Geometry
    implicit none
    INTEGER,INTENT(IN)  :: sites_p_len ! sites in my domain vv[#Node, #neighbor + rank_of_neighbor]
    INTEGER, INTENT(IN)   :: p_id
    DOUBLE COMPLEX, INTENT(OUT):: PsiIn_p(sites_p_len)
    
    ! internal variables
    INTEGER :: cadd, tmp, locId
    double precision :: val1, val2
    character(len=1024) :: filename
    
        write (filename, "(A7,I2,A9,I2,A4)") "test/psi",p_id,".dat"
        open(42, file=filename)
        
        ! save onsite energies
        write(*,*) "Save random phase state at ", filename
        do cadd=1,sites_p_len
            read(42,*) tmp, locId, val1, val2
            PsiIn_p(locId) = CMPLX(val1, val2)
        enddo
    
        close(42)
    
end subroutine read_randomphasestate



subroutine setOpticalDipole(p_id, filename, state, sites_p_len, sites_p, config, norm)  
    ! Reads transition matrix elements (optical dipoles) from file and sets initial state
    use Geometry
    use PARFILE_mod
    IMPLICIT NONE
    INTEGER i, j, ALLOC_ERR, stat, NumMatEl, p_id, sites_p_len, direction_index
    double complex, intent(out) :: state(sites_p_len)
    integer :: sites_p(sites_p_len), globId, Sx,Sy,Sz,m,l
    INTEGER, allocatable, DIMENSION(:,:) :: opticalDipoleMap
    double precision, allocatable, dimension(:,:) :: opticalDipole
    type(PARFILE), intent(in) :: config
    character(len=*),intent(in):: filename
    double complex, intent(out) :: norm
    character (len=500) :: line
    LOGICAL LFI_ONS


    if (config%DIRECTION.eq."x") then 
        direction_index = 1
    else if (config%DIRECTION.eq."y") then
        direction_index = 2
    else if (config%DIRECTION.eq."z") then
        direction_index = 3
    else 
        stop "Wrong cartesian direction for optical dipole moment"
    end if


    if (p_id.eq.0) then
        write(*,*) "Read file: ", filename
    endif

    inquire(FILE=filename,exist=LFI_ONS)
    if(.not.LFI_ONS) stop 'TRANSITION file does not exist.'

    open(114,FILE=filename,ACTION='READ')
    read(114,*) NumMatEl


    allocate(opticalDipoleMap(5,NumMatEl), STAT = ALLOC_ERR )
    IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION opticalDipoleMap !!!"

    allocate(opticalDipole(3, NumMatEl), STAT = ALLOC_ERR )
    IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION opticalDipole !!!"

    i = 1
    read_loop: do
        read (114, '(A)', iostat=stat)  line

        if ( stat .ne. 0 )  exit read_loop             ! check end of file or other errors
        if ( line .eq. ' ') cycle read_loop          ! skip empty lines
        if (index (line, "#").ne. 0)  cycle read_loop  ! skip comment lines
        if (i.gt.NumMatEl) STOP "More elements in TRANSITION file than expected. Please check the first line of the file!"

        ! store data in arrays
        read (line, *, iostat=stat) opticalDipoleMap(1,i),opticalDipoleMap(2,i),&
            & opticalDipoleMap(3,i), opticalDipoleMap(4,i),opticalDipoleMap(5,i),&
            & opticalDipole(1,i), opticalDipole(2,i), opticalDipole(3,i)

        if ( stat .ne. 0 ) then
            write(*,*) "line ", i
            STOP "Error while reading transition matrix element from line."
        endif

        i=i+1
    end do read_loop
    close(114)

    if ((i-1).ne.NumMatEl) STOP "Different number of elements in TRANSITION file than expected. &
        &Please check the first line of the file!"

    !if (p_id.eq.0) write(*,*) "Found ", i-1, "matrix elements. Done reading TRANSITION file."

    if (config%DEBUGGING) then
        if (p_id.eq.0) then
            write(*,*) "TRANSITION matrix elements (from file):"
            do i=1, min(NumMatEl, 10)
                write(*,*) opticalDipoleMap(1,i),opticalDipoleMap(2,i),opticalDipoleMap(3,i),&
                    & opticalDipoleMap(4,i),opticalDipoleMap(5,i), opticalDipole(1,i),&
                    & opticalDipole(2,i), opticalDipole(3,i)
            enddo
            if (NumMatEl.gt.10) write(*,*) " ... "
            write(*,*) "#elements = ", NumMatEl
        endif
    endif

    ! check if matrix elements are unique
    do i=1, NumMatEl
        do j=1, NumMatEl

            if (i.eq.j) cycle

            if ( (opticalDipoleMap(1,i).eq.opticalDipoleMap(1,j)).and.&
                &(opticalDipoleMap(2,i).eq.opticalDipoleMap(2,j)).and.&
                &(opticalDipoleMap(3,i).eq.opticalDipoleMap(3,j)).and.&
                &(opticalDipoleMap(4,i).eq.opticalDipoleMap(4,j)).and.&
                &(opticalDipoleMap(5,i).eq.opticalDipoleMap(5,j)) ) then

                write(*,*) i, j
                STOP "TRANSITION matrix elements are not unique!"
            endif
        enddo
    enddo

    ! normalize optical dipole state
    norm = cmplx(0.0,0.0)
    do j=1, NumMatEl
        norm = norm + cmplx(opticalDipole(direction_index,j), 0.0)**2.0d0  ! we cast to complex to get better accuracy
    enddo
    if (p_id.eq.0) then
        write(*,*) "Normalize optical dipole state... (initial norm = ", REAL(norm),")"
    endif


    ! create inital state from optical dipoles
    do i=1, sites_p_len
        globId = sites_p(i)
        call getExcitonIndex_shifted(globId, Sx,Sy,Sz,m,l)

        state(i) = cmplx(0.0,0.0)

        do j=1, NumMatEl
            if ((opticalDipoleMap(1,j).eq.l).and.(opticalDipoleMap(2,j).eq.m).and.&
                &(opticalDipoleMap(3,j).eq.Sx).and.(opticalDipoleMap(4,j).eq.Sy).and.&
                &(opticalDipoleMap(5,j).eq.Sz)) THEN


                state(i) = cmplx(opticalDipole(direction_index,j), 0.0) / dsqrt(REAL(norm))
                exit
            endif

        enddo
    enddo

    deallocate(opticalDipoleMap)
    deallocate(opticalDipole)

end subroutine setOpticalDipole


subroutine parallelTiming(slaveTime, maxTime, timeSum)
    !Takes execution time for part of the code at each process and computes
    !maximum and sum of these values, stored by master process
    use MPI
    implicit none

    double precision, intent(IN) :: slaveTime
    double precision, intent(OUT) :: maxTime, timeSum
    double precision :: tmp1d
    integer :: ierr 

    call MPI_REDUCE(slaveTime, tmp1d, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    maxTime = tmp1d

    call MPI_REDUCE(slaveTime, tmp1d, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    timeSum = tmp1d

end subroutine parallelTiming
