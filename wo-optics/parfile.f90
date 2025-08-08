
! Reads and writes the central input file (PARFILE)
! The PARFILE contains all parameters and flags that control the behaviour
! of the calculation. The file is written in the *.ini format.
!
! Author: Konrad Merkel


module PARFILE_mod
    use cfgio_mod
    implicit none
    !PRIVATE

    public :: PARFILE, parse_parfile, init_default_parfile


    type PARFILE
        integer :: nT, NRECURS, NPOL, N_ENERGIES, TWRITE, NNEIGH
        double precision :: EPS, DOSEMIN, DOSEMAX, T, DISORDER_STRENGTH
        ! LOGICAL :: saveWavPacks
        LOGICAL :: DEBUGGING, FORCE_HERMITIAN, CHECK_DUPLICATES
        character(len=1) :: DIRECTION
        !character(LEN=200) :: SEED_FILE
        character(LEN=20) :: DISORDER_MODEL
        character(LEN=20) :: DOS_PROJECTION

        contains
        procedure, public :: print_summary
        procedure, public :: save

    end type
    

    CONTAINS

    function parse_parfile(filename) result(param)
        character(len=*),intent(in):: filename
        type(PARFILE) :: param
        type(cfg_t):: cfg
        ! CHARACTER(len=100) :: homedir
        logical :: exists

        inquire(file=filename, exist=exists)
        if (.not.exists) then
            write(*,*) "Cannot parse PARFILE because file does not exists!"
            write(*,*) "filename: ", filename
            stop
        end if

        !read parameters from parfile
        cfg=parse_cfg(filename)
        call cfg%get("general","nT",param%nT)
        call cfg%get("general","NRECURS",param%NRECURS)
        call cfg%get("general","NPOL",param%NPOL)
        call cfg%get("general","N_ENERGIES",param%N_ENERGIES)
        call cfg%get("general","TWRITE",param%TWRITE)
        call cfg%get("general","T",param%T)
        call cfg%get("general","EPS",param%EPS)
        call cfg%get("general","DOSEMIN",param%DOSEMIN)
        call cfg%get("general","DOSEMAX",param%DOSEMAX)
        ! call cfg%get("general","saveWavPacks",param%saveWavPacks)
        call cfg%get("general","NNEIGH",param%NNEIGH)
        !call cfg%get("general","SEED_FILE",param%SEED_FILE)
        call cfg%get("general","DIRECTION",param%DIRECTION)

        if (cfg%has_key("general","DISORDER_MODEL")) then
            call cfg%get("general","DISORDER_MODEL",param%DISORDER_MODEL)
            call cfg%get("general","DISORDER_STRENGTH",param%DISORDER_STRENGTH)
        else
            param%DISORDER_MODEL="none"
            param%DISORDER_STRENGTH = 0.0
        endif

        if (cfg%has_key("general","DEBUGGING")) then
            call cfg%get("general","DEBUGGING",param%DEBUGGING)
        else
            param%DEBUGGING=.False.
        endif

        if (cfg%has_key("general","FORCE_HERMITIAN")) then
            call cfg%get("general","FORCE_HERMITIAN",param%FORCE_HERMITIAN)
        else
            param%FORCE_HERMITIAN=.FALSE.
        endif

        if (cfg%has_key("general","CHECK_DUPLICATES")) then
            call cfg%get("general","CHECK_DUPLICATES",param%CHECK_DUPLICATES)
        else
            param%CHECK_DUPLICATES=.FALSE.
        endif
        


        if (cfg%has_key("general","DOS_PROJECTION")) then
            call cfg%get("general","DOS_PROJECTION",param%DOS_PROJECTION)

            if (param%DOS_PROJECTION.ne."none" .and. param%DOS_PROJECTION.ne."x" .and. &
                &param%DOS_PROJECTION.ne."y" .and. param%DOS_PROJECTION.ne."z") THEN
                    stop "Invalid value for DOS_PROJECTION. Valid values are: none, x, y, z"
            endif
        else
            param%DOS_PROJECTION="none"
        endif


        ! param%SEED_FILE = TRIM(param%SEED_FILE)
        ! if (param%SEED_FILE(1:1).eq."~") THEN
        !     CALL get_environment_variable("HOME", homedir)
        !     param%SEED_FILE = TRIM(TRIM(homedir)//param%SEED_FILE(2:200))
        ! endif
    end function


    function init_default_parfile() result(param)
        type(PARFILE) :: param

        param%nT=1000
        param%NRECURS = 1000
        param%NPOL=20
        param%N_ENERGIES=32000
        param%TWRITE=1
        param%T= 0.2d0
        param%EPS = 0.000001d0
        param%DOSEMIN=-10.0d0
        param%DOSEMAX=10.0d0
        ! param%saveWavPacks = .False.
        param%NNEIGH = -1
        ! param%SEED_FILE = "~/Code/WannierOptics/seed/seed.0.dat"
        param%DIRECTION = "x"
        param%DISORDER_MODEL = "none"
        param%DISORDER_STRENGTH = 0.0
        param%DEBUGGING=.False.
        param%DOS_PROJECTION="none"
        param%FORCE_HERMITIAN = .FALSE.
        param%CHECK_DUPLICATES = .FALSE.
    end function


    subroutine print_summary(this)
        class(PARFILE), intent(in) :: this

        write(*,*) ""
        write(*,*) "Parameter from PARFILE:"
        write(*,*) "nT                  = ", this%nT
        write(*,*) "NRECURS             = ", this%NRECURS
        write(*,*) "NPOL                = ", this%NPOL
        write(*,*) "N_ENERGIES          = ", this%N_ENERGIES
        write(*,*) "TWRITE              = ", this%TWRITE
        write(*,*) "T                   = ", this%T
        write(*,*) "EPS                 = ", this%EPS
        write(*,*) "DOSEMIN             = ", this%DOSEMIN
        write(*,*) "DOSEMAX             = ", this%DOSEMAX
        write(*,*) "DOS_PROJECTION      = ", this%DOS_PROJECTION
        ! write(*,*) "saveWavPacks        = ", this%saveWavPacks
        write(*,*) "NNEIGH              = ", this%NNEIGH
        ! write(*,*) "SEED_FILE           = ", TRIM(this%SEED_FILE)
        write(*,*) "DIRECTION           = ", this%DIRECTION
        write(*,*) "DISORDER_MODEL      = ", this%DISORDER_MODEL
        write(*,*) "DISORDER_STRENGTH   = ", this%DISORDER_STRENGTH
        write(*,*) "DEBUGGING           = ", this%DEBUGGING
        write(*,*) "FORCE_HERMITIAN     = ", this%FORCE_HERMITIAN
        write(*,*) "CHECK_DUPLICATES    = ", this%CHECK_DUPLICATES
        write(*,*) ""

    end subroutine print_summary


    subroutine save(this, filename, overwrite)
        class(PARFILE), intent(in) :: this
        character(len=*),intent(in):: filename
        logical, intent(in) :: overwrite
        logical :: exists

        if (.not.overwrite) then
            inquire(file=filename, exist=exists)
            if (exists) then
                write(*,*) "Cannot save PARFILE because file already exists"
                write(*,*) "filename: ", filename
                return
            end if
        end if

        open(1, file = filename, status = 'REPLACE')
        write(1,*) "[general]"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) "#                              TECHNICAL PARAMETERS"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) ""
        write(1,*) "# Number of neighboring connections for every site:"
        write(1,*) "# values < 0 : use heuristics (might not work or is inefficient)"
        write(1,*) "NNEIGH       = ", this%NNEIGH
        write(1,*) ""
        ! write(1,*) "# Path to seed file:"
        ! write(1,*) "SEED_FILE    = ", TRIM(this%SEED_FILE)
        ! write(1,*) ""
        write(1,*) "# Additional tests and output? (T/F)"
        write(1,*) "DEBUGGING    = ", this%DEBUGGING
        write(1,*) ""
        write(1,*) "# Enforce Hamiltonian to be hermitian? (T/F)"
        write(1,*) "FORCE_HERMITIAN    = ", this%FORCE_HERMITIAN
        write(1,*) ""
        write(1,*) "# Check matrix elements of the Hamiltonian for duplicatesS? (T/F)"
        write(1,*) "# (Could be very slow!)"
        write(1,*) "CHECK_DUPLICATES    = ", this%CHECK_DUPLICATES
        write(1,*) ""
        write(1,*) ""
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) "#                                 TIME EVOLUTION"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) ""
        write(1,*) "# Length of time step (in internal units = largest transfer integral):"
        write(1,*) "T            = ", this%T
        write(1,*) ""
        write(1,*) "# Number of time steps:"
        write(1,*) "nT           = ", this%nT
        write(1,*) ""
        write(1,*) "# Number of Chebyshev polynomials:"
        write(1,*) "NPOL         = ", this%NPOL
        write(1,*) ""
        write(1,*) "# After which time steps you want to write output?"
        write(1,*) "TWRITE       = ", this%TWRITE
        write(1,*) ""
        write(1,*) ""
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) "#                              DENSITY OF STATES (DOS)"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) ""
        write(1,*) "# Number of Lanczos recursions"
        write(1,*) "NRECURS      = ", this%NRECURS
        write(1,*) ""
        write(1,*) "# Lorentzian broadening of the peaks"
        write(1,*) "EPS          = ", this%EPS
        write(1,*) ""
        write(1,*) "# Energy steps and boundarys for DOS output file:"
        write(1,*) "N_ENERGIES   = ", this%N_ENERGIES
        write(1,*) "DOSEMIN      = ", this%DOSEMIN
        write(1,*) "DOSEMAX      = ", this%DOSEMAX
        write(1,*) ""
        write(1,*) "# Use projected DOS? (values = none,x,y,z)"
        write(1,*) "DOS_PROJECTION = ", this%DOS_PROJECTION
        write(1,*) ""
        write(1,*) ""
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) "#                                      OPTIC"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) ""
        write(1,*) "# Polarization direction (values = x,y,z)"
        write(1,*) "DIRECTION    = ", this%DIRECTION
        write(1,*) ""
        write(1,*) ""
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) "#                                     DISORDER"
        write(1,*) "#------------------------------------------------------------------------------"
        write(1,*) ""
        write(1,*) "# Using additional disorder for diagonal parts of the Hamiltonian."
        write(1,*) "# This is only experimental! There are no derivations or extensive tests for"
        write(1,*) "# exciton Hamiltonians!!!!"
        write(1,*) "# Available disorder models = none, ANDERSON"
        write(1,*) "DISORDER_MODEL    = ", this%DISORDER_MODEL
        write(1,*) "DISORDER_STRENGTH = ", this%DISORDER_STRENGTH

        close(1) 

    end subroutine

end module PARFILE_mod
