
! This is a software test that tests various functions of the PARFILE_mod module
!
! Author: Konrad Merkel

program test_parfile

        ! use MPI
        use PARFILE_mod
        
        implicit none
        type(PARFILE) :: param, param2
        ! integer :: ierr = 0
        ! integer :: psize,p_id

        ! Initialize MPI once
        ! call MPI_init(ierr)
        ! call MPI_Comm_rank(MPI_COMM_WORLD, p_id, ierr)
        ! call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
        ! if(ierr /= 0) stop 'ERROR: cannot initialize MPI'

        param2 = init_default_parfile()
        call param2%save("PARFILE_default", .False.)

        param = parse_parfile("PARFILE")
        call param%print_summary()
        call param%save("PARFILE_new", .True.)

        param2 = parse_parfile("PARFILE_new")
        call param2%print_summary()

        


        ! call MPI_Finalize(ierr)
        ! if (ierr /= 0) then
        !     write(*,*) ' ERROR: finalize'
        !     stop
        ! endif


    end program
