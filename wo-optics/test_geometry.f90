
! This is a software test that tests various functions of the Geometry module
!
! Author: Konrad Merkel

program test_geometry

        use MPI
        use Geometry
        use test_geometry_mod
        
        implicit none
        integer :: p_id, ierr, psize
        integer, dimension(:), allocatable, target :: sites_p
        integer :: sites_p_len

        ! Initialize MPI once
        call MPI_init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, p_id, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
        if(ierr /= 0) stop 'ERROR: cannot initialize MPI'


        call parse_posfile(p_id, "POSFILE")
        call print_geometry()


        ! domain decomposition
        call getNumSiteInLocalDomain(p_id, sites_p_len)
        allocate(sites_p(sites_p_len), STAT=ierr)

        ! fill sites_p array = mapping of local and gobal indexes for local domain
        call getMappingLocalToGlobal(p_id, sites_p,sites_p_len)


        ! synchronize all processes
        call flush(6)  ! flush standard output for better appearance
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"

        ! perform several tests of the coordinate transforms and process mapping
        call test_indexes(p_id, sites_p, sites_p_len)
        call testProcessMapping(p_id, sites_p, sites_p_len,psize)
        call test_getLocalId(p_id, sites_p, sites_p_len)

        ! synchronize all processes
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (ierr /= 0) STOP "Error in MPI_BARRIER"




        call MPI_Finalize(ierr)
        if (ierr /= 0) stop ' ERROR: finalize'


    end program
