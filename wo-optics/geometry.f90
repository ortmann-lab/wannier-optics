
! The geometry module contains all parameters, functions and subroutines
! that are related to the construction of the supercell, calculation of positions,
! domain decomposition and conversion of indexes. 
!
! Author: Konrad Merkel

module Geometry
    use MPI

    implicit none
    Save
    ! functions to generate / read geometry
    public :: parse_posfile, print_geometry, reshapesamplesize
    
    ! funcitons to calculate indexes
    public :: getExcitonIndex, getExcitonIndex_shifted, getGlobIdFromExcitonIndex_shifted
    public :: getProcessId, getLocalId, getGlobIdFromExcitonIndex

    ! getter
    public :: get_supercell_dimensions, getCartCoordinates_val, getCartCoordinates_con
    public :: get_N_val, get_N_con, getNumSiteInLocalDomain, getMappingLocalToGlobal, get_total_num_sites
    public :: get_volume_unitcell

    ! setter
    public :: set_geometry_debugging

    ! parameters for supercell and unitcell
    PRIVATE
    integer     :: n1,n2,n3                 ! dimensions of the supercell
    integer     :: pbc1, pbc2, pbc3         ! boundary conditions
    integer     :: N_val, N_con             ! number of valence and conduction WF
    integer     :: N_ex                     ! = N_val*N_con number of excitons in a unit cell
    integer*4   :: nsites                   ! total number of sites in the supercell
    integer     :: nxparts,nyparts,nzparts  ! decomposition in local domains
    double precision, dimension (3,3) :: unitcell
    double precision, allocatable, dimension(:,:) :: pos_val, pos_con ! cartesian coordinates of singel WF
    logical     :: isReshaped               ! if true: dimensions of the supercell are changed to make better domain decomposition
    logical     :: DEBUGGING                ! enables extra output

    ! outer coordinates = valence
    ! innter coordinates = conduction


contains

    subroutine set_geometry_debugging(flag)
        logical, intent(in) :: flag
        DEBUGGING = flag
    end subroutine

    subroutine getNumSiteInLocalDomain(p_id, sites_p_len2)
        ! Determines the number of sites in the local domain
        implicit none
        integer, intent(OUT) :: sites_p_len2
        integer, intent(in) :: p_id
        integer :: m, part_M, psize, ierr

        call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
        if(ierr /= 0) stop 'ERROR: cannot get number of processes from MPI'

        
        if(nxparts*nyparts*nzparts.ne.psize) then
            write(*,*) 'Problem in getNumSiteInLocalDomain STOP',nxparts,nyparts,nzparts,psize
            STOP
        endif
        
        sites_p_len2=0    
        do m=1,nsites
            part_M = getProcessId(m)
            if (part_M == p_id+1) then
                sites_p_len2 = sites_p_len2 + 1
            endif
        enddo

        if (DEBUGGING) then
            write(*,*) "MPI Process: ", p_id+1, " has ", sites_p_len2, " sites in its local domain."
        endif

    end subroutine getNumSiteInLocalDomain

    subroutine getMappingLocalToGlobal(p_id, sites_p,sites_p_len2)
        ! creates mapping of localId to globalId that is stored in sites_p(...)
        implicit none
        integer, intent(IN)  :: sites_p_len2
        integer, intent(OUT) :: sites_p(sites_p_len2)
        integer :: m,padd, ierr !loop variables
        integer :: part_M,psize,p_id!,nxparts,nyparts,nzparts to be generalized to z-coordinate
        
        call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
        if(ierr /= 0) stop 'ERROR: cannot get number of processes from MPI'


        if(nxparts*nyparts*nzparts.ne.psize) then
            write(*,*) 'Problem in getMappingLocalToGlobal STOP',nxparts,nyparts,nzparts,psize
            STOP
        endif
        
        padd=0
        do m=1,nsites
            part_M = getProcessId(m)
            if (part_M == p_id+1) then
                padd = padd + 1
                sites_p(padd) = m
            endif
        enddo
        
    end subroutine getMappingLocalToGlobal

    subroutine get_supercell_dimensions(d1,d2,d3)
        ! returns the actual supercell dimensions
        integer, intent(out) :: d1, d2, d3
        d1 = n1; d2 = n2; d3 = n3;
    end subroutine

    function get_N_val() result(N_val_)
        ! returns the number of valence functions per unit cell
        integer :: N_val_
        N_val_ = N_val
    endfunction

    function get_total_num_sites() result(N)
        ! returns the total number of electron-hole pairs in the supercell
        ! N = N_ex *n1*n2*n3
        integer*4 :: N
        N = nsites
    endfunction

    function get_N_con() result(N_con_)
        ! returns the number of conduction functions per unit cell
        integer :: N_con_
        N_con_ = N_con
    endfunction

    function get_volume_unitcell() result(VOLUME)
        double precision :: VOLUME, H(3)

        H(1)=unitcell(2,1)*unitcell(3,2)-unitcell(3,1)*unitcell(2,2)
        H(2)=unitcell(3,1)*unitcell(1,2)-unitcell(1,1)*unitcell(3,2)
        H(3)=unitcell(1,1)*unitcell(2,2)-unitcell(2,1)*unitcell(1,2)
        VOLUME=H(1)*unitcell(1,3)+H(2)*unitcell(2,3)+H(3)*unitcell(3,3)
    endfunction

    subroutine parse_posfile(p_id, filename)
        ! Reads the POSFILE, parse everything and saves it as private variabes of the module
        ! Supercell dimensions will be reshaped for an optimal domain decomposition
        character(len=*),intent(in):: filename
        logical :: exists
        character*1:: charac,c1,c2,c3
        integer :: i, p_id
        double precision, allocatable, dimension(:,:) :: pos_help

        isReshaped = .False.
        nxparts = 0
        nyparts = 0
        nzparts = 0

        inquire(file=filename, exist=exists)
        if (.not.exists) then
            write(*,*) "Cannot parse POSFILE because file does not exists!"
            write(*,*) "filename: ", filename
            stop
        end if
        OPEN(UNIT=11,FILE=filename,STATUS='OLD')
        
        ! read supercell dimensions and boundary conditions
        ! set periodic boundaries         
        pbc1=0; pbc2=0; pbc3=0
        READ(11,*) n1, n2, n3, c1,c2,c3

        ! set boundary conditions
        if(c1.eq.'y') pbc1=1
        if(c2.eq.'y') pbc2=1
        if(c3.eq.'y') pbc3=1

        ! read unit cell
        READ(11,*) charac   ! comment
        DO i=1,3
            READ(11,*) unitcell(1,i),unitcell(2,i),unitcell(3,i) ! unit cell (here we also take the transpose)
        ENDDO

        ! ----------  read valence WF  --------------
        READ(11,*) N_val
        allocate (pos_val(3,N_val), pos_help(3,N_val))

        READ(11,*) charac   ! C or D or K
        IF(charac.NE.'K'.AND.charac.NE.'k'.AND.charac.NE.'C'.AND. &
            &   charac.NE.'c'.AND.charac.NE.'d'.AND.charac.NE.'D') THEN
            WRITE(*,*) 'POSFILE: Error reading coordinates for valence Wannier functions!!',charac
            STOP
        ENDIF

        DO i=1,N_val         ! read all positions
            READ(11,*) pos_help(1,i),pos_help(2,i),pos_help(3,i)
        ENDDO
        ! transform positions in cartesian coordiantes
        pos_val=0
        call CartesianCoordinates(p_id, charac, unitcell, N_val, pos_help, pos_val)
        deallocate(pos_help)



        ! ----------  read conduction WF  --------------
        READ(11,*) N_con
        allocate (pos_con(3,N_con), pos_help(3,N_con))

        READ(11,*) charac   ! C or D or K
        IF(charac.NE.'K'.AND.charac.NE.'k'.AND.charac.NE.'C'.AND. &
            &   charac.NE.'c'.AND.charac.NE.'d'.AND.charac.NE.'D') THEN
            WRITE(*,*) 'POSFILE: Error reading coordinates for conduction Wannier functions!!',charac
            STOP
        ENDIF

        DO i=1,N_con         ! read all positions
            READ(11,*) pos_help(1,i),pos_help(2,i),pos_help(3,i)
        ENDDO
        ! transform positions in cartesian coordiantes
        pos_con=0
        call CartesianCoordinates(p_id, charac, unitcell, N_con, pos_help, pos_con)
        deallocate(pos_help)


        ! reshape sample size for better domain decomposition
        call reshapesamplesize(p_id)
        N_ex=N_val*N_con
        Nsites=N_ex*n1*n2*n3

    end subroutine

    subroutine print_geometry()
        ! prints summary of the geometry
        integer :: i

        write(*,*) "Geometry of the supercell:"
        
        
        if (isReshaped) then
            write(*,*) "Dimensions (reshaped)           :", n1, n2, n3
            write(*,*) "Decomposition into local domains:", nxparts, nyparts, nzparts
        else
            write(*,*) "Dimensions (from POSFILE)       :", n1, n2, n3
            write(*,*) "Supercell is not reshaped yet!"
        endif
        
        write(*,*) "Boundary conditions             :", pbc1, pbc2, pbc3
        write(*,*) ""
        write(*,*) "Unit cell:"
        DO i=1,3
            write(*,*) unitcell(1,i),unitcell(2,i),unitcell(3,i)
        ENDDO
        write(*,*) ""
        write(*,*) "Positions of valence WF coordinates (in cartesian basis):"
        DO i=1,N_val
            write(*,*) i, ": ", pos_val(1,i),pos_val(2,i),pos_val(3,i)
        ENDDO
        write(*,*) ""


        write(*,*) "Positions of conduction WF coordinates (in cartesian basis):"
        DO i=1,N_con
            write(*,*) i, ": ", pos_con(1,i),pos_con(2,i),pos_con(3,i)
        ENDDO
        write(*,*) ""

    end subroutine

    subroutine reshapesamplesize(p_id)
        implicit none
        integer, intent(in) :: p_id

        integer :: psize,psizex,psizey,dn1,dn2,dn3
        integer :: ierr
        double precision :: psizex3D,psizey3D,psizez3D


        ! get number of processes
        call MPI_Comm_size(MPI_COMM_WORLD, psize, ierr)
        
        
        psizex3D=(dble(psize)*dble(n1)**2/dble(n2)/dble(n3))**(1.0d0/3.0d0)   
        psizey3D=(dble(psize)*dble(n2)**2/dble(n1)/dble(n3))**(1.0d0/3.0d0)
        psizez3D=(dble(psize)*dble(n3)**2/dble(n1)/dble(n2))**(1.0d0/3.0d0)

        nxparts=huge(nxparts)  ! largest possible value for that data type
        nyparts=huge(nyparts)
        nzparts=huge(nzparts)

        !  Decomposition of sample

        if ((psizex3D.lt.1) .or. ((n1.eq.1).and.(pbc1.eq.0))) then  ! check if monolayer in x-direction
        psizex3D=1
        psizey3D=sqrt(dble(psize)*dble(n2)/dble(n3))
        psizez3D=sqrt(dble(psize)*dble(n3)/dble(n2))
        endif

        if ((psizey3D.lt.1) .or. ((n2.eq.1).and.(pbc2.eq.0))) then  ! check if monolayer in y-direction
        psizey3D=1
        psizex3D=sqrt(dble(psize)*dble(n1)/dble(n3))
        psizez3D=sqrt(dble(psize)*dble(n3)/dble(n1))
        endif

        if ((psizez3D.lt.1) .or. ((n3.eq.1).and.(pbc3.eq.0))) then  ! check if monolayer in z-direction
        psizez3D=1
        psizex3D=sqrt(dble(psize)*dble(n1)/dble(n2))
        psizey3D=sqrt(dble(psize)*dble(n2)/dble(n1))
        endif
        

        if ((n1.eq.1).and.(pbc1.eq.0)) then 
        nxparts = 1
        else
        do psizex=1,psize
            if(mod(psize,psizex).eq.0) then
                if(abs(psizex3D-psizex).lt.abs(psizex3D-nxparts)) then
                    nxparts=psizex
                endif
            endif
        enddo
        
        psize=psize/nxparts
        endif

        if ((n2.eq.1).and.(pbc2.eq.0)) then
        nyparts=1
        else
        do psizey=1,psize
            if(mod(psize,psizey).eq.0) then
                if(abs(psizey3D-psizey).lt.abs(psizey3D-nyparts)) then
                    nyparts=psizey
                endif
            endif
        enddo

        psize=psize/nyparts
        endif
        
        
        if ((n3.eq.1).and.(pbc3.eq.0)) then
        nzparts = 1
        
        psize=psize*nyparts  ! undo changes of psize
        
        nyparts = psize
        psize=psize/nyparts
        else
        nzparts=psize
        psize=psize/nzparts
        endif
        
        
        
        if(psize.eq.1)then
        !  write(*,*) "decomposition of psize into : ",nxparts,nyparts,nzparts
        else
        write(*,*) "ERROR in reshapesamplesize", psize
        endif
        
        psize=nxparts*nyparts*nzparts

        do dn1=1,nxparts
        do dn2=1,nyparts
            do dn3=1,nzparts
                if(mod(n1,nxparts).ne.0.and.mod(n2,nyparts).ne.0.and.mod(n3,nzparts).ne.0) then
                    n1=n1+1
                    n2=n2+1
                    n3=n3+1
                elseif (mod(n1,nxparts).eq.0.and.mod(n2,nyparts).ne.0.and.mod(n3,nzparts).ne.0) then
                    n2=n2+1
                    n3=n3+1
                elseif (mod(n1,nxparts).ne.0.and.mod(n2,nyparts).eq.0.and.mod(n3,nzparts).ne.0) then
                    n1=n1+1
                    n3=n3+1
                elseif (mod(n1,nxparts).ne.0.and.mod(n2,nyparts).ne.0.and.mod(n3,nzparts).eq.0) then
                    n1=n1+1
                    n2=n2+1
                elseif (mod(n1,nxparts).eq.0.and.mod(n2,nyparts).eq.0.and.mod(n3,nzparts).ne.0) then
                    n3=n3+1
                elseif (mod(n1,nxparts).eq.0.and.mod(n2,nyparts).ne.0.and.mod(n3,nzparts).eq.0) then
                    n2=n2+1
                elseif (mod(n1,nxparts).ne.0.and.mod(n2,nyparts).eq.0.and.mod(n3,nzparts).eq.0) then
                    n1=n1+1
                endif
            enddo
        enddo
        enddo

        ! n1, n2 and n3 should be odd numbers if we have periodic boundary conditions
        if((mod(n1,2).eq.0) .and. (pbc1.gt.0)) n1=n1+1
        if((mod(n2,2).eq.0) .and. (pbc2.gt.0)) n2=n2+1
        if((mod(n3,2).eq.0) .and. (pbc3.gt.0)) n3=n3+1

        isReshaped = .True.  ! indicate that reshape was performed
        N_ex=N_val*N_con
        Nsites=N_ex*n1*n2*n3

    end subroutine reshapesamplesize

    subroutine CartesianCoordinates(p_id, charac, ABC, NBAT, bashelp, bas)
        ! Private helper routine to transform site positions in cartesian coordinates
        implicit none
        double precision, allocatable, dimension(:,:) :: bas, bashelp
        double precision, dimension (3,3) :: ABC
        character*1:: charac
        integer :: NBAT, p_id

        IF(charac.NE.'K'.AND.charac.NE.'k'.AND.charac.NE.'C'.AND. &
            &   charac.NE.'c'.AND.charac.NE.'d'.AND.charac.NE.'D') THEN
            WRITE(*,*) 'POSFILE reading error for conduction (inner) coordinates!!',charac
            STOP
        ENDIF

        IF (charac.EQ.'K'.OR.charac.EQ.'k'.OR.charac.EQ.'C'.OR.charac.EQ.'c') THEN
            if (p_id == 0) then
                WRITE(*,*) 'Positions given in cartesian coordinates.'
            endif
            bas(1,:)=bashelp(1,:)
            bas(2,:)=bashelp(2,:)
            bas(3,:)=bashelp(3,:)
        ELSE
            if (p_id == 0) then
                WRITE(*,*) 'Positions given in direct coordinates.'
            endif
            bas(1,:)=ABC(1,1)*bashelp(1,:)+ABC(1,2)*bashelp(2,:)+ABC(1,3)*bashelp(3,:)
            bas(2,:)=ABC(2,1)*bashelp(1,:)+ABC(2,2)*bashelp(2,:)+ABC(2,3)*bashelp(3,:)
            bas(3,:)=ABC(3,1)*bashelp(1,:)+ABC(3,2)*bashelp(2,:)+ABC(3,3)*bashelp(3,:)
        ENDIF
    end subroutine


    subroutine ijkl(iadd,i,j,k,l)  ! for inner cooridnates
        ! Calculates the unit cell indexes for inner coordinates from the index iadd
        ! of the (inner) site (equivalent to the global index in the transport
        ! simulation)
        ! For excitons use getExcitonIndex instead!!!
        !
        ! Input :   iadd   = global (inner) index
        ! Output:   i,j,k  = unit cell index in x,y,z direction
        !           l      = site index within the unit cell

        implicit none
        integer i,j,k,l,iadd
        intent(in) iadd
        intent(out) i,j,k,l

        l=mod(iadd-1,N_ex)+1  ! N_ex = number of sites in the unit cell
        k=mod((iadd-l)/N_ex,n3)+1
        j=mod(((iadd-l)/N_ex-(k-1))/n3,n2)+1
        i=(((iadd-l)/N_ex-(k-1))/n3-(j-1))/n2+1
        
        if(l.ge.1.and.l.le.N_ex.and.k.ge.1.and.k.le.n3.and. &
        & j.ge.1.and.j.le.n2.and.i.ge.1.and.i.le.n1) then
        !OK
        else
        write(*,*)iadd,i,j,k,l, 'ERROR in ijkl subroutine'
        STOP
        endif
    end subroutine ijkl


    subroutine getExcitonIndex(globId, ox,oy,oz,ol,il)
        ! Calculates the unit cell and exciton subindices from the global index of the site
        ! (globId)

        ! Input :   globId   = global index
        ! Output:   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell for outer coordinates (valence)
        !           il       = site index within the unit cell for inner coordinates (conduction)

        implicit none
        integer globId ,ox,oy,oz, ol, il, m
        intent(in) globId
        intent(out) ox,oy,oz,ol, il

        if ((globId.lt.1).or.(globId.gt.nsites)) stop "Invalid global id in getExcitonIndex"

        call ijkl(globId,ox,oy,oz,m)        ! returns lattice index and exciton index

        il = mod(m-1, N_con) + 1
        ol = (m-il) / N_con  + 1

    end subroutine getExcitonIndex

    subroutine getExcitonIndex_shifted(globId, ox,oy,oz,ol,il)
        ! Calculates the unit cell indexes from the global index of the site
        ! (globId) and shifts the unitcell vector such that (0,0,0) is in the
        ! middle of the supercell

        ! Input :   globId   = global index
        ! Output:   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell for outer coordinates (valence)
        !           il       = site index within the unit cell for inner coordinates (conduction)

        implicit none
        integer globId ,ox,oy,oz, ol, il, m
        intent(in) globId
        intent(out) ox,oy,oz,ol, il

        call ijkl(globId,ox,oy,oz,m)        ! inner cordinates

        il = mod(m-1, N_con) + 1
        ol = (m-il) / N_con  + 1

        ! shift unit cell vector
        ox = ox - int(n1/2) -1
        oy = oy - int(n2/2) -1
        oz = oz - int(n3/2) -1

    end subroutine getExcitonIndex_shifted


    subroutine getGlobIdFromExcitonIndex_shifted(ox,oy,oz,ol,il, globId)
        ! Calculates the globla index from shifted unit cell indexes for inner (i)
        ! and outer (o) coordinates.
        ! If ox,oy,oz,ix,iy,iz are not within the (global) supercell it uses
        ! the periodic boundary conditions to wrap it around
        ! Input :   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell for outer coordinates
        !           il       = site index within the unit cell for inner coordinates
        ! Output:   globId   = global index
        implicit none
        integer globId ,ox,oy,oz, ol, il, ox2,oy2,oz2
        intent(out) globId
        intent(in) ox,oy,oz,ol, il

        ! shift unit cell vector
        ox2 = ox + int(n1/2) +1
        oy2 = oy + int(n2/2) +1
        oz2 = oz + int(n3/2) +1

        call getGlobIdFromExcitonIndex(ox2,oy2,oz2,ol,il, globId)

    end subroutine getGlobIdFromExcitonIndex_shifted


    subroutine getGlobIdFromExcitonIndex(ox,oy,oz,ol,il, globId)
        ! Calculates the globla index from unit cell indexes for inner (i)
        ! and outer (o) coordinates.
        ! If ox,oy,oz,ix,iy,iz are not within the (global) supercell it uses
        ! the periodic boundary conditions to wrap it around
        ! Input :   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell for outer coordinates
        !           il       = site index within the unit cell for inner coordinates
        ! Output:   globId   = global index
        implicit none
        integer globId ,ox,oy,oz, ol, il
        integer oxw,oyw,ozw  ! warped around coordinates
        intent(out) globId
        intent(in) ox,oy,oz,ol, il

        ! use periodic boundary conditions to wrap around if necessary
        oxw=ox     ! x-direction
        if(oxw.lt.1) oxw=ox+n1
        if(oxw.gt.n1) oxw=ox-n1

        oyw=oy     ! y-direction
        if(oyw.lt.1) oyw=oy+n2
        if(oyw.gt.n2) oyw=oy-n2

        ozw=oz     ! z-direction
        if(ozw.lt.1) ozw=oz+n3
        if(ozw.gt.n3) ozw=oz-n3

        globId = il + N_con*(ol-1 +N_val*(ozw-1+n3*(oyw-1+n2*(oxw-1))))

        if (globId.gt.nsites) STOP "globId is larger than total number of sites!!!"
        if (globId.lt.1) STOP "globId is smaller than one!!!"

        ! check if periodic boundary conditions should be used
        if ((pbc1.lt.1).and.(ox.ne.oxw)) globId=-1
        if ((pbc2.lt.1).and.(oy.ne.oyw)) globId=-1
        if ((pbc3.lt.1).and.(oz.ne.ozw)) globId=-1


        ! if ((ox.ne.oxw).or.(oy.ne.oyw).or.(oz.ne.ozw)) write(*,*) "boundary: ", globId, oxw, oyw, ozw
        ! if ((globId.eq.1).and.((ox.ne.oxw).or.(oy.ne.oyw).or.(oz.ne.ozw))) write(*,*) "boundary: ", &
        !    &globId, oxw, oyw, ozw, ox, oy, oz

    end subroutine getGlobIdFromExcitonIndex

    subroutine getCartCoordinates_val(Rx,Ry,Rz,l, x,y,z)
        ! returns the cartesian coordinates (x,y,z) for the index
        ! (ox,oy,oz,ol) that corresponds to outer coordinates
        ! Input :   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell
        ! Output:   x,y,z    = cartesian coordinates
        implicit none
        integer, intent(in) :: Rx,Ry,Rz,l
        DOUBLE PRECISION, intent(out) :: x,y,z

        x = pos_val(1,l) +(Rx-1)*unitcell(1,1)+(Ry-1)*unitcell(1,2)+(Rz-1)*unitcell(1,3)
        y = pos_val(2,l) +(Rx-1)*unitcell(2,1)+(Ry-1)*unitcell(2,2)+(Rz-1)*unitcell(2,3)
        z = pos_val(3,l) +(Rx-1)*unitcell(3,1)+(Ry-1)*unitcell(3,2)+(Rz-1)*unitcell(3,3)

    end subroutine getCartCoordinates_val

    subroutine getCartCoordinates_con(Rx,Ry,Rz,l, x,y,z)
        ! returns the cartesian coordinates (x,y,z) for the index
        ! (ox,oy,oz,ol) that corresponds to outer coordinates
        ! Input :   ox,oy,oz = unit cell index in x,y,z direction
        !           ol       = site index within the unit cell
        ! Output:   x,y,z    = cartesian coordinates
        implicit none
        integer, intent(in) :: Rx,Ry,Rz,l
        DOUBLE PRECISION, intent(out) :: x,y,z

        x = pos_con(1,l) +(Rx-1)*unitcell(1,1)+(Ry-1)*unitcell(1,2)+(Rz-1)*unitcell(1,3)
        y = pos_con(2,l) +(Rx-1)*unitcell(2,1)+(Ry-1)*unitcell(2,2)+(Rz-1)*unitcell(2,3)
        z = pos_con(3,l) +(Rx-1)*unitcell(3,1)+(Ry-1)*unitcell(3,2)+(Rz-1)*unitcell(3,3)

    end subroutine getCartCoordinates_con

    function getProcessId(globSiteId) result(processId)
        ! Returns the process id of a given site in the TB-model
        ! Input : globSiteId  = global index of the site
        ! Output: processId = id of the process that contains this
        !                       site in the local domain
        integer, intent (in) :: globSiteId ! input
        integer              :: processId ! output
        integer :: Rx1,Ry1,Rz1,m1, l1
        integer :: nxlen, nylen, nzlen, dd, ff, gg

        ! get unit cell indexes
        call getExcitonIndex(globSiteId, Rx1,Ry1,Rz1,m1,l1)
        
        if (mod(n1,nxparts).eq.0) then
        nxlen=n1/nxparts
        else
        if (mod(n1-1,nxparts).ne.0) STOP 'Cannot find process for site because n1-1 is not multiple of nxparts!'
        nxlen=(n1-1)/nxparts
        if (Rx1.eq.n1) Rx1=Rx1-1
        endif

        if (mod(n2,nyparts).eq.0) then
        nylen=n2/nyparts
        else
        if (mod(n2-1,nyparts).ne.0) STOP 'Cannot find process for site because n2-1 is not multiple of nyparts!'
        nylen=(n2-1)/nyparts
        if (Ry1.eq.n2) Ry1=Ry1-1
        endif

        if (mod(n3,nzparts).eq.0) then
        nzlen=n3/nzparts
        else
        if (mod(n3-1,nzparts).ne.0) STOP 'Cannot find process for site because n3 is not multiple of nzparts!'
        nzlen=(n3-1)/nzparts
        if (Rz1.eq.n3) Rz1=Rz1-1
        endif

        dd=(Rx1-1)/nxlen
        ff=(Ry1-1)/nylen
        gg=(Rz1-1)/nzlen

        ! write(*,*) "nxlen        ", nxlen, nylen, nzlen
        ! write(*,*) "ExcitonIndex ", globSiteId, Rx1,Ry1,Rz1,m1,l1
        ! write(*,*) "dd           ", dd, ff, gg, 1+gg+ff*nzparts+dd*nyparts*nzparts

        processId = 1 + gg + nzparts*( ff + nyparts * dd)
        
    end function

    
    function getLocalId(globSiteId) result(localId)
        ! Returns the local id of a given site in the TB-model from the
        ! gloal id of the supercell.
        ! Input : globSiteId  = global index of the site
        ! Output: localId     = local index of the domain that contains this
        !                       site
        integer, intent (in) :: globSiteId ! input
        integer              :: localId    ! output
        integer :: Rx,Ry,Rz,il,ol
        integer :: nxlen, nylen, nzlen, processId, px,py,pz
        integer :: Lx, Ly, Lz
#ifdef DEBUG
        integer :: globid2, pId
#endif

        ! get unit cell indexes
        call getExcitonIndex(globSiteId, Rx,Ry,Rz,ol,il)
        processId = getProcessId(globSiteId)

        ! length of local domain without any additional boundary cells
        nxlen=n1/nxparts
        nylen=n2/nyparts
        nzlen=n3/nzparts

        ! domain coordinate in x-, y- and z-direction
        ! px,py and pz start at zero
        pz = mod(processId-1, nzparts)
        py = mod((processId-1-pz)/nzparts, nyparts)
        px = ((processId-1-pz)/nzparts-py) / nyparts

#ifdef DEBUG
        pId = 1+pz+ nzparts*(py + nyparts*px)
        if (pId.ne.processId) then
            write(*,*) nxparts, nyparts, nzparts
            write(*,*) processId, "-->",px,py,pz, "-->", pId
            STOP "Error while calcualting px,py,pz"
        endif
#endif

        ! first cell of the local domain for that process
        Lx = px*nxlen+1
        Ly = py*nylen+1
        Lz = pz*nzlen+1

#ifdef DEBUG
        call getGlobIdFromExcitonIndex(Lx,Ly,Lz,ol,il, globId2)
        pId = getProcessId(globId2)
        if (pId.ne.processId) then
            write(*,*) "nxlen:",nxlen, nylen, nzlen
            write(*,*) "processIds:", processId, "-->",px,py,pz, "-->", pId
            write(*,*) "Cell Ids:",Lx,Ly,Lz, "-->", pId, "?=", processId
            STOP "Error while calcualting first cell of domain"
        endif
#endif

        ! difference (unit cell) vector
        Rx = Rx-Lx
        Ry = Ry-Ly
        Rz = Rz-Lz

        ! from now on we need actual length of the local domain in each direction
        ! including additional boundary cells
        if (((px+1).eq.nxparts).and.(mod(n1,nxparts).ne.0)) nxlen = nxlen+1
        if (((py+1).eq.nyparts).and.(mod(n2,nyparts).ne.0)) nylen = nylen+1
        if (((pz+1).eq.nzparts).and.(mod(n3,nzparts).ne.0)) nzlen = nzlen+1

        localId = il + N_con*(ol-1 +N_val*(Rz + nzlen*( Ry + nylen*Rx )))  ! maybe without il and ol (only Rx,Ry,Rz???)

#ifdef DEBUG
        if ((localId.gt.N_con*N_val*nxlen*nylen*nzlen).or.(localId.lt.1)) then
            write(*,*) processId, localId, N_con*N_val*nxlen*nylen*nzlen
            write(*,*) Rx, Ry, Rz
            write(*,*) px,py,pz, "-->", nxlen,nylen,nzlen
            STOP "LocalId is out of range"
        endif
#endif
        
    end function

    

end module Geometry