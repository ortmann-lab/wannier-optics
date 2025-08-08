
! These file contains additional potentials that are needed for the construction
! of the Hamiltonian.
!
! Author: Konrad Merkel


    subroutine generatepot_anderson(sites_p_len, p_id, enList, WW)
        ! generate random Anderson onsite-disorder
        USE Geometry
        IMPLICIT NONE
        INTEGER j, cadd, p_id, sites_p_len,r,imax
        !INTEGER sites_p(sites_p_len)
        DOUBLE PRECISION WW
        DOUBLE COMPLEX enList(sites_p_len)
        REAL randf

        if (p_id == 0) then
            write(*,*) 'Generate Anderson potential: disorder strength = ', WW
        endif

        imax=100   ! if imax=1 corresponds to uniform random distribution of onsite energies
        do j=1,(p_id+1)*100
            call random_number(randf)
        enddo
        do cadd=1,sites_p_len
            do r=1,imax
            call random_number(randf)
                enList(cadd)=sqrt(12.d0/imax)*WW*(randf-0.5)+enList(cadd)
            enddo
        enddo
    end subroutine generatepot_anderson


    subroutine coulomb_density_density(sites_p_len, p_id, sites_p, enList, gamma)  ! TODO: update
        ! Generates coupling between valence and conduction sites (electron-hole interaction).
        ! Only density-density interaction that is read from the COULOMB file is used. Density-density matrix
        ! elements that are further away are calculated using monopole-monopole interaction.
        ! Parameter: interaction strength WW (screening/amplification)
        USE Geometry
        USE MPI
        IMPLICIT NONE
        INTEGER i,j,k, ii, jj, kk,m, cadd, p_id, sites_p_len, ALLOC_ERR
        INTEGER Sx,Sy,Sz, stat, num_diagonal, ierr
        INTEGER RD(3), Rc(3), Rv(3), S1(3), S2(3), c1,c2,v1,v2
        INTEGER sites_p(sites_p_len), NumMatEl
        INTEGER, allocatable, DIMENSION(:,:) :: CoulombMap
        character (len=500) :: line
        DOUBLE PRECISION WW,dist,deltax,deltay,deltaz
        DOUBLE PRECISION max_dist_coulomb, val1, val2, gamma, ARB_TO_EV
        DOUBLE PRECISION xo, yo, zo, xi, yi, zi
        DOUBLE COMPLEX enList(sites_p_len)
        DOUBLE COMPLEX, allocatable :: CoulombOnsite(:)
        LOGICAL LFI_ONS, already_counted

        ARB_TO_EV = 14.39964535  ! factor e^2/(4*pi*eps_0) in units of eV*Angstrom

        ! READ COULOMB file
        inquire(FILE='COULOMB',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'COULOMB does not exist. Stop.'
            stop
        endif
        open(114,FILE='COULOMB',ACTION='READ')
        read(114,*) NumMatEl, max_dist_coulomb, WW

        if (abs(WW).lt.1e-10) then
            IF ( p_id .eq. 0 ) write(*,*) "WW = 0 --> calculation without coulomb interaction"
            close(114)
            return
        endif

        if(p_id.eq.0) then
            write(*,*) ''
            write(*,*) 'Generate Coulomb potential: WW = ', WW   
        endif

        IF ( p_id .eq. 0 ) write(*,*) "--> Use Monopole-Monopole interaction for distances > ", max_dist_coulomb
        ! IF ( p_id .eq. 0 ) write(*,*) "Expect ", NumMatEl, " Coulomb matrix elements in file"

        allocate(CoulombMap(NumMatEl, 5), STAT = ALLOC_ERR )
        IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION CoulombMap !!!"
        CoulombMap = 0

        allocate(CoulombOnsite(NumMatEl), STAT = ALLOC_ERR )
        IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION CoulombOnsite !!!"
        CoulombOnsite = (0.d0, 0.d0)

        i = 1  ! counter for lines in file
        j = 1  ! counter for onsite terms
        read_loop: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
            if (i.gt.NumMatEl) STOP "More elements in COULOMB file than expected. Please check the first line of the file!"

            ! store data in arrays
            read (line, *, iostat=stat) c1, c2, v1, v2, (RD(ii), ii=1,3), (Rc(jj), jj=1,3), (Rv(kk), kk=1,3), val1,val2
   
            if ( stat /= 0 ) then
                  write(*,*) "line ", i
                  STOP "Error while reading local field effect matrix element from line."
            endif
            i=i+1  ! only count lines (not used matrix elements like for LFE)

            do ii=1,3
               S1(ii) = -RD(ii)
               S2(ii) = Rc(ii)-RD(ii)-Rv(ii)
            enddo

            ! filter for onsite elements
            if ((c1.eq.c2).and.(v1.eq.v2).and.(S1(1).eq.S2(1)).and.(S1(2).eq.S2(2)).and.(S1(3).eq.S2(3))) THEN
                CoulombMap(j, 1) = c1
                CoulombMap(j, 2) = v1
                CoulombMap(j, 3) = S1(1)
                CoulombMap(j, 4) = S1(2)
                CoulombMap(j, 5) = S1(3)
                CoulombOnsite(j) = - cmplx(val1, val2)* WW / gamma  ! minus because of hamiltonian

                ! check for hermiticity
                if (ABS(IMAG(CoulombOnsite(j))).gt.1e-8) then
                    write(*,*) "Entry in file:", c1, c2, v1, v2, (RD(ii), ii=1,3), (Rc(jj), jj=1,3), (Rv(kk), kk=1,3)
                    write(*,*) "Resulting matrix element:", CoulombMap(j,:), CoulombOnsite(j)
                    stop "Density-density-Coulomb matrix element has non-vanishing imaginary part!"
                endif

                j = j+1
            endif

        end do read_loop
        close(114)
        num_diagonal = j-1

        if ((i-1).ne.NumMatEl) STOP "Different number of elements in COULOMB file than expected. &
            &Please check the first line of the file!"

        if (p_id.eq.0) then
            write(*,*) "--> Total number of COULOMB matrix elements   = ", i-1
            write(*,*) "--> Number of diagonal matrix elements        = ", num_diagonal
        endif

! #ifdef DEBUG
!         if (p_id.eq.0) then
!             write(*,*) "density-density Coulomb matrix elements (from file, might contain duplicates):"
!             do i=1, num_diagonal
!                 write(*,*) (CoulombMap(i, ii), ii=1,5), CoulombOnsite(i)
!             enddo
!         endif
! #endif

        i=0  ! number of matrix elements that are actually used
        ii=0 ! number of matrix elements that are actually used (without duplicates)
        j=0  ! number of elements with monopole-monopole interaction

        ! set Coulomb matrix elements in the Hamiltonian (as onsite energies)
        do cadd=1,sites_p_len  ! go through all sites in the domain

            m = sites_p(cadd)  ! get global site index from local site index
            call getExcitonIndex_shifted(m, Sx,Sy,Sz,v1,c1)

            ! S is shift vector of unit cells; shift applies to the outer (valence) WF
            ! Attention RD=-S For the calculation of Monopole-Monopole we need RD
            ! call getCartCoordinates(-Sx+1,-Sy+1,-Sz+1,v1, ABC,bas_o,NBATouter, xo,yo,zo)  ! valence
            call getCartCoordinates_val(-Sx+1,-Sy+1,-Sz+1,v1, xo,yo,zo)  ! valence
            ! call getCartCoordinates(1,1,1,c1, ABC,bas_i,NBATinner, xi,yi,zi)  ! conduction
            call getCartCoordinates_con(1,1,1,c1, xi,yi,zi)  ! conduction

            ! outer index --> valence band
            ! inner index --> conduction band

            deltax = xi - xo
            deltay = yi - yo
            deltaz = zi - zo
            dist=sqrt(deltax**2+deltay**2+deltaz**2)  ! distance beween conduction WF and shifted valence WF

            if(abs(dist).gt.max_dist_coulomb) then  ! large distance is approximated by monopole-monopole interaction

                enList(cadd)=-1*WW/gamma /dist * ARB_TO_EV + enList(cadd) ! comensurable with PB
                j=j+1

                ! if (dist.lt.100.) write(*,*) c1,c1,v1,v1,-Sx,-Sy,-Sz,0,0,0,0,0,0,ARB_TO_EV/dist, 0.0, " Monopole_manuell"

            else ! look up value from COULOMB array
                already_counted = .false.  ! makes sure we only count ones even if the list might contain duplicates
                do k=1,num_diagonal  ! go through list of onsite coulomb matrix elements (might contain duplicates --> do not use cycle!)

                    ! site index are equal ?
                    if ((CoulombMap(k,1) .eq. c1) .and. (CoulombMap(k,2) .eq. v1)) then

                        ! shift vector is the same ?
                        if ((CoulombMap(k,3) .eq. Sx) .and. (CoulombMap(k,4) .eq. Sy) .and. &
                            & (CoulombMap(k,5) .eq. Sz)) then

                            ! write(*,*) "set onsite:", CoulombMap(k,3), CoulombMap(k,4), CoulombMap(k,5),&
                            !     & CoulombOnsite(k), -1*WW/gamma /dist * ARB_TO_EV

                            enList(cadd)=enList(cadd) + REAL(CoulombOnsite(k))  ! signs, WW and gamma are already contained
                            
                            i=i+1  ! number of matrix elements that are actually used
                            if (.not.already_counted) then
                                ii=ii+1  ! number of matrix elements that are actually used
                                already_counted = .true.
                            endif

                        endif
                    endif
                enddo
            endif
        enddo

        ! some final checks
        if (ii+j .ne. sites_p_len) THEN
            STOP "Some sites are without any Coulomb interaction (density-density)!!!&
                & Please make sure that the maximal distance is in agreement with the given matrix elements in the COULOMB file."
        endif

        ! check if all onsite energies are used
        call MPI_Reduce(i, jj, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (p_id.eq.0) THEN
            if (jj.ne.num_diagonal) THEN
                write(*,*) ""
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) jj , "of", num_diagonal, "excitonic onsite energies (COULOMB) are used!!!"
                write(*,*) "Please make sure that all known onsite energies are used. "
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) ""
            else
                write(*,*) "All Coulomb onsite energies are used :D"
                write(*,*) ""
            endif
        endif

        deallocate(CoulombMap)
        deallocate(CoulombOnsite)

    end subroutine coulomb_density_density


    subroutine localFieldEffects_onsite(sites_p_len, p_id, sites_p, enList,  gamma)
        ! Adds local field effects to the Hamiltonian.
        ! Only matrix elements for S=S', c1=c2 and v1=v2 are used (only onsite energy is modified).
        ! This is not necessarily a good approximation for local field effects since the terms
        ! c1 != c2 and v1 != v2 might contribute substantially (maybe even more than the diagonal ones!)
        USE Geometry
        USE MPI
        IMPLICIT NONE
        INTEGER cadd, p_id, sites_p_len,i, m, k, ALLOC_ERR
        INTEGER Sx,Sy,Sz,lo,li, stat, num_diagonal, jj, ierr
        INTEGER sites_p(sites_p_len), NumMatEl
        INTEGER, allocatable, DIMENSION(:,:) :: CoulombMap
        character (len=500) :: line
        !DOUBLE PRECISION dist,deltax,deltay,deltaz
        DOUBLE PRECISION WW
        DOUBLE PRECISION max_dist_coulomb, val1, val2, gamma
        DOUBLE COMPLEX enList(sites_p_len)
        DOUBLE COMPLEX, allocatable :: CoulombTI(:)
        LOGICAL LFI_ONS

        ! READ LOCALFIELDEFFECTS file
        inquire(FILE='LOCALFIELDEFFECTS',exist=LFI_ONS)
        if(.not.LFI_ONS) then
            write(*,*) 'LOCALFIELDEFFECTS does not exist. Stop.'
            stop
        endif
        open(114,FILE='LOCALFIELDEFFECTS',ACTION='READ')
        read(114,*) NumMatEl, max_dist_coulomb, WW

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

        allocate(CoulombMap(13,NumMatEl), STAT = ALLOC_ERR )
        IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION CoulombMap !!!"

        allocate(CoulombTI(NumMatEl), STAT = ALLOC_ERR )
        IF ( ALLOC_ERR .NE. 0 ) STOP "ERROR - ALLOCATION CoulombTI !!!"

        i = 1
        read_loop: do
            read (114, '(A)', iostat=stat)  line

            if ( stat /= 0 )  exit read_loop             ! check end of file or other errors
            if ( line .eq. ' ') cycle read_loop          ! skip empty lines
            if (index (line, "#")/= 0)  cycle read_loop  ! skip comment lines
            if (i.gt.NumMatEl) STOP "More elements in LOCALFIELDEFFECTS file than expected. "&
                &"Please check the first line of the file!"

            ! store data in arrays
            read (line, *, iostat=stat) CoulombMap(1,i),CoulombMap(2,i),CoulombMap(3,i), CoulombMap(4,i),CoulombMap(5,i),&
                & CoulombMap(6,i), CoulombMap(7,i), CoulombMap(8,i), CoulombMap(9,i), CoulombMap(10,i), &
                & CoulombMap(11,i), CoulombMap(12,i), CoulombMap(13,i), val1,val2

            if ( stat /= 0 ) then
                write(*,*) "line ", i
                STOP "Error while reading local field effect matrix element from line."
            endif

            ! local field effects have prefactor of (+2) that is not contained in the file
            CoulombTI(i) = 2* cmplx(val1, val2)* WW / gamma
            i=i+1

        end do read_loop
        close(114)

        if ((i-1).ne.NumMatEl) STOP "Different number of elements in LOCALFIELDEFFECTS file than expected. &
            &Please check the first line of the file!"

        if (p_id.eq.0) write(*,*) "--> Total number of LFE matrix elements = ", i-1

! #ifdef DEBUG
!         if (p_id.eq.0) then
!             write(*,*) "Local field effect matrix elements (from file):"
!             do i=1, NumMatEl
!                 write(*,*) CoulombMap(1,i),CoulombMap(2,i),CoulombMap(3,i), CoulombMap(4,i),CoulombMap(5,i) ,&
!                     & CoulombMap(6,i), CoulombMap(7,i), CoulombMap(8,i), CoulombMap(9,i), CoulombMap(10,i), &
!                     & CoulombMap(11,i), CoulombMap(12,i), CoulombMap(13,i), CoulombTI(i)
!             enddo
!         endif
! #endif

        ! count diagonal elements (c1=c2 & v1=v2 & S1=S2)
        i = 0
        do k=1,NumMatEl  ! go through list of coulomb matrix elements

            ! only diagonal matrix elements
            if (CoulombMap(1,k) .ne. CoulombMap(2,k)) then      ! c1=c2
                cycle
            endif
            if (CoulombMap(3,k) .ne. CoulombMap(4,k)) then      ! v1=v2
                cycle
            endif

            ! RD == 0 (this is just an artefact because we use the same format as for the Coulomb matrix elements)
            if ((CoulombMap(5,k) .ne. 0) .or. (CoulombMap(6,k) .ne. 0) .or. (CoulombMap(7,k) .ne. 0)) then
                write(*,*) "[Error] line", k, "in LOCALFIELDEFFECTS:"
                stop "RD has to be zero for local field effects."
            endif

            ! S1 == S2
            if ((CoulombMap(8,k) .ne. CoulombMap(11,k)) .or. (CoulombMap(9,k) .ne. CoulombMap(12,k)) .or. &
                & (CoulombMap(10,k) .ne. CoulombMap(13,k))) then
                cycle
            endif

            ! if (p_id.eq.0) write(*,*) "diagonal:", CoulombMap(1,k), CoulombMap(3,k), CoulombMap(8,k), &
            !     & CoulombMap(9,k), CoulombMap(10,k), " = ", CoulombTI(k) * gamma /WW /2

            i=i+1
        enddo
        num_diagonal = i
        if (p_id.eq.0) write(*,*) "--> Number of diagonal matrix elements  = ", num_diagonal

        i=0  ! number of matrix elements that are actually used

        ! set matrix elements in the Hamiltonian (as onsite energies)
        do cadd=1,sites_p_len  ! go through all sites in the domain

            m = sites_p(cadd)  ! get global site index from local site index
            call getExcitonIndex_shifted(m, Sx,Sy,Sz,lo,li)

            do k=1,NumMatEl  ! go through list of coulomb matrix elements

                ! only density-density matrix elements
                if (CoulombMap(1,k) .ne. CoulombMap(2,k)) then      ! c1=c2
                    cycle
                endif
                if (CoulombMap(3,k) .ne. CoulombMap(4,k)) then      ! v1=v2
                    cycle
                endif

                ! S1 == S2
                if ((CoulombMap(8,k) .ne. CoulombMap(11,k)) .or. (CoulombMap(9,k) .ne. CoulombMap(12,k)) .or. &
                    & (CoulombMap(10,k) .ne. CoulombMap(13,k))) then
                    cycle
                endif

                if (ABS(IMAG(CoulombTI(k))).gt.1e-8) then
                    write(*,*) "Diagonal matrix elements for local field effects has non-vanishing imaginary part: ", &
                        & CoulombTI(k)
                    stop
                endif

                ! look for matrix elements of the current site
                if ((CoulombMap(1,k) .eq. li) .and. (CoulombMap(3,k) .eq. lo)) then ! site index are equal ?
                    ! outer index --> valence band
                    ! inner index --> conduction band

                    ! shift vector is the same ?
                    if ((CoulombMap(8,k) .eq. Sx) .and. (CoulombMap(9,k) .eq. Sy) .and. &
                        & (CoulombMap(10,k) .eq. Sz)) then                           ! shift vector is the same ?

                        ! write(*,*) p_id, "set FE element: ", CoulombMap(1,k), CoulombMap(3,k), CoulombMap(8,k), &
                        !     & CoulombMap(9,k), CoulombMap(10,k), " = ", CoulombTI(k)

                        ! set onsite energy
                        enList(cadd)=enList(cadd) + CoulombTI(k)  ! prefactor and sign are already contained

                        i=i+1  ! number of matrix elements that are actually used
                    endif
                endif
            enddo
            !endif
        enddo

        ! check if all onsite energies are used
        call MPI_Reduce(i, jj, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (p_id.eq.0) THEN
            if (jj.ne.num_diagonal) THEN
                write(*,*) ""
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) jj , "of", num_diagonal, "excitonic onsite energies (LFE) are used!!!"
                write(*,*) "Please make sure that all known onsite energies are used. "
                write(*,*) ""
                write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                write(*,*) ""
                write(*,*) ""
            else
                write(*,*) "All LFE onsite energies are used :D"
                write(*,*) ""
            endif
        endif

        deallocate(CoulombMap)
        deallocate(CoulombTI)

    end subroutine localFieldEffects_onsite


    subroutine rand_gauss(rand_f1, rand_f2)
        ! create normal distribued random numbers with zero mean and
        ! variance one. For the calculation we use the polar methode
        ! to get a normal distribution out of a uniform distribution.
        !
        ! return two random numbers (rand_f1, rand_f2)
        implicit none
        DOUBLE PRECISION, intent(out)     :: rand_f1, rand_f2
        
        !internal
        DOUBLE PRECISION :: u, v, q, p

        do 
            call random_number(u) ! uniform random numbers
            call random_number(v)
            
            u = 2.0 * u -1
            v = 2.0 * v -1
            q = u*u + v*v
            
            if ((q.gt.1e-16).and.(q.le.1.0)) exit
        end do
            
        p = sqrt(-2 * log(q)/q)
        rand_f1 = u * p
        rand_f2 = v * p

    end subroutine rand_gauss
