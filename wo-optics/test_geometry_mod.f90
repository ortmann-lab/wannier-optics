
! This module contains test routines for the geometry.
!
! Author: Konrad Merkel


module test_geometry_mod

    use Geometry
    implicit none
    PRIVATE

    public :: test_getLocalId, test_indexes, testProcessMapping
    
contains

    subroutine test_getLocalId(p_id, sites_p, sites_p_len)
        ! This is a test if subroutines getExcitonIndex and getGlobIdFromExcitonIndex
        ! are compatible with each other
        IMPLICIT NONE
        INTEGER p_id,sites_p_len, sites_p(sites_p_len)
        INTEGER cadd, globId, localId

        write(*,*) "Test globalId to localId mapping", p_id, sites_p_len

        do cadd=1, sites_p_len
            !write(*,*) "cadd ", cadd
            globId = sites_p(cadd)
            localId = getLocalId(globId)

            if (cadd.ne.localId) then
                write(*,*) cadd, localId
                STOP "Wrong localId for the own domain."
            endif
        enddo

    end subroutine test_getLocalId


    subroutine test_indexes(p_id, sites_p, sites_p_len)
        ! This is a test if subroutines getExcitonIndex and getGlobIdFromExcitonIndex
        ! are compatible with each other
        IMPLICIT NONE
        INTEGER p_id,sites_p_len, sites_p(sites_p_len), i, j, k, li, lo
        INTEGER cadd, globId, globId2

        write(*,*) "Test coordinate systems and unit cell mapping ", p_id, sites_p_len

        do cadd=1, sites_p_len
            !write(*,*) "cadd ", cadd
            globId = sites_p(cadd)
            !write(*,*) "globId ", globId

            call getExcitonIndex(globId, i,j,k,lo,li)
            call getGlobIdFromExcitonIndex(i,j,k,lo,li, globId2)

            if (globId .ne. globId2) then
                write(*,*) globId, globId2, " : ", i,j,k,lo,li
                stop "Error in index calculation: "
            endif

            ! test shifted version
            call getExcitonIndex_shifted(globId, i,j,k,lo,li)
            call getGlobIdFromExcitonIndex_shifted(i,j,k,lo,li, globId2)

            if (globId .ne. globId2) then
                write(*,*) globId, globId2, " : ", i,j,k,lo,li
                stop "Error in (shifted) index calculation!"
            endif
        enddo

    end subroutine test_indexes


    subroutine testProcessMapping(p_id, sites_p, sites_p_len, psize)
        ! tests if the getProcess function gives the correct process
        ! for every site
        IMPLICIT NONE
        integer :: p_id, locId, globId, sites_p_len, procId, psize
        integer :: sites_p(sites_p_len), ox,oy,oz,ol,il, n1,n2,n3
        intent(in) p_id,sites_p,sites_p_len, psize

        call get_supercell_dimensions(n1,n2,n3)

        write(*,*) "Test process mapping", p_id, sites_p_len
        do locId=1,sites_p_len
            globId = sites_p(locId)
            procId = getProcessId(globId)
            ! write(*,*) "test: ", p_id+1, procId, globId, locId
            
            if (p_id+1 .ne. procId) then
                write(*,*) "Process for own domain is not valid: ", p_id+1, procId
                stop
            endif
        enddo


        if (p_id.eq.0) then  ! print decomposition details
            do globId=1,n1*n2*n3
                procId = getProcessId(globId)
                ! write(*,*) "glob test: ", globId, procId
                
                if (procId .gt. psize) then
                    call getExcitonIndex(globId, ox,oy,oz,ol,il)
                    write(*,*) globId, ox,oy,oz,ol,il
                    write(*,*) "ProcessId larger than number of processes: ", globId, procId
                    stop
                endif
            enddo
        endif

    end subroutine testProcessMapping



end module test_geometry_mod