!*******************************************************************************
!  MGPoisson2D - A simple Poisson solver with semicoarsening.
!  Developed by Ed Santilli
!  Copyright (C) 2015 Philadelphia University
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
!  USA
!
!  For up-to-date contact information, please visit the repository homepage,
!  https://github.com/EdSantilli.
!*******************************************************************************

program test
    use ArrayUtils
    implicit none

    integer         :: ilo, ihi, jlo, jhi
    real(dp)        :: dx, dy
    type(box)       :: valid
    type(box_data)  :: state
    type(bdry_data) :: diri_bc, neum_bc

    integer :: i, j

    ! Set up domain
    ilo = 1
    ihi = 32
    jlo = 1
    jhi = 32
    dx = one
    dy = one
    call define_box (valid, ilo, ihi, jlo, jhi, dx, dy)

    ! Set up field
    call define_box_data (state, valid, 1, 1)
    state%data = three


    ! TEST 1: Dirichlet BCs ----------------------------------------------------
    print*, 'Testing ArrayUtils::fill_ghosts with Dirichlet BCs...'

    ! Set up BCs
    call define_bdry_data (diri_bc, valid, &
                           BCTYPE_DIRI, &   ! xlo
                           BCTYPE_DIRI, &   ! xhi
                           BCTYPE_DIRI, &   ! ylo
                           BCTYPE_DIRI, &   ! yhi
                           BCMODE_UNIFORM, &    ! xlo
                           BCMODE_UNIFORM, &    ! xhi
                           BCMODE_UNIFORM, &    ! ylo
                           BCMODE_UNIFORM)      ! yhi
    diri_bc%data_xlo(1) = five
    diri_bc%data_xhi(1) = six
    diri_bc%data_ylo(1) = seven
    diri_bc%data_yhi(1) = eight

    ! Test BCs
    call fill_ghosts (state, diri_bc, .false.)

    do j = jlo, jhi
        do i = ilo, ihi
            if (state%data(i,j) .ne. three) then
                print*, 'i = ', i
                print*, 'j = ', j
                print*, 'state%data(i,j) = ', state%data(i,j)
                stop
            endif
        enddo
    enddo

    i = ilo-1
    do j = jlo, jhi
        if (state%data(i,j) .ne. seven) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    i = ihi+1
    do j = jlo, jhi
        if (state%data(i,j) .ne. nine) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    j = jlo-1
    do i = ilo, ihi
        if (state%data(i,j) .ne. eleven) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    j = jhi+1
    do i = ilo, ihi
        if (state%data(i,j) .ne. thirteen) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    print*, '...passed.'


    ! TEST 2: Neumann BCs ----------------------------------------------------
    print*, 'Testing ArrayUtils::fill_ghosts with Neumann BCs...'

    ! Set up BCs
    call define_bdry_data (neum_bc, valid, &
                           BCTYPE_NEUM, &   ! xlo
                           BCTYPE_NEUM, &   ! xhi
                           BCTYPE_NEUM, &   ! ylo
                           BCTYPE_NEUM, &   ! yhi
                           BCMODE_UNIFORM, &    ! xlo
                           BCMODE_UNIFORM, &    ! xhi
                           BCMODE_UNIFORM, &    ! ylo
                           BCMODE_UNIFORM)      ! yhi
    neum_bc%data_xlo(1) = five
    neum_bc%data_xhi(1) = six
    neum_bc%data_ylo(1) = seven
    neum_bc%data_yhi(1) = eight

    ! Test BCs
    call fill_ghosts (state, neum_bc, .true.)

    do j = jlo, jhi
        do i = ilo, ihi
            if (state%data(i,j) .ne. three) then
                print*, 'i = ', i
                print*, 'j = ', j
                print*, 'state%data(i,j) = ', state%data(i,j)
                stop
            endif
        enddo
    enddo

    i = ilo-1
    do j = jlo, jhi
        if (state%data(i,j) .ne. three) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    i = ihi+1
    do j = jlo, jhi
        if (state%data(i,j) .ne. three) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    j = jlo-1
    do i = ilo, ihi
        if (state%data(i,j) .ne. three) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    j = jhi+1
    do i = ilo, ihi
        if (state%data(i,j) .ne. three) then
            print*, 'i = ', i
            print*, 'j = ', j
            print*, 'state%data(i,j) = ', state%data(i,j)
            stop
        endif
    enddo

    print*, '...passed.'


    ! Free memory
    call undefine_bdry_data (diri_bc)
    call undefine_box_data (state)

end program test
