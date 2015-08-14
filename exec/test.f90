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

    integer, parameter          :: base_ilo = 1
    integer, parameter          :: base_ihi = 32
    integer, parameter          :: base_jlo = 1
    integer, parameter          :: base_jhi = 32

    real(8), parameter          :: L = one
    real(8), parameter          :: H = one

    type(box)                   :: valid

    integer                     :: r             ! Current refinement level
    integer, parameter          :: maxr = 7      ! Max refinement level
    real(dp), dimension(maxr)   :: errnorm       ! Error norm at each level
    real(dp), dimension(maxr-1) :: rate          ! Convergence rates

    integer, parameter          :: norm_type = 2


    ! Test 1: Non-uniform Dirichlet BCs
    errnorm = bogus_val
    do r = 1, maxr
        call define_domain (valid, r)
        errnorm(r) = test_nonuniform_diri_bcs (valid)
    enddo
    call compute_conv_rate (rate, errnorm)
    print*, 'Test 1: Non-uniform Dirichlet BCs'
    print*, 'Error norm                rate'
    print*, errnorm(1)
    do r = 2, maxr
        print*, errnorm(r), rate(r-1)
    enddo
    print*, ''


    ! Test 2: Non-uniform Neumann BCs
    errnorm = bogus_val
    do r = 1, maxr
        call define_domain (valid, r)
        errnorm(r) = test_nonuniform_neum_bcs (valid)
    enddo
    call compute_conv_rate (rate, errnorm)
    print*, 'Test 2: Non-uniform Neumann BCs'
    print*, 'Error norm                rate'
    print*, errnorm(1)
    do r = 2, maxr
        print*, errnorm(r), rate(r-1)
    enddo


contains

    ! --------------------------------------------------------------------------
    ! Define the valid region at the current level of refinement.
    ! --------------------------------------------------------------------------
    subroutine define_domain (valid, r)
        implicit none

        type(box), intent(out) :: valid
        integer, intent(in)    :: r

        integer                :: ilo, ihi, jlo, jhi
        integer                :: nx, ny
        real(dp)               :: dx, dy

        ilo = (base_ilo-1) * two**(r-1) + 1
        ihi = base_ihi * two**(r-1)
        jlo = (base_jlo-1) * two**(r-1) + 1
        jhi = base_jhi * two**(r-1)

        nx = ihi - ilo + 1
        ny = jhi - jlo + 1

        dx = L / nx
        dy = H / ny

        call define_box (valid, ilo, ihi, jlo, jhi, dx, dy)
    end subroutine define_domain


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine compute_conv_rate (rate, errnorm)
        implicit none
        real(dp), intent(out), dimension(:) :: rate
        real(dp), intent(in), dimension(:)  :: errnorm
        integer                             :: r, maxr

        maxr = size(errnorm)
        if (size(rate) .lt. maxr-1) then
            print*, 'compute_conv_rate: rate needs at least size(errnorm)-1 elements'
        endif

        rate = bogus_val
        do r = 1, maxr-1
            rate(r) = log(errnorm(r) / errnorm(r+1)) / log(two)
        enddo

    end subroutine compute_conv_rate


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_nonuniform_diri_bcs (valid) result (res)
        implicit none

        real(dp)              :: res
        type(box), intent(in) :: valid

        integer               :: i,j
        integer               :: ilo, ihi, jlo, jhi
        real(dp)              :: dx, dy
        real(dp)              :: x, y

        type(box_data)        :: soln, state
        type(bdry_data)       :: diri_bc
        real(dp)              :: val


        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Set up field
        call define_box_data (soln, valid, 1, 1)
        soln%data = three

        call define_box_data (state, valid, 1, 1)
        do j = jlo-1, jhi+1
            y = (j + half) * dy
            do i = ilo-1, ihi+1
                x = (i + half) * dx
                soln%data(i,j) = sin((half + four*y/H)*pi*x/L)
            enddo
        enddo
        state%data = bogus_val
        state%data(ilo:ihi, jlo:jhi) = soln%data(ilo:ihi, jlo:jhi)

        ! Define BCs
        call define_bdry_data (diri_bc, valid, &
                               BCTYPE_DIRI, &   ! xlo
                               BCTYPE_DIRI, &   ! xhi
                               BCTYPE_DIRI, &   ! ylo
                               BCTYPE_DIRI, &   ! yhi
                               BCMODE_NONUNIFORM, &    ! xlo
                               BCMODE_NONUNIFORM, &    ! xhi
                               BCMODE_NONUNIFORM, &    ! ylo
                               BCMODE_NONUNIFORM)      ! yhi

        x = ilo * dx
        do j = jlo, jhi
            y = (j + half) * dy
            diri_bc%data_xlo(j) = sin((half + four*y/H)*pi*x/L)
        enddo

        x = (ihi + 1) * dx
        do j = jlo, jhi
            y = (j + half) * dy
            diri_bc%data_xhi(j) = sin((half + four*y/H)*pi*x/L)
        enddo

        y = jlo * dy
        do i = ilo, ihi
            x = (i + half) * dx
            diri_bc%data_ylo(i) = sin((half + four*y/H)*pi*x/L)
        enddo

        y = (jhi + 1) * dy
        do i = ilo, ihi
            x = (i + half) * dx
            diri_bc%data_yhi(i) = sin((half + four*y/H)*pi*x/L)
        enddo

        ! Set BCs
        call fill_ghosts (state, diri_bc, .false.)
        state%data(ilo-1,jlo-1) = soln%data(ilo-1,jlo-1)
        state%data(ilo-1,jhi+1) = soln%data(ilo-1,jhi+1)
        state%data(ihi+1,jlo-1) = soln%data(ihi+1,jlo-1)
        state%data(ihi+1,jhi+1) = soln%data(ihi+1,jhi+1)

        ! state%data(ilo-1, jlo:jhi) = soln%data(ilo-1, jlo:jhi)
        ! state%data(ihi+1, jlo:jhi) = soln%data(ihi+1, jlo:jhi)
        ! state%data(ilo:ihi, jlo-1) = soln%data(ilo:ihi, jlo-1)
        ! state%data(ilo:ihi, jhi+1) = soln%data(ilo:ihi, jhi+1)

        ! Test valid cells
        do j = jlo, jhi
            do i = ilo, ihi
                val = soln%data(i,j)
                if (state%data(i,j) .ne. val) then
                    print*, 'i = ', i
                    print*, 'j = ', j
                    print*, 'state%data(i,j) = ', state%data(i,j)
                    print*, 'state%data(i,j) = ', val
                    stop
                endif
            enddo
        enddo

        ! Compute norm over ghosts
        state%data = state%data - soln%data
        res = pnorm (state, state%bx, norm_type)

        ! Free memory
        call undefine_bdry_data (diri_bc)
        call undefine_box_data (state)
        call undefine_box_data (soln)

    end function test_nonuniform_diri_bcs


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_nonuniform_neum_bcs (valid) result (res)
        implicit none

        real(dp)              :: res
        type(box), intent(in) :: valid

        integer               :: i,j
        integer               :: ilo, ihi, jlo, jhi
        real(dp)              :: dx, dy
        real(dp)              :: x, y

        type(box_data)        :: soln, state
        type(bdry_data)       :: neum_bc
        real(dp)              :: val


        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Set up field
        call define_box_data (soln, valid, 1, 1)
        soln%data = three

        call define_box_data (state, valid, 1, 1)
        do j = jlo-1, jhi+1
            y = (j + half) * dy
            do i = ilo-1, ihi+1
                x = (i + half) * dx
                soln%data(i,j) = sin((half + four*y/H)*pi*x/L)
            enddo
        enddo
        state%data = bogus_val
        state%data(ilo:ihi, jlo:jhi) = soln%data(ilo:ihi, jlo:jhi)

        ! Define BCs
        call define_bdry_data (neum_bc, valid, &
                               BCTYPE_NEUM, &   ! xlo
                               BCTYPE_NEUM, &   ! xhi
                               BCTYPE_NEUM, &   ! ylo
                               BCTYPE_NEUM, &   ! yhi
                               BCMODE_NONUNIFORM, &    ! xlo
                               BCMODE_NONUNIFORM, &    ! xhi
                               BCMODE_NONUNIFORM, &    ! ylo
                               BCMODE_NONUNIFORM)      ! yhi
        x = ilo * dx
        do j = jlo, jhi
            y = (j + half) * dy
            neum_bc%data_xlo(j) = (half + four*y/H) * (pi/L) * cos((half + four*y/H)*pi*x/L)
        enddo

        x = (ihi + 1) * dx
        do j = jlo, jhi
            y = (j + half) * dy
            neum_bc%data_xhi(j) = (half + four*y/H) * (pi/L) * cos((half + four*y/H)*pi*x/L)
        enddo

        y = jlo * dy
        do i = ilo, ihi
            x = (i + half) * dx
            neum_bc%data_ylo(i) = (four*pi*x/L/H) * cos((half + four*y/H)*pi*x/L)
        enddo

        y = (jhi + 1) * dy
        do i = ilo, ihi
            x = (i + half) * dx
            neum_bc%data_yhi(i) = (four*pi*x/L/H) * cos((half + four*y/H)*pi*x/L)
        enddo

        ! Set BCs
        call fill_ghosts (state, neum_bc, .false.)
        state%data(ilo-1,jlo-1) = soln%data(ilo-1,jlo-1)
        state%data(ilo-1,jhi+1) = soln%data(ilo-1,jhi+1)
        state%data(ihi+1,jlo-1) = soln%data(ihi+1,jlo-1)
        state%data(ihi+1,jhi+1) = soln%data(ihi+1,jhi+1)

        ! state%data(ilo-1, jlo:jhi) = soln%data(ilo-1, jlo:jhi)
        ! state%data(ihi+1, jlo:jhi) = soln%data(ihi+1, jlo:jhi)
        ! state%data(ilo:ihi, jlo-1) = soln%data(ilo:ihi, jlo-1)
        ! state%data(ilo:ihi, jhi+1) = soln%data(ilo:ihi, jhi+1)

        ! Test valid cells
        do j = jlo, jhi
            do i = ilo, ihi
                val = soln%data(i,j)
                if (state%data(i,j) .ne. val) then
                    print*, 'i = ', i
                    print*, 'j = ', j
                    print*, 'state%data(i,j) = ', state%data(i,j)
                    print*, 'state%data(i,j) = ', val
                    stop
                endif
            enddo
        enddo

        ! Compute norm over ghosts
        state%data = state%data - soln%data
        res = pnorm (state, state%bx, norm_type)

        ! Free memory
        call undefine_bdry_data (neum_bc)
        call undefine_box_data (state)
        call undefine_box_data (soln)

    end function test_nonuniform_neum_bcs

end program test
