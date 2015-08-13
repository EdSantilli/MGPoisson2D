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
    type(box_data)              :: soln

    integer                     :: r             ! Current refinement level
    integer, parameter          :: maxr = 3      ! Max refinement level
    real(dp), dimension(maxr)   :: errnorm       ! Error norm at each level
    real(dp), dimension(maxr-1) :: rate          ! Convergence rates


    ! Test 1: Uniform Dirichlet BCs
    errnorm = bogus_val
    do r = 1, maxr
        call define_domain (valid, r)
        call define_box_data (soln, valid, 1, 1)
        call fill_soln (soln)
        errnorm(r) = test_uniform_diri_bcs (soln)
        call undefine_box_data (soln)
    enddo
    call compute_conv_rate (rate, errnorm)
    print*, 'Test 1: Uniform Dirichlet BCs Convergence rates = '
    print*, rate

    ! integer         :: ilo, ihi, jlo, jhi
    ! real(dp)        :: dx, dy
    ! type(box)       :: valid
    ! type(box_data)  :: state
    ! type(bdry_data) :: diri_bc, neum_bc

    ! integer :: nx, ny, i, j
    ! real(dp) :: L, H, val
    ! real(dp), dimension(:), allocatable :: x, y

    ! ! Set up domain
    ! ilo = 1
    ! ihi = 32
    ! jlo = 1
    ! jhi = 32
    ! dx = one
    ! dy = one
    ! call define_box (valid, ilo, ihi, jlo, jhi, dx, dy)

    ! nx = ihi-ilo+1
    ! ny = jhi-jlo+1
    ! L = nx*dx
    ! H = ny*dy

    ! ! Define coordinates
    ! allocate(x(ilo-1:ihi+1))
    ! x(ilo-1:ihi+1) = (/ ((i + half) * dx, i = ilo-1, ihi+1) /)

    ! allocate(y(jlo-1:jhi+1))
    ! y(jlo-1:jhi+1) = (/ ((j + half) * dy, j = jlo-1, jhi+1) /)

    ! ! Set up field
    ! call define_box_data (state, valid, 1, 1)
    ! state%data = three


    ! ! TEST 1: Dirichlet BCs ----------------------------------------------------
    ! print*, 'Testing ArrayUtils::fill_ghosts with Dirichlet BCs...'

    ! ! Define BCs
    ! call define_bdry_data (diri_bc, valid, &
    !                        BCTYPE_DIRI, &   ! xlo
    !                        BCTYPE_DIRI, &   ! xhi
    !                        BCTYPE_DIRI, &   ! ylo
    !                        BCTYPE_DIRI, &   ! yhi
    !                        BCMODE_UNIFORM, &    ! xlo
    !                        BCMODE_UNIFORM, &    ! xhi
    !                        BCMODE_UNIFORM, &    ! ylo
    !                        BCMODE_UNIFORM)      ! yhi
    ! diri_bc%data_xlo(1) = five
    ! diri_bc%data_xhi(1) = six
    ! diri_bc%data_ylo(1) = seven
    ! diri_bc%data_yhi(1) = eight

    ! ! Set BCs
    ! call fill_ghosts (state, diri_bc, .false.)

    ! ! Test valid cells
    ! do j = jlo, jhi
    !     do i = ilo, ihi
    !         if (state%data(i,j) .ne. three) then
    !             print*, 'i = ', i
    !             print*, 'j = ', j
    !             print*, 'state%data(i,j) = ', state%data(i,j)
    !             stop
    !         endif
    !     enddo
    ! enddo

    ! ! Test ghost cells
    ! i = ilo-1
    ! do j = jlo, jhi
    !     if (state%data(i,j) .ne. seven) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! i = ihi+1
    ! do j = jlo, jhi
    !     if (state%data(i,j) .ne. nine) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jlo-1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. eleven) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jhi+1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. thirteen) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! print*, '...passed.'


    ! ! TEST 2: Neumann BCs ----------------------------------------------------
    ! print*, 'Testing ArrayUtils::fill_ghosts with Neumann BCs...'

    ! ! Define BCs
    ! call define_bdry_data (neum_bc, valid, &
    !                        BCTYPE_NEUM, &   ! xlo
    !                        BCTYPE_NEUM, &   ! xhi
    !                        BCTYPE_NEUM, &   ! ylo
    !                        BCTYPE_NEUM, &   ! yhi
    !                        BCMODE_UNIFORM, &    ! xlo
    !                        BCMODE_UNIFORM, &    ! xhi
    !                        BCMODE_UNIFORM, &    ! ylo
    !                        BCMODE_UNIFORM)      ! yhi
    ! neum_bc%data_xlo(1) = five
    ! neum_bc%data_xhi(1) = six
    ! neum_bc%data_ylo(1) = seven
    ! neum_bc%data_yhi(1) = eight

    ! ! Set BCs
    ! call fill_ghosts (state, neum_bc, .true.)

    ! ! Test valid cells
    ! do j = jlo, jhi
    !     do i = ilo, ihi
    !         if (state%data(i,j) .ne. three) then
    !             print*, 'i = ', i
    !             print*, 'j = ', j
    !             print*, 'state%data(i,j) = ', state%data(i,j)
    !             stop
    !         endif
    !     enddo
    ! enddo

    ! ! Test ghost cells
    ! i = ilo-1
    ! do j = jlo, jhi
    !     if (state%data(i,j) .ne. three) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! i = ihi+1
    ! do j = jlo, jhi
    !     if (state%data(i,j) .ne. three) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jlo-1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. three) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jhi+1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. three) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! print*, '...passed.'


    ! ! TEST 3: Non-uniform Dirichlet BCs ----------------------------------------
    ! print*, 'Testing ArrayUtils::fill_ghosts with non-uniform Dirichlet BCs...'

    ! ! Set valid region
    ! do j = jlo, jhi
    !     do i = ilo, ihi
    !         state%data(i,j) = sin((half + four*y(j)/H)*pi*x(i)/L)
    !     enddo
    ! enddo

    ! ! Define BCs
    ! call define_bdry_data (diri_bc, valid, &
    !                        BCTYPE_DIRI, &   ! xlo
    !                        BCTYPE_DIRI, &   ! xhi
    !                        BCTYPE_DIRI, &   ! ylo
    !                        BCTYPE_DIRI, &   ! yhi
    !                        BCMODE_NONUNIFORM, &    ! xlo
    !                        BCMODE_NONUNIFORM, &    ! xhi
    !                        BCMODE_NONUNIFORM, &    ! ylo
    !                        BCMODE_NONUNIFORM)      ! yhi
    ! diri_bc%data_xlo(jlo:jhi) = (/ (sin((half + four*y(j)/H)*pi*x(ilo-1)/L), j = jlo, jhi) /)
    ! diri_bc%data_xhi(jlo:jhi) = (/ (sin((half + four*y(j)/H)*pi*x(ihi+1)/L), j = jlo, jhi) /)
    ! diri_bc%data_ylo(ilo:ihi) = (/ (sin((half + four*y(jlo-1)/H)*pi*x(i)/L), i = ilo, ihi) /)
    ! diri_bc%data_yhi(ilo:ihi) = (/ (sin((half + four*y(jhi+1)/H)*pi*x(i)/L), i = ilo, ihi) /)

    ! ! Set BCs
    ! call fill_ghosts (state, diri_bc, .false.)

    ! ! Test valid cells
    ! do j = jlo, jhi
    !     do i = ilo, ihi
    !         val = sin((half + four*y(j)/H)*pi*x(i)/L)
    !         if (state%data(i,j) .ne. val) then
    !             print*, 'i = ', i
    !             print*, 'j = ', j
    !             print*, 'state%data(i,j) = ', state%data(i,j)
    !             stop
    !         endif
    !     enddo
    ! enddo

    ! ! Test ghost cells
    ! i = ilo-1
    ! do j = jlo, jhi
    !     val = sin((half + four*y(j)/H)*pi*x(i)/L)
    !     if (state%data(i,j) .ne. val) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         print*, 'expected ', val
    !         stop
    !     endif
    ! enddo

    ! i = ihi+1
    ! do j = jlo, jhi
    !     if (state%data(i,j) .ne. sin((half + four*y(j)/H)*pi*x(i)/L)) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jlo-1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. sin((half + four*y(j)/H)*pi*x(i)/L)) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! j = jhi+1
    ! do i = ilo, ihi
    !     if (state%data(i,j) .ne. sin((half + four*y(j)/H)*pi*x(i)/L)) then
    !         print*, 'i = ', i
    !         print*, 'j = ', j
    !         print*, 'state%data(i,j) = ', state%data(i,j)
    !         stop
    !     endif
    ! enddo

    ! print*, '...passed.'


    ! ! Free memory
    ! call undefine_bdry_data (diri_bc)
    ! call undefine_box_data (state)
    ! if (allocated(y)) deallocate(y)
    ! if (allocated(x)) deallocate(x)


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
    ! This fills soln's valid _and_ ghost cells.
    ! --------------------------------------------------------------------------
    subroutine fill_soln (soln)
        implicit none
        type(box_data), intent(inout) :: soln
        integer                       :: ilo, ihi, jlo, jhi
        real(dp)                      :: dx, dy
        integer                       :: i, j
        real(dp)                      :: x, y

        ilo = soln%bx%ilo
        ihi = soln%bx%ihi
        jlo = soln%bx%jlo
        jhi = soln%bx%jhi

        dx = soln%valid%dx
        dy = soln%valid%dy

        do j = jlo, jhi
            y = (j + half) * dy
            do i = ilo, ihi
                x = (i + half) * dx
                soln%data(i,j) = sin((half + four*y/H)*pi*x/L)
            enddo
        enddo
    end subroutine fill_soln


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
            rate(r) = log(errnorm(r+1) / errnorm(r)) / log(two)
        enddo

    end subroutine compute_conv_rate


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_uniform_diri_bcs (soln) result (res)
        implicit none

        real(dp)                   :: res
        type(box_data), intent(in) :: soln

        integer               :: i,j
        integer               :: ilo, ihi, jlo, jhi
        real(dp)              :: dx, dy
        real(dp), dimension(:), allocatable :: x, y

        type(box_data)        :: state
        type(bdry_data)       :: diri_bc, neum_bc


        ilo = soln%valid%ilo
        ihi = soln%valid%ihi
        jlo = soln%valid%jlo
        jhi = soln%valid%jhi

        dx = soln%valid%dx
        dy = soln%valid%dy


        ! ! Define coordinates
        ! allocate(x(ilo-1:ihi+1))
        ! x(ilo-1:ihi+1) = (/ ((i + half) * dx, i = ilo-1, ihi+1) /)

        ! allocate(y(jlo-1:jhi+1))
        ! y(jlo-1:jhi+1) = (/ ((j + half) * dy, j = jlo-1, jhi+1) /)

        ! ! Define BCs
        ! call define_bdry_data (diri_bc, valid, &
        !                        BCTYPE_DIRI, &   ! xlo
        !                        BCTYPE_DIRI, &   ! xhi
        !                        BCTYPE_DIRI, &   ! ylo
        !                        BCTYPE_DIRI, &   ! yhi
        !                        BCMODE_UNIFORM, &    ! xlo
        !                        BCMODE_UNIFORM, &    ! xhi
        !                        BCMODE_UNIFORM, &    ! ylo
        !                        BCMODE_UNIFORM)      ! yhi
        ! diri_bc%data_xlo(1) = five
        ! diri_bc%data_xhi(1) = six
        ! diri_bc%data_ylo(1) = seven
        ! diri_bc%data_yhi(1) = eight

        ! ! Set BCs
        ! call fill_ghosts (state, diri_bc, .false.)

        ! ! ! Test valid cells
        ! do j = jlo, jhi
        !     do i = ilo, ihi
        !         if (state%data(i,j) .ne. three) then
        !             print*, 'i = ', i
        !             print*, 'j = ', j
        !             print*, 'state%data(i,j) = ', state%data(i,j)
        !             stop
        !         endif
        !     enddo
        ! enddo

        ! ! Test ghost cells
        ! i = ilo-1
        ! do j = jlo, jhi
        !     if (state%data(i,j) .ne. seven) then
        !         print*, 'i = ', i
        !         print*, 'j = ', j
        !         print*, 'state%data(i,j) = ', state%data(i,j)
        !         stop
        !     endif
        ! enddo

        ! i = ihi+1
        ! do j = jlo, jhi
        !     if (state%data(i,j) .ne. nine) then
        !         print*, 'i = ', i
        !         print*, 'j = ', j
        !         print*, 'state%data(i,j) = ', state%data(i,j)
        !         stop
        !     endif
        ! enddo

        ! j = jlo-1
        ! do i = ilo, ihi
        !     if (state%data(i,j) .ne. eleven) then
        !         print*, 'i = ', i
        !         print*, 'j = ', j
        !         print*, 'state%data(i,j) = ', state%data(i,j)
        !         stop
        !     endif
        ! enddo

        ! j = jhi+1
        ! do i = ilo, ihi
        !     if (state%data(i,j) .ne. thirteen) then
        !         print*, 'i = ', i
        !         print*, 'j = ', j
        !         print*, 'state%data(i,j) = ', state%data(i,j)
        !         stop
        !     endif
        ! enddo



        ! ! Free memory
        ! ! call undefine_bdry_data (diri_bc)
        ! ! call undefine_box_data (state)
        ! if (allocated(y)) deallocate(y)
        ! if (allocated(x)) deallocate(x)

        res = zero
    end function test_uniform_diri_bcs

end program test
