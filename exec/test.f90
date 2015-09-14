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

    integer, parameter          :: base_ilo = 0
    integer, parameter          :: base_ihi = 31
    integer, parameter          :: base_jlo = 0
    integer, parameter          :: base_jhi = 31

    real(8), parameter          :: L = one
    real(8), parameter          :: H = one

    integer                     :: r             ! Current refinement level
    integer, parameter          :: maxr = 8      ! Max refinement level
    real(dp), dimension(maxr)   :: errnorm       ! Error norm at each level
    real(dp), dimension(maxr-1) :: rate          ! Convergence rates

    integer, parameter          :: norm_type = 2

    ! type(box_data) :: x, y
    type(box)                       :: valid
    type(geo_data), dimension(maxr) :: geo
    real(dp) :: t1, t2


    ! Geometry setup
    call cpu_time (t1)
    do r = 1, maxr
        call define_domain (valid, r)
        call define_geometry (geo(r), valid)
    enddo
    call cpu_time (t2)
    print*, 'Geometry setup: elapsed time (s) = ', t2-t1
    print*


    ! ! Test 1: Geometry
    ! errnorm = bogus_val
    ! do r = 1, maxr
    !     errnorm(r) = test_geometry (geo(r))
    ! enddo
    ! call compute_conv_rate (rate, errnorm)
    ! print*, 'Test 1: geometry'
    ! print*, 'Error norm                rate'
    ! print*, errnorm(1)
    ! do r = 2, maxr
    !     print*, errnorm(r), rate(r-1)
    ! enddo
    ! print*


    ! ! Test 2: Non-uniform Dirichlet BCs
    ! errnorm = bogus_val
    ! do r = 1, maxr
    !     errnorm(r) = test_nonuniform_diri_bcs (geo(r))
    ! enddo
    ! call compute_conv_rate (rate, errnorm)
    ! print*, 'Test 2: Non-uniform Dirichlet BCs'
    ! print*, 'Error norm                rate'
    ! print*, errnorm(1)
    ! do r = 2, maxr
    !     print*, errnorm(r), rate(r-1)
    ! enddo
    ! print*


    ! ! Test 3: Non-uniform Neumann BCs
    ! errnorm = bogus_val
    ! do r = 1, maxr
    !     errnorm(r) = test_nonuniform_neum_bcs (geo(r))
    ! enddo
    ! call compute_conv_rate (rate, errnorm)
    ! print*, 'Test 3: Non-uniform Neumann BCs'
    ! print*, 'Error norm                rate'
    ! print*, errnorm(1)
    ! do r = 2, maxr
    !     print*, errnorm(r), rate(r-1)
    ! enddo

    ! ! Test 4: Staggered partial derivatives
    ! errnorm = bogus_val
    ! do r = 1, maxr
    !     errnorm(r) = test_derivatives (geo(r))
    ! enddo
    ! call compute_conv_rate (rate, errnorm)
    ! print*, 'Test 4: Staggered partial derivatives'
    ! print*, 'Error norm                rate'
    ! print*, errnorm(1)
    ! do r = 2, maxr
    !     print*, errnorm(r), rate(r-1)
    ! enddo
    ! print*

    ! Test 5: Solver test
    ! errnorm = bogus_val
    ! do r = 1, maxr
    !     errnorm(r) = test_solver (geo(r))
    ! enddo
    ! call compute_conv_rate (rate, errnorm)
    ! print*, 'Test 5: Solver test'
    ! print*, 'Error norm                rate'
    ! print*, errnorm(1)
    ! do r = 2, maxr
    !     print*, errnorm(r), rate(r-1)
    ! enddo
    ! print*
    print*, 'Test 5: Solver test on ', geo(maxr)%J%valid%nx, ' x ', geo(maxr)%J%valid%ny
    errnorm(maxr) = test_solver (geo(maxr))
    print*, 'Error norm = ', errnorm(maxr)
    print*



    ! Prints x and y coordinates to the terminal.
    ! call define_domain (valid, 1)
    ! call define_box_data (x, valid, 0, 0, BD_NODE, BD_NODE)
    ! call define_box_data (y, valid, 0, 0, BD_NODE, BD_NODE)

    ! call fill_x (x)
    ! call fill_y (y)

    ! print*, 'x = ', x%data(:,x%bx%jlo)
    ! print*, 'y = ', y%data(:,y%bx%jlo)

    ! call undefine_box_data (x)
    ! call undefine_box_data (y)


    ! Free memory
    do r = 1, maxr
        call undefine_geometry (geo(r))
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

        ilo = base_ilo * two**(r-1)
        ihi = (base_ihi+1) * two**(r-1) - 1
        jlo = base_jlo * two**(r-1)
        jhi = (base_jhi+1) * two**(r-1) - 1

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
    ! Allocate and compute various geometric quantities.
    ! --------------------------------------------------------------------------
    subroutine define_geometry (geo, valid)
        implicit none

        type(geo_data), intent(out) :: geo
        type(box), intent(in)       :: valid

        type(box)                   :: ccbx
        integer, parameter          :: ng = 1

        ! Allocate space
        call define_box_data (geo%J, valid, ng, ng, BD_CELL, BD_CELL)
        call define_box_data (geo%Jgup_xx, valid, ng, ng, BD_NODE, BD_CELL)
        call define_box_data (geo%Jgup_xy, valid, ng, ng, BD_NODE, BD_CELL)
        call define_box_data (geo%Jgup_yx, valid, ng, ng, BD_CELL, BD_NODE)
        call define_box_data (geo%Jgup_yy, valid, ng, ng, BD_CELL, BD_NODE)

        ! Fill with data
        call fill_J (geo%J)
        call fill_Jgup (geo%Jgup_xx, 1, 1)
        call fill_Jgup (geo%Jgup_xy, 1, 2)
        call fill_Jgup (geo%Jgup_yx, 2, 1)
        call fill_Jgup (geo%Jgup_yy, 2, 2)
        geo%dx = valid%dx
        geo%dy = valid%dy

    end subroutine define_geometry


    ! --------------------------------------------------------------------------
    ! Frees memory used by a geo_data object.
    ! --------------------------------------------------------------------------
    subroutine undefine_geometry (geo)
        implicit none
        type(geo_data), intent(inout) :: geo

        call undefine_box_data (geo%J)
        call undefine_box_data (geo%Jgup_xx)
        call undefine_box_data (geo%Jgup_xy)
        call undefine_box_data (geo%Jgup_yx)
        call undefine_box_data (geo%Jgup_yy)
        geo%dx = bogus_val
        geo%dy = bogus_val

    end subroutine undefine_geometry


    ! --------------------------------------------------------------------------
    ! Fills a box_data with the Cartesian x cell locations.
    ! --------------------------------------------------------------------------
    subroutine fill_x (x)
        implicit none

        type(box_data), intent(inout) :: x
        real(dp)                      :: dxi1, offx
        integer                       :: ilo, ihi, i
        integer                       :: jlo, jhi, j

        ilo = x%bx%ilo
        ihi = x%bx%ihi
        jlo = x%bx%jlo
        jhi = x%bx%jhi
        dxi1 = x%bx%dx
        offx = x%offx

        ! Identity map
        x%data(:,jlo) = (/ ((i + offx) * dxi1, i = ilo, ihi) /)
        do j = jlo+1, jhi
            x%data(:,j) = x%data(:,jlo)
        enddo
    end subroutine fill_x


    ! --------------------------------------------------------------------------
    ! Fills a box_data with the Cartesian y cell locations.
    ! --------------------------------------------------------------------------
    subroutine fill_y (y)
        implicit none

        type(box_data), intent(inout) :: y

        type(box)                   :: xbox
        type(box_data)              :: x
        real(dp)                    :: dxi1, dxi2, offx, offy, L, H, el
        real(dp)                    :: yfrac
        integer                     :: ilo, ihi, i
        integer                     :: jlo, jhi, j

        ilo = y%bx%ilo
        ihi = y%bx%ihi
        jlo = y%bx%jlo
        jhi = y%bx%jhi
        dxi1 = y%bx%dx
        dxi2 = y%bx%dy
        offx = y%offx
        offy = y%offy
        L = y%L
        H = y%H

        ! Horizontal positions at bottom boundary
        call define_box (xbox, y%valid%ilo, y%valid%ihi, y%valid%jlo, y%valid%jlo, dxi1, dxi2)
        call define_box_data (x, xbox, y%ngx, 0, y%offi, BD_NODE)
        call fill_x (x)

        ! Bottom elevation
        call fill_elevation (y%data(:,jlo), x%data(:,x%valid%jlo), ilo, ihi, L, H)

        ! Map vertically. Do this in reverse order so we don't clobber
        ! the elevation data.
        do j = jhi, jlo, -1
            yfrac = (j + offy) * dxi2 / H
            do i = ilo, ihi
                el = y%data(i,jlo)
                y%data(i,j) = el + (H-el)*yfrac
            enddo
        enddo

        ! Free memory
        call undefine_box_data (x)

    end subroutine fill_y


    ! --------------------------------------------------------------------------
    ! Fills an array with bottom elevations.
    !
    !                         ///^\\\
    !                        /       \
    !                       /         \
    !                      /           \
    !                     /             \
    ! _________________///               \\\_________________
    ! |---------------|--|----|--|--|----|--|---------------|
    ! -LxMax        -C1 -C2  -C3 0  C4  C5  C6            LxMax
    !
    ! /// Represents smoothed regions.
    ! For a symmetric hill,
    !  C1 = -l*cos(angle) - 3P
    !  C2 = -l*cos(angle) - P
    !  C3 = -B
    !  C4 = B
    !  C5 = l*cos(angle) + P
    !  C6 = l*cos(angle) + 3P
    !  Where B and P are the sizes of the smoothed regions.
    ! LxMax = Lx/2
    !
    ! The total width of the ridge is then [2*lp*cos(angle) + 4*Bp + 2*Pp]*Lx
    !  where Lx is the domain length and Hz is the domain height
    ! --------------------------------------------------------------------------
    subroutine fill_elevation (el, x, ilo, ihi, Lx, Hz)
        implicit none

        integer, intent(in)                       :: ilo, ihi
        real(dp), intent(out), dimension(ilo:ihi) :: el
        real(dp), intent(in), dimension(ilo:ihi)  :: x
        real(dp), intent(in)                      :: Lx, Hz

        ! Masoud's lab-scale ridge (total ridge width ~ 3.8m with Lx = 40.5m)
        real(dp), parameter :: angle = 19.2877 * pi / 180.0
        real(dp), parameter :: lp = 0.009714_dp   ! Length of along-slope critical region   (fraction of the domain width)
        real(dp), parameter :: Bp = 0.01173_dp    ! Length of smoothed region at ridge base (fraction of the domain width)
        real(dp), parameter :: Pp = 0.0183542_dp  ! Half-width of smoothed region at ridge peak (fraction of the domain width)

        real(dp) :: l, B, P, F, scale
        real(dp) :: C1, C2, C3, C4, C5, C6
        real(dp) :: sa, ca, ta
        real(dp) :: LxMax, lstar
        real(dp) :: b0, b1, b2    ! Coefficents of smoothed regions at ridge base
        real(dp) :: p0, p2        ! Coefficents of smoothed regions at ridge peak

        integer  :: i
        real(dp) :: xx

        ! These will be used again and again. Cache the values.
        sa = sin(angle)
        ca = cos(angle)
        ta = tan(angle)

        ! lp, Bp, and Pp are fractions of the domain width.
        ! We need to bring them to scale
        l = lp*Lx
        B = Bp*Lx
        P = Pp*Lx

        ! Length of the critical slope extended to make a triangle shaped ridge
        lstar = l + (B+P)/ca

        ! Locations of functional changes
        C1 = -lstar*ca - B
        C2 = -lstar*ca + B
        C3 = -P
        C4 = P
        C5 = lstar*ca - B
        C6 = lstar*ca + B

        ! Cubic spline coefficients
        b0 = fourth*ta*(B+lstar*ca)*(B+lstar*ca)/B
        b1 = -half*ta*(B+lstar*ca)/B
        b2 = fourth*ta/B

        p0 = lstar*sa - half*ta*P
        p2 = -half*ta/P

        do i = ilo, ihi
            xx = x(i)

            ! el(i) = zero ! Temporary, flat bottom.
            ! el(i) = half * xx ! Sloped bottom
            el(i) = Hz*half*(one-xx/Lx) ! Sloped bottom (Use this one)
            ! el(i) = fourth * sin(pi*xx/Lx)   ! Temporary, sinusoidal bottom.

            ! if (xx .le. C1) then
            !     ! Left flat region
            !     el(i) = zero

            ! else if ((C1 .lt. xx) .and. (xx .lt. C2)) then
            !     ! Curving upwards (left ridge base)
            !     el(i) = b2*xx*xx - b1*xx + b0

            ! else if ((C2 .le. xx) .and. (xx .le. C3)) then
            !     ! Left critical region
            !     el(i) = lstar*sa + ta*xx

            ! else if ((C3 .lt. xx) .and. (xx .lt. C4)) then
            !     ! Curving downwards (ridge peak)
            !     el(i) = p2*xx*xx + p0

            ! else if ((C4 .le. xx) .and. (xx .le. C5)) then
            !     ! Right critical region
            !     el(i) = lstar*sa - ta*xx

            ! else if ((C5 .lt. xx) .and. (xx .lt. C6)) then
            !     ! Curving upwards (right ridge base)
            !     el(i) = b2*xx*xx + b1*xx + b0

            ! else
            !     ! Right flat region
            !     el(i) = zero

            ! endif
        enddo
    end subroutine fill_elevation


    ! --------------------------------------------------------------------------
    ! Fills a box_data with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    ! --------------------------------------------------------------------------
    subroutine fill_dxdXi (dest, mu, nu)
        implicit none

        type(box_data), intent(inout) :: dest
        integer, intent(in)           :: mu, nu

        type(box_data) :: xmu
        integer        :: ngx, ngy, ilo, ihi, jlo, jhi, i, j, offi, offj
        real(dp)       :: scale, xi, eta
        type(box)      :: cc_valid

        ilo = dest%bx%ilo
        ihi = dest%bx%ihi
        jlo = dest%bx%jlo
        jhi = dest%bx%jhi
        ngx = dest%ngx
        ngy = dest%ngy
        offi = dest%offi
        offj = dest%offj

        ! ! Shortcuts...
        ! if ((mu .eq. 1) .and. (nu .eq. 1)) then
        !     dest%data = one
        ! else if ((mu .eq. 1) .and. (nu .eq. 2)) then
        !     dest%data = zero
        ! endif

        ! Construct the cell-centered box needed by the define functions.
        cc_valid = dest%valid
        if (dest%offi .eq. BD_NODE) then
            cc_valid%ihi = cc_valid%ihi - 1
            cc_valid%nx = cc_valid%nx - 1
        endif
        if (dest%offj .eq. BD_NODE) then
            cc_valid%jhi = cc_valid%jhi - 1
            cc_valid%ny = cc_valid%ny - 1
        endif

        if (nu .eq. 1) then
            ! Allocate
            call define_box_data (xmu, cc_valid, ngx+1, ngy, stagger(offi), offj)

            ! Fill with x or y
            if (mu .eq. 1) then
                call fill_x (xmu)
            else
                call fill_y (xmu)
            endif

            ! Differentiate
            scale = one / dest%valid%dx
            if (offi .eq. BD_CELL) then
                dest%data(ilo:ihi,:) = (xmu%data(ilo+1:ihi+1,:) - xmu%data(ilo:ihi,:)) * scale
            else
                dest%data(ilo:ihi,:) = (xmu%data(ilo:ihi,:) - xmu%data(ilo-1:ihi-1,:)) * scale
            endif

        else
            ! Allocate
            call define_box_data (xmu, cc_valid, ngx, ngy+1, offi, stagger(offj))

            ! Fill with x or y
            if (mu .eq. 1) then
                call fill_x (xmu)
            else
                call fill_y (xmu)
            endif

            ! Differentiate
            scale = one / dest%valid%dy
            if (offj .eq. BD_CELL) then
                dest%data(:,jlo:jhi) = (xmu%data(:,jlo+1:jhi+1) - xmu%data(:,jlo:jhi)) * scale
            else
                dest%data(:,jlo:jhi) = (xmu%data(:,jlo:jhi) - xmu%data(:,jlo-1:jhi-1)) * scale
            endif
        endif

        ! Free memory
        call undefine_box_data (xmu)

    end subroutine fill_dxdXi


    ! --------------------------------------------------------------------------
    ! Fills a box_data with J = det[Jacobian] = (dy/dEta)
    ! --------------------------------------------------------------------------
    subroutine fill_J (dest)
        implicit none
        type(box_data), intent(inout) :: dest
        call fill_dxdXi (dest, 2, 2)
    end subroutine fill_J


    ! --------------------------------------------------------------------------
    ! Fills a box_data with the contravariant metric elements
    ! gup^{mu,nu} = Sum over rho [ dXi^{mu}/dx^{rho} * dXi^{nu}/dx^{rho} ]
    ! --------------------------------------------------------------------------
    subroutine fill_gup (dest, mu, nu)
        implicit none

        type(box_data), intent(inout) :: dest
        integer, intent(in)           :: mu, nu
        type(box_data)                :: tmp

        ! Don't do anything special here. Just remove the J from Jgup.
        call define_box_data (tmp, dest)

        call fill_J (tmp)
        call fill_Jgup (dest, mu, nu)
        dest%data = dest%data / tmp%data

        call undefine_box_data (tmp)

    end subroutine fill_gup


    ! --------------------------------------------------------------------------
    ! Fills a box_data with detJ * gup
    ! --------------------------------------------------------------------------
    subroutine fill_Jgup (dest, mu, nu)
        implicit none

        type(box_data), intent(inout) :: dest
        integer, intent(in)           :: mu, nu
        type(box_data)                :: tmp

        if (mu .eq. 1) then
            if (nu .eq. 1) then
                ! Jgup^{1,1} = y_eta
                call fill_dxdxi (dest, 2, 2)
            else
                ! Jgup^{1,2} = -y_xi
                call fill_dxdxi (dest, 2, 1)
                dest%data = -dest%data
            endif
        else
            if (nu .eq. 1) then
                ! Jgup^{2,1} = -y_xi
                call fill_dxdxi (dest, 2, 1)
                dest%data = -dest%data
            else
                ! Jgup^{2,2} = (1 + y_xi**2) / y_eta
                call define_box_data (tmp, dest)
                call fill_dxdxi (dest, 2, 1)
                call fill_dxdxi (tmp, 2, 2)
                dest%data = (one + dest%data**2) / tmp%data
                call undefine_box_data (tmp)
            endif
        endif

    end subroutine fill_Jgup


    ! --------------------------------------------------------------------------
    ! This only works with a constant-slope bottom.
    ! --------------------------------------------------------------------------
    function test_geometry (geo) result (res)
        implicit none

        real(dp)                   :: res
        type(geo_data), intent(in) :: geo

        type(box)                  :: valid
        integer                    :: i,j
        integer                    :: ilo, ihi, jlo, jhi
        real(dp)                   :: dx, dy, L, H
        type(box_data), target     :: bdx, bdy
        type(box_data), target     :: bdx_x, bdx_y
        type(box_data), target     :: bdy_x, bdy_y
        real(dp), dimension(:,:), pointer :: xp, yp

        type(box_data)             :: state, soln, soln_x, soln_y

        valid = geo%J%valid

        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        L = geo%J%L
        H = geo%J%H

        ! Compute Cartesian locations
        call define_box_data (bdx, valid, 0, 0, BD_CELL, BD_CELL)
        call define_box_data (bdx_x, valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (bdx_y, valid, 0, 0, BD_CELL, BD_NODE)

        call define_box_data (bdy, bdx)
        call define_box_data (bdy_x, valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (bdy_y, valid, 0, 0, BD_CELL, BD_NODE)

        bdx%data(:,jhi) = (/ ((i+half)*dx, i=ilo,ihi) /)
        bdx_x%data(:,jhi) = (/ (i*dx, i=ilo,ihi+1) /)
        bdx_y%data(:,jhi+1) = bdx%data(:,jhi)
        do j = jlo, jhi
            bdx%data(:,j) = bdx%data(:,jhi)
            bdx_x%data(:,j) = bdx_x%data(:,jhi)
            bdx_y%data(:,j) = bdx_y%data(:,jhi+1)
        enddo

        bdy%data(ihi,:) = (/ ((j+half)*dy, j=jlo,jhi) /)
        bdy_x%data(ihi+1,:) = bdy%data(ihi,:)
        bdy_y%data(ihi,:) = (/ (j*dy, j=jlo,jhi+1) /)
        do i = ilo, ihi
            bdy%data(i,:) = bdy%data(ihi,:)
            bdy_x%data(i,:) = bdy_x%data(ihi+1,:)
            bdy_y%data(i,:) = bdy_y%data(ihi,:)
        enddo

        ! Allocate the solution spaces
        call define_box_data (soln, valid, 0, 0, BD_CELL, BD_CELL)
        call define_box_data (soln_x, valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (soln_y, valid, 0, 0, BD_CELL, BD_NODE)

        ! CC tests...

        ! ! Test x
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_x (state)
        ! xp => bdx%data
        ! ! yp => bdy%data
        ! soln%data = xp
        ! nullify (xp)
        ! ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test y
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_y (state)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = yp + half*H*(one-xp/L)*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dxdXi
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_dxdxi (state, 1, 1)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = one
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dxdEta
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_dxdxi (state, 1, 2)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = zero
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dydXi
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_dxdxi (state, 2, 1)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = (-H*half/L)*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dy/dEta
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_dxdxi (state, 2, 2)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = half + half*xp/L
        ! nullify (xp)
        ! nullify (yp)
        ! ! do i = ilo,ihi
        ! !     do j = jlo,jhi
        ! !         print*, soln%data(i,j), state%data(i,j)
        ! !     enddo
        ! ! enddo
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test J = (L+Xi) / (2L)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = half + half*xp/L
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data(ilo:ihi,jlo:jhi) = geo%J%data(ilo:ihi,jlo:jhi) - soln%data(ilo:ihi,jlo:jhi)
        ! res = pnorm (soln, valid, norm_type)

        ! ! Test Jgup^{1,1}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_Jgup (state, 1, 1)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = half*(one+xp/L)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{1,2}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_Jgup (state, 1, 2)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{2,1}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_Jgup (state, 2, 1)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{2,2}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_CELL)
        ! call fill_Jgup (state, 2, 2)
        ! xp => bdx%data
        ! yp => bdy%data
        ! soln%data = (one + (half*H/L*(one-yp/H))**2) / (half+half*xp/L)
        ! nullify (xp)
        ! nullify (yp)
        ! soln%data = state%data - soln%data
        ! res = pnorm (soln, valid, norm_type)
        ! call undefine_box_data (state)


        ! FC tests...

        ! ! Test dxdXi
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_dxdxi (state, 1, 1)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = one
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dxdEta
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_dxdxi (state, 1, 2)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = zero
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dydXi
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_dxdxi (state, 2, 1)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = (-H*half/L)*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test dy/dEta
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_dxdxi (state, 2, 2)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = half + half*xp/L
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test J = (L+Xi) / (2L)
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_J (state)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = half + half*xp/L
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test J = (L+Xi) / (2L)
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_NODE)
        ! call fill_J (state)
        ! xp => bdx_y%data
        ! yp => bdy_y%data
        ! soln_y%data = half + half*xp/L
        ! nullify (xp)
        ! nullify (yp)
        ! soln_y%data = state%data - soln_y%data
        ! res = pnorm (soln_y, valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{1,1}
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = half*(one+xp/L)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data(ilo:ihi+1,jlo:jhi) = geo%Jgup_xx%data(ilo:ihi+1,jlo:jhi) - soln_x%data(ilo:ihi+1,jlo:jhi)
        ! res = pnorm (soln_x, soln_x%valid, norm_type)

        ! ! Test Jgup^{1,2}
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data(ilo:ihi+1,jlo:jhi) = geo%Jgup_xy%data(ilo:ihi+1,jlo:jhi) - soln_x%data(ilo:ihi+1,jlo:jhi)
        ! res = pnorm (soln_x, soln_x%valid, norm_type)

        ! ! Test Jgup^{2,1}
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_Jgup (state, 2, 1)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{2,2}
        ! call define_box_data (state, valid, 0, 0, BD_NODE, BD_CELL)
        ! call fill_Jgup (state, 2, 2)
        ! xp => bdx_x%data
        ! yp => bdy_x%data
        ! soln_x%data = (one + (half*H/L*(one-yp/H))**2) / (half+half*xp/L)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_x%data = state%data - soln_x%data
        ! res = pnorm (soln_x, soln_x%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{1,1}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_NODE)
        ! call fill_Jgup (state, 1, 1)
        ! xp => bdx_y%data
        ! yp => bdy_y%data
        ! soln_y%data = half*(one+xp/L)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_y%data = state%data - soln_y%data
        ! res = pnorm (soln_y, soln_y%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{1,2}
        ! call define_box_data (state, valid, 0, 0, BD_CELL, BD_NODE)
        ! call fill_Jgup (state, 1, 2)
        ! xp => bdx_y%data
        ! yp => bdy_y%data
        ! soln_y%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_y%data = state%data - soln_y%data
        ! res = pnorm (soln_y, soln_y%valid, norm_type)
        ! call undefine_box_data (state)

        ! ! Test Jgup^{2,1}
        ! xp => bdx_y%data
        ! yp => bdy_y%data
        ! soln_y%data = half*H/L*(one-yp/H)
        ! nullify (xp)
        ! nullify (yp)
        ! soln_y%data(ilo:ihi,jlo:jhi+1) = geo%Jgup_yx%data(ilo:ihi,jlo:jhi+1) - soln_y%data(ilo:ihi,jlo:jhi+1)
        ! res = pnorm (soln_y, soln_y%valid, norm_type)

        ! Test Jgup^{2,2}
        xp => bdx_y%data
        yp => bdy_y%data
        soln_y%data = (one + (half*H/L*(one-yp/H))**2) / (half+half*xp/L)
        nullify (xp)
        nullify (yp)
        soln_y%data(ilo:ihi,jlo:jhi+1) = geo%Jgup_yy%data(ilo:ihi,jlo:jhi+1) - soln_y%data(ilo:ihi,jlo:jhi+1)
        res = pnorm (soln_y, soln_y%valid, norm_type)



        ! Free memory
        call undefine_box_data (soln)
        call undefine_box_data (soln_x)
        call undefine_box_data (soln_y)

        call undefine_box_data (bdx)
        call undefine_box_data (bdx_x)
        call undefine_box_data (bdx_y)

        call undefine_box_data (bdy)
        call undefine_box_data (bdy_x)
        call undefine_box_data (bdy_y)

    end function test_geometry


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_nonuniform_diri_bcs (geo) result (res)
        implicit none

        real(dp)                   :: res
        type(geo_data), intent(in) :: geo

        type(box)                  :: valid
        integer                    :: i,j
        integer                    :: ilo, ihi, jlo, jhi
        real(dp)                   :: dx, dy
        type(box_data)             :: bdx, bdy
        type(box_data), target     :: bdx_xlo, bdx_xhi, bdx_ylo, bdx_yhi
        type(box_data), target     :: bdy_xlo, bdy_xhi, bdy_ylo, bdy_yhi
        real(dp), dimension(:), pointer :: xp, yp

        type(box_data)             :: soln, state
        type(bdry_data)            :: diri_bc
        real(dp)                   :: val

        valid = geo%J%valid

        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Compute Cartesian locations
        call define_box_data (bdx, valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data_bdry (bdx_xlo, bdx, 1, SIDE_LO)
        call define_box_data_bdry (bdx_xhi, bdx, 1, SIDE_HI)
        call define_box_data_bdry (bdx_ylo, bdx, 2, SIDE_LO)
        call define_box_data_bdry (bdx_yhi, bdx, 2, SIDE_HI)

        call define_box_data (bdy, bdx)
        call define_box_data_bdry (bdy_xlo, bdy, 1, SIDE_LO)
        call define_box_data_bdry (bdy_xhi, bdy, 1, SIDE_HI)
        call define_box_data_bdry (bdy_ylo, bdy, 2, SIDE_LO)
        call define_box_data_bdry (bdy_yhi, bdy, 2, SIDE_HI)

        call fill_x (bdx)
        call fill_x (bdx_xlo)
        call fill_x (bdx_xhi)
        call fill_x (bdx_ylo)
        call fill_x (bdx_yhi)

        call fill_y (bdy)
        call fill_y (bdy_xlo)
        call fill_y (bdy_xhi)
        call fill_y (bdy_ylo)
        call fill_y (bdy_yhi)

        ! Set up soln
        call define_box_data (soln, bdx)
        soln%data = sin((half + four*bdy%data/H)*pi*bdx%data/L)

        ! Set up state with the true solution in the interior (valid) cells
        ! and bogus values in the ghost cells.
        call define_box_data (state, bdx)
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

        xp => bdx_xlo%data(ilo,:)
        yp => bdy_xlo%data(ilo,:)
        diri_bc%data_xlo = sin((half + four*yp/H)*pi*xp/L)

        xp => bdx_xhi%data(ihi+1,:)
        yp => bdy_xhi%data(ihi+1,:)
        diri_bc%data_xhi = sin((half + four*yp/H)*pi*xp/L)

        xp => bdx_ylo%data(:,jlo)
        yp => bdy_ylo%data(:,jlo)
        diri_bc%data_ylo = sin((half + four*yp/H)*pi*xp/L)

        xp => bdx_yhi%data(:,jhi+1)
        yp => bdy_yhi%data(:,jhi+1)
        diri_bc%data_yhi = sin((half + four*yp/H)*pi*xp/L)

        nullify(xp)
        nullify(yp)

        ! Fill state's ghost cells.
        call fill_ghosts (state, diri_bc, geo, .false.)

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
        res = gpnorm (state, norm_type)

        ! Free memory
        call undefine_bdry_data (diri_bc)
        call undefine_box_data (state)
        call undefine_box_data (soln)

        call undefine_box_data (bdx)
        call undefine_box_data (bdx_xlo)
        call undefine_box_data (bdx_xhi)
        call undefine_box_data (bdx_ylo)
        call undefine_box_data (bdx_yhi)

        call undefine_box_data (bdy)
        call undefine_box_data (bdy_xlo)
        call undefine_box_data (bdy_xhi)
        call undefine_box_data (bdy_ylo)
        call undefine_box_data (bdy_yhi)

    end function test_nonuniform_diri_bcs


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_nonuniform_neum_bcs (geo) result (res)
        implicit none

        real(dp)                   :: res
        type(geo_data), intent(in) :: geo

        type(box)                  :: valid
        integer                    :: i,j
        integer                    :: ilo, ihi, jlo, jhi
        real(dp)                   :: dx, dy
        type(box_data)             :: bdx, bdy
        type(box_data), target     :: bdx_xlo, bdx_xhi, bdx_ylo, bdx_yhi
        type(box_data), target     :: bdy_xlo, bdy_xhi, bdy_ylo, bdy_yhi
        real(dp), dimension(:), pointer :: xp, yp

        type(box_data)             :: soln, state
        type(bdry_data)            :: neum_bc
        real(dp)                   :: val

        valid = geo%J%valid

        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Compute Cartesian locations
        call define_box_data (bdx, valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data_bdry (bdx_xlo, bdx, 1, SIDE_LO)
        call define_box_data_bdry (bdx_xhi, bdx, 1, SIDE_HI)
        call define_box_data_bdry (bdx_ylo, bdx, 2, SIDE_LO)
        call define_box_data_bdry (bdx_yhi, bdx, 2, SIDE_HI)

        call define_box_data (bdy, bdx)
        call define_box_data_bdry (bdy_xlo, bdy, 1, SIDE_LO)
        call define_box_data_bdry (bdy_xhi, bdy, 1, SIDE_HI)
        call define_box_data_bdry (bdy_ylo, bdy, 2, SIDE_LO)
        call define_box_data_bdry (bdy_yhi, bdy, 2, SIDE_HI)

        call fill_x (bdx)
        call fill_x (bdx_xlo)
        call fill_x (bdx_xhi)
        call fill_x (bdx_ylo)
        call fill_x (bdx_yhi)

        call fill_y (bdy)
        call fill_y (bdy_xlo)
        call fill_y (bdy_xhi)
        call fill_y (bdy_ylo)
        call fill_y (bdy_yhi)

        ! Set up soln
        call define_box_data (soln, bdx)
        soln%data = (bdx%data**3) * (bdy%data**3)

        ! Set up state with the true solution in the interior (valid) cells
        ! and bogus values in the ghost cells.
        call define_box_data (state, bdx)
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

        xp => bdx_xlo%data(ilo,:)
        yp => bdy_xlo%data(ilo,:)
        neum_bc%data_xlo = (three*xp**2*yp**3) * geo%Jgup_xx%data(ilo,:) &
                         + (xp**3*three*yp**2) * geo%Jgup_xy%data(ilo,:)

        xp => bdx_xhi%data(ihi+1,:)
        yp => bdy_xhi%data(ihi+1,:)
        neum_bc%data_xhi = (three*xp**2*yp**3) * geo%Jgup_xx%data(ihi+1,:) &
                         + (xp**3*three*yp**2) * geo%Jgup_xy%data(ihi+1,:)

        xp => bdx_ylo%data(:,jlo)
        yp => bdy_ylo%data(:,jlo)
        neum_bc%data_ylo = ((three*xp**2)*(yp**3)) * geo%Jgup_yx%data(:,jlo) &
                         + ((xp**3)*(three*yp**2)) * geo%Jgup_yy%data(:,jlo)

        xp => bdx_yhi%data(:,jhi+1)
        yp => bdy_yhi%data(:,jhi+1)
        neum_bc%data_yhi = (three*xp**2*yp**3) * geo%Jgup_yx%data(:,jhi+1) &
                         + (xp**3*three*yp**2) * geo%Jgup_yy%data(:,jhi+1)

        nullify(xp)
        nullify(yp)

        ! Fill state's ghost cells.
        call fill_ghosts (state, neum_bc, geo, .false.)

        state%data(ilo-1, jlo:jhi) = soln%data(ilo-1, jlo:jhi)
        state%data(ihi+1, jlo:jhi) = soln%data(ihi+1, jlo:jhi)
        ! state%data(ilo:ihi, jlo-1) = soln%data(ilo:ihi, jlo-1)
        state%data(ilo:ihi, jhi+1) = soln%data(ilo:ihi, jhi+1)

        ! Compute norm over ghosts
        state%data = state%data - soln%data
        res = gpnorm (state, norm_type)

        ! Free memory
        call undefine_bdry_data (neum_bc)
        call undefine_box_data (state)
        call undefine_box_data (soln)

        call undefine_box_data (bdx)
        call undefine_box_data (bdx_xlo)
        call undefine_box_data (bdx_xhi)
        call undefine_box_data (bdx_ylo)
        call undefine_box_data (bdx_yhi)

        call undefine_box_data (bdy)
        call undefine_box_data (bdy_xlo)
        call undefine_box_data (bdy_xhi)
        call undefine_box_data (bdy_ylo)
        call undefine_box_data (bdy_yhi)

    end function test_nonuniform_neum_bcs


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_derivatives (geo) result (res)
        use MGPoisson2D
        implicit none

        real(dp)                   :: res
        type(geo_data), intent(in) :: geo

        type(box)                  :: valid
        integer                    :: ilo, ihi, jlo, jhi, i, j
        real(dp)                   :: dx, dy
        type(box_data)             :: bdx, bdy

        type(box_data)             :: soln, state, phi

        valid = geo%J%valid

        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Compute cell-centered Cartesian locations
        call define_box_data (bdx, valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data (bdy, bdx)
        do j = jlo-1, jhi+1
            do i = ilo-1, ihi+1
                bdx%data(i,j) = (i + half) * dx
                bdy%data(i,j) = (j + half) * dy
            enddo
        enddo

        ! Define phi
        call define_box_data (phi, bdx)
        ! phi%data = (bdx%data**3) * (bdy%data**3)
        phi%data = bogus_val
        ! phi%data(:,jlo:jhi) = (bdx%data(:,jlo:jhi)**3) * (bdy%data(:,jlo:jhi)**3)
        phi%data(ilo:ihi,:) = (bdx%data(ilo:ihi,:)**3) * (bdy%data(ilo:ihi,:)**3)

        ! Compute x derivative
        call define_box_data (state, valid, 0, 0, BD_CELL, BD_NODE)
        state%data = bogus_val
        call compute_pd (state, phi, 2)
        ! state%data(ilo:ihi+1,jlo:jhi) = (phi%data(ilo:ihi+1,jlo:jhi) - phi%data(ilo-1:ihi,jlo:jhi)) / dx
        ! state%data(ilo:ihi,jlo:jhi+1) = (  phi%data(ilo+1:ihi+1,jlo:jhi+1) - phi%data(ilo-1:ihi-1,jlo:jhi+1) &
        !                                  + phi%data(ilo+1:ihi+1,jlo-1:jhi) - phi%data(ilo-1:ihi-1,jlo-1:jhi)  ) * fourth / dx

        ! Re-center Cartesian locations
        call undefine_box_data (bdx)
        call undefine_box_data (bdy)
        call define_box_data (bdx, state)
        call define_box_data (bdy, state)
        do j = jlo, jhi + 1 - state%offj
            do i = ilo, ihi + 1 - state%offi
                bdx%data(i,j) = (i + state%offx) * dx
                bdy%data(i,j) = (j + state%offy) * dy
            enddo
        enddo

        ! Define true soln
        call define_box_data (soln, state)
        soln%data = (bdx%data**3) * three*(bdy%data**2)

        ! Compute norm
        state%data = state%data - soln%data
        res = pnorm (state, state%valid, norm_type)

        ! Free memory
        call undefine_box_data (state)
        call undefine_box_data (soln)
        call undefine_box_data (phi)
        call undefine_box_data (bdx)
        call undefine_box_data (bdy)

    end function test_derivatives


    ! ! --------------------------------------------------------------------------
    ! ! --------------------------------------------------------------------------
    ! function test_solver (geo) result (res)
    !     use MGPoisson2D
    !     implicit none

    !     real(dp)                   :: res
    !     type(geo_data), intent(in) :: geo

    !     type(box)                  :: valid
    !     integer                    :: ilo, ihi, jlo, jhi, i, j
    !     real(dp)                   :: dx, dy
    !     real(dp)                   :: kx, ky
    !     logical, parameter         :: homog = .true.

    !     type(box_data)             :: bdx, bdy
    !     type(box_data)             :: soln
    !     type(bdry_data)            :: bc
    !     type(box_data)             :: lphi
    !     type(box_data)             :: invdiags
    !     type(box_data)             :: phi

    !     valid = geo%J%valid

    !     ilo = valid%ilo
    !     ihi = valid%ihi
    !     jlo = valid%jlo
    !     jhi = valid%jhi

    !     dx = valid%dx
    !     dy = valid%dy

    !     ! Compute cell-centered Cartesian locations
    !     call define_box_data (bdx, valid, 1, 1, BD_CELL, BD_CELL)
    !     call define_box_data (bdy, bdx)
    !     call fill_x (bdx)
    !     call fill_y (bdy)

    !     ! Set BCs
    !     call define_bdry_data (bc, valid, &
    !                            BCTYPE_NEUM, &   ! xlo
    !                            BCTYPE_NEUM, &   ! xhi
    !                            BCTYPE_NEUM, &   ! ylo
    !                            BCTYPE_NEUM, &   ! yhi
    !                            BCMODE_UNIFORM, &    ! xlo
    !                            BCMODE_UNIFORM, &    ! xhi
    !                            BCMODE_UNIFORM, &    ! ylo
    !                            BCMODE_UNIFORM)      ! yhi
    !     bc%data_xlo = zero
    !     bc%data_xhi = zero
    !     bc%data_ylo = zero
    !     bc%data_yhi = zero

    !     ! Set up solution
    !     kx = (valid%nx/2-2) * pi / L
    !     ky = (valid%ny/2-2) * pi / L
    !     call define_box_data (soln, valid, 1, 1, BD_CELL, BD_CELL)
    !     soln%data = cos(kx*bdx%data) * cos(ky*bdy%data)

    !     ! Set up RHS
    !     call define_box_data (lphi, valid, 0, 0, BD_CELL, BD_CELL)
    !     call compute_laplacian (lphi, soln, geo, bc, homog)
    !     ! lphi%data = -(kx**2 + ky**2) * soln%data(ilo:ihi,jlo:jhi)

    !     ! Set up invdiags
    !     call define_box_data (invdiags, lphi)
    !     call compute_inverse_diags (invdiags, geo)

    !     ! Solve
    !     call define_box_data (phi, soln)
    !     ! call relax_jacobi (phi, lphi, geo, bc, homog, invdiags, &
    !     !                    one,     & ! omega
    !     !                    1.0d-6,  & ! tol
    !     !                    100,     & ! maxiters
    !     !                    .true.)    ! zerophi

    !     call relax_gs (phi, lphi, geo, bc, homog, invdiags, &
    !                    one,     & ! omega
    !                    1.0d-6,  & ! tol
    !                    25,       & ! maxiters
    !                    .false., & ! redblack
    !                    .true.)    ! zerophi

    !     ! Compute norm
    !     phi%data = phi%data - soln%data
    !     res = pnorm (phi, phi%valid, norm_type)

    !     ! Free memory
    !     ! call undefine_box_data (r)
    !     call undefine_box_data (phi)
    !     call undefine_box_data (invdiags)
    !     call undefine_box_data (lphi)
    !     call undefine_box_data (soln)
    !     call undefine_box_data (bdx)
    !     call undefine_box_data (bdy)
    !     call undefine_bdry_data (bc)

    ! end function test_solver


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    function test_solver (geo) result (res)
        use MGPoisson2D
        implicit none

        real(dp)                   :: res
        type(geo_data), intent(in) :: geo

        type(box)                  :: valid
        integer                    :: ilo, ihi, jlo, jhi, i, j
        real(dp)                   :: dx, dy
        real(dp)                   :: kx, ky
        logical, parameter         :: homog = .true.
        integer, parameter         :: verbosity = 3

        real(dp), dimension(:), pointer :: xp, yp
        type(box_data),target      :: bdx, bdx_x, bdx_y
        type(box_data),target      :: bdy, bdy_x, bdy_y
        type(box_data)             :: soln
        type(bdry_data)            :: bc
        type(box_data)             :: lphi
        type(box_data)             :: invdiags
        type(box_data)             :: phi

        valid = geo%J%valid

        ilo = valid%ilo
        ihi = valid%ihi
        jlo = valid%jlo
        jhi = valid%jhi

        dx = valid%dx
        dy = valid%dy

        ! Compute cell-centered Cartesian locations
        call define_box_data (bdx  , valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data (bdx_x, valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (bdx_y, valid, 0, 0, BD_CELL, BD_NODE)
        call fill_x (bdx)
        call fill_x (bdx_x)
        call fill_x (bdx_y)

        call define_box_data (bdy  , valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data (bdy_x, valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (bdy_y, valid, 0, 0, BD_CELL, BD_NODE)
        call fill_y (bdy)
        call fill_y (bdy_x)
        call fill_y (bdy_y)

        ! ! Set BCs
        ! call define_bdry_data (bc, valid, &
        !                        BCTYPE_DIRI, &   ! xlo
        !                        BCTYPE_DIRI, &   ! xhi
        !                        BCTYPE_DIRI, &   ! ylo
        !                        BCTYPE_DIRI, &   ! yhi
        !                        BCMODE_NONUNIFORM, &    ! xlo
        !                        BCMODE_NONUNIFORM, &    ! xhi
        !                        BCMODE_NONUNIFORM, &    ! ylo
        !                        BCMODE_NONUNIFORM)      ! yhi

        ! kx = (valid%nx/2-2) * pi / L
        ! ky = (valid%ny/2-2) * pi / L

        ! xp => bdx_x%data(ilo,:)
        ! yp => bdy_x%data(ilo,:)
        ! bc%data_xlo = cos(kx*xp) * cos(ky*yp)

        ! xp => bdx_x%data(ihi+1,:)
        ! yp => bdy_x%data(ihi+1,:)
        ! bc%data_xhi = cos(kx*xp) * cos(ky*yp)

        ! xp => bdx_y%data(:,jlo)
        ! yp => bdy_y%data(:,jlo)
        ! bc%data_ylo = cos(kx*xp) * cos(ky*yp)

        ! xp => bdx_y%data(:,jhi+1)
        ! yp => bdy_y%data(:,jhi+1)
        ! bc%data_yhi = cos(kx*xp) * cos(ky*yp)

        ! Set BCs
        call define_bdry_data (bc, valid, &
                               BCTYPE_NEUM, &   ! xlo
                               BCTYPE_NEUM, &   ! xhi
                               BCTYPE_NEUM, &   ! ylo
                               BCTYPE_NEUM, &   ! yhi
                               BCMODE_NONUNIFORM, &    ! xlo
                               BCMODE_NONUNIFORM, &    ! xhi
                               BCMODE_NONUNIFORM, &    ! ylo
                               BCMODE_NONUNIFORM)      ! yhi

        kx = (valid%nx/2-2) * pi / L
        ky = (valid%ny/2-2) * pi / L

        xp => bdx_x%data(ilo,:)
        yp => bdy_x%data(ilo,:)
        bc%data_xlo = -kx * sin(kx*xp) * cos(ky*yp)

        xp => bdx_x%data(ihi+1,:)
        yp => bdy_x%data(ihi+1,:)
        bc%data_xhi = -kx * sin(kx*xp) * cos(ky*yp)

        xp => bdx_y%data(:,jlo)
        yp => bdy_y%data(:,jlo)
        bc%data_ylo = -ky * cos(kx*xp) * sin(ky*yp)

        xp => bdx_y%data(:,jhi+1)
        yp => bdy_y%data(:,jhi+1)
        bc%data_yhi = -ky * cos(kx*xp) * sin(ky*yp)

        nullify(xp)
        nullify(yp)

        ! Set up solution
        call define_box_data (soln, valid, 1, 1, BD_CELL, BD_CELL)
        soln%data = cos(kx*bdx%data) * cos(ky*bdy%data)

        ! Set up RHS
        call define_box_data (lphi, valid, 0, 0, BD_CELL, BD_CELL)
        call compute_laplacian (lphi, soln, geo, bc, homog)
        ! lphi%data = -(kx**2 + ky**2) * soln%data(ilo:ihi,jlo:jhi)

        ! Set up invdiags
        call define_box_data (invdiags, lphi)
        call compute_inverse_diags (invdiags, geo)

        ! Solve
        call define_box_data (phi, soln)
        ! call relax_jacobi (phi, lphi, geo, bc, homog, invdiags, &
        !                    one,     & ! omega
        !                    1.0d-6,  & ! tol
        !                    20,      & ! maxiters
        !                    .true.,  & ! zerophi
        !                    verbosity)
        ! call relax_gs (phi, lphi, geo, bc, homog, invdiags, &
        !                one,     & ! omega
        !                1.0d-6,  & ! tol
        !                25,      & ! maxiters
        !                .true.,  & ! redblack
        !                .true.,  & ! zerophi
        !                verbosity)
        ! call solve_bicgstab (phi, lphi, geo, bc, homog, &
        !                1.0d-6,  & ! tol
        !                80,      & ! max iters
        !                5,       & ! max restarts
        !                .true.,  & ! zerophi
        !                verbosity)
        call vcycle (phi, lphi, geo, bc, homog, 0, 0, &
                     1.0d-6,  & ! tol
                     5,       & ! max iters
                     -1,      & ! max depth
                     1,       & ! num cycles
                     4,       & ! smooth down
                     4,       & ! smooth up
                     2,       & ! smooth bottom
                     .true.,  & ! zerophi
                     verbosity)

        ! Compute norm
        phi%data = phi%data - soln%data
        res = pnorm (phi, phi%valid, norm_type)

        ! Free memory
        call undefine_box_data (phi)
        call undefine_box_data (invdiags)
        call undefine_box_data (lphi)
        call undefine_box_data (soln)
        call undefine_box_data (bdx)
        call undefine_box_data (bdx_x)
        call undefine_box_data (bdx_y)
        call undefine_box_data (bdy)
        call undefine_box_data (bdy_x)
        call undefine_box_data (bdy_y)
        call undefine_bdry_data (bc)

    end function test_solver

end program test
