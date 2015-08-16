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
    integer, parameter          :: maxr = 7      ! Max refinement level
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


    ! Test 1: Non-uniform Dirichlet BCs
    errnorm = bogus_val
    do r = 1, maxr
        errnorm(r) = test_nonuniform_diri_bcs (geo(r))
    enddo
    call compute_conv_rate (rate, errnorm)
    print*, 'Test 1: Non-uniform Dirichlet BCs'
    print*, 'Error norm                rate'
    print*, errnorm(1)
    do r = 2, maxr
        print*, errnorm(r), rate(r-1)
    enddo
    print*


    ! Test 2: Non-uniform Neumann BCs
    errnorm = bogus_val
    do r = 1, maxr
        errnorm(r) = test_nonuniform_neum_bcs (geo(r))
    enddo
    call compute_conv_rate (rate, errnorm)
    print*, 'Test 2: Non-uniform Neumann BCs'
    print*, 'Error norm                rate'
    print*, errnorm(1)
    do r = 2, maxr
        print*, errnorm(r), rate(r-1)
    enddo


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
        call fill_elevation (y%data(:,jlo), x%data(:,x%valid%jlo), ilo, ihi, L)

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
    !  where Lx is the domain length
    ! --------------------------------------------------------------------------
    subroutine fill_elevation (el, x, ilo, ihi, Lx)
        implicit none

        integer, intent(in)                       :: ilo, ihi
        real(dp), intent(out), dimension(ilo:ihi) :: el
        real(dp), intent(in), dimension(ilo:ihi)  :: x
        real(dp), intent(in)                      :: Lx

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

            ! el(i) = fourth * sin(pi*xx/Lx)
            el(i) = zero

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
        integer        :: ngx, ngy, ilo, ihi, jlo, jhi
        real(dp)       :: scale

        ilo = dest%bx%ilo
        ihi = dest%bx%ihi
        jlo = dest%bx%jlo
        jhi = dest%bx%jhi
        ngx = dest%ngx
        ngy = dest%ngy

        ! We need xmu to surround the dest region.
        if (dest%offi .eq. BD_CELL) then
            call define_box_data (xmu, dest%valid, ngx, ngy, BD_NODE, BD_NODE)

            if (mu .eq. 1) then
                call fill_x (xmu)
            else
                call fill_y (xmu)
            endif

            if (nu .eq. 1) then
                scale = one / dest%valid%dx
                dest%data(ilo:ihi,:) = (xmu%data(ilo+1:ihi+1,:) - xmu%data(ilo:ihi,:)) * scale
            else
                scale = one / dest%valid%dy
                dest%data(:,jlo:jhi) = (xmu%data(:,jlo+1:jhi+1) - xmu%data(:,jlo:jhi)) * scale
            endif
        else
            call define_box_data (xmu, dest%valid, ngx+1, ngy+1, BD_CELL, BD_CELL)

            if (mu .eq. 1) then
                call fill_x (xmu)
            else
                call fill_y (xmu)
            endif

            if (nu .eq. 1) then
                scale = one / dest%valid%dx
                dest%data(ilo:ihi,:) = (xmu%data(ilo:ihi,:) - xmu%data(ilo-1:ihi-1,:)) * scale
            else
                scale = one / dest%valid%dy
                dest%data(:,jlo:jhi) = (xmu%data(:,jlo:jhi) - xmu%data(:,jlo-1:jhi-1)) * scale
            endif
        endif

        ! Free memory
        call undefine_box_data (xmu)

    end subroutine fill_dxdXi


    ! --------------------------------------------------------------------------
    ! Fills a box_data with J = det[Jacobian] = (dx/dXi) * (dy/dEta)
    ! --------------------------------------------------------------------------
    subroutine fill_J (dest)
        implicit none

        type(box_data), intent(inout) :: dest
        type(box_data)                :: dxdxi

        ! Initialize with dx/dXi
        call fill_dxdxi (dest, 1, 1)

        ! Multiply by dy/dEta
        call define_box_data (dxdxi, dest)
        call fill_dxdXi (dxdxi, 2, 2)
        dest%data = dest%data * dxdxi%data

        ! Free memory
        call undefine_box_data(dxdxi)

    end subroutine fill_J


    ! --------------------------------------------------------------------------
    ! Fills a box_data with the contravariant metric elements
    ! gup^{mu,nu} = Sum over rho [ dXi^{mu}/dx^{rho} * dXi^{nu}/dx^{rho} ]
    ! --------------------------------------------------------------------------
    subroutine fill_gup (dest, mu, nu)
        implicit none

        type(box_data), intent(inout) :: dest
        integer, intent(in)           :: mu, nu

        type(box_data)                :: dxdxi_mu_rho
        type(box_data)                :: dxdxi_nu_rho
        integer                       :: rho

        ! Create space for Jacobian matrix elements
        call define_box_data (dxdxi_mu_rho, dest)
        call define_box_data (dxdxi_nu_rho, dest)

        ! rho = 1 (x)...
        call fill_dxdxi(dxdxi_mu_rho, mu, 1)
        call fill_dxdxi(dxdxi_nu_rho, nu, 1)
        dest%data = dxdxi_mu_rho%data * dxdxi_nu_rho%data

        ! rho = 2 (y)...
        call fill_dxdxi(dxdxi_mu_rho, mu, 2)
        call fill_dxdxi(dxdxi_nu_rho, nu, 2)
        dest%data = dest%data + dxdxi_mu_rho%data * dxdxi_nu_rho%data

        ! Free memory
        call undefine_box_data(dxdxi_mu_rho)
        call undefine_box_data(dxdxi_nu_rho)

    end subroutine fill_gup


    ! --------------------------------------------------------------------------
    ! Fills a box_data with detJ * gup
    ! --------------------------------------------------------------------------
    subroutine fill_Jgup (dest, mu, nu)
        implicit none

        type(box_data), intent(inout) :: dest
        integer, intent(in)           :: mu, nu

        type(box_data)                :: J

        ! Initialize to gup
        call fill_gup (dest, mu, nu)

        ! Multiply by J
        call define_box_data (J, dest)
        dest%data = dest%data * J%data

        ! Free memory
        call undefine_box_data (J)

    end subroutine fill_Jgup

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
        real(dp)                   :: x, y
        type(box_data)             :: bdx, bdy, xlo, xhi, ylo, yhi

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

        ! Compute Cartesian locations of cell-centers
        call define_box_data (bdx, valid, 1, 1, BD_CELL, BD_CELL)
        call define_box_data (bdy, bdx)

        ! TODO: Compute x and y arrays at boundaries and use them to fill bdry_data
        ! call define_box_data (xlo, bdry_box(valid, BD_CELL, 1, SIDE_LO), 0, 0, BD_NODE, )

        call fill_x (bdx)
        call fill_y (bdy)

        ! Set up soln
        call define_box_data (soln, bdx)
        soln%data = sin((half + four*bdy%data/H)*pi*bdx%data/L)
        ! do j = jlo-1, jhi+1
        !     y = (j + half) * dy
        !     do i = ilo-1, ihi+1
        !         x = (i + half) * dx
        !         soln%data(i,j) = sin((half + four*y/H)*pi*x/L)
        !     enddo
        ! enddo

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

        ! Fill state's ghost cells.
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
        call undefine_box_data (bdy)
        call undefine_box_data (bdx)

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
        real(dp)                   :: x, y

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

        ! Set up field
        call define_box_data (soln, valid, 1, 1, BD_CELL, BD_CELL)
        soln%data = three

        call define_box_data (state, valid, 1, 1, BD_CELL, BD_CELL)
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
