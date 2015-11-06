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


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module arrayutils
    use precision
    implicit none

    ! --------------------------------------------------------------------------
    ! Array index extents and size info.
    ! --------------------------------------------------------------------------
    type box
        integer  :: ilo, ihi        ! Lower and upper array indices in x dir
        integer  :: jlo, jhi        ! Lower and upper array indices in y dir
        integer  :: nx, ny          ! Total number of cells/nodes in each dir
        real(dp) :: dx, dy          ! The grid spacing
    end type box


    ! --------------------------------------------------------------------------
    ! The empty box definition.
    ! --------------------------------------------------------------------------
    type(box), parameter :: empty_box = box(0, -1, 0, -1, 0, 0, zero, zero)


    ! --------------------------------------------------------------------------
    ! These are the BC types.
    ! --------------------------------------------------------------------------
    integer, parameter :: BCTYPE_UNDEFINED = -2
    integer, parameter :: BCTYPE_NONE = -1
    integer, parameter :: BCTYPE_NEUM = 0
    integer, parameter :: BCTYPE_DIRI = 1
    integer, parameter :: BCTYPE_PERIODIC = 2
    integer, parameter :: BCTYPE_CF = 3
    integer, parameter :: BCTYPE_EXTRAP0 = 4
    integer, parameter :: BCTYPE_EXTRAP1 = 5
    integer, parameter :: BCTYPE_EXTRAP2 = 6


    ! --------------------------------------------------------------------------
    ! These are the BC modes.
    ! BCMODE_UNIFORM = the entire boundary uses the same BC value. For example,
    ! a homogeneous BC would use a uniform value of 0.
    ! BCMODE_NONUNIFORM = the BC values vary over the bounding face. For
    ! example, a value of sin(k*y) along an x boundary face is nonuniform.
    ! BCMODE_LAH = Low-accuracy, homogenous BCs. These are used in deferred-
    ! correction schemes.
    ! --------------------------------------------------------------------------
    integer, parameter :: BCMODE_UNIFORM = 0
    integer, parameter :: BCMODE_NONUNIFORM = 1
    integer, parameter :: BCMODE_LAH = 2


    ! --------------------------------------------------------------------------
    ! Used to hold BC types and values at the domain boundary faces.
    ! --------------------------------------------------------------------------
    type bdry_data
        ! Valid (not ghost) box location and size info.
        type(box) :: valid

        ! The types of BCs that the data represent.
        integer :: type_xlo, type_xhi
        integer :: type_ylo, type_yhi

        ! The BC modes. See the BCMODE_ comments for details.
        integer :: mode_xlo, mode_xhi
        integer :: mode_ylo, mode_yhi

        ! The boundary data that surrounds the domain.
        real(dp), dimension(:), allocatable :: data_xlo, data_xhi
        real(dp), dimension(:), allocatable :: data_ylo, data_yhi
    end type bdry_data


    ! --------------------------------------------------------------------------
    ! When creating a box_data, these constants make the centering clear.
    ! --------------------------------------------------------------------------
    integer, parameter :: BD_NODE = 0
    integer, parameter :: BD_CELL = 1


    ! --------------------------------------------------------------------------
    ! These objects hold CC state data along with its metadata.
    ! --------------------------------------------------------------------------
    type box_data
        type(box) :: bx         ! The data box. Includes ghosts.
        type(box) :: valid      ! The valid region only. Does not include ghosts.
        integer   :: ngx, ngy   ! Nomber of ghosts per side in each dir

        integer   :: offi, offj ! Index offsets.
        real(dp)  :: offx, offy ! Data offsets. Ex: offx=half, offy=zero gives
                                ! a face-centered box_data in y.

        real(dp)  :: L, H       ! The dimensions of the valid region

        ! The ghost and valid data array.
        real(dp), dimension(:,:), allocatable :: data
    end type box_data

    ! --------------------------------------------------------------------------
    ! Side enumeration.
    ! --------------------------------------------------------------------------
    integer, parameter :: SIDE_LO = 0
    integer, parameter :: SIDE_HI = 1

    ! --------------------------------------------------------------------------
    ! Contains:
    !   J = the CC Jacobian determinant
    !   Jgup = the FC inverse metric tensor scaled by J
    ! --------------------------------------------------------------------------
    type geo_data
        type(box_data) :: J             ! Cell-centered.
        type(box_data) :: Jgup_xx       ! Face-centered in x-direction.
        type(box_data) :: Jgup_xy       ! Face-centered in x-direction.
        type(box_data) :: Jgup_yx       ! Face-centered in y-direction.
        type(box_data) :: Jgup_yy       ! Face-centered in y-direction.
        real(dp)       :: dx, dy
    end type geo_data


    ! --------------------------------------------------------------------------
    ! This makes defining a box_data a one line operation.
    ! --------------------------------------------------------------------------
    interface define_box_data
        module procedure define_box_data_itemize
        module procedure define_box_data_dup
    end interface

contains

    ! --------------------------------------------------------------------------
    ! This makes setting up a box object a one line operation.
    ! --------------------------------------------------------------------------
    pure subroutine define_box (bx, ilo, ihi, jlo, jhi, dx, dy)
        type(box), intent(out) :: bx
        integer, intent(in)    :: ilo, ihi
        integer, intent(in)    :: jlo, jhi
        real(dp), intent(in)   :: dx, dy

        bx%ilo = ilo
        bx%ihi = ihi

        bx%jlo = jlo
        bx%jhi = jhi

        bx%nx = ihi - ilo + 1
        bx%ny = jhi - jlo + 1

        bx%dx = dx
        bx%dy = dy
    end subroutine define_box


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    pure subroutine coarsen_box (bx, refx, refy)
        type(box), intent(inout) :: bx
        integer, intent(in)      :: refx, refy

        ! Coarsen in x dir
        bx%ilo = bx%ilo / refx
        bx%ihi = bx%ihi / refx
        bx%nx = bx%ihi - bx%ilo + 1
        bx%dx = bx%dx * refx

        ! Coarsen in y dir
        bx%jlo = bx%jlo / refy
        bx%jhi = bx%jhi / refy
        bx%ny = bx%jhi - bx%jlo + 1
        bx%dy = bx%dy * refy

    end subroutine coarsen_box


    ! --------------------------------------------------------------------------
    ! Currently not used.
    ! --------------------------------------------------------------------------
    pure subroutine grow_box (bx, growx, growy)
        type(box), intent(inout) :: bx
        integer, intent(in)      :: growx, growy

        bx%ilo = bx%ilo - growx
        bx%ihi = bx%ihi + growx
        bx%jlo = bx%jlo - growy
        bx%jhi = bx%jhi + growy
    end subroutine grow_box


    ! --------------------------------------------------------------------------
    ! Returns a face-centered box at the boundary of src.
    ! off = src's centering in dir.
    ! --------------------------------------------------------------------------
    function bdry_box (src, off, dir, side) result (bx)
        type(box), intent(in) :: src
        integer, intent(in)   :: off, dir, side
        type(box)             :: bx
        integer               :: ilo, ihi, jlo, jhi
        real(dp)              :: dx, dy

        ilo = src%ilo
        ihi = src%ihi
        jlo = src%jlo
        jhi = src%jhi
        dx = src%dx
        dy = src%dy

        if (dir .eq. 1) then
            if (side .eq. SIDE_LO) then
                call define_box (bx, ilo, ilo, jlo, jhi, dx, dy)

            else if (side .eq. SIDE_HI) then
                call define_box (bx, ihi+off, ihi+off, jlo, jhi, dx, dy)

            else
                print*, 'bdry_box: Bad side.'
                stop
            endif

        else if (dir .eq. 2) then
            if (side .eq. SIDE_LO) then
                call define_box (bx, ilo, ihi, jlo, jlo, dx, dy)

            else if (side .eq. SIDE_HI) then
                call define_box (bx, ilo, ihi, jhi+off, jhi+off, dx, dy)

            else
                print*, 'bdry_box: Bad side.'
                stop
            endif

        else
            print*, 'bdry_box: Bad dir.'
            stop
        endif
    end function bdry_box


    ! --------------------------------------------------------------------------
    ! This makes setting up a bdry_data object a one line operation.
    ! --------------------------------------------------------------------------
    subroutine define_bdry_data (bcd, valid, &
                                 type_xlo, type_xhi, type_ylo, type_yhi, &
                                 mode_xlo, mode_xhi, mode_ylo, mode_yhi)
        type(bdry_data), intent(out) :: bcd
        type(box), intent(in)        :: valid
        integer, intent(in)          :: type_xlo, type_xhi, type_ylo, type_yhi
        integer, intent(in)          :: mode_xlo, mode_xhi, mode_ylo, mode_yhi

        integer                      :: ierr

        bcd%valid = valid

        bcd%type_xlo = type_xlo
        bcd%type_xhi = type_xhi
        bcd%type_ylo = type_ylo
        bcd%type_yhi = type_yhi

        bcd%mode_xlo = mode_xlo
        bcd%mode_xhi = mode_xhi
        bcd%mode_ylo = mode_ylo
        bcd%mode_yhi = mode_yhi

        ! xlo
        select case (mode_xlo)
            case (BCMODE_UNIFORM)
                allocate (bcd%data_xlo (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_NONUNIFORM)
                allocate (bcd%data_xlo (bcd%valid%jlo : bcd%valid%jhi), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_LAH)
                allocate (bcd%data_xlo (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case default
                print*, 'define_bdry_data: Bad BCMODE'
                stop
        end select

        ! xhi
        select case (mode_xhi)
            case (BCMODE_UNIFORM)
                allocate (bcd%data_xhi (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_NONUNIFORM)
                allocate (bcd%data_xhi (bcd%valid%jlo : bcd%valid%jhi), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_LAH)
                allocate (bcd%data_xhi (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case default
                print*, 'define_bdry_data: Bad BCMODE'
                stop
        end select

        ! ylo
        select case (mode_ylo)
            case (BCMODE_UNIFORM)
                allocate (bcd%data_ylo (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_NONUNIFORM)
                allocate (bcd%data_ylo (bcd%valid%ilo : bcd%valid%ihi), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_LAH)
                allocate (bcd%data_ylo (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case default
                print*, 'define_bdry_data: Bad BCMODE'
                stop
        end select

        ! xhi
        select case (mode_yhi)
            case (BCMODE_UNIFORM)
                allocate (bcd%data_yhi (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_NONUNIFORM)
                allocate (bcd%data_yhi (bcd%valid%ilo : bcd%valid%ihi), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case (BCMODE_LAH)
                allocate (bcd%data_yhi (1:1), stat=ierr)
                if (ierr .ne. 0) then
                    print*, 'define_bdry_data: Out of memory'
                    stop
                endif
            case default
                print*, 'define_bdry_data: Bad BCMODE'
                stop
        end select
    end subroutine define_bdry_data


    ! --------------------------------------------------------------------------
    ! Frees memory used by a bdry_data object.
    ! --------------------------------------------------------------------------
    pure subroutine undefine_bdry_data (bcd)
        type(bdry_data), intent(inout) :: bcd

        bcd%valid = empty_box

        bcd%type_xlo = BCTYPE_UNDEFINED
        bcd%type_xhi = BCTYPE_UNDEFINED
        bcd%type_ylo = BCTYPE_UNDEFINED
        bcd%type_yhi = BCTYPE_UNDEFINED

        if (allocated(bcd%data_xlo)) deallocate(bcd%data_xlo)
        if (allocated(bcd%data_xhi)) deallocate(bcd%data_xhi)
        if (allocated(bcd%data_ylo)) deallocate(bcd%data_ylo)
        if (allocated(bcd%data_yhi)) deallocate(bcd%data_yhi)
    end subroutine undefine_bdry_data


    ! --------------------------------------------------------------------------
    ! This makes setting up a box_data object a one line operation.
    ! ioff = BD_NODE = 0 produces node-centered data in x,
    ! ioff = BD_CELL = 1 produces cell-centered data in x,
    ! etc...
    ! --------------------------------------------------------------------------
    subroutine define_box_data_itemize (bd, cc_valid, ngx, ngy, offi, offj)
        type(box_data), intent(out) :: bd
        type(box), intent(in)       :: cc_valid
        integer, intent(in)         :: ngx, ngy
        integer, intent(in)         :: offi, offj

        type(box) :: bx, valid
        integer   :: ierr
        real(dp)  :: offx, offy

        ! Define valid box with the correct staggering.
        valid = cc_valid

        select case (offi)
            case (BD_NODE)
                valid%ihi = valid%ihi + 1
                valid%nx = valid%nx + 1
                offx = zero
            case (BD_CELL)
                offx = half
            case default
                print*, 'define_box_data: Bad offi'
        end select

        select case (offj)
            case (BD_NODE)
                valid%jhi = valid%jhi + 1
                valid%ny = valid%ny + 1
                offy = zero
            case (BD_CELL)
                offy = half
            case default
                print*, 'define_box_data: Bad offj'
        end select

        ! Define the data box (bx) with the correct number of ghosts.
        bx = valid
        bx%ilo = bx%ilo - ngx
        bx%ihi = bx%ihi + ngx
        bx%jlo = bx%jlo - ngy
        bx%jhi = bx%jhi + ngy

        ! Define box_data object.
        bd%valid = valid
        bd%bx = bx
        bd%ngx = ngx
        bd%ngy = ngy

        bd%offx = offx
        bd%offy = offy
        bd%offi = offi
        bd%offj = offj

        bd%L = cc_valid%nx * cc_valid%dx
        bd%H = cc_valid%ny * cc_valid%dy

        if (allocated(bd%data)) deallocate (bd%data)
        allocate (bd%data (bx%ilo : bx%ihi, bx%jlo : bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'define_box_data: Out of memory'
            stop
        endif
    end subroutine define_box_data_itemize


    ! --------------------------------------------------------------------------
    ! Copies all metadata from src to dest. This defines dest over the exact
    ! same box as src, but doesn't copy src's data.
    ! --------------------------------------------------------------------------
    subroutine define_box_data_dup (dest, src)
        type(box_data), intent(inout) :: dest
        type(box_data), intent(in)    :: src
        type(box)                     :: cc_valid

        ! Construct the cell-centered box needed by the define function.
        cc_valid = src%valid
        if (src%offi .eq. BD_NODE) then
            cc_valid%ihi = cc_valid%ihi - 1
            cc_valid%nx = cc_valid%nx - 1
        endif
        if (src%offj .eq. BD_NODE) then
            cc_valid%jhi = cc_valid%jhi - 1
            cc_valid%ny = cc_valid%ny - 1
        endif

        call define_box_data_itemize (dest, cc_valid, src%ngx, src%ngy, src%offi, src%offj)
    end subroutine define_box_data_dup


    ! --------------------------------------------------------------------------
    ! Defines bd at the boundary faces of src.
    ! --------------------------------------------------------------------------
    subroutine define_box_data_bdry (bd, src, dir, side)
        type(box_data), intent(inout) :: bd
        type(box_data), intent(in)    :: src
        integer, intent(in)           :: dir, side
        integer                       :: ierr

        if (dir .eq. 1) then
            bd%valid = bdry_box (src%valid, src%offi, dir, side)
            bd%bx = bd%valid
            bd%ngx = 0
            bd%ngy = 0

            bd%offi = BD_NODE
            bd%offx = zero
            bd%offj = src%offj
            bd%offy = src%offy

            bd%L = src%L
            bd%H = src%H

        else if (dir .eq. 2) then
            bd%valid = bdry_box (src%valid, src%offj, dir, side)
            bd%bx = bd%valid
            bd%ngx = 0
            bd%ngy = 0

            bd%offi = src%offi
            bd%offx = src%offx
            bd%offj = BD_NODE
            bd%offy = zero

            bd%L = src%L
            bd%H = src%H
        else
            print*, 'define_box_data_bdry: Bad dir'
        endif

        if (allocated(bd%data)) deallocate (bd%data)
        allocate (bd%data (bd%bx%ilo : bd%bx%ihi, bd%bx%jlo : bd%bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'define_box_data: Out of memory'
            stop
        endif
    end subroutine define_box_data_bdry


    ! --------------------------------------------------------------------------
    ! Currently not used.
    ! --------------------------------------------------------------------------
    subroutine define_box_data_coarsened (crse, fine, refx, refy)
        type(box_data), intent(out) :: crse
        type(box_data), intent(in)  :: fine
        integer, intent(in)         :: refx, refy

        type(box)                   :: crse_valid

        crse_valid%ilo = fine%valid%ilo / refx
        crse_valid%ihi = fine%valid%ihi / refx
        crse_valid%jlo = fine%valid%jlo / refy
        crse_valid%jhi = fine%valid%jhi / refy

        call define_box_data_itemize (crse, crse_valid, fine%ngx, fine%ngy, fine%offi, fine%offj)

    end subroutine define_box_data_coarsened


    ! --------------------------------------------------------------------------
    ! Frees memory used by a box_data object.
    ! --------------------------------------------------------------------------
    pure subroutine undefine_box_data (bd)
        type(box_data), intent(inout) :: bd

        bd%valid = empty_box
        bd%bx = empty_box
        bd%ngx = 0
        bd%ngy = 0
        bd%offx = bogus_val
        bd%offy = bogus_val
        bd%offi = 0
        bd%offj = 0
        bd%L = bogus_val
        bd%H = bogus_val

        if (allocated(bd%data)) deallocate(bd%data)
    end subroutine undefine_box_data


    ! --------------------------------------------------------------------------
    ! Returns the opposite centering
    ! Only used in fill_dxdXi.
    ! --------------------------------------------------------------------------
    pure function stagger (min) result (mout)
        integer, intent(in) :: min
        integer             :: mout
        mout = 1 - min
    end function stagger


    ! --------------------------------------------------------------------------
    ! Returns true if all members of box1 and box2 are equal.
    ! Currently not used.
    ! --------------------------------------------------------------------------
    pure function compatible_boxes (bx1, bx2) result (is_same)
        type(box), intent(in) :: bx1, bx2
        logical :: is_same

        is_same = .true.

        if (bx1%ilo .ne. bx2%ilo) is_same = .false.
        if (bx1%ihi .ne. bx2%ihi) is_same = .false.
        if (bx1%jlo .ne. bx2%jlo) is_same = .false.
        if (bx1%jhi .ne. bx2%jhi) is_same = .false.

        if (bx1%dx .ne. bx2%dx) is_same = .false.
        if (bx1%dy .ne. bx2%dy) is_same = .false.
    end function compatible_boxes


    ! --------------------------------------------------------------------------
    ! Sets all valid data to val.
    ! Currently not used.
    ! --------------------------------------------------------------------------
    pure subroutine setval_valid (bd, val)
        type(box_data), intent(inout) :: bd
        real(dp), intent(in)          :: val

        bd%data (bd%valid%ilo : bd%valid%ihi, bd%valid%jlo : bd%valid%jhi) = val
    end subroutine setval_valid


    ! --------------------------------------------------------------------------
    ! Sets all ghosts to val.
    ! Currently not used.
    ! --------------------------------------------------------------------------
    pure subroutine setval_ghosts (bd, val)
        type(box_data), intent(inout) :: bd
        real(dp), intent(in)          :: val

        ! The ghosts along the x-normal boundaries.
        if (bd%ngx .ge. 1) then
            bd%data (bd%bx%ilo : bd%valid%ilo-1, bd%bx%jlo : bd%bx%jhi) = val
            bd%data (bd%valid%ihi+1 : bd%bx%ihi, bd%bx%jlo : bd%bx%jhi) = val
        endif

        ! The ghosts along the y-normal boundaries.
        if (bd%ngy .ge. 1) then
            bd%data (bd%bx%ilo : bd%bx%ihi, bd%bx%jlo : bd%valid%jlo-1) = val
            bd%data (bd%bx%ilo : bd%bx%ihi, bd%valid%jhi+1 : bd%bx%jhi) = val
        endif
    end subroutine setval_ghosts


    ! --------------------------------------------------------------------------
    ! Computes ip = (bd1, bd2)
    ! bd1 and bd2 are box_data objects that must be exactly the same size.
    ! --------------------------------------------------------------------------
    function inner_prod (bd1, bd2) result (ip)
        type(box_data), intent(in) :: bd1, bd2
        real(dp)                   :: ip

        integer                    :: ilo, ihi, jlo, jhi, j
        real(dp)                   :: vol

        ilo = max(bd1%valid%ilo, bd2%valid%ilo)
        ihi = min(bd1%valid%ihi, bd2%valid%ihi)
        jlo = max(bd1%valid%jlo, bd2%valid%jlo)
        jhi = min(bd1%valid%jhi, bd2%valid%jhi)

        vol = bd1%valid%dx * bd1%valid%dy

        ! Bounds checks
        if ((ilo .gt. ihi) .or. (jlo .gt. jhi)) then
            print*, 'ERROR: inner_prod_box_data: given box_data do not overlap.'
            stop
        endif

        ip = zero
        do j = jlo, jhi
            ip = ip + vol * dot_product (bd1%data(ilo:ihi,j), bd2%data(ilo:ihi,j))
        enddo
    end function inner_prod


    ! --------------------------------------------------------------------------
    ! Computes sum = Sum_{i,j} bd*J*dx*dy
    ! If opt_jscale is true (default), dx*dy will be scaled by J.
    ! bd and bx must be cell-centered.
    ! WARNING: This function does not check that bd contains bx.
    ! --------------------------------------------------------------------------
    function integrate2d (bd, bx, geo, opt_jscale) result (sum)
        type(box_data), intent(in) :: bd
        type(box), intent(in)      :: bx
        type(geo_data), intent(in) :: geo
        logical, optional          :: opt_jscale
        logical                    :: jscale
        real(dp)                   :: sum

        integer :: i, j
        real(dp) :: volScale

        if ((bd%offi .ne. BD_CELL) .or. (bd%offj .ne. BD_CELL)) then
            print*, 'pnorm: Only works with cell-centered data for now.'
            stop
        endif

        ! Are we scaling by J?
        jscale = .true.
        if (present(opt_jscale)) then
            jscale = opt_jscale
        endif

        volScale = geo%dx * geo%dy
        sum = zero

        if (jscale) then
            ! Integrate with J scaling.
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    sum = sum + volScale * geo%J%data(i,j) * bd%data(i,j)
                enddo
            enddo
        else
            ! Integrate without J scaling.
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    sum = sum + volScale * bd%data(i,j)
                enddo
            enddo
        endif
    end function integrate2d


    ! --------------------------------------------------------------------------
    ! Computes sum of outward normal fluxes.
    ! It is assumed that xflux and yflux are already scaled by J.
    ! bx should be cell-centered.
    ! xflux should be node-centered in x.
    ! yflux should be node-centered in y.
    ! WARNING: This function does not check that xflux or yflux overlap bx.
    ! --------------------------------------------------------------------------
    function integrate2d_bdry (xflux, yflux, bx) result (sum)
        type(box_data), intent(in) :: xflux, yflux
        type(box), intent(in)      :: bx
        real(dp)                   :: sum

        integer :: ilo, ihi, i
        integer :: jlo, jhi, j
        real(dp) :: dx, dy

        ilo = bx%ilo
        ihi = bx%ihi
        jlo = bx%jlo
        jhi = bx%jhi
        dx = bx%dx
        dy = bx%dy

        ! Compute integral.
        sum = zero
        do j = jlo, jhi
            sum = sum + dx * (xflux%data(ihi+1,j) - xflux%data(ilo,j))
        enddo
        do i = ilo, ihi
            sum = sum + dy * (yflux%data(i,jhi+1) - yflux%data(i,jlo))
        enddo

    end function integrate2d_bdry


    ! --------------------------------------------------------------------------
    ! Computes pnorm = |bd|_p = [Sum_{i,j} bd^p]^(1/p) / (nx*ny)
    ! WARNING: This function does not check that bd contains bx.
    ! --------------------------------------------------------------------------
    function pnorm (bd, bx, p) result (pn)
        type(box_data), intent(in) :: bd
        type(box), intent(in)      :: bx
        integer, intent(in)        :: p
        real(dp)                   :: pn

        integer :: i, j
        real(dp) :: volScale

        volScale = bd%bx%dx * bd%bx%dy
        pn = zero

        if (p .eq. 0) then
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    pn = max(pn, abs(bd%data(i,j)))
                enddo
            enddo

        else if (p .eq. 1) then
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    pn = pn + volScale * abs(bd%data(i,j))
                enddo
            enddo

        else if (p .eq. 2) then
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    pn = pn + volScale * bd%data(i,j)**2
                enddo
            enddo
            pn = sqrt(pn)

        else
            do j = bx%jlo, bx%jhi
                do i = bx%ilo, bx%ihi
                    pn = pn + volScale * abs(bd%data(i,j))**p
                enddo
            enddo
            pn = pn**(one/p)
        endif

    end function pnorm


    ! --------------------------------------------------------------------------
    ! Computes the pnorm over ghost cells only, skipping corner ghosts.
    ! --------------------------------------------------------------------------
    function gpnorm (bd, p) result (pn)
        type(box_data), intent(in) :: bd
        integer, intent(in)        :: p
        real(dp)                   :: pn
        type(box)                  :: bx

        integer :: i, j
        real(dp) :: dx, dy

        if ((bd%offi .ne. BD_CELL) .or. (bd%offj .ne. BD_CELL)) then
            print*, 'gpnorm: Only works with cell-centered data for now.'
            stop
        endif

        bx = bd%valid
        dx = bx%dx
        dy = bx%dy

        pn = zero

        if (p .eq. 0) then
            ! Lower x-boundary
            i = bx%jlo-1
            do j = bx%jlo, bx%jhi
                pn = max(pn, abs(bd%data(i,j)))
            enddo
            ! Upper x-boundary
            i = bx%jhi+1
            do j = bx%jlo, bx%jhi
                pn = max(pn, abs(bd%data(i,j)))
            enddo
            ! Lower y-boundary
            j = bx%ilo-1
            do i = bx%ilo, bx%ihi
                pn = max(pn, abs(bd%data(i,j)))
            enddo
            ! Upper y-boundary
            j = bx%ihi+1
            do i = bx%ilo, bx%ihi
                pn = max(pn, abs(bd%data(i,j)))
            enddo
        else
            ! Lower x-boundary
            i = bx%jlo-1
            do j = bx%jlo, bx%jhi
                pn = pn + dy * bd%data(i,j)**p
            enddo
            ! Upper x-boundary
            i = bx%jhi+1
            do j = bx%jlo, bx%jhi
                pn = pn + dy * bd%data(i,j)**p
            enddo
            ! Lower y-boundary
            j = bx%ilo-1
            do i = bx%ilo, bx%ihi
                pn = pn + dx * bd%data(i,j)**p
            enddo
            ! Upper y-boundary
            j = bx%ihi+1
            do i = bx%ilo, bx%ihi
                pn = pn + dx * bd%data(i,j)**p
            enddo
            pn = pn**(one/p)
        endif
    end function gpnorm


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine fill_ghosts_neum_dir_side (phi, bc, geo, homog, dir, side)
        type(box_data), intent(inout)       :: phi
        type(bdry_data), intent(in), target :: bc
        type(geo_data), intent(in)          :: geo
        logical, intent(in)                 :: homog
        integer, intent(in)                 :: dir
        integer, intent(in)                 :: side

        real(dp), dimension(:), pointer     :: bcd
        integer                             :: ilo, ihi, i
        integer                             :: jlo, jhi, j
        real(dp)                            :: dx, dy
        real(dp)                            :: cross, dpn, dpf
        integer                             :: bcmode, e

        ! Right now, we can only handle cell-centered data
        if ((phi%offi .ne. BD_CELL) .or. (phi%offj .ne. BD_CELL)) then
            print*, 'fill_ghosts: Only works with cell-centered data for now.'
            stop
        endif

        ! Right now, we can only handle 1 ghost layer.
        if ((phi%ngx .gt. 1) .or. (phi%ngy .gt. 1)) then
            print*, 'fill_ghosts: Can only handle 1 ghost layer max.'
            stop
        endif

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        dx = phi%valid%dx
        dy = phi%valid%dy

        if (dir .eq. 1) then
            if (side .lt. 1) then
                bcd => bc%data_xlo
                bcmode = bc%mode_xlo
                e = -1
                i = ilo
                dx = -dx
            else
                bcd => bc%data_xhi
                bcmode = bc%mode_xhi
                e = 1
                i = ihi
            endif
        else
            if (side .lt. 1) then
                bcd => bc%data_ylo
                bcmode = bc%mode_ylo
                e = -1
                j = jlo
                dy = -dy
            else
                bcd => bc%data_yhi
                bcmode = bc%mode_yhi
                e = 1
                j = jhi
            endif
        endif

        if (dir .eq. 1) then
            if (bcmode .eq. BCMODE_LAH) then
                phi%data(i+e,jlo:jhi) = phi%data(i,jlo:jhi)

            else if (homog) then
                j = jlo
                dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (-cross)*dx/geo%Jgup_xx%data(i,j)

                do j = jlo+1, jhi-1
                    dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                    dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                    phi%data(i+e,j) = phi%data(i,j) - (-cross)*dx/geo%Jgup_xx%data(i,j)
                enddo

                j = jhi
                dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (-cross)*dx/geo%Jgup_xx%data(i,j)

            else if (bcmode .eq. BCMODE_UNIFORM) then
                j = jlo
                dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (bcd(1)-cross)*dx/geo%Jgup_xx%data(i,j)

                do j = jlo, jhi
                    dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                    dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                    phi%data(i+e,j) = phi%data(i,j) - (bcd(1)-cross)*dx/geo%Jgup_xx%data(i,j)
                enddo

                j = jhi
                dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (bcd(1)-cross)*dx/geo%Jgup_xx%data(i,j)

            else
                j = jlo
                dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (bcd(j)-cross)*dx/geo%Jgup_xx%data(i,j)

                do j = jlo, jhi
                    dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                    dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                    phi%data(i+e,j) = phi%data(i,j) - (bcd(j)-cross)*dx/geo%Jgup_xx%data(i,j)
                enddo

                j = jhi
                dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                phi%data(i+e,j) = phi%data(i,j) - (bcd(j)-cross)*dx/geo%Jgup_xx%data(i,j)

            endif
        else
            if (bcmode .eq. BCMODE_LAH) then
                phi%data(ilo:ihi,j+e) = phi%data(ilo:ihi,j)

            else if (homog) then
                i = ilo
                dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (-cross)*dy/geo%Jgup_yy%data(i,j)

                do i = ilo, ihi
                    dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                    dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                    phi%data(i,j+e) = phi%data(i,j) - (-cross)*dy/geo%Jgup_yy%data(i,j)
                enddo

                i = ihi
                dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (-cross)*dy/geo%Jgup_yy%data(i,j)

            else if (bcmode .eq. BCMODE_UNIFORM) then
                i = ilo
                dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (bcd(1)-cross)*dy/geo%Jgup_yy%data(i,j)

                do i = ilo, ihi
                    dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                    dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                    phi%data(i,j+e) = phi%data(i,j) - (bcd(1)-cross)*dy/geo%Jgup_yy%data(i,j)
                enddo

                i = ihi
                dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (bcd(1)-cross)*dy/geo%Jgup_yy%data(i,j)

            else
                i = ilo
                dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (bcd(i)-cross)*dy/geo%Jgup_yy%data(i,j)

                do i = ilo, ihi
                    dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                    dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                    cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                    phi%data(i,j+e) = phi%data(i,j) - (bcd(i)-cross)*dy/geo%Jgup_yy%data(i,j)
                enddo

                i = ihi
                dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                phi%data(i,j+e) = phi%data(i,j) - (bcd(i)-cross)*dy/geo%Jgup_yy%data(i,j)
            endif
        endif
    end subroutine fill_ghosts_neum_dir_side


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine fill_ghosts_diri_dir_side (phi, bc, geo, homog, dir, side, order)
        type(box_data), intent(inout)       :: phi
        type(bdry_data), intent(in), target :: bc
        type(geo_data), intent(in)          :: geo
        logical, intent(in)                 :: homog
        integer, intent(in)                 :: dir
        integer, intent(in)                 :: side
        integer, intent(in)                 :: order

        real(dp), dimension(:), pointer     :: bcd
        integer                             :: ilo, ihi, i
        integer                             :: jlo, jhi, j
        integer                             :: bcmode, e

        real(dp), parameter                 :: eightthird     = eight / three
        real(dp), parameter                 :: sixteenfifth   = sixteen / five
        real(dp), parameter                 :: thirtyfifth    = one / (twenty+fifteen)
        real(dp), parameter                 :: onetwentyeight = 128.0E0_dp
        real(dp), parameter                 :: oneforty       = 140.0E0_dp
        real(dp), parameter                 :: seventy        = 70.0E0_dp
        real(dp), parameter                 :: twentyeight    = 28.0E0_dp

        ! Right now, we can only handle cell-centered data
        if ((phi%offi .ne. BD_CELL) .or. (phi%offj .ne. BD_CELL)) then
            print*, 'fill_ghosts: Only works with cell-centered data for now.'
            stop
        endif

        ! Right now, we can only handle 1 ghost layer.
        if ((phi%ngx .gt. 1) .or. (phi%ngy .gt. 1)) then
            print*, 'fill_ghosts: Can only handle 1 ghost layer max.'
            stop
        endif

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        if (dir .eq. 1) then
            if (side .lt. 1) then
                bcmode = bc%mode_xlo
                bcd => bc%data_xlo
                e = -1
                i = ilo
            else
                bcd => bc%data_xhi
                bcmode = bc%mode_xhi
                e = 1
                i = ihi
            endif
        else
            if (side .lt. 1) then
                bcd => bc%data_ylo
                bcmode = bc%mode_ylo
                e = -1
                j = jlo
            else
                bcd => bc%data_yhi
                bcmode = bc%mode_yhi
                e = 1
                j = jhi
            endif
        endif

        if (dir .eq. 1) then
            if (bcmode .eq. BCMODE_LAH) then
                phi%data(i+e, jlo:jhi) = -phi%data(i, jlo:jhi)
            else
                select case (order)
                    case (1)
                        if (homog) then
                            phi%data(i+e, jlo:jhi) = -phi%data(i, jlo:jhi)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(i+e, jlo:jhi) = two*bcd(1) - phi%data(i, jlo:jhi)
                        else
                            phi%data(i+e,jlo:jhi) = two*bcd(jlo:jhi) - phi%data(i,jlo:jhi)
                        endif
                    case (2)
                        if (homog) then
                            phi%data(i+e, jlo:jhi) =                                   &
                                                   -        two*phi%data(i  , jlo:jhi) &
                                                   +      third*phi%data(i-e, jlo:jhi)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(i+e, jlo:jhi) = eightthird*bcd(1)                 &
                                                   -        two*phi%data(i  , jlo:jhi) &
                                                   +      third*phi%data(i-e, jlo:jhi)
                        else
                            phi%data(i+e, jlo:jhi) = eightthird*bcd(jlo:jhi)           &
                                                   -        two*phi%data(i  , jlo:jhi) &
                                                   +      third*phi%data(i-e, jlo:jhi)
                        endif
                    case (3)
                        if (homog) then
                            phi%data(i+e, jlo:jhi) =                                       &
                                                   -        three*phi%data(i    , jlo:jhi) &
                                                   +              phi%data(i-  e, jlo:jhi) &
                                                   -        fifth*phi%data(i-2*e, jlo:jhi)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(i+e, jlo:jhi) = sixteenfifth*bcd(1)                   &
                                                   -        three*phi%data(i    , jlo:jhi) &
                                                   +              phi%data(i-  e, jlo:jhi) &
                                                   -        fifth*phi%data(i-2*e, jlo:jhi)
                        else
                            phi%data(i+e, jlo:jhi) = sixteenfifth*bcd(jlo:jhi)             &
                                                   -        three*phi%data(i    , jlo:jhi) &
                                                   +              phi%data(i-  e, jlo:jhi) &
                                                   -        fifth*phi%data(i-2*e, jlo:jhi)
                        endif
                    case (4)
                        if (homog) then
                            phi%data(i+e, jlo:jhi) = thirtyfifth*(-    oneforty*phi%data(i    , jlo:jhi) &
                                                                  +     seventy*phi%data(i-  e, jlo:jhi) &
                                                                  - twentyeight*phi%data(i-2*e, jlo:jhi) &
                                                                  +        five*phi%data(i-3*e, jlo:jhi))
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(i+e, jlo:jhi) = thirtyfifth*(onetwentyeight*bcd(1)                   &
                                                                  -     oneforty*phi%data(i    , jlo:jhi) &
                                                                  +      seventy*phi%data(i-  e, jlo:jhi) &
                                                                  -  twentyeight*phi%data(i-2*e, jlo:jhi) &
                                                                  +         five*phi%data(i-3*e, jlo:jhi))
                        else
                            phi%data(i+e, jlo:jhi) = thirtyfifth*(onetwentyeight*bcd(jlo:jhi)             &
                                                                  -     oneforty*phi%data(i    , jlo:jhi) &
                                                                  +      seventy*phi%data(i-  e, jlo:jhi) &
                                                                  -  twentyeight*phi%data(i-2*e, jlo:jhi) &
                                                                  +         five*phi%data(i-3*e, jlo:jhi))
                        endif
                    case default
                        print*, 'fill_ghosts: Bad order'
                end select
            endif
        else
            if (bcmode .eq. BCMODE_LAH) then
                phi%data(ilo:ihi, j+e) = -phi%data(ilo:ihi, j)
            else
                select case (order)
                    case (1)
                        if (homog) then
                            phi%data(ilo:ihi, j+e) = -phi%data(ilo:ihi, j)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(ilo:ihi, j+e) = two*bcd(1) - phi%data(ilo:ihi, j)
                        else
                            phi%data(ilo:ihi, j+e) = two*bcd(ilo:ihi) - phi%data(ilo:ihi, j)
                        endif
                    case (2)
                        if (homog) then
                            phi%data(ilo:ihi, j+e) =                                   &
                                                   -        two*phi%data(ilo:ihi, j  ) &
                                                   +      third*phi%data(ilo:ihi, j-e)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(ilo:ihi, j+e) = eightthird*bcd(1)                 &
                                                   -        two*phi%data(ilo:ihi, j  ) &
                                                   +      third*phi%data(ilo:ihi, j-e)
                        else
                            phi%data(ilo:ihi, j+e) = eightthird*bcd(jlo:jhi)           &
                                                   -        two*phi%data(ilo:ihi, j  ) &
                                                   +      third*phi%data(ilo:ihi, j-e)
                        endif
                    case (3)
                        if (homog) then
                            phi%data(ilo:ihi, j+e) =                                       &
                                                   -        three*phi%data(ilo:ihi, j    ) &
                                                   +              phi%data(ilo:ihi, j-  e) &
                                                   -        fifth*phi%data(ilo:ihi, j-2*e)
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(ilo:ihi, j+e) = sixteenfifth*bcd(1)                   &
                                                   -        three*phi%data(ilo:ihi, j    ) &
                                                   +              phi%data(ilo:ihi, j-  e) &
                                                   -        fifth*phi%data(ilo:ihi, j-2*e)
                        else
                            phi%data(ilo:ihi, j+e) = sixteenfifth*bcd(jlo:jhi)             &
                                                   -        three*phi%data(ilo:ihi, j    ) &
                                                   +              phi%data(ilo:ihi, j-  e) &
                                                   -        fifth*phi%data(ilo:ihi, j-2*e)
                        endif
                    case (4)
                        if (homog) then
                            phi%data(ilo:ihi, j+e) = thirtyfifth*(-    oneforty*phi%data(ilo:ihi, j    ) &
                                                                  +     seventy*phi%data(ilo:ihi, j-  e) &
                                                                  - twentyeight*phi%data(ilo:ihi, j-2*e) &
                                                                  +        five*phi%data(ilo:ihi, j-3*e))
                        else if (bcmode .eq. BCMODE_UNIFORM) then
                            phi%data(ilo:ihi, j+e) = thirtyfifth*(onetwentyeight*bcd(1)                   &
                                                                  -     oneforty*phi%data(ilo:ihi, j    ) &
                                                                  +      seventy*phi%data(ilo:ihi, j-  e) &
                                                                  -  twentyeight*phi%data(ilo:ihi, j-2*e) &
                                                                  +         five*phi%data(ilo:ihi, j-3*e))
                        else
                            phi%data(ilo:ihi, j+e) = thirtyfifth*(onetwentyeight*bcd(jlo:jhi)             &
                                                                  -     oneforty*phi%data(ilo:ihi, j    ) &
                                                                  +      seventy*phi%data(ilo:ihi, j-  e) &
                                                                  -  twentyeight*phi%data(ilo:ihi, j-2*e) &
                                                                  +         five*phi%data(ilo:ihi, j-3*e))
                        endif
                    case default
                        print*, 'fill_ghosts: Bad order'
                end select
            endif
        endif

        nullify(bcd)

    end subroutine fill_ghosts_diri_dir_side


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine fill_ghosts (phi, bc, geo, homog, do_neum_opt)
        type(box_data), intent(inout) :: phi
        type(bdry_data), intent(in)   :: bc
        type(geo_data), intent(in)    :: geo
        logical, intent(in)           :: homog
        logical, intent(in), optional :: do_neum_opt

        integer, parameter            :: diri_order_xlo = 3
        integer, parameter            :: diri_order_xhi = 3
        integer, parameter            :: diri_order_ylo = 3
        integer, parameter            :: diri_order_yhi = 3

        integer  :: ilo, ihi, jlo, jhi
        real(dp) :: dx, dy
        logical  :: do_neum

        real(dp) :: refx, c1x, c2x
        real(dp) :: refy, c1y, c2y

        refx = four                 ! TODO: Make this an input parameter
        c1x = 2*(refx-1)/(refx+1)
        c1x =  -(refx-1)/(refx+3)

        refy = four                 ! TODO: Make this an input parameter
        c1y = 2*(refy-1)/(refy+1)
        c1y =  -(refy-1)/(refy+3)

        if ((phi%offi .ne. BD_CELL) .or. (phi%offj .ne. BD_CELL)) then
            print*, 'fill_ghosts: Only works with cell-centered data for now.'
            stop
        endif

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        dx = phi%valid%dx
        dy = phi%valid%dy

        ! By default, we apply all BCs.
        if (.not. present(do_neum_opt)) then
            do_neum = .true.
        else
            do_neum = do_neum_opt
        endif

        ! Right now, we can only handle 1 ghost layer.
        if ((phi%ngx .gt. 1) .or. (phi%ngy .gt. 1)) then
            print*, 'fill_ghosts: Can only handle 1 ghost layer max.'
            stop
        endif

        ! Lower x ghosts
        select case (bc%type_xlo)
            case (BCTYPE_NEUM)
                if (do_neum) then
                    call fill_ghosts_neum_dir_side (phi, bc, geo, homog, 1, -1)
                endif

            case (BCTYPE_DIRI)
                call fill_ghosts_diri_dir_side (phi, bc, geo, homog, 1, -1, diri_order_xlo)

            case (BCTYPE_PERIODIC)
                phi%data(ilo-1, jlo:jhi) = phi%data(ihi, jlo:jhi)

            case (BCTYPE_CF)
                phi%data(ilo-1, jlo:jhi) = c1x*phi%data(ilo, jlo:jhi) + c2x*phi%data(ilo+1, jlo:jhi)

            case (BCTYPE_EXTRAP0)
                phi%data(ilo-1,:) = phi%data(ilo,:)

            case (BCTYPE_EXTRAP1)
                phi%data(ilo-1,:) = two*phi%data(ilo,:) - phi%data(ilo+1,:)

            case (BCTYPE_EXTRAP2)
                phi%data(ilo-1,:) = three*(phi%data(ilo,:) - phi%data(ilo+1,:)) + phi%data(ilo+2,:)

            case (BCTYPE_NONE)
                ! Do nothing

            case default
                print*, 'fill_ghosts: invalid BCTYPE_'
                stop
        end select

        ! Upper x ghosts
        select case (bc%type_xhi)
            case (BCTYPE_NEUM)
                if (do_neum) then
                    call fill_ghosts_neum_dir_side (phi, bc, geo, homog, 1, 1)
                endif

            case (BCTYPE_DIRI)
                call fill_ghosts_diri_dir_side (phi, bc, geo, homog, 1, 1, diri_order_xhi)

            case (BCTYPE_PERIODIC)
                phi%data(ihi+1, jlo:jhi) = phi%data(ilo, jlo:jhi)

            case (BCTYPE_CF)
                phi%data(ihi+1, jlo:jhi) = c1x*phi%data(ihi, jlo:jhi) + c2x*phi%data(ihi-1, jlo:jhi)

            case (BCTYPE_EXTRAP0)
                phi%data(ihi+1,:) = phi%data(ihi,:)

            case (BCTYPE_EXTRAP1)
                phi%data(ihi+1,:) = two*phi%data(ihi,:) - phi%data(ihi-1,:)

            case (BCTYPE_EXTRAP2)
                phi%data(ihi+1,:) = three*(phi%data(ihi,:) - phi%data(ihi-1,:)) + phi%data(ihi-2,:)

            case (BCTYPE_NONE)
                ! Do nothing

            case default
                print*, 'fill_ghosts: invalid BCTYPE_'
                stop
        end select

        ! Lower y ghosts
        select case (bc%type_ylo)
            case (BCTYPE_NEUM)
                if (do_neum) then
                    call fill_ghosts_neum_dir_side (phi, bc, geo, homog, 2, -1)
                endif

            case (BCTYPE_DIRI)
                call fill_ghosts_diri_dir_side (phi, bc, geo, homog, 2, -1, diri_order_ylo)

            case (BCTYPE_PERIODIC)
                phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jhi)

            case (BCTYPE_CF)
                phi%data(ilo:ihi, jlo-1) = c1y*phi%data(ilo:ihi, jlo) + c2y*phi%data(ilo:ihi, jlo+1)

            case (BCTYPE_EXTRAP0)
                phi%data(:,jlo-1) = phi%data(:,jlo)

            case (BCTYPE_EXTRAP1)
                phi%data(:,jlo-1) = two*phi%data(:,jlo) - phi%data(:,jlo+1)

            case (BCTYPE_EXTRAP2)
                phi%data(:,jlo-1) = three*(phi%data(:,jlo) - phi%data(:,jlo+1)) + phi%data(:,jlo+2)

            case (BCTYPE_NONE)
                ! Do nothing

            case default
                print*, 'fill_ghosts: invalid BCTYPE_'
                stop
        end select

        ! Upper y ghosts
        select case (bc%type_yhi)
            case (BCTYPE_NEUM)
                if (do_neum) then
                    call fill_ghosts_neum_dir_side (phi, bc, geo, homog, 2, 1)
                endif

            case (BCTYPE_DIRI)
                call fill_ghosts_diri_dir_side (phi, bc, geo, homog, 2, 1, diri_order_yhi)

            case (BCTYPE_PERIODIC)
                phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jlo)

            case (BCTYPE_CF)
                phi%data(ilo:ihi, jhi+1) = c1y*phi%data(ilo:ihi, jhi) + c2y*phi%data(ilo:ihi, jhi-1)

            case (BCTYPE_EXTRAP0)
                phi%data(:,jhi+1) = phi%data(:,jhi)

            case (BCTYPE_EXTRAP1)
                phi%data(:,jhi+1) = two*phi%data(:,jhi) - phi%data(:,jhi-1)

            case (BCTYPE_EXTRAP2)
                phi%data(:,jhi+1) = three*(phi%data(:,jhi) - phi%data(:,jhi-1)) + phi%data(:,jhi-2)

            case (BCTYPE_NONE)
                ! Do nothing

            case default
                print*, 'fill_ghosts: invalid BCTYPE_'
                stop
        end select

        ! Fill corner ghosts.
        if ((bc%type_xlo .ne. BCTYPE_PERIODIC) .and. (bc%type_ylo .ne. BCTYPE_PERIODIC)) then
            ! Put nans in the corner ghosts
            phi%data(ilo-1,jlo-1) = bogus_val
            phi%data(ilo-1,jhi+1) = bogus_val
            phi%data(ihi+1,jlo-1) = bogus_val
            phi%data(ihi+1,jhi+1) = bogus_val

        else if ((bc%type_xlo .eq. BCTYPE_PERIODIC) .and. (bc%type_ylo .eq. BCTYPE_PERIODIC)) then
            ! Both directions periodic
            phi%data(ilo-1,jlo-1) = phi%data(ihi,jhi)
            phi%data(ilo-1,jhi+1) = phi%data(ihi,jlo)
            phi%data(ihi+1,jlo-1) = phi%data(ilo,jhi)
            phi%data(ihi+1,jhi+1) = phi%data(ilo,jlo)

        else if (bc%type_xlo .eq. BCTYPE_PERIODIC) then
            ! Only x-dir is periodic
            phi%data(ilo-1,jlo-1) = phi%data(ihi,jlo-1)
            phi%data(ilo-1,jhi+1) = phi%data(ihi,jhi+1)
            phi%data(ihi+1,jlo-1) = phi%data(ilo,jlo-1)
            phi%data(ihi+1,jhi+1) = phi%data(ilo,jhi+1)

        else if (bc%type_ylo .eq. BCTYPE_PERIODIC) then
            ! Only y-dir is periodic
            phi%data(ilo-1,jlo-1) = phi%data(ilo-1,jhi)
            phi%data(ihi+1,jlo-1) = phi%data(ihi+1,jhi)
            phi%data(ilo-1,jhi+1) = phi%data(ilo-1,jlo)
            phi%data(ihi+1,jhi+1) = phi%data(ihi+1,jlo)
        endif

    end subroutine fill_ghosts


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine extrapolate_ghosts (phi, order)
        type(box_data), intent(inout) :: phi
        integer, intent(in)           :: order

        integer  :: ilo, ihi, jlo, jhi

        ! Only works with CC data.
        if ((phi%offi .ne. BD_CELL) .or. (phi%offj .ne. BD_CELL)) then
            print*, 'extrapolate_ghosts: Only works with cell-centered data.'
            stop
        endif

        ! Right now, we can only handle 1 ghost layer at most.
        if ((phi%ngx .gt. 1) .or. (phi%ngy .gt. 1)) then
            print*, 'fill_ghosts: Can only handle 1 ghost layer max.'
            stop
        endif

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        select case (order)
            case (0)
                if (phi%ngx .gt. 0) then
                    phi%data(ilo-1,:) = phi%data(ilo,:)
                    phi%data(ihi+1,:) = phi%data(ihi,:)
                endif

                if (phi%ngy .gt. 0) then
                    phi%data(:,jlo-1) = phi%data(:,jlo)
                    phi%data(:,jhi+1) = phi%data(:,jhi)
                endif
            case (1)
                if (phi%ngx .gt. 0) then
                    phi%data(ilo-1,:) = two*phi%data(ilo,:) - phi%data(ilo+1,:)
                    phi%data(ihi+1,:) = two*phi%data(ihi,:) - phi%data(ihi-1,:)
                endif

                if (phi%ngy .gt. 0) then
                    phi%data(:,jlo-1) = two*phi%data(:,jlo) - phi%data(:,jlo+1)
                    phi%data(:,jhi+1) = two*phi%data(:,jhi) - phi%data(:,jhi-1)
                endif
            case (2)
                if (phi%ngx .gt. 0) then
                    phi%data(ilo-1,:) = three*(phi%data(ilo,:) - phi%data(ilo+1,:)) + phi%data(ilo+2,:)
                    phi%data(ihi+1,:) = three*(phi%data(ihi,:) - phi%data(ihi-1,:)) + phi%data(ihi-2,:)
                endif

                if (phi%ngy .gt. 0) then
                    phi%data(:,jlo-1) = three*(phi%data(:,jlo) - phi%data(:,jlo+1)) + phi%data(:,jlo+2)
                    phi%data(:,jhi+1) = three*(phi%data(:,jhi) - phi%data(:,jhi-1)) + phi%data(:,jhi-2)
                endif
            case (3)
                if (phi%ngx .gt. 0) then
                    phi%data(ilo-1,:) = four*(phi%data(ilo,:) + phi%data(ilo+2,:)) - six*phi%data(ilo+1,:) - phi%data(ilo+3,:)
                    phi%data(ihi+1,:) = four*(phi%data(ihi,:) + phi%data(ihi-2,:)) - six*phi%data(ihi-1,:) - phi%data(ihi-3,:)
                endif

                if (phi%ngy .gt. 0) then
                    phi%data(:,jlo-1) = four*(phi%data(:,jlo) + phi%data(:,jlo+2)) - six*phi%data(:,jlo+1) - phi%data(:,jlo+3)
                    phi%data(:,jhi+1) = four*(phi%data(:,jhi) + phi%data(:,jhi-2)) - six*phi%data(:,jhi-1) - phi%data(:,jhi-3)
                endif
            case default
                print*, 'extrapolate_ghosts: order = ', order, ' is not supported.'
                stop
        end select

    end subroutine extrapolate_ghosts


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine fill_boundary_fluxes (xflux, yflux, bc, homog)
        type(box_data), intent(inout) :: xflux, yflux
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        if (homog) then
            if (bc%type_xlo .eq. BCTYPE_NEUM) xflux%data(xflux%valid%ilo,:) = zero
            if (bc%type_xhi .eq. BCTYPE_NEUM) xflux%data(xflux%valid%ihi,:) = zero
            if (bc%type_ylo .eq. BCTYPE_NEUM) yflux%data(:,yflux%valid%jlo) = zero
            if (bc%type_yhi .eq. BCTYPE_NEUM) yflux%data(:,yflux%valid%jhi) = zero
        else
            if (bc%type_xlo .eq. BCTYPE_NEUM) then
                if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                    xflux%data(xflux%valid%ilo,:) = bc%data_xlo(1)
                else if (bc%mode_xlo .eq. BCMODE_NONUNIFORM) then
                    xflux%data(xflux%valid%ilo,:) = bc%data_xlo
                else
                    xflux%data(xflux%valid%ilo,:) = zero
                endif
            endif

            if (bc%type_xhi .eq. BCTYPE_NEUM) then
                if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                    xflux%data(xflux%valid%ihi,:) = bc%data_xhi(1)
                else if (bc%mode_xhi .eq. BCMODE_NONUNIFORM) then
                    xflux%data(xflux%valid%ihi,:) = bc%data_xhi
                else
                    xflux%data(xflux%valid%ihi,:) = zero
                endif
            endif

            if (bc%type_ylo .eq. BCTYPE_NEUM) then
                if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                    yflux%data(:,yflux%valid%jlo) = bc%data_ylo(1)
                else if (bc%mode_ylo .eq. BCMODE_NONUNIFORM) then
                    yflux%data(:,yflux%valid%jlo) = bc%data_ylo
                else
                    yflux%data(:,yflux%valid%jlo) = zero
                endif
            endif

            if (bc%type_yhi .eq. BCTYPE_NEUM) then
                if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                    yflux%data(:,yflux%valid%jhi) = bc%data_yhi(1)
                else if (bc%mode_yhi .eq. BCMODE_NONUNIFORM) then
                    yflux%data(:,yflux%valid%jhi) = bc%data_yhi
                else
                    yflux%data(:,yflux%valid%jhi) = zero
                endif
            endif
        endif
    end subroutine fill_boundary_fluxes

end module ArrayUtils



! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module Poisson2D
    use ArrayUtils
    implicit none

    save

contains

    ! --------------------------------------------------------------------------
    ! Computes 1.0 / the Laplacian's diagonal matrix elements.
    ! invdiags must be prepared (allocated and box set) prior to call.
    ! NOTE: This assumes the Laplacian is not scaling by 1/J.
    ! --------------------------------------------------------------------------
    subroutine compute_invdiags (invdiags, geo, bc, homog)
        type(box_data), intent(inout) :: invdiags
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        integer                       :: ilo, ihi, inner_ilo, inner_ihi, i
        integer                       :: jlo, jhi, inner_jlo, inner_jhi, j
        real(dp)                      :: ee,en,ew,es
        real(dp)                      :: ene,enw,esw,ese
        real(dp)                      :: xfe,xfw,yfn,yfs
        real(dp)                      :: idx, idy

        real(dp), parameter           :: p = one
        real(dp), parameter           :: pe = zero
        real(dp), parameter           :: pn = zero
        real(dp), parameter           :: pw = zero
        real(dp), parameter           :: ps = zero

        if ((bc%mode_xlo .eq. BCMODE_LAH) .or. (bc%mode_xhi .eq. BCMODE_LAH) .or. &
            (bc%mode_ylo .eq. BCMODE_LAH) .or. (bc%mode_yhi .eq. BCMODE_LAH)) then
            print*, 'Poisson2d::compute_invdiags received BCMODE_LAH.'
            stop
        endif

        ilo = invdiags%valid%ilo
        ihi = invdiags%valid%ihi
        jlo = invdiags%valid%jlo
        jhi = invdiags%valid%jhi

        if (bc%type_xlo .eq. BCTYPE_PERIODIC) then
            inner_ilo = ilo
            inner_ihi = ihi
        else
            inner_ilo = ilo+1
            inner_ihi = ihi-1
        endif

        if (bc%type_ylo .eq. BCTYPE_PERIODIC) then
            inner_jlo = jlo
            inner_jhi = jhi
        else
            inner_jlo = jlo+1
            inner_jhi = jhi-1
        endif

        idx = one / geo%dx
        idy = one / geo%dy

        ! Lower x boundary (avoid west), lower y boundary (avoid south)
        if (bc%type_ylo .ne. BCTYPE_PERIODIC) then
            j = jlo
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                ee = zero
                en = zero
                ew = three*(p-ee)+zero !two*p-ee
                es = three*(p-en)+zero !two*p-en

                ene = zero
                enw = three*(en-ene)+zero !two*en-ene
                ese = three*(ee-ene)+zero !two*ee-ene
                esw = half*((three*(es-ese)+zero)+(three*(ew-enw)+zero)) !half*((two*es-ese)+(two*ew-enw))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif
                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x, lower y boundary (avoid south)
            do i = inner_ilo, inner_ihi
                ee = zero
                ew = zero
                en = zero
                es = three*(p-en)+zero

                ene = zero
                enw = zero
                ese = three*(ee-ene)+zero !two*ee-ene
                esw = three*(ew-enw)+zero !two*ew-enw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), lower y boundary (avoid south)
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                ew = zero
                en = zero
                ee = three*(p-ew)+zero !two*p-ew
                es = three*(p-en)+zero !two*p-en

                enw = zero
                ene = three*(pn-enw)+zero !two*pn-enw
                esw = three*(pw-enw)+zero !two*pw-enw
                ese = half*((three*(es-esw)+zero)+(three*(ee-ene)+zero)) !half*((two*es-esw)+(two*ee-ene))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif
                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        endif ! not periodic in y

        ! Interior to y...
        do j = inner_jlo, inner_jhi
            ! Lower x boundary (avoid west), interior to y
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                ee = zero
                en = zero
                es = zero
                ew = three*(p-ee)+zero !two*p-ee

                ene = zero
                ese = zero
                enw = three*(en-ene)+zero !two*en-ene
                esw = three*(es-ese)+zero !two*es-ese

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x and y
            do i = inner_ilo, inner_ihi
                ee = zero
                ew = zero
                en = zero
                es = zero

                ene = zero
                enw = zero
                ese = zero
                esw = zero

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), interior to y
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                ew = zero
                en = zero
                es = zero
                ee = three*(p-ew)+zero !two*p-ew

                enw = zero
                esw = zero
                ene = three*(en-enw)+zero !two*en-enw
                ese = three*(es-esw)+zero !two*es-esw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        enddo !j


        ! Upper y boundary...
        if (bc%type_yhi .ne. BCTYPE_PERIODIC) then
            j = jhi

            ! Lower x boundary (avoid west), upper y boundary (avoid north)
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                ee = zero
                es = zero
                ew = three*(p-ee)+zero !two*p-ee
                en = three*(p-es)+zero !two*p-es

                ese = zero
                ene = three*(ee-ese)+zero !two*ee-ese
                esw = three*(es-ese)+zero !two*es-ese
                enw = half*((three*(ew-esw)+zero)+(three*(en-ene)+zero)) !half*((two*ew-esw)+(two*en-ene))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif
                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x, upper y boundary (avoid north)
            do i = inner_ilo, inner_ihi
                ee = zero
                ew = zero
                es = zero
                en = three*(p-es)+zero !two*p-es

                ese = zero
                esw = zero
                ene = three*(ee-ese)+zero !two*ee-ese
                enw = three*(ew-esw)+zero !two*ew-esw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), upper y boundary (avoid north)
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                ew = zero
                es = zero
                ee = three*(p-ew)+zero !two*p-ew
                en = three*(p-es)+zero !two*p-es

                esw = zero
                enw = three*(ew-esw)+zero !two*ew-esw
                ese = three*(es-esw)+zero !two*es-esw
                ene = half*((three*(en-enw)+zero)+(three*(ee-ese)+zero)) !half*((two*en-enw)+(two*ee-ese))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif
                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                invdiags%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        endif ! not periodic in y

        ! Invert the diag values
        invdiags%data = one / invdiags%data

    end subroutine compute_invdiags


    ! ------------------------------------------------------------------------------
    ! Computes partial derivatives.
    ! phi is expected to be cell-centered with ghosts filled if needed.
    ! pd is expected to be face-centered. No ghosts will be filled.
    !
    ! This function uses the EXACT same stencils as compute_laplacian if
    ! simple_bdry_stencil = .false., otherwise the stencils differ by O(dx^2)
    ! at the boundaries of transverse derivatives.
    !
    ! NOTE: This function only uses ghosts if we are differentiating in the
    !       nodal direction.
    ! NOTE: pd must be nodal in exactly one direction.
    ! ------------------------------------------------------------------------------
    subroutine compute_pd (pd, phi, dir)
        type(box_data), intent(inout) :: pd
        type(box_data), intent(in)    :: phi
        integer, intent(in)           :: dir

        real(dp)                      :: scale
        integer                       :: ilo, ihi, jlo, jhi
        integer                       :: nodedir

        integer                       :: i, j
        real(dp)                      :: p,pe,pn,pw,ps
        real(dp)                      :: ee,en,ew,es
        real(dp)                      :: ene,enw,esw,ese

        logical, parameter            :: simple_bdry_stencil = .false.

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        ! Find the nodal direction
        if ((pd%offi .ne. 0) .and. (pd%offj .ne. 0)) then
            print*, 'compute_pd cannot handle totally cell-centered data.'
            stop
        endif

        if ((pd%offi .eq. 0) .and. (pd%offj .eq. 0)) then
            print*, 'compute_pd cannot handle totally nodal-centered data.'
            stop
        endif

        if (pd%offi .eq. 0) then
            nodedir = 1
        else
            nodedir = 2
        endif


        ! Compute xflux...

        if (dir .eq. 1) then
            if (nodedir .eq. dir) then
                scale = one / phi%valid%dx
                pd%data(ilo:ihi+1,jlo:jhi) = scale * (  phi%data(ilo:ihi+1,jlo:jhi) &
                                                      - phi%data(ilo-1:ihi,jlo:jhi))
            else
                scale = fourth / phi%valid%dx

                ! Away from y boundaries...
                ! Interior
                pd%data(ilo+1:ihi-1,jlo+1:jhi) = scale * (  phi%data(ilo+2:ihi,jlo+1:jhi) &
                                                          - phi%data(ilo:ihi-2,jlo+1:jhi) &
                                                          + phi%data(ilo+2:ihi,jlo:jhi-1) &
                                                          - phi%data(ilo:ihi-2,jlo:jhi-1)  )

                ! Lower x boundary
                pd%data(ilo,jlo+1:jhi) = -scale * (  three*phi%data(ilo,jlo+1:jhi)   &
                                                   -  four*phi%data(ilo+1,jlo+1:jhi) &
                                                   +       phi%data(ilo+2,jlo+1:jhi) &
                                                   + three*phi%data(ilo,jlo:jhi-1)   &
                                                   -  four*phi%data(ilo+1,jlo:jhi-1) &
                                                   +       phi%data(ilo+2,jlo:jhi-1)  )

                ! Upper x boundary
                pd%data(ihi,jlo+1:jhi) = scale * (  three*phi%data(ihi,jlo+1:jhi)   &
                                                  -  four*phi%data(ihi-1,jlo+1:jhi) &
                                                  +       phi%data(ihi-2,jlo+1:jhi) &
                                                  + three*phi%data(ihi,jlo:jhi-1)   &
                                                  -  four*phi%data(ihi-1,jlo:jhi-1) &
                                                  +       phi%data(ihi-2,jlo:jhi-1)  )

                ! At y boundaries...
                if (simple_bdry_stencil) then
                    ! Lower y boundary
                    pd%data(ilo:ihi,jlo) = three*(  pd%data(ilo:ihi,jlo+1)    &
                                                  - pd%data(ilo:ihi,jlo+2)  ) &
                                         + pd%data(ilo:ihi,jlo+3)

                    ! Upper y boundary
                    pd%data(ilo:ihi,jhi+1) = three*(  pd%data(ilo:ihi,jhi)      &
                                                    - pd%data(ilo:ihi,jhi-1)  ) &
                                           + pd%data(ilo:ihi,jhi-2)
                else
                    ! Lower x boundary (avoid west), lower y boundary (avoid south)
                    j = jlo
                    i = ilo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    en = phi%data(i  ,j+1)
                    ew = three*(p-ee)+phi%data(i+2,j  ) !two*p-ee
                    es = three*(p-en)+phi%data(i  ,j+2) !two*p-en

                    ene = phi%data(i+1,j+1)
                    enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                    ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                    esw = half*((three*(es-ese)+phi%data(i+2,j-1))+(three*(ew-enw)+phi%data(i-1,j+2))) !half*((two*es-ese)+(two*ew-enw))

                    pd%data(i,j  ) = scale*(ee-ew+ese-esw)

                    ! Interior to x, lower y boundary (avoid south)
                    do i = ilo+1, ihi-1
                        p  = phi%data(i  ,j  )
                        pe = phi%data(i+1,j  )
                        pw = phi%data(i-1,j  )
                        pn = phi%data(i  ,j+1)
                        ps = phi%data(i  ,j-1)

                        ee = phi%data(i+1,j  )
                        ew = phi%data(i-1,j  )
                        en = phi%data(i  ,j+1)
                        es = three*(p-en)+phi%data(i,j+2) !two*p-en

                        ene = phi%data(i+1,j+1)
                        enw = phi%data(i-1,j+1)
                        ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                        esw = three*(ew-enw)+phi%data(i-1,j+2) !two*ew-enw

                        pd%data(i,j  ) = scale*(ee-ew+ese-esw)
                    enddo !i

                    ! Upper x boundary (avoid east), lower y boundary (avoid south)
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    es = three*(p-en)+phi%data(i,j+2) !two*p-en

                    enw = phi%data(i-1,j+1)
                    ene = three*(pn-enw)+phi%data(i-2,j+1) !two*pn-enw
                    esw = three*(pw-enw)+phi%data(i-1,j+2) !two*pw-enw
                    ese = half*((three*(es-esw)+phi%data(i-2,j-1))+(three*(ee-ene)+phi%data(i+1,j+2))) !half*((two*es-esw)+(two*ee-ene))

                    pd%data(i,j  ) = scale*(ee-ew+ese-esw)


                    ! Upper y boundary...
                    j = jhi

                    ! Lower x boundary (avoid west), upper y boundary (avoid north)
                    i = ilo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    es = phi%data(i  ,j-1)
                    ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    ese = phi%data(i+1,j-1)
                    ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                    esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese
                    enw = half*((three*(ew-esw)+phi%data(i-1,j-2))+(three*(en-ene)+phi%data(i+2,j+1))) !half*((two*ew-esw)+(two*en-ene))

                    pd%data(i,j+1) = scale*(ene-enw+ee-ew)

                    ! Interior to x, upper y boundary (avoid north)
                    do i = ilo+1, ihi-1
                        p  = phi%data(i  ,j  )
                        pe = phi%data(i+1,j  )
                        pw = phi%data(i-1,j  )
                        pn = phi%data(i  ,j+1)
                        ps = phi%data(i  ,j-1)

                        ee = phi%data(i+1,j  )
                        ew = phi%data(i-1,j  )
                        es = phi%data(i  ,j-1)
                        en = three*(p-es)+phi%data(i,j-2) !two*p-es

                        ese = phi%data(i+1,j-1)
                        esw = phi%data(i-1,j-1)
                        ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                        enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw

                        pd%data(i,j+1) = scale*(ene-enw+ee-ew)
                    enddo !i

                    ! Upper x boundary (avoid east), upper y boundary (avoid north)
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    es = phi%data(i  ,j-1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    esw = phi%data(i-1,j-1)
                    enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw
                    ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw
                    ene = half*((three*(en-enw)+phi%data(i-2,j+1))+(three*(ee-ese)+phi%data(i+1,j-2))) !half*((two*en-enw)+(two*ee-ese))

                    pd%data(i,j+1) = scale*(ene-enw+ee-ew)
                endif
            endif
        else
            if (nodedir .eq. dir) then
                scale = one / phi%valid%dy
                pd%data(ilo:ihi,jlo:jhi+1) = scale * (  phi%data(ilo:ihi,jlo:jhi+1) &
                                                      - phi%data(ilo:ihi,jlo-1:jhi)  )
            else
                scale = fourth / phi%valid%dy

                ! Away from x boundaries...
                ! Interior
                pd%data(ilo+1:ihi,jlo+1:jhi-1) = scale * (  phi%data(ilo+1:ihi,jlo+2:jhi) &
                                                          - phi%data(ilo+1:ihi,jlo:jhi-2) &
                                                          + phi%data(ilo:ihi-1,jlo+2:jhi) &
                                                          - phi%data(ilo:ihi-1,jlo:jhi-2)  )

                ! Lower y boundary
                pd%data(ilo+1:ihi,jlo) = -scale * (  three*phi%data(ilo+1:ihi,jlo)   &
                                                   -  four*phi%data(ilo+1:ihi,jlo+1) &
                                                   +       phi%data(ilo+1:ihi,jlo+2) &
                                                   + three*phi%data(ilo:ihi-1,jlo)   &
                                                   -  four*phi%data(ilo:ihi-1,jlo+1) &
                                                   +       phi%data(ilo:ihi-1,jlo+2)  )

                ! Upper y boundary
                pd%data(ilo+1:ihi,jhi) = scale * (  three*phi%data(ilo+1:ihi,jhi)   &
                                                  -  four*phi%data(ilo+1:ihi,jhi-1) &
                                                  +       phi%data(ilo+1:ihi,jhi-2) &
                                                  + three*phi%data(ilo:ihi-1,jhi)   &
                                                  -  four*phi%data(ilo:ihi-1,jhi-1) &
                                                  +       phi%data(ilo:ihi-1,jhi-2)  )

                ! At x boundaries...
                if (simple_bdry_stencil) then
                    ! Lower x boundary
                    pd%data(ilo,jlo:jhi) = three*(  pd%data(ilo+1,jlo:jhi)    &
                                                  - pd%data(ilo+2,jlo:jhi)  ) &
                                         + pd%data(ilo+3,jlo:jhi)

                    ! Upper x boundary
                    pd%data(ihi+1,jlo:jhi) = three*(  pd%data(ihi,jlo:jhi)      &
                                                    - pd%data(ihi-1,jlo:jhi)  ) &
                                           + pd%data(ihi-2,jlo:jhi)
                else
                    ! Lower x boundary (avoid west), lower y boundary (avoid south)
                    i = ilo
                    j = jlo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    en = phi%data(i  ,j+1)
                    ew = three*(p-ee)+phi%data(i+2,j  ) !two*p-ee
                    es = three*(p-en)+phi%data(i  ,j+2) !two*p-en

                    ene = phi%data(i+1,j+1)
                    enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                    ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                    esw = half*((three*(es-ese)+phi%data(i+2,j-1))+(three*(ew-enw)+phi%data(i-1,j+2))) !half*((two*es-ese)+(two*ew-enw))

                    pd%data(i  ,j) = scale*(en-es+enw-esw)

                    do j = jlo+1, jhi-1
                        ! Lower x boundary (avoid west), interior to y
                        p  = phi%data(i  ,j  )
                        pe = phi%data(i+1,j  )
                        pw = phi%data(i-1,j  )
                        pn = phi%data(i  ,j+1)
                        ps = phi%data(i  ,j-1)

                        ee = phi%data(i+1,j  )
                        en = phi%data(i  ,j+1)
                        es = phi%data(i  ,j-1)
                        ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee

                        ene = phi%data(i+1,j+1)
                        ese = phi%data(i+1,j-1)
                        enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                        esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese

                        pd%data(i  ,j) = scale*(en-es+enw-esw)
                    enddo

                    ! Lower x boundary (avoid west), upper y boundary (avoid north)
                    j = jhi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    es = phi%data(i  ,j-1)
                    ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    ese = phi%data(i+1,j-1)
                    ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                    esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese
                    enw = half*((three*(ew-esw)+phi%data(i-1,j-2))+(three*(en-ene)+phi%data(i+2,j+1))) !half*((two*ew-esw)+(two*en-ene))

                    pd%data(i  ,j) = scale*(en-es+enw-esw)


                    ! Upper x boundary (avoid east), lower y boundary (avoid south)
                    i = ihi
                    j = jlo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    es = three*(p-en)+phi%data(i,j+2) !two*p-en

                    enw = phi%data(i-1,j+1)
                    ene = three*(pn-enw)+phi%data(i-2,j+1) !two*pn-enw
                    esw = three*(pw-enw)+phi%data(i-1,j+2) !two*pw-enw
                    ese = half*((three*(es-esw)+phi%data(i-2,j-1))+(three*(ee-ene)+phi%data(i+1,j+2))) !half*((two*es-esw)+(two*ee-ene))

                    pd%data(i+1,j) = scale*(ene-ese+en-es)


                    ! Upper x boundary (avoid east), interior to y
                    do j = jlo+1,jhi-1
                        p  = phi%data(i  ,j  )
                        pe = phi%data(i+1,j  )
                        pw = phi%data(i-1,j  )
                        pn = phi%data(i  ,j+1)
                        ps = phi%data(i  ,j-1)

                        ew = phi%data(i-1,j  )
                        en = phi%data(i  ,j+1)
                        es = phi%data(i  ,j-1)
                        ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew

                        enw = phi%data(i-1,j+1)
                        esw = phi%data(i-1,j-1)
                        ene = three*(en-enw)+phi%data(i-2,j+1) !two*en-enw
                        ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw

                        pd%data(i+1,j) = scale*(ene-ese+en-es)
                    enddo !j

                    ! Upper x boundary (avoid east), upper y boundary (avoid north)
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    es = phi%data(i  ,j-1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    esw = phi%data(i-1,j-1)
                    enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw
                    ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw
                    ene = half*((three*(en-enw)+phi%data(i-2,j+1))+(three*(ee-ese)+phi%data(i+1,j-2))) !half*((two*en-enw)+(two*ee-ese))

                    pd%data(i+1,j) = scale*(ene-ese+en-es)
                endif
            endif
        endif
    end subroutine compute_pd


    ! ------------------------------------------------------------------------------
    ! Computes J*Grad[Phi].
    ! phi is expected to be cell-centered and have at least 1 ghost layer.
    ! *flux is expected to be face-centered. No ghosts will be filled.
    ! bc* = 0 for Neum, 1 for Diri, 2 for Periodic, 3 for CF.
    ! xwk is workspace that must be node-centered in x and cell-centered in y.
    ! ywk is workspace that must be node-centered in y and cell-centered in x.
    !
    ! This function uses the EXACT same stencils as compute_laplacian as long as
    ! simple_bdry_stencil = .false. in compute_pd.
    ! ------------------------------------------------------------------------------
    subroutine compute_grad (xflux, yflux, phi, geo, bc, homog, xwk, ywk)
        type(box_data), intent(inout) :: xflux, yflux
        type(box_data), intent(inout) :: phi
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog
        type(box_data), intent(inout) :: xwk, ywk

        integer                       :: ilo, ihi, jlo, jhi, i, j

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        ! Fill ghosts
        call fill_ghosts (phi, bc, geo, homog, .false.)

        ! Compute xflux...
        call compute_pd (xflux, phi, 1)
        call compute_pd (xwk, phi, 2)
        do j = jlo, jhi
            do i = ilo, ihi+1
                xflux%data(i,j) = geo%Jgup_xx%data(i,j) * xflux%data(i,j) &
                                + geo%Jgup_xy%data(i,j) *   xwk%data(i,j)
            enddo
        enddo

        ! Compute yflux...
        call compute_pd (ywk, phi, 1)
        call compute_pd (yflux, phi, 2)
        do j = jlo, jhi+1
            do i = ilo, ihi
                yflux%data(i,j) = geo%Jgup_yx%data(i,j) *   ywk%data(i,j) &
                                + geo%Jgup_yy%data(i,j) * yflux%data(i,j)
            enddo
        enddo

        ! Set boundary fluxes
        call fill_boundary_fluxes (xflux, yflux, bc, homog)

    end subroutine compute_grad


    ! ------------------------------------------------------------------------------
    ! Computes J*Div[Flux].
    ! We multiply by J (that is, don't divide by J) because the rhs is expected to
    ! be scaled by J.
    ! *flux is expected to be face-centered with BCs computed.
    ! div is expected to be cell-centered with no ghosts.
    ! ------------------------------------------------------------------------------
    subroutine compute_div (div, xflux, yflux)
        type(box_data), intent(inout) :: div
        type(box_data), intent(in)    :: xflux, yflux

        real(dp)                      :: invdx, invdy
        integer                       :: ilo, ihi, jlo, jhi

        invdx = one / div%valid%dx
        invdy = one / div%valid%dy

        ilo = div%valid%ilo
        ihi = div%valid%ihi
        jlo = div%valid%jlo
        jhi = div%valid%jhi

        div%data(ilo:ihi,jlo:jhi) = &
              invdx * (xflux%data(ilo+1:ihi+1,jlo:jhi) - xflux%data(ilo:ihi,jlo:jhi)) &
            + invdy * (yflux%data(ilo:ihi,jlo+1:jhi+1) - yflux%data(ilo:ihi,jlo:jhi))

    end subroutine compute_div


    ! ------------------------------------------------------------------------------
    ! If opt_jscale is false (default), this function will not scale the result
    ! by 1/J. Note that a true Laplacian should scale by 1/J.
    ! ------------------------------------------------------------------------------
    subroutine compute_laplacian (lap, phi, geo, bc, homog, opt_jscale)
        type(box_data), intent(inout) :: lap
        type(box_data), intent(inout) :: phi
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog
        logical, intent(in), optional :: opt_jscale

        logical, parameter            :: neum_ghosts = .false.

        integer                       :: ilo, ihi, inner_ilo, inner_ihi, i
        integer                       :: jlo, jhi, inner_jlo, inner_jhi, j
        real(dp)                      :: p,pe,pn,pw,ps
        real(dp)                      :: ee,en,ew,es
        real(dp)                      :: ene,enw,esw,ese
        real(dp)                      :: xfe,xfw,yfn,yfs
        real(dp)                      :: idx, idy

        if ((bc%mode_xlo .eq. BCMODE_LAH) .or. (bc%mode_xhi .eq. BCMODE_LAH) .or. &
            (bc%mode_ylo .eq. BCMODE_LAH) .or. (bc%mode_yhi .eq. BCMODE_LAH)) then
            print*, 'Poisson2d::compute_invdiags received BCMODE_LAH.'
            stop
        endif

        ilo = lap%valid%ilo
        ihi = lap%valid%ihi
        jlo = lap%valid%jlo
        jhi = lap%valid%jhi

        if (bc%type_xlo .eq. BCTYPE_PERIODIC) then
            inner_ilo = ilo
            inner_ihi = ihi
        else
            inner_ilo = ilo+1
            inner_ihi = ihi-1
        endif

        if (bc%type_ylo .eq. BCTYPE_PERIODIC) then
            inner_jlo = jlo
            inner_jhi = jhi
        else
            inner_jlo = jlo+1
            inner_jhi = jhi-1
        endif

        idx = one / geo%dx
        idy = one / geo%dy

        ! Fill ghost cells (except at Neum BCs)
        call fill_ghosts (phi, bc, geo, homog, neum_ghosts)

        ! Lower x boundary (avoid west), lower y boundary (avoid south)
        if (bc%type_ylo .ne. BCTYPE_PERIODIC) then
            j = jlo
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                en = phi%data(i  ,j+1)
                ew = three*(p-ee)+phi%data(i+2,j  ) !two*p-ee
                es = three*(p-en)+phi%data(i  ,j+2) !two*p-en

                ene = phi%data(i+1,j+1)
                enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                esw = half*((three*(es-ese)+phi%data(i+2,j-1))+(three*(ew-enw)+phi%data(i-1,j+2))) !half*((two*es-ese)+(two*ew-enw))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif
                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x, lower y boundary (avoid south)
            do i = inner_ilo, inner_ihi
                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                ew = phi%data(i-1,j  )
                en = phi%data(i  ,j+1)
                es = three*(p-en)+phi%data(i,j+2) !two*p-en

                ene = phi%data(i+1,j+1)
                enw = phi%data(i-1,j+1)
                ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                esw = three*(ew-enw)+phi%data(i-1,j+2) !two*ew-enw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), lower y boundary (avoid south)
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ew = phi%data(i-1,j  )
                en = phi%data(i  ,j+1)
                ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                es = three*(p-en)+phi%data(i,j+2) !two*p-en

                enw = phi%data(i-1,j+1)
                ene = three*(pn-enw)+phi%data(i-2,j+1) !two*pn-enw
                esw = three*(pw-enw)+phi%data(i-1,j+2) !two*pw-enw
                ese = half*((three*(es-esw)+phi%data(i-2,j-1))+(three*(ee-ene)+phi%data(i+1,j+2))) !half*((two*es-esw)+(two*ee-ene))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif
                if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                    if (homog) then
                        yfs = zero
                    else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                        yfs = bc%data_ylo(1)
                    else
                        yfs = bc%data_ylo(i)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        endif ! not periodic in y

        ! Interior to y...
        do j = inner_jlo, inner_jhi
            ! Lower x boundary (avoid west), interior to y
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                en = phi%data(i  ,j+1)
                es = phi%data(i  ,j-1)
                ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee

                ene = phi%data(i+1,j+1)
                ese = phi%data(i+1,j-1)
                enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x and y
            do i = inner_ilo, inner_ihi
                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                ew = phi%data(i-1,j  )
                en = phi%data(i  ,j+1)
                es = phi%data(i  ,j-1)

                ene = phi%data(i+1,j+1)
                enw = phi%data(i-1,j+1)
                ese = phi%data(i+1,j-1)
                esw = phi%data(i-1,j-1)

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), interior to y
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ew = phi%data(i-1,j  )
                en = phi%data(i  ,j+1)
                es = phi%data(i  ,j-1)
                ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew

                enw = phi%data(i-1,j+1)
                esw = phi%data(i-1,j-1)
                ene = three*(en-enw)+phi%data(i-2,j+1) !two*en-enw
                ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        enddo !j


        ! Upper y boundary...
        if (bc%type_yhi .ne. BCTYPE_PERIODIC) then
            j = jhi

            ! Lower x boundary (avoid west), upper y boundary (avoid north)
            if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                i = ilo

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                es = phi%data(i  ,j-1)
                ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee
                en = three*(p-es)+phi%data(i,j-2) !two*p-es

                ese = phi%data(i+1,j-1)
                ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese
                enw = half*((three*(ew-esw)+phi%data(i-1,j-2))+(three*(en-ene)+phi%data(i+2,j+1))) !half*((two*ew-esw)+(two*en-ene))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                    if (homog) then
                        xfw = zero
                    else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                        xfw = bc%data_xlo(1)
                    else
                        xfw = bc%data_xlo(j)
                    endif
                endif
                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x

            ! Interior to x, upper y boundary (avoid north)
            do i = inner_ilo, inner_ihi
                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ee = phi%data(i+1,j  )
                ew = phi%data(i-1,j  )
                es = phi%data(i  ,j-1)
                en = three*(p-es)+phi%data(i,j-2) !two*p-es

                ese = phi%data(i+1,j-1)
                esw = phi%data(i-1,j-1)
                ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo !i

            ! Upper x boundary (avoid east), upper y boundary (avoid north)
            if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                i = ihi

                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                ew = phi%data(i-1,j  )
                es = phi%data(i  ,j-1)
                ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                en = three*(p-es)+phi%data(i,j-2) !two*p-es

                esw = phi%data(i-1,j-1)
                enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw
                ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw
                ene = half*((three*(en-enw)+phi%data(i-2,j+1))+(three*(ee-ese)+phi%data(i+1,j-2))) !half*((two*en-enw)+(two*ee-ese))

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                    if (homog) then
                        xfe = zero
                    else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                        xfe = bc%data_xhi(1)
                    else
                        xfe = bc%data_xhi(j)
                    endif
                endif
                if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                    if (homog) then
                        yfn = zero
                    else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                        yfn = bc%data_yhi(1)
                    else
                        yfn = bc%data_yhi(j)
                    endif
                endif

                lap%data(i,j) = (xfe-xfw)*idx + (yfn-yfs)*idy
            endif ! not periodic in x
        endif ! not periodic in y

        ! Scale by 1/J if necessary.
        if (present(opt_jscale)) then
            if (opt_jscale) then
                lap%data(ilo:ihi,jlo:jhi) = lap%data(ilo:ihi,jlo:jhi) &
                                          / geo%J%data(ilo:ihi,jlo:jhi)
            endif
        endif

    end subroutine compute_laplacian


    ! ------------------------------------------------------------------------------
    ! Compute res = rhs - L[phi].
    ! If opt_jscale is false (default), this function will not scale the result
    ! by 1/J. Note that a true Laplacian should do this.
    ! ------------------------------------------------------------------------------
    subroutine compute_residual (res, rhs, phi, geo, bc, homog, opt_jscale)
        type(box_data), intent(inout) :: res, phi
        type(box_data), intent(in)    :: rhs
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        logical, intent(in), optional :: opt_jscale
        logical                       :: jscale

        integer                       :: ilo, ihi, jlo, jhi

        if (present(opt_jscale)) then
            jscale = opt_jscale
        else
            jscale = .false.
        endif

        call compute_laplacian (res, phi, geo, bc, homog, jscale)

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi
        res%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) - res%data(ilo:ihi,jlo:jhi)

    end subroutine compute_residual


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine relax_jacobi (phi, rhs, geo, bc, homog, invdiags, &
                             omega, tol, maxiters, zerophi, verbosity)
        type(box_data), intent(inout)   :: phi
        type(box_data), intent(in)      :: rhs
        type(geo_data), intent(in)      :: geo
        type(bdry_data), intent(in)     :: bc
        logical, intent(in)             :: homog
        type(box_data), intent(in)      :: invdiags
        real(dp), intent(in)            :: omega
        real(dp), intent(in)            :: tol
        integer, intent(in)             :: maxiters
        logical, intent(in)             :: zerophi
        integer, intent(in)             :: verbosity

        type(box_data)                  :: r
        real(dp)                        :: rscale
        real(dp), dimension(0:maxiters) :: relres
        real(dp)                        :: newphi
        integer                         :: iter
        integer                         :: ilo, ihi, i
        integer                         :: jlo, jhi, j
        real(dp)                        :: sum

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        ! Do we even need to be here?
        if (maxIters .eq. 0) then
            return
        endif

        ! Allocate workspace
        call define_box_data (r, rhs)

        ! Initialize phi to zero if necessary
        if (zerophi) then
            phi%data = zero
        endif

        ! Compute initial residual
        call compute_residual (r, rhs, phi, geo, bc, homog)

        if (tol .gt. zero) then
            rscale = pnorm (r, r%valid, 2)
            if (rscale .eq. zero) then
                relres(0) = rscale
                rscale = one
            else
                relres(0) = one
            endif

            if (verbosity .ge. 3) then
                print*, 'scale |res| = ', rscale
                print*, 'iter ', 0, ': rel |res| = ', relres(0)
            endif
        endif

        ! Iterate
        do iter = 1, maxIters
            if (verbosity .ge. 5) then
                sum = integrate2d (r, r%valid, geo, .false.)
                print*, ' sum rhs = ', sum
            endif

            ! Update phi
            do j = jlo, jhi
                do i = ilo, ihi
                    newphi = phi%data(i,j) + r%data(i,j)*invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                enddo
            enddo

            ! Compute new residual
            call compute_residual (r, rhs, phi, geo, bc, homog)

            if (tol .gt. zero) then
                relres(iter) = pnorm (r, r%valid, 2) / rscale
                if (verbosity .ge. 3) then
                    print*, 'iter ', iter, ': rel |res| = ', relres(iter)
                endif

                ! Did we converge?
                if (relres(iter) .le. tol) then
                    if (verbosity .ge. 3) then
                        print*, "Converged."
                    endif
                    exit
                endif
            endif
        enddo

    end subroutine relax_jacobi


    ! ------------------------------------------------------------------------------
    !             color = 0
    !             do j = jlo, jhi
    !                 imin = ilo + mod(color+j, 2) ! Removed +1
    !                 do i = imin, ihi, 2
    ! ------------------------------------------------------------------------------
    subroutine relax_gs (phi, rhs, geo, bc, homog, invdiags, &
                         omega, tol, maxiters, zerophi, verbosity)
        type(box_data), intent(inout)   :: phi
        type(box_data), intent(in)      :: rhs
        type(geo_data), intent(in)      :: geo
        type(bdry_data), intent(in)     :: bc
        logical, intent(in)             :: homog
        type(box_data), intent(in)      :: invdiags
        real(dp), intent(in)            :: omega
        real(dp), intent(in)            :: tol
        integer, intent(in)             :: maxiters
        logical, intent(in)             :: zerophi
        integer, intent(in)             :: verbosity

        type(box_data)                  :: r
        real(dp)                        :: rscale
        real(dp), dimension(0:maxiters) :: relres
        integer                         :: iter

        integer                         :: ilo, ihi, inner_ilo, inner_ihi, i
        integer                         :: jlo, jhi, inner_jlo, inner_jhi, j
        real(dp)                        :: p,pe,pn,pw,ps
        real(dp)                        :: ee,en,ew,es
        real(dp)                        :: ene,enw,esw,ese
        real(dp)                        :: xfe,xfw,yfn,yfs
        real(dp)                        :: idx, idy
        real(dp)                        :: lphi, newphi, sum

        ! Do we even need to be here?
        if (maxiters .eq. 0) then
            return
        endif

        if ((bc%mode_xlo .eq. BCMODE_LAH) .or. (bc%mode_xhi .eq. BCMODE_LAH) .or. &
            (bc%mode_ylo .eq. BCMODE_LAH) .or. (bc%mode_yhi .eq. BCMODE_LAH)) then
            print*, 'Poisson2d::compute_invdiags received BCMODE_LAH.'
            stop
        endif

        ! Define some useful quantities
        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        if (bc%type_xlo .eq. BCTYPE_PERIODIC) then
            inner_ilo = ilo
            inner_ihi = ihi
        else
            inner_ilo = ilo+1
            inner_ihi = ihi-1
        endif

        if (bc%type_ylo .eq. BCTYPE_PERIODIC) then
            inner_jlo = jlo
            inner_jhi = jhi
        else
            inner_jlo = jlo+1
            inner_jhi = jhi-1
        endif

        idx = one / geo%dx
        idy = one / geo%dy

        ! Initialize phi to zero
        if (zerophi) then
            phi%data = zero
        endif

        ! Compute initial residual
        if (tol .gt. zero) then
            call define_box_data (r, rhs)
            call compute_residual (r, rhs, phi, geo, bc, homog)

            rscale = pnorm (r, r%valid, 2)
            if (rscale .eq. zero) then
                relres(0) = rscale
                rscale = one
            else
                relres(0) = one
            endif

            if (verbosity .ge. 3) then
                print*, 'scale |res| = ', rscale
                print*, 'iter ', 0, ': rel |res| = ', relres(0)
            endif
        endif

        ! Iterate
        do iter = 1, maxiters
            if (verbosity .ge. 5) then
                sum = integrate2d (r, r%valid, geo, .false.)
                print*, ' sum rhs = ', sum
            endif

            ! Update phi via standard Gauss-Seidel
            call fill_ghosts (phi, bc, geo, homog, .false.)

            ! The main comutation...

            ! Lower x boundary (avoid west), lower y boundary (avoid south)
            if (bc%type_ylo .ne. BCTYPE_PERIODIC) then
                j = jlo
                if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                    i = ilo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    en = phi%data(i  ,j+1)
                    ew = three*(p-ee)+phi%data(i+2,j  ) !two*p-ee
                    es = three*(p-en)+phi%data(i  ,j+2) !two*p-en

                    ene = phi%data(i+1,j+1)
                    enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                    ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                    esw = half*((three*(es-ese)+phi%data(i+2,j-1))+(three*(ew-enw)+phi%data(i-1,j+2))) !half*((two*es-ese)+(two*ew-enw))

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                        if (homog) then
                            xfw = zero
                        else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                            xfw = bc%data_xlo(1)
                        else
                            xfw = bc%data_xlo(j)
                        endif
                    endif
                    if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                        if (homog) then
                            yfs = zero
                        else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                            yfs = bc%data_ylo(1)
                        else
                            yfs = bc%data_ylo(i)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x

                ! Interior to x, lower y boundary (avoid south)
                do i = inner_ilo, inner_ihi
                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    es = three*(p-en)+phi%data(i,j+2) !two*p-en

                    ene = phi%data(i+1,j+1)
                    enw = phi%data(i-1,j+1)
                    ese = three*(ee-ene)+phi%data(i+1,j+2) !two*ee-ene
                    esw = three*(ew-enw)+phi%data(i-1,j+2) !two*ew-enw

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                        if (homog) then
                            yfs = zero
                        else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                            yfs = bc%data_ylo(1)
                        else
                            yfs = bc%data_ylo(i)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                enddo !i

                ! Upper x boundary (avoid east), lower y boundary (avoid south)
                if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    es = three*(p-en)+phi%data(i,j+2) !two*p-en

                    enw = phi%data(i-1,j+1)
                    ene = three*(pn-enw)+phi%data(i-2,j+1) !two*pn-enw
                    esw = three*(pw-enw)+phi%data(i-1,j+2) !two*pw-enw
                    ese = half*((three*(es-esw)+phi%data(i-2,j-1))+(three*(ee-ene)+phi%data(i+1,j+2))) !half*((two*es-esw)+(two*ee-ene))

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                        if (homog) then
                            xfe = zero
                        else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                            xfe = bc%data_xhi(1)
                        else
                            xfe = bc%data_xhi(j)
                        endif
                    endif
                    if (bc%type_ylo .eq. BCTYPE_NEUM) then ! south terms are zero
                        if (homog) then
                            yfs = zero
                        else if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                            yfs = bc%data_ylo(1)
                        else
                            yfs = bc%data_ylo(i)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x
            endif ! not periodic in y

            ! Interior to y...
            do j = inner_jlo, inner_jhi
                ! Lower x boundary (avoid west), interior to y
                if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                    i = ilo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    en = phi%data(i  ,j+1)
                    es = phi%data(i  ,j-1)
                    ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee

                    ene = phi%data(i+1,j+1)
                    ese = phi%data(i+1,j-1)
                    enw = three*(en-ene)+phi%data(i+2,j+1) !two*en-ene
                    esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                        if (homog) then
                            xfw = zero
                        else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                            xfw = bc%data_xlo(1)
                        else
                            xfw = bc%data_xlo(j)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x

                ! Interior to x and y
                do i = inner_ilo, inner_ihi
                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    es = phi%data(i  ,j-1)

                    ene = phi%data(i+1,j+1)
                    enw = phi%data(i-1,j+1)
                    ese = phi%data(i+1,j-1)
                    esw = phi%data(i-1,j-1)

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                enddo !i

                ! Upper x boundary (avoid east), interior to y
                if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    en = phi%data(i  ,j+1)
                    es = phi%data(i  ,j-1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew

                    enw = phi%data(i-1,j+1)
                    esw = phi%data(i-1,j-1)
                    ene = three*(en-enw)+phi%data(i-2,j+1) !two*en-enw
                    ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                        if (homog) then
                            xfe = zero
                        else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                            xfe = bc%data_xhi(1)
                        else
                            xfe = bc%data_xhi(j)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x
            enddo !j


            ! Upper y boundary...
            if (bc%type_yhi .ne. BCTYPE_PERIODIC) then
                j = jhi

                ! Lower x boundary (avoid west), upper y boundary (avoid north)
                if (bc%type_xlo .ne. BCTYPE_PERIODIC) then
                    i = ilo

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    es = phi%data(i  ,j-1)
                    ew = three*(p-ee)+phi%data(i+2,j) !two*p-ee
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    ese = phi%data(i+1,j-1)
                    ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                    esw = three*(es-ese)+phi%data(i+2,j-1) !two*es-ese
                    enw = half*((three*(ew-esw)+phi%data(i-1,j-2))+(three*(en-ene)+phi%data(i+2,j+1))) !half*((two*ew-esw)+(two*en-ene))

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xlo .eq. BCTYPE_NEUM) then ! west terms are zero
                        if (homog) then
                            xfw = zero
                        else if (bc%mode_xlo .eq. BCMODE_UNIFORM) then
                            xfw = bc%data_xlo(1)
                        else
                            xfw = bc%data_xlo(j)
                        endif
                    endif
                    if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                        if (homog) then
                            yfn = zero
                        else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                            yfn = bc%data_yhi(1)
                        else
                            yfn = bc%data_yhi(j)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x

                ! Interior to x, upper y boundary (avoid north)
                do i = inner_ilo, inner_ihi
                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ee = phi%data(i+1,j  )
                    ew = phi%data(i-1,j  )
                    es = phi%data(i  ,j-1)
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    ese = phi%data(i+1,j-1)
                    esw = phi%data(i-1,j-1)
                    ene = three*(ee-ese)+phi%data(i+1,j-2) !two*ee-ese
                    enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                        if (homog) then
                            yfn = zero
                        else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                            yfn = bc%data_yhi(1)
                        else
                            yfn = bc%data_yhi(j)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                enddo !i

                ! Upper x boundary (avoid east), upper y boundary (avoid north)
                if (bc%type_xhi .ne. BCTYPE_PERIODIC) then
                    i = ihi

                    p  = phi%data(i  ,j  )
                    pe = phi%data(i+1,j  )
                    pw = phi%data(i-1,j  )
                    pn = phi%data(i  ,j+1)
                    ps = phi%data(i  ,j-1)

                    ew = phi%data(i-1,j  )
                    es = phi%data(i  ,j-1)
                    ee = three*(p-ew)+phi%data(i-2,j) !two*p-ew
                    en = three*(p-es)+phi%data(i,j-2) !two*p-es

                    esw = phi%data(i-1,j-1)
                    enw = three*(ew-esw)+phi%data(i-1,j-2) !two*ew-esw
                    ese = three*(es-esw)+phi%data(i-2,j-1) !two*es-esw
                    ene = half*((three*(en-enw)+phi%data(i-2,j+1))+(three*(ee-ese)+phi%data(i+1,j-2))) !half*((two*en-enw)+(two*ee-ese))

                    xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx + geo%Jgup_xy%data(i+1,j)*fourth*(ene-ese+en-es)*idy
                    xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx + geo%Jgup_xy%data(i  ,j)*fourth*(en-es+enw-esw)*idy
                    yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy + geo%Jgup_yx%data(i,j+1)*fourth*(ene-enw+ee-ew)*idx
                    yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy + geo%Jgup_yx%data(i,j  )*fourth*(ee-ew+ese-esw)*idx

                    if (bc%type_xhi .eq. BCTYPE_NEUM) then ! east terms are zero
                        if (homog) then
                            xfe = zero
                        else if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                            xfe = bc%data_xhi(1)
                        else
                            xfe = bc%data_xhi(j)
                        endif
                    endif
                    if (bc%type_yhi .eq. BCTYPE_NEUM) then ! north terms are zero
                        if (homog) then
                            yfn = zero
                        else if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                            yfn = bc%data_yhi(1)
                        else
                            yfn = bc%data_yhi(j)
                        endif
                    endif

                    lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                    newphi = phi%data(i,j) + (rhs%data(i,j) - lphi) * invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                endif ! not periodic in x
            endif ! not periodic in y

            ! Diagnostics
            if (tol .gt. zero) then
                ! Compute new residual
                ! TODO: Once this function has been corrected,
                ! this can move inside the i(tol>0) block.
                call compute_residual (r, rhs, phi, geo, bc, homog)

                relres(iter) = pnorm (r, r%valid, 2) / rscale
                if (verbosity .ge. 3) then
                    print*, 'iter ', iter, ': rel |res| = ', relres(iter)
                endif

                ! Did we converge?
                if (relres(iter) .le. tol) then
                    if (verbosity .ge. 3) then
                        print*, "Converged."
                    endif
                    exit
                endif
            endif
        enddo

        ! Free memory
        if (tol .gt. zero) then
            call undefine_box_data (r)
        endif

    end subroutine relax_gs


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine solve_bicgstab (phi, rhs, geo, bc, homog, invdiags, &
                               tol, max_iters, max_restarts, zerophi, verbosity)
        type(box_data), intent(inout)                 :: phi
        type(box_data), intent(in)                    :: rhs
        type(geo_data), intent(in)                    :: geo
        type(bdry_data), intent(in)                   :: bc
        logical, intent(in)                           :: homog
        type(box_data), intent(in)                    :: invdiags
        real(dp), intent(in)                          :: tol
        integer, intent(in)                           :: max_iters, max_restarts
        logical, intent(in)                           :: zerophi
        integer, intent(in)                           :: verbosity

        integer                                       :: ilo, ihi
        integer                                       :: jlo, jhi
        type(box_data)                                :: r, r0, nu, p, t, y, z
        real(dp)                                      :: rscale, sum
        real(dp), dimension(0:max_iters+max_restarts) :: rho, omega, relres
        real(dp)                                      :: alpha, beta, lastres
        integer                                       :: iter, i, num_restarts
        logical                                       :: is_restart

        real(dp), parameter                           :: hang = 1.0E-8_dp
        integer, parameter                            :: norm_type = 0

        ! Do we even need to be here?
        if (max_iters .eq. 0) then
            return
        endif

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        ! Allocate workspace
        call define_box_data (r, phi)
        call define_box_data (r0, phi)
        call define_box_data (nu, phi)
        call define_box_data (p, phi)
        call define_box_data (t, phi)
        call define_box_data (y, phi)
        call define_box_data (z, phi)

        ! Initialize phi to zero
        if (zerophi) then
            phi%data = zero
        endif

        is_restart = .false.
        i = 0
        relres = zero

        ! Compute initial residual
        call compute_residual (r, rhs, phi, geo, bc, homog)
        r0%data = r%data
        rscale = pnorm (r0, r0%valid, norm_type)
        relres(0) = one

        if (rscale .eq. zero) then
            ! Free memory and exit
            call undefine_box_data (r)
            call undefine_box_data (r0)
            call undefine_box_data (nu)
            call undefine_box_data (p)
            call undefine_box_data (t)
            call undefine_box_data (y)
            call undefine_box_data (z)
            return
        endif

        if (verbosity .ge. 3) then
            print*, 'scale |res| = ', rscale
            print*, 'iter ', 0, ': rel |res| = ', relres(0)
        endif

        ! Initialize all other workspace variables
        alpha = one
        rho(i) = one
        omega(i) = one
        nu%data = zero
        p%data = zero
        num_restarts = 0

        ! Iterate...
        do iter = 1, max_iters
            ! Increment index for bookkeeping vars.
            i = i + 1

            if (verbosity .ge. 5) then
                sum = integrate2d (r, r%valid, geo, .false.)
                print*, ' sum rhs = ', sum
            endif

            rho(i) = inner_prod (r0, r)
            beta = (rho(i) / rho(i-1)) * (alpha / omega(i-1))
            if (beta .eq. zero) then
                print*, 'solve_bicgstab: beta is zero. Probably due to a zero rhs.'
                stop
            endif
            p%data = beta*p%data
            p%data = p%data                  &
                   - beta*omega(i-1)*nu%data &
                   + r%data

            ! Preconditioner
            y%data(ilo:ihi,jlo:jhi) = p%data(ilo:ihi,jlo:jhi) * invdiags%data(ilo:ihi,jlo:jhi)
            call relax_gs (y, p, geo, bc, homog, invdiags, &
                           one-third, & ! omega
                           -one,      & ! tol
                           2,         & ! maxiters
                           .false.,   & ! zerophi
                           0)           ! verbosity

            call compute_laplacian (nu, y, geo, bc, homog)
            alpha = inner_prod (r0, nu)
            alpha = rho(i) / alpha
            r%data = r%data - alpha*nu%data

            ! If |r| is small, set phi = phi + alpha*p and quit

            ! A Preconditioner
            z%data(ilo:ihi,jlo:jhi) = r%data(ilo:ihi,jlo:jhi) * invdiags%data(ilo:ihi,jlo:jhi)
            call relax_gs (z, r, geo, bc, homog, invdiags, &
                           one-third, & ! omega
                           -one,      & ! tol
                           2,         & ! maxiters
                           .false.,   & ! zerophi
                           0)           ! verbosity

            call compute_laplacian (t, z, geo, bc, homog)
            omega(i) = inner_prod (t, r) / inner_prod (t, t)

            ! This would also change with a preconditioner
            phi%data = phi%data         &
                     + alpha * y%data   &
                     + omega(i) * z%data

            ! Compute new residual
            r%data = r%data - omega(i)*t%data

            ! If this is a restart, we expect the residual to rise.
            ! Don't let this stop the solver from proceeding.
            if (.not. is_restart) then
                lastres = relres(i-1)
            else
                lastres = 1.0E200_dp
            endif

            ! Check if we are at tol
            relres(i) = pnorm (r, r%valid, norm_type) / rscale
            if (verbosity .ge. 3) then
                print*, 'iter ', iter, ': rel |res| = ', relres(i)
            endif

            ! Did we converge?
            if (relres(i) .le. tol) then
                if (verbosity .ge. 3) then
                    print*, "Converged."
                endif
                exit
            endif

            ! Are we hanging?
            if (abs(relres(i) - lastres) .lt. tol*hang) then
                if (num_restarts .lt. max_restarts) then
                    ! The act of restarting will produce a new residual which we
                    ! would like to include in our bookkeeping, so we increase i,
                    ! recompute the residual, and reset all other bookkeeping vars.

                    ! Increment
                    num_restarts = num_restarts + 1
                    i = i + 1

                    ! Compute new residual
                    call compute_residual (r, rhs, phi, geo, bc, homog)
                    r0%data = r%data
                    relres(i) = pnorm (r, r%valid, norm_type) / rscale
                    if (verbosity .ge. 3) then
                        print*, "Hanging, restart number ", num_restarts, ', current rel |res| = ', relres(i)
                    endif

                    ! Reset bookkeeping variables
                    alpha = one
                    rho(i) = one
                    omega(i) = one
                    nu%data = zero
                    p%data = zero

                    ! Start new iteration
                    is_restart = .true.
                    cycle
                else
                    if (verbosity .ge. 3) then
                        print*, "Hanging, max restarts reached."
                    endif
                    exit
                endif
            endif

            ! Are we diverging?
            if (relres(i) .gt. lastres) then
                if (num_restarts .lt. max_restarts) then
                    ! The act of restarting will produce a new residual which we
                    ! would like to include in our bookkeeping, so we increase i,
                    ! recompute the residual, and reset all other bookkeeping vars.

                    ! Increment
                    num_restarts = num_restarts + 1
                    i = i + 1

                    ! Compute new residual
                    call compute_residual (r, rhs, phi, geo, bc, homog)
                    r0%data = r%data
                    relres(i) = pnorm (r0, r0%valid, norm_type) / rscale
                    if (verbosity .ge. 3) then
                        print*, "Hanging, restart number ", num_restarts, ', current rel |res| = ', relres(i)
                    endif

                    ! Reset bookkeeping variables
                    alpha = one
                    rho(i) = one
                    omega(i) = one
                    nu%data = zero
                    p%data = zero

                    ! Start new iteration
                    is_restart = .true.
                    cycle
                else
                    if (verbosity .ge. 3) then
                        print*, 'Diverging.'
                    endif

                    ! Undo last correction
                    ! TODO: It would be better to remember the best solution
                    r%data = r%data + omega(i)*t%data
                    phi%data = phi%data         &
                             - alpha * y%data   &
                             - omega(i) * z%data
                    exit
                endif
            endif

            is_restart = .false.
        enddo

        ! Free memory
        call undefine_box_data (r)
        call undefine_box_data (r0)
        call undefine_box_data (nu)
        call undefine_box_data (p)
        call undefine_box_data (t)
        call undefine_box_data (y)
        call undefine_box_data (z)
    end subroutine solve_bicgstab


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine solve_bicgstab2 (phi, rhs, geo, bc, homog, invdiags, &
                                tol, max_iters, max_restarts, zerophi, verbosity)
        type(box_data), intent(inout)                 :: phi
        type(box_data), intent(in)                    :: rhs
        type(box_data), intent(in)                    :: invdiags
        type(geo_data), intent(in)                    :: geo
        type(bdry_data), intent(in)                   :: bc
        logical, intent(in)                           :: homog
        real(dp), intent(in)                          :: tol
        integer, intent(in)                           :: max_iters, max_restarts
        logical, intent(in)                           :: zerophi
        integer, intent(in)                           :: verbosity

        integer                                       :: ilo, ihi
        integer                                       :: jlo, jhi
        type(box_data)                                :: r, r_tilde, e, p, p_tilde, s_tilde, t, v
        real(dp), dimension(0:3)                      :: rho
        real(dp), dimension(0:1)                      :: norm, alpha, beta, omega
        real(dp)                                      :: initial_norm, initial_rnorm, m
        integer                                       :: recount, i, restarts
        logical                                       :: init

        real(dp), parameter                           :: reps = 1.0E-12_dp
        real(dp), parameter                           :: hang = 1.0E-7_dp
        real(dp), parameter                           :: small = 1.0E-30_dp
        integer, parameter                            :: norm_type = 2

        ! Do we even need to be here?
        if (max_iters .eq. 0) then
            return
        endif

        ! Initialize phi to zero
        if (zerophi) then
            phi%data = zero
        endif

        ! Allocate workspace
        call define_box_data (e, phi)
        call define_box_data (p_tilde, phi)
        call define_box_data (s_tilde, phi)
        call define_box_data (r, rhs)
        call define_box_data (r_tilde, rhs)
        call define_box_data (p, rhs)
        call define_box_data (t, rhs)
        call define_box_data (v, rhs)

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        recount = 0

        call compute_residual (r, rhs, phi, geo, bc, homog)

        r_tilde%data = r%data
        e%data = zero
        p_tilde%data = zero
        s_tilde%data = zero

        i = 0

        norm(0) = pnorm(r, r%valid, norm_type)
        initial_norm = norm(0)
        initial_rnorm = norm(0)
        norm(1) = norm(0)

        rho = zero
        alpha = zero
        beta = zero
        omega = zero

        init = .true.
        restarts = 0

        if (verbosity .ge. 3) then
            ! print*, 'scale |res| = ', rscale
            print*, 'iter ', 0, ': rel |res| = ', initial_norm
        endif

        do while (((i .lt. max_iters) .and. (norm(0) .gt. tol*norm(1)))  .and. (norm(1) .gt. zero))
            i = i + 1

            norm(1) = norm(0)
            alpha(1) = alpha(0)
            beta(1) = beta(0)
            omega(1) = omega(0)

            rho(3) = rho(2)
            rho(2) = rho(1)
            rho(1) = inner_prod (r_tilde, r)

            if (rho(1) .eq. zero) then
                phi%data = phi%data + e%data

                if (verbosity .ge. 3) then
                    print*, 'rho = 0, returning.'
                    print*, 'residual norm = ', norm(0)
                endif

                ! Free memory
                call undefine_box_data (e)
                call undefine_box_data (p_tilde)
                call undefine_box_data (s_tilde)
                call undefine_box_data (r)
                call undefine_box_data (r_tilde)
                call undefine_box_data (p)
                call undefine_box_data (t)
                call undefine_box_data (v)

                return
            endif

            if (init) then
                p%data = r%data
                init = .false.
            else
                beta(1) = (rho(1)/rho(2))*(alpha(1)/omega(1))
                p%data = beta(1)*(p%data-omega(1)*v%data) + r%data
            endif

            ! Precond (p_rilde, p)
            p_tilde%data(ilo:ihi,jlo:jhi) = p%data(ilo:ihi,jlo:jhi) * invdiags%data(ilo:ihi,jlo:jhi)
            call relax_gs (p_tilde, p, geo, bc, .true., invdiags, &
                           one-third,     &  ! omega
                           -one,    &  ! tol
                           2,       &  ! maxiters
                           .false.,  & ! zerophi
                           0)          ! verbosity

            call compute_laplacian (v, p_tilde, geo, bc, .true.)
            m = inner_prod (r_tilde, v)
            alpha(0) = rho(1)/m

            if (abs(m) .gt. small*abs(rho(1))) then
                r%data = r%data - alpha(0)*v%data
                norm(0) = pnorm (r, r%valid, norm_type)
                e%data = e%data + alpha(0)*p_tilde%data
            else
                r%data = zero
                norm(0) = zero
            endif

            if ((norm(0) .gt. tol*initial_norm) .and. (norm(0) .gt. reps*initial_rnorm)) then
                ! Precond (s_tilde, r)
                s_tilde%data(ilo:ihi,jlo:jhi) = r%data(ilo:ihi,jlo:jhi) * invdiags%data(ilo:ihi,jlo:jhi)
                call relax_gs (s_tilde, r, geo, bc, .true., invdiags, &
                               one-third,     &  ! omega
                               -one,    &  ! tol
                               2,       &  ! maxiters
                               .false.,  & ! zerophi
                               0)          ! verbosity

                call compute_laplacian (t, s_tilde, geo, bc, .true.)
                omega(0) = inner_prod (t, r) / inner_prod (t, t)
                e%data = e%data + omega(0)*s_tilde%data
                r%data = r%data - omega(0)*t%data
                norm(0) = pnorm(r, r%valid, norm_type)
            endif

            if (verbosity .ge. 3) then
                print*, 'iter ', i, ': rel |res| = ', norm(0)
            endif

            if ((norm(0) .le. tol*initial_norm) .or. (norm(0) .le. reps*initial_rnorm)) then
                if (verbosity .ge. 3) then
                    print*, "Converged."
                endif
                exit
            endif

            if ((omega(0) .eq. zero) .or. (norm(0) .gt. (one-hang)*norm(1))) then
                if (recount .eq. 0) then
                    recount = 1
                else
                    recount = 0
                    phi%data = phi%data + e%data

                    if (restarts .eq. max_restarts) then
                        if (verbosity .ge. 3) then
                            print*, 'Max restarts reached.'
                        endif
                        exit
                    endif

                    call compute_residual (r, rhs, phi, geo, bc, homog)
                    norm(0) = pnorm(r, r%valid, norm_type)
                    ! rho(0) = zero
                    rho(1) = zero
                    rho(2) = zero
                    rho(3) = zero
                    alpha(0) = zero
                    beta(0) = zero
                    omega(0) = zero

                    r_tilde%data = r%data
                    e%data = zero
                    restarts = restarts + 1

                    if (verbosity .ge. 3) then
                        print*, 'Restart number ', restarts
                    endif

                    init = .true.
                endif
            endif
        enddo

        if (verbosity .ge. 3) then
            print*, i, ' iterations. final rel |res| = ', norm(0)
        endif

        phi%data = phi%data + e%data

        ! Free memory
        call undefine_box_data (e)
        call undefine_box_data (p_tilde)
        call undefine_box_data (s_tilde)
        call undefine_box_data (r)
        call undefine_box_data (r_tilde)
        call undefine_box_data (p)
        call undefine_box_data (t)
        call undefine_box_data (v)
    end subroutine solve_bicgstab2

end module Poisson2D


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module MGPoisson2D
    use Poisson2D
    implicit none

    save

contains

    ! --------------------------------------------------------------------------
    ! This will overshoot the true maxDepth.
    ! --------------------------------------------------------------------------
    subroutine estimate_maxdepth (maxdepth, valid)
        integer, intent(inout) :: maxdepth
        type(box), intent(in)  :: valid

        integer                :: mgnx, mgny, d

        ! Prepare for first iter
        d = 0
        mgnx = valid%nx
        mgny = valid%ny

        do d = 0, maxdepth-1
            ! Can we coarsen this grid in any direction?
            if (mod(mgnx,2) .eq. 0) then
                mgnx = mgnx / 2
            else if (mod(mgny,2) .eq. 0) then
                mgny = mgny / 2
            else
                exit
            endif
        enddo

        maxdepth = min(maxdepth, d+1)

    end subroutine estimate_maxdepth


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine create_mgschedule (refx, refy, maxdepth, valid)
        integer, intent(inout)                        :: maxdepth
        integer, intent(out), dimension(0:maxdepth-1) :: refx, refy
        type(box), intent(in)                         :: valid

        integer                                       :: mgnx, mgny, d
        real(dp)                                      :: mgdx, mgdy, max_dx
        integer, parameter                            :: min_nx = 4
        integer, parameter                            :: min_ny = 4

        ! Initialize to bogus values
        refx = 0
        refy = 0

        ! Prepare for first iter
        d = 0
        mgnx = valid%nx
        mgny = valid%ny
        mgdx = valid%dx
        mgdy = valid%dy

        do d = 0, maxdepth-1
            ! Compute coarsening factors
            if ((mod(mgnx,2) .ne. 0) .or. (mod(mgny,2) .ne. 0)) then
                exit
            else
                ! Which directions should we coarsen to promote isotropy?
                max_dx = max(mgdx, mgdy)

                if (mgdx .le. max_dx / 2) then
                    refx(d) = 2
                else
                    refx(d) = 1
                endif

                if (mgdy .le. max_dx / 2) then
                    refy(d) = 2
                else
                    refy(d) = 1
                endif

                ! Is anisotropic coarsening worth it?
                if (refx(d)*refy(d) .eq. 1) then
                    refx(d) = 2
                    refy(d) = 2
                endif

                ! Is this amount of coarsening allowed?
                mgnx = mgnx / refx(d)
                mgny = mgny / refy(d)
                if ((mgnx .lt. min_nx) .or. (mgny .lt. min_ny)) then
                    exit
                endif

                ! We have chosen our coarsening factors. Move on to the next depth.
                mgdx = mgdx * real(refx(d),8)
                mgdy = mgdy * real(refy(d),8)
            endif
        enddo

        maxdepth = d

    end subroutine create_mgschedule


    ! ------------------------------------------------------------------------------
    ! fine and crse can have NO ghosts!
    ! ------------------------------------------------------------------------------
    subroutine restrict (fine, crse)
        type(box_data), intent(in)    :: fine
        type(box_data), intent(inout) :: crse

        logical, parameter            :: full_weighting = .false.
        integer                       :: refx, refy
        integer                       :: fi, fj, ci, cj
        real(dp)                      :: scale


        if ((fine%offi .eq. BD_CELL) .and. (fine%offj .eq. BD_CELL)) then
            ! Cell-centered in all directions.
            refx = fine%valid%nx / crse%valid%nx
            refy = fine%valid%ny / crse%valid%ny
            scale = one / (refx*refy)

            crse%data = zero
            do fj = fine%valid%jlo, fine%valid%jhi
                cj = fj/refy
                do fi = fine%valid%ilo, fine%valid%ihi
                    ci = fi/refx
                    crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                enddo
            enddo

        else if ((fine%offi .eq. BD_NODE) .and. (fine%offj .eq. BD_CELL)) then
            ! Nodal in x
            refx = (fine%valid%nx-1) / (crse%valid%nx-1)
            refy = fine%valid%ny / crse%valid%ny
            scale = one / refy
            crse%data = zero

            if (full_weighting) then
                do fj = fine%valid%jlo, fine%valid%jhi
                    cj = fj/refy
                    ci = crse%valid%ilo
                    fi = ci*refx
                    crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                enddo

                do fj = fine%valid%jlo, fine%valid%jhi
                    cj = fj/refy
                    ci = crse%valid%ihi
                    fi = ci*refx
                    crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                enddo

                scale = fourth / refy
                do fj = fine%valid%jlo, fine%valid%jhi
                    cj = fj/refy
                    do ci = crse%valid%ilo+1, crse%valid%ihi-1
                        fi = ci*refx
                        crse%data(ci,cj) = crse%data(ci,cj) &
                                         + scale*(fine%data(fi-1,fj) + two*fine%data(fi,fj) + fine%data(fi+1,fj))
                    enddo
                enddo

            else
                do fj = fine%valid%jlo, fine%valid%jhi
                    cj = fj/refy
                    do ci = crse%valid%ilo, crse%valid%ihi
                        fi = ci*refx
                        crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                    enddo
                enddo
            endif

        else if ((fine%offi .eq. BD_CELL) .and. (fine%offj .eq. BD_NODE)) then
            ! Nodal in y
            refx = fine%valid%nx / crse%valid%nx
            refy = (fine%valid%ny-1) / (crse%valid%ny-1)
            scale = one / refx
            crse%data = zero

            if (full_weighting) then
                cj = crse%valid%jlo
                fj = cj*refy
                do fi = fine%valid%ilo, fine%valid%ihi
                    ci = fi/refx
                    crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                enddo

                cj = crse%valid%jhi
                fj = cj*refy
                do fi = fine%valid%ilo, fine%valid%ihi
                    ci = fi/refx
                    crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                enddo

                scale = fourth / refx
                do cj = crse%valid%jlo+1, crse%valid%jhi-1
                    fj = cj*refy
                    do fi = fine%valid%ilo, fine%valid%ihi
                        ci = fi/refx
                        crse%data(ci,cj) = crse%data(ci,cj) &
                                         + scale*(fine%data(fi,fj-1) + two*fine%data(fi,fj) + fine%data(fi,fj+1))
                    enddo
                enddo

            else
                do cj = crse%valid%jlo, crse%valid%jhi
                    fj = cj*refy
                    do fi = fine%valid%ilo, fine%valid%ihi
                        ci = fi/refx
                        crse%data(ci,cj) = crse%data(ci,cj) + scale*fine%data(fi,fj)
                    enddo
                enddo
            endif
        else
            ! Nodal in all directions. Just copy nodes that don't vanish.
            refx = (fine%valid%nx-1) / (crse%valid%nx-1)
            refy = (fine%valid%ny-1) / (crse%valid%ny-1)

            do cj = crse%valid%jlo, crse%valid%jhi
                fj = cj*refy
                do ci = crse%valid%ilo, crse%valid%ihi
                    fi = ci*refx
                    crse%data(ci,cj) = fine%data(fi,fj)
                enddo
            enddo
        endif

    end subroutine restrict


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine prolong (fine, crse, fine_geo, crse_geo, crse_bc, order)
        type(box_data), intent(inout) :: fine
        type(box_data), intent(inout) :: crse
        type(geo_data), intent(in)    :: fine_geo
        type(geo_data), intent(in)    :: crse_geo
        type(bdry_data), intent(in)   :: crse_bc
        integer, intent(in)           :: order

        logical, parameter            :: homog = .true.
        logical, parameter            :: do_neum = .true.

        integer                       :: refx, refy
        integer                       :: fi, fj, ci, cj
        real(dp)                      :: mx, my, ml, mr, mc, eta
        logical                       :: remove_avg, rax, ray
        real(dp)                      :: dxprod, vol, dvol, sum

        ! Should we remove the average?
        rax = (crse_bc%type_xlo .eq. BCTYPE_PERIODIC) .and. (crse_bc%type_xhi .eq. BCTYPE_PERIODIC)
        rax = rax .or. &
              (crse_bc%type_xlo .eq. BCTYPE_NEUM) .and. (crse_bc%type_xhi .eq. BCTYPE_NEUM)

        ray = (crse_bc%type_ylo .eq. BCTYPE_PERIODIC) .and. (crse_bc%type_yhi .eq. BCTYPE_PERIODIC)
        ray = ray .or. &
              (crse_bc%type_ylo .eq. BCTYPE_NEUM) .and. (crse_bc%type_yhi .eq. BCTYPE_NEUM)

        remove_avg = rax .and. ray

        ! Make sure the geos weren't passed in reverse order.
        dxprod = fine_geo%dx * fine_geo%dy
        if (dxprod .gt. (crse_geo%dx * crse_geo%dy)) then
            print*, 'prolong: crse and fine geos passed in wrong order'
            stop
        endif

        if ((fine%offi .eq. BD_CELL) .and. (fine%offj .eq. BD_CELL)) then
            ! Cell-centered in all directions.
            refx = fine%valid%nx / crse%valid%nx
            refy = fine%valid%ny / crse%valid%ny

            if ((refx .eq. 2) .and. (refy .eq. 2)) then
                ! Constant interp
                do cj = crse%valid%jlo, crse%valid%jhi
                    fj = refy*cj
                    do ci = crse%valid%ilo, crse%valid%ihi
                        fi = refx*ci

                        fine%data(fi  ,fj  ) = fine%data(fi  ,fj  ) + crse%data(ci,cj)
                        fine%data(fi+1,fj  ) = fine%data(fi+1,fj  ) + crse%data(ci,cj)
                        fine%data(fi  ,fj+1) = fine%data(fi  ,fj+1) + crse%data(ci,cj)
                        fine%data(fi+1,fj+1) = fine%data(fi+1,fj+1) + crse%data(ci,cj)
                    enddo
                enddo

                if (order .ge. 1) then
                    ! Upgrade to linear interp
                    call fill_ghosts (crse, crse_bc, crse_geo, homog, do_neum)
                    do cj = crse%valid%jlo, crse%valid%jhi
                        fj = refy*cj
                        do ci = crse%valid%ilo, crse%valid%ihi
                            fi = refx*ci

                            mx = half * (crse%data(ci+1,cj) - crse%data(ci-1,cj))
                            my = half * (crse%data(ci,cj+1) - crse%data(ci,cj-1))

                            fine%data(fi  ,fj  ) = fine%data(fi  ,fj  ) + fourth*(-mx - my)
                            fine%data(fi+1,fj  ) = fine%data(fi+1,fj  ) + fourth*( mx - my)
                            fine%data(fi  ,fj+1) = fine%data(fi  ,fj+1) + fourth*(-mx + my)
                            fine%data(fi+1,fj+1) = fine%data(fi+1,fj+1) + fourth*( mx + my)
                        enddo
                    enddo
                endif

                if (order .ge. 2) then
                    ! Upgrade to quadratic interp
                    do cj = crse%valid%jlo, crse%valid%jhi
                        fj = refy*cj
                        do ci = crse%valid%ilo, crse%valid%ihi
                            fi = refx*ci

                            mx = fourth * (crse%data(ci+1,cj) - two*crse%data(ci,cj) + crse%data(ci-1,cj))
                            my = fourth * (crse%data(ci,cj+1) - two*crse%data(ci,cj) + crse%data(ci,cj-1))

                            fine%data(fi  ,fj  ) = fine%data(fi  ,fj  ) + eighth*(mx + my)
                            fine%data(fi+1,fj  ) = fine%data(fi+1,fj  ) + eighth*(mx + my)
                            fine%data(fi  ,fj+1) = fine%data(fi  ,fj+1) + eighth*(mx + my)
                            fine%data(fi+1,fj+1) = fine%data(fi+1,fj+1) + eighth*(mx + my)
                        enddo
                    enddo

                    call extrapolate_ghosts (crse, 3)
                    do cj = crse%valid%jlo, crse%valid%jhi
                        fj = refy*cj
                        do ci = crse%valid%ilo, crse%valid%ihi
                            fi = refx*ci

                            mc = sixteenth * (  crse%data(ci+1,cj+1) &
                                              - crse%data(ci+1,cj-1) &
                                              - crse%data(ci-1,cj+1) &
                                              + crse%data(ci-1,cj-1)  )

                            fine%data(fi  ,fj  ) = fine%data(fi  ,fj  ) + fourth*mc
                            fine%data(fi+1,fj  ) = fine%data(fi+1,fj  ) - fourth*mc
                            fine%data(fi  ,fj+1) = fine%data(fi  ,fj+1) - fourth*mc
                            fine%data(fi+1,fj+1) = fine%data(fi+1,fj+1) + fourth*mc
                        enddo
                    enddo
                endif

                ! Remove average if necessary
                if (remove_avg) then
                    sum = zero
                    vol = zero
                    do fj = fine%valid%jlo, fine%valid%jhi
                        do fi = fine%valid%ilo, fine%valid%ihi
                            dvol = dxprod * fine_geo%J%data(fi,fj)
                            sum = sum + dvol * fine%data(fi,fj)
                            vol = vol + dvol
                        enddo
                    enddo
                    fine%data = fine%data - sum/vol
                endif

            else
                print*, "prolong: Needs updatin'"
                stop
            endif

        else
            print*, 'prolong: Cannot handle nodal data yet.'
            stop
        endif
    end subroutine prolong


    ! ------------------------------------------------------------------------------
    ! rhs will be scaled with J, but restored before the function exits.
    ! ------------------------------------------------------------------------------
    subroutine vcycle (phi, rhs, geo, bc, homog, amrrefx, amrrefy, &
                       tol, maxiters, maxdepth_user, numcycles, &
                       smooth_down, smooth_up, smooth_bottom, &
                       zerophi, verbosity)
        type(box_data), intent(inout) :: phi
        type(box_data), intent(inout) :: rhs
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog
        integer, intent(in)           :: amrrefx, amrrefy
        real(dp), intent(in)          :: tol
        integer, intent(in)           :: maxiters
        integer, intent(in)           :: maxdepth_user
        integer, intent(in)           :: numcycles
        integer, intent(in)           :: smooth_down, smooth_up, smooth_bottom
        logical, intent(in)           :: zerophi
        integer, intent(in)           :: verbosity

        type(box)                                  :: valid
        integer                                    :: maxdepth
        integer, dimension(:), allocatable         :: refx, refy
        integer                                    :: ierr, d, iter
        real(dp)                                   :: mgdx, mgdy
        real(dp), dimension(:), allocatable        :: relres
        type(box_data), dimension(:), allocatable  :: e, r, invdiags, work1
        type(geo_data), dimension(:), allocatable  :: mggeo
        type(bdry_data), dimension(:), allocatable :: mgbc
        integer                                    :: ilo, ihi, jlo, jhi
        real(dp)                                   :: rscale, sum

        ! Estimate size of scheduling vectors
        valid = rhs%valid
        if (maxdepth_user .lt. 0) then
            maxdepth = 1000
        else
            maxdepth = maxdepth_user
        endif
        call estimate_maxdepth (maxdepth, valid)

        ! Allocate schedule vectors
        allocate (refx(0:maxdepth-1), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'vcycle: Out of memory'
            stop
        endif

        allocate (refy(0:maxdepth-1), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'vcycle: Out of memory'
            stop
        endif

        ! Create coarsening schedule
        call create_mgschedule(refx, refy, maxdepth, valid)
        if (verbosity .ge. 1) then
            print*, 'max MG depth = ', maxdepth
        endif
        if (verbosity .ge. 2) then
            mgdx = valid%dx
            mgdy = valid%dy
            do d = 0, maxdepth-1
                ! print*, 'ref(', d, ') = ', refx(d), ', ', refy(d)
                print*, d, ': dx = ', mgdx, ', dy = ', mgdy
                mgdx = mgdx * real(refx(d),8)
                mgdy = mgdy * real(refy(d),8)
            enddo
            print*, d, ': dx = ', mgdx, ', dy = ', mgdy
        endif

        ! Allocate and set up workspace
        allocate (relres(0:maxiters), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (e(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (r(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (invdiags(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (mggeo(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (work1(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (mgbc(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        do d = 0, maxdepth
            ! print*, '*** valid = ', valid

            ! All BCs will be homogeneous after first residual calculation.
            call define_bdry_data (mgbc(d), valid, &
                                   bc%type_xlo, bc%type_xhi, bc%type_ylo, bc%type_yhi, &
                                   BCMODE_UNIFORM, BCMODE_UNIFORM, BCMODE_UNIFORM, BCMODE_UNIFORM)
            mgbc(d)%data_xlo = zero
            mgbc(d)%data_xhi = zero
            mgbc(d)%data_ylo = zero
            mgbc(d)%data_yhi = zero

            call define_box_data (e(d), valid, 1, 1, BD_CELL, BD_CELL)
            call define_box_data (r(d), valid, 0, 0, BD_CELL, BD_CELL)
            call define_box_data (work1(d), valid, 0, 0, BD_CELL, BD_CELL)

            call define_box_data (mggeo(d)%J, valid, geo%J%ngx, geo%J%ngy, BD_CELL, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_xx, valid, geo%Jgup_xx%ngx, geo%Jgup_xx%ngy, BD_NODE, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_xy, valid, geo%Jgup_xy%ngx, geo%Jgup_xy%ngy, BD_NODE, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_yx, valid, geo%Jgup_yx%ngx, geo%Jgup_yx%ngy, BD_CELL, BD_NODE)
            call define_box_data (mggeo(d)%Jgup_yy, valid, geo%Jgup_yy%ngx, geo%Jgup_yy%ngy, BD_CELL, BD_NODE)
            call define_box_data (invdiags(d), valid, 0, 0, BD_CELL, BD_CELL)

            if (d .eq. 0) then
                ! For now, just copy the data. It would be more economical to
                ! point to the data that already exists.
                mggeo(d)%dx = valid%dx
                mggeo(d)%dy = valid%dy
                mggeo(d)%J%data = geo%J%data
                mggeo(d)%Jgup_xx%data = geo%Jgup_xx%data
                mggeo(d)%Jgup_xy%data = geo%Jgup_xy%data
                mggeo(d)%Jgup_yx%data = geo%Jgup_yx%data
                mggeo(d)%Jgup_yy%data = geo%Jgup_yy%data
                call compute_invdiags (invdiags(d), mggeo(d), mgbc(d), .true.)
            else
                ! Coarsen the data from MG level d-1
                mggeo(d)%dx = valid%dx
                mggeo(d)%dy = valid%dy
                call restrict (mggeo(d-1)%J, mggeo(d)%J)
                call restrict (mggeo(d-1)%Jgup_xx, mggeo(d)%Jgup_xx)
                call restrict (mggeo(d-1)%Jgup_xy, mggeo(d)%Jgup_xy)
                call restrict (mggeo(d-1)%Jgup_yx, mggeo(d)%Jgup_yx)
                call restrict (mggeo(d-1)%Jgup_yy, mggeo(d)%Jgup_yy)
                call compute_invdiags (invdiags(d), mggeo(d), mgbc(d), .true.)
            endif

            ! Move on to next depth
            if (d .lt. maxDepth) then
                call coarsen_box (valid, refx(d), refy(d))
            endif
        enddo

        ! Initialize phi to zero if necessary
        if (zerophi) then
            phi%data = zero
        endif

        ! Scale rhs by J
        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi
        ! rhs%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) &
        !                           * geo%J%data(ilo:ihi,jlo:jhi)


        ! Set up residual equation
        call compute_residual (r(0), rhs, phi, geo, bc, homog)
        rscale = pnorm (r(0), r(0)%valid, 2)

        relres = zero
        if (rscale .eq. zero) then
            relres(0) = rscale
            rscale = one
        else
            relres(0) = one
        endif

        if (verbosity .ge. 1) then
            sum = integrate2d (r(0), r(0)%valid, mggeo(0), .false.)
            print*, 'scale |res| = ', rscale, ', sum res = ', sum
            print*, 'V-Cycle iter ', 0, ': rel |res| = ', relres(0)
        endif

        ! The main iteration loop.
        do iter = 1, maxiters
            ! Set the initial guess
            ! e(0)%data = zero
            e(0)%data(ilo:ihi,jlo:jhi) = r(0)%data(ilo:ihi,jlo:jhi) * invdiags(0)%data(ilo:ihi,jlo:jhi)

            ! Solve for phi's correction.
            call vcycle_noinit (e, r, invdiags, work1, mggeo, mgbc, &
                                refx, refy, &
                                tol, maxdepth, 0, &
                                numcycles, &
                                smooth_down, smooth_up, smooth_bottom, &
                                verbosity)

            ! Apply correction
            phi%data = phi%data + e(0)%data

            ! Compute residual
            call compute_residual (r(0), rhs, phi, mggeo(0), bc, homog)
            relres(iter) = pnorm (r(0), r(0)%valid, 2) / rscale
            if (verbosity .ge. 1) then
                print*, 'V-Cycle iter ', iter, ': rel |res| = ', relres(iter)
            endif

            ! Did we converge?
            if (relres(iter) .le. tol) then
                if (verbosity .ge. 1) then
                    print*, "Converged."
                endif

                exit
            endif

            ! Are we diverging?
            if (relres(iter) .gt. relres(iter-1)) then
                if (verbosity .ge. 1) then
                    print*, 'Diverging.'
                endif

                ! Undo last correction
                phi%data = phi%data - e(0)%data

                exit
            endif
        enddo

        ! ! Restore rhs scaling.
        ! rhs%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) &
        !                           / geo%J%data(ilo:ihi,jlo:jhi)

        ! Free memory
        do d = 0, maxDepth
            call undefine_bdry_data (mgbc(d))
            call undefine_box_data (invdiags(d))
            call undefine_box_data (mggeo(d)%J)
            call undefine_box_data (mggeo(d)%Jgup_xx)
            call undefine_box_data (mggeo(d)%Jgup_xy)
            call undefine_box_data (mggeo(d)%Jgup_yx)
            call undefine_box_data (mggeo(d)%Jgup_yy)
            call undefine_box_data (work1(d))
            call undefine_box_data (r(d))
            call undefine_box_data (e(d))
        enddo

        if (allocated(mgbc)) deallocate(mgbc)
        if (allocated(invDiags)) deallocate(invDiags)
        if (allocated(mggeo)) deallocate(mggeo)
        if (allocated(work1)) deallocate(work1)
        if (allocated(r)) deallocate(r)
        if (allocated(e)) deallocate(e)
        if (allocated(relres)) deallocate(relres)
        if (allocated(refx)) deallocate(refx)
        if (allocated(refy)) deallocate(refy)

    end subroutine vcycle


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    recursive subroutine vcycle_noinit (e, r, invdiags, work1, geo, bc, &
                                        refx, refy, &
                                        tol, maxdepth, depth, &
                                        numcycles, &
                                        smooth_down, smooth_up, smooth_bottom, &
                                        verbosity)
        integer, intent(in)                                  :: maxdepth
        type(box_data), intent(inout), dimension(0:maxdepth) :: e, r, work1
        type(box_data), intent(in), dimension(0:maxdepth)    :: invdiags
        type(geo_data), intent(in), dimension(0:maxdepth)    :: geo
        type(bdry_data), intent(in), dimension(0:maxdepth)   :: bc
        integer, intent(in), dimension(0:maxdepth-1)         :: refx, refy
        real(dp), intent(in)                                 :: tol
        integer, intent(in)                                  :: depth
        integer, intent(in)                                  :: smooth_down, smooth_up, smooth_bottom
        integer, intent(in)                                  :: numcycles
        integer, intent(in)                                  :: verbosity

        character*2                                          :: indent = '  '
        logical, parameter                                   :: homog = .true.
        integer                                              :: ilo, ihi, jlo, jhi
        integer                                              :: cilo, cihi, cjlo, cjhi
        integer                                              :: curcycle
        real(dp)                                             :: norm, sum
        type(box_data)                                       :: tmp

        integer, parameter                                   :: prolong_order = 1

        ! Relaxation params
        real(dp), parameter                                  :: relax_tol = -one
        real(dp), parameter                                  :: relax_omega = one !-third ! Was 1.33
        ! logical, parameter                                   :: relax_redblack = .true.
        integer, parameter                                   :: relax_verbosity = 0

        ! Bottom solver params
        real(dp), parameter                                  :: bottom_tol = 1.0E-6_dp
        integer, parameter                                   :: bottom_maxiters = 80
        integer, parameter                                   :: bottom_maxrestarts = 5
        integer, parameter                                   :: bottom_verbosity = 0

        if (verbosity .ge. 7) then
            print*, repeat(indent,depth), 'MG depth = ', depth
        endif

        ! Compute the rhs integral
        if (verbosity .ge. 8) then
            call define_box_data (tmp, r(depth));

            norm = pnorm (r(depth), r(depth)%valid, 2)
            sum = integrate2d (r(depth), r(depth)%valid, geo(depth), .false.)
            print*, repeat(indent,depth), 'sq |rhs| = ', norm, ', sum rhs = ', sum
        endif

        ilo = r(depth)%valid%ilo
        ihi = r(depth)%valid%ihi
        jlo = r(depth)%valid%jlo
        jhi = r(depth)%valid%jhi

        if (depth .eq. maxdepth) then
            ! Use bottom solver...

            ! --- Bottom relaxation ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Bottom relax'
            endif
            call relax_gs (e(depth), r(depth), geo(depth), bc(depth), homog, &
                           invdiags(depth), relax_omega, relax_tol, &
                           smooth_bottom, &
                           .false., & ! zero phi?
                           relax_verbosity)

            ! Diagnostics
            if (verbosity .ge. 8) then
                call compute_residual (tmp, r(depth), e(depth), geo(depth), bc(depth), homog)
                norm = pnorm (tmp, tmp%valid, 2)
                print*, repeat(indent,depth), 'sq |rhs| = ', norm
            endif

            ! --- Bottom solver ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Bottom solver'
            endif
            call solve_bicgstab2 (e(depth), r(depth), geo(depth), bc(depth), homog, invdiags(depth), &
                                  bottom_tol, bottom_maxiters, bottom_maxrestarts, &
                                  .false., & ! zero phi?
                                  bottom_verbosity)

            ! Diagnostics
            if (verbosity .ge. 8) then
                call compute_residual (tmp, r(depth), e(depth), geo(depth), bc(depth), homog)
                norm = pnorm (tmp, tmp%valid, 2)
                print*, repeat(indent,depth), 'sq |rhs| = ', norm
            endif

        else
            ! V-Cycle...

            ! --- Downward relaxation ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Smooth down'
            endif
            call relax_gs (e(depth), r(depth), geo(depth), bc(depth), homog, &
                           invdiags(depth), relax_omega, relax_tol, &
                           smooth_down, &
                           .false., & ! zero phi?
                           relax_verbosity)

            ! Diagnostics
            if (verbosity .ge. 8) then
                call compute_residual (tmp, r(depth), e(depth), geo(depth), bc(depth), homog)
                norm = pnorm (tmp, tmp%valid, 2)
                print*, repeat(indent,depth), 'sq |rhs| = ', norm
            endif

            ! --- Restrict residual ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Restrict resudual'
            endif
            call compute_residual (work1(depth), r(depth), e(depth), &
                                   geo(depth), bc(depth), homog)
            call restrict (work1(depth), r(depth+1))


            ! --- Coarse level solve ---
            ! Set the initial guess
            cilo = r(depth+1)%valid%ilo
            cihi = r(depth+1)%valid%ihi
            cjlo = r(depth+1)%valid%jlo
            cjhi = r(depth+1)%valid%jhi
            ! e(depth)%data = zero
            e(depth+1)%data(cilo:cihi,cjlo:cjhi) = r(depth+1)%data(cilo:cihi,cjlo:cjhi) &
                                                 * invdiags(depth+1)%data(cilo:cihi,cjlo:cjhi)

            ! Solve
            do curcycle = 1, numcycles
                call vcycle_noinit (e, r, invdiags, work1, geo, bc, &
                                    refx, refy, &
                                    tol, maxdepth, depth+1, &
                                    numcycles, &
                                    smooth_down, smooth_up, smooth_bottom, &
                                    verbosity)
            enddo

            ! --- Prolong correction ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Prolong and add correction'
            endif
            call prolong (e(depth), e(depth+1), geo(depth), geo(depth+1), bc(depth+1), prolong_order)

            ! Diagnostics
            if (verbosity .ge. 8) then
                call compute_residual (tmp, r(depth), e(depth), geo(depth), bc(depth), homog)
                norm = pnorm (tmp, tmp%valid, 2)
                print*, repeat(indent,depth), 'sq |rhs| = ', norm
            endif

            ! --- Upward relaxation ---
            if (verbosity .ge. 7) then
                print*, repeat(indent,depth), 'Smooth up'
            endif
            call relax_gs (e(depth), r(depth), geo(depth), bc(depth), homog, &
                           invdiags(depth), relax_omega, relax_tol, &
                           smooth_up, &
                           .false., & ! zero phi?
                           relax_verbosity)

            ! Diagnostics
            if (verbosity .ge. 8) then
                call compute_residual (tmp, r(depth), e(depth), geo(depth), bc(depth), homog)
                norm = pnorm (tmp, tmp%valid, 2)
                print*, repeat(indent,depth), 'sq |rhs| = ', norm
            endif
        endif

        ! Free memory
        if (verbosity .ge. 8) then
            call undefine_box_data (tmp);
        endif

    end subroutine vcycle_noinit


    ! ------------------------------------------------------------------------------
    ! rhs will be scaled with J, but restored before the function exits.
    ! ------------------------------------------------------------------------------
    subroutine fmg (phi, rhs, geo, bc, homog, amrrefx, amrrefy, &
                    tol, maxiters, maxdepth_user, numcycles, &
                    smooth_down, smooth_up, smooth_bottom, &
                    zerophi, verbosity)
        type(box_data), intent(inout) :: phi
        type(box_data), intent(inout) :: rhs
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog
        integer, intent(in)           :: amrrefx, amrrefy
        real(dp), intent(in)          :: tol
        integer, intent(in)           :: maxiters
        integer, intent(in)           :: maxdepth_user
        integer, intent(in)           :: numcycles
        integer, intent(in)           :: smooth_down, smooth_up, smooth_bottom
        logical, intent(in)           :: zerophi
        integer, intent(in)           :: verbosity

        type(box)                                  :: valid
        integer                                    :: maxdepth
        integer, dimension(:), allocatable         :: refx, refy
        integer                                    :: ierr, d, iter
        real(dp)                                   :: mgdx, mgdy
        real(dp), dimension(:), allocatable        :: relres
        type(box_data), dimension(:), allocatable  :: e, r, invdiags, work1
        type(geo_data), dimension(:), allocatable  :: mggeo
        type(bdry_data), dimension(:), allocatable :: mgbc
        integer                                    :: ilo, ihi, jlo, jhi
        real(dp)                                   :: rscale, sum

        ! Estimate size of scheduling vectors
        valid = rhs%valid
        if (maxdepth_user .lt. 0) then
            maxdepth = 1000
        else
            maxdepth = maxdepth_user
        endif
        call estimate_maxdepth (maxdepth, valid)

        ! Allocate schedule vectors
        allocate (refx(0:maxdepth-1), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'vcycle: Out of memory'
            stop
        endif

        allocate (refy(0:maxdepth-1), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'vcycle: Out of memory'
            stop
        endif

        ! Create coarsening schedule
        call create_mgschedule(refx, refy, maxdepth, valid)
        if (verbosity .ge. 1) then
            print*, 'max MG depth = ', maxdepth
        endif
        ! if (verbosity .ge. 2) then
        !     mgdx = valid%dx
        !     mgdy = valid%dy
        !     do d = 0, maxdepth-1
        !         ! print*, 'ref(', d, ') = ', refx(d), ', ', refy(d)
        !         print*, d, ': dx = ', mgdx, ', dy = ', mgdy
        !         mgdx = mgdx * real(refx(d),8)
        !         mgdy = mgdy * real(refy(d),8)
        !     enddo
        !     print*, d, ': dx = ', mgdx, ', dy = ', mgdy
        ! endif

        ! Allocate and set up workspace
        allocate (relres(0:maxiters), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (e(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (r(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (invdiags(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (mggeo(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (work1(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        allocate (mgbc(0:maxdepth), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'Out of memory'
            stop
        endif

        do d = 0, maxdepth
            ! print*, '*** valid = ', valid

            ! All BCs will be homogeneous after first residual calculation.
            call define_bdry_data (mgbc(d), valid, &
                                   bc%type_xlo, bc%type_xhi, bc%type_ylo, bc%type_yhi, &
                                   BCMODE_UNIFORM, BCMODE_UNIFORM, BCMODE_UNIFORM, BCMODE_UNIFORM)
            mgbc(d)%data_xlo = zero
            mgbc(d)%data_xhi = zero
            mgbc(d)%data_ylo = zero
            mgbc(d)%data_yhi = zero

            call define_box_data (e(d), valid, 1, 1, BD_CELL, BD_CELL)
            call define_box_data (r(d), valid, 0, 0, BD_CELL, BD_CELL)
            call define_box_data (work1(d), valid, 0, 0, BD_CELL, BD_CELL)

            call define_box_data (mggeo(d)%J, valid, geo%J%ngx, geo%J%ngy, BD_CELL, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_xx, valid, geo%Jgup_xx%ngx, geo%Jgup_xx%ngy, BD_NODE, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_xy, valid, geo%Jgup_xy%ngx, geo%Jgup_xy%ngy, BD_NODE, BD_CELL)
            call define_box_data (mggeo(d)%Jgup_yx, valid, geo%Jgup_yx%ngx, geo%Jgup_yx%ngy, BD_CELL, BD_NODE)
            call define_box_data (mggeo(d)%Jgup_yy, valid, geo%Jgup_yy%ngx, geo%Jgup_yy%ngy, BD_CELL, BD_NODE)
            call define_box_data (invdiags(d), valid, 0, 0, BD_CELL, BD_CELL)

            if (d .eq. 0) then
                ! For now, just copy the data. It would be more economical to
                ! point to the data that already exists.
                mggeo(d)%dx = valid%dx
                mggeo(d)%dy = valid%dy
                mggeo(d)%J%data = geo%J%data
                mggeo(d)%Jgup_xx%data = geo%Jgup_xx%data
                mggeo(d)%Jgup_xy%data = geo%Jgup_xy%data
                mggeo(d)%Jgup_yx%data = geo%Jgup_yx%data
                mggeo(d)%Jgup_yy%data = geo%Jgup_yy%data
                call compute_invdiags (invdiags(d), mggeo(d), mgbc(d), .true.)
            else
                ! Coarsen the data from MG level d-1
                mggeo(d)%dx = valid%dx
                mggeo(d)%dy = valid%dy
                call restrict (mggeo(d-1)%J, mggeo(d)%J)
                call restrict (mggeo(d-1)%Jgup_xx, mggeo(d)%Jgup_xx)
                call restrict (mggeo(d-1)%Jgup_xy, mggeo(d)%Jgup_xy)
                call restrict (mggeo(d-1)%Jgup_yx, mggeo(d)%Jgup_yx)
                call restrict (mggeo(d-1)%Jgup_yy, mggeo(d)%Jgup_yy)
                call compute_invdiags (invdiags(d), mggeo(d), mgbc(d), .true.)
            endif

            ! Move on to next depth
            if (d .lt. maxDepth) then
                call coarsen_box (valid, refx(d), refy(d))
            endif
        enddo

        ! Initialize phi to zero if necessary
        if (zerophi) then
            phi%data = zero
        endif

        ! Scale rhs by J
        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi
        ! rhs%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) &
        !                           * geo%J%data(ilo:ihi,jlo:jhi)


        ! Set up residual equation
        call compute_residual (r(0), rhs, phi, geo, bc, homog)
        rscale = pnorm (r(0), r(0)%valid, 2)

        relres = zero
        if (rscale .eq. zero) then
            relres(0) = rscale
            rscale = one
        else
            relres(0) = one
        endif

        if (verbosity .ge. 1) then
            sum = integrate2d (r(0), r(0)%valid, mggeo(0), .false.)
            print*, 'scale |res| = ', rscale, ', sum res = ', sum
            print*, 'FMG-Cycle iter ', 0, ': rel |res| = ', relres(0)
        endif

        do iter = 1, maxiters
            ! Solve!
            call fmg_noinit (e, r, invdiags, work1, mggeo, mgbc, &
                             refx, refy, &
                             tol, maxdepth, 0, &
                             numcycles, &
                             smooth_down, smooth_up, smooth_bottom, &
                             verbosity)

            ! Add correction to phi
            phi%data = phi%data + e(0)%data

            ! Compute new residual
            call compute_residual (r(0), rhs, phi, mggeo(0), bc, homog)
            relres(iter) = pnorm (r(0), r(0)%valid, 2) / rscale

            if (verbosity .ge. 1) then
                print*, 'FMG-Cycle iter ', iter, ': rel |res| = ', relres(iter)
            endif

            ! Did we converge?
            if (relres(iter) .le. tol) then
                if (verbosity .ge. 1) then
                    print*, "Converged."
                endif

                exit
            endif

            ! Are we diverging?
            if (relres(iter) .gt. relres(iter-1)) then
                if (verbosity .ge. 1) then
                    print*, 'Diverging.'
                endif

                ! Undo last correction
                phi%data = phi%data - e(0)%data

                exit
            endif
        enddo

        ! Free memory
        do d = 0, maxDepth
            call undefine_bdry_data (mgbc(d))
            call undefine_box_data (invdiags(d))
            call undefine_box_data (mggeo(d)%J)
            call undefine_box_data (mggeo(d)%Jgup_xx)
            call undefine_box_data (mggeo(d)%Jgup_xy)
            call undefine_box_data (mggeo(d)%Jgup_yx)
            call undefine_box_data (mggeo(d)%Jgup_yy)
            call undefine_box_data (work1(d))
            call undefine_box_data (r(d))
            call undefine_box_data (e(d))
        enddo

        if (allocated(mgbc)) deallocate(mgbc)
        if (allocated(invDiags)) deallocate(invDiags)
        if (allocated(mggeo)) deallocate(mggeo)
        if (allocated(work1)) deallocate(work1)
        if (allocated(r)) deallocate(r)
        if (allocated(e)) deallocate(e)
        if (allocated(relres)) deallocate(relres)
        if (allocated(refx)) deallocate(refx)
        if (allocated(refy)) deallocate(refy)

    end subroutine fmg


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    recursive subroutine fmg_noinit (e, r, invdiags, work1, geo, bc, &
                                     refx, refy, &
                                     tol, maxdepth, depth, &
                                     numcycles, &
                                     smooth_down, smooth_up, smooth_bottom, &
                                     verbosity)
        integer, intent(in)                                  :: maxdepth
        type(box_data), intent(inout), dimension(0:maxdepth) :: e, r, work1
        type(box_data), intent(in), dimension(0:maxdepth)    :: invdiags
        type(geo_data), intent(in), dimension(0:maxdepth)    :: geo
        type(bdry_data), intent(in), dimension(0:maxdepth)   :: bc
        integer, intent(in), dimension(0:maxdepth-1)         :: refx, refy
        real(dp), intent(in)                                 :: tol
        integer, intent(in)                                  :: depth
        integer, intent(in)                                  :: smooth_down, smooth_up, smooth_bottom
        integer, intent(in)                                  :: numcycles
        integer, intent(in)                                  :: verbosity

        integer                                              :: ilo, ihi, jlo, jhi

        character*2                                          :: indent = '  '
        integer                                              :: curcycle
        integer, parameter                                   :: prolong_order = 2

        if (verbosity .ge. 7) then
            print*, repeat(indent,depth), 'MG depth = ', depth
        endif

        ! Initialize solution
        e(depth)%data = zero
        ! ilo = r(depth)%valid%ilo
        ! ihi = r(depth)%valid%ihi
        ! jlo = r(depth)%valid%jlo
        ! jhi = r(depth)%valid%jhi
        ! e(depth)%data(ilo:ihi,jlo:jhi) = r(depth)%data(ilo:ihi,jlo:jhi) * invdiags(depth)%data(ilo:ihi,jlo:jhi)

        ! Go to the coarser MG level if possible.
        if (depth .lt. maxdepth) then
            ! Coarsen the residual
            call restrict (r(depth), r(depth+1))

            ! Solve
            call fmg_noinit (e, r, invdiags, work1, geo, bc, &
                             refx, refy, &
                             tol, maxdepth, depth+1, &
                             numcycles, &
                             smooth_down, smooth_up, smooth_bottom, &
                             verbosity)

            ! Refine solution
            call prolong (e(depth), e(depth+1), geo(depth), geo(depth+1), bc(depth+1), prolong_order)
        endif

        ! Solve with a V-Cycle.
        do curcycle = 1, numcycles
            call vcycle_noinit (e, r, invdiags, work1, geo, bc, &
                                refx, refy, &
                                tol, maxdepth, depth, &
                                numcycles, &
                                smooth_down, smooth_up, smooth_bottom, &
                                verbosity)
        enddo
    end subroutine fmg_noinit


end module MGPoisson2D


! NOTES ............................
! 1. inner_prod should scale by J
! 2. check if restrict or prolong need J scaling too.
! 3. Test compute_pd
! 4. Refinement ratio is hard-coded at 4,4 in fill_ghosts.
! 5. If rscale = zero, manually set to one.
! 6. Handle beta = zero in bicgstab.
