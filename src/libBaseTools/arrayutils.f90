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
