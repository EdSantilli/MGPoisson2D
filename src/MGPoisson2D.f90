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
! Defines the precision of all numbers that we will use.
! ------------------------------------------------------------------------------
module precision
    integer, parameter  :: dp = kind(0.d0)

    real(dp), parameter :: zero      =  0.0_dp
    real(dp), parameter :: one       =  1.0_dp
    real(dp), parameter :: two       =  2.0_dp
    real(dp), parameter :: three     =  3.0_dp
    real(dp), parameter :: four      =  4.0_dp
    real(dp), parameter :: five      =  5.0_dp
    real(dp), parameter :: six       =  6.0_dp
    real(dp), parameter :: seven     =  7.0_dp
    real(dp), parameter :: eight     =  8.0_dp
    real(dp), parameter :: nine      =  9.0_dp
    real(dp), parameter :: ten       = 10.0_dp
    real(dp), parameter :: eleven    = 11.0_dp
    real(dp), parameter :: twelve    = 12.0_dp
    real(dp), parameter :: thirteen  = 13.0_dp
    real(dp), parameter :: fourteen  = 14.0_dp
    real(dp), parameter :: fifteen   = 15.0_dp
    real(dp), parameter :: sixteen   = 16.0_dp
    real(dp), parameter :: seventeen = 17.0_dp
    real(dp), parameter :: eighteen  = 18.0_dp
    real(dp), parameter :: nineteen  = 19.0_dp
    real(dp), parameter :: twenty    = 20.0_dp

    real(dp), parameter :: half      = one / two
    real(dp), parameter :: third     = one / three
    real(dp), parameter :: fourth    = one / four

    real(dp), parameter :: threehalves  = three / two
    real(dp), parameter :: threefourths = three / four

    real(dp), parameter :: pi        = four * atan(one)
    real(dp), parameter :: bogus_val = 1.2345E300_dp

end module precision


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


    ! --------------------------------------------------------------------------
    ! These are the BC modes.
    ! BCMODE_UNIFORM = the entire boundary uses the same BC value. For example,
    ! a homogeneous BC would use a uniform value of 0.
    ! BCMODE_NONUNIFORM = the BC values vary over the bounding face. For
    ! example, a value of sin(k*y) along an x boundary face is nonuniform.
    ! --------------------------------------------------------------------------
    integer, parameter :: BCMODE_UNIFORM = 0
    integer, parameter :: BCMODE_NONUNIFORM = 1


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


    ! --------------------------------------------------------------------------
    ! Computes inner product of data arrays for norms.
    ! --------------------------------------------------------------------------
    interface inner_prod
        module procedure inner_prod_array
        module procedure inner_prod_box_data
    end interface

contains

    ! --------------------------------------------------------------------------
    ! Sends a box_data to gnuplot
    ! --------------------------------------------------------------------------
    subroutine plot (x, y, data)
        type(box_data), intent(in) :: x, y, data

        ! integer :: plot_type                                 ! 1 for linear plot, 2 for log plot, 3 for log-log plot
        ! character(len=*) :: xlabel,ylabel,title1,title2      ! plot axis labels and title
        ! !---------------
        ! integer :: i
        ! integer :: ret
        ! !---------------

        ! ! write data on two separate files
        ! open (10, access='SEQUENTIAL', file='xydata1.dat')
        ! do i = 1,n1
        !    write(10,*) xydata1(i,1),xydata1(i,2)
        ! enddo
        ! CLOSE(10,STATUS='KEEP')

        ! OPEN(10,ACCESS='SEQUENTIAL',FILE='xydata2.dat')
        ! do i=1,n2
        !    write(10,*) xydata2(i,1),xydata2(i,2)
        ! enddo
        ! CLOSE(10,STATUS='KEEP')

        ! ! create gnuplot command file
        ! OPEN(10,ACCESS='SEQUENTIAL',FILE='gp.txt')
        ! write(10,*) 'set terminal postscript'
        ! write(10,*) 'set output "plot.ps"'
        ! write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
        ! write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
        ! if (plot_type==2) write(10,*) 'set log y'
        ! if (plot_type==3) then
        !    write(10,*) 'set log x'
        !    write(10,*) 'set log y'
        ! endif

        ! if(n1>0.AND.n2>0) then
        !     write(10,*) 'plot "xydata1.dat" using 1:2 with lines title "'//TRIM(title1)//&
        !                 &'", "xydata2.dat" using 1:2 with lines title "'//TRIM(title2)//'"'
        ! endif

        ! if(n1>0.AND.n2==0) write(10,*) 'plot "xydata1.dat" using 1:2 with lines title "'//TRIM(title1)//'"'

        ! if(n2>0.AND.n1==0) write(10,*) 'plot "xydata2.dat" using 1:2 with lines title "'//TRIM(title2)//'"'

        ! CLOSE(10,STATUS='KEEP')

        ! ! plot curve with gnuplot and cleanup files
        ! ret=SYSTEM('gnuplot gp.txt')
        ! ret=SYSTEM('rm gp.txt')
        ! ret=SYSTEM('rm xydata1.dat')
        ! ret=SYSTEM('rm xydata2.dat')

    end subroutine plot


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
        integer               :: ilo, ihi, jlo, jhi, offi, offj
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
        real(dp)  :: L, H

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

        call define_box_data (dest, cc_valid, src%ngx, src%ngy, src%offi, src%offj)
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

        allocate (bd%data (bd%bx%ilo : bd%bx%ihi, bd%bx%jlo : bd%bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'define_box_data: Out of memory'
            stop
        endif
    end subroutine define_box_data_bdry


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
    ! --------------------------------------------------------------------------
    pure function stagger (min) result (mout)
        integer, intent(in) :: min
        integer             :: mout
        mout = 1 - min
    end function stagger


    ! --------------------------------------------------------------------------
    ! Returns true if all members of box1 and box2 are equal.
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
    ! --------------------------------------------------------------------------
    pure subroutine setval_valid (bd, val)
        type(box_data), intent(inout) :: bd
        real(dp), intent(in)          :: val

        bd%data (bd%valid%ilo : bd%valid%ihi, bd%valid%jlo : bd%valid%jhi) = val
    end subroutine setval_valid


    ! --------------------------------------------------------------------------
    ! Sets all ghosts to val.
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
    ! Computes ip = (arr1, arr2)
    ! arr1 and arr2 must be exactly the same size and rank.
    ! --------------------------------------------------------------------------
    pure function inner_prod_array (arr1, arr2, nx, ny) result (ip)
        integer, intent(in)                      :: nx, ny
        real(dp), intent(in), dimension(1:nx*ny) :: arr1, arr2
        real(dp)                                 :: ip

        ip = dot_product(arr1, arr2)
    end function inner_prod_array


    ! --------------------------------------------------------------------------
    ! Computes ip = (bd1, bd2)
    ! bd1 and bd2 are box_data objects that must be exactly the same size.
    ! --------------------------------------------------------------------------
    function inner_prod_box_data (bd1, bd2) result (ip)
        type(box_data), intent(in) :: bd1, bd2
        real(dp)                   :: ip

        integer                    :: ilo, ihi, jlo, jhi, j

        ilo = max(bd1%valid%ilo, bd2%valid%ilo)
        ihi = min(bd1%valid%ihi, bd2%valid%ihi)
        jlo = max(bd1%valid%jlo, bd2%valid%jlo)
        jhi = min(bd1%valid%jhi, bd2%valid%jhi)

        ! Bounds checks
        if ((ilo .gt. ihi) .or. (jlo .gt. jhi)) then
            print*, 'ERROR: inner_prod_box_data: given box_data do not overlap.'
            stop
        endif

        ip = zero
        do j = jlo, jhi
            ip = ip + dot_product (bd1%data(ilo:ihi,j), bd2%data(ilo:ihi,j))
        enddo
    end function inner_prod_box_data


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

        ! if ((bd%offi .ne. BD_CELL) .or. (bd%offj .ne. BD_CELL)) then
        !     print*, 'pnorm: Only works with cell-centered data for now.'
        !     stop
        ! endif

        volScale = bd%bx%dx * bd%bx%dy

        pn = zero
        do j = bx%jlo, bx%jhi
            do i = bx%ilo, bx%ihi
                pn = pn + volScale * bd%data(i,j)**p
            enddo
        enddo
        pn = pn**(one/p)
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
    end function gpnorm


    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    subroutine fill_ghosts (phi, bcd, geo, homog, do_neum_opt)
        type(box_data), intent(inout) :: phi
        type(bdry_data), intent(in)   :: bcd
        type(geo_data), intent(in)    :: geo
        logical, intent(in)           :: homog
        logical, intent(in), optional :: do_neum_opt

        integer  :: xlo, xhi, ylo, yhi
        integer  :: ilo, ihi, jlo, jhi
        integer  :: i, j, e
        real(dp) :: dx, dy, bcval
        real(dp) :: cross, dpn, dpf
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

        xlo = bcd%type_xlo
        xhi = bcd%type_xhi
        ylo = bcd%type_ylo
        yhi = bcd%type_yhi

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

        if (phi%ngx .gt. 0) then
            ! Lower x ghosts
            select case (xlo)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            e = -1
                            j = jlo
                            i = ilo
                            bcval = zero

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                        else if (bcd%mode_xlo .eq. BCMODE_UNIFORM) then
                            e = -1
                            j = jlo
                            i = ilo
                            bcval = bcd%data_xlo(1)

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                        else
                            e = -1
                            j = jlo
                            i = ilo

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xlo(j)-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xlo(j)-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xlo(j)-cross)*dx/geo%Jgup_xx%data(i,j)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo-1, jlo:jhi) = -phi%data(ilo, jlo:jhi)
                    else if (bcd%mode_xlo .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_xlo(1)
                        phi%data(ilo-1, jlo:jhi) = two*bcval - phi%data(ilo, jlo:jhi)
                    else
                        phi%data(ilo-1,jlo:jhi) = two*bcd%data_xlo(jlo:jhi) - phi%data(ilo,jlo:jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo-1, jlo:jhi) = phi%data(ihi, jlo:jhi)

                case (BCTYPE_CF)
                    phi%data(ilo-1, jlo:jhi) = c1x*phi%data(ilo, jlo:jhi) + c2x*phi%data(ilo+1, jlo:jhi)

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select

            ! Upper x ghosts
            select case (xhi)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            e = 1
                            j = jlo
                            i = ilo
                            bcval = zero

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                        else if (bcd%mode_xhi .eq. BCMODE_UNIFORM) then
                            e = 1
                            j = jlo
                            i = ilo
                            bcval = bcd%data_xhi(1)

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcval-cross)*dx/geo%Jgup_xx%data(i,j)

                        else
                            e = 1
                            j = jlo
                            i = ilo

                            dpn = -(three*phi%data(i  ,j) - four*phi%data(i  ,j+1) + phi%data(i  ,j+2)) * half/dy
                            dpf = -(three*phi%data(i-e,j) - four*phi%data(i-e,j+1) + phi%data(i-e,j+2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xhi(j)-cross)*dx/geo%Jgup_xx%data(i,j)

                            do j = jlo+1, jhi-1
                                dpn = (phi%data(i  ,j+1) - phi%data(i  ,j-1)) * half/dy
                                dpf = (phi%data(i-e,j+1) - phi%data(i-e,j-1)) * half/dy
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_xy%data(i,j)
                                phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xhi(j)-cross)*dx/geo%Jgup_xx%data(i,j)
                            enddo

                            j = jhi
                            dpn = (three*phi%data(i  ,j) - four*phi%data(i  ,j-1) + phi%data(i  ,j-2)) * half/dy
                            dpf = (three*phi%data(i-e,j) - four*phi%data(i-e,j-1) + phi%data(i-e,j-2)) * half/dy
                            cross = half*(three*dpn - dpf) * geo%Jgup_xy%data(i,j)
                            phi%data(i+e,j) = phi%data(i,j) - (bcd%data_xhi(j)-cross)*dx/geo%Jgup_xx%data(i,j)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ihi+1, jlo:jhi) = -phi%data(ihi, jlo:jhi)
                    else if (bcd%mode_xhi .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_xhi(1)
                        phi%data(ihi+1, jlo:jhi) = two*bcval - phi%data(ihi, jlo:jhi)
                    else
                        phi%data(ihi+1, jlo:jhi) = two*bcd%data_xhi(jlo:jhi) - phi%data(ihi, jlo:jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ihi+1, jlo:jhi) = phi%data(ilo, jlo:jhi)

                case (BCTYPE_CF)
                    phi%data(ihi+1, jlo:jhi) = c1x*phi%data(ihi, jlo:jhi) + c2x*phi%data(ihi-1, jlo:jhi)

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select
        endif

        if (phi%ngy .gt. 0) then
            ! Lower y ghosts
            select case (ylo)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            e = -1
                            j = jlo
                            i = ilo
                            bcval = zero

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                        else if (bcd%mode_ylo .eq. BCMODE_UNIFORM) then
                            e = -1
                            j = jlo
                            i = ilo
                            bcval = bcd%data_ylo(1)

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                        else
                            e = -1
                            j = jlo
                            i = ilo

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcd%data_ylo(i)-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcd%data_ylo(i)-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcd%data_ylo(i)-cross)*dy/geo%Jgup_yy%data(i,j)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo:ihi, jlo-1) = -phi%data(ilo:ihi, jlo)
                    else if (bcd%mode_ylo .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_ylo(1)
                        phi%data(ilo:ihi, jlo-1) = two*bcval - phi%data(ilo:ihi, jlo)
                    else
                        phi%data(ilo:ihi, jlo-1) = two*bcd%data_ylo(ilo:ihi) - phi%data(ilo:ihi, jlo)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jhi)

                case (BCTYPE_CF)
                    phi%data(ilo:ihi, jlo-1) = c1x*phi%data(ilo:ihi, jlo) + c2x*phi%data(ilo:ihi, jlo+1)

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select

            ! Upper y ghosts
            select case (yhi)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            e = 1
                            j = jhi
                            i = ilo
                            bcval = zero

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                        else if (bcd%mode_yhi .eq. BCMODE_UNIFORM) then
                            e = 1
                            j = jhi
                            i = ilo
                            bcval = bcd%data_yhi(1)

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcval-cross)*dy/geo%Jgup_yy%data(i,j)

                        else
                            e = 1
                            j = jhi
                            i = ilo

                            dpn = -(three*phi%data(i,j  ) - four*phi%data(i+1,j  ) + phi%data(i+2,j  )) * half/dx
                            dpf = -(three*phi%data(i,j-e) - four*phi%data(i+1,j-e) + phi%data(i+2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcd%data_yhi(i)-cross)*dy/geo%Jgup_yy%data(i,j)

                            do i = ilo+1, ihi-1
                                dpn = (phi%data(i+1,j  ) - phi%data(i-1,j  )) * half/dx
                                dpf = (phi%data(i+1,j-e) - phi%data(i-1,j-e)) * half/dx
                                cross = (threehalves*dpn - half*dpf) * geo%Jgup_yx%data(i,j)
                                phi%data(i,j+e) = phi%data(i,j) - (bcd%data_yhi(i)-cross)*dy/geo%Jgup_yy%data(i,j)
                            enddo

                            i = ihi
                            dpn = (three*phi%data(i,j  ) - four*phi%data(i-1,j  ) + phi%data(i-2,j  )) * half/dx
                            dpf = (three*phi%data(i,j-e) - four*phi%data(i-1,j-e) + phi%data(i-2,j-e)) * half/dx
                            cross = half*(three*dpn - dpf) * geo%Jgup_yx%data(i,j)
                            phi%data(i,j+e) = phi%data(i,j) - (bcd%data_yhi(i)-cross)*dy/geo%Jgup_yy%data(i,j)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo:ihi, jhi+1) = -phi%data(ilo:ihi, jhi)
                    else if (bcd%mode_yhi .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_yhi(1)
                        phi%data(ilo:ihi, jhi+1) = two*bcval - phi%data(ilo:ihi, jhi)
                    else
                        phi%data(ilo:ihi, jhi+1) = two*bcd%data_yhi(ilo:ihi) - phi%data(ilo:ihi, jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jlo)

                case (BCTYPE_CF)
                    phi%data(ilo:ihi, jhi+1) = c1x*phi%data(ilo:ihi, jhi) + c2x*phi%data(ilo:ihi, jhi-1)

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select
        endif

        ! Put nans in the corner ghosts
        phi%data(ilo-1,jlo-1) = bogus_val
        phi%data(ilo-1,jhi+1) = bogus_val
        phi%data(ihi+1,jlo-1) = bogus_val
        phi%data(ihi+1,jhi+1) = bogus_val

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
                else
                    xflux%data(xflux%valid%ilo,:) = bc%data_xlo
                endif
            endif

            if (bc%type_xhi .eq. BCTYPE_NEUM) then
                if (bc%mode_xhi .eq. BCMODE_UNIFORM) then
                    xflux%data(xflux%valid%ihi,:) = bc%data_xhi(1)
                else
                    xflux%data(xflux%valid%ihi,:) = bc%data_xhi
                endif
            endif

            if (bc%type_ylo .eq. BCTYPE_NEUM) then
                if (bc%mode_ylo .eq. BCMODE_UNIFORM) then
                    yflux%data(:,yflux%valid%jlo) = bc%data_ylo(1)
                else
                    yflux%data(:,yflux%valid%jlo) = bc%data_ylo
                endif
            endif

            if (bc%type_yhi .eq. BCTYPE_NEUM) then
                if (bc%mode_yhi .eq. BCMODE_UNIFORM) then
                    yflux%data(:,yflux%valid%jhi) = bc%data_yhi(1)
                else
                    yflux%data(:,yflux%valid%jhi) = bc%data_yhi
                endif
            endif
        endif

    end subroutine fill_boundary_fluxes

end module ArrayUtils



! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module MGPoisson2D
    use ArrayUtils
    implicit none

    save

contains

    ! --------------------------------------------------------------------------
    ! Computes 1.0 / the Laplacian's diagonal matrix elements.
    ! invDiags must be prepared (allocated and box set) prior to call.
    ! --------------------------------------------------------------------------
    pure subroutine compute_inverse_diags (idiags, geo)
        type(box_data), intent(inout) :: idiags
        type(geo_data), intent(in)    :: geo

        real(dp)                      :: invdxsq, invdysq
        integer                       :: ilo, ihi, jlo, jhi

        invdxsq = -one / (idiags%valid%dx**2)
        invdysq = -one / (idiags%valid%dy**2)

        ilo = idiags%valid%ilo
        ihi = idiags%valid%ihi
        jlo = idiags%valid%jlo
        jhi = idiags%valid%jhi

        idiags%data(ilo:ihi,jlo:jhi) = &
            ( (geo%Jgup_xx%data(ilo+1:ihi+1,jlo:jhi) + geo%Jgup_xx%data(ilo:ihi,jlo:jhi)) * invdxsq    &
             +(geo%Jgup_yy%data(ilo:ihi,jlo+1:jhi+1) + geo%Jgup_yy%data(ilo:ihi,jlo:jhi)) * invdysq  ) &
            / geo%J%data(ilo:ihi,jlo:jhi)

        idiags%data = one / idiags%data

    end subroutine compute_inverse_diags


    ! ------------------------------------------------------------------------------
    ! Computes partial derivatives.
    ! phi is expected to be cell-centered with ghosts filled if needed.
    ! pd is expected to be face-centered with no ghosts.
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

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        ! Find the nodal direction
        nodedir = 1
        if (pd%offj .eq. 0) nodedir = 2

        ! Compute xflux...

        if (dir .eq. 1) then
            if (nodedir .eq. dir) then
                ! TODO: Need to fill ghosts.
                scale = one / phi%valid%dx
                pd%data(ilo:ihi+1,jlo:jhi) = scale * (phi%data(ilo:ihi+1,jlo:jhi) - phi%data(ilo-1:ihi,jlo:jhi))
            else
                scale = fourth / phi%valid%dx

                ! Away from y boundaries...
                ! Interior
                pd%data(ilo+1:ihi-1,jlo+1:jhi) = scale * (  phi%data(ilo+2:ihi,jlo+1:jhi) - phi%data(ilo:ihi-2,jlo+1:jhi) &
                                                          + phi%data(ilo+2:ihi,jlo:jhi-1) - phi%data(ilo:ihi-2,jlo:jhi-1)  )
                ! Upper x boundary
                pd%data(ihi,jlo+1:jhi) = scale * (  three*phi%data(ihi,jlo+1:jhi) - four*phi%data(ihi-1,jlo+1:jhi) + phi%data(ihi-2,jlo+1:jhi) &
                                                  + three*phi%data(ihi,jlo:jhi-1) - four*phi%data(ihi-1,jlo:jhi-1) + phi%data(ihi-2,jlo:jhi-1)  )

                ! Lower x boundary
                pd%data(ilo,jlo+1:jhi) = -scale * (  three*phi%data(ilo,jlo+1:jhi) - four*phi%data(ilo+1,jlo+1:jhi) + phi%data(ilo+2,jlo+1:jhi) &
                                                   + three*phi%data(ilo,jlo:jhi-1) - four*phi%data(ilo+1,jlo:jhi-1) + phi%data(ilo+2,jlo:jhi-1)  )

                ! At y boundaries...
                ! Lower y boundary
                pd%data(ilo:ihi,jlo) = two*pd%data(ilo:ihi,jlo+1) - pd%data(ilo:ihi,jlo+2)

                ! Upper y boundary
                pd%data(ilo:ihi,jhi+1) = two*pd%data(ilo:ihi,jhi) - pd%data(ilo:ihi,jhi-1)
            endif
        else
            if (nodedir .eq. dir) then
                ! TODO: Need to fill ghosts.
                scale = one / phi%valid%dy
                pd%data(ilo:ihi,jlo:jhi+1) = scale * (phi%data(ilo:ihi,jlo:jhi+1) - phi%data(ilo:ihi,jlo-1:jhi))
            else
                scale = fourth / phi%valid%dy

                ! Away from x boundaries...
                ! Interior
                pd%data(ilo+1:ihi,jlo+1:jhi-1) = scale * (  phi%data(ilo+1:ihi,jlo+2:jhi) - phi%data(ilo+1:ihi,jlo:jhi-2) &
                                                          + phi%data(ilo:ihi-1,jlo+2:jhi) - phi%data(ilo:ihi-1,jlo:jhi-2)  )
                ! Upper y boundary
                pd%data(ilo+1:ihi,jhi) = scale * (  three*phi%data(ilo+1:ihi,jhi) - four*phi%data(ilo+1:ihi,jhi-1) + phi%data(ilo+1:ihi,jhi-2) &
                                                  + three*phi%data(ilo:ihi-1,jhi) - four*phi%data(ilo:ihi-1,jhi-1) + phi%data(ilo:ihi-1,jhi-2)  )

                ! Lower y boundary
                pd%data(ilo+1:ihi,jlo) = -scale * (  three*phi%data(ilo+1:ihi,jlo) - four*phi%data(ilo+1:ihi,jlo+1) + phi%data(ilo+1:ihi,jlo+2) &
                                                   + three*phi%data(ilo:ihi-1,jlo) - four*phi%data(ilo:ihi-1,jlo+1) + phi%data(ilo:ihi-1,jlo+2)  )

                ! At x boundaries...
                ! Lower x boundary
                pd%data(ilo,jlo:jhi) = two*pd%data(ilo+1,jlo:jhi) - pd%data(ilo+2,jlo:jhi)

                ! Upper x boundary
                pd%data(ihi+1,jlo:jhi) = two*pd%data(ihi,jlo:jhi) - pd%data(ihi-1,jlo:jhi)
            endif
        endif

    end subroutine compute_pd


    ! ------------------------------------------------------------------------------
    ! Computes J*Grad[Phi].
    ! phi is expected to be cell-centered and have at least 1 ghost layer.
    ! *flux is expected to be face-centered with no ghosts.
    ! bc* = 0 for Neum, 1 for Diri, 2 for Periodic, 3 for CF.
    ! xwk is workspace that must be node-centered in x and cell-centered in y.
    ! ywk is workspace that must be node-centered in y and cell-centered in x.
    ! ------------------------------------------------------------------------------
    subroutine compute_grad (xflux, yflux, phi, geo, bc, homog, xwk, ywk)
        type(box_data), intent(inout) :: xflux, yflux
        type(box_data), intent(inout) :: phi
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog
        type(box_data), intent(inout) :: xwk, ywk

        integer                       :: ilo, ihi, jlo, jhi

        ilo = phi%valid%ilo
        ihi = phi%valid%ihi
        jlo = phi%valid%jlo
        jhi = phi%valid%jhi

        ! Fill ghosts
        call fill_ghosts (phi, bc, geo, homog, .false.)

        ! Compute xflux...
        call compute_pd (xflux, phi, 1)
        call compute_pd (xwk, phi, 2)
        xflux%data(ilo:ihi+1,jlo:jhi) = geo%Jgup_xx%data(ilo:ihi+1,jlo:jhi) * xflux%data(ilo:ihi+1,jlo:jhi) &
                                      + geo%Jgup_xy%data(ilo:ihi+1,jlo:jhi) *   xwk%data(ilo:ihi+1,jlo:jhi)

        ! Compute yflux...
        call compute_pd (ywk, phi, 1)
        call compute_pd (yflux, phi, 2)
        yflux%data(ilo:ihi,jlo:jhi+1) = geo%Jgup_yx%data(ilo:ihi,jlo:jhi+1) *   ywk%data(ilo:ihi,jlo:jhi+1) &
                                      + geo%Jgup_yy%data(ilo:ihi,jlo:jhi+1) * yflux%data(ilo:ihi,jlo:jhi+1)

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
    ! ------------------------------------------------------------------------------
    subroutine compute_laplacian (lap, phi, geo, bc, homog)
        type(box_data), intent(inout) :: lap
        type(box_data), intent(inout) :: phi
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        type(box_data)                :: xflux, yflux
        type(box_data)                :: xwk, ywk

        ! Allocate scratch space
        call define_box_data (xflux, phi%valid, 0, 0, BD_NODE, BD_CELL)
        call define_box_data (yflux, phi%valid, 0, 0, BD_CELL, BD_NODE)
        call define_box_data (xwk, xflux)
        call define_box_data (ywk, yflux)

        ! Compute Div[Grad[phi]]
        call compute_grad (xflux, yflux, phi, geo, bc, homog, xwk, ywk)
        call compute_div (lap, xflux, yflux)

        ! Free memory
        call undefine_box_data (xflux)
        call undefine_box_data (yflux)
        call undefine_box_data (xwk)
        call undefine_box_data (ywk)

    end subroutine compute_laplacian


    ! ------------------------------------------------------------------------------
    ! Compute res = rhs - L[phi].
    ! ------------------------------------------------------------------------------
    subroutine compute_residual (res, rhs, phi, geo, bc, homog)
        type(box_data), intent(inout) :: res, phi
        type(box_data), intent(in)    :: rhs
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        integer                       :: ilo, ihi
        integer                       :: jlo, jhi

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        call compute_laplacian (res, phi, geo, bc, homog)
        res%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) - res%data(ilo:ihi,jlo:jhi)

    end subroutine compute_residual


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine relax_jacobi (phi, rhs, geo, bc, homog, &
                             invdiags, omega, tol, maxiters, zerophi)
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

        type(box_data)                  :: r
        real(dp)                        :: rscale
        real(dp), dimension(0:maxiters) :: relres
        real(dp)                        :: newphi
        integer                         :: iter
        integer                         :: ilo, ihi, i
        integer                         :: jlo, jhi, j

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
        rscale = inner_prod (r, r)
        relres(0) = one
        print*, 'scale sq res = ', rscale
        print*, 'iter ', 0, ': sq res = ', relres(0)

        ! Iterate
        do iter = 1, maxIters
            ! Update phi
            do j = jlo, jhi
                do i = ilo, ihi
                    newphi = phi%data(i,j) + r%data(i,j)*invdiags%data(i,j)
                    phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                enddo
            enddo

            ! Compute new residual
            call compute_residual (r, rhs, phi, geo, bc, homog)
            relres(iter) = inner_prod(r, r) / rscale
            print*, 'iter ', iter, ': sq res = ', relres(iter)

            ! Did we converge?
            if (relres(iter) .le. tol) then
                print*, "Converged."
                exit
            endif
        enddo

    end subroutine relax_jacobi


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine relax_gs (phi, rhs, geo, bc, homog, &
                         invdiags, omega, tol, maxiters, redblack, zerophi)
        type(box_data), intent(inout)   :: phi
        type(box_data), intent(in)      :: rhs
        type(geo_data), intent(in)      :: geo
        type(bdry_data), intent(in)     :: bc
        logical, intent(in)             :: homog
        type(box_data), intent(in)      :: invdiags
        real(dp), intent(in)            :: omega
        real(dp), intent(in)            :: tol
        integer, intent(in)             :: maxiters
        logical, intent(in)             :: redblack
        logical, intent(in)             :: zerophi

        type(box_data)                  :: r
        real(dp)                        :: rscale
        real(dp), dimension(0:maxiters) :: relres
        integer                         :: iter, color
        real(dp)                        :: newphi
        integer                         :: ilo, ihi, i, imin
        integer                         :: jlo, jhi, j

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        ! Do we even need to be here?
        if (maxiters .eq. 0) then
            return
        endif

        ! Allocate workspace
        call define_box_data (r, rhs)

        ! Initialize phi to zero
        if (zerophi .ne. 0) then
            phi%data = zero
        endif

        ! Compute initial residual
        call compute_residual (r, rhs, phi, geo, bc, homog)
        rscale = inner_prod (r, r)
        relres(0) = one
        print*, 'scale sq res = ', rscale
        print*, 'iter ', 0, ': sq res = ', relres(0)

        ! Iterate
        do iter = 1, maxiters
            if (redblack) then
                ! Update phi via Red-Black Gauss-Seidel
                if (omega .eq. one) then
                    color = 0
                    do j = jlo, jhi
                        imin = ilo + mod(color+j, 2) ! Removed +1
                        phi%data(imin:ihi:2,j) = phi%data(imin:ihi:2,j) &
                                               + r%data(imin:ihi:2,j) * invdiags%data(imin:ihi:2,j)
                    enddo

                    color = 1
                    call compute_residual (r, rhs, phi, geo, bc, homog)
                    do j = jlo, jhi
                        imin = ilo + mod(color+j, 2) ! Removed +1
                        phi%data(imin:ihi:2,j) = phi%data(imin:ihi:2,j) &
                                               + r%data(imin:ihi:2,j) * invdiags%data(imin:ihi:2,j)
                    enddo
                else
                    color = 0
                    do j = jlo, jhi
                        imin = ilo + mod(color+j, 2) ! Removed +1
                        do i = imin, ihi, 2
                            newphi = phi%data(i,j) + r%data(i,j) * invdiags%data(i,j)
                            phi%data(i,j) = (one-omega)*phi%data(i,j) + omega*newphi
                        enddo
                    enddo

                    color = 1
                    call compute_residual (r, rhs, phi, geo, bc, homog)
                    do j = jlo, jhi
                        imin = ilo + mod(color+j, 2) ! Removed +1
                        do i = imin, ihi, 2
                            newphi = phi%data(i,j) + r%data(i,j) * invdiags%data(i,j)
                            phi%data(i,j) = (one-omega)*phi%data(i,j) + omega*newphi
                        enddo
                    enddo
                endif
            else
                ! Update phi via standard Gauss-Seidel
                if (omega .eq. one) then
                    phi%data(ilo:ihi,jlo:jhi) = phi%data(ilo:ihi,jlo:jhi) &
                                              + r%data(ilo:ihi,jlo:jhi) * invdiags%data(ilo:ihi,jlo:jhi)
                else
                    do j = jlo, jhi
                        do i = ilo, ihi
                            newphi = phi%data(i,j) + r%data(i,j) * invdiags%data(i,j)
                            phi%data(i,j) = (one-omega)*phi%data(i,j) + omega*newphi
                        enddo
                    enddo
                endif
            endif

            ! Compute new residual
            call compute_residual (r, rhs, phi, geo, bc, .true.)
            relres(iter) = inner_prod (r, r) / rscale
            print*, 'iter ', iter, ': sq res = ', relres(iter)

            ! Did we converge?
            if (relres(iter) .le. tol) then
                print*, "Converged."
                exit
            endif
        enddo

        ! Free memory
        call undefine_box_data (r)

    end subroutine relax_gs


    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    subroutine solve_bicgstab (phi, rhs, geo, bc, homog, tol, max_iters, max_restarts, zerophi)
        type(box_data), intent(inout)                 :: phi
        type(box_data), intent(in)                    :: rhs
        type(geo_data), intent(in)                    :: geo
        type(bdry_data), intent(in)                   :: bc
        logical, intent(in)                           :: homog
        real(dp), intent(in)                          :: tol
        integer, intent(in)                           :: max_iters, max_restarts
        logical, intent(in)                           :: zerophi

        integer                                       :: ilo, ihi
        integer                                       :: jlo, jhi
        type(box_data)                                :: r, r0, nu, p, t
        real(dp)                                      :: rscale
        real(dp), dimension(0:max_iters+max_restarts) :: rho, omega, relres
        real(dp)                                      :: alpha, beta, lastres
        integer                                       :: iter, i, num_restarts
        logical                                       :: is_restart

        real(dp), parameter                           :: hang = 1.0E-7_dp

        ! Do we even need to be here?
        if (max_iters .eq. 0) then
            return
        endif

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        ! Allocate workspace
        call define_box_data (r, rhs)
        call define_box_data (r0, rhs)
        call define_box_data (nu, phi)
        call define_box_data (p, rhs)
        call define_box_data (t, rhs)

        ! Initialize phi to zero
        if (zerophi .ne. 0) then
            phi%data = zero
        endif

        is_restart = .false.
        i = 0
        relres = zero

        ! Compute initial residual
        call compute_residual (r, rhs, phi, geo, bc, homog)
        r0%data = r%data
        rscale = inner_prod (r0, r0)
        relres(0) = one
        print*, 'scale sq res = ', rscale
        print*, 'iter ', 0, ': sq res = ', relres(0)

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

            rho(i) = inner_prod (r0, r)
            beta = (rho(i) / rho(i-1)) * (alpha / omega(i-1))
            p%data = beta*p%data
            p%data(ilo:ihi,jlo:jhi) = p%data(ilo:ihi,jlo:jhi)                  &
                                    - beta*omega(i-1)*nu%data(ilo:ihi,jlo:jhi) &
                                    + r%data(ilo:ihi,jlo:jhi)

            ! A preconditioner would go here
            call compute_laplacian (nu, p, geo, bc, homog)
            ! alpha = inner_prod (r0, nu)
            ! alpha = rho(i) / alpha
            ! r%data(ilo:ihi,jlo:jhi) = r%data(ilo:ihi,jlo:jhi) - alpha*nu%data(ilo:ihi,jlo:jhi)

        !     ! A preconditioner would go here
        !     call compute_laplacian (t, r, geo, bc, homog)
        !     omega(i) = inner_prod (t, r) / inner_prod (t, t)

        !     ! This would also change with a preconditioner
        !     phi%data(ilo:ihi,jlo:jhi) = phi%data(ilo:ihi,jlo:jhi)         &
        !                               + alpha * p%data(ilo:ihi,jlo:jhi)   &
        !                               + omega(i) * r%data(ilo:ihi,jlo:jhi)

        !     ! Compute new residual
        !     r%data = r%data - omega(i)*t%data

        !     ! If this is a restart, we expect the residual to rise.
        !     ! Don't let this stop the solver from proceeding.
        !     if (is_restart .eq. .false.) then
        !         lastres = relres(i-1)
        !     else
        !         lastres = 1.0E200_dp
        !     endif

        !     ! Check if we are at tol
        !     relres(i) = inner_prod (r, r) / rscale
        !     print*, 'iter ', iter, ': sq res = ', relres(i)

        !     ! Did we converge?
        !     if (relres(i) .le. tol) then
        !         print*, "Converged."
        !         exit
        !     endif

        !     ! Are we hanging?
        !     if (abs(relres(i) - lastres) .lt. hang) then
        !         if (num_restarts .lt. max_restarts) then
        !             ! The act of restarting will produce a new residual which we
        !             ! would like to include in our bookkeeping, so we increase i,
        !             ! recompute the residual, and reset all other bookkeeping vars.

        !             ! Increment
        !             num_restarts = num_restarts + 1
        !             i = i + 1

        !             ! Compute new residual
        !             call compute_residual (r, rhs, phi, geo, bc, homog)
        !             r0%data = r%data
        !             relres(i) = inner_prod (r0, r0) / rscale
        !             print*, "Hanging, restart number ", num_restarts, ', current sq res = ', relres(i)

        !             ! Reset bookkeeping variables
        !             alpha = one
        !             rho(i) = one
        !             omega(i) = one
        !             nu%data = zero
        !             p%data = zero

        !             ! Start new iteration
        !             is_restart = .true.
        !             cycle
        !         else
        !             print*, "Hanging, max restarts reached."
        !             exit
        !         endif
        !     endif

        !     ! Are we diverging?
        !     if (relres(i) .gt. lastres) then
        !         print*, 'Diverging.'

        !         ! Undo last correction
        !         ! TODO: It would be better to remember the best solution
        !         r%data = r%data + omega(i)*t%data
        !         phi%data(ilo:ihi,jlo:jhi) = phi%data(ilo:ihi,jlo:jhi)         &
        !                                   - alpha * p%data(ilo:ihi,jlo:jhi)   &
        !                                   - omega(i) * r%data(ilo:ihi,jlo:jhi)
        !         exit
        !     endif

            is_restart = .false.
        enddo

        ! Free memory
        call undefine_box_data (r)
        call undefine_box_data (r0)
        call undefine_box_data (nu)
        call undefine_box_data (p)
        call undefine_box_data (t)

    end subroutine solve_bicgstab

end module MGPoisson2D

