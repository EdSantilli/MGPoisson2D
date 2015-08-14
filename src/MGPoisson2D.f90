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

    real(dp), parameter :: half      = 0.5_dp

    real(dp), parameter :: pi        = 4*atan(one)
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
    ! These objects hold CC state data along with its metadata.
    ! --------------------------------------------------------------------------
    type box_data
        type(box) :: bx         ! The data box. Includes ghosts.
        type(box) :: valid      ! The valid region only. Does not include ghosts.
        integer   :: ngx, ngy   ! Nomber of ghosts per side in each dir

        ! The ghost and valid data array.
        real(dp), dimension(:,:), allocatable :: data
    end type box_data


    ! --------------------------------------------------------------------------
    ! These objects hold FC state data along with its metadata.
    ! --------------------------------------------------------------------------
    type face_data
        type(box) :: bx         ! The FC data box. Includes ghosts.
        type(box) :: valid      ! The FC valid region only. Does not include ghosts.
        integer   :: ngx, ngy   ! Nomber of ghosts per side in each dir
        integer   :: facedir    ! The FC direction

        ! The ghost and valid data array.
        real(dp), dimension(:,:), allocatable :: data
    end type face_data


    ! --------------------------------------------------------------------------
    ! Contains:
    !   J = the CC Jacobian determinant
    !   Jgup = the FC inverse metric tensor scaled by J
    ! --------------------------------------------------------------------------
    type geo_data
        type(box_data)  :: J
        type(face_data) :: Jgup_x
        type(face_data) :: Jgup_y
        real(dp)        :: dx, dy
    end type geo_data


    ! --------------------------------------------------------------------------
    ! Computes inner product of data arrays for norms.
    ! --------------------------------------------------------------------------
    interface inner_prod
        module procedure inner_prod_array
        module procedure inner_prod_box_data
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
    ! --------------------------------------------------------------------------
    subroutine define_box_data (bd, valid, ngx, ngy)
        type(box_data), intent(out) :: bd
        type(box), intent(in)       :: valid
        integer, intent(in)         :: ngx, ngy

        type(box) :: bx
        integer   :: ierr

        call define_box(bx, &
                        valid%ilo-ngx, valid%ihi+ngx, &
                        valid%jlo-ngy, valid%jhi+ngy, &
                        valid%dx, valid%dy)

        bd%valid = valid
        bd%bx = bx
        bd%ngx = ngx
        bd%ngy = ngy

        allocate (bd%data (bx%ilo : bx%ihi, bx%jlo : bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'define_box_data: Out of memory'
            stop
        endif
    end subroutine define_box_data


    ! --------------------------------------------------------------------------
    ! Frees memory used by a box_data object.
    ! --------------------------------------------------------------------------
    pure subroutine undefine_box_data (bd)
        type(box_data), intent(inout) :: bd

        bd%valid = empty_box
        bd%bx = empty_box
        bd%ngx = 0
        bd%ngy = 0

        if (allocated(bd%data)) deallocate(bd%data)
    end subroutine undefine_box_data


    ! --------------------------------------------------------------------------
    ! This makes setting up a face_data object a one line operation.
    ! --------------------------------------------------------------------------
    subroutine define_face_data (fd, valid, facedir, ngx, ngy)
        type(face_data), intent(out) :: fd
        type(box), intent(in)        :: valid
        integer, intent(in)          :: facedir
        integer, intent(in)          :: ngx, ngy

        type(box) :: bx, fcvalid
        integer   :: ierr

        if (facedir .eq. 1) then
            fcvalid = valid
            fcvalid%ihi = fcvalid%ihi + 1
            call define_box(bx, &
                            valid%ilo-ngx, valid%ihi+ngx+1, &
                            valid%jlo-ngy, valid%jhi+ngy, &
                            valid%dx, valid%dy)

        else if (facedir .eq. 2) then
            fcvalid = valid
            fcvalid%jhi = fcvalid%jhi + 1
            call define_box(bx, &
                            valid%ilo-ngx, valid%ihi+ngx, &
                            valid%jlo-ngy, valid%jhi+ngy+1, &
                            valid%dx, valid%dy)

        else
            print*, 'define_face_data: Bad facedir'
            stop
        endif

        fd%valid = fcvalid
        fd%bx = bx
        fd%ngx = ngx
        fd%ngy = ngy
        fd%facedir = facedir

        allocate (fd%data (bx%ilo : bx%ihi, bx%jlo : bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'define_box_data: Out of memory'
            stop
        endif
    end subroutine define_face_data


    ! --------------------------------------------------------------------------
    ! Frees memory used by a face_data object.
    ! --------------------------------------------------------------------------
    pure subroutine undefine_face_data (fd)
        type(face_data), intent(inout) :: fd

        fd%valid = empty_box
        fd%bx = empty_box
        fd%ngx = 0
        fd%ngy = 0
        fd%facedir = 0

        if (allocated(fd%data)) deallocate(fd%data)
    end subroutine undefine_face_data


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

        ! Bounds checks
        if (.not.compatible_boxes(bd1%valid, bd2%valid)) then
            print*, 'ERROR: inner_prod_box_data: given box_data are incompatible.'
            stop
        endif

        ! Compute over entire region (not just valid region).
        ip = inner_prod_array(bd1%data, bd2%data, bd1%bx%nx, bd1%bx%ny)
    end function inner_prod_box_data


    ! --------------------------------------------------------------------------
    ! Computes pnorm = |bd|_p = [Sum_{i,j} bd^p]^(1/p) / (nx*ny)
    ! WARNING: This function does not check that bd contains bx.
    ! --------------------------------------------------------------------------
    pure function pnorm (bd, bx, p) result (pn)
        type(box_data), intent(in) :: bd
        type(box), intent(in)      :: bx
        integer, intent(in)        :: p
        real(dp)                   :: pn

        integer :: i, j
        real(dp) :: volScale

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
    ! --------------------------------------------------------------------------
    subroutine fill_ghosts (phi, bcd, homog, do_neum_opt)
        type(box_data), intent(inout) :: phi
        type(bdry_data), intent(in)   :: bcd
        logical, intent(in)           :: homog
        logical, intent(in), optional :: do_neum_opt

        integer  :: xlo, xhi, ylo, yhi
        integer  :: ilo, ihi, jlo, jhi
        real(dp) :: dx, dy, bcval
        logical  :: do_neum

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
                            phi%data(ilo-1, jlo:jhi) = phi%data(ilo, jlo:jhi)
                        else if (bcd%mode_xlo .eq. BCMODE_UNIFORM) then
                            bcval = dx * bcd%data_xlo(1)
                            phi%data(ilo-1, jlo:jhi) = phi%data(ilo, jlo:jhi) - bcval
                        else
                            phi%data(ilo-1, jlo:jhi) = phi%data(ilo, jlo:jhi) - dx * bcd%data_xlo(jlo:jhi)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo-1, jlo:jhi) = -phi%data(ilo, jlo:jhi)
                    else if (bcd%mode_xlo .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_xlo(1)
                        phi%data(ilo-1, jlo:jhi) = bcval - phi%data(ilo, jlo:jhi)
                    else
                        phi%data(ilo-1, jlo:jhi) = two * bcd%data_xlo(jlo:jhi) - phi%data(ilo, jlo:jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo-1, jlo:jhi) = phi%data(ihi, jlo:jhi)

                case (BCTYPE_CF)
                    ! TODO
                    print*, 'fill_ghosts: cannot handle BCTYPE_CF yet.'
                    stop

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select

            ! Upper x ghosts
            select case (xhi)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            phi%data(ihi+1, jlo:jhi) = phi%data(ihi, jlo:jhi)
                        else if (bcd%mode_xhi .eq. BCMODE_UNIFORM) then
                            bcval = dx * bcd%data_xhi(1)
                            phi%data(ihi+1, jlo:jhi) = phi%data(ihi, jlo:jhi) + bcval
                        else
                            phi%data(ihi+1, jlo:jhi) = phi%data(ihi, jlo:jhi) + dx * bcd%data_xhi(jlo:jhi)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ihi+1, jlo:jhi) = -phi%data(ihi, jlo:jhi)
                    else if (bcd%mode_xhi .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_xhi(1)
                        phi%data(ihi+1, jlo:jhi) = bcval - phi%data(ihi, jlo:jhi)
                    else
                        phi%data(ihi+1, jlo:jhi) = two * bcd%data_xhi(jlo:jhi) - phi%data(ihi, jlo:jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ihi+1, jlo:jhi) = phi%data(ilo, jlo:jhi)

                case (BCTYPE_CF)
                    ! TODO
                    print*, 'fill_ghosts: cannot handle BCTYPE_CF yet.'
                    stop

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
                            phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jlo)
                        else if (bcd%mode_ylo .eq. BCMODE_UNIFORM) then
                            bcval = dx * bcd%data_ylo(1)
                            phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jlo) - bcval
                        else
                            phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jlo) - dx * bcd%data_ylo(ilo:ihi)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo:ihi, jlo-1) = -phi%data(ilo:ihi, jlo)
                    else if (bcd%mode_ylo .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_ylo(1)
                        phi%data(ilo:ihi, jlo-1) = bcval - phi%data(ilo:ihi, jlo)
                    else
                        phi%data(ilo:ihi, jlo-1) = two * bcd%data_ylo(ilo:ihi) - phi%data(ilo:ihi, jlo)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo:ihi, jlo-1) = phi%data(ilo:ihi, jhi)

                case (BCTYPE_CF)
                    ! TODO
                    print*, 'fill_ghosts: cannot handle BCTYPE_CF yet.'
                    stop

                case default
                    print*, 'fill_ghosts: invalid BCTYPE_'
                    stop
            end select

            ! Upper y ghosts
            select case (yhi)
                case (BCTYPE_NEUM)
                    if (do_neum) then
                        if (homog) then
                            phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jhi)
                        else if (bcd%mode_yhi .eq. BCMODE_UNIFORM) then
                            bcval = dx * bcd%data_yhi(1)
                            phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jhi) + bcval
                        else
                            phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jhi) + dx * bcd%data_yhi(ilo:ihi)
                        endif
                    endif

                case (BCTYPE_DIRI)
                    if (homog) then
                        phi%data(ilo:ihi, jhi+1) = -phi%data(ilo:ihi, jhi)
                    else if (bcd%mode_yhi .eq. BCMODE_UNIFORM) then
                        bcval = two * bcd%data_yhi(1)
                        phi%data(ilo:ihi, jhi+1) = bcval - phi%data(ilo:ihi, jhi)
                    else
                        phi%data(ilo:ihi, jhi+1) = two * bcd%data_yhi(ilo:ihi) - phi%data(ilo:ihi, jhi)
                    endif

                case (BCTYPE_PERIODIC)
                    phi%data(ilo:ihi, jhi+1) = phi%data(ilo:ihi, jlo)

                case (BCTYPE_CF)
                    ! TODO
                    print*, 'fill_ghosts: cannot handle BCTYPE_CF yet.'
                    stop

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

end module ArrayUtils



! ! ------------------------------------------------------------------------------
! ! ------------------------------------------------------------------------------
! module MGPoisson2D
!     use ArrayUtils
!     implicit none
!     private

!     save

! contains

!     ! --------------------------------------------------------------------------
!     ! Computes 1.0 / the Laplacian's diagonal matrix elements.
!     ! invDiags must be prepared (allocated and box set) prior to call.
!     ! --------------------------------------------------------------------------
!     subroutine ComputeInvDiags (invDiags)
!         implicit none
!         type(BoxData), intent(inout) :: invDiags
!         real(dp)                     :: dxScale, dyScale

!         dxScale = 2.0d0 / (invDiags%box%dx**2)
!         dyScale = 2.0d0 / (invDiags%box%dy**2)
!         invDiags%data = -1.0d0 / (dxScale + dyScale)
!     end subroutine ComputeInvDiags


! end module MGPoisson2D

