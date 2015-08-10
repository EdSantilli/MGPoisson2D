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
    integer, parameter:: dp = kind(0.d0)
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
    ! These are the valid types of BCs.
    ! --------------------------------------------------------------------------
    integer, parameter :: BC_UNDEFINED = -2
    integer, parameter :: BC_NONE = -1
    integer, parameter :: BC_NEUM = 0
    integer, parameter :: BC_DIRI = 1
    integer, parameter :: BC_PERIODIC = 2
    integer, parameter :: BC_CF = 3


    ! --------------------------------------------------------------------------
    ! Used to hold BC types and values at the domain boundary faces.
    ! --------------------------------------------------------------------------
    type bdry_data
        ! Valid (not ghost) box location and size info.
        type(box) :: valid

        ! The types of BCs that the data represent.
        integer :: type_xlo, type_xhi
        integer :: type_ylo, type_yhi

        ! The boundary data that surrounds the domain.
        real(dp), dimension(:), allocatable :: data_xlo, data_xhi
        real(dp), dimension(:), allocatable :: data_ylo, data_yhi
    end type bdry_data


    ! --------------------------------------------------------------------------
    ! These objects hold the data along with its metadata.
    ! --------------------------------------------------------------------------
    type box_data
        type(box) :: bx         ! The data box. Includes ghosts.
        type(box) :: valid      ! The valid region only. Does not include ghosts.
        integer :: ngx, ngy     ! Nomber of ghosts per side in each dir

        ! The ghost and valid data array.
        real(dp), dimension(:,:), allocatable :: data

        ! BC info
        type(bdry_data) :: bc
    end type box_data


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
        real(8), intent(in)    :: dx, dy

        bx%ilo = ilo
        bx%ihi = ihi

        bx%jlo = jlo
        bx%jhi = jhi

        bx%dx = dx
        bx%dy = dy
    end subroutine define_box


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
    ! Returns true if all members of box1 and box2 are equal.
    ! --------------------------------------------------------------------------
    pure function compatible_boxes (bx1, bx2) result (is_same)
        implicit none
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
        implicit none
        type(box_data), intent(inout) :: bd
        real(dp), intent(in)          :: val

        bd%data (bd%valid%ilo : bd%valid%ihi, bd%valid%jlo : bd%valid%jhi) = val
    end subroutine setval_valid


    ! --------------------------------------------------------------------------
    ! Sets all ghosts to val.
    ! --------------------------------------------------------------------------
    pure subroutine setval_ghosts (bd, val)
        implicit none
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
        implicit none
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
        implicit none
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


    ! ! --------------------------------------------------------------------------
    ! ! --------------------------------------------------------------------------
    ! subroutine fillGhosts (bd)

    ! end subroutine fillGhosts

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

