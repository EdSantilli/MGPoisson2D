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
module Precision
    integer, parameter:: dp = kind(0.d0)
end module Precision


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module ArrayUtils
    use Precision
    implicit none

    ! --------------------------------------------------------------------------
    ! Array index extents and size info.
    ! --------------------------------------------------------------------------
    type Box
        integer  :: ilo, ihi        ! Lower and upper array indices in x dir
        integer  :: jlo, jhi        ! Lower and upper array indices in y dir
        integer  :: nx, ny          ! Total number of cells/nodes in each dir
        real(dp) :: dx, dy          ! The grid spacing
    end type Box


    ! --------------------------------------------------------------------------
    ! Used to identify the type of BC at each boundary.
    ! --------------------------------------------------------------------------
    type BCTypes
        integer :: xlo, xhi
        integer :: ylo, yhi
    end type BCTypes


    ! --------------------------------------------------------------------------
    ! These are the valid types of BCs.
    ! --------------------------------------------------------------------------
    integer, parameter :: UNDEFINED = -2
    integer, parameter :: NONE = -1
    integer, parameter :: NEUM = 0
    integer, parameter :: DIRI = 1
    integer, parameter :: PERIODIC = 2
    integer, parameter :: CF = 3


    ! --------------------------------------------------------------------------
    ! Used to hold data at the faces of each boundary.
    ! --------------------------------------------------------------------------
    type BoundaryData
        ! Valid (not ghost) box location and size info.
        type(Box) :: valid

        ! The boundary data that surrounds the domain.
        real(dp), dimension(:), allocatable :: xloData
        real(dp), dimension(:), allocatable :: xhiData
        real(dp), dimension(:), allocatable :: yloData
        real(dp), dimension(:), allocatable :: yhiData
    end type BoundaryData


    ! --------------------------------------------------------------------------
    ! These objects hold the data along with its metadata.
    ! --------------------------------------------------------------------------
    type BoxData
        type(Box) :: bx         ! The data box. Includes ghosts.
        type(Box) :: valid      ! The valid region only. Does not include ghosts.
        integer :: ngx, ngy     ! Nomber of ghosts per side in each dir

        ! The ghost and valid data array.
        real(dp), dimension(:,:), allocatable :: data

        ! Boundary condition information
        type(BCTypes) :: bcTypes
        type(BoundaryData) :: bdryData
    end type BoxData


    ! --------------------------------------------------------------------------
    ! Computes inner product of data arrays for norms.
    ! --------------------------------------------------------------------------
    interface InnerProd
        module procedure InnerProd_Array
        module procedure InnerProd_BoxData
    end interface

contains

    ! --------------------------------------------------------------------------
    ! This makes setting up a Box object a one line operation.
    ! --------------------------------------------------------------------------
    pure subroutine DefineBox (bx, ilo, ihi, jlo, jhi, dx, dy)
        type(Box), intent(out) :: bx
        integer, intent(in)    :: ilo, ihi
        integer, intent(in)    :: jlo, jhi
        real(8), intent(in)    :: dx, dy

        bx%ilo = ilo
        bx%ihi = ihi

        bx%jlo = jlo
        bx%jhi = jhi

        bx%dx = dx
        bx%dy = dy
    end subroutine DefineBox


    ! --------------------------------------------------------------------------
    ! This makes setting up a BoxData object a one line operation.
    ! --------------------------------------------------------------------------
    subroutine DefineBoxData (bd, valid, ngx, ngy)
        type(BoxData), intent(out) :: bd
        type(Box), intent(in)      :: valid
        integer, intent(in)        :: ngx, ngy

        type(Box) :: bx
        integer   :: ierr

        call DefineBox(bx, &
                       valid%ilo-ngx, valid%ihi+ngx, &
                       valid%jlo-ngy, valid%jhi+ngy, &
                       valid%dx, valid%dy)

        bd%valid = valid
        bd%bx = bx
        bd%ngx = ngx
        bd%ngy = ngy

        allocate (bd%data (bx%ilo : bx%ihi, bx%jlo : bx%jhi), stat=ierr)
        if (ierr .ne. 0) then
            print*, 'DefineDataBox: Out of memory'
            stop
        endif
    end subroutine DefineBoxData


    ! --------------------------------------------------------------------------
    ! Returns true if all members of box1 and box2 are equal.
    ! --------------------------------------------------------------------------
    pure function Compatible (bx1, bx2) result (isSame)
        implicit none
        type(Box), intent(in) :: bx1, bx2
        logical :: isSame

        isSame = .true.

        if (bx1%ilo .ne. bx2%ilo) isSame = .false.
        if (bx1%ihi .ne. bx2%ihi) isSame = .false.
        if (bx1%jlo .ne. bx2%jlo) isSame = .false.
        if (bx1%jhi .ne. bx2%jhi) isSame = .false.

        if (bx1%dx .ne. bx2%dx) isSame = .false.
        if (bx1%dy .ne. bx2%dy) isSame = .false.
    end function Compatible


    ! --------------------------------------------------------------------------
    ! Sets all valid data to val.
    ! --------------------------------------------------------------------------
    pure subroutine SetValValid (bd, val)
        implicit none
        type(BoxData), intent(inout) :: bd
        real(dp), intent(in)         :: val

        bd%data (bd%valid%ilo : bd%valid%ihi, bd%valid%jlo : bd%valid%jhi) = val
    end subroutine SetValValid


    ! --------------------------------------------------------------------------
    ! Sets all ghosts to val.
    ! --------------------------------------------------------------------------
    pure subroutine SetValGhosts (bd, val)
        implicit none
        type(BoxData), intent(inout) :: bd
        real(dp), intent(in)         :: val

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
    end subroutine SetValGhosts


    ! --------------------------------------------------------------------------
    ! Computes ip = (arr1, arr2)
    ! arr1 and arr2 must be exactly the same size and rank.
    ! --------------------------------------------------------------------------
    pure function InnerProd_Array (arr1, arr2, nx, ny) result (ip)
        implicit none
        integer, intent(in)                      :: nx, ny
        real(dp), intent(in), dimension(1:nx*ny) :: arr1, arr2
        real(dp)                                 :: ip

        ip = dot_product(arr1, arr2)
    end function InnerProd_Array


    ! --------------------------------------------------------------------------
    ! Computes ip = (bd1, bd2)
    ! bd1 and bd2 are BoxData objects that must be exactly the same size.
    ! --------------------------------------------------------------------------
    function InnerProd_BoxData (bd1, bd2) result (ip)
        implicit none
        type(BoxData), intent(in) :: bd1, bd2
        real(dp)                  :: ip

        ! Bounds checks
        if (.not.Compatible(bd1%valid, bd2%valid)) then
            print*, 'ERROR: InnerProd_BoxData: given BoxData are incompatible.'
            stop
        endif

        ! Compute over entire region (not just valid region).
        ip = InnerProd_Array(bd1%data, bd2%data, bd1%bx%nx, bd1%bx%ny)

    end function InnerProd_BoxData
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

