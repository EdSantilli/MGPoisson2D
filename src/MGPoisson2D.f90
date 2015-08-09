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

module BoundaryUtils
    implicit none

    ! BC types
    integer, parameter :: UNDEFINED = -2
    integer, parameter :: NONE = -1
    integer, parameter :: NEUM = 0
    integer, parameter :: DIRI = 1
    integer, parameter :: PERIODIC = 2
    integer, parameter :: CF = 3

    type BCTypes
        integer :: xlo, xhi
        integer :: ylo, yhi
    end type BCTypes

    type BoundaryData
        ! Valid (not ghost) box location and size info.
        integer :: ilo, ihi         ! Lower and upper array indices in x dir
        integer :: jlo, jhi         ! Lower and upper array indices in y dir
        integer :: nx, ny           ! Total number of cells/nodes in each dir

        ! The boundary data that surrounds the valid box.
        real(8), dimension(:), allocatable :: xloData
        real(8), dimension(:), allocatable :: xhiData
        real(8), dimension(:), allocatable :: yloData
        real(8), dimension(:), allocatable :: yhiData
    end type BoundaryData

contains

end module BoundaryUtils


module ArrayUtils
    use BoundaryUtils
    implicit none

    type BoxData
        ! Box location and size info.
        integer :: ilo, ihi         ! Lower and upper array indices in x dir
        integer :: jlo, jhi         ! Lower and upper array indices in y dir
        integer :: nx, ny           ! Total number of cells/nodes in each dir
        integer :: nvx, nvy         ! Number of valid (not ghost) cells/nodes in each dir
        integer :: ngx, ngy         ! Nomber of ghosts per side in each dir

        ! The ghost and valid data array.
        real(8), dimension(:,:), allocatable :: data

        ! Boundary condition information
        type(BCTypes) :: bcTypes
        type(BoundaryData) :: bdryData
    end type BoxData

contains

end module ArrayUtils


module MGPoisson2D
    use ArrayUtils
    implicit none
    private

    save

contains

end module MGPoisson2D

