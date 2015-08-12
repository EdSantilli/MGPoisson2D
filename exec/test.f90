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

    type(box)       :: valid
    type(box_data)  :: state
    type(bdry_data) :: bc


    print*, 'Testing ArrayUtils::fill_ghosts in Cartesian...'

    ! Set up domain
    call define_box (valid, 1, 32, 1, 32, one, one)

    ! Set up field
    call define_box_data (state, valid, 1, 1)
    state%data = -three

    ! Set up BCs
    call define_bdry_data (bc, valid, &
                           BCTYPE_DIRI, &   ! xlo
                           BCTYPE_DIRI, &   ! xhi
                           BCTYPE_DIRI, &   ! ylo
                           BCTYPE_DIRI, &   ! yhi
                           BCMODE_UNIFORM, &    ! xlo
                           BCMODE_UNIFORM, &    ! xhi
                           BCMODE_UNIFORM, &    ! ylo
                           BCMODE_UNIFORM)      ! yhi
    bc%data_xlo(1) = five
    bc%data_ylo(1) = six
    bc%data_xlo(1) = seven
    bc%data_ylo(1) = eight

    ! Fill ghosts
    call fill_ghosts (state, bc, .true.)

    ! Check field
    ! TODO

    ! Free memory
    call undefine_bdry_data (bc)
    call undefine_box_data (state)

end program test
