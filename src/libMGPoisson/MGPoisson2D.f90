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
