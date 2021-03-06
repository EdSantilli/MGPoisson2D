! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module DeferredPoisson2D
    use ArrayUtils
    implicit none

    save

contains

    ! --------------------------------------------------------------------------
    ! Computes 1.0 / the Laplacian's diagonal matrix elements.
    ! invdiags must be prepared (allocated and box set) prior to call.
    ! NOTE: This assumes the Laplacian is not scaling by 1/J.
    ! --------------------------------------------------------------------------
    pure subroutine compute_invdiags (invdiags, geo)
        type(box_data), intent(inout) :: invdiags
        type(geo_data), intent(in)    :: geo

        integer                       :: ilo, ihi
        integer                       :: jlo, jhi
        real(dp)                      :: idx2, idy2

        ilo = invdiags%valid%ilo
        ihi = invdiags%valid%ihi
        jlo = invdiags%valid%jlo
        jhi = invdiags%valid%jhi

        idx2 = one / geo%dx**2
        idy2 = one / geo%dy**2

        ! Compute -diagonals.
        invdiags%data(ilo:ihi,jlo:jhi) = (geo%Jgup_xx%data(ilo+1:ihi+1,jlo:jhi) + geo%Jgup_xx%data(ilo:ihi,jlo:jhi)) * idx2 &
                                       + (geo%Jgup_yy%data(ilo:ihi,jlo+1:jhi+1) + geo%Jgup_yy%data(ilo:ihi,jlo:jhi)) * idy2

        ! Invert the diag values.
        invdiags%data = -one / invdiags%data

    end subroutine compute_invdiags


    ! ------------------------------------------------------------------------------
    ! This function will not scale the result by 1/J So, it is not a true Laplacian.
    ! ------------------------------------------------------------------------------
    subroutine compute_laplacian (lap, phi, geo, bc, homog)
        type(box_data), intent(inout) :: lap
        type(box_data), intent(inout) :: phi
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        logical, parameter            :: neum_ghosts = .true.

        integer                       :: ilo, ihi, i
        integer                       :: jlo, jhi, j
        real(dp)                      :: p,pe,pn,pw,ps
        real(dp)                      :: pne, pnw, pse, psw
        real(dp)                      :: gxye, gxyw, gyxn, gyxs
        real(dp)                      :: xfe, xfw, yfn, yfs
        real(dp)                      :: idx, idy, xyscale

        integer, parameter            :: extrap_order = 0

        ilo = lap%valid%ilo
        ihi = lap%valid%ihi
        jlo = lap%valid%jlo
        jhi = lap%valid%jhi

        idx = one / geo%dx
        idy = one / geo%dy
        xyscale = fourth*idx*idy

        ! Fill ghost cells, including corners
        call extrapolate_ghosts (phi, extrap_order)

        ! Compute the cross derivatives
        do j = jlo, jhi
            do i = ilo, ihi
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                pne = phi%data(i+1,j+1)
                pnw = phi%data(i-1,j+1)
                psw = phi%data(i-1,j-1)
                pse = phi%data(i+1,j-1)

                gxye = geo%Jgup_xy%data(i+1,j  )*(pne-pse+pn-ps)
                gxyw = geo%Jgup_xy%data(i  ,j  )*(pn-ps+pnw-psw)
                gyxn = geo%Jgup_yx%data(i  ,j+1)*(pne-pnw+pe-pw)
                gyxs = geo%Jgup_yx%data(i  ,j  )*(pe-pw+pse-psw)

                lap%data(i,j) = (gxye-gxyw+gyxn-gyxs)*xyscale
            enddo
        enddo

        ! Fill ghost cells
        call fill_ghosts (phi, bc, geo, homog, neum_ghosts)

        ! Add the dominant terms
        do j = jlo, jhi
            do i = ilo, ihi
                p  = phi%data(i  ,j  )
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx
                xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx
                yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy
                yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy

                lap%data(i,j) = lap%data(i,j) + (xfe-xfw)*idx + (yfn-yfs)*idy
            enddo
        enddo

    end subroutine compute_laplacian


    ! ------------------------------------------------------------------------------
    ! Compute res = rhs - L[phi].
    ! If opt_jscale is false (default), this function will not scale the result
    ! by 1/J. Note that a true Laplacian should do this.
    ! ------------------------------------------------------------------------------
    subroutine compute_residual (res, rhs, phi, geo, bc, homog)
        type(box_data), intent(inout) :: res, phi
        type(box_data), intent(in)    :: rhs
        type(geo_data), intent(in)    :: geo
        type(bdry_data), intent(in)   :: bc
        logical, intent(in)           :: homog

        integer                       :: ilo, ihi, jlo, jhi

        call compute_laplacian (res, phi, geo, bc, homog)

        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi
        res%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) - res%data(ilo:ihi,jlo:jhi)

    end subroutine compute_residual


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

        type(box_data)                  :: r, modrhs
        real(dp)                        :: rscale
        real(dp), dimension(0:maxiters) :: relres
        integer                         :: iter

        integer                         :: ilo, ihi, i
        integer                         :: jlo, jhi, j
        real(dp)                        :: p,pe,pn,pw,ps
        real(dp)                        :: xfe,xfw,yfn,yfs
        real(dp)                        :: pne, pnw, pse, psw
        real(dp)                        :: gxye, gxyw, gyxn, gyxs
        real(dp)                        :: idx, idy, xyscale
        real(dp)                        :: lphi, newphi, sum

        integer                         :: color, imin

        logical, parameter              :: neum_ghosts = .true.
        integer, parameter              :: extrap_order = 0
        logical, parameter              :: red_black = .false.

        ! Do we even need to be here?
        if (maxiters .eq. 0) then
            return
        endif

        ! Define some useful quantities
        ilo = rhs%valid%ilo
        ihi = rhs%valid%ihi
        jlo = rhs%valid%jlo
        jhi = rhs%valid%jhi

        idx = one / geo%dx
        idy = one / geo%dy
        xyscale = half*idx*idy

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

        ! Add cross derivative terms to r
        call define_box_data (modrhs, rhs)
        modrhs%data = rhs%data

        call extrapolate_ghosts (phi, extrap_order)
        do j = jlo, jhi
            do i = ilo, ihi
                pe = phi%data(i+1,j  )
                pw = phi%data(i-1,j  )
                pn = phi%data(i  ,j+1)
                ps = phi%data(i  ,j-1)

                pne = phi%data(i+1,j+1)
                pnw = phi%data(i-1,j+1)
                psw = phi%data(i-1,j-1)
                pse = phi%data(i+1,j-1)

                gxye = geo%Jgup_xy%data(i+1,j  )*(pne-pse)
                gxyw = geo%Jgup_xy%data(i  ,j  )*(pnw-psw)
                gyxn = geo%Jgup_yx%data(i  ,j+1)*(pne-pnw)
                gyxs = geo%Jgup_yx%data(i  ,j  )*(pse-psw)

                modrhs%data(i,j) = modrhs%data(i,j) - (gxye-gxyw+gyxn-gyxs)*xyscale
            enddo
        enddo

        ! Iterate
        do iter = 1, maxiters
            if (verbosity .ge. 5) then
                sum = integrate2d (r, r%valid, geo, .false.)
                print*, ' sum rhs = ', sum
            endif

            ! Update phi via standard Gauss-Seidel
            call fill_ghosts (phi, bc, geo, homog, neum_ghosts)

            if (red_black) then
                ! Red-black iteration.
                do color = 0, 1
                    do j = jlo, jhi
                        imin = ilo + mod(color+j, 2) ! Removed +1
                        do i = imin, ihi, 2
                            p  = phi%data(i  ,j  )
                            pe = phi%data(i+1,j  )
                            pw = phi%data(i-1,j  )
                            pn = phi%data(i  ,j+1)
                            ps = phi%data(i  ,j-1)

                            xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx
                            xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx
                            yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy
                            yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy

                            lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                            newphi = phi%data(i,j) + (modrhs%data(i,j) - lphi) * invdiags%data(i,j)
                            phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                        enddo
                    enddo
                enddo
            else
                ! Standard iteration.
                do j = jlo, jhi
                    do i = ilo, ihi
                        p  = phi%data(i  ,j  )
                        pe = phi%data(i+1,j  )
                        pw = phi%data(i-1,j  )
                        pn = phi%data(i  ,j+1)
                        ps = phi%data(i  ,j-1)

                        xfe = geo%Jgup_xx%data(i+1,j)*(pe-p)*idx
                        xfw = geo%Jgup_xx%data(i  ,j)*(p-pw)*idx
                        yfn = geo%Jgup_yy%data(i,j+1)*(pn-p)*idy
                        yfs = geo%Jgup_yy%data(i,j  )*(p-ps)*idy

                        lphi = (xfe-xfw)*idx + (yfn-yfs)*idy

                        newphi = phi%data(i,j) + (modrhs%data(i,j) - lphi) * invdiags%data(i,j)
                        phi%data(i,j) = omega*newphi + (one-omega)*phi%data(i,j)
                    enddo
                enddo
            endif

            ! Diagnostics
            if (tol .gt. zero) then
                ! Compute new residual
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
        call undefine_box_data (modrhs)
        if (tol .gt. zero) then
            call undefine_box_data (r)
        endif
    end subroutine relax_gs


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
        integer, parameter                            :: norm_type = 0

        real(dp), parameter                           :: relax_omega = four/five
        integer, parameter                            :: relax_iters = 2

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
                           relax_omega, &
                           -one,    &  ! tol
                           relax_iters, &
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
                               relax_omega, &
                               -one,    &  ! tol
                               relax_iters, &
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

end module DeferredPoisson2D



! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
module DeferredMGPoisson2D
    use ArrayUtils
    use MGPoisson2D, only: estimate_maxdepth, &
                           create_mgschedule, &
                           restrict, &
                           prolong
    implicit none

    save

contains

    ! ------------------------------------------------------------------------------
    ! rhs will be scaled with J, but restored before the function exits.
    ! ------------------------------------------------------------------------------
    subroutine vcycle (phi, rhs, geo, bc, homog, amrrefx, amrrefy, &
                       tol, maxiters, maxdepth_user, numcycles, &
                       smooth_down, smooth_up, smooth_bottom, &
                       zerophi, verbosity)

        use Poisson2D, only: compute_residual, compute_invdiags, relax_gs, compute_pd
        use DeferredPoisson2D, only: compute_deferred_invdiags => compute_invdiags

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

        character*2                                          :: indent = '  '
        real(dp)                                             :: norm


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
        if (verbosity .ge. 4) then
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
                                   BCMODE_LAH, BCMODE_LAH, BCMODE_LAH, BCMODE_LAH)
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
                ! call compute_invdiags (invdiags(d), geo, bc, .true.)
                call compute_deferred_invdiags (invdiags(d), mggeo(d))
            else
                ! Coarsen the data from MG level d-1
                mggeo(d)%dx = valid%dx
                mggeo(d)%dy = valid%dy
                call restrict (mggeo(d-1)%J, mggeo(d)%J)
                call restrict (mggeo(d-1)%Jgup_xx, mggeo(d)%Jgup_xx)
                call restrict (mggeo(d-1)%Jgup_xy, mggeo(d)%Jgup_xy)
                call restrict (mggeo(d-1)%Jgup_yx, mggeo(d)%Jgup_yx)
                call restrict (mggeo(d-1)%Jgup_yy, mggeo(d)%Jgup_yy)
                call compute_deferred_invdiags (invdiags(d), mggeo(d))
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
        rhs%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) &
                                  * geo%J%data(ilo:ihi,jlo:jhi)


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
            print*, 'Deferred correction V-Cycle iter ', 0, ': rel |res| = ', relres(0)
        endif

        ! The main iteration loop.
        do iter = 1, maxiters
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
            call compute_residual (r(0), rhs, phi, geo, bc, homog)
            relres(iter) = pnorm (r(0), r(0)%valid, 2) / rscale
            if (verbosity .ge. 1) then
                print*, 'Deferred correction V-Cycle iter ', iter, ': rel |res| = ', relres(iter)
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

        ! Restore rhs scaling.
        rhs%data(ilo:ihi,jlo:jhi) = rhs%data(ilo:ihi,jlo:jhi) &
                                  / geo%J%data(ilo:ihi,jlo:jhi)

        ! Free memory
        do d = 0, maxDepth
            call undefine_bdry_data (mgbc(d))
            call undefine_box_data (invdiags(d))
            call undefine_box_data (mggeo(d)%J)
            call undefine_box_data (mggeo(d)%Jgup_xx)
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
        use DeferredPoisson2D

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
        integer                                              :: curcycle
        real(dp)                                             :: norm, sum
        type(box_data)                                       :: tmp

        integer, parameter                                   :: prolong_order = 1

        ! Relaxation params
        real(dp), parameter                                  :: relax_tol = -one
        real(dp), parameter                                  :: relax_omega = four/five
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

            ! Set the initial guess
            ! e(depth)%data = zero
            e(depth)%data(ilo:ihi,jlo:jhi) = r(depth)%data(ilo:ihi,jlo:jhi) * invdiags(depth)%data(ilo:ihi,jlo:jhi)

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

            ! Set the initial guess
            ! e(depth)%data = zero
            e(depth)%data(ilo:ihi,jlo:jhi) = r(depth)%data(ilo:ihi,jlo:jhi) * invdiags(depth)%data(ilo:ihi,jlo:jhi)

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


end module DeferredMGPoisson2D


! NOTES ............................
! 1. inner_prod should scale by J
! 2. check if restrict or prolong need J scaling too.
! 3. Test compute_pd
! 4. Refinement ratio is hard-coded at 4,4 in fill_ghosts.
! 5. If rscale = zero, manually set to one.
! 6. Handle beta = zero in bicgstab.
