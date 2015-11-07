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
