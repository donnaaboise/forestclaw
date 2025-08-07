!!  \file
!! 
!!  Routines described here average are used to fill coarse grid ghost
!!  cells, and average sibling grids onto a parent grid.  Indices
!!  for cells at block boundaries are transformed using encodings
!!  stored in `transform_cptr`.
!! 
!! -------------------------------------------------------------------
!!  @brief @copybrief ::clawpatch_fort_average_face_t
!! 
!!  Implementation for clawpack 4.6
!! 
!!  @details @copydetails ::clawpatch_fort_average_face_t
!! -------------------------------------------------------------------
subroutine sphere_fort_average_face(mx,my,mbc,meqn, mcapa, mbathy,&
           qcoarse,qfine,auxcoarse, auxfine, maux, &
           idir,iface_coarse,igrid, transform_cptr)
    implicit none

    integer mx,my,mbc,meqn,igrid,idir,iface_coarse, maux
    integer mbathy, mcapa

    integer*8 transform_cptr
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer ic, ibc
    integer jc, jbc
    integer m

    double precision capac, etasum, hsum, husum, hvsum, hwsum
    double precision bf, hf, huf, hvf, hwf, etaf, capa
    double precision etaavg, havg, hc, huc, hvc, hwc

    integer i2(0:3),j2(0:3)

    logical fclaw2d_clawpatch_is_valid_average, skip_this_grid

    !! # 'iface' is relative to the coarse grid

    !! # Average fine grid onto coarse grid
    if (idir .eq. 0) then
        do jc = 1,my
            do ibc = 1,mbc
                !! # ibc = 1 corresponds to first layer of ghost cells, and
                !! # ibc = 2 corresponds to the second layer

                if (iface_coarse .eq. 0) then
                    ic = 1-ibc
                elseif (iface_coarse .eq. 1) then
                    ic = mx+ibc
                endif

                call  fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2, transform_cptr)
                !! # ---------------------------------------------
                !! # Two 'half-size' neighbors will be passed into
                !! # this routine.  Only half of the coarse grid ghost
                !! # indices will be valid for the particular grid
                !! # passed in.  We skip those ghost cells that will
                !! # have to be filled in by the other half-size
                !! # grid.
                !! # ---------------------------------------------
                skip_this_grid = .false.
                do m = 0,3
                    if (.not. fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my)) then
                        skip_this_grid = .true.
                        return
                    endif
                end do

                if (.not. skip_this_grid) then
                    capac = auxcoarse(ic,jc,mcapa)

                    etasum = 0
                    hsum   = 0
                    husum  = 0
                    hvsum  = 0
                    hwsum  = 0

                    do m = 0,3
                        capa = auxfine(i2(m),j2(m),mcapa)
                        bf = auxfine(i2(m),j2(m),mbathy)*capa

                        hf   = qfine(i2(m),j2(m),1)*capa
                        huf  = qfine(i2(m),j2(m),2)*capa
                        hvf  = qfine(i2(m),j2(m),3)*capa
                        hwf  = qfine(i2(m),j2(m),4)*capa
                        etaf = hf+bf

                        hsum   = hsum + hf
                        husum  = husum + huf
                        hvsum  = hvsum + hvf
                        hwsum  = hwsum + hwf
                        etasum = etasum + etaf
                    end do
                    etaavg = etasum/4
                    havg = hsum/4
                    hc = min(havg,(max(etaavg - auxcoarse(ic,jc,mbathy)*capac,0.d0)))
                    huc = (hc/hsum)*husum
                    hvc = (hc/hsum)*hvsum
                    hwc = (hc/hsum)*hwsum

                    qcoarse(ic,jc,1) = hc / capac
                    qcoarse(ic,jc,2) = huc / capac
                    qcoarse(ic,jc,3) = hvc / capac
                    qcoarse(ic,jc,4) = hwc / capac
                endif
            enddo
        enddo
    else
        !! dir == 1
        !! # idir = 1 (faces 2,3)
        do jbc = 1,mbc
            do ic = 1,mx

                if (iface_coarse .eq. 2) then
                    jc = 1-jbc
                else if (iface_coarse .eq. 3) then
                    jc = my+jbc
                endif

                call fclaw2d_clawpatch_transform_face_half(ic,jc,i2,j2, transform_cptr)
                skip_this_grid = .false.
                do m = 0,3
                    if (.not. fclaw2d_clawpatch_is_valid_average(i2(m),j2(m),mx,my)) then
                        skip_this_grid = .true.
                    endif
                end do
                if (.not. skip_this_grid) then
                    capac = auxcoarse(ic,jc,mcapa)

                    etasum = 0
                    hsum   = 0
                    husum  = 0
                    hvsum  = 0

                    do m = 0,3
                        capa = auxfine(i2(m),j2(m),mcapa)
                        bf = auxfine(i2(m),j2(m),mbathy)*capa

                        hf = qfine(i2(m),j2(m),1)*capa
                        huf= qfine(i2(m),j2(m),2)*capa
                        hvf= qfine(i2(m),j2(m),3)*capa
                        hwf= qfine(i2(m),j2(m),4)*capa
                        etaf = hf+bf

                        hsum   = hsum + hf
                        husum  = husum + huf
                        hvsum  = hvsum + hvf
                        hwsum = hwsum + hwf
                        etasum = etasum + etaf
                    end do
                    etaavg = etasum/4
                    havg = hsum/4
                    hc = min(havg,(max(etaavg- & 
                            auxcoarse(ic,jc,mbathy)*capac,0.d0)))
                    huc = (hc/hsum)*husum
                    hvc = (hc/hsum)*hvsum
                    hwc = (hc/hsum)*hwsum

                    qcoarse(ic,jc,1) = hc / capac
                    qcoarse(ic,jc,2) = huc / capac
                    qcoarse(ic,jc,3) = hvc / capac
                    qcoarse(ic,jc,4) = hwc / capac
                endif      !! skip grid loop
            end do
        end do
    endif
    
end subroutine sphere_fort_average_face


!! -------------------------------------------------------------------
!!  @brief @copybrief ::clawpatch_fort_average_corner_t
!! 
!!  Implementation for clawpack 4.6
!! 
!!  @details @copydetails ::clawpatch_fort_average_corner_t
!! -------------------------------------------------------------------
subroutine sphere_fort_average_corner(mx,my,mbc,meqn, mcapa, mbathy, &
           qcoarse,qfine,auxcoarse, auxfine, maux, &
           icorner_coarse,transform_cptr)
    implicit none

    integer mx,my,mbc,meqn,icorner_coarse,maux
    integer mbathy, mcapa

    integer*8 transform_cptr
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    double precision capac, etasum, hsum, husum, hvsum, hwsum
    double precision bf, hf, huf, hvf, hwf, etaf, capa
    double precision etaavg, havg, hc, huc, hvc, hwc

    !! # This should be refratio*refratio.
    integer i1,j1,m, ibc, jbc
    integer i2(0:3),j2(0:3)

    !! # Loop over four corner cells on coarse grid
    do ibc = 1,mbc
        do jbc = 1,mbc
            !! # Average fine grid corners onto coarse grid ghost corners
            if (icorner_coarse .eq. 0) then
               i1 = 1-ibc
               j1 = 1-jbc
            elseif (icorner_coarse .eq. 1) then
               i1 = mx+ibc
               j1 = 1-jbc
            elseif (icorner_coarse .eq. 2) then
               i1 = 1-ibc
               j1 = my+jbc
            elseif (icorner_coarse .eq. 3) then
               i1 = mx+ibc
               j1 = my+jbc
            endif

            !! # Again, a fake routine until the real one is
            !! # available (be sure to pass in (i1,j1)
            call fclaw2d_clawpatch_transform_corner_half(i1,j1,i2,j2,transform_cptr)


            capac = auxcoarse(i1,j1,mcapa)
            etasum = 0.d0
            hsum   = 0.d0
            husum  = 0.d0
            hvsum  = 0.d0

            do m = 0,3
                capa = auxfine(i2(m),j2(m),mcapa)
                bf = auxfine(i2(m),j2(m),mbathy)*capa

                hf  = qfine(i2(m),j2(m),1)*capa
                huf = qfine(i2(m),j2(m),2)*capa
                hvf = qfine(i2(m),j2(m),3)*capa
                hwf = qfine(i2(m),j2(m),4)*capa
                etaf = hf + bf

                hsum   = hsum + hf
                husum  = husum + huf
                hvsum  = hvsum + hvf
                hwsum  = hwsum + hwf
                etasum = etasum + etaf
            end do
            etaavg = etasum/4
            havg = hsum/4
            hc = min(havg,(max(etaavg - auxcoarse(i1,j1,mbathy)*capac,0.d0)))
            huc = (hc/hsum)*husum
            hvc = (hc/hsum)*hvsum
            hwc = (hc/hsum)*hwsum

            qcoarse(i1,j1,1) = hc / capac
            qcoarse(i1,j1,2) = huc / capac
            qcoarse(i1,j1,3) = hvc / capac
            qcoarse(i1,j1,4) = hwc / capac            
        end do
    end do

end subroutine sphere_fort_average_corner


!! -------------------------------------------------------------------
!!  @brief @copybrief ::clawpatch_fort_average2coarse_t
!! 
!!  Implementation for clawpack 4.6
!! 
!!  @details @copydetails ::clawpatch_fort_average2coarse_t
!! -------------------------------------------------------------------
subroutine sphere_fort_average2coarse(mx,my,mbc,meqn, mcapa,mbathy,& 
           qcoarse,qfine, auxcoarse, auxfine, maux, igrid)
    implicit none

    integer mx,my,mbc,meqn,igrid, maux,mcapa,mbathy
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer i,j, ig, jg, ic_add, jc_add, ii, jj

    double precision capac, etasum, hsum, husum, hvsum, hwsum
    double precision bf, hf, huf, hvf, hwf, etaf, capa
    double precision etaavg, havg, hc, huc, hvc, hwc


    !! # This should be refratio*refratio.
    integer i1,j1, m
    integer i2(0:3),j2(0:3)

    !! # 'iface' is relative to the coarse grid

    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,2)
    jg = (igrid-ig)/2

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/2
    jc_add = jg*mx/2

    do j = 1,my/2
        do i = 1,mx/2
            i1 = i + ic_add
            j1 = j + jc_add
            m = 0
            do jj = 1,2
                do ii = 1,2
                    i2(m) = (i-1)*2 + ii
                    j2(m) = (j-1)*2 + jj
                    m = m + 1
                end do
            end do
            etasum = 0
            hsum   = 0
            husum  = 0
            hvsum  = 0
            hwsum  = 0

            do m = 0,3
                capa = auxfine(i2(m),j2(m),mcapa)
                bf = auxfine(i2(m),j2(m),mbathy)*capa

                hf   = qfine(i2(m),j2(m),1)*capa
                huf  = qfine(i2(m),j2(m),2)*capa
                hvf  = qfine(i2(m),j2(m),3)*capa
                hwf  = qfine(i2(m),j2(m),4)*capa
                etaf = hf+bf

                hsum   = hsum + hf
                husum  = husum + huf
                hvsum  = hvsum + hvf
                hwsum = hwsum + hwf
                etasum = etasum + etaf
            end do
            capac = auxcoarse(i1,j1,mcapa)

            etaavg = etasum/4
            havg = hsum/4
            hc = min(havg,(max(etaavg - auxcoarse(i1,j1,mbathy)*capac,0.d0)))
            huc = (hc/hsum)*husum
            hvc = (hc/hsum)*hvsum
            hwc = (hc/hsum)*hwsum

            qcoarse(i1,j1,1) = hc / capac
            qcoarse(i1,j1,2) = huc / capac
            qcoarse(i1,j1,3) = hvc / capac
            qcoarse(i1,j1,4) = hwc / capac
        end do
    end do
end subroutine sphere_fort_average2coarse


