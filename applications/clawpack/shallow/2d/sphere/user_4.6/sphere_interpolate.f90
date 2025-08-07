!! ----------------------------------------------------------
!!  @file
!!  Interpolation routines - (i,j,mq) ordering
!! ----------------------------------------------------------
!! interpolate_face_ghost
!! interpolate_corner_ghost
!! interpolate_to_fine_patch
!! Other routines :
!! fclaw2d_clawpatch_compute_slopes 
!! (for limited function reconstruction)
!! fixcapaq (to preserve conservation)
!! Note that fixcapaq is only used when regridding;  ghost
!! cell interpolation is not conservative in the mapped case.
!! (Should it be?  We are going to correct the flux mixmatch
!! anyhow, so maybe the accuracy of the ghost cell values is
!! more important.)
!! ----------------------------------------------------------


!! -------------------------------------------------------------------
!!  @brief @copybrief ::clawpatch_fort_interpolate_face_t
!! 
!!  Implementation for clawpack 4.6.
!!  This routine is used for both mapped and non-mapped cases.
!! 
!!  @details @copydetails ::clawpatch_fort_interpolate_face_t
!! -------------------------------------------------------------------
subroutine sphere_fort_interpolate_face & 
           (mx,my,mbc,meqn,qcoarse,qfine,  & 
            auxcoarse, auxfine, maux, mbathy, & 
            idir,iface_coarse, igrid, transform_ptr)
    implicit none
    integer mx,my,mbc,meqn,maux, igrid,idir,iface_coarse, mbathy
    integer*8 transform_ptr
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer mq,m
    integer ibc, i1
    integer jbc, j1
    integer ic, jc, mth
    double precision gradx, grady, qc, sl, sr, b, value
    double precision fclaw2d_clawpatch_compute_slopes

    !! # This should be refratio*refratio.
    integer i2(0:3),j2(0:3)
    logical fclaw2d_clawpatch_is_valid_interp
    logical skip_this_grid

    integer a(2,2), f(2)
    integer ii,jj,dc(2),df(2,0:3),iff,jff
    double precision shiftx(0:3),shifty(0:3)
    double precision qc_quad(-1:1,-1:1,meqn)

    integer sweep_dir

    mth = 5

    call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)

    !! # This needs to be written for refratios .ne. 2.
    m = 0
    do jj = 0,1
        do ii = 0,1
            !! # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

            !! # Direction on fine grid (converted using metric). Divide
            !! # by refratio=2 to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/2
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/2

            !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5)/2
            shifty(m) = (jj-0.5)/2
            m = m + 1
        end do
    end do
    !! Create map :


    do mq = 1,meqn
       if (idir .eq. 0) then
            !! # this ensures that we get 'hanging' corners

            sweep_dir = 0
            !! call sphere_fix_corners(mx,my,mbc,meqn,qcoarse, sweep_dir)

            do ibc = 1,mbc/2
                if (iface_coarse .eq. 0) then
                    ic = ibc
                elseif (iface_coarse .eq. 1) then
                    ic = mx - ibc + 1
                else
                    write(6,*) 'interpolate : Problem with iface_coarse'
                    write(6,*) 'iface_coarse = ', iface_coarse
                    stop               
                endif
                do jc = 1,mx
                    i1 = ic
                    j1 = jc
                    call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2, transform_ptr)
                    skip_this_grid = .false.
                    do m = 0,3
                        if (.not. fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc)) then
                             skip_this_grid = .true.
                             return
                         endif
                     end do
                     if (.not. skip_this_grid) then
                         !! Calculate surface elevation eta using dry limiting
                         do ii = -1, 1
                             do jj = -1, 1
                                 qc_quad(ii,jj,mq) = qcoarse(ic+ii,jc+jj,mq)
                                 if (mq .eq. 1) then
                                     b = auxcoarse(ic+ii,jc+jj,mbathy)
                                    qc_quad(ii,jj,mq) = qc_quad(ii,jj,mq) + b
                                endif 
                            end do
                        end do

                        qc = qc_quad(0,0,mq)

                        sl = qc - qc_quad(-1,0,mq)
                        sr = qc_quad(1,0,mq) - qc
                        gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                        sl = qc - qc_quad(0,-1,mq)
                        sr = qc_quad(0,1,mq) - qc
                        grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)
                        do m = 0,3
                            iff = i2(0) + df(1,m)
                            jff = j2(0) + df(2,m)
                            value = qc + gradx*shiftx(m) + grady*shifty(m)
                            if (mq .eq. 1) then 
                                value = max(value-auxfine(iff,jff,mbathy),0.0)
                            endif
                            qfine(iff,jff,mq) = value
                        end do
                    endif  !! end skip grid
                end do
            end do
        else
            !! idir = 1
            sweep_dir = 1
            !!call sphere_fix_corners(mx,my,mbc,meqn,qcoarse,sweep_dir)

            do jbc = 1,mbc/2
                if (iface_coarse .eq. 2) then
                    jc = jbc
                elseif (iface_coarse .eq. 3) then
                    !! iface_coarse = 3
                    jc = my - jbc + 1
                else
                    write(6,*) 'interpolate : Problem with iface_coarse'
                    write(6,*) 'iface_coarse = ', iface_coarse
                    stop
                endif
                do ic = 1,mx
                    i1 = ic
                    j1 = jc
                    call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2,transform_ptr)
                    !!  # ---------------------------------------------
                    !!  # Two 'half-size' neighbors will be passed into
                    !!  # this routine.  Only half of the coarse grid ghost
                    !!  # indices will be valid for the particular grid
                    !!  # passed in.  We skip those ghost cells that will
                    !!  # have to be filled in by the other half-size
                    !!  # grid.
                    !!  # ---------------------------------------------
                    skip_this_grid = .false.
                    do m = 0,1
                        if (.not. fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc)) then
                            skip_this_grid = .true.
                            return
                        endif
                     end do
                     if (.not. skip_this_grid) then
                         !! Calculate surface elevation eta using dry limiting
                         do ii = -1, 1
                             do jj = -1, 1
                                qc_quad(ii,jj,mq) = qcoarse(ic+ii,jc+jj,mq)
                                if (mq .eq. 1) then
                                    b = auxcoarse(ic+ii,jc+jj,mbathy)
                                    qc_quad(ii,jj,mq) = qc_quad(ii,jj,mq) + b
                                endif 
                            end do
                        end do

                        qc = qc_quad(0,0,mq)

                        sl = (qc - qc_quad(-1,0,mq))
                        sr = (qc_quad(1,0,mq) - qc)
                        gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                        sl = (qc - qc_quad(0,-1,mq))
                        sr = (qc_quad(0,1,mq) - qc)
                        grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                        do m = 0,3
                            iff = i2(0) + df(1,m)
                            jff = j2(0) + df(2,m)
                            value = qc + gradx*shiftx(m) + grady*shifty(m)
                            if (mq .eq. 1) then 
                                value = max(value-auxfine(iff,jff,mbathy),0.0)
                            endif
                            qfine(iff,jff,mq) = value
                       end do
                    endif                    !! Don't skip this grid
                end do                       !! i loop
            end do                       !! end of jbc loop
        endif                          !! end idir branch
    end do                             !! endo mq loop

end subroutine sphere_fort_interpolate_face

!! # -------------------------------------------------------------------
!! #  @brief @copybrief ::clawpatch_fort_interpolate_corner_t
!! # 
!! #  Implementation for clawpack 4.6.
!! # 
!! #  @details @copydetails ::clawpatch_fort_interpolate_corner_t
!! # -------------------------------------------------------------------
subroutine sphere_fort_interpolate_corner(mx,my,mbc,meqn,& 
     qcoarse,qfine,auxcoarse, auxfine, maux, mbathy, & 
     icorner_coarse,transform_ptr)
    implicit none

    integer mx,my,mbc,meqn,icorner_coarse,mbathy, maux
    integer*8 transform_ptr
    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer ic, jc, mq, ibc,jbc, mth
    double precision qc, sl, sr, gradx, grady
    double precision fclaw2d_clawpatch_compute_slopes, value
    double precision qc_quad(-1:1,-1:1,meqn), b

    !! # This should be refratio*refratio.
    integer i1,j1,m
    integer i2(0:3),j2(0:3)

    integer a(2,2), f(2)
    integer ii,jj,iff,jff,dc(2),df(2,0:3)
    double precision shiftx(0:3), shifty(0:3)

    call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)

    m = 0
    do jj = 0,1
        do ii = 0,1
            !! # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

            !! # Direction on fine grid (converted using metric). Divide
            !! # by 2 (refratio) to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/2
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/2

            !! # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = (ii-0.5d0)/2.d0
            shifty(m) = (jj-0.5d0)/2.d0
            m = m + 1
        end do
    end do


    mth = 5

    do ibc = 1,mbc/2
        do jbc = 1,mbc/2
            if (icorner_coarse .eq. 0) then
               ic = ibc
               jc = jbc
            elseif (icorner_coarse .eq. 1) then
               ic = mx - ibc + 1
               jc = jbc
            elseif (icorner_coarse .eq. 2) then
               ic = ibc
               jc = my - jbc + 1
            elseif (icorner_coarse .eq. 3) then
               ic = mx - ibc + 1
               jc = my - jbc + 1
           else
               write(6,*) "interpolate : Problem with icorner_coarse"
               write(6,*) "icorner_coarse = ", icorner_coarse
               stop
           endif

            !! # Interpolate coarse grid corners to fine grid corner ghost cells
            i1 = ic
            j1 = jc
            call fclaw2d_clawpatch_transform_corner_half(i1,j1,i2,j2,transform_ptr)

            do mq = 1,meqn
                !! Calculate surface elevation eta using dry limiting
                do ii = -1, 1
                    do jj = -1, 1
                        qc_quad(ii,jj,mq) = qcoarse(ic+ii,jc+jj,mq)
                        if (mq .eq. 1) then
                            b = auxcoarse(ic+ii,jc+jj,mbathy)
                            qc_quad(ii,jj,mq) = qc_quad(ii,jj,mq) + b
                        endif 
                    end do
                end do
                qc = qc_quad(0,0,mq)

                !! # Compute limited slopes in both x and y. Note we are not
                !! # really computing slopes, but rather just differences.
                !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                sl = (qc - qc_quad(-1,0,mq))
                sr = (qc_quad(1,0,mq) - qc)
                gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                sl = (qc - qc_quad(0,-1,mq))
                sr = (qc_quad(0,1,mq) - qc)
                grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                do m = 0,3
                    iff = i2(0) + df(1,m)
                    jff = j2(0) + df(2,m)
                    value = qc + gradx*shiftx(m) + grady*shifty(m)
                    if (mq .eq. 1) then 
                        value = max(value-auxfine(iff,jff,mbathy),0.0)
                    endif
                    qfine(iff,jff,mq) = value
                end do
            end do
        end do 
    end do
end subroutine sphere_fort_interpolate_corner
 

!! # -------------------------------------------------------------------
!! #  @brief @copybrief ::clawpatch_fort_interpolate2fine_t
!! # 
!! #  Implementation for clawpack 4.6.
!! # 
!! #  @details @copydetails ::clawpatch_fort_interpolate2fine_t
!! # -------------------------------------------------------------------
subroutine sphere_fort_interpolate2fine & 
          (mx,my,mbc,meqn,qcoarse, qfine, areacoarse,  & 
           areafine, igrid, manifold)
    implicit none

    integer mx,my,mbc,meqn
    integer igrid, manifold

    double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer ii, jj, i,j, i1, i2, j1, j2, ig, jg, mq, mth
    integer ic,jc,ic_add, jc_add
    double precision qc, shiftx, shifty, sl, sr, gradx, grady
    double precision fclaw2d_clawpatch_compute_slopes

    integer p4est_refineFactor,refratio

    p4est_refineFactor = 2
    refratio = 2

    !! # Use limiting done in AMRClaw.
    mth = 5

    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    i1 = 1-ig
    i2 = mx/p4est_refineFactor + (1-ig)
    ic_add = ig*mx/p4est_refineFactor

    j1 = 1-jg
    j2 = my/p4est_refineFactor + (1-jg)
    jc_add = jg*my/p4est_refineFactor

    do mq = 1,meqn
        do j = j1,j2
            do i = i1,i2
                ic = i + ic_add
                jc = j + jc_add
                qc = qcoarse(ic,jc,mq)

                !! # Compute limited slopes in both x and y. Note we are not
                !! # really computing slopes, but rather just differences.
                !! # Scaling is accounted for in 'shiftx' and 'shifty', below.
                sl = (qc - qcoarse(ic-1,jc,mq))
                sr = (qcoarse(ic+1,jc,mq) - qc)
                gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                sl = (qc - qcoarse(ic,jc-1,mq))
                sr = (qcoarse(ic,jc+1,mq) - qc)
                grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                !! # Fill in refined values on coarse grid cell (ic,jc)
                do ii = 1,refratio
                    do jj = 1,refratio
                        shiftx = (ii - refratio/2.d0 - 0.5d0)/refratio
                        shifty = (jj - refratio/2.d0 - 0.5d0)/refratio
                        qfine((i-1)*refratio + ii,(j-1)*refratio + jj,mq) & 
                        = qc + shiftx*gradx + shifty*grady
                    end do
                end do
            end do
        end do
    end do

    if (manifold .ne. 0) then
!!         call sphere_fort_fixcapaq2(mx,my,mbc,meqn,
!!     &         qcoarse,qfine, areacoarse,areafine,igrid)
    endif


end subroutine sphere_fort_interpolate2fine


!! logical function fclaw2d_clawpatch_is_valid_interp(i,j,mx,my,mbc)
!!     implicit none
!!     integer i,j,mx, my, mbc
!! 
!!     logical i1, j1
!! 
!!     i1 = 1-mbc .le. i .and. i .le. mx+mbc
!!     j1 = 1-mbc .le. j .and. j .le. my+mbc
!! 
!!     sphere_is_valid_interp = i1 .and. j1
!! 
!! end function sphere_is_valid_interp


!! # -------------------------------------------------------------------
!! #  @brief ensures interpolated fine grid has same mass as coarse grid
!! # 
!! #  So far, this is only used by the interpolation from
!! #  coarse to fine when regridding.  But maybe it should
!! #  be used by the ghost cell routines as well?
!! # 
!! #  @param[in] mx, my the number of cells in the x and y directions
!! #  @param[in] mbc the number of ghost cells
!! #  @param[in] meqn the number of equations
!! #  @param[in] qcoarse the coarse solution
!! #  @param[in,out] qfine the interpolated fine solution, update for 
!! #                 conservation
!! #  @param[in] areacoarse, areafine the areas of the fine and coarse grids
!! #  @param[in] igrid the index of the fine grid in the child array
!! # -------------------------------------------------------------------
subroutine sphere_fort_fixcapaq2(mx,my,mbc,meqn, & 
           qcoarse,qfine, areacoarse,areafine,igrid)
  implicit none

  integer mx,my,mbc,meqn, refratio, igrid
  integer p4est_refineFactor

  double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
  double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
  double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
  double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

  integer i,j,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
  double precision kf, kc, r2, sum, cons_diff, qf, qc


  p4est_refineFactor = 2
  refratio = 2

    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/p4est_refineFactor
    jc_add = jg*my/p4est_refineFactor

    !! # ------------------------------------------------------
    !! # This routine ensures that the interpolated solution
    !! # has the same mass as the coarse grid solution
    !! # -------------------------------------------------------


    r2 = refratio*refratio
    do m = 1,meqn
        do i = 1,mx/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
                sum = 0.d0
                do ii = 1,refratio
                    do jj = 1,refratio
                        ifine = (i-1)*refratio + ii
                        jfine = (j-1)*refratio + jj
                        kf = areafine(ifine,jfine)
                        qf = qfine(ifine,jfine,m)
                        sum = sum + kf*qf
                    end do
                end do

                kc = areacoarse(i+ic_add,j+jc_add)
                qc = qcoarse(i+ic_add, j+jc_add,m)
                cons_diff = (qc*kc - sum)/r2

               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine  = (i-1)*refratio + ii
                     jfine  = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     qfine(ifine,jfine,m) = qfine(ifine,jfine,m) + cons_diff/kf
                  end do
               end do
            end do  !! end of meqn
        end do
    end do

end subroutine sphere_fort_fixcapaq2


subroutine sphere_set_corner_count(corner_count)
    implicit none
    integer corner_count(0:3)

    integer block_corner_count(0:3)
    common /com_interp/ block_corner_count

    integer k

    do k = 0,3
        block_corner_count(k) = corner_count(k)
    end do

end subroutine sphere_set_corner_count


!! #  See 'cubed_sphere_corners.ipynb'
subroutine sphere_fix_corners(mx,my,mbc,meqn,q, sweep_dir)
    implicit none

    integer :: mx,my,mbc,meqn,sweep_dir
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    !!double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer :: k,m,idata,jdata
    double precision :: ihat(0:3),jhat(0:3)
    integer :: i1, j1, ibc, jbc
    logical :: use_b

    integer block_corner_count(0:3)
    common /com_interp/ block_corner_count

    !! # Lower left corner
    ihat(0) = 0.5
    jhat(0) = 0.5

    !! # Lower right corner
    ihat(1) = mx+0.5
    jhat(1) = 0.5

    !! # Upper left corner
    ihat(2) = 0.5
    jhat(2) = my+0.5

    !! # Upper right corner
    ihat(3) = mx+0.5
    jhat(3) = my+0.5

    do k = 0,3
        if (block_corner_count(k) .ne. 3) then
            cycle
        endif
        use_b = sweep_dir .eq. 0 .and. (k .eq. 0 .or. k .eq. 3) & 
             .or.  sweep_dir .eq. 1 .and. (k .eq. 1 .or. k .eq. 2)
        do ibc = 1,mbc
            do jbc = 1,mbc
                !! # Average fine grid corners onto coarse grid ghost corners
                if (k .eq. 0) then
                    i1 = 1-ibc
                    j1 = 1-jbc
                elseif (k .eq. 1) then
                    i1 = mx+ibc
                    j1 = 1-jbc
                elseif (k .eq. 2) then
                    i1 = 1-ibc
                    j1 = my+jbc
                elseif (k .eq. 3) then
                    i1 = mx+ibc
                    j1 = my+jbc
                endif

                if (use_b) then
                    !! # Transform involves B                
                    idata =  j1 + int(ihat(k) - jhat(k))
                    jdata = -i1 + int(ihat(k) + jhat(k))
                else
                    !! # Transform involves B.transpose()             
                    idata = -j1 + int(ihat(k) + jhat(k))
                    jdata =  i1 - int(ihat(k) - jhat(k))
                endif 
                do m = 1,meqn
                    q(i1,j1,m) = q(idata,jdata,m)
                end do           
!!               do m = 1,maux
!!                  aux(i1,j1,m) = aux(idata,jdata,m)
!!               end do           
            end do              !! ibc
        end do                 !! jbc
    end do                    !! corner 'k' loop
end subroutine sphere_fix_corners

