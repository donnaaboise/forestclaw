double precision function fdisc(blockno,xc,yc,zc)
    implicit none

    double precision xc,yc, zc
    integer blockno

    double precision pi, pi2, deg2rad
    common /compi/ pi, pi2, deg2rad

    integer init_cond
    common /latlong_initcond/ init_cond

    double precision ring_inner, ring_outer, center(2)
    common /latlong_initcond_parms2/ ring_inner, ring_outer, center

    integer ring_units
    common /latlong_initcond_parms3/ ring_units

    integer*8 cont, fclaw_map_get_context

    double precision xp, yp, zp, rp, a, m, d, c(3), d1
    double precision phi, ri, ro, th, tp, xcc, ycc, zcc
    integer i

    cont = fclaw_map_get_context()

    call fclaw_map_3d_c2m(cont,blockno,xc,yc,zc,xp,yp,zp)

    rp = sqrt(xp**2 + yp**2 + zp**2)

    if (init_cond .eq. 0) then
        th = atan2(yp,xp)
        tp = abs(th)
        fdisc = tp - pi/4.d0
    else
        ri = ring_inner
        ro = ring_outer
        if (ring_units .gt. 0) then
            if (ring_units .eq. 1) then
                !! Convert degrees to radians
                ri = ri*deg2rad
                ro = ro*deg2rad
            elseif (ring_units .eq. 2) then
                !! Compute Gaussian distance to radians
                ri = ri/rp
                ro = ro/rp
            endif
        endif
        !! Initialize smoke in the center of lat/long region
        xcc = 0.5d0
        ycc = 0.5d0
        zcc = zc
        call fclaw_map_3d_c2m(cont,blockno,xcc,ycc,zcc,c(1),c(2),c(3))        
        d1 = sqrt(c(1)**2 + c(2)**2 + c(3)**2)
        do i = 1,3
            c(i) = c(i)/d1
        end do
        d = c(1)*xp + c(2)*yp + c(3)*zp
        if (d/rp .gt. 1) then
            if (d/rp - 1 .lt. 1e-10) then
                d = rp
            else
                write(6,*) 'd/rp = ', d/rp
                write(6,*) 'fdisc : problem taking acos()'
                stop
            endif
        endif
        phi = acos(d/rp)  
        a = abs(ro - ri)/2.0   !! Half-width of the ring
        m = (ro + ri)/2.0   !! Midpoint of the ring

        !! Negative is "inside" the domain (w = 1)
        fdisc = abs(phi-m) - a
    endif
end



