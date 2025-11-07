subroutine sphere_fort_conservation_check(blockno,  & 
          mx,my,mbc,mfields,xlower, ylower, & 
          dx,dy,area,q,sum,c_kahan)
    implicit none

    integer blockno, mx,my,mbc,mfields
    double precision xlower, ylower, dx, dy
    double precision sum(mfields), c_kahan(mfields)
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision t, y, area_ij
    double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer i,j,m
    double precision xc, yc, hu(3), cq(mfields), xvec(3)

    integer*8 cont, fclaw_map_get_context

    logical use_kahan, angular

    use_kahan = .true.

    angular = .true.

    cont = fclaw_map_get_context()

    do j = 1,my
        do i = 1,mx
            if (angular) then
                cq(1) = q(i,j,1)

                !! Compute angular momentum
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                call fclaw_map_2d_c2m(cont,blockno,xc,yc,xvec(1),xvec(2),xvec(3))

                do m = 1,3
                    hu(m) = q(i,j,m+1)
                end do
                call sphere_compute_cross(xvec,hu,cq(2:4))
            else
                do m = 1,mfields
                    cq(m) = q(i,j,m)
                end do
            endif

            area_ij = area(i,j)  !! Area varies
            do m = 1,mfields
                if (use_kahan) then
                    y = cq(m)*area_ij - c_kahan(m)
                    t = sum(m) + y
                    c_kahan(m) = (t-sum(m)) - y
                    sum(m) = t
                else
                    sum(m) = sum(m) + cq(m)*area_ij
                endif
            end do
        enddo
    enddo

end subroutine sphere_fort_conservation_check

subroutine sphere_compute_cross(u,v,w)
    implicit none

    double precision u(3), v(3), w(3)

    !! From sympy
    !! (u_{2}*v_{3} - u_{3}*v_{2})*R.i + 
    !! (-u_{1}*v_{3} + u_{3}*v_{1})*R.j + 
    !! (u_{1}*v_{2} - u_{2}*v_{1})*R.k

    w(1) =  u(2)*v(3) - u(3)*v(2)
    w(2) = -u(1)*v(3) + u(3)*v(1)
    w(3) =  u(1)*v(2) - u(2)*v(1)

end subroutine sphere_compute_cross