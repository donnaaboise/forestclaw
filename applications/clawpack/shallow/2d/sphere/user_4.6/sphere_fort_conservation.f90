subroutine sphere_fort_conservation_check  & 
          (mx,my,mbc,mfields,dx,dy,area,q,sum,c_kahan)
    implicit none

    integer mx,my,mbc,mfields
    double precision dx, dy
    double precision sum(mfields), c_kahan(mfields)
    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision t, y, area_ij
    double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer i,j,m
    integer*8 cont, fclaw_map_get_context
    logical fclaw_map_is_used

    logical use_kahan

    use_kahan = .true.

    cont = fclaw_map_get_context()

    area_ij = dx*dy  !! Area in each mesh cell is constant
    do m = 1,mfields
        do j = 1,my
            do i = 1,mx
                if (fclaw_map_is_used(cont)) then
                    area_ij = area(i,j)  !! Area varies
                endif
                if (use_kahan) then
                    y = q(i,j,m)*area_ij - c_kahan(m)
                    t = sum(m) + y
                    c_kahan(m) = (t-sum(m)) - y
                    sum(m) = t
                else
                    sum(m) = sum(m) + q(i,j,m)*area_ij
                endif
            enddo
        enddo
    enddo

end subroutine sphere_fort_conservation_check