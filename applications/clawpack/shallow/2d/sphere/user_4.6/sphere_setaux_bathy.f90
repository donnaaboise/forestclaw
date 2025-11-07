!! The aux array has the following elements:
!! 1  kappa = ratio of cell area to dxc*dyc
!! 2  enx = x-component of normal vector to left edge in tangent plane
!! 3  eny = y-component of normal vector to left edge in tangent plane
!! 4  enz = z-component of normal vector to left edge in tangent plane
!! 5  etx = x-component of tangent vector to left edge in tangent plane
!! 6  ety = y-component of tangent vector to left edge in tangent plane
!! 7  etz = z-component of tangent vector to left edge in tangent plane
!! 8  enx = x-component of normal vector to bottom edge in tangent plane
!! 9  eny = y-component of normal vector to bottom edge in tangent plane
!! 10  enz = z-component of normal vector to bottom edge in tangent plane
!! 11  etx = x-component of tangent vector to bottom edge in tangent plane
!! 12  ety = y-component of tangent vector to bottom edge in tangent plane
!! 13  etz = z-component of tangent vector to bottom edge in tangent plane
!! 14  erx = x-component of unit vector in radial direction at cell ctr
!! 15  ery = y-component of unit vector in radial direction at cell ctr
!! 16  erz = z-component of unit vector in radial direction at cell ctr

!! (not used here)
!! 17  xc at cell center
!! 18  yc at cell center
!! 19  bathymetry - averaged over all possible finer cells


subroutine sphere_setaux(blockno, mx,my,mbc,xlower,ylower, & 
    dx,dy,area,xnormals,ynormals, & 
    xtangents,ytangents,edgelengths, surfnormals,aux,maux)

    use lambwave_module, only : integrate_bathy
      
    implicit none
    integer mbc, mx,my, maux, blockno
    double precision xlower, ylower, dx,dy
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    !! # Edge based quantities
    double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
    double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
    double precision    xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
    double precision    ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
    double precision  edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
   
    double precision eta_in, eta_out, height_top, height_temp
    common /swe_initcond_parms4/ eta_in, eta_out, height_top, height_temp

    !! Cell centered quantities
    double precision        area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
    double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

    integer i,j,m
    double precision dxdy

    dxdy = dx*dy
    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc

            !! (1) Capacity (mcapa)
            aux(i,j,1) = area(i,j)/dxdy
            do m = 1,3
                !! (2-7) normals and tangents at left x-face
                aux(i,j,1+m) = xnormals(i,j,m)
                aux(i,j,4+m) = xtangents(i,j,m)

                !! (8-13) normals and tangents at bottom y-face
                aux(i,j,7+m) = ynormals(i,j,m)
                aux(i,j,10+m) = ytangents(i,j,m)

                !! (14-16) surface normal at cell centers
                aux(i,j,13+m) = surfnormals(i,j,m)

            enddo

            aux(i,j,17) = edgelengths(i,j,1)/dy      
            aux(i,j,18) = edgelengths(i,j,2)/dx

            !! Positive value : eta = h + height_temp
            !!aux(i,j,19) = height_temp
        enddo
    enddo

    if (integrate_bathy .eq. 0) then
        call sphere_set_bathy(blockno,mx,my,mbc, & 
                              xlower,ylower, dx,dy,aux,maux) 
    else
        call sphere_set_bathy_new(blockno,mx,my,mbc, & 
                                  xlower,ylower, dx,dy,aux,maux) 
    endif
end


subroutine sphere_set_bathy(blockno, mx,my,mbc,xlower,ylower, & 
    dx,dy,aux,maux)

    use lambwave_module, only : pi, earth_radius, mbathy
      
    implicit none

    integer mbc, mx,my, maux, blockno
    double precision xlower, ylower, dx,dy
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    integer*8 cont, fclaw_map_get_context

    logical is_corner


    double precision xc,yc, xp,yp,zp, r, lat, long

    integer i,j

    double precision bathymetry_interp

    cont = fclaw_map_get_context()

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc

           if (is_corner(i,j,mx,my)) then
              aux(i,j,mbathy) = 10000
!                cycle
               
            endif

            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            cont = fclaw_map_get_context()
            call fclaw_map_2d_c2m(cont,blockno,xc,yc,xp,yp,zp)

            r = sqrt(xp**2 + yp**2 + zp**2)

            if (abs(r-earth_radius)/earth_radius > 1d-8) then
                write(6,*) 'setaux : Something went wrong with radius calculation'
                write(6,*) 'radius = ', r
                write(6,*) 'earth_radius = ', earth_radius
                write(6,*) 'diff =', abs(r-earth_radius)/earth_radius
                write(6,*) xp,yp,zp
                write(6,*) 'xc,yc : ',xc,yc
                write(6,*) 'blockno ', blockno
                write(6,*) i,j
                stop
            endif 


            if (r .ne. r) then
                write(6,*) 'setaux : r is NAN'
                stop
            endif

            if (abs(zp/r) > 1.0) then
                write(6,*) 'setaux : zp/r > 1'
                write(6,'(A,F24.16)') 'zp/r = ', zp/r
                stop
            endif
            lat = -asin(zp/r)  !! value in [-pi/2,pi/2]
            long = atan2(yp, xp)  !! [-pi,pi]


            if (long < 0) then
                long = long + (2*pi)
            endif

            lat = lat*180.d0/pi
            long = long*180.d0/pi
!            write(6,*) ' '
!            write(6,'(A,2F24.16)') "lat/long ", lat, long
!            write(6,*) ' '

            
            !! Evaluate bathymetry at cell centers
            aux(i,j,mbathy) = bathymetry_interp(lat,long)
!!          write(6,*) ' '
!!          write(6,*) i,j
!!          write(6,'(A,2F24.16)') "lat/long ", lat, long
!!          write(6, '(A,2F24.16)') "aux: ", aux(i,j,mbathy)
!!          write(6,*) 

        enddo
    enddo 


    return 

end

!! This uses the exact area 
subroutine sphere_set_bathy_new(blockno, mx,my,mbc,xlower,ylower, & 
    dx,dy,aux,maux)
    use lambwave_module, only : mbathy, mcapa
      
    implicit none

    integer mbc, mx,my, maux, blockno
    double precision xlower, ylower, dx,dy
    double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    double precision area

    !!integer*8 cont, fclaw_map_get_context

    logical is_corner

    double precision xll, yll, int_BJ

    integer i,j

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc

            !! not sure why we are doing this
            if (is_corner(i,j,mx,my)) then
                aux(i,j,mbathy) = 10000
            endif


            !! Four corners of the mesh cell (i,j) in an aux array
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

            !! Use Gauss-Lobatto quadrature to integrate
            !! B(x,y)*J(x,y) over [xll,xur] x [yll, yur]

            !! This integrates the piecewise bilinear bathymetry surface
            !! over the mesh cell (i,j)

            !!int_BJ = bathy_integral(xll,yll,dx,dy)

            int_BJ = dx*dy  !! TODO : Include correct integral value

            !! Multiply capacity by dx*dy to get area
            area = aux(i,j,mcapa)*dx*dy

            !! Average value of bathymetry
            aux(i,j,mbathy) = int_BJ/area

        enddo
    enddo 

end






