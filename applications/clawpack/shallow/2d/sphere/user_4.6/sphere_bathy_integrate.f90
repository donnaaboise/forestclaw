#if 0
!! Fix syntax highlighting
#endif

module gauss_labatto
    implicit none

    double precision, dimension(4) :: X4 = &
        (/-1.d0, & 
         -0.447213595499957939281834733746d0, & 
          0.447213595499957939281834733746d0, & 
          1.d0/)

    double precision, dimension(4) :: W4 = & 
        (/0.166666666666666666666666666667d0, & 
        0.833333333333333333333333333333d0, & 
        0.833333333333333333333333333333d0, & 
        0.166666666666666666666666666667d0/)

    double precision, dimension(6) :: X6 = &
        (/-1.d0, & 
        -0.765055323929464692851002973959d0, & 
        -0.285231516480645096314150994041d0, & 
        0.285231516480645096314150994041d0, & 
        0.765055323929464692851002973959d0, & 
        1.d0/)

    double precision, dimension(6) :: W6 = & 
        (/0.0666666666666666666666666666667d0, & 
        0.378474956297846980316612808212d0, & 
        0.554858377035486353016720525121d0, & 
        0.554858377035486353016720525121d0, & 
        0.378474956297846980316612808212d0, & 
        0.0666666666666666666666666666667d0/)

        !! Add subroutine that computes the points above. 

end


!! This creates a new bathymetry data for each of the six panels of the
!! cubed sphere.   The size is hardwired in the lambwave module. 
subroutine sphere_bathy_project()

    use lambwave_module, only : pi, earth_radius, bathy_project, & 
        mx_bathy, my_bathy, dx_bathy, dy_bathy
      
    implicit none

    !! Local variables
    integer*8 cont, fclaw_map_get_context

    !!logical is_corner

    double precision xc,yc, xp,yp,zp, r, lat, long

    integer blockno
    integer i,j

    double precision bathymetry_interp

    cont = fclaw_map_get_context()

    !! Loop over nodes of each panel of the cubed sphere;  interpolate
    !! from original bathymetry to this new six-panel bathymetry. 
    do blockno = 0,5
        do i = 0,mx_bathy+1
            do j = 0,my_bathy +1

                !! computational coordinates are [0,1]x[0,1] for 
                !! each panel
                xc = i*dx_bathy
                yc = j*dy_bathy

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
                    write(6,*) 'bathy_project : r is NAN'
                    stop
                endif

                if (abs(zp/r) > 1.0) then
                    write(6,*) 'bathy_project : zp/r > 1'
                    write(6,'(A,F24.16)') 'zp/r = ', zp/r
                    stop
                endif
                lat = -asin(zp/r)     !! value in [-pi/2,pi/2]
                long = atan2(yp, xp)  !! [-pi,pi]


                if (long < 0) then
                    long = long + (2*pi)
                endif

                lat = lat*180.d0/pi
                long = long*180.d0/pi

                bathy_project(i,j,blockno) = bathymetry_interp(lat,long)                
            enddo
        enddo 
    enddo 

    return 

end


double precision function sphere_pw_bilinear_bathy(blockno, xc,yc)

    use lambwave_module, only : bathy_project, & 
        mx_bathy, my_bathy, dx_bathy, dy_bathy

    implicit none

    integer blockno
    double precision xc,yc

    integer i1, j1
    double precision b00, b01, b10, b11
    double precision a00, a01, a10, a11
    double precision :: b

    i1 = int(floor(xc/dx_bathy))
    j1 = int(floor(yc/dy_bathy))

    if (i1 .lt. 0 .or. i1 .gt. mx_bathy) then
        write(6,*) 'sphere_bathy_project_eval : '
        write(6,*) 'i1 < 0 or i1 > 1 : ', i1
        stop
    endif

    if (j1 .lt. 0 .or. j1 .gt. my_bathy) then
        write(6,*) 'sphere_bathy_project_eval : '
        write(6,*) 'j1 < 0 or j1 > 1 : ', j1
        stop
    endif

!!            b00 = bmount(xc[i1],yc[j1],**kwargs)
!!            b01 = bmount(xc[i1+1],yc[j1],**kwargs)
!!            b10 = bmount(xc[i1],yc[j1+1],**kwargs)
!!            b11 = bmount(xc[i1+1],yc[j1+1],**kwargs)
!!            a00 = b00
!!            a01 = b01 - b00
!!            a10 = b10 - b00
!!            a11 = b11 - b01 - b10 + b00
!!            return a00 + a10*xi + a01*eta + a11*xi*eta

    !! get bathymetry values at four corners of mesh cell containing
    !! (xc,yc)
    b00 = bathy_project(i1,j1,blockno)
    b01 = bathy_project(i1+1,j1,blockno)
    b10 = bathy_project(i1,j1+1,blockno)
    b11 = bathy_project(i1+1,j1+1,blockno)

    !! Construct bilinear surface of bathymetry value
    a00 = b00
    a01 = b01 - b00
    a10 = b10 - b00
    a11 = b11 - b01 - b10 + b00
    b =  a00 + a10*xc + a01*yc + a11*xc*yc

    sphere_pw_bilinear_bathy = b

    return
end




logical function is_corner(i,j,mx,my)

    implicit none

    integer i,j,mx,my


    is_corner = .false.
    if (i .lt. 1 .or. i .gt. mx) then
        if (j .lt. 1 .or. j .gt. my) then 
            is_corner = .true.
            return
        endif
    else if (j .lt. 1 .or.  j .gt. my) then
        if (i .lt. 1 .or. i .gt. mx) then 
            is_corner = .true.
            return
        endif
    endif

    return

end function is_corner



double precision function bathy_integral(xi,yj,dx,dy)

!!###########################################################################
!! This code was taken from the Geoclaw library and adopted for use in 
!! ForestClaw
!!
!!     topointegral integrates a surface over a rectangular region
!!     that is the intersection with a Cartesion grid (of topography data)
!!     the surface integrated is defined by a piecewise bilinear through the
!!     nodes of the Cartesian grid.
!!
!!      The rectangular intersection (a computational grid cell) has coords:
!!
!!              xim <= x <= xip, yjm <= y <= yjp
!!
!!      The topography data grid has coords:
!!
!!           xlow_bathy <= x <= xhi_mbathy, ylow_bathy <= y <= yhi_mbathy, 
!!
!!      with grid cell size dx_mbathy by dy_mbathy and 
!!      mx_bathy by my_bathy cells.
!!
!!      For the cubedsphere mesh, we have xlow_mbathy = ylow_mbathy = 0
!!      and xhi_mbathy = yhi_mbathy = 1.  These are the computational
!!      coordinates for a single panel of the cubed sphere.
!!
!!                                                written by David L. George
!!                                                Seattle, WA 7/16/08
!!###########################################################################

    use lambwave_module, only : mx_bathy, my_bathy, dx_bathy, dy_bathy

    implicit none

    double precision xi,yj,dx,dy

    double precision x1, x2, y1, y2, xip, yjp, bint
    integer istart, iend, jstart, jend
    integer i,j


!! =============INTEGRATE PIECEWISE BILINEAR OVER RECTANGULAR REGION====
!!  # don't waste time looping through the entire grid
!!         just find indices that include the computational mesh cell

    xip = xi + dx
    yjp = yj + dy

    istart = int(floor(xi/dx_bathy))
    jstart = int(floor(yj/dy_bathy))

    iend = int(floor(xip/dx_bathy))
    jend = int(floor(yjp/dy_bathy))

    istart = max(istart,0)
    jstart = max(jstart,0)
    iend = min(iend,mx_bathy)
    jend = min(jend,my_bathy)


    !! # initialize:
    bint = 0.d0
    do j = jstart,jend
        y1 = j*dy_bathy
        y2 = y1 + dy_bathy

        do i = istart,iend
            x1 = i*dx_bathy
            x2 = x1 + dx_bathy

            !! Integrate over the intersection of the mesh cell with 
            !! the bathymetry cell (ii,jj)
            !! bint = bint + gll_integrate( )
            bint = dx*dy  !! TODO : Correct using gll integration
        enddo
    enddo


    bathy_integral = bint
    return
end function