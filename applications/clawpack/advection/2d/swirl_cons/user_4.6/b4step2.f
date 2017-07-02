      subroutine clawpack46_b4step2(maxmx, maxmy, mbc,mx,my,meqn,q,
     &      xlower,ylower,dx,dy,time,dt,maux,aux)
      implicit none

      integer mbc, mx, my, meqn, maux, maxmx, maxmy
      double precision xlower, ylower, dx, dy, time, dt
      double precision q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer i, j
      double precision tperiod, pi2, vt, xll, yll, psi, uec, vec

      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in
c        # setaux should be used for all time.
         return
      endif

C       vt = cos(pi2*(time+dt/2.d0)/tperiod)
      vt = 1.0

      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
C             aux(i,j,1) = (psi(xll, yll+dy) - psi(xll,yll)) / dy
C             aux(i,j,2) =  -(psi(xll+dx, yll) - psi(xll,yll)) / dx
            
c           # edge velocity
            aux(i,j,1) = uec(xll,yll)
            aux(i,j,2) = vec(xll,yll)

c           # multiply by time-factor:
            aux(i,j,1) = vt * aux(i,j,1)
            aux(i,j,2) = vt * aux(i,j,2)
         enddo
      enddo

      return
      end

      double precision function uec(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

C       u = -2*((sin(pi*xp))**2 * sin(pi*yp) * cos(pi*yp)) + 2
      if(xp .gt. 0) then
            uec = exp(-(xp-0.5)**2/(2*0.1**2))
      else
            uec = 0
      endif

C       if((xp .ge. 0.5 .and. xp .lt. 1) .or. (xp .lt. 0)) then
C             u = 1.0
C       else
C             u = 0.5
C       endif

      return
      end

      double precision function vec(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

C       v = 2*((sin(pi*yp))**2 * sin(pi*xp) * cos(pi*xp))
      vec = 0

      return
      end
