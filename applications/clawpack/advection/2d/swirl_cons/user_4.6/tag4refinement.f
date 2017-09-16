      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax, xc, yc

      tag_patch = 0


c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
c            qmin = min(q(i,j,mq),qmin)
c            qmax = max(q(i,j,mq),qmax)
            if (q(i,j,mq) .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
c            xc = xlower + (i-0.5)*dx
c            yc = ylower + (j-0.5)*dy
c            if (abs(yc-0.5d0) .lt. dy .and. abs(xc-0.5d0) .lt. dx) then
c               tag_patch = 1
c               return
c            endif

         enddo
      enddo

      end