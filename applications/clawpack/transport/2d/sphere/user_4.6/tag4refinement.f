      subroutine sphere_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,t, blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy, t
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer example
      common /example_comm/ example      

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer initchoice
      common /initchoice_comm/ initchoice

      integer*8 cont, get_context

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      logical refine

      cont = get_context()

c     # Refine based only on first variable in system.
      refine = .false.
      mq = 1
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            if (init choice .le. 3) then
                refine = q(i,j,mq) .gt.  tag_threshold              
            endif
            if (refine) then
                tag_patch = 1
                return
            endif
         enddo
      enddo


      end

