c
c
c     =====================================================
      subroutine clawpack_flux3fw(ixy,maxm,meqn,maux,mbc,mx,
     &      q1d,dtdx1d,aux1,aux2,aux3,
     &      faddm,faddp,gaddm,gaddp,cfl1d,fwave,s,
     &      amdq,apdq,cqxx,bmasdq,bpasdq,rpn2,rpt2,
     &      mwaves,mcapa,method,mthlim)
c     =====================================================
c
c     # clawpack routine ...  modified for AMRCLAW
c
c--------------------------------------------------------------------
c     # flux2fw is a modified version of flux2 to use fwave instead of wave.
c     # A modified Riemann solver rp2n must be used in conjunction with this
c     # routine, which returns fwave's instead of wave's.
c     # See http://amath.washington.edu/~claw/fwave.html
c
c     # Limiters are applied to the fwave's, and the only significant
c     # modification of this code is in the "do 119" loop, for the
c     # second order corrections.
c
c--------------------------------------------------------------------
c
c
c     # Compute the modification to fluxes f and g that are generated by
c     # all interfaces along a 1D slice of the 2D grid.
c     #    ixy = 1  if it is a slice in x
c     #          2  if it is a slice in y
c     # This value is passed into the Riemann solvers. The flux modifications
c     # go into the arrays fadd and gadd.  The notation is written assuming
c     # we are solving along a 1D slice in the x-direction.
c
c     # fadd(i,.) modifies F to the left of cell i
c     # gadd(i,.,1) modifies G below cell i
c     # gadd(i,.,2) modifies G above cell i
c
c     # The method used is specified by method(2:3):
c
c         method(2) = 1 if only first order increment waves are to be used.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.
c
c         method(3) = 0 if no transverse propagation is to be applied.
c                       Increment and perhaps correction waves are propagated
c                       normal to the interface.
c                   = 1 if transverse propagation of increment waves
c                       (but not correction waves, if any) is to be applied.
c                   = 2 if transverse propagation of correction waves is also
c                       to be included.
c
c     Note that if mcapa>0 then the capa array comes into the second
c     order correction terms, and is already included in dtdx1d:
c     If ixy = 1 then
c        dtdx1d(i) = dt/dx                      if mcapa= 0
c                  = dt/(dx*aux(i,jcom,mcapa))  if mcapa = 1
c     If ixy = 2 then
c        dtdx1d(j) = dt/dy                      if mcapa = 0
c                  = dt/(dy*aux(icom,j,mcapa))  if mcapa = 1
c
c     Notation:
c        The jump in q (q1d(i,:)-q1d(i-1,:))  is split by rpn2 into
c            amdq =  the left-going flux difference  A^- Delta q
c            apdq = the right-going flux difference  A^+ Delta q
c        Each of these is split by rpt2 into
c            bmasdq = the down-going transverse flux difference B^- A^* Delta q
c            bpasdq =   the up-going transverse flux difference B^+ A^* Delta q
c        where A^* represents either A^- or A^+.
c
c
      implicit none

      external rpn2, rpt2

      integer ixy,mbc,maxm,meqn,maux,mx
      integer mwaves, mcapa, method(7), mthlim(mwaves)

      double precision cfl1d
      double precision    q1d(1-mbc:maxm+mbc, meqn)
      double precision   amdq(1-mbc:maxm+mbc, meqn)
      double precision   apdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   cqxx(1-mbc:maxm+mbc, meqn)
      double precision   faddm(1-mbc:maxm+mbc, meqn)
      double precision   faddp(1-mbc:maxm+mbc, meqn)
      double precision   gaddm(1-mbc:maxm+mbc, meqn, 2)
      double precision   gaddp(1-mbc:maxm+mbc, meqn, 2)
      double precision dtdx1d(1-mbc:maxm+mbc)
      double precision aux1(1-mbc:maxm+mbc, maux)
      double precision aux2(1-mbc:maxm+mbc, maux)
      double precision aux3(1-mbc:maxm+mbc, maux)
c
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision fwave(1-mbc:maxm+mbc, meqn, mwaves)
c
      logical limit
      double precision dtcom,dxcom,dycom,tcom
      integer                                icom,jcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
      double precision gupdate
      integer mw,i,jside,m

      limit = .false.
      do 5 mw=1,mwaves
         if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
c     # initialize flux increments:
c     -----------------------------
c
      do 30 jside=1,2
         do 20 m=1,meqn
            do 10 i = 1-mbc, mx+mbc
               faddm(i,m) = 0.d0
               faddp(i,m) = 0.d0
               gaddm(i,m,jside) = 0.d0
               gaddp(i,m,jside) = 0.d0
   10          continue
   20       continue
   30    continue
c
c
c     # solve Riemann problem at each interface and compute Godunov updates
c     ---------------------------------------------------------------------
c
      call rpn2(ixy,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &          aux2,aux2,fwave,s,amdq,apdq)
c
c     # Set fadd for the donor-cell upwind method (Godunov)
      do 40 i=1,mx+1
         do 40 m=1,meqn
            faddp(i,m) = faddp(i,m) - apdq(i,m)
            faddm(i,m) = faddm(i,m) + amdq(i,m)
   40       continue
c
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do 50 mw=1,mwaves
         do 50 i=1,mx+1
c          # if s>0 use dtdx1d(i) to compute CFL,
c          # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(i,mw),
     &                          -dtdx1d(i-1)*s(i,mw))
   50       continue
c
      if (method(2).eq.1) go to 130
c
c     # modify F fluxes for second order q_{xx} correction terms:
c     -----------------------------------------------------------
c
c     # apply limiter to fwaves:
      if (limit) call clawpack_inlinelimiter(maxm,meqn,mwaves,mbc,mx,
     &      fwave,s,mthlim)
c
      do 120 i = 1, mx+1
c
c        # For correction terms below, need average of dtdx in cell
c        # i-1 and i.  Compute these and overwrite dtdx1d:
c
         dtdx1d(i-1) = 0.5d0 * (dtdx1d(i-1) + dtdx1d(i))
c
         do 120 m=1,meqn
            cqxx(i,m) = 0.d0
            do 119 mw=1,mwaves
c
c              # second order corrections:
               cqxx(i,m) = cqxx(i,m) + dsign(1.d0,s(i,mw))
     &            * (1.d0 - dabs(s(i,mw))*dtdx1d(i-1)) * fwave(i,m,mw)
c
  119          continue
            faddm(i,m) = faddm(i,m) + 0.5d0 * cqxx(i,m)
            faddp(i,m) = faddp(i,m) + 0.5d0 * cqxx(i,m)
  120       continue
c
c
  130  continue
c
       if (method(3).eq.0) go to 999   !# no transverse propagation
c
       if (method(3).eq.2) then
c         # incorporate cqxx into amdq and apdq so that it is split also.
          do 150 i = 1, mx+1
             do 150 m=1,meqn
                amdq(i,m) = amdq(i,m) + cqxx(i,m)
                apdq(i,m) = apdq(i,m) - cqxx(i,m)
  150           continue
          endif
c
c
c      # modify G fluxes for transverse propagation
c      --------------------------------------------
c
c
c     # split the left-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &          q1d,q1d,aux1,aux2,aux3,
     &          1,amdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
      do 160 m=1,meqn
          do 160 i = 1, mx+1
               gupdate = 0.5d0*dtdx1d(i-1) * bmasdq(i,m)
               gaddm(i-1,m,1) = gaddm(i-1,m,1) - gupdate
               gaddp(i-1,m,1) = gaddp(i-1,m,1) - gupdate
c
               gupdate = 0.5d0*dtdx1d(i-1) * bpasdq(i,m)
               gaddm(i-1,m,2) = gaddm(i-1,m,2) - gupdate
               gaddp(i-1,m,2) = gaddp(i-1,m,2) - gupdate
  160          continue
c
c     # split the right-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &          q1d,q1d,aux1,aux2,aux3,
     &          2,apdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
      do 180 m=1,meqn
          do 180 i = 1, mx+1
               gupdate = 0.5d0*dtdx1d(i-1) * bmasdq(i,m)
               gaddm(i,m,1) = gaddm(i,m,1) - gupdate
               gaddp(i,m,1) = gaddp(i,m,1) - gupdate
c
               gupdate = 0.5d0*dtdx1d(i-1) * bpasdq(i,m)
               gaddm(i,m,2) = gaddm(i,m,2) - gupdate
               gaddp(i,m,2) = gaddp(i,m,2) - gupdate
  180          continue
c
  999 continue
      return
      end
