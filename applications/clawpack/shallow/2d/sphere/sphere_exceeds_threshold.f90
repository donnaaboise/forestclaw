!! # check to see if value exceeds threshold

integer function user_exceeds_threshold(blockno,meqn,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,ivar_threshold, threshold, &
                                     init_flag, is_ghost)
    implicit none
    
    integer :: meqn,ivar_threshold
    double precision :: qval(meqn),qmin(meqn),qmax(meqn),threshold
    double precision :: quad(-1:1,-1:1,meqn)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag, refine
    logical(kind=4) :: is_ghost

   double precision hin, hout
   common /swe_initcond_parms4/  hin,hout

    integer :: mq
    double precision :: pi


    refine = 0

!!    if (blockno .eq. 4) then
!!        refine = 1
!!    endif

    block
        integer*8 cont, fclaw_map_get_context
        double precision xp,yp,zp

        cont = fclaw_map_get_context()

        call fclaw_map_2d_c2m(cont,blockno,xc,yc,xp,yp,zp)

        if (xp .ge. 0.1) then
            refine = 1
        endif
    end block

    user_exceeds_threshold = refine

    return

    pi = 4.d0*atan(1.d0)

    mq = ivar_threshold

    if (.not. is_ghost) then
        if (abs(qval(1)-hout) .gt. threshold) then
            refine = 1
        endif
    endif

    user_exceeds_threshold = refine

end function user_exceeds_threshold
