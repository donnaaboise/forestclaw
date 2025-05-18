subroutine fclaw_map_2d_c2m_cubedsphere(blockno,xc,yc,xp,yp,zp)
   implicit none

   double precision xc,yc,xp,yp,zp
   double precision xc1,yc1
   integer blockno,blockno1

   call csphere_handle_ghost(blockno,xc,yc,xc1,yc1,blockno1)

   if (blockno1 .eq. 0) then
      call csphere_basic(xc1,yc1,yp,xp,zp)
      zp = -zp
   elseif (blockno1 .eq. 1) then
      call csphere_basic(xc1,yc1,zp,xp,yp)
   elseif (blockno1 .eq. 2) then
      call csphere_basic(xc1,yc1,zp,yp,xp)
      xp = -xp
   elseif (blockno1 .eq. 3) then
      call csphere_basic(xc1,yc1,xp,yp,zp)
   elseif (blockno1 .eq. 4) then
      call csphere_basic(xc1,yc1,xp,zp,yp)
      yp = -yp
   elseif (blockno1 .eq. 5) then
      call csphere_basic(xc1,yc1,yp,zp,xp)
   endif

end subroutine fclaw_map_2d_c2m_cubedsphere


subroutine csphere_basic(xc,yc,xp,yp,zp)
   implicit none

   double precision xc,yc,xp,yp,zp
   double precision R, tan_xi, tan_eta

      
   double precision pi, pi2
   common /compi/ pi, pi2

   R = 1.d0
   tan_xi = tan(0.5d0*pi*(xc-0.5d0))
   tan_eta = tan(0.5d0*pi*(yc-0.5d0))
   zp = R/sqrt(tan_xi**2 + tan_eta**2 + 1.d0)
   xp = zp*tan_xi
   yp = zp*tan_eta

end subroutine csphere_basic


!!   T_left   = {edge.bottom : lambda x,y : (y,-x),  edge.right : lambda x,y : (x+1,y) }
!!   T_right  = {edge.left   : lambda x,y : (x-1,y), edge.top   : lambda x,y : (y,2-x) }
!!   T_bottom = {edge.left   : lambda x,y : (-y,x),  edge.top   : lambda x,y : (x,1+y) }
!!   T_top    = {edge.bottom : lambda x,y : (x,y-1), edge.right : lambda x,y : (2-y,x) }
!!    
!!   T = {edge.left : T_left, edge.right : T_right, edge.bottom : T_bottom, edge.top : T_top}


!! --------------- T_left -------------
subroutine T_left_bottom(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = y
   y1 = -x
   return
end subroutine T_left_bottom

subroutine T_left_right(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = x+1
   y1 = y
   return
end subroutine T_left_right

!! --------------- T_right -------------
subroutine T_right_left(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = x-1
   y1 = y
   return
end subroutine T_right_left

subroutine T_right_top(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = y
   y1 = 2-x
   return
end subroutine T_right_top

!! --------------- T_bottom -------------
subroutine T_bottom_left(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = -y
   y1 = x
   return
end subroutine T_bottom_left

subroutine T_bottom_top(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = x
   y1 = 1+y
   return
end subroutine T_bottom_top


!! --------------- T_top -------------
subroutine T_top_bottom(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = x
   y1 = y-1
   return
end subroutine T_top_bottom

subroutine T_top_right(x,y,x1,y1)
   implicit none
   double precision x,y,x1,y1
   x1 = 2-y
   y1 = x
   return
end subroutine T_top_right


subroutine csphere_handle_ghost(blockno,xc,yc,xc1,yc1,blockno1)
   implicit none

   double precision xc,yc,xc1,yc1
   integer blockno, blockno1

   logical x_in, y_in, interior, corner

   integer blocks(0:5,0:3,0:1), e0, e1

   !! (blockno, edge,0) = neighbor block (in [0,5])
   !! (blockno, edge,0) = neighbor block edge (in [0,3])

   !! Block 0
   blocks(0,0,0) = 4
   blocks(0,0,1) = 2
   blocks(0,1,0) = 1
   blocks(0,1,1) = 0
   blocks(0,2,0) = 2
   blocks(0,2,1) = 0
   blocks(0,3,0) = 5
   blocks(0,3,1) = 2

   !! Block 1
   blocks(1,0,0) = 0
   blocks(1,0,1) = 1
   blocks(1,1,0) = 3
   blocks(1,1,1) = 3
   blocks(1,2,0) = 2
   blocks(1,2,1) = 3
   blocks(1,3,0) = 5
   blocks(1,3,1) = 1

   !! Block 2
   blocks(2,0,0) = 0
   blocks(2,0,1) = 2
   blocks(2,1,0) = 3
   blocks(2,1,1) = 0
   blocks(2,2,0) = 4
   blocks(2,2,1) = 0
   blocks(2,3,0) = 1
   blocks(2,3,1) = 2

   !! Block 3
   blocks(3,0,0) = 2
   blocks(3,0,1) = 1
   blocks(3,1,0) = 5
   blocks(3,1,1) = 3
   blocks(3,2,0) = 4
   blocks(3,2,1) = 3
   blocks(3,3,0) = 1
   blocks(3,3,1) = 1

   !! Block 4
   blocks(4,0,0) = 2
   blocks(4,0,1) = 2
   blocks(4,1,0) = 5
   blocks(4,1,1) = 0
   blocks(4,2,0) = 0
   blocks(4,2,1) = 0
   blocks(4,3,0) = 3
   blocks(4,3,1) = 2

   !! Block 5
   blocks(5,0,0) = 4
   blocks(5,0,1) = 1
   blocks(5,1,0) = 1
   blocks(5,1,1) = 3
   blocks(5,2,0) = 0
   blocks(5,2,1) = 3
   blocks(5,3,0) = 3
   blocks(5,3,1) = 1

   x_in = (0 .le. xc) .and. (xc .le. 1)
   y_in = (0 .le. yc) .and. (yc .le. 1)

   interior = x_in .and. y_in

   corner = (.not. x_in) .and. (.not. y_in)

   if ((.not. interior) .and. (.not. corner)) then
        
      if (xc .lt. 0) then
         e0 = 0
      else if (xc .gt. 1) then
         e0 = 1
      elseif (yc .lt. 0) then
         e0 = 2
      else
         e0 = 3
      endif

      !! blockno1 = block[blockno][e0].blockno
      blockno1 = blocks(blockno,e0,0)
            
      !! e1 = block[blockno][e0].edge
      e1 = blocks(blockno,e0,1)   !! block[blockno][e0].edge
        
      !!xc1,yc1 = T[e0][e1](xc,yc)
      if (e0 .eq. 0) then
         if (e1 .eq. 2) then
            call T_left_bottom(xc,yc,xc1,yc1)
         else
            call T_left_right(xc,yc,xc1,yc1)
         endif
      else if (e0 .eq. 1) then
         if (e1 .eq. 0) then
            call T_right_left(xc,yc,xc1,yc1)
         else
            call T_right_top(xc,yc,xc1,yc1)
         endif
      else if (e0 .eq. 2) then
         if (e1 .eq. 0) then
            call T_bottom_left(xc,yc,xc1,yc1)
         else
            call T_bottom_top(xc,yc,xc1,yc1)
         endif
      else
         if (e1 .eq. 2) then
            call T_top_bottom(xc,yc,xc1,yc1)
         else
            call T_top_right(xc,yc,xc1,yc1)
         endif
      endif

   else
      blockno1 = blockno
      xc1 = xc
      yc1 = yc
   endif

   return
end subroutine csphere_handle_ghost
    
