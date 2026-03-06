c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================

      Subroutine MeanField(PtINI,PtEND,X,UMF)

      Implicit None

      Integer i,i0,i1,i2,i3
      Integer npt
      Parameter (npt = 4)
      Real*8 ptINI(2,npt) ! IN -> initial coord of corners
      Real*8 ptEND(2,npt) ! IN -> final coord of corners
      Real*8 X(2)         ! IN ->  initial coord of the point
      Real*8 UMF(2)       ! OUT -> computed meanfield displacement
      Real*8 xij,yij
      Real*8 x0,y0,x1,y1,x2,y2,x3,y3,xi,yi
      Real*8 s,t,s2,t2
      Real*8 xmean,ymean,xx0,xx1,xx2,xx3,yy0,yy1,yy2,yy3,xxi,yyi

c-----The 4 corners should be numbered like that
c      p2 --- p3
c      |      |
c    t |  p   |
c      |      |
c      p0 --- p1
c         s

c     -- mean value for X and Y coordinates
      xmean = 0.d0
      ymean = 0.d0
      Do i = 1,npt
         xmean = xmean + ptINI(1,i)
         ymean = ymean + ptINI(2,i)
      EndDo
      xmean = xmean / dfloat(npt)
      ymean = ymean / dfloat(npt)

c     -- the vector from the center to p0 gives (cos,sin) < 0
      Do i = 1,npt
         xij = ptINI(1,i) - xmean
         yij = ptINI(2,i) - ymean
         If ( (xij.LT.0.d0).AND.(yij.LT.0.d0) ) i0 = i
      EndDo
c     -- the vector from the center to p1 gives cos >0, sin <0
      Do i = 1,npt
         xij = ptINI(1,i) - xmean
         yij = ptINI(2,i) - ymean
         If ( (xij.GT.0.d0).AND.(yij.LT.0.d0) ) i1 = i
      EndDo
c     -- the vector from the center to p2 gives cos  < 0, sin  > 0
      Do i = 1,npt
         xij = ptINI(1,i) - xmean
         yij = ptINI(2,i) - ymean
         If ( (xij.LT.0.d0).AND.(yij.GT.0.d0) ) i2 = i
      EndDo
c     -- the vector from the center to p3 gives cos  > 0, sin  > 0
      Do i = 1,npt
         xij = ptINI(1,i) - xmean
         yij = ptINI(2,i) - ymean
         If ( (xij.GT.0.d0).AND.(yij.GT.0.d0) ) i3 = i
      EndDo

c-----Compute XMF(1) and XMF(2), the "meanfield" disp of X.
      x0 = ptINI(1,i0)
      y0 = ptINI(2,i0)
      x1 = ptINI(1,i1)
      y1 = ptINI(2,i1)
      x2 = ptINI(1,i2)
      y2 = ptINI(2,i2)
      x3 = ptINI(1,i3)
      y3 = ptINI(2,i3)
      xi = X(1)
      yi = X(2)
      Call inverseBilinear(x0,y0,x1,y1,x2,y2,x3,y3,xi,yi,s,t,s2,t2)

      xx0 = ptEND(1,i0) - ptINI(1,i0)
      yy0 = ptEND(2,i0) - ptINI(2,i0)
      xx1 = ptEND(1,i1) - ptINI(1,i1)
      yy1 = ptEND(2,i1) - ptINI(2,i1)
      xx2 = ptEND(1,i2) - ptINI(1,i2)
      yy2 = ptEND(2,i2) - ptINI(2,i2)
      xx3 = ptEND(1,i3) - ptINI(1,i3)
      yy3 = ptEND(2,i3) - ptINI(2,i3)
      Call bilinear(xx0,yy0,xx1,yy1,xx2,yy2,xx3,yy3,s,t,xxi,yyi)
      UMF(1) = xxi
      UMF(2) = yyi

      Return
      End
c=======================================================================

      Subroutine inverseBilinear(x0,y0,x1,y1,x2,y2,x3,y3,x,y,
     &                           sout,tout,s2out,t2out)

      Implicit Real*8 (a-h,o-z)
      Integer t_valid, t2_valid
      Real*8 a,b1,b2,c,b
      Real*8 s,s2,t,t2
      Real*8 am2bpc,sqrtbsqmac
      Real*8 tdenom_x,tdenom_y
      Integer num_valid_s
      Logical equals,in_range

      a  = cross2( x0-x, y0-y, x0-x2, y0-y2 )
      b1 = cross2( x0-x, y0-y, x1-x3, y1-y3 )
      b2 = cross2( x1-x, y1-y, x0-x2, y0-y2 )
      c  = cross2( x1-x, y1-y, x1-x3, y1-y3 )
      b  = 0.5d0 * (b1 + b2)
      am2bpc = a - 2.d0*b+c
      num_valid_s = 0

      If (equals(am2bpc,0.d0,1.d-10)) Then
         If (equals(a-c,0.d0,1.d-10)) Then
            Write(6,*) 'Looks like the input is a line'
            Write(6,*) 'You could set s=0.5 and solve for t'//
     &                 ' if you wanted to'
            Call Exit(2)
         EndIf
         s = a / (a-c)
         If (in_range(s,0.d0,1.d0,1.d-10)) num_valid_s = 1
      Else
         sqrtbsqmac = sqrt(b*b - a*c)
         s  = ((a-b) - sqrtbsqmac) / am2bpc
         s2 = ((a-b) + sqrtbsqmac) / am2bpc
         num_valid_s = 0
         If (in_range(s,0.d0,1.d0,1.d-10)) Then
            num_valid_s = num_valid_s + 1
            If (in_range(s2,0.d0,1.d0,1.d-10)) Then
               num_valid_s = num_valid_s + 1
            EndIf
         Else
            If (in_range(s2,0.d0,1.d0,1.d-10)) Then
               num_valid_s = num_valid_s + 1
               s = s2
            EndIf
         EndIf
      EndIf

      If (num_valid_s.EQ.0) Call Exit(2)
    
      t_valid = 0

      If (num_valid_s.GE.1) Then
         tdenom_x = (1.d0-s)*(x0-x2) + s*(x1-x3)
         tdenom_y = (1.d0-s)*(y0-y2) + s*(y1-y3)
         t_valid = 1
         If (equals(tdenom_x,0.d0,1.d-10).AND.
     &       equals(tdenom_y,0.d0,1.d-10)
     &      ) Then
            t_valid = 0
         Else 
c           -- Choose the more robust denominator
            If (abs(tdenom_x).GT.abs(tdenom_y)) Then
               t = ( (1.d0-s)*(x0-x) + s*(x1-x) ) / tdenom_x
            Else
               t = ( (1.d0-s)*(y0-y) + s*(y1-y) ) / tdenom_y;
            EndIf
            If (.NOT.in_range(t,0.d0,1.d0,1.d-10 )) t_valid = 0
         EndIf
      EndIf
c     --Same thing for s2 and t2
      If (num_valid_s.EQ.2) then
         tdenom_x = (1.d0-s2)*(x0-x2) + s2*(x1-x3)
         tdenom_y = (1.d0-s2)*(y0-y2) + s2*(y1-y3)
         t2_valid = 1
         If (equals(tdenom_x,0.d0,1.d-10).AND.
     &       equals(tdenom_y,0.d0,1.d-10)
     &      ) Then
            t2_valid = 0
         Else
c          -- Choose the more robust denominator
            If (abs(tdenom_x).GT.abs(tdenom_y)) Then
               t2 = ( (1.d0-s2)*(x0-x) + s2*(x1-x) ) / tdenom_x 
            Else
               t2 = ( (1.d0-s2)*(y0-y) + s2*(y1-y) ) / tdenom_y
            EndIf
            If (.NOT.in_range(t2,0.d0,1.d0,1.d-10)) t2_valid = 0
	 EndIf
      EndIf

c     -- Final cleanup
      If ((t2_valid.EQ.1).AND.(t_valid.NE.1)) Then
         s = s2
         t = t2
         t_valid = t2_valid
         t2_valid = 0
      EndIf
c     -- Output
      If (t_valid.EQ.1) Then
         sout = s
         tout = t
      EndIf

      If (t2_valid.EQ.1) Then
         s2out = s2
         t2out = t2
      EndIf

c	return t_valid + t2_valid;

      Return
      End
c=======================================================================

      Subroutine bilinear(x0,y0,x1,y1,x2,y2,x3,y3,s,t,x,y)

      Implicit Real*8 (a-h,o-z)

      x = t*(s*x3+(1.d0-s)*x2) + (1.d0-t)*(s*x1+(1.d0-s)*x0)
      y = t*(s*y3+(1.d0-s)*y2) + (1.d0-t)*(s*y1+(1.d0-s)*y0)

      Return
      End 
c=======================================================================
      Real*8 Function cross2(x0,y0,x1,y1)

      Implicit Real*8 (a-h,o-z)
      Real*8 x0,y0,x1,y1

      cross2 = x0*y1 - y0*x1
    
      Return
      End
c=======================================================================

      Logical Function equals(a,b,tolerance)

      Implicit Real*8 (a-h,o-z)
      Real*8 a,b,tolerance

      bmin = b - tolerance
      bmax = b + tolerance
      If ((a.LE.bmax).AND.(a.GE.bmin)) Then
         equals = .TRUE.
      Else
         equals = .FALSE.
      EndIf

      Return
      End
c=======================================================================
  
      Logical Function in_range(val,range_min,range_max,tol)

      Implicit Real*8 (a-h,o-z)
      Real*8 val,range_min,range_max,tol

      If ( ((val+tol).GE.range_min).AND.
     &     ((val-tol).LE.range_max) 
     &   ) Then
         in_range = .TRUE.
      Else
         in_range = .FALSE.
      EndIf

      Return
      End
c=======================================================================


