c2345678901234567890123456789012345678901234567890123456789012
c=============================================================
      Program KStest


      Parameter(Lmax=5000000)
      Common/tab/ Val(3,Lmax),nval
      Real*4 input(Lmax),fit(Lmax)
      Dimension p(4)
      Dimension xi(4,4)
      Character*50 Fichier

      Write(6,*) 'Nom de fichier ?'
      Read(5,*) Fichier
      
c-----Ouverture du fichier
      Open(1,file=Fichier,status='old',err=99)
c-----Lecture du fichier
      k = 0
1     Read(1,*,end=2) x1,x2
          k = k + 1
          val(1,k) = x1 
          val(2,k) = x2
          GoTo 1
2     Close(1)
      nval = k
      Write(6,*) 'Nbre de points :',nval

c-----echantillonnage
      xmin  = 1.e20
      xmax = -1.e20
      Do i = 1,nval
         xmin = min(xmin,val(1,i))
         xmax = max(xmax,val(1,i))
      EndDo
      Write(6,*) "xmin,xmax = ",xmin,xmax
      Write(6,*) "xrange ? (x1min,x1max)" 
      Read(5,*) x1min,x1max

      Open(1,file=Fichier,status='old',err=99)
      k = 0
10     Read(1,*,end=20) x1,x2,x3
        If ( (x1.GE.x1min).AND.(x1.LE.x1max)) Then
          k = k + 1
          val(1,k) = x1 
          val(2,k) = x2
          val(3,i) = x3
        EndIf
        GoTo 10
20     Close(1)
      nval = k
      Write(6,*) 'Nbre de points :',nval

c-----Valeurs initiales pour le powell
      n = 4
      np = 4
      p(1) = 0.0001
      p(2) = 0.14
      p(3) = 0.1
      p(4) = 1.64

      Do i = 1,n
         Do j = 1,n  
            xi(i,j) = 0.d0
         EndDo
      EndDo
      xi(1,1) = -0.01d0
      xi(2,2) = -0.01d0
      xi(3,3) = -0.01d0
      xi(4,4) = 0.01d0

      ftol = 1.d-15

      fret = func(p)
      Write(6,*) 'fret =  ',1-fret

      Call Powell(p,xi,n,np,ftol,iter,fret)
      Write(6,*) 'iter =  ',iter
      Write(6,*) 'fret =  ',1-fret

      Write(6,*) "x0 = ",p(1)
      Write(6,*) "a = ", p(2)
      Write(6,*) "b = ", p(3)
      Write(6,*) "q = ", p(4)
      Write(6,*) "f(x) = a*(1.-(1-q)*b*(x-x0)*(x-x0)"//
     &            ")**(1./(1.-q))" 
     
c-----Verification - coeff de correlation de pearson
      x0 = p(1)
      a = p(2)
      b = p(3)
      q = p(4)
      Do i = 1,nval
         input(i) = val(2,i)
         xxx = val(1,i)
         fit(i) = a*(1.-(1.-q)*b*(xxx-x0)*(xxx-x0))**(1./(1.-q))
         Write(64,*) input(i),fit(i)
      EndDo

      Call pearsn(input,fit,nval,r,prob,z)
      Write(6,*) "r = ",r
      Write(6,*) "prob = ", prob

     
      Stop

99    Write(6,*) 'File does not exist'
      Stop

      End

c=============================================================
      Real Function func(p)
      real p(*)
      Parameter(Lmax=5000000)
      Common/tab/ Val(3,Lmax),nval

      x0 = p(1)
      xa = p(2)
      xb = p(3)
      xq = p(4)

c      sum1 = 0.
c      Do i = 1,nval
c         xxx = Val(1,i)
c         yyy = Val(2,i)
c         qgauss = xa*(1.-(1.-xq)*xb*(xxx-x0)*(xxx-x0))**(1./(1.-xq))
c         sum1 = sum1 + (yyy - qgauss)**2. * qgauss
c      EndDo
c      func = sum1
c      GoTo 1


      xm = 0.
      ym = 0.
      Do i = 1,nval
         xxx = Val(1,i)
         yyy = Val(2,i)
         xm = xm + (yyy)
         qgauss = xa*(1.-(1.-xq)*xb*(xxx-x0)*(xxx-x0))**(1./(1.-xq))
         ym = ym + (qgauss)
      EndDo
      xm = xm / float(nval)
      ym = ym / float(nval)
      rnum = 0.
      rdem1 = 0.
      rdem2 = 0.
      Do i = 1,nval
         xxx = Val(1,i)
         yyy = (Val(2,i))
         yyy = (yyy)
         qgauss =  xa*(1.-(1.-xq)*xb*(xxx-x0)*(xxx-x0))**(1./(1.-xq))  
         qgauss = (qgauss)
         rnum = rnum + (yyy - xm)*(qgauss-ym)
         rdem1 = rdem1 + (yyy - xm)*(yyy - xm)
         rdem2 = rdem2 + (qgauss-ym)*(qgauss-ym)
      EndDo
      func = 1 - rnum / sqrt(rdem1*rdem2)




1      Return
      End
c=============================================================
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100000,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      REAL pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      REAL fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.e-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
      INTEGER iter,n,np,NMAX,ITMAX
      REAL fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=20,ITMAX=2000)
CU    USES func,linmin
      INTEGER i,ibig,j
      REAL del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2.-del*(fp-fptt)**2.
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      FUNCTION betacf(a,b,x)
      INTEGER MAXIT
      REAL betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER m,m2
      REAL aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a or b too big, or MAXIT too small in betacf'
1     betacf=h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.

      FUNCTION betai(a,b,x)
      REAL betai,a,b,x
CU    USES betacf,gammln
      REAL bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.


      SUBROUTINE pearsn(x,y,n,r,prob,z)
      INTEGER n
      REAL prob,r,z,x(n),y(n),TINY
      PARAMETER (TINY=1.e-20)
CU    USES betai
      INTEGER j
      REAL ax,ay,df,sxx,sxy,syy,t,xt,yt,betai
      ax=0.
      ay=0.
      do 11 j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
11    continue
      ax=ax/n
      ay=ay/n
      sxx=0.
      syy=0.
      sxy=0.
      do 12 j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
12    continue
      r=sxy/(sqrt(sxx*syy)+TINY)
      z=0.5*log(((1.+r)+TINY)/((1.-r)+TINY))
      df=n-2
      t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
      prob=betai(0.5*df,0.5,df/(df+t**2))
C     prob=erfcc(abs(z*sqrt(n-1.))/1.4142136)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
