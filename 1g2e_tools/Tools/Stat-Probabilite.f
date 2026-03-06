      Program statisticsONforces

      Implicit Real*8 (a-h,o-z)

      Parameter(Lmax=5000000)
      Character*50 Fichier,filename
      Character*6 lign
      Real*8 fn(Lmax)
      Real*8 hist(1000000,6)

      pi = 4.d0*atan(1.d0)

      ncpt = 0

      Write(6,*)'Nom du fichier ?'
      Read(5,*) Fichier
      Open(1,File=Fichier,status='old')

      Write(6,*) 'lecture des valeurs de fluct avec'
      Write(6,*) ' un seuil sur les deplacements en pixels'
      Write(6,*) ' la colonne 3 du fichier fluc contient les dep'
      write(6,*) ' en pixels'
      Write(6,*) 'depmin, depmax ?'
      Read(5,*) depmin,depmax
      
      ntotal = 0
      nval = 0
1     Read(1,*,end=2) z,zz,zz
          ntotal = ntotal + 1
          If ( (z.GE.depmin).AND.(z.LE.depmax) ) Then
             nval = nval + 1
          EndIf
          GoTo 1
2     Rewind(1)
c      Write(6,*) 'ntotal = ',ntotal
c      Write(6,*) 'nval = ',nval
      ilk = 0
      Write(6,*) 'Symetrisation des donnees (0=oui, 1=non)?'
      Read(5,*) irep
      Do il = 1,ntotal
         Read(1,*) z1,z2,z3
         If ( (z1.GE.depmin).AND.(z1.LE.depmax) ) Then
            ilk = ilk + 1
            fn(ilk) = (z1)
            If (irep.EQ.0) fn(ilk+nval) = -(z1)
         EndIf
      EndDo  
      ntotal = nval
      If (irep.EQ.0) ntotal = 2*nval
      Write(6,*) 'ntotal = ',ntotal
      
       
      Close(1)

      rmoy = 0.
      Do i = 1,ntotal
         rmoy = rmoy + fn(i)
      EndDo
      rmoy = rmoy / dfloat(ntotal)
      Write(6,*) 'rmoy = ',rmoy
c      Do i = 1,ntotal
c         fn(i) = fn(i) / rmoy
c      EndDo
         
      Call Sort(ntotal,fn)
      Do i = 1,ntotal
         Write(34,*) fn(i),float(i)/float(ntotal)
      EndDo

      Call moment(fn,ntotal,ave,adev,sdev,var,skew,curt)
      Write(6,*) 'skew,curt =  ',skew,curt
      Open(62,file='kurto',access='append')
      Write(62,*) curt
      Close(62)

11    Write(6,*) '--------------------------------------------'
      Write(6,*) 'Which way to compute PDF ?'
      Write(6,*) ' 1- number of beams is fixed '
      Write(6,*) ' 2- number of value per beam is fixed'
      Write(6,*) '- - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*) 'INFO: number of values: ',ntotal
      Write(6,*) '- - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*) 'answer 1 or answer 2 ?'
      Read(5,*) ians
      If ((ians.LT.1).OR.(ians.GT.2)) GoTo 11
      If (ians.EQ.2) Then
12       Write(6,*) 'number of values per beams ?'
         Read(5,*) nval
         nnn = ntotal - nval*(ntotal/nval)
         Write(6,*) 'number of values in the last beam: ',nnn
         Write(6,*) 'OK ? (OK = 0)'
         Read(5,*) irep 
         If (irep.NE.0) GoTo 12
         Call StatForces(ntotal,nval,fn,hist,ncl)
      Else 
13       Write(6,*) 'Number of beams ?'
         Read(5,*) nval
         Write(6,*) 'Min number of values per beams ?'
         Read(5,*) nvalbeam
         Call StatForces3(ntotal,nval,nvalbeam,fn,hist,ncl,iflag)
         If (iflag.EQ.0) GoTo 13
      EndIf
      Write(6,*) 'output file ?'
      ReaD(5,*) filename
      Open(66,file=filename)
      rv = 0.
      Do i = 1,ncl
         If (hist(i,1).GT.1) Then
         Write(66,*) hist(i,2),
     &           (hist(i,1)/hist(i,3)/dfloat(ntotal)),hist(i,4),
     &           hist(i,1)
         EndIf
      EndDo
      Close(66)
 
      End
	
c=======================================================================
      
      Subroutine StatForces2(nc,nval,xx,hist,ncl,iflag)
      Implicit Real*8 (a-h,o-z)
      Integer nc,nval,nclasse,ideb,ifin
      Real*8 xx(nc),hist(1000000,6)

c-----ncl = number of beams
c     in this subroutine, nval = ncl
c
      ncl = nval

      xmin = 1d10
      xmax = -1d10
      Do i = 1,nc
         xmin = min(xmin,xx(i))
         xmax = max(xmax,xx(i))
      EndDo
      dxx = xmax - xmin
      dx = dxx / dfloat(ncl)
      xfin = xmin
      Write(6,*) 'xmin = ',xmin
      Write(6,*) 'xmax = ',xmax
      Write(6,*) 'dx = ',dx
      Write(6,*) 'OK (1 = no, other = ok)?'
      Read(5,*) ii
      If (ii.EQ.1) Then
         iflag = 0
         Return
      EndIf
      
      Do il = 1,ncl
         xdeb = xfin
         xfin = xdeb + dx
         rmoy = 0.d0
         ncpt = 0
         i = ideb
         Do i = 1,nc
            If ((xx(i).GE.xdeb).AND.(xx(i).LT.xfin)) Then
               rmoy = rmoy + xx(i)
               ncpt = ncpt + 1
            EndIf
         EndDo
         rmoy = rmoy / dfloat(ncpt)
         rvar = 0.d0
         ncpt = 0
         Do i = 1,nc
            If ((xx(i).GE.xdeb).AND.(xx(i).LT.xfin)) Then
               rvar = rvar + (xx(i)-rmoy)*(xx(i)-rmoy)
               ncpt = ncpt + 1
            EndIf
         EndDo
         rvar = sqrt(rvar/dfloat(ncpt))         
         hist(il,1) = ncpt ! Nbre de valeurs.
         hist(il,2) = rmoy ! Centre de la classe.
         hist(il,3) = dx   ! Largeur de la classe.
         hist(il,4) = rvar ! Variance sur la classe.
         hist(il,5) = xdeb
         hist(il,6) = xfin
      EndDo

      iflag = 1

      Return
      End

c=======================================================================
      
      Subroutine StatForces3(nc,nval,nvalbeam,xx,hist,ncl,iflag)
      Implicit Real*8 (a-h,o-z)
      Integer nc,nval,nclasse,ideb,ifin
      Real*8 xx(nc),hist(1000000,6)

c-----ncl = number of beams
c     in this subroutine, nval = ncl
c
      ncl = nval

      xmin = 1d10
      xmax = -1d10
      Do i = 1,nc
         xmin = min(xmin,xx(i))
         xmax = max(xmax,xx(i))
      EndDo
      dxx = xmax - xmin
      dx = dxx / dfloat(ncl)
      dx0 = dx
      xfin = xmin
      Write(6,*) 'xmin = ',xmin
      Write(6,*) 'xmax = ',xmax
      Write(6,*) 'dx = ',dx
      Write(6,*) 'OK (1 = no, other = ok)?'
      Read(5,*) ii
      If (ii.EQ.1) Then
         iflag = 0
         Return
      EndIf
      
      ideb = 1
      Do il = 1,ncl
         dx = dx0
         xdeb = xfin
1        Continue
         xfin = xdeb + dx
         rmoy = 0.d0
         ncpt = 0
         i = ideb
c        -- on compte le nombre de valeur dans la classe
2        If (i.LE.nc) Then
            If ((xx(i).GE.xdeb).AND.(xx(i).LT.xfin)) Then
               ncpt = ncpt + 1
               i = i + 1
               GoTo 2
            EndIf
         EndIf
         If (ncpt.LT.nvalbeam) Then ! le nbre de val est trop faible
             dx = dx + dx0
             xxx = xdeb + dx
             If (xxx.LT.xmax) GoTo 1
         EndIf
         ncpt = 0
         i = ideb
3        If (i.LE.nc) Then
            If ((xx(i).GE.xdeb).AND.(xx(i).LT.xfin)) Then
               rmoy = rmoy + xx(i)
               ncpt = ncpt + 1
               i = i + 1
               GoTo 3
            EndIf
         EndIf
         rmoy = rmoy / dfloat(ncpt)
         rvar = 0.d0
         ncpt = 0
         i = ideb
4        If (i.LE.nc) Then
            If ((xx(i).GE.xdeb).AND.(xx(i).LT.xfin)) Then
               rvar = rvar + (xx(i)-rmoy)*(xx(i)-rmoy)
               ncpt = ncpt + 1
               i = i + 1
               GoTo 4
            EndIf
         EndIf
         ideb = i
         rvar = sqrt(rvar/dfloat(ncpt))         
         hist(il,1) = ncpt ! Nbre de valeurs.
         hist(il,2) = rmoy ! Centre de la classe.
         hist(il,3) = dx   ! Largeur de la classe.
         hist(il,4) = rvar ! Variance sur la classe.
         hist(il,5) = xdeb
         hist(il,6) = xfin
      EndDo

      iflag = 1

      Return
      End

c=======================================================================

      Subroutine StatForces(nc,nval,xx,hist,ncl)
      Implicit Real*8 (a-h,o-z)
      Integer nc,nval,nclasse,ideb,ifin
      Real*8 xx(nc),hist(1000000,6)

c-----On calcule un histogramme ou l'on veut nval valeurs par classes
      nclasse = nc/nval
      If (mod(nc,nval).NE.0) nclasse = nclasse + 1
      ideb = 1
      ifin = nval
      Do ncl = 1,nclasse-1
c        --Moyenne.
         rmoy = 0.
         ncpt = 0
         Do i = ideb,ifin
            rmoy = rmoy + xx(i)
            ncpt = ncpt + 1
         EndDo
         rmoy = rmoy/float(ncpt)
c        --Variance.
         var = 0.
         Do i = ideb,ifin
            var = var + (xx(i)-rmoy)*(xx(i)-rmoy)
         EndDo
         var = sqrt(var/float(ncpt))
         If (ideb.EQ.1) Then
            Pas = (xx(ifin)+xx(ifin+1))*0.5-xx(ideb) 
         Else
            Pas = (xx(ifin)+xx(ifin+1))*0.5-(xx(ideb)+xx(ideb-1))*0.5
         EndIf
c        --stockage des valeurs calculees.
         hist(ncl,1) = ncpt ! Nbre de valeurs.
         hist(ncl,2) = rmoy ! Centre de la classe.
         hist(ncl,3) = Pas ! Largeur de la classe.
         hist(ncl,4) = var ! Variance sur la classe.
         hist(ncl,5) = xx(ideb)
         hist(ncl,6) = xx(ifin)
         If (ncl.NE.nclasse-1) Then
            ideb = ifin + 1
            ifin = ideb + nval - 1
         EndIf
      EndDo
c
      ncl = nclasse
      ideb = ifin + 1
      ifin = nc
c     --Moyenne.
      rmoy = 0.
      ncpt = 0
      Do i = ideb,ifin
         rmoy = rmoy + xx(i)
         ncpt = ncpt + 1
      EndDo
      rmoy = rmoy/dfloat(ncpt)
c     --Variance.
      var = 0.
      Do i = ideb,ifin
         var = var + (xx(i)-rmoy)*(xx(i)-rmoy)
      EndDo
      var = sqrt(var/dfloat(ncpt))
      Pas = xx(ifin) - (xx(ideb)+xx(ideb-1))*0.5

c     --stockage des valeurs calculees.
      hist(ncl,1) = (ifin-ideb)+1
      hist(ncl,2) = rmoy
      hist(ncl,3) = Pas
      hist(ncl,4) = var
      hist(ncl,5) = xx(ideb)
      hist(ncl,6) = xx(ifin)

      Return
      End 

c=======================================================================
      Subroutine Sort(n,ra)

      Real*8 ra(n)

      l = n/2 + 1
      ir = n

10    Continue
      If (L.GT.1) Then
      l = l - 1
      rra = ra(L)
      Else
      rra = ra(ir)
      ra(ir) = ra(1)
      ir = ir-1
      If (ir.EQ.1) Then
      ra(1) = rra
      Return
      EndIf
      EndIf
      i = l
      j = l + l
20    If (j.LE.ir) Then
      If (j.LT.ir) Then
      If (ra(j).Lt.ra(j+1)) j = j + 1
      EndIf
      If (rra.LT.ra(j)) Then
      ra(i) = ra(j)
      i = j
      j = j + j
      Else
      j = ir+1
      EndIf
      GoTo 20
      EndIf
      ra(i) = rra
      GoTo 10	
      
      End	
c=======================================================================

      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
      Implicit Real (a-h,o-z)
      INTEGER n
      REAL*8 adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL*8 p,s,ep
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      stotal = 0.d0
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
        pause 'no skew or kurtosis when zero variance in moment'
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
c======================================================================
