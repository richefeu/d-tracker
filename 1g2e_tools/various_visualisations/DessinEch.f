c23456789012345678901234567890123456789012345678901234567890123456789012
c=======================================================================
      Program CourbePGPLOT

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      Integer PGOPEN

      IF (PGOPEN('?').LE.0) STOP

      Call LectureFichiers
      Call DessineEch
      Call Overlap

      CALL PGCLOS
      End
c=======================================================================
      Subroutine LectureFichiers

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      Character*20 Fichier

      Write(6,*) 'Fichier ?'
      Read(5,*) Fichier

      xrmax = 0.
      Open(1,file=Fichier)
      Read(1,*) npa
      Write(6,*)  'npa = ',npa
      Do i = 1,npa
         Read(1,*) x(1,i),x(2,i),xr(i),rota(i)
         xrmax = max(xrmax,xr(i))
      EndDo
      Read(1,*) z
c-----lecture des 4 coins
      Do i = 1,4
         Read(1,*) xp(i),yp(i)
      EndDo


      xmin = 0.
      xmax = 0.
      ymin = 0.
      ymax = 0.
      Do i = 1,npa
         xmin = min(xmin,xp(i)-xrmax)
         xmax = max(xmax,xp(i)+xrmax)
         ymin = min(ymin,yp(i)-xrmax)
         ymax = max(ymax,yp(i)+xrmax)
      EndDo

      CALL PGENV(xmin,xmax,ymin,ymax,1,0)

      Return
      End

c=======================================================================
      Subroutine DessineEch

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/par/ xp(4),yp(4),xnip(2,4)
      Real ppx(2),ppy(2)

      CALL PGSAVE
            
c-----Echantillon initial
      CALL PGSFS(1)
      CALL PGSCI(15)
      Do i = 1,npa
         Call PGCIRC(x(1,i),x(2,i),xr(i))
      EndDo
      CALL PGSFS(2)
      CALL PGSCI(1)
      Do i = 1,npa
         Call PGCIRC(x(1,i),x(2,i),xr(i))
      EndDo
c-----dessin du rayon
      Do i = 1,npa
         ppx(1) = x(1,i)
	 ppy(1) = x(2,i)
	 ppx(2) = x(1,i) + xr(i)*cos(rota(i))
	 ppy(2) = x(2,i) + xr(i)*sin(rota(i))
         CALL PGLINE(2,ppx,ppy)
      EndDo	 
c-----les 4 coins
      xrmax = 0.
      Do i = 1,npa
         xrmax = max(xrmax,xr(i))
      EndDo
      CALL PGSCI(2)
      Do i = 1,4
         ppx(1) = xp(i) - xrmax
         ppx(2) = xp(i) + xrmax
         ppy(1) = yp(i) !- xrmax
         ppy(2) = yp(i) !+ xrmax   
         CALL PGLINE(2,ppx,ppy)
         ppx(1) = xp(i) !- xrmax
         ppx(2) = xp(i) !+ xrmax
         ppy(1) = yp(i) - xrmax
         ppy(2) = yp(i) + xrmax   
         CALL PGLINE(2,ppx,ppy)
      EndDo	 


      Return
      End
c=======================================================================
      Subroutine Overlap

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Real ppx(2),ppy(2)

      ncont = 0
      violmin = 0.
      Do i = 1,npa
         Do j = i+1,npa
	    xij = x(1,j)-x(1,i)
	    yij = x(2,j)-x(2,i)
	    dij = sqrt(xij*xij+yij*yij)
	    hij = dij - xr(i) - xr(j)
	    If (hij.LT.0.) Then
	       ncont = ncont + 1
	       violmin = min(violmin,hij)
	    EndIf
	 EndDo      
      EndDo
      Write(6,*) "ncont, violmin = ",ncont,violmin 

      Return
      End
c=======================================================================
      Subroutine DessineForces

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      Integer ifn

      If (itype.EQ.0) Then
c--------Nbre de cellules pour la recherche des voisins
         Mlh = sqrt(dfloat(npa)/10.d0)
         If (Mlh.LT.3) Stop 'Pas assez de cellules'
         Ncell = Mlh*Mlh
         Write(6,*) 'Nbre de cel pour recherche voisins = ',Ncell
         Call Cellules_Voisines
         Call Recherche_Voisins3
         Call DetectForces
      EndIf
     
      fnmoy = 0.
      Do il = 1,ntotal
         fnmoy = fnmoy + fn(il)
      EndDo
      fnmoy = fnmoy / float(ntotal)
      Do il = 1,ntotal
         fn(il) = fn(il) / fnmoy
      EndDo

      Write(6,*) 'ncont,ntotal = ',ncont,ntotal
      Write(6,*) 'force moyenne = ',fnmoy

      fnmax = 0.
      fnmin = 100.
      Do il = 1,ntotal
         fnmax = max(fnmax,fn(il))
         fnmin = min(fnmin,fn(il))
      EndDo
      Write(6,*) 'fnmin, fnmax = ',fnmin,fnmax
     

c-----Chaine de forces
      CALL PGSAVE
      Call PGSCI(7)
      Do il = 1,ncont
         i = ior(il)
         j = iex(il)
         px(1) = x(1,i)
         px(2) = x(1,j)
         py(1) = x(2,i)
         py(2) = x(2,j)
         a = (20.-1.) / (fnmax - fnmin)
         b = 1. - a * fnmin
         ifn = a * fn(il) + b
         CALL PGSLW(ifn)
c         If (fn(il).GE.1) Then
c            Call PGSCI(7)
c         Else
            CALL PGSCI(2)
c         EndIf
         CALL PGLINE(2,px,py)
      EndDo
      Do il = ncont+1,ntotal
         i = ior(il)
         j = iex(il)
         px(1) = x(1,i)
         py(1) = x(2,i)
         px(2) = x(1,i) - xr(i)*xnip(1,j)
         py(2) = x(2,i) - xr(i)*xnip(2,j)
         a = (20.-1.) / (fnmax - fnmin)
         b = 1. - a * fnmin
         ifn = a * fn(il) + b
         CALL PGSLW(ifn)
c         If (fn(il).GE.1) Then
c            Call PGSCI(7)
c         Else
            CALL PGSCI(2)
c         EndIf
         CALL PGLINE(2,px,py)
      EndDo
         
         
      CALL PGUNSA

      Return
      End
c=======================================================================
      Subroutine Recherche_Voisins3

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      Real*8 x0(2,npam)
      Real*8 xd1,xd2,xmax,ymax,Celli

      rmax = 0.d0
      Do i = 1,npa
         rmax = max(rmax,xr(i))
      EndDo
           
      xd1 = 1.15d0*2.*rmax
      xd2 = xd1*xd1

      xmax = 0.d0
      ymax = 0.d0
      Do i = 1,npa
         xmax = max(xmax,abs(x(1,i)))
         ymax = max(ymax,abs(x(2,i)))
      EndDo
      Do i = 1,npa
         x0(1,i) = x(1,i) / xmax / 2.001d0 !position des particules
         x0(2,i) = x(2,i) / ymax / 2.001d0 !entre -0.5 et 0.5 (x et y)
      EndDo
      Do i = 1,Ncell
         head(i) = 0
      EndDo
      Celli = dfloat(Mlh)
      Do i = 1,npa
         icell = 1 + Int((x0(1,i) + 0.5d0)*Celli)
     &             + Int((x0(2,i) + 0.5d0)*Celli)*Mlh
         List(i) = Head(icell)
         Head(icell) = i
      EndDo

c-----Recherche des voisins entre les Particules
      npoint = 0
      Do icell = 1,Ncell
         i = Head(icell)
1        If (i.GT.0) Then
            j = List(i)
2           If (j.GT.0) Then
                dx = x(1,j)-x(1,i)
                dy = x(2,j)-x(2,i)
                dnk = dx*dx + dy*dy
                If (dnk.LT.xd2) Then
                   npoint = npoint + 1
                   Liste(1,npoint) = i
                   Liste(2,npoint) = j
                EndIf
                j = List(j)
                GoTo 2
            EndIf
            jcell0 = 4 * (icell-1)
            Do NaBord = 1,4
               jcell = Map(jcell0+NaBord)
               j = Head(jcell)
3              If (j.NE.0) Then
                  dx = x(1,j)-x(1,i)
                  dy = x(2,j)-x(2,i)
                  dnk = dx*dx + dy*dy
                  If (dnk.LT.xd2) Then
                     npoint = npoint + 1
                     Liste(1,npoint) = i
                     Liste(2,npoint) = j
                  EndIf
                  j = List(j)
                  GoTo 3
               EndIf
            EndDo
            i = List(i)
            GoTo 1
         EndIf
      EndDo
      nl = npoint
 
c-----Recherche des voisins entre echantillon et paroi
      Do i = 1,npa
         xi = x(1,i)
         yi = x(2,i)
         Do j = 1,4
            pos = xnip(1,j)*xp(j)+xnip(2,j)*yp(j)
            dij = abs(pos-xnip(1,j)*xi-xnip(2,j)*yi)
            If (dij.LE.xd1) Then
               npoint = npoint + 1
               Liste(1,npoint) = i
               Liste(2,npoint) = j
            EndIf
         EndDo
      EndDo
      nlt = npoint !nbre total de couples de voisins 
	
      Return
      End
c=======================================================================

      Subroutine Cellules_Voisines

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      Do iy = 1,Mlh
         Do ix = 1,Mlh
            imap = (NumCell(ix,iy,Mlh) - 1)*4
            map(imap+1) = NumCell(ix+1,iy  ,Mlh)
            map(imap+2) = NumCell(ix+1,iy+1,Mlh)
            map(imap+3) = NumCell(ix  ,iy+1,Mlh)
            map(imap+4) = NumCell(ix-1,iy+1,Mlh)
         EndDo
      EndDo

      Return
      End
c=======================================================================

      Integer Function NumCell(ik,jk,MM)

      Integer ik,jk,MM

      NumCell = 1 + Mod(ik - 1 + MM,MM) + Mod(jk - 1 + MM,MM)*MM

      Return
      End
c=======================================================================
      Subroutine DetectForces

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam,nlm=30*npam)
      Common/pos/ x(2,npam),xr(npam),rota(npam)
      Common/ori/ ior(nvois),iex(nvois),fn(nvois)
      Common/list/ Liste(2,nlm),itype
      Common/nbr/ npoint,npa,ncont,ntotal,nl,nlt
      Common/des/ px(2),py(2)
      Common/par/ xp(4),yp(4),xnip(2,4)
      Parameter(Mmax=100,NCellMax=Mmax*Mmax,MapSize=4*NCellMax)
      Integer List,Head,Map
      Common /Voisins1/ List(npam),Head(NCellMax),Map(MapSize)
      Common /Voisins2/ Mlh,Ncell

      kn = 1

      ncont = 0
      Do il = 1,nl
         i = Liste(1,il)
         j = Liste(2,il)
         xij = x(1,j) - x(1,i)
         yij = x(2,j) - x(2,i)
         dij = xij*xij + yij*yij
         xri = xr(i)
         xrj = xr(j)
         xr2 = (xri+xrj)*(xri+xrj)
         If (dij.LT.xr2) Then
            ncont = ncont + 1
            ior(ncont) = i
            iex(ncont) = j
            dij = sqrt(dij)
            hij = dij - xri - xrj
            fn(ncont) = - hij*kn
         EndIf
      EndDo
      ntotal = ncont
      Do il = nl+1,nlt
         i = Liste(1,il)
         j = Liste(2,il)
         xi = x(1,i)
         yi = x(2,i)
         pos = xnip(1,j)*xp(j)+xnip(2,j)*yp(j)
         hij = abs(pos-xnip(1,j)*xi-xnip(2,j)*yi) - xr(i)
         If (hij.LT.0.) Then
            ntotal = ntotal + 1
            ior(ntotal) = i
            iex(ntotal) = j
            fn(ntotal) = - hij*kn
         EndIf
      EndDo
      Write(6,*) nl,nlt

      Return
      End

c=======================================================================
c=======================================================================

      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END

c=======================================================================
