c=======================================================================
      Program Deformation_Locale_Triangulation_Delaunay
c     http://math.nist.gov/~JBernal/JBernal_Sft.html
      Implicit Real*8 (a-h,o-z)

      Character*20 ConfIni,ConfFin
      Parameter(npam=30000)
      Dimension lt(3,npam)
      Dimension ltg(npam,20),ntg(npam)
      Dimension r1(2,npam),r2(2,npam)
      Dimension rr1(npam),rot1(npam),rr2(npam),rot2(npam)
      Dimension rt1(2,npam),rt2(2,npam)
      Integer jgauch(npam),num(npam),num0(npam)
      Dimension vpn(2,npam)
      Integer PGOPEN
      
      Pi = 4.*atan(1.)
c-----Initialisation de la fenetre graphique PGPLOT
      If (PGOPEN('?') .LE. 0) STOP

c-----Lecture des Fichiers CONF INI et FIN
      Write(6,*) 'Fichier CONF initial ?'
      Read(5,*) ConfIni
      Write(6,*) 'Fichier CONF final ?'
      Read(5,*) ConfFin

      Open(1,file=ConfIni,status='old',err=901)
      Open(2,file=ConfFin,status='old',err=902)
      Read(1,*) npa !,nb
      Read(2,*) npa !,nb
      nb = npa
      Do i = 1,nb
         Read(1,*) r1(1,i),r1(2,i),rr1(i),rot1(i)
         Read(2,*) r2(1,i),r2(2,i),rr2(i),rot2(i)
      EndDo
      Close(1)
      Close(2)
c     -- recentrage des echantillons
      xmin = 1d7
      xmax = -1d7
      ymin = 1d7
      ymax = -1d7
      Do k = 1,nb
         xmin = min(xmin,r1(1,k))
         xmax = max(xmax,r1(1,k))
         ymin = min(ymin,r1(2,k))
         ymax = max(ymax,r1(2,k))
      EndDo
      dx = (xmin+xmax)*0.5
      dy = (ymin+ymax)*0.5
      Do k = 1,nb
         r1(1,k) = r1(1,k) - dx
         r1(2,k) = r1(2,k) - dy
      EndDo   
      xmin = 1d7
      xmax = -1d7
      ymin = 1d7
      ymax = -1d7
      Do k = 1,nb
         xmin = min(xmin,r2(1,k))
         xmax = max(xmax,r2(1,k))
         ymin = min(ymin,r2(2,k))
         ymax = max(ymax,r2(2,k))
      EndDo
      dx = (xmin+xmax)*0.5
      dy = (ymin+ymax)*0.5
      Do k = 1,nb
         r2(1,k) = r2(1,k) - dx
         r2(2,k) = r2(2,k) - dy
      EndDo   

      iajout = 1


c     --rectangle avec ajout de pts, iajout = 0
      If (iajout.EQ.0) Then
         Call AjoutePoints(nb,r1,rot1,np,rt1)
         Call AjoutePoints(nb,r2,rot2,np,rt2)
         Write(6,*) 'nb = ',nb
         Write(6,*) 'np = ',np           
c        --Triangulation de Delaunay
         Call Triangule(np,rt1,ivnxt,lt)
         Call CorrectionDimRectangle(nb,r1,rot1,np,rt1)
         Call CorrectionDimRectangle(nb,r2,rot2,np,rt2)
c        --Reduction du nombre triangles
         Call RetireTrianglesPlats(np,rt1,ivnxt,lt)
         Call RetireTrianglesRectangle(ivnxt,lt)
         Call CalculDef(np,rt1,rt2,ivnxt,lt,nb,r2,rot2)
      Else
c     --rectangles sans ajout de pts, iajout = 1
        np = nb
c        --Triangulation de Delaunay
        Do i = 1,np
	   r1(1,i) = r1(1,i)*100.
	   r1(2,i) = r1(2,i)*100.
	EndDo
        Call Triangule(np,r1,ivnxt,lt)
        Do i = 1,np
	   r1(1,i) = r1(1,i)/100.
	   r1(2,i) = r1(2,i)/100.
	EndDo
c        --Reduction du nombre triangles
        Call RetireTrianglesPlats(np,r1,ivnxt,lt)
        Call CalculDef(np,r1,r2,ivnxt,lt,nb,r2,rot2)
      EndIf
      Stop
901   Write(6,*) 'PB fichier '//ConfIni
      Stop
902   Write(6,*) 'PN fichier '//ConfFin

c-----Nombre et numero des triangles pour chaque grain
c      Call CompteTrianglesGrain(np,ltg,ntg,ivnxt,lt)

      End
c=======================================================================
      
      Subroutine AjoutePoints(npa,r,rot,nn,rt)
      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension r(2,npam),rt(2,npam),rot(npam)

c-----Les points ajoutes correspondent aux sommets des
c     rectangles.

      Pi = 4.d0*atan(1.d0)

c-----geometrie des rectangles
c       hyp : rectangle l = 0.02 m et L = 0.055 m 
      zl = 0.015
      zll = 0.050
      rl = sqrt((zl/2.)**2.+(zll/2.)**2.)
      ra = acos(zll/2./rl)
c-----On ajoute les nn points
      nn = 0
      Do i = 1,npa
         nn = nn + 1
         rt(1,nn) = r(1,i) !centre de masse X
         rt(2,nn) = r(2,i) !centre de masse Y
         rott = rot(i)*Pi/180.
         rt(1,nn+1) = r(1,i) + rl*cos(-ra+rott) ! coin 1 - X
         rt(2,nn+1) = r(2,i) + rl*sin(-ra+rott) ! coin 1 - Y
         rt(1,nn+2) = r(1,i) + rl*cos(ra+rott)  ! coin 2 - X
         rt(2,nn+2) = r(2,i) + rl*sin(ra+rott)  ! coin 2 - Y
         rt(1,nn+3) = r(1,i) + rl*cos(Pi-ra+rott) ! coin 3 - X
         rt(2,nn+3) = r(2,i) + rl*sin(Pi-ra+rott) ! coin 3 - Y
         rt(1,nn+4) = r(1,i) + rl*cos(Pi+ra+rott) ! coin 4 - X
         rt(2,nn+4) = r(2,i) + rl*sin(Pi+ra+rott) ! coin 4 - Y

         rt(1,nn+5) = (rt(1,nn+2)+rt(1,nn+3))*0.5 ! milieu grand cote 2-3 X
         rt(2,nn+5) = (rt(2,nn+2)+rt(2,nn+3))*0.5 ! milieu grand cote 2-3 Y

         rt(1,nn+6) = (rt(1,nn+4)+rt(1,nn+1))*0.5 ! milieu grand cote 4-1 X
         rt(2,nn+6) = (rt(2,nn+4)+rt(2,nn+1))*0.5 ! milieu grand cote 4-1 Y

         rt(1,nn+7) = (rt(1,nn+5)+rt(1,nn+2))*0.5
         rt(2,nn+7) = (rt(2,nn+5)+rt(2,nn+2))*0.5

         rt(1,nn+8) = (rt(1,nn+5)+rt(1,nn+3))*0.5
         rt(2,nn+8) = (rt(2,nn+5)+rt(2,nn+3))*0.5

         rt(1,nn+9) = (rt(1,nn+6)+rt(1,nn+4))*0.5
         rt(2,nn+9) = (rt(2,nn+6)+rt(2,nn+4))*0.5

         rt(1,nn+10) = (rt(1,nn+6)+rt(1,nn+1))*0.5
         rt(2,nn+10) = (rt(2,nn+6)+rt(2,nn+1))*0.5

         rt(1,nn+11) = (rt(1,nn+1)+rt(1,nn+2))*0.5
         rt(2,nn+11) = (rt(2,nn+1)+rt(2,nn+2))*0.5

         rt(1,nn+12) = (rt(1,nn+3)+rt(1,nn+4))*0.5
         rt(2,nn+12) = (rt(2,nn+3)+rt(2,nn+4))*0.5

         rt(1,nn+13) = (rt(1,nn)+rt(1,nn+11))*0.5
         rt(2,nn+13) = (rt(2,nn)+rt(2,nn+11))*0.5

         rt(1,nn+14) = (rt(1,nn)+rt(1,nn+12))*0.5
         rt(2,nn+14) = (rt(2,nn)+rt(2,nn+12))*0.5

         nn = nn + 14
      EndDo

      Return
      End
c=======================================================================
      
      Subroutine CorrectionDimRectangle(npa,r,rot,nn,rt)
      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension r(2,npam),rt(2,npam),rot(npam)

c-----Les points ajoutes correspondent aux sommets des
c     rectangles.

      Pi = 4.d0*atan(1.d0)

c-----geometrie des rectangles
c       hyp : rectangle l = 0.02 m et L = 0.055 m 
      zl = 0.02
      zll = 0.055
      rl = sqrt((zl/2.)**2.+(zll/2.)**2.)
      ra = acos(zll/2./rl)
c-----On ajoute les nn points
      nn = 0
      Do i = 1,npa
         nn = nn + 1
         rt(1,nn) = r(1,i) !centre de masse X
         rt(2,nn) = r(2,i) !centre de masse Y
         rott = rot(i)*Pi/180.
         rt(1,nn+1) = r(1,i) + rl*cos(-ra+rott) ! coin 1 - X
         rt(2,nn+1) = r(2,i) + rl*sin(-ra+rott) ! coin 1 - Y
         rt(1,nn+2) = r(1,i) + rl*cos(ra+rott)  ! coin 2 - X
         rt(2,nn+2) = r(2,i) + rl*sin(ra+rott)  ! coin 2 - Y
         rt(1,nn+3) = r(1,i) + rl*cos(Pi-ra+rott) ! coin 3 - X
         rt(2,nn+3) = r(2,i) + rl*sin(Pi-ra+rott) ! coin 3 - Y
         rt(1,nn+4) = r(1,i) + rl*cos(Pi+ra+rott) ! coin 4 - X
         rt(2,nn+4) = r(2,i) + rl*sin(Pi+ra+rott) ! coin 4 - Y

         rt(1,nn+5) = (rt(1,nn+2)+rt(1,nn+3))*0.5 ! milieu grand cote 2-3 X
         rt(2,nn+5) = (rt(2,nn+2)+rt(2,nn+3))*0.5 ! milieu grand cote 2-3 Y

         rt(1,nn+6) = (rt(1,nn+4)+rt(1,nn+1))*0.5 ! milieu grand cote 4-1 X
         rt(2,nn+6) = (rt(2,nn+4)+rt(2,nn+1))*0.5 ! milieu grand cote 4-1 Y

         rt(1,nn+7) = (rt(1,nn+5)+rt(1,nn+2))*0.5
         rt(2,nn+7) = (rt(2,nn+5)+rt(2,nn+2))*0.5

         rt(1,nn+8) = (rt(1,nn+5)+rt(1,nn+3))*0.5
         rt(2,nn+8) = (rt(2,nn+5)+rt(2,nn+3))*0.5

         rt(1,nn+9) = (rt(1,nn+6)+rt(1,nn+4))*0.5
         rt(2,nn+9) = (rt(2,nn+6)+rt(2,nn+4))*0.5

         rt(1,nn+10) = (rt(1,nn+6)+rt(1,nn+1))*0.5
         rt(2,nn+10) = (rt(2,nn+6)+rt(2,nn+1))*0.5

         rt(1,nn+11) = (rt(1,nn+1)+rt(1,nn+2))*0.5
         rt(2,nn+11) = (rt(2,nn+1)+rt(2,nn+2))*0.5

         rt(1,nn+12) = (rt(1,nn+3)+rt(1,nn+4))*0.5
         rt(2,nn+12) = (rt(2,nn+3)+rt(2,nn+4))*0.5

         rt(1,nn+13) = (rt(1,nn)+rt(1,nn+11))*0.5
         rt(2,nn+13) = (rt(2,nn)+rt(2,nn+11))*0.5

         rt(1,nn+14) = (rt(1,nn)+rt(1,nn+12))*0.5
         rt(2,nn+14) = (rt(2,nn)+rt(2,nn+12))*0.5

         nn = nn + 14
      EndDo

      Return
      End
c=======================================================================

      Subroutine CompteTrianglesGrain(npa,ltg,ntg,nt,lt)
      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension lt(3,npam)
      Dimension ltg(npam,20),ntg(npam)

      Do i = 1,npa
         ntg(i) = 0
      EndDo
      Do i = 1,nt
         i1 = lt(1,i)
         i2 = lt(2,i)
         i3 = lt(3,i)
         ntg(i1) = ntg(i1) + 1
         ntg(i2) = ntg(i2) + 1
         ntg(i3) = ntg(i3) + 1
         ltg(i1,ntg(i1)) = i
         ltg(i2,ntg(i2)) = i
         ltg(i3,ntg(i3)) = i
      EndDo
      
      Return
      End 
c=======================================================================

      Subroutine CalculDef(np,x1,x2,nt,lt,npa,r,rot)
      Implicit Real*8 (a-h,o-z)

      Parameter(npam=30000)
      Dimension lt(3,npam)
      Dimension x1(2,npam),x2(2,npam)
      Dimension u(2,npam)
      Dimension Exx(npam),Eyy(npam),Exy(npam),Ew(npam),Dd(npam)
      Dimension E1(npam),E2(npam),Dis(npam),Vol(npam)
      Dimension r(2,npam),rot(npam)
      Real*4 rx(4),ry(4)
      Real*8 Surf,nx,ny,SS(npam)
      Real*8 xmax,ymax,xmin,ymin
      Real*4 rxmax,rymax,rxmin,rymin,xp(3),yp(3),xx,yy,size
      Real*4 fDmin,fDmax
      Real*4 BRIGHT,CONTRA
      Character*50 FichierPictures
 
      Write(6,*) 'Calcul des deformations...'

      Do i = 1,np
         Exx(i) = 0.
         Eyy(i) = 0.
         Exy(i) = 0.
         Ew(i)  = 0.
         E1(i)  = 0.
         E2(i)  = 0.
         Dis(i) = 0.
         Vol(i) = 0.
         u(1,i) = x1(1,i) - x2(1,i)
         u(2,i) = x1(2,i) - x2(2,i)
         Write(95,*) x1(1,i),x1(2,i)
         Write(96,*) x2(1,i),x2(2,i)
      EndDo
 
      Do i = 1,nt
         i01 = lt(1,i)
         i02 = lt(2,i)
         i03 = lt(3,i)
c        --numerotation des sommets dans le sens inverse trigo
c        --i1 est le sommet completement a gauche
c        --i2 est le sommet le plus haut entre i2 et i3
c        --i3 est le sommet le plus bas entre i2 et i3
         xi1 = x1(1,i01)
         xi2 = x1(1,i02)
         xi3 = x1(1,i03)
         xmin = min(xi1,xi2,xi3)
         If (xmin.EQ.xi1) i1 = i01
         If (xmin.EQ.xi2) i1 = i02
         If (xmin.EQ.xi3) i1 = i03
         If (i1.EQ.i01) Then
            xi2 = x1(2,i02)
            xi3 = x1(2,i03)
            If (xi2.GT.xi3) Then
               i2 = i02
               i3 = i03
            Else
               i2 = i03
               i3 = i02
            EndIf
         EndIf
         If (i1.EQ.i02) Then
            xi2 = x1(2,i01)
            xi3 = x1(2,i03)
            If (xi2.GT.xi3) Then
               i2 = i01
               i3 = i03
            Else
               i2 = i01
               i3 = i03
            EndIf
         EndIf
         If (i1.EQ.i03) Then
            xi2 = x1(2,i01)
            xi3 = x1(2,i02)
            If (xi2.GT.xi3) Then
               i2 = i01
               i3 = i02
            Else
               i2 = i02
               i3 = i01
            EndIf
         EndIf
c        --Surface de la cellule consideree
         SS(i) = 0.d0
         SS(i) = Surf(x1,i1,i2,i3)
c         --Calcul des deformations locales
         axx = 0.d0
         ayy = 0.d0
         ayx = 0.d0
         axy = 0.d0

         dx = x1(1,i2)-x1(1,i1)
         dy = x1(2,i2)-x1(2,i1)
         nx = - dy  ! normale sortante
         ny = dx 
         Ux = (u(1,i2)+u(1,i1))*0.5
         Uy = (u(2,i2)+u(2,i1))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy

         dx = x1(1,i3)-x1(1,i2)
         dy = x1(2,i3)-x1(2,i2)
         nx = - dy  ! normale sortante
         ny = dx 
         Ux = (u(1,i3)+u(1,i2))*0.5
         Uy = (u(2,i3)+u(2,i2))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy

         dx = x1(1,i1)-x1(1,i3)
         dy = x1(2,i1)-x1(2,i3)
         nx = - dy  ! normale sortante
         ny = dx 
         Ux = (u(1,i3)+u(1,i1))*0.5
         Uy = (u(2,i3)+u(2,i1))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy

         Exx(i) = axx / SS(i)
         Eyy(i) = ayy / SS(i)
         Exy(i) = 0.5d0*(axy+ayx) / SS(i)
         Ew(i) = 0.5d0*(ayx-axy) / SS(i)
         Delta = ((Exx(i)-Eyy(i))*0.5)**2.d0 + (Exy(i))**2.d0
         E1(i) = (Exx(i)+Eyy(i))*0.5d0 + sqrt(Delta)
         E2(i) = (Exx(i)+Eyy(i))*0.5d0 - sqrt(Delta)
         Dis(i) = E1(i) - E2(i)
         Vol(i) = -(E1(i) + E2(i)) 
      EndDo

      rxmin = 1.e10
      rymin = 1.e10
      rxmax = -1.e10
      rymax = -1.e10
      Do i = 1,np
c         rxmax = max(rxmax,x1(1,i),x2(1,i))
c         rymax = max(rymax,x1(2,i),x2(2,i))
c         rxmin = min(rxmin,x1(1,i),x2(1,i))
c         rymin = min(rymin,x1(2,i),x2(2,i))
         rxmax = max(rxmax,x2(1,i))
         rymax = max(rymax,x2(2,i))
         rxmin = min(rxmin,x2(1,i))
         rymin = min(rymin,x2(2,i))
      EndDo
      Write(6,*) rxmin,rxmax,rymin,rymax
      rxmax = max(abs(rxmin),abs(rymin),
     &            abs(rxmax),abs(rymax)
     &           )
c      rxmax = rxmax + 0.05
c      CALL PGENV(-rxmax,rxmax,-rxmax,rxmax,1,-2)
      CALL PGENV(rxmin,rxmax,rymin,rymax,1,-2)
c      CALL PGQVP(0,rx(1),rx(2),ry(1),ry(2))
c      Write(6,*) rx(1),rx(2),ry(1),ry(2)

1     Continue
      Write(6,*) '    ____________________________________________'
      Write(6,*) '    |    Min   |   Max    |   Moy    |   Macro  '
      Write(6,*) '________________________________________________'
      Call MinMaxMoy2(nt,Exy,Dmin,Dmax,Dmoy,Dvar,SS)
      Write(6,60) 'Exy |',Dmin,'|',Dmax,'|',Dmoy,'|   -'
      ExyT = Dmoy
      ExyVAR = Dvar
      Call MinMaxMoy2(nt,Ew,Dmin,Dmax,Dmoy,Dvar,SS)
      Write(6,60) ' Ew |',Dmin,'|',Dmax,'|',Dmoy,'|   -'
      EwT = Dmoy
      EwVAR = Dvar
      Call MinMaxMoy2(nt,Exx,Dmin,Dmax,Dmoy,Dvar,SS)
      Write(6,60) 'Exx |',Dmin,'|',Dmax,'|',Dmoy,'|   -'
      ExxT = Dmoy
      ExxVAR = Dvar
      Call MinMaxMoy2(nt,Eyy,Dmin,Dmax,Dmoy,Dvar,SS)
      Write(6,60) 'Eyy |',Dmin,'|',Dmax,'|',Dmoy,'|   -'
      EyyT = Dmoy
      EyyVAR = Dvar
      Call MinMaxMoy(nt,E1,Dmin,Dmax,Dmoy,SS)
      Write(6,60) 'E1  |',Dmin,'|',Dmax,'|',Dmoy,'|   -'

      Call MinMaxMoy(nt,E2,Dmin,Dmax,Dmoy,SS)
      Write(6,60) 'E2  |',Dmin,'|',Dmax,'|',Dmoy,'|   -'

c      Call moment(Dis,nt,ave,adev,sdev,var,skew,curt,SS)
      Call MinMaxMoy2(nt,Dis,Dmin,Dmax,Dmoy,Dvar,SS)
c      Call IndicLoca(nt,Dis,SS,S2)
      Write(6,60) 'Dis |',Dmin,'|',Dmax,'|',Dmoy,'|   -'
      DisT = Dmoy
      DisVAR = Dvar

      Call MinMaxMoy(nt,Vol,Dmin,Dmax,Dmoy,SS)
      Write(6,60) 'Vol |',Dmin,'|',Dmax,'|',Dmoy,'|   -'


c      Open(15,file='fort.15',access='APPEND')
c      Write(15,'(10(e15.7,1x))') Eyym,Eyym-Exxm,ave,adev,
c     &                          sdev,var,skew,curt,DisT,S2
c      Write(15,'(3(e15.7,1x))') Eyym,DisT,DisVAR
c      Write(15,'(12(e15.7,1x))') Exym,ExyT,Ewm,EwT,Exxm,ExxT,Eyym,EyyT,
c     &                           ExyVAR,EwVAR,ExxVAR,EyyVAR


      Write(6,*) '________________________________________________'
60    Format(1x,a5,e10.3,a1,e10.3,a1,e10.3,a5)
61    Format(1x,a5,e10.3,a1,e10.3,a1,e10.3,a1,e10.3)

      Write(6,*) ' 1=Exx, 2=Eyy, 3=E1, 4=E2, 5=Dis, 6=Vol, 0=SORTIE'
      Write(6,*) ' 7=Exy, 8=Eww, 9 = echantillon+triangle'
      Read(5,*) irep1
      If ((irep1.LT.0).OR.(irep1.GT.9)) GoTo 1
    
      If (irep1.EQ.0) Then
         CALL PGCLOS
         STOP
      EndIf
      If (irep1.EQ.1) Then
         Do i = 1,nt
            Dd(i) = Exx(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.2) Then
         Do i = 1,nt
            Dd(i) = Eyy(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.3) Then
         Do i = 1,nt
            Dd(i) = E1(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.4) Then
         Do i = 1,nt
            Dd(i) = E2(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.5) Then
         Do i = 1,nt
            Dd(i) = Dis(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.6) Then
         Do i = 1,nt
            Dd(i) = Vol(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.7) Then
         Do i = 1,nt
            Dd(i) = Exy(i)*100.
         EndDo
      EndIf
      If (irep1.EQ.8) Then
         Do i = 1,nt
            Dd(i) = Ew(i)
         EndDo
      EndIf
      If (irep1.EQ.9) Then
         GoTo 3
      EndIf

c      Write(6,*) '---------------------------------------'

      Call MinMaxMoy2(nt,Dd,Dmin,Dmax,Dmoy,Dvar,SS)
      Write(6,*) 'Dmin,Dmax = ',Dmin,Dmax
      Write(6,*) 'Dmoy,Dvar = ',Dmoy,Dvar

c=====Le dessin des def avec PGPLOT

      CALL PGQCIR(iC1,iC2)
      NC = MAX(0,iC2-iC1+1)
c      WRITE(6,*) 'Nbre de couleur utilisees: ',NC
      BRIGHT = 0.5
      CONTRA  = 1.0
2     Write(6,*) '- - - - - - - - - - - - - - - - - - - -'
      Write(6,*) 'Type de couleurs ?'
      Write(6,*) ' 1 - niveaux de gris'
      Write(6,*) ' 2 - arc en ciel'
      Write(6,*) ' 3 - dominante jaune-rouge'
      Write(6,*) ' 4 - couleurs distincts (degradees)'
      Write(6,*) ' 5 - couleurs distincts'
      Write(6,*) ' 6 - symboles monochromes'
      Read(5,*) irep2
      If ((irep2.LT.1).OR.(irep2.GT.6)) GoTo 2
      If (irep2.EQ.6) GoTo 4 
      Write(6,*) 'luminosite = ',BRIGHT
      Write(6,*) 'Contraste  = ',CONTRA
      CALL PGERAS
      CALL PALETT(irep2, CONTRA, BRIGHT)     
      Do i = 1,nt
         a = (ic2-ic1) / (Dmax - Dmin)
         b = (ic1) - a * Dmin
         icoul = a*(Dd(i)) + b
         If (icoul.LT.ic1) icoul = ic1
         If (icoul.GT.ic2) icoul = ic2
         i1 = lt(1,i)
         i2 = lt(2,i)
         i3 = lt(3,i)
         xp(1) = x2(1,i1)
         yp(1) = x2(2,i1)
         xp(2) = x2(1,i2)
         yp(2) = x2(2,i2)
         xp(3) = x2(1,i3)
         yp(3) = x2(2,i3)
         CALL PGSFS(1)
         CALL PGSCI(icoul)
         CALL PGPOLY(3,xp,yp)
      EndDo
      CALL PGSCI(1)
      fDmin = Dmin
      fDmax = Dmax
      If (irep1.EQ.1) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,'Exx (%)')
      If (irep1.EQ.2) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,'Eyy (%)')
      If (irep1.EQ.3) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Deformation principale majeure (%)')
      If (irep1.EQ.4) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Deformation principale mineure (%)')
      If (irep1.EQ.5) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Intensite de cisaillement (%)')
      If (irep1.EQ.6) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Deformation volumique (%)')
      If (irep1.EQ.7) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Exy (%)')
      If (irep1.EQ.8) CALL PGWEDG('BI',0.,3.0,fDmin,fDmax,
     &                            'Rotation (rd)')

      GoTo 1

c-----dessin avec des symboles
4     Continue
      Write(6,*) 'symboles...'
      CALL PGERAS
      If (Dmin.LT.0) Then
c        --Il y a des valeurs <0
         aa = Max(abs(Dmin),abs(Dmax))
         binf = 3.*Dmin/aa
         bsup = 3.*Dmax/aa
         a = (bsup-binf) / (Dmax - Dmin)
         b = (binf) - a * Dmin
      Else
c        --toutes les valeurs sont >0
         binf = 0.001
         bsup = 3.
         a = (bsup-binf) / (Dmax - Dmin) 
         b = (binf) - a * Dmin
      EndIf
      Do i = 1,nt
         i1 = lt(1,i)
         i2 = lt(2,i)
         i3 = lt(3,i)
         xp(1) = x2(1,i1)
         yp(1) = x2(2,i1)
         xp(2) = x2(1,i2)
         yp(2) = x2(2,i2)
         xp(3) = x2(1,i3)
         yp(3) = x2(2,i3)
         CALL PGSFS(2)
         CALL PGSCI(1)
         size = a*Dd(i) + b
         size = abs(size)
         CALL PGSCH(size)
         Write(77,*) Dd(i),size
         xx = (xp(1)+xp(2)+xp(3))/3.
         yy = (yp(1)+yp(2)+yp(3))/3.
         If (Dd(i).GE.0.) Then
            CALL PGPT1(xx,yy,19)
         Else
            CALL PGPT1(xx,yy,24)
         EndIf
      EndDo

      GoTo 1

c-----Dessin de l'echantillon et de la triangulation
3     Continue
c      CALL PGERAS
      Pi = 4.d0*atan(1.d0)

c-----les triangles
      Do i = 1,nt
         i1 = lt(1,i)
         i2 = lt(2,i)
         i3 = lt(3,i)
         xp(1) = x2(1,i1)
         yp(1) = x2(2,i1)
         xp(2) = x2(1,i2)
         yp(2) = x2(2,i2)
         xp(3) = x2(1,i3)
         yp(3) = x2(2,i3)
         CALL PGSFS(2)
         CALL PGSCI(2)
         CALL PGPOLY(3,xp,yp)
      EndDo


      GoTo 1

      End 

c=======================================================================
      Subroutine MinMaxMoy(Nbre,Tab,rmin,rmax,rmoy,SurfT)

      Implicit Real*8 (a-h,o-z)
      Integer Nbre
      Dimension Tab(Nbre),SurfT(Nbre)

      rmin =  1e20
      rmax = -1e20
      rmoy = 0.d0
      stotal = 0.d0
      Do i = 1,Nbre
         rmin = min(rmin,Tab(i))
         rmax = max(rmax,Tab(i))
         rmoy = rmoy + Tab(i)*Surft(i)
         stotal = stotal + Surft(i)
      EndDo
      rmoy = rmoy / stotal

      Return
      End

c=======================================================================
      Subroutine MinMaxMoy2(Nbre,Tab,rmin,rmax,rmoy,rvar,SurfT)

      Implicit Real*8 (a-h,o-z)
      Integer Nbre
      Dimension Tab(Nbre),SurfT(Nbre)

      rmin =  1e20
      rmax = -1e20
      rmoy = 0.
      stotal = 0.
      Do i = 1,Nbre
         rmin = min(rmin,Tab(i))
         rmax = max(rmax,Tab(i))
         rmoy = rmoy + Tab(i)*Surft(i)
         stotal = stotal + Surft(i)
      EndDo
      rmoy = rmoy / stotal
      rvar = 0.d0
      Do i = 1,Nbre
         rvar = rvar + (rmoy-Tab(i))
     &                *(rmoy-Tab(i))
      EndDo
      rvar = sqrt(rvar/dfloat(Nbre))

      Return
      End

c=======================================================================
      Subroutine IndicLoca(Nbre,Tab,Surft,S2)

      Implicit Real*8 (a-h,o-z)
      Integer Nbre
      Dimension Tab(Nbre),SurfT(Nbre)

      stotal = 0.
      Do i = 1,Nbre
         stotal = stotal + Surft(i)
      EndDo

      Do i = 1,Nbre
         rNum = rNum + Tab(i)
         rVar = rVar + Tab(i)*Tab(i)
      EndDo
      rNum2 = rNum*rNum
      S2 = (rNum2/rVar)/float(Nbre)

      Write(77,*) S2

      Return
      End

c=======================================================================
      Real*8 Function Surf1(x,i1,i2,i3)
  
      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Integer i1,i2,i3
      Dimension x(2,npam)

      a1 = x(1,i2)-x(1,i1)
      b1 = x(2,i2)-x(2,i1)
      a2 = x(1,i3)-x(1,i1)
      b2 = x(2,i3)-x(2,i1)
      Surf = (a1*b2 - b1*a2)*0.5
      If (Surf.LT.0.d0) Surf = -Surf
   
      Return
      End
c=======================================================================
      Real*8 Function Surf(x,i1,i2,i3)
  
      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Integer i1,i2,i3
      Real*8 x(2,npam)
 
      ax = x(1,i2)-x(1,i1)
      ay = x(2,i2)-x(2,i1)
      a = sqrt(ax*ax+ay*ay)
      bx = x(1,i3)-x(1,i2)
      by = x(2,i3)-x(2,i2)
      b = sqrt(bx*bx+by*by)
      cx = x(1,i3)-x(1,i1)
      cy = x(2,i3)-x(2,i1)
      c = sqrt(cx*cx+cy*cy)
      p = a+b+c
      Surf = sqrt(p*(p-2*a)*(p-2*b)*(p-2*c)/16.)

      Return
      End
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
      Subroutine RetireTrianglesPlats(npa,r1,nt,lt)

      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension lt(3,npam),lt0(3,npam)
      Dimension r1(2,npam),x(2,npam)

      Do i = 1,npa
         x(1,i) = r1(1,i)
         x(2,i) = r1(2,i)
      EndDo
      Do i = 1,nt
         lt0(1,i) = lt(1,i)
         lt0(2,i) = lt(2,i)
         lt0(3,i) = lt(3,i)
      EndDo

      pi = 4.*atan(1.)

      ntt = 0
      Do i = 1,nt
         i0 = lt0(1,i)
         i1 = lt0(2,i)
         i2 = lt0(3,i)
         u01x = x(1,i1)-x(1,i0)
         u01y = x(2,i1)-x(2,i0)
         u12x = x(1,i2)-x(1,i1)
         u12y = x(2,i2)-x(2,i1)
         u20x = x(1,i2)-x(1,i0)
         u20y = x(2,i2)-x(2,i0)
         s0112 = u01x*u12x + u01y*u12y
         s1220 = u12x*u20x + u12y*u20y
         s2001 = u20x*u01x + u20y*u01y
         d01 = sqrt((x(1,i0)-x(1,i1))**2. + (x(2,i0)-x(2,i1))**2.)
         d12 = sqrt((x(1,i1)-x(1,i2))**2. + (x(2,i1)-x(2,i2))**2.)
         d20 = sqrt((x(1,i2)-x(1,i0))**2. + (x(2,i2)-x(2,i0))**2.)
         cos0112 = s0112/d01/d12
         cos1220 = s1220/d12/d20
         cos2001 = s2001/d20/d01
         a0112 = acos(cos0112)*180./pi
         a1220 = acos(cos1220)*180./pi
         a2001 = acos(cos2001)*180./pi 
         If (cos0112.LT.0) a0112 = a0112 - 180.
         If (cos1220.LT.0) a1220 = a1220 - 180.
         If (cos2001.LT.0) a2001 = a2001 - 180.

         S = Surf(x,i0,i1,i2)
         Srec = 0.!0.05*0.02/30.
 
         tol = 30.
         angletol = 10. !degres
         If ((max(d01,d12,d20).GT.tol).OR.
     &       (abs(a0112).LT.angletol).OR.
     &        (abs(a1220).LT.angletol).OR.
     &        (abs(a2001).LT.angletol).OR.
     &        (S.LT.Srec)
     &      ) Then
            Continue
         Else
            ntt = ntt + 1
            lt(1,ntt) = i0
            lt(2,ntt) = i1
            lt(3,ntt) = i2
         EndIf
      EndDo

      nt = ntt
   
      Write(6,*) 'Retrait des triangles plats :'
      Write(6,*) ' Nbre total de triangles apres elimination: ',nt

      Return
      End
c=======================================================================
      Subroutine RetireTrianglesRectangle(nt,lt)

      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension lt(3,npam),lt0(3,npam)

      Do i = 1,nt
         lt0(1,i) = lt(1,i)
         lt0(2,i) = lt(2,i)
         lt0(3,i) = lt(3,i)
      EndDo
      
      ntt = 0
      Do i = 1,nt
         i0 = lt0(1,i)
         i1 = lt0(2,i)
         i2 = lt0(3,i)
         k1 = (i0-1)/15
         k2 = (i1-1)/15
         k3 = (i2-1)/15
         If (k1.EQ.k2.AND.k1.EQ.k3) Then
            Continue
         Else
            ntt = ntt + 1
            lt(1,ntt) = i0
            lt(2,ntt) = i1
            lt(3,ntt) = i2
         EndIf
      EndDo       
      nt = ntt
      Write(6,*) 'Retrait des triangles des rectangles :'
      Write(6,*) ' Nbre total de triangles apres elimination: ',nt

      Return
      End
c=======================================================================
      Subroutine Triangule(npa,r1,ivnxt,lt)

      Implicit Real*8 (a-h,o-z)
      Parameter(npam=30000)
      Dimension lt(3,npam)
      Dimension r1(2,npam)
      Parameter (nmax=60000,nsqrt=250,nvmax=2.1*nmax,nicon=6)
      Parameter (nfmax=1.5*nsqrt,nhmax=2.*nsqrt)
      Parameter (nkmax=130,njmax=2*nkmax,npmax=2000)
      Real*8 x(nmax),y(nmax)
      Integer is(nmax),icon(nicon,nvmax),id(nvmax)
      Integer iave(njmax),iaze(nkmax)

c-----Initialisations
      nvmz = -nvmax - 1
      nvmn = nvmz - 1
      nvmp = -nvmn
      eps = 0.000001
      epz = 0.00003
      epw = 0.000003
      epv = 1.0

      nv = 0
      Do i = 1,npa
         nv = nv + 1
         xcor = r1(1,i)
         ycor = r1(2,i)
         x(nv) = xcor
         y(nv) = ycor
         If (nv.NE.1) Then
            If (xcor.GT.xmax) xmax = xcor
            If (xcor.LT.xmin) xmin = xcor
            If (ycor.GT.ymax) ymax = ycor
            If (ycor.LT.ymin) ymin = ycor
         Else
            xmax = xcor
            xmin = xcor
            ymax = ycor
            ymin = ycor
         EndIf
      EndDo

      If (nv.LT.3) Then
         Write(6,*) 'n < 3'
         stop
      EndIf
      Write(6,*) 'Nbre de points lus',nv
      n = nv
      iase = 1
      iaze(1) = n+1

      Write(6,*) 'xmin = ',xmin,' xmax = ',xmax
      Write(6,*) 'ymin = ',ymin,' ymax = ',ymax

c-----Triangulation
      Write(6,*) 'Calcul Triangulation...'
      Do i = 1,nv
         is(i) = 0
      EndDo
c     --Triangle initial
      xlen = (xmax-xmin)/4.0
      ylen = (ymax-ymin)/4.0
      If (xlen.LT.epv) epv = xlen
      If (ylen.LT.epv) epv = ylen

      i1 = 1
      Do i2 = i1+1, nv-1
          If (is(i2).EQ.0) Then
             xtemp = x(i2) - x(i1)
             ytemp = y(i2) - y(i1)
             If ( (xtemp.LE.-epv).OR.
     &            (xtemp.GE. epv).OR.
     &            (ytemp.LE.-epv).OR. 
     &            (ytemp.GE.epv)
     &          ) Then
                GoTo 320
             EndIf
          EndIf
      EndDo

  320 Continue
      If ((xtemp.GT.-eps).AND.(xtemp.LT.eps)) xtemp = 0.0
      If ((ytemp.GT.-eps).AND.(ytemp.LT.eps)) ytemp = 0.0
      def = sqrt(xtemp*xtemp + ytemp*ytemp)
      Do i3 = i2+1, nv
          If (is(i3).NE.0) GoTo 350
          delx = x(i3) - x(i2)
          dely = y(i3) - y(i2)
          If ( (delx.GT.-epv).AND.
     &         (delx.LT. epv).AND.
     &         (dely.GT.-epv).AND.
     &         (dely.LT. epv)
     &       ) Then
             GoTo 350
          EndIf
          If ((delx.GT.-eps).AND.(delx.LT.eps)) delx = 0.0
          If ((dely.GT.-eps).AND.(dely.LT.eps)) dely = 0.0
          daf = sqrt(delx**2 + dely**2)
          delx = y(i1) - y(i3)
          dely = x(i3) - x(i1)
          If ((delx.GT.-epv).AND.(delx.LT.epv).AND.
     *       (dely.GT.-epv).AND.(dely.LT.epv)) go to 350
          If ((delx.GT.-eps).AND.(delx.LT.eps)) delx = 0.0
          If ((dely.GT.-eps).AND.(dely.LT.eps)) dely = 0.0
          dif = sqrt(delx**2 + dely**2)
          If (dif.LT.epv) go to 350
          dmax = max(dif,daf,def)
          dot = (delx*xtemp + dely*ytemp)/dmax
          If (dot.LE.-epv) go to 400
          If (dot.GE.epv) go to 380
  350     Continue
      EndDo
      Stop

  380 Continue
      it = i2
      i2 = i3
      i3 = it

  400 Continue
      ivnxt = 1
      If (ivnxt.GT.nvmax) Then
          Write(6,*) 'nvmax atteint'
          Stop
      endif

      icon(1,1) = 0
      icon(2,1) = 0
      icon(3,1) = 0
      icon(4,1) = i1
      icon(5,1) = i2
      icon(6,1) = i3

      is(i1) = 1
      is(i2) = 1
      is(i3) = 1

c     --insert other points
      iabe=1
      iave(1)=ivnxt+1
      ileft=i1
      insrt=1
      if(n.eq.3) go to 480
      idupl=0
c
      do 450 isite = 2, n
         if(is(isite).ne.0) go to 450
         xdes=x(isite)
         ydes=y(isite)
         call pntins(x, y, is, icon, iave, iaze, eps, epz, epw, nicon,
     *            nvmz, ileft, iright, n, xdes, ydes, iabe, iase, isite,
     *            nmax, nvmax, ivnxt, insrt, itype)
         if(itype.eq.1.or.itype.eq.-1)idupl=idupl+1
c
         ileft = iright
c         if(isite.le.(isite/500)*500)
c     *      Write(6,*)'number of points inserted=',isite
  450 continue
      If (idupl.NE.0) Then
      Write(6,*)'number of points read that were duplications=',idupl
      EndIf
c
  480 continue

c-----Sauvegarde des num des grains des triangles
c      Write(20,890)((icon(i,j),i=4,6),j=1,ivnxt)
c  890 format(3i8)

c-----tests
      call contst(icon, is, id, nv, ivnxt, nicon, nvmz)
      call circri (icon, ivnxt, nicon, x, y, epz, eps)
c      Write(6,*)'probable number of vertices in triangulation: ',n
      Write(6,*)'Nbre total de triangle: ',ivnxt
      Write(6,*)''

      Do j = 1,ivnxt
         lt(1,j) = icon(4,j)
         lt(2,j) = icon(5,j)
         lt(3,j) = icon(6,j)
      EndDo

      Return
      End
c=======================================================================
c
c     subroutine pntins to -
c
c     insert a point into a Delaunay triangulation and obtain a
c     Delaunay triangulation that contains the inserted point
c
c     February 11, 1991
c
      subroutine pntins(x, y, is, icon, iave, iaze, eps, epz, epw,
     *                  nicon, nvmz, ileft, iright, n, xdes, ydes, iabe,
     *                  iase, isite, nmax, nvmax, ivnxt, insrt, itype)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer is(*), icon(nicon,*), iave(*), iaze(*)
c
c     find position of point in triangulation
c
      call gettri(x, y, is, icon, eps, epz, epw, nicon, nvmz, ileft,
     *            iright, n, itype, ivfnd, xdes, ydes, isite)
c
c     insert point
c
      if(itype .eq. -1 .or. itype .eq. 1) go to 400
      if(insrt .eq. 0) go to 400
      if(isite.ne.0) go to 200
      isite = iaze(iase)
      if(iase .eq. 1) then
          n = isite
          if(n. gt. nmax) then
              Write(6,*)' pntins: parameter nmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iaze(1) = n + 1
      else
          iase = iase - 1
      endif
      x(isite) = xdes
      y(isite) = ydes
  200 continue
      if(itype .eq. 2) then
          call edgint(icon, is, iave, nicon, ivnxt, nvmax, nvmz,
     *                ivfnd, isite, iright, iabe)
      elseif(itype .eq. 3) then
          call intins(icon, is, iave, nicon, ivnxt, nvmax, nvmz,
     *                ivfnd, isite, iabe)
      else
          call outhul(x, y, icon, is, iave, nicon, ivnxt, nvmax,
     *                nvmz, ivfnd, isite, iright, iabe, eps, epz)
      endif
      iright = isite
c
c     optimize triangulation with respect to isite
c
      ivcur = ivfnd
  300 continue
      ivaft = icon(2,ivcur)
      call opttri(x, y, icon, is, eps, epz, nicon, ivcur, isite, nvmz)
      if(ivaft .eq. 0 .or. ivaft .eq. nvmz) go to 400
      ivcur = iabs(ivaft)
      if(ivcur .ne. ivfnd) go to 300
c
  400 continue
c
      return
      end
c=======================================================================
c
c     subroutine gettri to -
c
c     obtain position in a triangulation for
c     a point (xdes,ydes) by moving through the triangulation from
c     a known vertex (ileft) in the triangulation to the point
c
c     the following holds on output: 
c
c     if itype equals -1 then (xdes,ydes) is exactly the vertex with
c     index ileft
c
c     if itype equals 0 then (xdes,ydes) is a point outside the convex
c     hull of the triangulation and iright is a vertex of the convex
c     hull visible from (xdes,ydes)
c
c     if itype equals 1 then (xdes,ydes) is exactly the vertex in the
c     triangulation with index iright and iright is different from
c     ileft
c
c     if itype equals 2 then (xdes,ydes) is a point in the interior of
c     a side of the triangle in the triangulation with index ivfnd and
c     iright is the first vertex of this triangle which is found when
c     moving from the point in a counterclockwise direction around the
c     boundary of the triangle
c
c     if itype equals 3 then (xdes,ydes) is a point in the interior of
c     the triangle in the triangulation with index ivfnd
c
c     February 6, 1991
c
      subroutine gettri(x, y, is, icon, eps, epz, epw, nicon, nvmz,
     *                ileft, iright, n, itype, ivfnd, xdes, ydes, isite)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer is(*), icon(nicon,*)
c
c     compute direction parameters
c
      if(ileft .lt. 1 .or. ileft .gt. n) then
          Write(6,*)' gettri: vertex out of range'
          Write(6,*)' program terminated'
          stop
      endif
      delx = y(ileft) - ydes
      dely = xdes - x(ileft)
      if(delx .gt. -eps .and. delx .lt. eps) delx = 0.0
      if(dely .gt. -eps .and. dely .lt. eps) dely = 0.0
      dif = sqrt(delx**2 + dely**2)
      itype = -1
      iright = ileft
      if(delx .gt. -epw .and. delx .lt. epw .and.
     *   dely .gt. -epw .and. dely .lt. epw) go to 2000
      if(dif .lt. epz) go to 2000
      isi2 = ileft
      itype = 0
c
  100 continue
      iloft = isi2
      ivini = is(iloft)
      ivcur = ivini
c
c     test current triangle
c
  200 continue
      if(icon(4,ivcur) .eq. iloft) then
          ivadj = icon(2,ivcur)
          ivfol = icon(1,ivcur)
          isi2 = icon(5,ivcur)
          isi1 = icon(6,ivcur)
      elseif(icon(5,ivcur) .eq. iloft) then
          ivadj = icon(3,ivcur)
          ivfol = icon(2,ivcur)
          isi2 = icon(6,ivcur)
          isi1 = icon(4,ivcur)
      elseif(icon(6,ivcur) .eq. iloft) then
          ivadj = icon(1,ivcur)
          ivfol = icon(3,ivcur)
          isi2 = icon(4,ivcur)
          isi1 = icon(5,ivcur)
      else
          Write(6,*)' gettri 1: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          WRITE(6,*)'ILEFT=',ILEFT,' ILOFT=',ILOFT,' IRIGHT=',IRIGHT
          stop
      endif
c
      xtemp = x(isi2) - xdes
      ytemp = y(isi2) - ydes
      if(xtemp .le. -epw .or. xtemp .ge. epw .or.
     *   ytemp .le. -epw .or. ytemp .ge. epw) go to 300
      itype = 1
      iright = isi2
      go to 2000
  300 continue
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dls2 = sqrt(xtemp**2 + ytemp**2)
      xtamp = x(isi2) - x(ileft)
      ytamp = y(isi2) - y(ileft)
      if(xtamp .gt. -eps .and. xtamp .lt. eps) xtamp = 0.0
      if(ytamp .gt. -eps .and. ytamp .lt. eps) ytamp = 0.0
      dlg2 = sqrt(xtamp**2 + ytamp**2)
      dmax2 = max(dif,dls2,dlg2)
      dot1 = (delx * xtamp + dely * ytamp)/dmax2
      if(dot1 .ge. epz) go to 500
      if(dot1 .le.-epz) go to 320
      xtimp = x(isi2) - x(iloft)
      ytimp = y(isi2) - y(iloft)
      if(xtimp .gt. -eps .and. xtimp .lt. eps) xtimp = 0.0
      if(ytimp .gt. -eps .and. ytimp .lt. eps) ytimp = 0.0
      duf = sqrt(xtimp**2 + ytimp**2)
      dmex = max(dif,duf)
      dut1 = (dely * xtimp - delx * ytimp)/dmex
      if(dut1 .gt. 0.0) go to 380
      dot2 = -1.0
      go to 500
  320 continue
      xtump = x(isi1) - xdes
      ytump = y(isi1) - ydes
      if(xtump .le. -epw .or. xtump .ge. epw .or.
     *   ytump .le. -epw .or. ytump .ge. epw) go to 350
      itype = 1
      iright = isi1
      go to 2000
  350 continue
      if(xtump .gt. -eps .and. xtump .lt. eps) xtump = 0.0
      if(ytump .gt. -eps .and. ytump .lt. eps) ytump = 0.0
      dls1 = sqrt(xtump**2 + ytump**2)
      xtamp = x(isi1) - x(ileft)
      ytamp = y(isi1) - y(ileft)
      if(xtamp .gt. -eps .and. xtamp .lt. eps) xtamp = 0.0
      if(ytamp .gt. -eps .and. ytamp .lt. eps) ytamp = 0.0
      dlg1 = sqrt(xtamp**2 + ytamp**2)
      dmax1 = max(dif,dls1,dlg1)
      dot2 = (delx * xtamp + dely * ytamp)/dmax1
      if(dot2 .lt. epz) go to 500
      if(dot1 .le. -epz) go to 900
  380 continue
      dot3 = (-dely * xtemp + delx * ytemp)/dmax2
      if(dot3 .ge. epz) go to 100
      if(dot3 .le. -epz) go to 400
      itype = 1
      iright = isi2
      go to 2000
  400 continue
      itype = 2
      iright = isi2
      ivfnd = ivcur
      go to 2000
c
  500 continue
      ivpre = ivcur
      ivcur = iabs(ivadj)
      if(ivcur .ne. ivini) go to 600
      Write(6,*)' gettri: can not direct segment'
      Write(6,*)' program terminated'
      stop
  600 continue
      if(ivadj .ne. 0 .and. ivadj .ne. nvmz) go to 200
      iright = iloft
      if(dot1 .ge. epz) go to 2000
      if(dot2 .le. -epz) go to 2000
      dot3 = (-dely * xtump + delx * ytump)/dmax1
      isi2 = isi1
      if(dot3 .ge. epz) go to 100
      if(dot3 .le. -epz) go to 800
      itype = 1
      iright = isi1
      go to 2000
  800 continue
      itype = 2
      ivfnd = ivpre
      go to 2000
c
  900 continue
      dalx = y(isi2) - y(isi1)
      daly = x(isi1) - x(isi2)
      if(dalx .gt. -eps .and. dalx .lt. eps) dalx = 0.0
      if(daly .gt. -eps .and. daly .lt. eps) daly = 0.0
      daf = sqrt(dalx**2 + daly**2)
      if(daf .lt. eps) then
          Write(6,*)' gettri 1: distinct sites appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      dmax3 = max(daf,dls1,dls2)
      dot3 = (dalx * xtemp + daly * ytemp)/dmax3
      if(dot3 .ge. epz) go to 1100
      if(dot3 .le. -epz) go to 1000
      itype = 2
      iright = isi1
      ivfnd = ivcur
      go to 2000
 1000 continue
      itype = 3
      ivfnd = ivcur
      go to 2000
c
c     test next triangle
c
 1100 continue
      iright = isi2
      if(ivfol .eq. 0 .or. ivfol .eq. nvmz) go to 2000
      ivcur = iabs(ivfol)
      if(icon(4,ivcur) .eq. isi2) then
          isi3 = icon(5,ivcur)
          ivadj = icon(1,ivcur)
          ivpre = icon(3,ivcur)
      elseif(icon(5,ivcur) .eq. isi2) then
          isi3 = icon(6,ivcur)
          ivadj = icon(2,ivcur)
          ivpre = icon(1,ivcur)
      elseif(icon(6,ivcur) .eq. isi2) then
          isi3 = icon(4,ivcur)
          ivadj = icon(3,ivcur)
          ivpre = icon(2,ivcur)
      else
          Write(6,*)' gettri 2: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      xtemp = x(isi3) - xdes
      ytemp = y(isi3) - ydes
      if(xtemp .le. -epw .or. xtemp .ge. epw .or.
     *   ytemp .le. -epw .or. ytemp .ge. epw) go to 1200
      itype = 1
      iright = isi3
      go to 2000
 1200 continue
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dls = sqrt(xtemp**2 + ytemp**2)
      xtamp = x(isi3) - x(ileft)
      ytamp = y(isi3) - y(ileft)
      if(xtamp .gt. -eps .and. xtamp .lt. eps) xtamp = 0.0
      if(ytamp .gt. -eps .and. ytamp .lt. eps) ytamp = 0.0
      dlg = sqrt(xtamp**2 + ytamp**2)
      dmax = max(dif,dls,dlg)
      dot = (delx * xtamp + dely * ytamp)/dmax
      if(dot .ge. epz) go to 1400
      if(dot .le. -epz) go to 1700
      dot = (-dely * xtemp + delx * ytemp)/dmax
      isi2 = isi3
      if(dot .ge. epz) go to 100
      if(dot .le. -epz) go to 1300
      itype = 1
      iright = isi3
      go to 2000
 1300 continue
      itype = 3
      ivfnd = ivcur
      go to 2000
c
 1400 continue
      dls1 = dls
      dalx = y(isi3) - y(isi2)
      daly = x(isi2) - x(isi3)
      if(dalx .gt. -eps .and. dalx .lt. eps) dalx = 0.0
      if(daly .gt. -eps .and. daly .lt. eps) daly = 0.0
      daf = sqrt(dalx**2 + daly**2)
      if(daf .lt. eps) then
          Write(6,*)' gettri 2: distinct sites appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      dmax3 = max(daf,dls1,dls2)
      dot = (dalx * xtemp + daly * ytemp)/dmax3
      if(dot .le. -epz) go to 1600
      if(dot .ge. epz) go to 1500
      itype = 2
      iright = isi3
      ivfnd = ivcur
      go to 2000
 1500 continue
      itype = 3
      ivfnd = ivcur
      go to 2000
 1600 continue
      isi1 = isi3
      ivfol = ivpre
      go to 1100
c
 1700 continue
      dls2 = dls
      dalx = y(isi3) - y(isi1)
      daly = x(isi1) - x(isi3)
      if(dalx .gt. -eps .and. dalx .lt. eps) dalx = 0.0
      if(daly .gt. -eps .and. daly .lt. eps) daly = 0.0
      daf = sqrt(dalx**2 + daly**2)
      if(daf .lt. eps) then
          Write(6,*)' gettri 3: distinct sites appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      dmax3 = max(daf,dls1,dls2)
      dot = (dalx * xtemp + daly * ytemp)/dmax3
      if(dot .ge. epz) go to 1900
      if(dot .le. -epz) go to 1800
      itype = 2
      iright = isi1
      ivfnd = ivcur
      go to 2000
 1800 continue
      itype = 3
      ivfnd = ivcur
      go to 2000
 1900 continue
      isi2 = isi3
      ivfol = ivadj
      go to 1100
c
 2000 continue
      return
      end
c=======================================================================
c
c     subroutine edgint to -
c
c     insert a point into current triangulation when point is in the
c     interior of an edge of a triangle
c
c     January 31, 1991
c
      subroutine edgint(icon, is, iave, nicon, ivnxt, nvmax, nvmz,
     *                  ivcur, isite, islt, iabe)
c
      Implicit Real*8 (a-h,o-z)
      integer icon(nicon,*), is(*), iave(*)
c
      it1 = ivcur
      if(icon(4,ivcur) .eq. islt) then
          isrt = icon(6,ivcur)
          isad = icon(5,ivcur)
          ivadj = icon(2,ivcur)
          iv1 = icon(3,ivcur)
          iv2 = icon(1,ivcur)
      elseif(icon(5,ivcur) .eq. islt) then
          isrt = icon(4,ivcur)
          isad = icon(6,ivcur)
          ivadj = icon(3,ivcur)
          iv1 = icon(1,ivcur)
          iv2 = icon(2,ivcur)
      elseif(icon(6,ivcur) .eq. islt) then
          isrt = icon(5,ivcur)
          isad = icon(4,ivcur)
          ivadj = icon(1,ivcur)
          iv1 = icon(2,ivcur)
          iv2 = icon(3,ivcur)
      else
          Write(6,*)' edgint 1: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
      it3 = iabs(ivadj)
      it4 = it3
      if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 100
      if(icon(4,it4) .eq. islt) then
          isop = icon(6,it4)
          iv3 = icon(1,it4)
          iv4 = icon(2,it4)
      elseif(icon(5,it4) .eq. islt) then
          isop = icon(4,it4)
          iv3 = icon(2,it4)
          iv4 = icon(3,it4)
      elseif(icon(6,it4) .eq. islt) then
          isop = icon(5,it4)
          iv3 = icon(3,it4)
          iv4 = icon(1,it4)
      else
          Write(6,*)' edgint 2: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
  100 continue
      it2 = iave(iabe)
      if(iabe .eq. 1) then
          ivnxt = it2
          if(ivnxt .gt. nvmax) then
              Write(6,*)' edgint 1: parameter nvmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(1) = ivnxt + 1
      else
          iabe = iabe - 1
      endif
      if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 200
      it3 = iave(iabe)
      if(iabe .eq. 1) then
          ivnxt = it3
          if(ivnxt .gt. nvmax) then
              Write(6,*)' edgint 2: parameter nvmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(1) = ivnxt + 1
      else
          iabe = iabe - 1
      endif
c
  200 continue
      is(isite) = it1
      if(is(isad) .eq. it1) is(isad) = it2
      if(is(isrt) .eq. it1) is(isrt) = it2
      if(is(isrt) .eq. it4) is(isrt) = it3
c
c     create new triangles
c
      icon(1,it1) = iv1
      icon(2,it1) = it2
      icon(3,it1) = it4
      icon(4,it1) = isite
      icon(5,it1) = islt
      icon(6,it1) = isad
c
      icon(1,it2) = iv2
      icon(2,it2) = it3
      icon(3,it2) = it1
      icon(4,it2) = isite
      icon(5,it2) = isad
      icon(6,it2) = isrt
c
      if(islt.eq.isad.or.isrt.eq.isad.or.isite.eq.islt.or.
     *   isite.eq.isrt.or.isite.eq.isad)then
         Write(6,*)' edgint 1: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' islt,isrt,isad:',islt,isrt,isad
         Write(6,*)' isite=',isite
         stop
      endif
c
      if(ivadj .lt. 0) then
          icon(3,it1) = -it4
          icon(2,it2) = -it3
      endif
c
      if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 300
      icon(1,it3) = iv3
      icon(2,it3) = it4
      icon(3,it3) = it2
      icon(4,it3) = isite
      icon(5,it3) = isrt
      icon(6,it3) = isop
c
      icon(1,it4) = iv4
      icon(2,it4) = it1
      icon(3,it4) = it3
      icon(4,it4) = isite
      icon(5,it4) = isop
      icon(6,it4) = islt
c
      if(islt.eq.isop.or.isrt.eq.isop.or.isite.eq.isop)then
         Write(6,*)' edgint 2: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' islt,isrt,isop:',islt,isrt,isop
         Write(6,*)' isite=',isite
         stop
      endif
c
      if(ivadj .lt. 0) then
          icon(3,it3) = -it2
          icon(2,it4) = -it1
      endif
c
c     update neighboring triangles
c
  300 continue
      if(iv2 .eq. 0 .or. iv2 .eq. nvmz) go to 400
      ivodj = it1
      if(iv2 .lt. 0) then
          iv2 = -iv2
          it2 = -it2
          ivodj = -ivodj
      endif
c
      if(icon(1,iv2) .eq. ivodj) then
          icon(1,iv2) = it2
      elseif(icon(2,iv2) .eq. ivodj) then
          icon(2,iv2) = it2
      elseif(icon(3,iv2) .eq. ivodj) then
          icon(3,iv2) = it2
      else
          Write(6,*)' edgint 1: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  400 continue
c
      if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 500
      if(iv3 .eq. 0 .or. iv3 .eq. nvmz) go to 500
      ivodj = it4
      if(iv3 .lt. 0) then
          iv3 = -iv3
          it3 = -it3
          ivodj = -ivodj
      endif
c
      if(icon(1,iv3) .eq. ivodj) then
          icon(1,iv3) = it3
      elseif(icon(2,iv3) .eq. ivodj) then
          icon(2,iv3) = it3
      elseif(icon(3,iv3) .eq. ivodj) then
          icon(3,iv3) = it3
      else
          Write(6,*)' edgint 2: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  500 continue
c
      return
      end
c=======================================================================
c
c     subroutine intins to -
c
c     insert a point into current triangulation when point is in
c     the interior of a triangle
c
c     January 31, 1991
c
      subroutine intins(icon, is, iave, nicon, ivnxt, nvmax, nvmz,
     *                  ivcur, isite, iabe)
c
      Implicit Real*8 (a-h,o-z)
      integer icon(nicon,*), is(*), iave(*)
c
      it1 = ivcur
      it2 = iave(iabe)
      if(iabe .eq. 1) then
          ivnxt = it2
          if(ivnxt .gt. nvmax) then
              Write(6,*)' intins 1: parameter nvmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(1) = ivnxt + 1
      else
          iabe = iabe - 1
      endif
      it3 = iave(iabe)
      if(iabe .eq. 1) then
          ivnxt = it3
          if(ivnxt .gt. nvmax) then
              Write(6,*)' intins 2: parameter nvmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(1) = ivnxt + 1
      else
          iabe = iabe - 1
      endif
c
      iv1 = icon(1,ivcur)
      iv2 = icon(2,ivcur)
      iv3 = icon(3,ivcur)
c
      is4 = icon(4,ivcur)
      is5 = icon(5,ivcur)
      is6 = icon(6,ivcur)
c
      if(is4.eq.is5.or.is4.eq.is6.or.is5.eq.is6.or.isite.eq.is4.or.
     *   isite.eq.is5.or.isite.eq.is6)then
         Write(6,*)' intins: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' is4,is5,is6:',is4,is5,is6
         Write(6,*)' isite=',isite,' ivcur=',ivcur
         stop
      endif
c
      is(isite) = it1
      if(is(is4) .eq. ivcur) is(is4) = it3
      if(is(is6) .eq. ivcur) is(is6) = it2
c
c     create new triangles
c
      icon(1,it1) = iv1
      icon(2,it1) = it2
      icon(3,it1) = it3
      icon(4,it1) = isite
      icon(5,it1) = is5
      icon(6,it1) = is6
c
      icon(1,it2) = iv2
      icon(2,it2) = it3
      icon(3,it2) = it1
      icon(4,it2) = isite
      icon(5,it2) = is6
      icon(6,it2) = is4
c
      icon(1,it3) = iv3
      icon(2,it3) = it1
      icon(3,it3) = it2
      icon(4,it3) = isite
      icon(5,it3) = is4
      icon(6,it3) = is5
c
c     update neighboring triangles
c
      if(iv2 .eq. 0 .or. iv2 .eq. nvmz) go to 100
      ivadj = ivcur
      if(iv2 .lt. 0) then
          iv2 = -iv2
          it2 = -it2
          ivadj = -ivadj
      endif
c
      if(icon(1,iv2) .eq. ivadj) then
          icon(1,iv2) = it2
      elseif(icon(2,iv2) .eq. ivadj) then
          icon(2,iv2) = it2
      elseif(icon(3,iv2) .eq. ivadj) then
          icon(3,iv2) = it2
      else
          Write(6,*)' intins 1: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  100 continue
c
      if(iv3 .eq. 0 .or. iv3 .eq. nvmz) go to 200
      ivadj = ivcur
      if(iv3 .lt. 0) then
          iv3 = -iv3
          it3 = -it3
          ivadj = -ivadj
      endif
c
      if(icon(1,iv3) .eq. ivadj) then
          icon(1,iv3) = it3
      elseif(icon(2,iv3) .eq. ivadj) then
          icon(2,iv3) = it3
      elseif(icon(3,iv3) .eq. ivadj) then
          icon(3,iv3) = it3
      else
          Write(6,*)' intins 2: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  200 continue
c
      return
      end
c=======================================================================
c
c     subroutine outhul to -
c
c     insert a point that lies outside current triangulation
c
c     February 20, 1991
c
      subroutine outhul(x, y, icon, is, iave, nicon, ivnxt, nvmax,
     *                  nvmz, ivfnd, isite, islt, iabe, eps, epz)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer icon(nicon,*), is(*), iave(*)
c
c     search for first visible site from isite
c
      ivcur = is(islt)
  100 continue
      if(icon(4,ivcur) .eq. islt) then
          isrt = icon(6,ivcur)
          ivadj = icon(2,ivcur)
      elseif(icon(5,ivcur) .eq. islt) then
          isrt = icon(4,ivcur)
          ivadj = icon(3,ivcur)
      elseif(icon(6,ivcur) .eq. islt) then
          isrt = icon(5,ivcur)
          ivadj = icon(1,ivcur)
      else
          Write(6,*)' outhul 1: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
      if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 200
      ivcur = iabs(ivadj)
      go to 100
c
  200 continue
      delx = y(isite) - y(isrt)
      dely = x(isrt) - x(isite)
      if(delx .gt. -eps .and. delx .lt. eps) delx = 0.0
      if(dely .gt. -eps .and. dely .lt. eps) dely = 0.0
      dif = sqrt(delx ** 2 + dely ** 2)
      if(dif .lt. eps) then
          Write(6,*)' outhul 1: unexpected identical sites'
          Write(6,*)' program terminated'
          stop
      endif
c
      xtemp = x(islt) - x(isite)
      ytemp = y(islt) - y(isite)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dlst=sqrt(xtemp**2+ytemp**2)
c
      xtamp = x(islt) - x(isrt)
      ytamp = y(islt) - y(isrt)
      if(xtamp .gt. -eps .and. xtamp .lt. eps) xtamp = 0.0
      if(ytamp .gt. -eps .and. ytamp .lt. eps) ytamp = 0.0
      dlgt=sqrt(xtamp**2+ytamp**2)
c
      dmax = max(dif,dlst,dlgt)
      dot = (delx * xtemp + dely * ytemp)/dmax
      if(dot .gt. -epz) go to 300
      islt = isrt
      go to 100
c
c     generate triangles with isite as a vertex
c
  300 continue
      isfn = islt
      ivlst = 0
  400 continue
      ivcur = is(islt)
      isrt = islt
      if(icon(4,ivcur) .eq. isrt) then
          indx = 3
          islt = icon(5,ivcur)
          ivadj = icon(3,ivcur)
      elseif(icon(5,ivcur) .eq. isrt) then
          indx = 1
          islt = icon(6,ivcur)
          ivadj = icon(1,ivcur)
      elseif(icon(6,ivcur) .eq. isrt) then
          indx = 2
          islt = icon(4,ivcur)
          ivadj = icon(2,ivcur)
      else
          Write(6,*)' outhul 2: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      if(ivadj .ne. 0 .and. ivadj .ne. nvmz) then
          Write(6,*)' outhul: boundary triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
      delx = y(isite) - y(isrt)
      dely = x(isrt) - x(isite)
      if(delx .gt. -eps .and. delx .lt. eps) delx = 0.0
      if(dely .gt. -eps .and. dely .lt. eps) dely = 0.0
      dif = sqrt(delx ** 2 + dely ** 2)
      if(dif .lt. eps) then
          Write(6,*)' outhul 2: unexpected identical sites'
          Write(6,*)' program terminated'
          stop
      endif
c
      xtemp = x(islt) - x(isite)
      ytemp = y(islt) - y(isite)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dlst=sqrt(xtemp**2+ytemp**2)
c
      xtamp = x(islt) - x(isrt)
      ytamp = y(islt) - y(isrt)
      if(xtamp .gt. -eps .and. xtamp .lt. eps) xtamp = 0.0
      if(ytamp .gt. -eps .and. ytamp .lt. eps) ytamp = 0.0
      dlgt=sqrt(xtamp**2+ytamp**2)
c
      dmax = max(dif,dlst,dlgt)
      dot = (delx * xtemp + dely * ytemp)/dmax
      if(dot .gt. -epz) go to 500
c
      ivfnd = iave(iabe)
      if(iabe .eq. 1) then
          ivnxt = ivfnd
          if(ivnxt .gt. nvmax) then
              Write(6,*)' outhul: parameter nvmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(1) = ivnxt + 1
      else
          iabe = iabe - 1
      endif
c
      ivcar = ivcur
      ivfad = ivfnd
      if(ivadj .eq. nvmz) then
          ivcar = -ivcar
          ivfad = -ivfad
      endif
      icon(indx,ivcur) = ivfad
      icon(1,ivfnd) = ivcar
      icon(2,ivfnd) = ivlst
      icon(4,ivfnd) = isite
      icon(5,ivfnd) = islt
      icon(6,ivfnd) = isrt
c
      if(islt.eq.isrt.or.isite.eq.islt.or.isite.eq.isrt)then
         Write(6,*)' outhul: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' islt,isrt:',islt,isrt
         Write(6,*)' isite=',isite
         stop
      endif
c
      if(ivlst .ne. 0) icon(3,ivlst) = ivfnd
      ivlst = ivfnd
      if(isrt .ne. isfn) then
          if(ivadj .eq. nvmz) go to 400
          ivfin = icon(2,ivfnd)
          ivnow = ivcur
  450     continue
          if(icon(4,ivnow) .eq. isrt) then
              ivout = icon(2,ivnow)
          elseif(icon(5,ivnow) .eq. isrt) then
              ivout = icon(3,ivnow)
          elseif(icon(6,ivnow) .eq. isrt) then
              ivout = icon(1,ivnow)
          else
              Write(6,*)' outhul: missing vertex for a triangle'
              Write(6,*)' program terminated'
              stop
          endif
          if(ivout .lt. 0) then
              is(isrt) = -ivout
              go to 400
          endif
          if(ivout .eq. ivfin) go to 400
          ivnow = ivout
          go to 450
      else
          is(isrt) = ivfnd
          go to 400
      endif
c
  500 continue
      if(ivlst .eq. 0) then
          Write(6,*)' outhul: can not insert point into triangulation'
          Write(6,*)' program terminated'
          stop
      endif
      icon(3,ivlst) = 0
      is(isite) = ivfnd
c
      return
      end
c=======================================================================
c
c     subroutine opttri to -
c
c     optimize triangulation by edge swapping with respect to a point,
c     in the direction of a triangle that has the point as a vertex
c
c     February 15, 1991
c
      subroutine opttri(x, y, icon, is, eps, epz, nicon, ivcur,
     *                  isite, nvmz)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer icon(nicon,*), is(*)
c
c     initialize
c
      isadj = icon(5,ivcur)
      iscur = icon(6,ivcur)
      isfin = isadj
c
      if(isite.eq.isadj.or.isite.eq.iscur.or.isadj.eq.iscur)then
         Write(6,*)' opttri 1: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' isite,isadj,iscur:',isite,isadj,iscur
         stop
      endif
c
  100 continue
      if(icon(1,ivcur) .le. 0) go to 200
      ivadj = icon(1,ivcur)
c
      if(icon(2,ivadj) .eq. ivcur) go to 150
      if(icon(1,ivadj) .eq. ivcur) then
          icon(1,ivadj) = icon(3,ivadj)
          icon(3,ivadj) = icon(2,ivadj)
          icon(5,ivadj) = icon(4,ivadj)
      elseif(icon(3,ivadj) .eq. ivcur) then
          icon(3,ivadj) = icon(1,ivadj)
          icon(1,ivadj) = icon(2,ivadj)
          icon(5,ivadj) = icon(6,ivadj)
      else
          Write(6,*)' opttri 1: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      icon(2,ivadj) = ivcur
      icon(4,ivadj) = isadj
      icon(6,ivadj) = iscur
  150 continue
      ishat = icon(5,ivadj)
c
      if(ishat.eq.isite.or.ishat.eq.isadj.or.ishat.eq.iscur)then
         Write(6,*)' opttri 2: unexpected equal sites'
         Write(6,*)' program terminated'
         Write(6,*)' ishat,isite,isadj,iscur:',ishat,isite,isadj,iscur
         stop
      endif
c
c     compute center of circumcircle for triangle isite-isadj-iscur
c
      delx2 = y(isadj) - y(iscur)
      dely2 = x(iscur) - x(isadj)
      if(delx2 .gt. -eps .and. delx2 .lt. eps) delx2 = 0.0
      if(dely2 .gt. -eps .and. dely2 .lt. eps) dely2 = 0.0
c
      delx3 = y(isadj) - y(isite)
      dely3 = x(isite) - x(isadj)
      if(delx3 .gt. -eps .and. delx3 .lt. eps) delx3 = 0.0
      if(dely3 .gt. -eps .and. dely3 .lt. eps) dely3 = 0.0
      dif3 = sqrt(delx3**2 + dely3**2)
      if(dif3 .lt. eps) then
          Write(6,*)' opttri 1: sites ', isadj, ' and ', isite
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
c
      delx1 = x(isite) - x(iscur)
      dely1 = y(isite) - y(iscur)
      if(delx1 .gt. -eps .and. delx1 .lt. eps) delx1 = 0.0
      if(dely1 .gt. -eps .and. dely1 .lt. eps) dely1 = 0.0
      dif1 = sqrt(delx1**2 + dely1**2)
      if(dif1 .lt. eps) then
          Write(6,*)' opttri 2: sites ', iscur, ' and ', isite
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
c
      dat = delx1*dely3 - dely1*delx3
      if(dif3 .ge. dif1) then
          dot = delx2*dely3 - dely2*delx3
          dot2 = dot/dif3
          dat2 = dat/dif3
      else
          dot = delx2*delx1 + dely2*dely1
          dot2 = dot/dif1
          dat2 = dat/dif1
      endif
c
      if(dot2 .ge. eps) then
          if(dot .ge. eps) then
              zlamb = dat/dot
          else
              zlamb = dat2/dot2
          endif
      else
          Write(6,*)' opttri: sites ',isadj,', ',iscur,' and ',isite
          Write(6,*)' appear collinear'
          Write(6,*)' program terminated'
          Write(6,*)' dot2=',dot2
          Write(6,*)' increasing tolerance epz might help'
          stop
      endif
      vx = .5*(x(isadj)+x(iscur) + zlamb*delx2)
      vy = .5*(y(isadj)+y(iscur) + zlamb*dely2)
c     dist = (vx - x(isite))**2 + (vy - y(isite))**2
c     dust = (vx - x(ishat))**2 + (vy - y(ishat))**2
c     if(dust .gt. dist-epz) go to 200
      call biscrc (x, y, ishat, isadj, epz, vx, vy, tdist)
      if(tdist.le.0.0) go to 200
      call biscrc (x, y, ishat, iscur, epz, vx, vy, tdist)
      if(tdist.le.0.0) go to 200
      call biscrc (x, y, ishat, isite, epz, vx, vy, tdist)
      if(tdist.le.0.0) go to 200
c
c     switch diagonals of quadrilateral
c
      icon(1,ivcur) = icon(3,ivadj)
      icon(3,ivadj) = ivcur
      icon(2,ivadj) = icon(2,ivcur)
      icon(2,ivcur) = ivadj
      icon(6,ivcur) = ishat
      icon(4,ivadj) = isite
      ivout = icon(1,ivcur)
      if(ivout .eq. 0 .or. ivout .eq. nvmz) go to 170
      ivcar = ivcur
      ivodj = ivadj
      if(ivout .lt. 0) then
          ivout = -ivout
          ivcar = -ivcar
          ivodj = -ivodj
      endif
      if(icon(1,ivout) .eq. ivodj) then
          icon(1,ivout) = ivcar
      elseif(icon(2,ivout) .eq. ivodj) then
          icon(2,ivout) = ivcar
      elseif(icon(3,ivout) .eq. ivodj) then
          icon(3,ivout) = ivcar
      else
          Write(6,*)' opttri 2: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  170 continue
      ivout = icon(2,ivadj)
      if(ivout .eq. 0 .or. ivout .eq. nvmz) go to 180
      ivodj = ivadj
      if(ivout .lt. 0) then
          ivout = -ivout
          ivodj = -ivodj
      endif
      icon(3,ivout) = ivodj
  180 continue
c
      if(is(isadj) .eq. ivadj) is(isadj) = ivcur
      if(is(iscur) .eq. ivcur) is(iscur) = ivadj
c
      ivcur = ivadj
      isadj = ishat
      go to 100
c
c     move to next triangle with isite as a vertex
c
  200 continue
      if(isadj .eq. isfin) go to 300
      ivadj = icon(3,ivcur)
      ivcur = ivadj
      iscur = isadj
      isadj = icon(5,ivcur)
      go to 100
c
  300 continue
      return
      end
c=======================================================================
c
c This subroutine will compute the distance from center of circle to
c bisector line between vertex of triangle and ishat
c
      subroutine biscrc (x, y, ishat, ivrt, epz, xctr, yctr, tdist)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*), norm
c
c     find midpoint of edge from ishat to ivrt
c
      xm = (x(ishat) + x(ivrt)) / 2.0
      ym = (y(ishat) + y(ivrt)) / 2.0
c
c     find vector from ivrt to ishat
c
      xu = x(ishat) - x(ivrt)
      yu = y(ishat) - y(ivrt)
c 
      norm = sqrt (xu**2 + yu**2)
      if (norm .lt. epz) then
         Write (6,*) 'biscrc: vertices ishat, ivrt appear identical'
         Write (6,*) 'program terminated'
         stop
      endif
c
      xd = xctr-xm
      yd = yctr-ym
      dif = sqrt(xd**2 + yd**2)
c
      dmax = max(norm,dif)
c
c     compute distance
c
      tdist = (xd*xu + yd*yu)/dmax
c
      return
      end
c=======================================================================
c
c     subroutine reconc to -
c
c     reconcile a line segment with current triangulation so as to
c     obtain a Delaunay triangulation constrained by the line segments
c     that have been reconciled with the triangulation so far
c
c     November 30, 1987
c
      subroutine reconc(x, y, is, icon, ifun, irun, isun, izun,
     *                  igun, iwun, ia, eps, epz, nicon, nfmax, nhmax,
     *                  nvmz, nvmp, nvmn, ileft, iright, n)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer is(*), icon(nicon,*), ia(*)
      integer ifun(*), irun(*), isun(*), izun(*)
      integer igun(*), iwun(*)
c
c     test nvmn and nvmp values
c
      if(nvmn .ne. nvmz - 1 .or. nvmp .ne. -nvmn) then
          Write(6,*)' reconc: wrong nvmn or nvmp value'
          Write(6,*)' program terminated'
          stop
      endif
c
c     compute endpoints parameters
c
      if(ileft.lt.1.or.iright.lt.1.or.ileft.gt.n.or.iright.gt.n)then
          Write(6,*)' reconc: endpoint out of range'
          Write(6,*)' program terminated'
          stop
      endif
      if(ileft .eq. iright) go to 600
      delx = y(ileft) - y(iright)
      dely = x(iright) - x(ileft)
      if(delx .gt. -eps .and. delx .lt. eps) delx = 0.0
      if(dely .gt. -eps .and. dely .lt. eps) dely = 0.0
      dif = sqrt(delx**2 + dely**2)
      if(dif .lt. eps) then
          Write(6,*)' reconc: distinct sites appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      delx = delx/dif
      dely = dely/dif
      iloft = ileft
c
   50 continue
      ivini = is(iloft)
      ivcur = ivini
c
c     test vertices of current triangle
c
  100 continue
      if(icon(4,ivcur) .eq. iloft) then
          ivadj = icon(2,ivcur)
          ivpre = icon(3,ivcur)
          ivfol = icon(1,ivcur)
          isi2 = icon(5,ivcur)
          isi1 = icon(6,ivcur)
          indx = 3
          indy = 2
      elseif(icon(5,ivcur) .eq. iloft) then
          ivadj = icon(3,ivcur)
          ivpre = icon(1,ivcur)
          ivfol = icon(2,ivcur)
          isi2 = icon(6,ivcur)
          isi1 = icon(4,ivcur)
          indx = 1
          indy = 3
      elseif(icon(6,ivcur) .eq. iloft) then
          ivadj = icon(1,ivcur)
          ivpre = icon(2,ivcur)
          ivfol = icon(3,ivcur)
          isi2 = icon(4,ivcur)
          isi1 = icon(5,ivcur)
          indx = 2
          indy = 1
      else
          Write(6,*)' reconc 1: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
      xtemp = x(isi2) - x(ileft)
      ytemp = y(isi2) - y(ileft)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dot1 = delx * xtemp + dely * ytemp
      if(dot1 .ge. epz) go to 200
      xtemp = x(isi1) - x(ileft)
      ytemp = y(isi1) - y(ileft)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dot2 = delx * xtemp + dely * ytemp
      if(dot2 .lt. epz) go to 200
      if(dot1 .le. -epz .and. isi2 .ne. iright) go to 150
      if(ivpre .lt. 0) go to 500
      if(ivpre .eq. 0) then
          icon(indx,ivcur) = nvmz
          go to 500
      endif
      icon(indx,ivcur) = -ivpre
      if(icon(1,ivpre) .eq. ivcur) then
          icon(1,ivpre) = -ivcur
      elseif(icon(2,ivpre) .eq. ivcur) then
          icon(2,ivpre) = -ivcur
      elseif(icon(3,ivpre) .eq. ivcur) then
          icon(3,ivpre) = -ivcur
      else
          Write(6,*)' reconc: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      go to 500
  150 continue
      if(2 .gt. nfmax) then
          Write(6,*)' reconc: nfmax < 2'
          Write(6,*)' program terminated'
          stop
      endif
      ifun(1) = iloft
      ifun(2) = isi2
      irun(1) = iloft
      irun(2) = isi1
      lf = 2
      lr = 2
      isun(1) = ivpre
      izun(1) = ivadj
      icon(4,ivcur) = -1
      ian = 1
      if(1 .gt. nhmax) then
          Write(6,*)' reconc: size of array ia < 1'
          Write(6,*)' program terminated'
          stop
      endif
      ia(1) = ivcur
      go to 300
c
  200 continue
      ivpro = ivcur
      ivcur = iabs(ivadj)
      if(ivcur .eq. ivini) then
          Write(6,*)' reconc: can not direct segment'
          Write(6,*)' program terminated'
          stop
      else
          if(ivadj .ne. 0 .and. ivadj .ne. nvmz) go to 100
          isi2 = isi1
          if(ivadj .eq. nvmz) go to 500
          icon(indy,ivpro) = nvmz
          go to 500
      endif
c
c     find vertex for next triangle
c
  300 continue
      if(ivfol .le. 0) then
          Write(6,*)' reconc: inappropriate crossing of line segments'
          Write(6,*)' endpoints of first segment are'
          Write(6,*)' ( ',x(isi1),' , ',y(isi1),' ) and ( ',
     *                    x(isi2),' , ',y(isi2),' )'
          Write(6,*)' endpoints of second segment are'
          Write(6,*)' ( ',x(ileft),' , ',y(ileft),' ) and ( ',
     *                    x(iright),' , ',y(iright),' )'
          Write(6,*)' program terminated'
          stop
      endif
      if(icon(4,ivfol) .eq. isi2) then
          isi3 = icon(5,ivfol)
          ivadj = icon(1,ivfol)
          ivpre = icon(3,ivfol)
      elseif(icon(5,ivfol) .eq. isi2) then
          isi3 = icon(6,ivfol)
          ivadj = icon(2,ivfol)
          ivpre = icon(1,ivfol)
      elseif(icon(6,ivfol) .eq. isi2) then
          isi3 = icon(4,ivfol)
          ivadj = icon(3,ivfol)
          ivpre = icon(2,ivfol)
      else
          Write(6,*)' reconc 2: vertex of triangle appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
      xtemp = x(isi3) - x(ileft)
      ytemp = y(isi3) - y(ileft)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dot = delx * xtemp + dely * ytemp
      icon(4,ivfol) = -1
      ian = ian + 1
      if(ian .gt. nhmax) then
          Write(6,*)' reconc: size of array ia < ian = ', ian
          Write(6,*)' program terminated'
          stop
      endif
      ia(ian) = ivfol
      if((dot .gt. -epz .and. dot .lt. epz) .or. (isi3 .eq. iright))
     * then
          lf = lf + 1
          lr = lr + 1
          if(lf .gt. nfmax) then
              Write(6,*)' reconc 1: nfmax < lf = ', lf
              Write(6,*)' program terminated'
              stop
          endif
          if(lr .gt. nfmax) then
              Write(6,*)' reconc 1: nfmax < lr = ', lr
              Write(6,*)' program terminated'
              stop
          endif
          ifun(lf) = isi3
          irun(lr) = isi3
          isun(lf-1) = ivpre
          izun(lr-1) = ivadj
          isi2 = isi3
          go to 400
      elseif(dot .ge. epz) then
          lr = lr + 1
          if(lr .gt. nfmax) then
              Write(6,*)' reconc 2: nfmax < lr = ', lr
              Write(6,*)' program terminated'
              stop
          endif
          irun(lr) = isi3
          izun(lr-1) = ivadj
          ivfol = ivpre
          isi1 = isi3
      else
          lf = lf + 1
          if(lf .gt. nfmax) then
              Write(6,*)' reconc 2: nfmax < lf = ', lf
              Write(6,*)' program terminated'
              stop
          endif
          ifun(lf) = isi3
          isun(lf-1) = ivpre
          ivfol = ivadj
          isi2 = isi3
      endif
      go to 300
c
c     retriangulate
c
  400 continue
c
c     test for eliminated adjacent triangles
c
      do 410 l = 1, lf - 1
          isam = isun(l)
          if(isam .eq. 0 .or. isam .eq. nvmz) go to 410
          if(icon(4,iabs(isam)) .ne. -1) go to 410
          if(isam .gt. 0) then
              isun(l) = nvmp
          else
              isun(l) = nvmn
          endif
  410 continue
      do 420 l = 1, lr - 1
          isam = izun(l)
          if(isam .eq. 0 .or. isam .eq. nvmz) go to 420
          if(icon(4,iabs(isam)) .ne. -1) go to 420
          if(isam .gt. 0) then
              izun(l) = nvmp
          else
              izun(l) = nvmn
          endif
  420 continue
c
c     triangulate 'right' edge star-shaped polygon
c
      mslop = 1
      mintr = 0
      iam = 0
      iopt = 3
c
      call edgstr(x, y, icon, ifun, igun, iwun, ia, eps,
     *            epz, nicon, lf, ian, iam, mslop, mintr, j, iopt)
c
      if(j .ne. 2) then
          Write(6,*)' reconc 1: unexpected j value'
          Write(6,*)' program terminated'
          stop
      endif
      ial = iwun(1)
      iap = ial
c
c     reconcile new triangles with rest of triangulation
c
      iflg = 0
      call mattri(is, icon, ifun, isun, nicon, ial, lf,
     *            nvmz, nvmp, nvmn, iflg)
c
c     triangulate 'left' edge star-shaped polygon
c
      lf = lr
      do 480 l = 1, lr - 1
          ifun(l) = irun(lr + 1 - l)
          isun(l) = izun(lr - l)
  480 continue
      ifun(lr) = irun(1)
      mslop = -1
      mintr = ian + iam + 1
c
      call edgstr(x, y, icon, ifun, igun, iwun, ia, eps,
     *            epz, nicon, lf, ian, iam, mslop, mintr, j, iopt)
c
      if(j .ne. 2) then
          Write(6,*)' reconc 2: unexpected j value'
          Write(6,*)' program terminated'
          stop
      endif
      ial = iwun(1)
      icon(2,ial) = -iap
      icon(2,iap) = -ial
c
c     reconcile new triangles with rest of triangulation
c
      iflg = 1
      call mattri(is, icon, ifun, isun, nicon, ial, lf,
     *            nvmz, nvmp, nvmn, iflg)
c
c     reinitialize to reconcile next portion of line segment
c
  500 continue
      if(isi2 .eq. iright) go to 600
      iloft = isi2
      go to 50
c
  600 continue
      return
      end
c=======================================================================
c
c     subroutine edgsrt to -
c
c     compute a (Delaunay) triangulation for a subset of the boundary
c     of one or more adjacent edge star-shaped simple polygons
c     constrained by their boundaries
c
c     December 4, 1987
c
      subroutine edgstr(x, y, icon, ifun, igun, iwun, ia, eps, epz,
     *                  nicon, lf, ian, iam, mslop, mintr, j, iopt)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer icon(nicon,*), ia(*)
      integer ifun(*), igun(*), iwun(*)
c
      iwun(1) = 0
      igun(1) = ifun(1)
      igun(2) = ifun(2)
      i = 2
      j = 2
c
  100 continue
      if(i .eq. lf) go to 400
      iwun(j) = 0
      i = i + 1
      if(i .eq. lf .and. iopt .eq. 4) iopt = 1
      j = j + 1
      igun(j) = ifun(i)
      iq1 = igun(j-2)
      iq2 = igun(j-1)
      iq3 = igun(j)
c
  200 continue
      if(i .eq. lf .and. iopt .gt. 0) go to 300
      dulx = y(iq1) - y(iq2)
      duly = x(iq2) - x(iq1)
      if(dulx .gt. -eps .and. dulx .lt. eps) dulx = 0.0
      if(duly .gt. -eps .and. duly .lt. eps) duly = 0.0
      dif = sqrt(dulx**2 + duly**2)
      if(dif .lt. eps) then
          Write(6,*)' edgstr: sites ', iq1, ' and ', iq2
          Write(6,*)' appear identical - program terminated'
          stop
      endif
      dulx = dulx/dif
      duly = duly/dif
      xtemp = x(iq3) - x(iq1)
      ytemp = y(iq3) - y(iq1)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dot = dulx * xtemp + duly * ytemp
c
      if(dot .lt. epz) go to 100
c
  300 continue
      iam = iam + 1
      if(iam .gt. ian) then
          Write(6,*)' edgstr: the Euler formula has been violated'
          Write(6,*)' program terminated'
          stop
      endif
      ial = ia(mslop*iam+mintr)
      icon(3,ial) = iwun(j-2)
      icon(1,ial) = iwun(j-1)
      icon(2,ial) = 0
      icon(4,ial) = iq1
      icon(5,ial) = iq2
      icon(6,ial) = iq3
      if(iwun(j-2) .ne. 0) icon(2,iwun(j-2)) = ial
      if(iwun(j-1) .ne. 0) icon(2,iwun(j-1)) = ial
c
c     update triangulation to obtain a Delaunay triangulation of the
c     union of the triangles obtained so far inside the polygon
c
      if(iopt .gt. 1)
     * call updtri(x, y, icon, eps, epz, nicon, ial, iq1, iq2, iq3)
c
      j = j - 1
      iwun(j-1) = ial
      igun(j) = iq3
c
      if(j .eq. 2) go to 100
      iq1 = igun(j-2)
      iq2 = igun(j-1)
      go to 200
c
  400 continue
c
      return
      end
c=======================================================================
c
c     subroutine mattri to -
c
c     reconcile triangulation of an edge star-shaped polygon with
c     rest of triangulation
c
c     December 8, 1987
c
      subroutine mattri(is, icon, ifun, isun, nicon, ial, lf,
     *                  nvmz, nvmp, nvmn, iflg)
c
      Implicit Real*8 (a-h,o-z)
      integer is(*), icon(nicon,*), ifun(*), isun(*)
c
c     reconcile triangulation with triangles outside polygon
c
      ivcur = ial
      do 400 l = 1, lf-1
          isite = ifun(l)
  100     continue
          if(icon(4,ivcur) .eq. isite) then
              ivadj = icon(3,ivcur)
              index = 3
          elseif(icon(5,ivcur) .eq. isite) then
              ivadj = icon(1,ivcur)
              index = 1
          elseif(icon(6,ivcur) .eq. isite) then
              ivadj = icon(2,ivcur)
              index = 2
          else
              Write(6,*)' mattri 1: vertex of triangle appears'
              Write(6,*)' otherwise -- program terminated'
              stop
          endif
          if(ivadj .eq. 0) go to 200
          ivcur = ivadj
          go to 100
c
  200     continue
          is(isite) = ivcur
          ivout = isun(l)
          if(ivout .eq. nvmp .or. ivout .eq. nvmn) go to 300
          isun(l) = 0
          icon(index,ivcur) = ivout
          if(ivout .eq. 0 .or. ivout .eq. nvmz) go to 400
          ivadj = ivcur
          if(ivout .lt. 0) ivadj = -ivcur
          ivabs = iabs(ivout)
          if(icon(4,ivabs) .eq. isite) then
              icon(2,ivabs) = ivadj
          elseif(icon(5,ivabs) .eq. isite) then
              icon(3,ivabs) = ivadj
          elseif(icon(6,ivabs) .eq. isite) then
              icon(1,ivabs) = ivadj
          else
              Write(6,*)' mattri 2: vertex of triangle appears'
              Write(6,*)' otherwise -- program terminated'
              stop
          endif
          go to 400
c
  300     continue
          isun(l) = ivcur
          if(ivout .eq. nvmn) isun(l) = -ivcur
c
  400 continue
c
c     reconcile triangulation with triangles inside polygon
c
      do 600 l = 1, lf-1
          if(isun(l) .eq. 0) go to 600
          ivout = isun(l)
          ivcur = iabs(ivout)
          isite = ifun(l)
          if(icon(4,ivcur) .eq. isite) then
              ivpre = icon(3,ivcur)
              ind1 = 3
          elseif(icon(5,ivcur) .eq. isite) then
              ivpre = icon(1,ivcur)
              ind1 = 1
          elseif(icon(6,ivcur) .eq. isite) then
              ivpre = icon(2,ivcur)
              ind1 = 2
          else
              Write(6,*)' mattri 3: vertex of triangle appears'
              Write(6,*)' otherwise -- program terminated'
              stop
          endif
          if(ivpre .ne. 0) then
              Write(6,*)' mattri: unexpected non-zero triangle'
              Write(6,*)' program terminated'
              stop
          endif
          if(ivcur .ne. is(isite)) go to 600
          isadj = ifun(l+1)
c
          ivadj = ivcur
  500     continue
          ivnow = iabs(ivadj)
          if(icon(4,ivnow) .eq. isite) then
              ivadj = icon(2,ivnow)
              ind2 = 2
              isnxt = icon(6,ivnow)
          elseif(icon(5,ivnow) .eq. isite) then
              ivadj = icon(3,ivnow)
              ind2 = 3
              isnxt = icon(4,ivnow)
          elseif(icon(6,ivnow) .eq. isite) then
              ivadj = icon(1,ivnow)
              ind2 = 1
              isnxt = icon(5,ivnow)
          else
              Write(6,*)' mattri 4: vertex of triangle appears'
              Write(6,*)' otherwise -- program terminated'
              stop
          endif
          if(ivadj .ne. 0) go to 500
          if(isadj .ne. isnxt) then
              Write(6,*)' mattri: unexpected non-matching vertices'
              Write(6,*)' program terminated'
              stop
          endif
          icon(ind1,ivcur) = ivnow
          if(ivout .lt. 0) icon(ind1,ivcur) = -ivnow
          icon(ind2,ivnow) = ivout
c
  600 continue
c
c     update array is
c
      if(iflg .eq. 0) then
          l1 = 2
          l2 = lf - 1
      else
          l1 = 1
          l2 = lf
      endif
      do 900 l = l1, l2
          isite = ifun(l)
          ivcur = is(isite)
          ivini = ivcur
  700     continue
          if(icon(4,ivcur) .eq. isite) then
              ivadj = icon(3,ivcur)
          elseif(icon(5,ivcur) .eq. isite) then
              ivadj = icon(1,ivcur)
          elseif(icon(6,ivcur) .eq. isite) then
              ivadj = icon(2,ivcur)
          else
              Write(6,*)' mattri: can not find vertex of triangle'
              Write(6,*)' program terminated'
              stop
          endif
          if(ivadj .eq. 0 .or. ivadj .eq. nvmz) go to 800
          if(ivadj .lt. 0) is(isite) = ivcur
          ivcur = iabs(ivadj)
          if(ivcur .eq. ivini) go to 900
          go to 700
  800     continue
          is(isite) = ivcur
  900 continue
c
      return
      end
c=======================================================================
c
c     subroutine updtri to -
c
c     update a triangulation of a generalized multiply-connected
c     polygonal region to obtain a Delaunay triangulation
c     of the region
c
c     December 10, 1987
c
      subroutine updtri(x, y, icon, eps, epz, nicon, ial,
     *                  iq1, iq2, iq3)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer icon(nicon,*)
c
c     initialize
c
      ivcur = ial
      isadj = iq1
      iscur = iq2
c
  100 continue
      if(icon(3,ivcur) .eq. 0) go to 200
      ivadj = icon(3,ivcur)
      ishat = icon(5,ivadj)
c
c     compute center of circumcircle for triangle isadj-iscur-iq3
c
      delx2 = y(isadj) - y(iscur)
      dely2 = x(iscur) - x(isadj)
      if(delx2 .gt. -eps .and. delx2 .lt. eps) delx2 = 0.0
      if(dely2 .gt. -eps .and. dely2 .lt. eps) dely2 = 0.0
      dif = sqrt(delx2**2 + dely2**2)
      if(dif .lt. eps) then
          Write(6,*)' updtri 1: sites ', isadj, ' and ', iscur
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
      delx2 = delx2/dif
      dely2 = dely2/dif
c
      delx3 = y(isadj) - y(iq3)
      dely3 = x(iq3) - x(isadj)
      if(delx3 .gt. -eps .and. delx3 .lt. eps) delx3 = 0.0
      if(dely3 .gt. -eps .and. dely3 .lt. eps) dely3 = 0.0
      dif = sqrt(delx3**2 + dely3**2)
      if(dif .lt. eps) then
          Write(6,*)' updtri 2: sites ', isadj, ' and ', iq3
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
      delx3 = delx3/dif
      dely3 = dely3/dif
c
      delx1 = x(iq3) - x(iscur)
      dely1 = y(iq3) - y(iscur)
      if(delx1 .gt. -eps .and. delx1 .lt. eps) delx1 = 0.0
      if(dely1 .gt. -eps .and. dely1 .lt. eps) dely1 = 0.0
      dat = .5*(delx1*dely3 - dely1*delx3)
      dot = delx2*dely3 - dely2*delx3
      if(dot .ge. eps .or. dot .le. -eps) then
          zlamb = dat/dot
      else
          Write(6,*)' updtri: sites ',isadj,', ',iscur,' and ',iq3
          Write(6,*)' appear colinear -- program terminated'
          stop
      endif
      vx = .5*(x(isadj)+x(iscur)) + zlamb*delx2
      vy = .5*(y(isadj)+y(iscur)) + zlamb*dely2
      dist = (vx - x(iq3))**2 + (vy - y(iq3))**2
      dust = (vx - x(ishat))**2 + (vy - y(ishat))**2
      if(dust .gt. dist-epz) go to 200
c
c     switch diagonals of quadrilateral
c
      icon(3,ivcur) = icon(3,ivadj)
      icon(3,ivadj) = icon(1,ivadj)
      icon(1,ivadj) = icon(1,ivcur)
      icon(1,ivcur) = ivadj
      icon(5,ivcur) = ishat
      icon(4,ivadj) = ishat
      icon(5,ivadj) = iscur
      icon(6,ivadj) = iq3
      if(icon(3,ivcur) .ne. 0) icon(2,icon(3,ivcur)) = ivcur
      if(icon(1,ivadj) .ne. 0) icon(2,icon(1,ivadj)) = ivadj
c
      ivcur = ivadj
      isadj = ishat
      go to 100
c
c     move to next triangle with iq3 as a vertex
c
  200 continue
      if(isadj .eq. iq1) go to 300
      ivadj = icon(2,ivcur)
      ivcur = ivadj
      iscur = isadj
      isadj = icon(4,ivcur)
      go to 100
c
  300 continue
c
      return
      end
c=======================================================================
c
c     subroutine deleli to -
c
c     eliminate a vertex from a Delaunay triangulation so as to obtain
c     a Delaunay triangulation that does not include the vertex
c     (vertex can not be in an inserted edge)
c
c     February 28, 1991
c
      subroutine deleli(x, y, is, icon, ifun, irun, isun, izun, igun,
     *                  iwun, ia, iave, iaze, eps, epz, nicon, nfmax,
     *                  nhmax, njmax, nkmax, nvmz, nvmp, nvmn, isite,
     *                  n, isucs, isnew, iabe, iase)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer is(*), icon(nicon,*), ia(*)
      integer ifun(*), irun(*), isun(*), izun(*)
      integer igun(*), iwun(*), iave(*), iaze(*)
c
c     test nvmn and nvmp values
c
      if(nvmn .ne. nvmz - 1 .or. nvmp .ne. -nvmn) then
          Write(6,*)' deleli: wrong nvmn or nvmp value'
          Write(6,*)' program terminated'
          stop
      endif
c
c     obtain 'right' side information
c
      isucs = 1
      if(n - iase + 1 .le. 3) then
          isucs = -2
          go to 2300
      endif
      if(isite .lt. 1 .or. isite .gt. n) then
          Write(6,*)' deleli: point to be eliminated is out of range'
          Write(6,*)' program terminated'
          stop
      endif
      ivini = is(isite)
      if(icon(4,ivini) .eq. isite) then
          ivpre = icon(1,ivini)
          ispre = icon(5,ivini)
          ivout = icon(3,ivini)
          ivadj = icon(2,ivini)
          isadj = icon(6,ivini)
      elseif(icon(5,ivini) .eq. isite) then
          ivpre = icon(2,ivini)
          ispre = icon(6,ivini)
          ivout = icon(1,ivini)
          ivadj = icon(3,ivini)
          isadj = icon(4,ivini)
      elseif(icon(6,ivini) .eq. isite) then
          ivpre = icon(3,ivini)
          ispre = icon(4,ivini)
          ivout = icon(2,ivini)
          ivadj = icon(1,ivini)
          isadj = icon(5,ivini)
      else
          Write(6,*)' deleli 1: missing triangle vertex'
          Write(6,*)' program terminated'
          stop
      endif
c
      if(ivout .lt. 0) then
          isucs = 0
          go to 2300
      endif
      if(ivout .eq. 0) go to 100
c
c     compute 'splitting' vector
c
      lf = 0
      delx = y(ispre) - y(isite)
      dely = x(isite) - x(ispre)
      if(delx .gt. -eps .and. delx .lt. eps) delx = 0.0
      if(dely .gt. -eps .and. dely .lt. eps) dely = 0.0
      dif = sqrt(delx ** 2 + dely ** 2)
      if(dif .lt. eps) then
          Write(6,*)' reconc: distinct sites appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      delx = delx/dif
      dely = dely/dif
c
c     eliminate first triangle
c
  100 continue
      lr = 2
      if(lr .gt. nfmax) then
          Write(6,*)' deleli 1: parameter nfmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      irun(1) = ispre
      irun(2) = isadj
      izun(1) = ivpre
      icon(4,ivini) = -1
      ian = 1
      if(ian .gt. nhmax) then
          Write(6,*)' deleli 1: parameter nhmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      ia(1) = ivini
c
c     analyze next triangle if any
c
  200 continue
      if(ivadj .lt. 0) then
          isucs = 0
          go to 2300
      endif
      if(ivadj .eq. 0) go to 700
      ivcur = ivadj
      if(icon(4,ivcur) .eq. isite) then
          ivpre = icon(1,ivcur)
          ivadj = icon(2,ivcur)
          isadj = icon(6,ivcur)
      elseif(icon(5,ivcur) .eq. isite) then
          ivpre = icon(2,ivcur)
          ivadj = icon(3,ivcur)
          isadj = icon(4,ivcur)
      elseif(icon(6,ivcur) .eq. isite) then
          ivpre = icon(3,ivcur)
          ivadj = icon(1,ivcur)
          isadj = icon(5,ivcur)
      else
          Write(6,*)' deleli 2: missing triangle vertex'
          Write(6,*)' program terminated'
          stop
      endif
      if(ivout .eq. 0) go to 300
      if(lf .ne. 0) go to 400
c
c     ascertain which side new vertex is on
c
      xtemp = x(isadj) - x(ispre)
      ytemp = y(isadj) - y(ispre)
      if(xtemp .gt. -eps .and. xtemp .lt. eps) xtemp = 0.0
      if(ytemp .gt. -eps .and. ytemp .lt. eps) ytemp = 0.0
      dot = delx * xtemp + dely * ytemp
      if(dot .ge. epz) go to 400
c
c     add next vertex for 'right' side
c
  300 continue
      lr = lr + 1
      if(lr .gt. nfmax) then
          Write(6,*)' deleli 2: parameter nfmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      irun(lr) = isadj
      izun(lr-1) = ivpre
      icon(4,ivcur) = -1
      ian = ian + 1
      if(ian .gt. nhmax) then
          Write(6,*)' deleli 2: parameter nhmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      ia(ian) = ivcur
      go to 200
c
c     obtain 'left' side information
c
  400 continue
      lf = lf + 1
      if(lf .gt. nfmax) then
          Write(6,*)' deleli 3: parameter nfmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      ifun(lf) = isadj
      isun(lf) = ivpre
      icon(4,ivcur) = -1
      ian = ian + 1
      if(ian .gt. nhmax) then
          Write(6,*)' deleli 3: parameter nhmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      ia(ian) = ivcur
      if(isadj .eq. ispre) go to 500
      go to 200
  500 continue
      if(lf .lt. 2) then
          Write(6,*)' deleli: not enough left side information'
          Write(6,*)' program terminated'
          stop
      endif
      ivtro = isun(1)
      do 600 l = 2, lf
          isun(l-1) = isun(l)
  600 continue
c
c     test for eliminated adjacent triangles
c
  700 continue
      do 800 l = 1, lr - 1
          isam = izun(l)
          if(isam .eq. 0 .or. isam .eq. nvmz) go to 800
          if(icon(4,iabs(isam)) .ne. -1) go to 800
          if(isam .gt. 0) then
              izun(l) = nvmp
          else
              izun(l) = nvmn
          endif
  800 continue
      mslop = 1
      mintr = 0
      iam = 0
      isnew = irun(lr)
      if(ivout .eq. 0) go to 1000
      do 900 l = 1, lf -1
          isam = isun(l)
          if(isam .eq. 0 .or. isam .eq. nvmz) go to 900
          if(icon(4,iabs(isam)) .ne. -1) go to 900
          if(isam .gt. 0) then
              isun(l) = nvmp
          else
              isun(l) = nvmn
          endif
  900 continue
      go to 1200
c
c     point to be eliminated is on convex hull
c
 1000 continue
      lr = lr + 1
      if(lr .gt. nfmax) then
          Write(6,*)' deleli 4: parameter nfmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      irun(lr) = isite
      izun(lr-1) = 0
c
c     triangulate 'right' side in a Delaunay
c
      iopt = 4
c
      call edgstr(x, y, icon, irun, igun, iwun, ia, eps, epz, nicon,
     *            lr, ian, iam, mslop, mintr, jr, iopt)
c
      ial = iwun(1)
      is(isite) = ial
c
c     reconcile with rest of triangulation
c
      iflg = 1
      call mattri(is, icon, irun, izun, nicon, ial, lr, nvmz, nvmp,
     *            nvmn, iflg)
c
c     eliminate superfluous triangles
c
 1100 continue
      ivcur = icon(3,ial)
      ivabs = iabs(ivcur)
      iscor = icon(5,ial)
      if(icon(4,ivabs) .eq. iscor) then
          indx = 3
      elseif(icon(5,ivabs) .eq. iscor) then
          indx = 1
      elseif(icon(6,ivabs) .eq. iscor) then
          indx = 2
      else
          Write(6,*)' deleli: triangle adjacency problem'
          Write(6,*)' program terminated'
          stop
      endif
      icon(indx,ivabs) = 0
      if(ivcur .lt. 0) icon(indx,ivabs) = nvmz
      is(iscor) = ivabs
      icon(4,ial) = -1
      iabe = iabe + 1
      if(iabe .gt. njmax) then
          Write(6,*)' deleli 1: parameter njmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      iave(iabe) = ial
      ial = icon(1,ial)
      if(ial .ne. 0) go to 1100
      go to 2200
c
c     point to be eliminated is not on convex hull
c
 1200 continue
c
c     triangulate 'right' side (not necessarily Delaunay)
c
      iopt = 0
c
      call edgstr(x, y, icon, irun, igun, iwun, ia, eps, epz, nicon,
     *            lr, ian, iam, mslop, mintr, jr, iopt)
c
c     update 'left' side information
c
      lg = lf
      do 1300 j = 2, jr
          lf = lf + 1
          if(lf .gt. nfmax) then
              Write(6,*)' deleli 5: parameter nfmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          ifun(lf) = igun(j)
          isun(lf-1) = iwun(j-1)
 1300 continue
c
c     triangulate 'left' side in a Delaunay fashion
c
      iopt = 2
c
      call edgstr(x, y, icon, ifun, igun, iwun, ia, eps, epz, nicon,
     *            lf, ian, iam, mslop, mintr, jf, iopt)
c
      if(jf .ne. 2) then
          Write(6,*)' deleli: unexpected jf value'
          Write(6,*)' program terminated'
          stop
      endif
      ial = iwun(1)
c
c     merge the two sides
c
      lf = lg
      do 1400 j = 1, jr - 1
          iwun(j) = isun(lf+j-1)
 1400 continue
c
      j = jr
      ivcur = ial
 1500 continue
      ivadj = icon(1,ivcur)
      if(ivadj .eq. 0) go to 1600
      ivcur = ivadj
      go to 1500
 1600 continue
      j = j - 1
      ivadj = iwun(j)
      icon(1,ivcur) = ivadj
      if(ivadj .ne. 0) icon(2,ivadj) = ivcur
      ivcur = icon(3,ivcur)
      if(j .ne. 1) go to 1500
c
c     add each triangle to 'left' polygon and optimize
c
      j = jr - 1
 1700 continue
      if(j .eq. 0) go to 1900
      ivcur = iwun(j)
      if(ivcur .eq. 0) go to 1800
      iwun(j) = icon(3,ivcur)
      j = j + 1
      iwun(j) = icon(1,ivcur)
      indx = 2
      isnow = icon(5,ivcur)
      isadj = icon(6,ivcur)
      iscur = icon(4,ivcur)
c
      call addtri(x, y, icon, eps, epz, nicon, ivcur, indx, isnow,
     *            isadj, iscur)
c
      go to 1700
 1800 continue
      j = j - 1
      go to 1700
 1900 continue
c
c     reconcile everything with rest of triangulation
c
      if(icon(4,ial) .ne. ifun(1) .or. icon(6,ial) .ne. isnew) then
          Write(6,*)' deleli: can not find initial triangle'
          Write(6,*)' program terminated'
          stop
      endif
      icon(2,ial) = ivtro
      if(ivtro .eq. 0 .or. ivtro .eq. nvmz) go to 2000
      iel = ial
      if(ivtro .lt. 0) then
          ivtro = -ivtro
          iel = -iel
      endif
      if(icon(4,ivtro) .eq. isnew) then
          icon(2,ivtro) = iel
      elseif(icon(5,ivtro) .eq. isnew) then
          icon(3,ivtro) = iel
      elseif(icon(6,ivtro) .eq. isnew) then
          icon(1,ivtro) = iel
      else
          Write(6,*)' deleli 3: missing triangle vertex'
          Write(6,*)' program terminated'
          stop
      endif
c
 2000 continue
      is(isnew) = ial
      do 2100 l = 2, lr
          lf = lf + 1
          if(lf .gt. nfmax) then
              Write(6,*)' deleli 6: parameter nfmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          ifun(lf) = irun(l)
          isun(lf-1) = izun(l-1)
 2100 continue
c
      iflg = 1
      call mattri(is, icon, ifun, isun, nicon, ial, lf, nvmz, nvmp,
     *            nvmn, iflg)
c
c     record eliminated triangles
c
      if(iam .ge. ian) then
          Write(6,*)' deleli: Euler formula has been violated'
          Write(6,*)' program terminated'
          stop
      endif
      iam = iam + 1
      do 2150 i = iam, ian
          iabe = iabe + 1
          if(iabe .gt. njmax) then
              Write(6,*)' deleli 2: parameter njmax exceeded'
              Write(6,*)' program terminated'
              stop
          endif
          iave(iabe) = ia(i)
 2150 continue
c
c     record eliminated point
c
 2200 continue
      is(isite) = 0
      iase = iase + 1
      if(iase .gt. nkmax) then
          Write(6,*)' deleli: parameter nkmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
      iaze(iase) = isite
c
 2300 continue
      return
      end
c=======================================================================
c
c     subroutine addtri to -
c
c     add a triangle to a Delaunay triangulation of a polygon and get
c     a new Delaunay triangulation of the new polygon if the new
c     triangle has only one side in common with the polygon
c
c     March 6, 1991
c
      subroutine addtri(x, y, icon, eps, epz, nicon, ivcur, indx,
     *                  isite, isadj, iscur)
c
      Implicit Real*8 (a-h,o-z)
      real*8 x(*), y(*)
      integer icon(nicon,*)
c
c     initialize
c
      isfin = isadj
c
  100 continue
      ivadj = icon(indx,ivcur)
      if(ivadj .eq. 0) go to 200
c
      if(icon(1,ivadj) .eq. ivcur) then
          indy = 4
      elseif(icon(2,ivadj) .eq. ivcur) then
          indy = 5
      elseif(icon(3,ivadj) .eq. ivcur) then
          indy = 6
      else
          Write(6,*)' addtri: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      ishat = icon(indy,ivadj)
c
c     compute center of circumcircle for triangle isite-isadj-iscur
c
      delx2 = y(isadj) - y(iscur)
      dely2 = x(iscur) - x(isadj)
      if(delx2 .gt. -eps .and. delx2 .lt. eps) delx2 = 0.0
      if(dely2 .gt. -eps .and. dely2 .lt. eps) dely2 = 0.0
      dif = sqrt(delx2**2 + dely2**2)
      if(dif .lt. eps) then
          Write(6,*)' addtri 1: sites ', isadj, ' and ', iscur
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
      delx2 = delx2/dif
      dely2 = dely2/dif
c
      delx3 = y(isadj) - y(isite)
      dely3 = x(isite) - x(isadj)
      if(delx3 .gt. -eps .and. delx3 .lt. eps) delx3 = 0.0
      if(dely3 .gt. -eps .and. dely3 .lt. eps) dely3 = 0.0
      dif = sqrt(delx3**2 + dely3**2)
      if(dif .lt. eps) then
          Write(6,*)' addtri 2: sites ', isadj, ' and ', isite
          Write(6,*)' appear identical -- program terminated'
          stop
      endif
      delx3 = delx3/dif
      dely3 = dely3/dif
c
      delx1 = x(isite) - x(iscur)
      dely1 = y(isite) - y(iscur)
      if(delx1 .gt. -eps .and. delx1 .lt. eps) delx1 = 0.0
      if(dely1 .gt. -eps .and. dely1 .lt. eps) dely1 = 0.0
      dat = .5*(delx1*dely3 - dely1*delx3)
      dot = delx2*dely3 - dely2*delx3
      if(dot .ge. eps .or. dot .le. -eps) then
          zlamb = dat/dot
      else
          Write(6,*)' addtri: sites ',isadj,', ',iscur,' and ',isite
          Write(6,*)' appear colinear -- program terminated'
          stop
      endif
      vx = .5*(x(isadj)+x(iscur)) + zlamb*delx2
      vy = .5*(y(isadj)+y(iscur)) + zlamb*dely2
      dist = (vx - x(isite))**2 + (vy - y(isite))**2
      dust = (vx - x(ishat))**2 + (vy - y(ishat))**2
      if(dust .gt. dist-epz) go to 200
c
c     switch diagonals of quadrilateral
c
      iv1 = ivadj
      iv2 = ivcur
      if(indy .eq. 5) then
          iv1 = ivcur
          iv2 = ivadj
      endif
      if(indx .ne. 1 .and. indy .ne. 4) then
          icon(3,iv1) = icon(3,iv2)
          icon(3,iv2) = icon(1,iv2)
          icon(1,iv2) = icon(1,iv1)
          icon(1,iv1) = iv2
          icon(5,iv1) = icon(5,iv2)
          icon(4,iv2) = icon(5,iv2)
          icon(5,iv2) = icon(6,iv2)
          icon(6,iv2) = icon(6,iv1)
          if(icon(3,iv1) .ne. 0) icon(2,icon(3,iv1)) = iv1
          if(icon(1,iv2) .ne. 0) icon(2,icon(1,iv2)) = iv2
          ivcur = ivadj
      else
          icon(1,iv1) = icon(1,iv2)
          icon(1,iv2) = icon(3,iv2)
          icon(3,iv2) = icon(3,iv1)
          icon(3,iv1) = iv2
          icon(5,iv1) = icon(5,iv2)
          icon(6,iv2) = icon(5,iv2)
          icon(5,iv2) = icon(4,iv2)
          icon(4,iv2) = icon(4,iv1)
          if(icon(1,iv1) .ne. 0) icon(2,icon(1,iv1)) = iv1
          if(icon(3,iv2) .ne. 0) icon(2,icon(3,iv2)) = iv2
      endif
      if(indy .eq. 4) indx = 3
      isadj = ishat
      go to 100
c
c     move to next triangle with isite as a vertex
c
  200 continue
      if(isadj .eq. isfin) go to 300
      if(indx .eq. 1) then
          ivadj = icon(3,ivcur)
      elseif(indx .eq. 2) then
          ivadj = icon(1,ivcur)
      else
          ivadj = icon(2,ivcur)
      endif
      ivcur = ivadj
      iscur = isadj
      if(icon(4,ivcur) .eq. isite) then
          indx = 1
          isadj = icon(5,ivcur)
      elseif(icon(5,ivcur) .eq. isite) then
          indx = 2
          isadj = icon(6,ivcur)
      elseif(icon(6,ivcur) .eq. isite) then
          indx = 3
          isadj = icon(4,ivcur)
      else
          Write(6,*)' addtri: missing vertex of triangle'
          Write(6,*)' program terminated'
          stop
      endif
      go to 100
c
  300 continue
      return
      end
c=======================================================================
c
c     subroutine comprs to -
c
c     compress array icon so that there are no gaps in it
c
c     March 15, 1991
c
      subroutine comprs(is, icon, id, n, ivnxt, nvmax, nicon, nvmz)
c
      Implicit Real*8 (a-h,o-z)
      integer is(*), icon(nicon,*), id(*)
c
c     test n, ivnxt, nvmz
c
      if(n .gt. nvmax .or. ivnxt .gt. nvmax) then
          Write(6,*)' comprs: parameter nvmax exceeded'
          Write(6,*)' program terminated'
          stop
      endif
c
      if(nvmz .ne. -nvmax - 1) then
          Write(6,*)' comprs: wrong nvmz value'
          Write(6,*)' program terminated'
          stop
      endif
c
c     compress icon
c
      ivnot = 0
      do 200 i = 1, ivnxt
          if(icon(4,i) .eq. -1) go to 200
          ivnot = ivnot + 1
          id(i) = ivnot
          do 100 j = 1, 6
              icon(j,ivnot) = icon(j,i)
  100     continue
  200 continue
      ivnxt = ivnot
c
c     update icon for triangles
c
      do 400 i = 1, ivnxt
          do 300 j = 1, 3
              ivcur = icon(j,i)
              if(ivcur .eq. 0 .or. ivcur .eq. nvmz) go to 300
              if(ivcur .gt. 0) then
                  icon(j,i) = id(ivcur)
              else
                  ivcur = -ivcur
                  icon(j,i) = -id(ivcur)
              endif
  300     continue
  400 continue
c
c     update is for triangles
c
      do 500 i = 1, n
          if(is(i) .le. 0) go to 500
          is(i) = id(is(i))
  500 continue
c
      return
      end
c=======================================================================
c
c     subroutine contst to -
c
c     test consistency of triangulation
c
c     November 9, 1993
c
      subroutine contst(icon, is, id, nv, ivnxt, nicon, nvmz)
c
      Implicit Real*8 (a-h,o-z)
      integer icon(nicon,*), is(*), id(*)
      integer site1, site2
c
c     test initial triangle for each site
c
      do 50 i = 1, nv
          iscur = is(i)
          if(iscur .le. 0) go to 50
          if(icon(4,iscur) .ne. i .and. icon(5,iscur) .ne. i .and.
     *       icon(6,iscur) .ne. i) then
              Write(6,*)' contst: unexpected triangle vertex'
              Write(6,*)' program terminated'
              stop
          endif
   50 continue
c
c     initialize
c
      do 60 i = 1, ivnxt
          id(i) = 1
   60 continue
c
      do 70 i = 1, nv
          if(is(i) .gt. 0) go to 80
   70 continue
      Write(6,*)' contst: no initial triangle'
      Write(6,*)' program terminated'
      stop
   80 continue
      ivini = is(i)
      itrct = 1
      id(ivini) = 0
      ivcur = ivini
  100 continue
      ivcur2 = icon(2,ivcur)
      if(ivcur2 .eq. 0 .or. ivcur2 .eq. nvmz) go to 200
      ivabs2 = iabs(ivcur2)
      ivlst = ivcur
      if(ivcur2 .lt. 0) ivlst = -ivcur
      site1 = icon(6,ivcur)
      site2 = icon(4,ivcur)
      if(id(ivabs2) .eq. 1) go to 500
      if(icon(2,ivabs2) .eq. ivlst) then
          isit1 = icon(4,ivabs2)
          isit2 = icon(6,ivabs2)
      elseif(icon(3,ivabs2) .eq. ivlst) then
          isit1 = icon(5,ivabs2)
          isit2 = icon(4,ivabs2)
      elseif(icon(1,ivabs2) .eq. ivlst .and. ivabs2 .eq. ivini) then
          isit1 = icon(6,ivabs2)
          isit2 = icon(5,ivabs2)
      else
          Write(6,*)' contst 1: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      if(site1 .ne. isit1 .or. site2 .ne. isit2) then
          Write(6,*)' contst 1: triangle vertex appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  200 continue
      ivcur2 = icon(3,ivcur)
      if(ivcur2 .eq. 0 .or. ivcur2 .eq. nvmz) go to 300
      ivabs2 = iabs(ivcur2)
      ivlst = ivcur
      if(ivcur2 .lt. 0) ivlst = -ivcur
      site1 = icon(4,ivcur)
      site2 = icon(5,ivcur)
      if(id(ivabs2) .eq. 1) go to 500
      if(icon(2,ivabs2) .eq. ivlst) then
          isit1 = icon(4,ivabs2)
          isit2 = icon(6,ivabs2)
      elseif(icon(3,ivabs2) .eq. ivlst) then
          isit1 = icon(5,ivabs2)
          isit2 = icon(4,ivabs2)
      elseif(icon(1,ivabs2) .eq. ivlst .and. ivabs2 .eq. ivini) then
          isit1 = icon(6,ivabs2)
          isit2 = icon(5,ivabs2)
      else
          Write(6,*)' contst 2: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      if(site1 .ne. isit1 .or. site2 .ne. isit2) then
          Write(6,*)' contst 2: triangle vertex appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
  300 continue
      ivcur2 = icon(1,ivcur)
      if(ivcur .eq. ivini) go to 400
      if(ivcur2 .eq. 0 .or. ivcur2 .eq. nvmz) then
          Write(6,*)' contst: unexpected adjacent triangle'
          Write(6,*)' program terminated'
          stop
      endif
      ivabs2 = iabs(ivcur2)
      ivlst = ivcur
      if(ivcur2 .lt. 0) ivlst = -ivcur
      if(icon(3,ivabs2) .eq. ivlst) then
          ivcur = ivabs2
          go to 300
      elseif(icon(2,ivabs2) .eq. ivlst) then
          ivcur = ivabs2
          go to 200
      elseif(icon(1,ivabs2) .eq. ivlst .and. ivabs2 .eq. ivini) then
          go to 900
      else
          Write(6,*)' contst: unexpected adjacency'
          Write(6,*)' program terminated'
          stop
      endif
  400 continue
      if(ivcur2 .eq. 0 .or. ivcur2 .eq. nvmz) go to 900
      ivabs2 = iabs(ivcur2)
      ivlst = ivcur
      if(ivcur2 .lt. 0) ivlst = -ivcur
      site1 = icon(5,ivcur)
      site2 = icon(6,ivcur)
      if(id(ivabs2) .eq. 1) go to 500
      if(icon(2,ivabs2) .eq. ivlst) then
          isit1 = icon(4,ivabs2)
          isit2 = icon(6,ivabs2)
      elseif(icon(3,ivabs2) .eq. ivlst) then
          isit1 = icon(5,ivabs2)
          isit2 = icon(4,ivabs2)
      else
          Write(6,*)' contst 3: adjacent triangles appear otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      if(site1 .ne. isit1 .or. site2 .ne. isit2) then
          Write(6,*)' contst 3: triangle vertex appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
      go to 900
c
c     shift icon for ivabs2
c
  500 continue
      if(icon(1,ivabs2) .eq. ivlst) go to 600
      if(icon(2,ivabs2) .eq. ivlst) then
          itemp = icon(1,ivabs2)
          icon(1,ivabs2) = icon(2,ivabs2)
          icon(2,ivabs2) = icon(3,ivabs2)
          icon(3,ivabs2) = itemp
          itemp = icon(4,ivabs2)
          icon(4,ivabs2) = icon(5,ivabs2)
          icon(5,ivabs2) = icon(6,ivabs2)
          icon(6,ivabs2) = itemp
      elseif(icon(3,ivabs2) .eq. ivlst) then
          itemp = icon(1,ivabs2)
          icon(1,ivabs2) = icon(3,ivabs2)
          icon(3,ivabs2) = icon(2,ivabs2)
          icon(2,ivabs2) = itemp
          itemp = icon(4,ivabs2)
          icon(4,ivabs2) = icon(6,ivabs2)
          icon(6,ivabs2) = icon(5,ivabs2)
          icon(5,ivabs2) = itemp
      else
          Write(6,*)' contst: can not find triangle'
          Write(6,*)' program terminated'
          stop
      endif
c
      if(site1 .ne. icon(6,ivabs2) .or.
     *   site2 .ne. icon(5,ivabs2)) then
          Write(6,*)' contst 4: triangle vertex appears otherwise'
          Write(6,*)' program terminated'
          stop
      endif
c
  600 continue
      ivcur = ivabs2
      id(ivcur) = 0
      itrct = itrct + 1
      go to 100
c
  900 continue
      if(itrct .ne. ivnxt) then
          Write(6,*)' contst: unexpected number of triangles'
          Write(6,*)' program terminated'
          stop
      endif
c
      do 1000 i = 1, ivnxt
          if(id(i) .ne. 0) then
              Write(6,*)' contst: unexpected triangle mark'
              Write(6,*)' program terminated'
              stop
          endif
 1000 continue
c
c      Write(6,*)'****************************************'
c      Write(6,*)'consistency check satisfied'
c      Write(6,*)'****************************************'
      return
      end
c=======================================================================
c
c     this subroutine tests how well the circle criterion is satisfied
c     by the triangles
c
      subroutine circri (icon, ivnxt, nicon, x, y, epz, eps)
c
      Implicit Real*8 (a-h,o-z)
      integer icon(nicon,*)
      real*8 x(*), y(*)
c
c     initialize
c
      dmax = 0.0
c
c     test all triangles
c
      do 200 i = 1, ivnxt
         isite = icon(4,i)
         isadj = icon(5,i)
         iscur = icon(6,i)
c
c    compute center for current triangle
c
         delx2 = y(isadj) - y(iscur)
         dely2 = x(iscur) - x(isadj)
         if(delx2 .gt. -eps .and. delx2 .lt. eps) delx2 = 0.0
         if(dely2 .gt. -eps .and. dely2 .lt. eps) dely2 = 0.0
c
         delx3 = y(isadj) - y(isite)
         dely3 = x(isite) - x(isadj)
         if(delx3 .gt. -eps .and. delx3 .lt. eps) delx3 = 0.0
         if(dely3 .gt. -eps .and. dely3 .lt. eps) dely3 = 0.0
         dif3 = sqrt(delx3**2 + dely3**2)
         if(dif3 .lt. eps) then
            Write(6,*)' circri 1: sites ', isadj, ' and ', isite
            Write(6,*)' appear identical -- program terminated'
            stop
         endif
c
         delx1 = x(isite) - x(iscur)
         dely1 = y(isite) - y(iscur)
         if(delx1 .gt. -eps .and. delx1 .lt. eps) delx1 = 0.0
         if(dely1 .gt. -eps .and. dely1 .lt. eps) dely1 = 0.0
         dif1 = sqrt(delx1**2 + dely1**2)
         if(dif1 .lt. eps) then
            Write(6,*)' circri 2: sites ', iscur, ' and ', isite
            Write(6,*)' appear identical -- program terminated'
            stop
         endif
c
         dat = delx1*dely3 - dely1*delx3
         if(dif3 .ge. dif1) then
            dot = delx2*dely3 - dely2*delx3
            dot2 = dot/dif3
            dat2 = dat/dif3
         else
            dot = delx2*delx1 + dely2*dely1
            dot2 = dot/dif1
            dat2 = dat/dif1
         endif
c
         if(dot2 .ge. eps) then
            if(dot .ge. eps) then
               zlamb = dat/dot
            else
               zlamb = dat2/dot2
            endif
         else
            Write(6,*)' circri: sites ',isadj,', ',iscur,' and ',isite
            Write(6,*)' appear collinear'
            Write(6,*)' program terminated'
            Write(6,*)' dot2=',dot2
            Write(6,*)' increasing tolerance epz might help'
            stop
         endif
         vx = .5*(x(isadj)+x(iscur) + zlamb*delx2)
         vy = .5*(y(isadj)+y(iscur) + zlamb*dely2)
c
         do 100 j = 1, 3
            ivadj = icon(j,i)
            if(ivadj.le.0) go to 100
            do 50 k = 1, 3
                 if(icon(k,ivadj).eq.i) go to 75
   50       continue
            Write(6,*)' circri: unexpected adjacent triangles'
            Write(6,*)' program terminated'
            stop
   75       continue
            ishat = icon(k+3,ivadj)
c
            call biscrc(x,y,ishat,isadj,epz,vx,vy,tdist1)
            if(tdist1.le.0.0) go to 100
            call biscrc(x,y,ishat,iscur,epz,vx,vy,tdist2)
            if(tdist2.le.0.0) go to 100
            call biscrc(x,y,ishat,isite,epz,vx,vy,tdist3)
            if(tdist3.le.0.0) go to 100
            tdist=min(tdist1,tdist2,tdist3)
            if(tdist.gt.dmax) dmax=tdist
  100    continue
  200 continue
c
c      Write(6,*) '******************************************'
c      Write(6,*) 'Circle criterion error = ', dmax
c      Write(6,*) '(0.0 is the desired error)'
c      Write(6,*) '******************************************'
c      Write(6,*) ''
c
      return
      end
c=======================================================================
      subroutine indexi(n,iarrin,indx)
      dimension iarrin(n),indx(n)
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=iarrin(indxt)
        else
          indxt=indx(ir)
          q=iarrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(iarrin(indx(j)).lt.iarrin(indx(j+1)))j=j+1
          endif
          if(q.lt.iarrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
      end
c======================================================================
      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt,SS)
      Implicit Real*8 (a-h,o-z)
      INTEGER n
      REAL*8 adev,ave,curt,sdev,skew,var,data(n),SS(n)
      INTEGER j
      REAL*8 p,s,ep
      if(n.le.1) Then
        Write(6,*)'n must be at least 2 in moment'
	Read(5,*)
      endif	
      s=0.
      stotal = 0.d0
      do 11 j=1,n
        s=s+data(j)*SS(n)
        stotal = stotal + SS(n)
11    continue
      ave=s/stotal
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
        Write(6,*) 'no skew or kurtosis when zero variance in moment'
	Read(5,*) 
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
c======================================================================
