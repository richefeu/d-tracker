c23456789012345678901234567890123456789012345678901234567890123456789012
c=======================================================================
c Gael COMBE, 5 avril 2020
c
c a drawing program that needs pgplot to be compiled
c
c
c To compile (it can work if pgplot is installed)
c   gfortran -O3 draw-real8-t.f -o draw -L/usr/local/lib -lpgplot -fbackslash -ffixed-line-length-132
c
c-----------------------------------------------------------------------
  
      Subroutine Output_Help
      
      Write(6,*) " _____                      _____ ____  _   _ ______ " 
      Write(6,*) "|  __ \                    / ____/ __ \| \ | |  ____|"
      Write(6,*) "| |  | |_ __ __ ___      _| |   | |  | |  \| | |__   "
      Write(6,*) "| |  | | '__/ _` \ \ /\ / / |   | |  | | . ` |  __|  "
      Write(6,*) "| |__| | | | (_| |\ V  V /| |___| |__| | |\  | |     "
      Write(6,*) "|_____/|_|  \__,_| \_/\_/  \_____\____/|_| \_|_|     "
      Write(6,*) "  " 
      Write(6,*) " Usages :"
      Write(6,*) "-------------" 
      Write(6,*) " " 
      Write(6,*) " CASE 1 : ./draw 34" 
      Write(6,*) "   --> will open a graphical windows" 
      Write(6,*) "       load the file named conf34"
      Write(6,*) "       use + to go to next conf file"
      Write(6,*) "       use - to go to previous conf file"
      Write(6,*) " "   
      Write(6,*) " CASE 2 : ./draw 34 45" 
      Write(6,*) "   --> will open a graphical windows" 
      Write(6,*) "       load the files named conf34 and conf45"
      Write(6,*) "       use + to go to next conf file"
      Write(6,*) "       use - to go to previous conf file"
      Write(6,*) " "
      Write(6,*) " CASE 3 : ./draw 34 45 10" 
      Write(6,*) "   --> will open a graphical windows" 
      Write(6,*) "       load the files named conf34 and conf45"
      Write(6,*) "       use + to add 10 to conf numbers"
      Write(6,*) "       use - to remove 10 to conf numbers"
      Write(6,*) " "
      Write(6,*) " WITH the mouse pointer in the window :"
      Write(6,*) "       use G or g to plot the Grains" 
      Write(6,*) "       use C or c to plot Contact forces" 
      Write(6,*) "       left click to zoom in"
      Write(6,*) "       right click to zoom out" 
      Write(6,*) "  --> with Case 2 or 3 :"
      Write(6,*) "       use D or d to plot the grain Displacements" 
      Write(6,*) "       use F or f to plot the grain Fluctuation" 
      Write(6,*) "       use T or t to plot a Delaunay triangulation"
      Write(6,*) "       use S or s to plot local shear strains"
      Write(6,*) " " 
      Write(6,*) " ESC to exit draw" 
      Stop

      End
c=======================================================================

      Program draw_main
      
      Implicit None      

      Character*50 FichierConf(2)
      Integer iCharLen5
      Integer ierr_read
c      Integer i
      Integer npa,npam
      Parameter (npam=100000)
      Real*8 r1(3,npam)
      Real*8 r2(3,npam)
      Integer ntotal,ntotalm
      Parameter(ntotalm=50*npam)
      Integer ior(ntotalm),iex(ntotalm)
      Real*8 hxx,hxy,hyx,hyy
      Real*8 hxx2,hxy2,hyx2,hyy2
      Real*8 fn(ntotalm),ft(ntotalm)
c     -- decla des arguments
      Character*5 arg,cc,ch_arg(3)
      Integer nc
      Integer i,iarg,iconf(3)
c     -- pgplot      
      Real*4 xmax
      Character ch
      Real*4 posx,posy
      integer ier,pgcurs
      Integer PGOPEN
      Integer idisc,iforce,idep,ifluct,itriangle,istrain
      Real*4 x1,x2,y1,y2
      Real*8 CoinsIni(2,4),CoinsFin(2,4)

c     -- lecture des arguments 
      If (iargc().EQ.0) Call Output_Help
      Do iarg = 1,iargc()
         Call getarg(iarg,arg)
         ch_arg(iarg) = arg
         Read(arg(1:5),'(i5)') iconf(iarg) ! convert string to integer
         Write(6,*) iconf(iarg)
      EndDo
      iarg = iarg - 1 !En sortie du Do, le nbre d argument est de 1 trop grand

      Call ModifString(ch_arg(1))
      FichierConf(1) = 'CONF'//ch_arg(1)
      If (iarg.GE.2) Then
         Call ModifString(ch_arg(2))
         FichierConf(2) = 'CONF'//ch_arg(2)
      EndIf 

      If (PGOPEN('/xwin').LE.0) Then 
         Write(6,*) " pb to open xwindow" 
         Stop
      EndIf 
                    
      Call Read_Conf_File(FichierConf(1),hxx,hxy,hyx,hyy,npa,r1,CoinsIni,ntotal,ior,iex,fn,ft,ierr_read)
      If (ierr_read.EQ.0) Stop ! le fichier conf n existe pas

c     -- PGPLOT drawing      
c      Call PGSCRN(0,'DarkSalmon',ier)
      Call PGSCRN(0,'navy',ier)

c     -- lecture positions et rayon des particules
      xmax = -1e20
      Do i = 1,npa
         xmax = max(xmax,r1(1,i)+r1(3,i))
      EndDo
      xmax = xmax*1.1
      Call PGENV(-xmax,xmax,-xmax,xmax,1,-2) ! axis system definition

c --- dessin avec un seul fichier conf
      If (iarg.EQ.1) Then 
         idisc = 1
         iforce = 0
1        Call PGERAS
         If (idisc.EQ.1) Call Plot_Disc(npa,hxx,hxy,hyx,hyy,r1)
         If (iforce.EQ.1) Call Plot_Force(hxx,hxy,hyx,hyy,npa,r1,ntotal,ior,iex,fn)
         ier = PGCURS(posx,posy,ch)
c        -- Pour sortir (mais pourquoi sortir ? )
         If (ch.EQ.char(27)) Stop
c        -- dessin des grains
         If (ch.EQ."G".OR.ch.EQ."g") Then
            If (idisc.EQ.0) Then 
                idisc = 1 
            Else
                idisc = 0
            EndIf
         EndIf
c        -- dessin des forces
         If (ch.EQ."C".OR.ch.EQ."c") Then
            If (iforce.EQ.0) Then 
                iforce = 1 
            Else
                iforce = 0
            EndIf
         EndIf
c        -- passage au fichier conf suivant
         If (ch.EQ."+") Then
            iconf(1) = iconf(1)+1
            Write(arg,'(i5)') iconf(1)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(1) = 'CONF'//cc
            Write(6,*) FichierConf(1)
            Call Read_Conf_File(FichierConf(1),hxx,hxy,hyx,hyy,npa,r1,CoinsIni,ntotal,ior,iex,fn,ft,ierr_read)
         EndIf 
c        -- passage au fichier conf precedent 
         If (ch.EQ."-") Then
            iconf(1) = iconf(1)-1
            If (iconf(1).LT.0) iconf(1) = 0
            Write(arg,'(i5)') iconf(1)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(1) = 'CONF'//cc
            Write(6,*) FichierConf(1)
            Call Read_Conf_File(FichierConf(1),hxx,hxy,hyx,hyy,npa,r1,CoinsIni,ntotal,ior,iex,fn,ft,ierr_read)
         EndIf  
c        -- zoom out : A = click gauche
         If (ch.EQ."A") Then
            xmax = xmax * 0.9
            Call PGSWIN(-xmax,xmax,-xmax,xmax)
c        -- zoom in : X = click droit
         ElseIf (ch.EQ."X") Then
            xmax = xmax * 1.1 
            Call PGSWIN(-xmax,xmax,-xmax,xmax)
         EndIf
         GoTo 1
      Else
c --- dessin avec deux fichiers conf
c        -- lecture du second fichier conf
         Call Read_Conf_File(FichierConf(2),hxx2,hxy2,hyx2,hyy2,npa,r2,CoinsFin,ntotal,ior,iex,fn,ft,ierr_read)
         If (ierr_read.EQ.0) Then
            Write(6,*) "no second conf file"
            Stop
         EndIf
         idisc = 1
         iforce = 0
         idep = 0
         ifluct = 0
         itriangle = 0
         istrain = 0
2        Call PGERAS
         If (idisc.EQ.1) Call Plot_Disc(npa,hxx,hxy,hyx,hyy,r1)
         If (iforce.EQ.1) Call Plot_Force(hxx,hxy,hyx,hyy,npa,r1,ntotal,ior,iex,fn)
c         If (ifluct.EQ.1) Call Plot_Fluct(npa,r1,r2,hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2)
         If (ifluct.EQ.1) Call Plot_Fluct2(npa,r1,r2,coinsIni,coinsFin)
         If (idep.EQ.1) Call Plot_Dep(npa,r1,r2,hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2)
         If (itriangle.EQ.1) Call Plot_Triangle(npa,r1,hxx,hxy,hyx,hyy)
         If (istrain.EQ.1) Call Plot_Deformation(npa,r1,r2,hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2)
         ier = PGCURS(posx,posy,ch)
c        -- Pour sortir (mais pourquoi sortir ? )
         If (ch.EQ.char(27)) Stop
c        -- dessin des grains
         If (ch.EQ."G".OR.ch.EQ."g") Then
            If (idisc.EQ.0) Then 
                idisc = 1 
            Else
                idisc = 0
            EndIf
         EndIf
c        -- dessin des deplacements
         If (ch.EQ."D".OR.ch.EQ."d") Then
            If (idep.EQ.0) Then 
                idep = 1 
            Else
                idep = 0
            EndIf
         EndIf
c        -- dessin des forces
         If (ch.EQ."C".OR.ch.EQ."c") Then
            If (iforce.EQ.0) Then 
                iforce = 1 
            Else
                iforce = 0
            EndIf
         EndIf
c        -- dessin des fluctuations
         If (ch.EQ."F".OR.ch.EQ."f") Then
            If (ifluct.EQ.0) Then 
                ifluct = 1 
            Else
                ifluct = 0
            EndIf
         EndIf
c        -- dessin des triangle
         If (ch.EQ."T".OR.ch.EQ."t") Then
            If (itriangle.EQ.0) Then 
                itriangle = 1 
            Else
                itriangle = 0
            EndIf
         EndIf
c        -- dessin des deformation
         If (ch.EQ."S".OR.ch.EQ."s") Then
            If (istrain.EQ.0) Then 
                istrain = 1 
            Else
                istrain = 0
            EndIf
         EndIf
c        -- passage au couple de fichiers conf suivant
         If (ch.EQ."+") Then
            If (iarg.EQ.3) Then
               iconf(1) = iconf(1) + iconf(3)
            Else 
               iconf(1) = iconf(1)+1
            EndIf
            Write(arg,'(i5)') iconf(1)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(1) = 'CONF'//cc
            Write(6,*) FichierConf(1)
            Call Read_Conf_File(FichierConf(1),hxx,hxy,hyx,hyy,npa,r1,CoinsIni,ntotal,ior,iex,fn,ft,ierr_read)
            If (iarg.EQ.3) Then
               iconf(2) = iconf(2) + iconf(3)
            Else 
               iconf(2) = iconf(2)+1
            EndIf
            Write(arg,'(i5)') iconf(2)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(2) = 'CONF'//cc
            Write(6,*) FichierConf(2)
            Call Read_Conf_File(FichierConf(2),hxx2,hxy2,hyx2,hyy2,npa,r2,CoinsFin,ntotal,ior,iex,fn,ft,ierr_read)
         EndIf 
c        -- passage au couple de fichiers conf precedent
         If (ch.EQ."-") Then
            If (iarg.EQ.3) Then
               iconf(1) = iconf(1) - iconf(3)
            Else 
               iconf(1) = iconf(1)- 1
            EndIf
            Write(arg,'(i5)') iconf(1)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(1) = 'CONF'//cc
            Write(6,*) FichierConf(1)
            Call Read_Conf_File(FichierConf(1),hxx,hxy,hyx,hyy,npa,r1,CoinsIni,ntotal,ior,iex,fn,ft,ierr_read)
            If (iarg.EQ.3) Then
               iconf(2) = iconf(2) - iconf(3)
            Else 
               iconf(2) = iconf(2) - 1
            EndIf
            Write(arg,'(i5)') iconf(2)
            nc = iCharLen5(arg) - 1
            cc = arg(5-nc:5)
            Call ModifString(cc)
            FichierConf(2) = 'CONF'//cc
            Write(6,*) FichierConf(2)
            Call Read_Conf_File(FichierConf(2),hxx2,hxy2,hyx2,hyy2,npa,r2,CoinsFin,ntotal,ior,iex,fn,ft,ierr_read)
         EndIf 
c        -- zoom out : A = click gauche
         If (ch.EQ."A") Then
            xmax = xmax * 0.9
            Call PGSWIN(-xmax,xmax,-xmax,xmax)
c        -- zoom in : X = click droit
         ElseIf (ch.EQ."X") Then
            xmax = xmax * 1.1 
            Call PGSWIN(-xmax,xmax,-xmax,xmax)
         EndIf
         GoTo 2
      EndIf
              

      CALL PGCLOS      
     
      
      End
c=======================================================================
      Subroutine Plot_Force(hxx,hxy,hyx,hyy,npa,r,ntotal,ior,iex,fn)
      
      Implicit None

      Integer i,j,npa
      Integer il,ntotal
      Real*8 hxx,hxy,hyx,hyy
      Integer npam
      Parameter (npam=100000)
      Real*8 r(3,npam)
      Integer ntotalm
      Parameter(ntotalm=40*npam)
      Integer ior(ntotalm),iex(ntotalm)
      Real*8 fn(ntotalm),fn2(ntotalm)
      Real*8 s(2,npam)     
      Real*8 den
      Real*8 fnmoy,fnmin,fnmax
      Real*8 sxij,syij
      Real*8 a,b
c     -- pgplot
      Integer width
      Real*4 px(2),py(2)
      
      If (ntotal.LT.1) Then
         Write(6,*) " No forces" 
         Return
      EndIf
      
      den = hxx*hyy-hxy*hyx      
      Do i = 1,npa
         s(1,i) = (hyy*r(1,i)-hxy*r(2,i))/den
         s(2,i) = (hxx*r(2,i)-hyx*r(1,i))/den
      EndDo
      fnmoy = 0.
      Do il = 1,ntotal
         fnmoy = fnmoy + fn(il)
      EndDo
      fnmoy = fnmoy / float(ntotal)
      Do il = 1,ntotal
         fn2(il) = fn(il) / fnmoy
      EndDo
      Write(6,*) 'ncont = ',ntotal
      Write(6,*) 'Mean force = ',fnmoy      
      fnmax = 0.
      fnmin = 100.
      Do il = 1,ntotal
         fnmax = max(fnmax,fn2(il))
         fnmin = min(fnmin,fn2(il))
      EndDo
      Write(6,*) 'fnmin, fnmax = ',fnmin,fnmax  
      If (fnmax.EQ.fnmin) Then
         a = 1
         b = 0
      Else  
         a = (20.-1.) / (fnmax-fnmin)
         b = 1. - a*fnmin
      EndIf
      Call PGSCI(2)  
      Do il = 1,ntotal
         i = ior(il)
         j = iex(il)
         sxij = s(1,j)-s(1,i)
         sxij = sxij - int(sxij+sxij)
         syij = s(2,j)-s(2,i)
         syij = syij - int(syij+syij)
         px(1) = (s(1,j)-sxij)*hxx + hxy*(s(2,j)-syij)
         px(2) = s(1,j)*hxx + hxy*s(2,j)
         py(1) = (s(2,j)-syij)*hyy + hyx*(s(1,j)-sxij)
         py(2) = s(2,j)*hyy + s(1,j)*hyx
         width = int(a * fn2(il) + b)
         Call PGSLW(width)
         Call PGLINE(2,px,py)
      EndDo
            
      End 
c=======================================================================
      
      Subroutine Plot_Disc(npa,hxx,hxy,hyx,hyy,r)
      
      Implicit None
      
      Integer npam
      Parameter (npam=100000)
      Integer npa,i
      Integer m,mm,n,nn
      Real*8 r(3,npam),s(2,npam)
      Real*8 hxx,hxy,hyx,hyy,den
      Real*8 xxx,yyy
c     -- pgplot
      Real*4 rx,ry,xray 
      
      den = hxx*hyy-hxy*hyx      
      Do i = 1,npa
         s(1,i) = (hyy*r(1,i)-hxy*r(2,i))/den
         s(2,i) = (hxx*r(2,i)-hyx*r(1,i))/den
      EndDo


      Call PGSLW(1)

      CALL PGSFS(1)
      CALL PGSCI(15)
      Do i = 1,npa
         xxx = s(1,i)
         yyy = s(2,i)
         rx = hxx*xxx + hxy*yyy
         ry = hyx*xxx + hyy*yyy
         xray = r(3,i)
         Call PGCIRC(rx,ry,xray)
      EndDo

      CALL PGSFS(2)
      CALL PGSCI(1)
      Do i = 1,npa
         xxx = s(1,i)
         yyy = s(2,i)
         rx = hxx*xxx + hxy*yyy
         ry = hyx*xxx + hyy*yyy
         xray = r(3,i)
         Call PGCIRC(rx,ry,xray)
      EndDo

      End  

c=======================================================================
      
      Subroutine Plot_Dep(npa,r,r2,hxxIni,hxyIni,hyxIni,hyyIni,hxxFin,hxyFin,hyxFin,hyyFin)
      
      Implicit None
      
      Integer npam
      Parameter (npam=100000)
      Integer i,k
      Integer npa
      Real*8 r(3,npam)
      Real*8 hxxIni,hxyIni,hyxIni,hyyIni
      Real*8 r2(3,npam)
      Real*8 hxxFin,hxyFin,hyxFin,hyyFin
      real*8 dhxx,dhyy,dhxy,dhyx
      Real*8 xIni(npam),yIni(npam)
      Real*8 xFin(npam),yFin(npa)
      Real*8 sxIni(npam),syIni(npam)
      Real*8 sxFin(npam),syFin(npam)
      Real*8 depX(npam),depY(npam),dep(npam)
      Real*8 dmax,dmin,den
      Real*8 dsxx,dsyy
      Real*8 xmin,xmax,ymin,ymax
      Real*8 fact,a,b
c     -- pgplot
      Real*4 fleche
      Real*4 x1(2,npam),x2(2,npam)

      Do k = 1,npa
         xIni(k) = r(1,k)
         yIni(k) = r(2,k)
      EndDo
      Do k = 1,npa
         xFin(k) = r2(1,k)
         yFin(k) = r2(2,k)
      EndDo
      dmax = -1
      dmin = 1e10
      dhxx = hxxFin-hxxIni
      dhyy = hyyFin-hyyIni
      dhxy = hxyFin-hxyIni
      dhyx = hyxFin-hyxIni
      Do k = 1,npa
         dsxx = xFin(k) - xIni(k)
         dsyy = yFin(k) - yIni(k)
         depX(k) = dsxx 
         depY(k) = dsyy
         dep(k) = sqrt(depX(k)*depX(k) + depY(k)*depY(k))
         dmax = max(dmax,dep(k))
         dmin = min(dmin,dep(k))
      EndDo
      Write(6,*) 'Deplacement max = ',dmax
      If (dmax.EQ.0.) Return

      xmin = 1e7
      xmax = -1e7
      ymin = 1e7
      ymax = -1e7
      Do i = 1,npa
         x1(1,i) = xIni(i)
         x1(2,i) = yIni(i)
         x2(1,i) = xIni(i)+depX(i)
         x2(2,i) = yIni(i)+depY(i)
         xmin = min(xmin,x1(1,i),x2(1,i))
         xmax = max(xmax,x1(1,i),x2(1,i))
         ymin = min(ymin,x1(2,i),x2(2,i))
         ymax = max(ymax,x1(2,i),x2(2,i))
      EndDo
c      Write(6,*) xmin,ymin,xmax,ymax
      xmin = 1.05*xmin
      xmax = 1.05*xmax
      ymin = 1.05*ymin
      ymax = 1.05*ymax

      fact = abs(xmax-xmin)/5/dmax
      
c      Write(6,*) 'fact = ',fact
      Do i = 1,npa
         x2(1,i) = x1(1,i)+depX(i)*fact
         x2(2,i) = x1(2,i)+depY(i)*fact
      EndDo
      Call PGSCI(6)
      CALL PGSAH(1.,45.,0.5)
      Do i = 1,npa
         a = (1.-0.01)/(dmax-dmin)
         b = 0.01 - a*dmin
         fleche = a*dep(i)+b
         CALL PGSCH(fleche) !taille de la tete de fleche
         CALL PGARRO(x1(1,i),x1(2,i),x2(1,i),x2(2,i))
      EndDo

      End

c=======================================================================
      Subroutine Plot_Fluct(npa,r1,r2,hxxIni,hxyIni,hyxIni,hyyIni,hxxFin,hxyFin,hyxFin,hyyFin)
      
      Implicit None
      
      Integer npam
      Parameter (npam=100000)
      Integer k
      Integer npa
      Real*8 r1(3,npam)
      Real*8 hxxIni,hxyIni,hyxIni,hyyIni
      Real*8 r2(3,npam)
      Real*8 hxxFin,hxyFin,hyxFin,hyyFin
      Real*8 xIni(npam),yIni(npam)
      Real*8 xFin(npam),yFin(npam)
      Real*8 sxIni(npam),syIni(npam)
      Real*8 sxFin(npam),syFin(npam)
      Real*8 depX(npam),depY(npam),dep(npam)
      Real*8 dmax,dmin,den
      Real*8 dsxx,dsyy
      Real*8 xmin,xmax,ymin,ymax
      Real*8 fact,a,b
c     -- pgplot
      Real*4 x1(2,npam),x2(2,npam)
      Real*4 fleche

      den = hxxIni*hyyIni - hyxIni*hxyIni
      Do k = 1,npa
         xIni(k) = r1(1,k)
         yIni(k) = r1(2,k)
         sxIni(k) = (hyyIni*xIni(k) - hxyIni*yIni(k))/den
         syIni(k) = (hxxIni*yIni(k) - hyxIni*xIni(k))/den
      EndDo
      den = hxxFin*hyyFin - hyxFin*hxyFin
      Do k = 1,npa
         xFin(k) = r2(1,k)
         yFin(k) = r2(2,k)
         sxFin(k) = (hyyFin*xFin(k) - hxyFin*yFin(k))/den
         syFin(k) = (hxxFin*yFin(k) - hyxFin*xFin(k))/den 
      EndDo
      dmax = -1
      dmin = 1e10
      Do k = 1,npa
         dsxx = sxFin(k) - sxIni(k)
         dsyy = syFin(k) - syIni(k)
         dsxx = dsxx - int(dsxx+dsxx)
         dsyy = dsyy - int(dsyy+dsyy)
         depX(k) = dsxx * hxxFin + dsyy * hxyFin
         depY(k) = dsxx * hyxFin + dsyy * hyyFin
         dep(k) = sqrt(depX(k)*depX(k) + depY(k)*depY(k))
         dmax = max(dmax,dep(k))
         dmin = min(dmin,dep(k))
      EndDo
      Write(6,*) 'Fluctuation max = ',dmax
      If (dmax.EQ.0.) Return

      xmin = 1e7
      xmax = -1e7
      ymin = 1e7
      ymax = -1e7
      Do k = 1,npa
         x1(1,k) = r1(1,k)
         x1(2,k) = r1(2,k)
         x2(1,k) = x1(1,k)+depX(k)
         x2(2,k) = x1(2,k)+depY(k)
         xmin = min(xmin,x1(1,k),x2(1,k))
         xmax = max(xmax,x1(1,k),x2(1,k))
         ymin = min(ymin,x1(2,k),x2(2,k))
         ymax = max(ymax,x1(2,k),x2(2,k))
      EndDo
c      Write(6,*) xmin,ymin,xmax,ymax
      xmin = 1.05*xmin
      xmax = 1.05*xmax
      ymin = 1.05*ymin
      ymax = 1.05*ymax

      fact = abs(xmax-xmin)/5/dmax
      Call PGSCI(7)
      CALL PGSAH(1.,45.,0.5)
      a = (1.-0.01)/(dmax-dmin)
      b = 0.01 - a*dmin
      Do k = 1,npa
         fleche = a*dep(k)+b
         CALL PGSCH(fleche) !taille de la tete de fleche
         x2(1,k) = x1(1,k)+depX(k)*fact
         x2(2,k) = x1(2,k)+depY(k)*fact
         CALL PGARRO(x1(1,k),x1(2,k),x2(1,k),x2(2,k))
      EndDo

      End

c=======================================================================

      Subroutine Plot_Fluct2(npa,r1,r2,coinsIni,coinsFin)
      
      Implicit None
      
      Integer npam
      Parameter (npam=100000)
      Integer k
      Integer npa
      Real*8 coinsIni(2,4),coinsFin(2,4)
      Real*8 XX(2),UU(2)
      Real*8 r1(3,npam)
      Real*8 r2(3,npam)
      Real*8 xIni(npam),yIni(npam)
      Real*8 xFin(npam),yFin(npam)
      Real*8 sxIni(npam),syIni(npam)
      Real*8 sxFin(npam),syFin(npam)
      Real*8 depX(npam),depY(npam),dep(npam)
      Real*8 dmax,dmin,den
      Real*8 dsxx,dsyy
      Real*8 xmin,xmax,ymin,ymax
      Real*8 fact,a,b
c     -- pgplot
      Real*4 x1(2,npam),x2(2,npam)
      Real*4 fleche

      dmax = 0.
      dmin = 1d10
      Do k = 1,npa
         XX(1) = r1(1,k)
         XX(2) = r1(2,k)
         Call MeanField(coinsIni,coinsFin,XX,UU)
         depX(k) = (r2(1,k)-r1(1,k)) - UU(1)
         depY(k) = (r2(2,k)-r1(2,k)) - UU(2)
         dep(k) = sqrt(depX(k)*depX(k) + depY(k)*depY(k))
         dmax = max(dmax,dep(k))
         dmin = min(dmin,dep(k))
      EndDo

      Write(6,*) 'Fluctuation max = ',dmax
      If (dmax.EQ.0.) Return

      xmin = 1e7
      xmax = -1e7
      ymin = 1e7
      ymax = -1e7
      Do k = 1,npa
         x1(1,k) = r1(1,k)
         x1(2,k) = r1(2,k)
         x2(1,k) = x1(1,k)+depX(k)
         x2(2,k) = x1(2,k)+depY(k)
         xmin = min(xmin,x1(1,k),x2(1,k))
         xmax = max(xmax,x1(1,k),x2(1,k))
         ymin = min(ymin,x1(2,k),x2(2,k))
         ymax = max(ymax,x1(2,k),x2(2,k))
      EndDo
c      Write(6,*) xmin,ymin,xmax,ymax
      xmin = 1.05*xmin
      xmax = 1.05*xmax
      ymin = 1.05*ymin
      ymax = 1.05*ymax

      fact = abs(xmax-xmin)/5/dmax
      Call PGSCI(7)
      CALL PGSAH(1.,45.,0.5)
      a = (1.-0.01)/(dmax-dmin)
      b = 0.01 - a*dmin
      Do k = 1,npa
         fleche = a*dep(k)+b
         CALL PGSCH(fleche) !taille de la tete de fleche
         x2(1,k) = x1(1,k)+depX(k)*fact
         x2(2,k) = x1(2,k)+depY(k)*fact
         CALL PGARRO(x1(1,k),x1(2,k),x2(1,k),x2(2,k))
      EndDo

      End

c=======================================================================
      Subroutine Read_Conf_File(filename,hxx,hxy,hyx,hyy,npa,r,coin,ntotal,ior,iex,fn,ft,ierr_read)
      
      Implicit None      

      Integer npam
      Parameter (npam=100000)
      Character*50 filename
      Integer i,j,npa
      Real*8 z
      Real*8 hxx,hxy,hyx,hyy,x1,x2,x3,x4,y1,y2,y3,y4
      Real*8 r(3,npam),xmin,xmax,ymin,ymax,xdep,ydep
      Real*8 xij,yij,hij
      Integer ntotal
      Real*8 fn(ntotal),ft(ntotal)
      Integer ior(npam),iex(npam)
      Logical lval
      Integer ierr_read
      Real*8 coin(2,4)
      
c------Lecture du fichier conf   
      ierr_read = 1
      Inquire(file=filename, exist=lval)  
      If (.NOT.lval) Then
         Write(6,*) " None existing file: ",filename
         ierr_read = 0
         Return
      EndIf
 
      Write(6,*) "Reading ",filename

c     -- ouverture du fichier CONF
      Open(1,file=filename,status='old')
c     -- lecture du nbre de particules
      Read(1,*) npa
c     -- lecture positions et rayon des particules
      xmin = 1e20
      xmax = -1e20
      ymin = 1e20
      ymax = -1e20
      Do i = 1,npa
         Read(1,*) r(1,i),r(2,i),r(3,i)
         xmin = min(xmin,r(1,i)-r(3,i))
         xmax = max(xmax,r(1,i)+r(3,i))
         ymin = min(ymin,r(2,i)-r(3,i))
         ymax = max(ymax,r(2,i)+r(3,i))
      EndDo
c     -- on recentre l'echantillon
      xdep = (xmax+xmin)*0.5d0
      ydep = (ymax+ymin)*0.5d0
      Write(6,*) xdep,ydep
      Do i = 1,npa
         r(1,i) = r(1,i) - xdep
         r(2,i) = r(2,i) - ydep
      EndDo
c     -- on saute une ligne
      Read(1,*) z
c     -- lecture des coord. des 4 coins
      Read(1,*) x1,y1 ! en bas a  droite
      Read(1,*) x2,y2 ! en haut a droite, pas necessaire
      Read(1,*) x3,y3 ! en haut a  gauche
      Read(1,*) x4,y4 ! en bas a  gauche
      coin(1,1) = x1 - xdep
      coin(2,1) = y1 - ydep
      coin(1,2) = x2 - xdep
      coin(2,2) = y2 - ydep
      coin(1,3) = x3 - xdep
      coin(2,3) = y3 - ydep
      coin(1,4) = x4 - xdep
      coin(2,4) = y4 - ydep

      Close(1)
c     -- calcul de la matrice h equivalente
      hxx = x1-x4
      hyx = y1-y4
      hxy = x3-x4
      hyy = y3-y4  
      
c     -- les interactions
      ntotal = 0
      Do i = 1,npa
         Do j = i+1,npa
            xij = r(1,j) - r(1,i)
            yij = r(2,j) - r(2,i)
            hij = sqrt(xij*xij + yij*yij) - r(3,i) - r(3,j)
            If (hij.LE.0) Then
               ntotal = ntotal + 1
               ior(ntotal) = i
               iex(ntotal) = j
               fn(ntotal) = 1. ! force fictive
               ft(ntotal) = 0. ! force fictive
            EndIf
         EndDo
      EndDo
      Write(6,*) "npa, ntotal = ", npa,ntotal

      End
c=======================================================================      

      Integer Function iCharLen5(text)

      Integer i
      Character*(*) text
      Character*1 cc

      i = 5 + 1
1     i = i - 1
      cc = text(i:i)
      If (cc.NE.'') GoTo 1

      iCharLen5 = 5 - i

      Return
      End

c======================================================================      

      Subroutine Plot_Triangle(npa,r1,hxx,hxy,hyx,hyy)

      Implicit None

      Integer npam
      Parameter (npam=100000)
      Integer i,j,k,nt,npa
      Integer ivnxt
      Real*8 r1(3,npam),r2(3,npam)
      Real*8 s(2,npam)
      Real*8 s_fake(2,npam)
      Real*8 hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2
      Real*8 den,z
c     -- triangles
      Integer lt(3,npam*5)
c     -- pgplot
      Real*4 px(3),py(3)

      den = hxx*hyy-hxy*hyx      
      Do i = 1,npa
         s(1,i) = (hyy*r1(1,i)-hxy*r1(2,i))/den
         s(2,i) = (hxx*r1(2,i)-hyx*r1(1,i))/den
      EndDo

      Call Triangule(npa,s,ivnxt,lt) 

      den = hxx*hyy-hxy*hyx      
      Do i = 1,npa
         s(1,i) = (hyy*r1(1,i)-hxy*r1(2,i))/den
         s(2,i) = (hxx*r1(2,i)-hyx*r1(1,i))/den
      EndDo
      Call RetireTrianglesPlats(npa,s,ivnxt,lt)        

      Write(6,*) "Nbre triangle = ",ivnxt
      If (ivnxt.EQ.0) Return

c     -- dessin pgplot
      Call PGSCI(2)  
      Call PGSLW(1)
      Do nt = 1,ivnxt
         i = lt(1,nt)
         j = lt(2,nt)   
         k = lt(3,nt)   
         px(1) = r1(1,i)
         px(2) = r1(1,j)
         px(3) = r1(1,k)
         py(1) = r1(2,i)
         py(2) = r1(2,j)
         py(3) = r1(2,k)
         Call PGLINE(3,px,py)
      EndDo
 
      End

c=======================================================================

      Subroutine Plot_Deformation(npa,r1,r2,hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2)

      Implicit None

      Integer npam
      Parameter (npam=100000)
      Integer i,j,k,nt,npa
      Integer i01,i02,i03
      Integer i1,i2,i3
      Integer ivnxt
      Real*8 r1(3,npam),r2(3,npam)
      Real*8 s(2,npam)
      Real*8 hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2
      Real*8 den
c     -- triangles
      Integer lt(3,npam*5)
c     -- Deformation
      Real*8 strain(npam*5)
      Real*8 dmin,dmax,dmoy
      real*8 Moment
c     -- pgplot
      Real*4 px(3),py(3),size,xx,yy
      real*4 a,b

      den = hxx*hyy-hxy*hyx      
      Do i = 1,npa
         s(1,i) = (hyy*r1(1,i)-hxy*r1(2,i))/den
         s(2,i) = (hxx*r1(2,i)-hyx*r1(1,i))/den
      EndDo

      Call Triangule(npa,s,ivnxt,lt)
      Call RetireTrianglesPlats(npa,s,ivnxt,lt)
c     -- les sommets sont mis dans le sens inverse trigo
      Do nt = 1,ivnxt
         i01 = lt(1,nt)
         i02 = lt(2,nt)
         i03 = lt(3,nt)
         Call Renum(i01,i02,i03,s,i1,i2,i3)
         lt(1,nt) = i1
         lt(2,nt) = i2
         lt(3,nt) = i3
      EndDo
         
      Write(6,*) "Nbre triangle = ",ivnxt
      If (ivnxt.EQ.0) Return

      Call CalculDeformation(npa,r1,r2,hxx,hxy,hyx,hyy,hxx2,hxy2,hyx2,hyy2,ivnxt,lt,strain)
      dmin = Moment(strain,ivnxt,1)
      dmax = Moment(strain,ivnxt,2)
      dmoy = Moment(strain,ivnxt,3)
      Write(6,*) "Deformatio min,max,moy ",dmin,dmax,dmoy
      a = (5.-0.1) / (dmax)!-dmin)
      b = 0.1 - a*dmin
      Call PGSCI(1)
      Do i = 1,ivnxt
         i1 = lt(1,i)
         i2 = lt(2,i)
         i3 = lt(3,i)
         size = a*strain(i)+b
         CALL PGSCH(size)
         px(1) = r1(1,i1) 
         py(1) = r1(2,i1)
         px(2) = r1(1,i2)
         py(2) = r1(2,i2)
         px(3) = r1(1,i3)
         py(3) = r1(2,i3)
         xx = (px(1)+px(2)+px(3))/3.
         yy = (py(1)+py(2)+py(3))/3.
         CALL PGPT1(xx,yy,19)
      EndDo
 
      End

c=======================================================================
      Real*8 Function Moment(data,n,im)

      Implicit None 

      Real*8 data(*)
      Integer i,n
      Integer im ! 1 = valeur minimum
      Real*8 datamin,datamax,datamoy

      If (im.EQ.1) Then
         datamin = 1e10
         Do i = 1,n
            datamin = min(datamin,data(i))
         EndDo
         Moment = datamin
      Else If (im.EQ.2) Then
         datamax = -1e10
         Do i = 1,n
            datamax = max(datamax,data(i))
         EndDo
         Moment = datamax
      Else If (im.EQ.3) Then
         datamoy = 0.d0
         Do i = 1,n
            datamoy = datamoy + data(i)
         EndDo
         datamoy = datamoy / dfloat(n)
         Moment = datamoy
      EndIf

      Return
      End

c=======================================================================

      Subroutine CalculDeformation(npa,r1,r2,hxx1,hxy1,hyx1,hyy1,hxx2,hxy2,hyx2,hyy2,ivnxt,lt,strain)

      Implicit None

      Integer npam
      Parameter (npam=100000)
      Integer i,j,k,nt,npa
      Integer i1,i2,i3
      Integer ivnxt
      Real*8 r1(3,npam),r2(3,npam)
      Real*8 s1(2,npam),s2(2,npam)
      Real*8 hxx1,hxy1,hyx1,hyy1,hxx2,hxy2,hyx2,hyy2
      Real*8 dhxx,dhxy,dhyx,dhyy
      Real*8 den1,den2
c     -- triangles
      Integer lt(3,npam*5)
c     -- Deformation
      Real*8 u(2,npam)
      Real*8 strain(npam*5)
      Real*8 Surf,SS
      Real*8 dsx,dx,dxx,dsy,dy,dyy
      Real*8 axx,axy,ayx,ayy
      Real*8 nx,ny,Ux,Uy
      Real*8 Exx,Eyy,Exy,Ew,Delta,E1,E2,Dis,Vol

c     -- coordonnees reduites
      den1 = hxx1*hyy1-hxy1*hyx1      
      Do i = 1,npa
         s1(1,i) = (hyy1*r1(1,i)-hxy1*r1(2,i))/den1
         s1(2,i) = (hxx1*r1(2,i)-hyx1*r1(1,i))/den1
      EndDo
      den2 = hxx2*hyy2-hxy2*hyx2      
      Do i = 1,npa
         s2(1,i) = (hyy2*r2(1,i)-hxy2*r2(2,i))/den2
         s2(2,i) = (hxx2*r2(2,i)-hyx2*r2(1,i))/den2
      EndDo
c     -- increment de H
      dhxx = hxx2 - hxx1
      dhxy = hxy2 - hxy1
      dhyx = hyx2 - hyx1
      dhyy = hyy2 - hyy1  
c     -- calcul des deplacements
      Do i = 1,npa
         dsx = s2(1,i)-s1(1,i)
         dsx = dsx - int(dsx+dsx)
         dsy = s2(2,i)-s1(2,i)
         dsy = dsy - int(dsy+dsy)
         u(1,i) = hxx1*dsx+hxy1*dsy + dhxx*s2(1,i)+dhxy*s2(2,i)
         u(2,i) = hyx1*dsx+hyy1*dsy + dhyx*s2(1,i)+dhyy*s2(2,i)
      EndDo
c     -- calcul des deformations
      Do nt = 1,ivnxt
         i1 = lt(1,nt)
         i2 = lt(2,nt)
         i3 = lt(3,nt)
c         --Calcul des deformations locales
         axx = 0.d0
         ayy = 0.d0
         ayx = 0.d0
         axy = 0.d0
c
         dx = s1(1,i2)-s1(1,i1)
	 dx = dx - int(dx+dx)
         dy = s1(2,i2)-s1(2,i1)
	 dy = dy - int(dy+dy)
         dxx = hxx1*dx+hxy1*dy
         dyy = hyx1*dx+hyy1*dy
         nx = - dyy  ! normale sortante
         ny = dxx 
         Ux = (u(1,i2)+u(1,i1))*0.5
         Uy = (u(2,i2)+u(2,i1))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy
c
         dx = s1(1,i3)-s1(1,i2)
	 dx = dx - int(dx+dx)
         dy = s1(2,i3)-s1(2,i2)
	 dy = dy - int(dy+dy)
         dxx = hxx1*dx+hxy1*dy
         dyy = hyx1*dx+hyy1*dy
         nx = - dyy  ! normale sortante
         ny = dxx 
         Ux = (u(1,i3)+u(1,i2))*0.5
         Uy = (u(2,i3)+u(2,i2))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy
c
         dx = s1(1,i1)-s1(1,i3)
	 dx = dx - int(dx+dx)
         dy = s1(2,i1)-s1(2,i3)
	 dy = dy - int(dy+dy)
         dxx = hxx1*dx+hxy1*dy
         dyy = hyx1*dx+hyy1*dy
         nx = - dyy  ! normale sortante
         ny = dxx 
         Ux = (u(1,i1)+u(1,i3))*0.5
         Uy = (u(2,i1)+u(2,i3))*0.5
         axx = axx + nx*Ux
         ayy = ayy + ny*Uy
         ayx = ayx + ny*Ux
         axy = axy + nx*Uy
c
         SS = Surf(r1,i1,i2,i3)
         Exx = axx / SS
         Eyy = ayy / SS
         Exy = 0.5d0*(axy+ayx) / SS
         Ew = 0.5d0*(ayx-axy) / SS
         Delta = ((Exx-Eyy)*0.5)**2.d0 + (Exy)**2.d0
         E1 = (Exx+Eyy)*0.5d0 + sqrt(Delta)
         E2 = (Exx+Eyy)*0.5d0 - sqrt(Delta)
         Dis = E1 - E2
         Vol = E1 + E2 
         strain(nt) = Dis
c         Write(67,*) i1,i2,i3,Dis,strain(nt)
      EndDo

      End
c=======================================================================

      Real*8 Function Surf(x,i1,i2,i3)
  
      Implicit None

      Integer i1,i2,i3
      Real*8 x(3,*)
      Real*8 a1,b1,a2,b2

      a1 = x(1,i2)-x(1,i1)
      b1 = x(2,i2)-x(2,i1)
      a2 = x(1,i3)-x(1,i1)
      b2 = x(2,i3)-x(2,i1)
      Surf = abs(a1*b2 - b1*a2)*0.5
   
      Return
      End
c=======================================================================
 
      Subroutine Renum(i01,i02,i03,s,i1,i2,i3)

c        --numerotation des sommets dans le sens inverse trigo
c        --i1 est le sommet completement a gauche
c        --i2 est le sommet le plus haut entre i2 et i3
c        --i3 est le sommet le plus bas entre i2 et i3

      Implicit None

      Integer npam
      Parameter (npam=100000)
c     -- entree
      Integer i01,i02,i03
      Real*8 s(2,npam)
c     -- sortie
      Integer i1,i2,i3
c     -- interne
      Real*8 xi1,xi2,xi3
      Real*8 xmin

      xi1 = s(1,i01)
      xi2 = s(1,i02)
      xi3 = s(1,i03)
      xmin = min(xi1,xi2,xi3)
      If (xmin.EQ.xi1) i1 = i01
      If (xmin.EQ.xi2) i1 = i02
      If (xmin.EQ.xi3) i1 = i03
      If (i1.EQ.i01) Then
         xi2 = s(2,i02)
         xi3 = s(2,i03)
         If (xi2.GT.xi3) Then
            i2 = i02
            i3 = i03
         Else
            i2 = i03
            i3 = i02
         EndIf
      EndIf
      If (i1.EQ.i02) Then
         xi2 = s(2,i01)
         xi3 = s(2,i03)
         If (xi2.GT.xi3) Then
            i2 = i01
            i3 = i03
         Else
            i2 = i01
            i3 = i03
         EndIf
      EndIf
      If (i1.EQ.i03) Then
         xi2 = s(2,i01)
         xi3 = s(2,i02)
         If (xi2.GT.xi3) Then
            i2 = i01
            i3 = i02
         Else
            i2 = i02
            i3 = i01
         EndIf
      EndIf

      End
c=======================================================================
      
      Subroutine Triangule(npa,s,ivnxt,lt)

      Implicit Real*8 (a-h,o-z)
      Integer npam
      Parameter (npam=100000)
      Integer lt(3,npam*5)
      Real*8 s(2,npam)
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
         xcor = s(1,i)
         ycor = s(2,i)
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

c      Write(6,*) 'xmin = ',xmin,' xmax = ',xmax
c      Write(6,*) 'ymin = ',ymin,' ymax = ',ymax

c-----Triangulation
c      Write(6,*) 'Calcul Triangulation...'
      Do i = 1,nv
         is(i) = 0
      EndDo
c     --Triangle initial
      xlen = (xmax-xmin)/4.d0
      ylen = (ymax-ymin)/4.d0
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
      If ((xtemp.GT.-eps).AND.(xtemp.LT.eps)) xtemp = 0.d0
      If ((ytemp.GT.-eps).AND.(ytemp.LT.eps)) ytemp = 0.d0
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
          If ((delx.GT.-eps).AND.(delx.LT.eps)) delx = 0.d0
          If ((dely.GT.-eps).AND.(dely.LT.eps)) dely = 0.d0
          daf = sqrt(delx**2 + dely**2)
          delx = y(i1) - y(i3)
          dely = x(i3) - x(i1)
          If ((delx.GT.-epv).AND.(delx.LT.epv).AND.
     &       (dely.GT.-epv).AND.(dely.LT.epv)) go to 350
          If ((delx.GT.-eps).AND.(delx.LT.eps)) delx = 0.d0
          If ((dely.GT.-eps).AND.(dely.LT.eps)) dely = 0.d0
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

  450 continue
      If (idupl.NE.0) Then
      Write(6,*)'number of points read that were duplications=',idupl
      EndIf
c
  480 continue

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
      Real*8 x(*), y(*)
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
      Real*8 x(*), y(*)
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
      Real*8 x(*), y(*)
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
      Real*8 x(*), y(*)
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
      Real*8 x(*), y(*), norm
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
      Subroutine RetireTrianglesPlats(npa,x,nt,lt)

      Implicit Real*4 (a-h,o-z)
      Integer npam
      Parameter (npam=100000)
      Integer i,npa,nt,ntt
      Integer i0,i1,i2
      Integer lt(3,npam*5),lt0(3,npam*5)
      Real*8 x(2,npam)
      Real*4 Surf1,u01x,u01y,u02x,u02y,u10x,u10y
      Real*4 u21x,u21y,u20x,u20y
      Real*4 s0112,s1220,s2001
      Real*4 d01,d12,d20 
      Real*4 cos0112,cos1220,cos2001
      Real*4 ac0112,ac1220,ac2001
      Real*4 tol,angletol,Srec,pi2deg,S

      Do i = 1,nt
         lt0(1,i) = lt(1,i)
         lt0(2,i) = lt(2,i)
         lt0(3,i) = lt(3,i)
      EndDo

      pi2deg = 180./(4.*atan(1.))

      tol = 0.1
      angletol = 10. !degres

      ntt = 0
      Do i = 1,nt
         i0 = lt0(1,i)
         i1 = lt0(2,i)
         i2 = lt0(3,i)

         u01x = x(1,i1)-x(1,i0)
         u01y = x(2,i1)-x(2,i0)
         u02x = x(1,i2)-x(1,i0)
         u02y = x(2,i2)-x(2,i0)
         u21x = x(1,i1)-x(1,i2)
         u21y = x(2,i1)-x(2,i2)

         d01 = sqrt(u01x*u01x + u01y*u01y)
         d12 = sqrt(u21x*u21x + u21y*u21y)
         d20 = sqrt(u02x*u02x + u02y*u02y)

         s0102 = (u01x*u02x + u01y*u02y)
         cos0102 = (s0102/d20/d01)
         ac0102 = acos(cos0102)*pi2deg 

         u10x = x(1,i0)-x(1,i1)
         u10y = x(2,i0)-x(2,i1)
         u21x = x(1,i2)-x(1,i1)
         u21y = x(2,i2)-x(2,i1)
         s1021= (u10x*u21x + u10y*u21y)
         cos1021 = (s1021/d01/d12)
         ac1021 = acos(cos1021)*pi2deg

         u20x = x(1,i0)-x(1,i2)
         u20y = x(2,i0)-x(2,i2)
         s2021= (u20x*u21x + u20y*u21y)
         cos2021 = (s2021/d20/d12)
         ac2021 = acos(cos2021)*pi2deg

c         Write(67,*) ac0102,ac1021,ac2021,ac0102+ac1021+ac2021
 
         If ((max(d01,d12,d20).GT.tol).OR.
     &       ((ac0102).LT.angletol).OR.
     &        ((ac1021).LT.angletol).OR.
     &        ((ac2021).LT.angletol) 
     &      ) Then
            GoTo 1
         Else
            ntt = ntt + 1
            lt(1,ntt) = i0
            lt(2,ntt) = i1
            lt(3,ntt) = i2
         EndIf
1        Continue
      EndDo

      nt = ntt
   
      Write(6,*) 'Retrait des triangles plats :'
      Write(6,*) ' Nbre total de triangles apres elimination: ',nt

      Return
      End
c=======================================================================
      Subroutine RetireTrianglesPlats2(npa,r1,nt,lt)

      Implicit None
      Integer npam
      Parameter (npam=100000)
      Integer i,npa,nt,ntt
      Integer i0,i1,i2
      Integer lt(3,npam*5),lt0(3,npam*5)
      Real*8 r1(3,npam)
      Real*8 x(2,npam)
      Real*4 Surf1,u01x,u01y,u12x,u12y,u20x,u20y
      Real*4 s0112,s1220,s2001
      Real*4 d01,d12,d20 
      Real*4 cos0112,cos1220,cos2001
      Real*4 ac0112,ac1220,ac2001
      Real*4 tol,angletol,Srec,pi2deg,S

      Do i = 1,npa
         x(1,i) = (r1(1,i))
         x(2,i) = (r1(2,i))
      EndDo
      Do i = 1,nt
         lt0(1,i) = lt(1,i)
         lt0(2,i) = lt(2,i)
         lt0(3,i) = lt(3,i)
      EndDo

      pi2deg = 180./(4.*atan(1.))

      tol = 2.
      angletol = 0. !degres
      Srec = 0.

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
         s0112 = (u01x*u12x + u01y*u12y)
         s1220 = (u12x*u20x + u12y*u20y)
         s2001 = (u20x*u01x + u20y*u01y)
         d01 = sqrt((x(1,i0)-x(1,i1))**2. + (x(2,i0)-x(2,i1))**2.)
         d12 = sqrt((x(1,i1)-x(1,i2))**2. + (x(2,i1)-x(2,i2))**2.)
         d20 = sqrt((x(1,i2)-x(1,i0))**2. + (x(2,i2)-x(2,i0))**2.)
         cos0112 = (s0112/d01/d12)
         cos1220 = (s1220/d12/d20)
         cos2001 = (s2001/d20/d01)
         ac0112 = acos(cos0112)*pi2deg
         ac1220 = acos(cos1220)*pi2deg
         ac2001 = acos(cos2001)*pi2deg 

c         S = Surf1(x,i0,i1,i2)
         Write(67,*) ac0112,ac1220,ac2001,d01,d12,d20
 
         If ((max(d01,d12,d20).GT.tol).OR.
     &       ((ac0112).LT.angletol).OR.
     &        ((ac1220).LT.angletol).OR.
     &        ((ac2001).LT.angletol).OR.
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
      Real*4 Function Surf1(x,i1,i2,i3)
  
      Implicit Real*4 (a-h,o-z)

      Integer npam
      Parameter(npam=300000)
      Integer i1,i2,i3
      Dimension x(2,npam)

      a1 = x(1,i2)-x(1,i1)
      b1 = x(2,i2)-x(2,i1)
      a2 = x(1,i3)-x(1,i1)
      b2 = x(2,i3)-x(2,i1)
      Surf1 = abs(a1*b2 - b1*a2)*0.5
   
      Return
      End
c=======================================================================
     
      Subroutine ModifString(cc)

      Character*5 cc
      Character*1 cc1
      Character*2 cc2
      Character*3 cc3
      Character*4 cc4
      Integer iconf

      Read(cc,*) iconf ! convert string to integer
      If (iconf.LT.10) Then
         Write(cc1,'(i1)') iconf ! convert interer to string
         cc = "000"//cc1
      ElseIf (iconf.LT.100) Then
         Write(cc2,'(i2)') iconf ! convert interer to string
         cc = "00"//cc2
      ElseIf (iconf.LT.1000) Then
         Write(cc3,'(i3)') iconf ! convert interer to string
         cc = "0"//cc3
      EndIf
                   
      Return
      End
c=======================================================================

      Subroutine MeanField(PtINI,PtEND,X,UMF)

      Implicit None

      Integer k,i,i0,i1,i2,i3
      Integer npt
      Parameter (npt = 4)
      Real*8 ptINI(2,npt) ! IN -> initial coord of corners
      Real*8 ptEND(2,npt) ! IN -> final coord of corners
      Real*8 X(2)         ! IN ->  initial coord of the point
      Real*8 UMF(2)       ! OUT -> computed meanfield displacement
      Real*8 xij,yij
      Real*8 x0,y0,x1,y1,x2,y2,x3,y3,xi,yi
      Real*8 s,t,s2,t2
      Real*8 xmean,ymean,xx0,xx1,xx2,xx3,yy0,yy1,yy2,yy3,xxi,yyi,z

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

      If (num_valid_s.EQ.0) Then
         Write(6,*) 'num_valid_s = ',num_valid_s
         Call Exit(2)
      EndIf

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

