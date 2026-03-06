c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c
c  Post process dic_out_%d.txt files comming from Tracker program
c
c  generate CONF%d files with:
c
c  number of grains
c  xi yi radii rotation
c  corners positions
c  corresponding deformations
c  scale factor
c
c
c  vincent.richefeu@hmg.inpg.fr
c  gael.combe@hmg.inpg.fr (*)
c
c January 2012
c
c
cmodifications in PostProcessDIC prog: 2012/01/10
c  - in dic file, corners must have a radii equal to 1
c    4 corners are needed.
c
c  - in dic file, fixed points must have a radii equal to 2.
c  the two first fixed point must be the two "down" points. All
c  the other fixed points can be anywhere.
c
c
c  in dic file, corners and fixed point can be anywhere.
c  there is no order for corners because they are renumbered
c  (clockwise - must be checked) by the program
c
c  April 2012 : use distorsion parameters to correct coordinates
c 
c
c


c=======================================================================
      Program PostProcessDIC

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Character*1 a1num
      Character*2 a2num
      Character*3 a3num
      Character*4 a4num
      Character*50 FichDIC,FichConf
      Character*50 NomFichDIC
      Character*50 PATH
      Integer iDICdeb,iDICfin,npath
      Dimension FichDIC(NbrePhotosMax)
      Dimension FichConf(NbrePhotosMax)
      Dimension r_grain(4,5000)
      Dimension real_coord(2,5000)
      Dimension r(2,ncoins,NbrePhotosMax) 
      Dimension x1(2,ncoins),x2(2,ncoins)
      Dimension u(2,ncoins)
      Dimension Gxx(NbrePhotosMax),Gyy(NbrePhotosMax)
      Dimension Gxy(NbrePhotosMax),Gyx(NbrePhotosMax)
      Dimension Evol(NbrePhotosMax)
      Dimension Exx(NbrePhotosMax),Eyy(NbrePhotosMax)
      Dimension Exy(NbrePhotosMax),Erot(NbrePhotosMax)
      Dimension Orienta(NbrePhotosMax,5000)
      Real*8 nx,ny
      Dimension ptfix0(2,200),ptfix1(2,200)
      Dimension transf(3,NbrePhotosMax)
      Dimension DistoParam(8)

      Pi = 4.d0*atan(1.d0)

      Write(6,*) '--------------------------'
      Write(6,*) 'Names of files to be read :'
      Write(6,*) ' NOM = dic_out_'
      Write(6,*) ' - - - - -'
      NomFichDIC = 'dic_out_'
      n1 = iCharlen(NomFichDIC)
      Write(6,*) 'NUMBERS : First, End'
      Read(5,*) iDICdeb,iDICfin



      PATH = '../'


      npath = iCharlen(PATH)
      Write(6,*) "PATH = ",PATH(1:npath)
      Write(6,*) " " 
      k = 0
      Do i = iDICdeb,iDICfin
         k = k + 1
         If (i.LT.10) Then
            Write(a1num,'(i1)') i
            FichDIC(k) = PATH(1:npath)//NomFichDIC(1:n1)//a1num//'.txt'
         Else
            If (i.LT.100) Then
               Write(a2num,'(i2)') i
               FichDIC(k) = PATH(1:npath)//NomFichDIC(1:n1)//
     &                      a2num//'.txt'
            Else
	       If (i.LT.1000) Then
                  Write(a3num,'(i3)') i
                  FichDIC(k) = PATH(1:npath)//NomFichDIC(1:n1)//
     &                         a3num//'.txt'
	       Else
                  Write(a4num,'(i4)') i
                  FichDIC(k) = PATH(1:npath)//NomFichDIC(1:n1)//
     &            a4num//'.txt'
	       EndIf
            EndIf
         EndIf
      EndDo
 
      NbrePhotos = k
      Write(6,*) 'Numbers of pictures : ',NbrePhotos 

      Write(6,*) 'Correction of distorsion ? (0 = no, other = yes)'
      Read(5,*) idisto
      If (idisto.GT.0) Then
c         Write(6,*) 'Parameters for distorsion correction ?'
c         Write(6,*) 'center of distorsion: xc,yc ?'
c         Read(5,*) DistoParam(1),DistoParam(2)
c         Write(6,*) 'For radial correction : K1,K2,K3 ?'
c         Read(5,*) DistoParam(3),DistoParam(4),DistoParam(5)
c         Write(6,*) 'For tangential correction : P1,P2,P3 ?'
c         Read(5,*) DistoParam(6),DistoParam(7),DistoParam(8)
          DistoParam(1) = 3024  
          DistoParam(2) = 2016
          DistoParam(3) = -5.22618E-10
          DistoParam(4) = 5.65538E-20
          DistoParam(5) = 0.0
          DistoParam(6) = -7.55309E-8
          DistoParam(7) = -2.0303E-7
          DistoParam(8) = 0.0
      EndIf





      fact_echelle = 1.0

c==================================
c     ==Depouillement des coord  ==
c==================================
c
      Do i = 1,9
         Write(a1num,'(i1)') i
         FichConf(i) ='CONF000'//a1num
      EndDo
      Do i = 10,99
         Write(a2num,'(i2)') i
         FichConf(i) ='CONF00'//a2num
      EndDo
      Do i = 100,999
         Write(a3num,'(i3)') i
         FichConf(i) ='CONF0'//a3num
      EndDo
      Do i = 1000,NbrePhotos
         Write(a4num,'(i4)') i
         FichConf(i) ='CONF'//a4num
      EndDo
c
c     --lecture du nombre de points
      i = 1
      Open(1,file=FichDIC(i),status='old')
      Read(1,*) ngrains
      Close(1) 
      ngrains = ngrains
      Write(6,*) 'nbre de grains =',ngrains
c
      Do i = 1,NbrePhotos
c        -- lecture du fichier
         Open(1,file=FichDIC(i),status='old')
	 Read(1,*) ng
         ngrains = 0
         Do j = 1,ng
            Read(1,*) xini,yini,z,radii,dx,dy,drot
            If (radii.NE.2) Then
               ngrains = ngrains + 1
               xd = xini + dx
               xu = xd
               yd = yini + dy
               yu = yd
               If (idisto.GT.0) Call Undisto(DistoParam,xd,yd,xu,yu)
               xj = xu
               yj = yu 
               xi = xj
               yi = yj
               r_grain(1,ngrains) = xi * fact_echelle
               r_grain(2,ngrains) = yi * fact_echelle
               r_grain(3,ngrains) = drot
               r_grain(4,ngrains) = radii * fact_echelle
            EndIf
         EndDo
         Close(1)
c        -- ecriture dans le fichier CONF
         If ((i.EQ.1).OR.(mod(i,100).EQ.0).OR.(i.EQ.NbrePhotos)) Then
            Write(6,*) 'Ecriture du fichier ',FichConf(i)
         EndIf
         Open(1,file=FichConf(i))
c         Write(1,11) ngrains
         Do j = 1,ngrains
            Write(1,12) r_grain(1,j),r_grain(2,j),
     &                  r_grain(3,j),r_grain(4,j)
         EndDo
      EndDo
      Call Exit(2)

11    Format(5(i5,1x))
12    Format(9(e15.8,1x))
13    Format(2(e15.8,1x))

999   Write(6,*) 'Fichier inexistant'
      Write(6,*) FichDIC(num_photo)
      Write(6,*) 'Prog STOPPE !'
      

      End
c=======================================================================

      Real*8 Function Surf1(x,i1,i2,i3)
  
      Implicit Real*8 (a-h,o-z)
      Integer i1,i2,i3
      Dimension x(2,*)

      a1 = x(1,i2)-x(1,i1)
      b1 = x(2,i2)-x(2,i1)
      a2 = x(1,i3)-x(1,i1)
      b2 = x(2,i3)-x(2,i1)
      Surf1 = abs(a1*b2 - b1*a2)*0.5d0
   
      Return
      End
c=======================================================================

      Integer Function iCharLen(text)
      
      Integer i
      Character*50 text
      Character*1 cc

      i = 0
1     i = i + 1
      cc = text(i:i)
      If (cc.NE.'') GoTo 1
      
      iCharLen = i - 1

      Return
      End

c======================================================================

      Subroutine Undisto(DistoParam,xd,yd,xu,yu)

      Implicit Real*8 (a-h,o-z)
      Real*8 DistoParam(*)
      Real*8 xd,yd,xu,yu

      ftheta = 1.d0 !no correction for parallaxe

      dx = xd - DistoParam(1)
      dy = yd - DistoParam(2)
      r = sqrt(dx*dx+dy*dy)
      r2 = r*r
      xu = xd + ftheta*dx*(DistoParam(3)*r2 + DistoParam(4)*r2*r2 + 
     &                     DistoParam(5)*r2*r2*r2) +
     &    (DistoParam(6)*(r2 + 2.0*dx*dx) + 2.0*Distoparam(7)*dx*dy)*
     &    (1.0 + DistoParam(8)*r2)
      yu = yd + ftheta*dy*(DistoParam(3)*r2 + DistoParam(4)*r2*r2 +
     &     DistoParam(5)*r2*r2*r2) + (DistoParam(7)*(r2 + 2.0*dy*dy) + 
     &     2.0*DistoParam(6)*dx*dy)*(1.0 + DistoParam(8)*r2)

      Return
      End

c======================================================================

