c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c
c  Post process conf files files comming from PBC2D program
c
c  generate CONF%d files with:
c
c  number of grains
c  xi yi radii rotation
c  corners positions
c  corresponding deformations
c  scale factor
cs
c
c
c=======================================================================
      Include 'string.inc'

      Program PostProcessDIC

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Character*50 FichDIC
      Character*50 FichConf1,FichConf2
      Dimension FichDIC(NbrePhotosMax)
      Dimension r_grains(6,50000)
      Dimension r_coins(2,ncoins,NbrePhotosMax) 
      Dimension Fxx(NbrePhotosMax),Fyy(NbrePhotosMax)
      Dimension Fxy(NbrePhotosMax),Fyx(NbrePhotosMax)
      Dimension Exx(NbrePhotosMax),Eyy(NbrePhotosMax)
      Dimension Exy(NbrePhotosMax),Evol(NbrePhotosMax)
      Dimension Evide(NbrePhotosMax)
      Dimension Strains(8)
      Dimension pt_fix(3,200,5000) ! 2 coord + 1 rot, 100 ptfix, 5000 photos
      Dimension transf(3,NbrePhotosMax)
      Dimension DistoParam(8)
      Real*8 fact_echelle(NbrePhotosMax)
      Logical lgDICfile(NbrePhotosMax)
      Logical Equal_real8    

      Pi = 4.d0 * atan(1.d0)

c-----Show a title on the current output
      Call Title
      Call Init(NbrePhotos,FichDIC,idisto,DistoParam)

c==============================================================
c     ==Changement de repere et correction translation-rotation
c==============================================================
c
c     -- les points fixes sont les points avec un rayon de 2
c     -- le point 1 est le point de reference
c     -- les deux premiers pt fixes sont les 2 points du bas.
c     -- les autres points fixes sont éventuellement la pour le plaisir !

c-----Nbre de points fixes ?
c      num_photo = 1 ! First DIC file
c      Open(1,file=FichDIC(num_photo),status='old',err=999)
c      Read(1,*) nligne

      npt_fix1 = 0
      npt_fix2 = 0

c      Do i = 1,nligne
c         Read(1,*) z,z,z,rad
c         If (Equal_real8(rad,2.1d0)) Then
c            npt_fix1 = npt_fix1 + 1
c         EndIf 
c         If (Equal_real8(rad,2.2d0)) Then
c            npt_fix2 = npt_fix2 + 1
c         EndIf 
c      EndDo
c      Close(1)
c      Write(6,*) 'Number of points N1 for re-scaling = ', npt_fix1
c      Write(6,*) 'Number of points N2 for undesired movement correction = ', npt_fix2
c      If (npt_fix1.NE.2) Then
c         Write(6,*) 'N1 is not 2! The points for computation of pixel size are missing'
c         Write(6,*) 'Scaling will NOT be done!'
cc         Stop
c      EndIf
c      If (npt_fix2.LT.1) Then
c         Write(6,*) 'N2 < 1! Parasite motions cannot be corrected!!!!'
cc         Stop
c      EndIf 
c      npt_fix = npt_fix2

c-----Calcul du facteur d'echelle
c      If (npt_fix1.EQ.2) Then
c         Call ScaleFactor(NbrePhotos,FichDIC,idisto,DistoParam,
c     &                    fact_echelle)
c      Else
c         Do num_photo = 1,NbrePhotos
c         fact_echelle(num_photo) = 1.d0
c      EndDo
c      EndIf
	
      Do num_photo = 1,NbrePhotos
         fact_echelle(num_photo) = 1.d0
      EndDo

       

c-----Lecture des coord de tous les pt_fix de correction de translation
c      Do num_photo = 1,NbrePhotos
c         Open(1,file=FichDIC(num_photo),status='old')
c         Read(1,*) nligne
c         npt = 0
c         Do i = 1,nligne
c            Read(1,*) xx,yy,rot,rad,dx,dy,drot
c            If (Equal_real8(rad,2.2d0)) Then
c               npt = npt + 1
c               xd = xx + dx  ! distorded X new pos
c               yd = yy + dy  ! distorded Y new pos
c               xu = xd  ! to be undistorded X new pos
c               yu = yd  ! to be undistorded Y new pos
c               If (idisto.GT.0) Call Undisto(DistoParam,xd,yd,xu,yu)
c               pt_fix(1,npt,num_photo) = xu
c               pt_fix(2,npt,num_photo) = yu
c               pt_fix(3,npt,num_photo) = drot
c               If (npt.EQ.npt_fix) GoTo 1
c            EndIf
c         EndDo
c1        Close(1)
c      EndDo

c-----calcul de la moyenne des pos et rot de tous les pts fixes 
c      Open(30,file='checktranslationPTFIX.txt',status='unknown')
c      Do num_photo = 1,NbrePhotos
c         depx = 0.d0
c         depy = 0.d0
c         rot = 0.d0
c         Do iptfix = 1,npt_fix
c            depx = depx + pt_fix(1,iptfix,num_photo)
c            depy = depy + pt_fix(2,iptfix,num_photo)
c            rot = rot + pt_fix(3,iptfix,num_photo)
c         EndDo
c         transf(1,num_photo) = depx / dfloat(npt_fix)
c         transf(2,num_photo) = depy / dfloat(npt_fix)
c         transf(3,num_photo) = rot / dfloat(npt_fix) 
c         Write(30,*) transf(1,num_photo),transf(2,num_photo),
c     &               transf(3,num_photo)
c      EndDo
c      Close(30)
               
c-----lecture des coordonnees des 4 coins
      Call ReadCorners(NbrePhotos,FichDIC,r_coins,
     &                 idisto,DistoParam)
c      Call TransfCorners(NbrePhotos,transf,r_coins,fact_echelle)

c-----Calcul des deformations a partir des 4 coins (hyp: grandes def)
      Write(6,*) 'Strains computation from 1g2e corners...'
      Call ComputeDef(NbrePhotos,r_coins,Fxx,Fxy,Fyx,Fyy,Exx,Exy,
     &                Eyy,Evol,Evide)

c-----Creation des fichiers CONF
      FichConf1 = 'CONF%04d'
      Do num_photo = 1,NbrePhotos
c        -- lecture des pos des grains
         Call ReadGrains(FichDIC(num_photo),idisto,DistoParam,
     &        ngrains,r_grains)
c        -- correction mouvement camera et mise a l'echelle
c         Call TransfGrains(num_photo,transf,ngrains,r_grains,
c     &        fact_echelle)

c        -- Indice des vides
         Vsolid = 0.d0
         Do i = 1,ngrains
            Vsolid = Vsolid + Pi * r_grains(3,i)**2.d0
         EndDo
	 Write(6,*) Vsolid
         Evide(num_photo) = Evide(num_photo) / Vsolid - 1.d0

c        -- Ecriture du fichier CONF correspondant         
         FichConf2 = FichConf1
         Call NumberingString(FichConf2,num_photo) 
         Strains(1) = Fxx(num_photo)
         Strains(2) = Fyy(num_photo)
         Strains(3) = Fxy(num_photo)
         Strains(4) = Fyx(num_photo)
         Strains(5) = Exx(num_photo)
         Strains(6) = Eyy(num_photo)
         Strains(7) = Exy(num_photo)
         Strains(8) = Evol(num_photo)
         Call SaveConfFile(num_photo,FichConf2,ngrains,r_grains,
     &                     r_coins,Strains,fact_echelle)
      EndDo

      Open(15,file='Deformation1G2E.txt')
      Write(15,*) 'Fxx          Fyy          Fxy          '//
     &            'Fyx          Exx          Eyy          '//
     &            'Exy          Evol         Evide'
      Do num_photo = 1,NbrePhotos
         Write(15,'(9(e12.5,1x))') Fxx(num_photo),Fyy(num_photo),
     &                             Fxy(num_photo),Fyx(num_photo),
     &                             Exx(num_photo),Eyy(num_photo),
     &                             Exy(num_photo),Evol(num_photo),
     &                             Evide(num_photo)
      EndDo

      Call Exit(0)

999   Write(6,*) 'File not found'
      Write(6,*) FichDIC(num_photo)
      Write(6,*) 'Application STOPPED!'
      Call Exit(1)      

      End
c=======================================================================
      Subroutine SaveConfFile(num_photo,FileName,ngrains,r_grains,
     &                        r_coins,Strains,fact_echelle)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Character*50 FileName
      Dimension r_grains(6,50000)
      Dimension r_coins(2,ncoins,NbrePhotosMax) 
      Dimension Strains(8)
      Real*8 fact_echelle(NbrePhotosMax)

      Write(6,*) 'Writing file ',FileName
      Open(1,file=FileName,status='unknown')
c
      Write(1,11) ngrains
      Do j = 1,ngrains
         Write(1,12) r_grains(1,j),r_grains(2,j),
     &               r_grains(3,j),r_grains(4,j),
     &               r_grains(5,j),r_grains(6,j)
      EndDo
      Write(1,'(a37)') '0 - corners coordinates - counterclockwise'
      Do k = 1,ncoins
         Write(1,13) r_coins(1,k,num_photo),r_coins(2,k,num_photo)
      EndDo
      Write(1,'(a68)') '1 - Associated macroscopic deformation Fxx,'//
     &                 'Fyy,Fxy,Fyx,Exx,Eyy,Exy,Evol :'
      Write(1,14) (Strains(j),j=1,8)
      Write(1,'(a32)') '2 - scale factor (meters/pixels)'
      Write(1,'(e20.13)') fact_echelle(num_photo)
c
      Close(1)     

11    Format(5(i5,1x))
12    Format(4(e20.13,1x),e10.3,1x,e10.3)
13    Format(2(e15.8,1x))
14    Format(8(e15.8,1x))
      Return
      End
c=======================================================================

      Subroutine Undisto(DistoParam,xd,yd,xu,yu)

      Implicit Real*8 (a-h,o-z)
      Real*8 DistoParam(*)
      Real*8 xd,yd,xu,yu
      Real*8 Xc,Yc ! centre of the distorsion
      Real*8 K1,K2,K3 ! coeff for radial disto
      Real*8 P1,P2,P3 ! coeef for ortho-radial disto
      
      Xc = DistoParam(1)
      Yc = DistoParam(2)
      K1 = DistoParam(3)
      K2 = DistoParam(4)
      K3 = DistoParam(5)
      P1 = DistoParam(6)
      P2 = DistoParam(7)
      P3 = DistoParam(8)

      dx = xd - DistoParam(1)
      dy = yd - DistoParam(2)
      r = sqrt(dx*dx+dy*dy)
      r2 = r*r
      r4 = r2 * r2
      r6 = r2 * r2 * r2
      xu = Xc + dx*(1.d0 + K1*r2 + K2*r4 + K3*r6) +
     &    (P1*(r2+2.0*dx*dx) + 2.0*P2*dx*dy)*(1.0+P3*r2)
      yu = Yc + dy*(1.d0 + K1*r2 + K2*r4 + K3*r6) + 
     &    (P1*(r2+2.0*dy*dy) + 2.0*P2*dx*dy)*(1.0+P3*r2)

      Return
      End
c=======================================================================

      Subroutine TransfCorners(NbrePhotos,transf,r_coins,fact_echelle)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Real*8 r_coins(2,ncoins,NbrePhotosMax) 
      Dimension transf(3,NbrePhotosMax)
      Real*8 fact_echelle(NbrePhotosMax)

      Do num_photo = 1,NbrePhotos
         Do i = 1,ncoins
            r_coins(1,i,num_photo) = (r_coins(1,i,num_photo) - 
     &                               transf(1,num_photo))*
     &                               fact_echelle(num_photo)
            r_coins(2,i,num_photo) = (-(r_coins(2,i,num_photo) -
     &                               transf(2,num_photo)))*
     &                               fact_echelle(num_photo)
         EndDo
      EndDo

      Return
      End
c=======================================================================

      Subroutine TransfGrains(num_photo,transf,ngrains,r_grains
     &                       ,fact_echelle)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Real*8 r_grains(6,50000) 
      Dimension transf(3,NbrePhotosMax)
      Real*8 fact_echelle(NbrePhotosMax)

      Do i = 1,ngrains
         r_grains(1,i) =   r_grains(1,i) - transf(1,num_photo)
         r_grains(2,i) = -(r_grains(2,i) - transf(2,num_photo)) ! on inverse l'axe des y

      r_grains(1,i) = r_grains(1,i)*fact_echelle(num_photo)
      r_grains(2,i) = r_grains(2,i)*fact_echelle(num_photo)
      r_grains(3,i) = r_grains(3,i)*fact_echelle(num_photo)
      EndDo

      Return
      End
c=======================================================================

      Subroutine ReadCorners(NbrePhotos,FichDIC,r_coins,
     &                       idisto,DistoParam)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Character*50 FichDIC(NbrePhotosMax)
      Real*8 r_coins(2,ncoins,NbrePhotosMax) 
      Dimension DistoParam(*)
      Character*2 KeyWord_h

      Write(6,*) 'Reading h components...'
      Do i = 1,NbrePhotos
         Open(1,file=FichDIC(i),status='old')
         KeyWord_h = ''
         Do While (KeyWord_h.NE.'h ')
            Read(1,'(a2)') KeyWord_h
         EndDo
         Backspace(1)
         Read(1,*) KeyWord_h,hxx,hxy,hyx,hyy
         Close(1)
         r_coins(1,1,i) = hxx
         r_coins(1,2,i) = hxy
         r_coins(1,3,i) = hyx
         r_coins(1,4,i) = hyy
      EndDo

      Return
      End
c=======================================================================

      Subroutine ComputeDef(NbrePhotos,r_coins,Fxx,Fxy,Fyx,Fyy,
     &                Exx,Exy,Eyy,Evol,Evide)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Real*8 r_coins(2,ncoins,NbrePhotosMax) 
      Real*8 Fxx,Fyy,Fxy,Fyx
      Real*8 Exx,Eyy,Exy,Evol,Evide
      Real*8 hxx0,hyy0,hxy0,hyx0
      Real*8 hxx,hyy,hxy,hyx
      Dimension Fxx(NbrePhotosMax),Fyy(NbrePhotosMax)
      Dimension Fxy(NbrePhotosMax),Fyx(NbrePhotosMax)
      Dimension Exx(NbrePhotosMax),Eyy(NbrePhotosMax)
      Dimension Exy(NbrePhotosMax),Evol(NbrePhotosMax)
      Dimension Evide(NbrePhotosMax)
      Integer i1,i2,i3,i4

      num_photo = 1

c-----numerotation des sommets dans le sens trigo
c     --i1 est le coins an bas a gauche 
c     -- mean value for X and Y coordinates
c      xmean = 0.d0
c      ymean = 0.d0
c      Do i = 1,ncoins
c         xmean = xmean + r_coins(1,i,num_photo)
c         ymean = ymean + r_coins(2,i,num_photo)
c      EndDo
c      xmean = xmean / dfloat(ncoins)
c      ymean = ymean / dfloat(ncoins)
c     -- the vector from the center to i1 gives (cos,sin) < 0
c      Do i = 1,ncoins
c         xij = r_coins(1,i,num_photo) - xmean
c         yij = r_coins(2,i,num_photo) - ymean
c         If ( (xij.LT.0.d0).AND.(yij.LT.0.d0) ) i1 = i
c         If ( (xij.GT.0.d0).AND.(yij.LT.0.d0) ) i2 = i
c         If ( (xij.GT.0.d0).AND.(yij.GT.0.d0) ) i3 = i
c         If ( (xij.LT.0.d0).AND.(yij.GT.0.d0) ) i4 = i
c      EndDo

c-----Forme intiale du trapeze et determinant, num_photo = 1
      hxx0 = r_coins(1,1,num_photo) 
      hxy0 = r_coins(1,2,num_photo) 
      hyx0 = r_coins(1,3,num_photo) 
      hyy0 = r_coins(1,4,num_photo) 
      deth = hxx0 * hyy0 - hyx0 * hxy0

c-----Boucle sur toutes les photos
      Do num_photo = 1,NbrePhotos
         hxx = r_coins(1,1,num_photo) 
         hxy = r_coins(1,2,num_photo) 
         hyx = r_coins(1,3,num_photo) 
         hyy = r_coins(1,4,num_photo) 
         Evide(num_photo) = hxx*hyy-hyx*hxy
c        -- gradient de transformation
         Fxx(num_photo) = (hxx * hyy0 - hxy * hyx0) / deth
         Fxy(num_photo) = (hxy * hxx0 - hxx * hxy0) / deth
         Fyx(num_photo) = (hyx * hyy0 - hyy * hyx0) / deth
         Fyy(num_photo) = (hyy * hxx0 - hyx * hxy0) / deth
c           -- tenseur de def (Green-Lagrange)
         Exx(num_photo) = 0.5d0 * (Fxx(num_photo) * Fxx(num_photo) +
     &                           Fyx(num_photo) * Fyx(num_photo) - 1.d0)
         Exy(num_photo) = 0.5d0 * (Fxx(num_photo) * Fxy(num_photo) +
     &                           Fyy(num_photo) * Fyx(num_photo))
         Eyy(num_photo) = 0.5d0 * (Fyy(num_photo) * Fyy(num_photo) +
     &                           Fxy(num_photo) * Fxy(num_photo) - 1.d0)
c        -- deformation volumique
         Evol(num_photo) = Fxx(num_photo) * Fyy(num_photo) -
     &                     Fxy(num_photo) * Fyx(num_photo) - 1.d0
      EndDo

      Return
      End
c======================================================================
      Subroutine ScaleFactor(NbrePhotos,FichDIC,idisto,DistoParam,
     &                       fact_echelle)

      Implicit Real*8 (a-h,o-z)

      Integer NbrePhotosMax
      Parameter (NbrePhotosMax = 5000) 
      Character*50 FichDIC(NbrePhotosMax)
      Real*8 ptfix0(2,2)
      Real*8 DistoParam(8)
      Integer num_photo,npt_fix1,idisto,NbrePhotos
      Integer i,nligne
      Real*8 vrai_dist
      Real*8 xx,yy,rot,rad,dx,dy
      Real*8 xd,yd,xu,yu
      Real*8 vec0_x,vec0_y,vec0
      Real*8 fact_echelle(NbrePhotosMax)
      Logical Equal_Real8

      vrai_dist = 0.3355d0 ! metres
      
      Do num_photo = 1,NbrePhotos
         Open(1,file=FichDIC(num_photo),status='old')
         Read(1,*) nligne
         npt_fix1 = 0
         Do i = 1,nligne
            Read(1,*) xx,yy,rot,rad,dx,dy
            If (Equal_Real8(rad,2.1d0)) Then 
               xd = xx + dx
               yd = yy + dy
               xu = xd
               yu = yd
               If (idisto.GT.0) Call Undisto(DistoParam,xd,yd,xu,yu)
               npt_fix1 = npt_fix1 + 1
               ptfix0(1,npt_fix1) = xu
               ptfix0(2,npt_fix1) = yu
            EndIf 
            If (npt_fix1.EQ.2) GoTo 1
         EndDo
1        Close(1)
c        -- vecteur joignant les deux pt fixes sur la photo initiale
         vec0_x = ptfix0(1,2) - ptfix0(1,1)
         vec0_y = ptfix0(2,2) - ptfix0(2,1)
         vec0 = sqrt(vec0_x * vec0_x + vec0_y * vec0_y)
         fact_echelle(num_photo) = vrai_dist / vec0 ! (m/pixels)
      EndDo

      Return 
      End
c=======================================================================

      Subroutine ReadGrains(FileName,idisto,DistoParam,
     &           ngrains,r_grains)

      Implicit Real*8 (a-h,o-z)

      Character*50 FileName
      Dimension r_grains(6,50000)
      Dimension DistoParam(8)
      Character*9 KeyWord_Particles

      Open(1,file=FileName,status='old')

      Rewind(1)
      KeyWord_Particles = ''
      Do While (KeyWord_Particles.NE.'Particles')
         Read(1,'(a9)') KeyWord_Particles
      EndDo
      Backspace(1)
      Read(1,*) KeyWord_Particles,nligne

      ngrains = 0
      Do i = 1,nligne
         Read(1,*) xu,yu,z,z,z,z,rot,z,z,radii
         ngrains = ngrains + 1
         r_grains(1,ngrains) = xu
         r_grains(2,ngrains) = yu
         r_grains(3,ngrains) = radii
         r_grains(4,ngrains) = rot
         r_grains(5,ngrains) = 1.
         r_grains(6,ngrains) = 1.
      EndDo

      Close(1)

      Return
      End
c======================================================================

      Subroutine Init(NbrePhotos,FichDIC,idisto,DistoParam)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Character*50 FichDIC,fichTMP
      Dimension FichDIC(NbrePhotosMax)
      Character*50 NomFichDIC1
      Dimension DistoParam(8)
      Character*50 cdisto
      Logical lgYES

c-----Input DIC files
      Write(6,*) '-------------------------------'
      Write(6,*) ' conf to be readed '
      Write(6,*) ' - - - - - - - - - - - - - - -'
      Write(6,*) ' Ex: ../conf%4d.txt' 
      Read(5,'(a50)') NomFichDIC1 
      Write(6,*) 'NUMBERS: First, Last dic file number'
      Write(6,*) ' --> they must exist!' 
      Read(5,*) iDICdeb,iDICfin
      
c     -- we check that the first and the last DIC files exists.
c     ------ first DIC
      fichTMP = NomFichDIC1
      Call NumberingString(fichTMP,iDICdeb) 
      lgYES = .FALSE.
      Inquire(file=fichTMP,exist=lgYES)
      If (.NOT.lgYES) Then
         Write(6,*) 'The file ',fichTMP,' does not exist'
         Write(6,*) 'Program STOPPED'
         Call Exit(1)
      EndIf
c     ------ last DIC
      fichTMP = NomFichDIC1
      Call NumberingString(fichTMP,iDICfin) 
      lgYES = .FALSE.
      Inquire(file=fichTMP,exist=lgYES)
      If (.NOT.lgYES) Then
         Write(6,*) 'The file ',fichTMP,' does not exist'
         Write(6,*) 'Program STOPPED'
         Call Exit(1)
      EndIf
      
c     -- we make the list of exsiting DIC files
      k = 0
      Do i = iDICdeb,iDICfin
      fichTMP = NomFichDIC1
         Call NumberingString(fichTMP, i)
         Inquire(file=fichTMP,exist=lgYES)
      If (lgYES) Then
        k = k + 1
        FichDIC(k) = fichTMP
      EndIf
      EndDo
      NbrePhotos = k
      Write(6,*) 'Total number of DICfiles: ', NbrePhotos 

c-----Parameters needed to undistor the coordinates
c      Write(6,*) 'Correction of distorsion? (n,N,no,No, y,Y,yes,Yes)'
c      Write(6,*) ' if yes, fake_undistor only!!!'
c      Read(5,*) cdisto
      cdisto = 'n'
      If ( (cdisto.EQ.'y').OR.(cdisto.EQ.'Y')
     &     .OR.
     &     (cdisto.EQ.'yes').OR.(cdisto.EQ.'Yes')
     &   ) Then
         Write(6,*) 'Lens distortion WILL be corrected'
         idisto = 1
      Else
         Write(6,*) 'Lens distortion WILL NOT be corrected'
         idisto = 0
      EndIf
      If (idisto.GT.0) Then
         Write(6,*) 'Parameters for distortion correction'
         Write(6,*) 'Centre of distortion: xc,yc >'
         Read(5,*) DistoParam(1),DistoParam(2)
         Write(6,*) 'For radial correction: K1,K2,K3 >'
         Read(5,*) DistoParam(3),DistoParam(4),DistoParam(5)
         Write(6,*) 'For tangential correction: P1,P2,P3 >'
         Read(5,*) DistoParam(6),DistoParam(7),DistoParam(8)
      EndIf

      Return
      End
c=======================================================================

      Subroutine Title

      Write(6,*) '               __    _______  _______  _______ '
      Write(6,*) '              /  \  (  ____ \/ ___   )(  ____ \'
      Write(6,*) '              \/) ) | (    \/\/   )  || (    \/'
      Write(6,*) '                | | | |          /   )| (__    '
      Write(6,*) '                | | | | ____   _/   / |  __)   '
      Write(6,*) '                | | | | \_  ) /   _/  | (      '
      Write(6,*) '              __) (_| (___) |(   (__/\| (____/\'
      Write(6,*) '              \____/(_______)\_______/(_______/'
      Write(6,*) ' '
      Write(6,*) '  __                     '
      Write(6,*) ' /  |           /        '
      Write(6,*) '(___| ___  ___ (___  ___ '
      Write(6,*) '|    |   )|___ |         '
      Write(6,*) '|    |__/  __/ |__       '
      Write(6,*) '  '
      Write(6,*) '                  __	                             '
      Write(6,*) '                 /  |				     '
      Write(6,*) '          ___   (___| ___  ___  ___  ___  ___  ___ '
      Write(6,*) '                |    |   )|   )|    |___)|___ |___ '
      Write(6,*) '                |    |    |__/ |__  |__   __/  __/ ' 
      Write(6,*) ' '
      Write(6,*) 'ver du 13 Mars 2015, gael'
      Write(6,*) '    '

      Return
      End
c======================================================================

      Logical Function Equal_Real8(x, y)

      Real*8 x
      Real*8 y
      Real*8 diff
      real*8 epsilon
      
      epsilon = 1.e-14
      diff = abs(x - y)
      If (diff.LT.epsilon) Then
         Equal_Real8 = .TRUE. 
      Else
         Equal_Real8 = .FALSE. 
      EndIf
            
      End
c======================================================================
     


