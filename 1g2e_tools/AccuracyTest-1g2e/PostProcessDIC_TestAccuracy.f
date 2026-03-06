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
c  gael.combe@grenoble-inp.fr
c
c Fevrier 2015
c Mars 2015
c   - les ZNCC sont maintenant ecrits dans les fichiers CONF
c   - correction bug : les calculs de def ne fonctionnaient pas
c                      en l'absence de correction de disto
c   - in dic file, corners must have a radii equal to 1.
c       4 corners are needed.
c   - in dic file, fixed points must have a radii equal to 
c       2.1 for the two fixed points used to rescale the coord
c       2.2 for all the fixed points used correct the translation+rot
c
c   ATTENTION : la rotation de la camera n'est pas corrigee!
c
c Octobre 2018
c   Quelques modifs mineures: 
c     Les DIC peuvent avoir des nums qui ne se suivent pas
c Novembre 2019
c   Quelques modifications de forme (anglais + ajout d'une mini-doc)
c
c Octobre 2021
c     prise en compte de fake_undisto 0 et 1
c     deformation en HPP -- mais calcul du gradiant de tranformation 
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
      Dimension transf(5,NbrePhotosMax)
      Dimension DistoParam(8)
      Real*8 fact_echelle(NbrePhotosMax)
c      Logical lgDICfile(NbrePhotosMax)
      Logical Equal_real8    
      Dimension seg0(50000),seg1(50000)
      Dimension segx_0(50000),segy_0(50000)
      Dimension segx_1(50000),segy_1(50000)
      Dimension seg2(50000),rlon(50000),CC(50000)
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
      num_photo = 1 ! First DIC file
      Open(1,file=FichDIC(num_photo),status='old',err=999)
      Read(1,*) nligne
      npt_fix1 = 0
      npt_fix2 = 0
      Do i = 1,nligne
         Read(1,*) z,z,z,rad
         If (Equal_real8(rad,2.1d0)) Then
            npt_fix1 = npt_fix1 + 1
         EndIf 
         If (Equal_real8(rad,2.2d0)) Then
            npt_fix2 = npt_fix2 + 1
         EndIf 
      EndDo
      Close(1)
      Write(6,*) 'Number of points N1 for re-scaling = ', npt_fix1
      Write(6,*) 'Number of points N2 for undesired movement '//
     &           'correction = ', npt_fix2
      If (npt_fix1.NE.2) Then
         Write(6,*) 'N1 is not 2! The points for computation '//
     &              'of pixel size are missing'
         Write(6,*) 'Scaling will NOT be done!'
c         Stop
      EndIf
      If (npt_fix2.LT.1) Then
         Write(6,*) 'N2 < 1! Parasite motions cannot be corrected!!!!'
c         Stop
      EndIf 
      npt_fix = npt_fix2

c-----Calcul du facteur d'echelle
      If (npt_fix1.EQ.2) Then
         Call ScaleFactor(NbrePhotos,FichDIC,idisto,DistoParam,
     &                    fact_echelle)
      Else
         Do num_photo = 1,NbrePhotos
         fact_echelle(num_photo) = 1.d0
      EndDo
      EndIf

c-----Lecture des coord de tous les pt_fix de correction de translation
      Do num_photo = 1,NbrePhotos
         Open(1,file=FichDIC(num_photo),status='old')
         Read(1,*) nligne
         npt = 0
         Do i = 1,nligne
            Read(1,*) xx,yy,rot,rad,dx,dy,drot
            If (Equal_real8(rad,2.2d0)) Then
               npt = npt + 1
               xd = xx + dx  ! distorded X new pos
               yd = yy + dy  ! distorded Y new pos
               xu = xd  ! to be undistorded X new pos
               yu = yd  ! to be undistorded Y new pos
               Call UnDisto(DistoParam,xd,yd,xu,yu,idisto)
               pt_fix(1,npt,num_photo) = xu
               pt_fix(2,npt,num_photo) = yu
               pt_fix(3,npt,num_photo) = rot + drot
               If (npt.EQ.npt_fix) GoTo 1
            EndIf
         EndDo
1        Close(1)
      EndDo

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
c         transf(4,num_photo) = 0.d0
c         transf(5,num_photo) = 0.d0
c         Write(30,*) transf(1,num_photo),transf(2,num_photo),
c     &               transf(3,num_photo)
c      EndDo
c      Close(30)

c     --calcul pt de rotation instantannee
c     ----Coord du point A a t0
c      num_photo = 1
c      Do iptfix = 1,npt_fix
c         transf(4,num_photo) = transf(4,num_photo) + 
c     &                         pt_fix(1,iptfix,num_photo)
c         transf(5,num_photo) = transf(5,num_photo) + 
c     &                         pt_fix(2,iptfix,num_photo)
c      EndDo
c      transf(4,num_photo) = transf(4,num_photo) / dfloat(npt_fix)
c      transf(5,num_photo) = transf(5,num_photo) / dfloat(npt_fix)
c      Write(6,*) "Coord du pt A a t0 = ",
c     &            transf(4,num_photo),transf(5,num_photo)
c      Do num_photo = 2,NbrePhotos
cc        -- coord du point A a t
c         Do iptfix = 1,npt_fix
c            transf(4,num_photo) = transf(4,num_photo) + 
c     &                            pt_fix(1,iptfix,num_photo)
c            transf(5,num_photo) = transf(5,num_photo) + 
c     &                            pt_fix(2,iptfix,num_photo)
c         EndDo
c         transf(4,num_photo) = transf(4,num_photo) / dfloat(npt_fix)
c         transf(5,num_photo) = transf(5,num_photo) / dfloat(npt_fix)
c        -- coord du point B a t (CIR)
c         Bx0 = - transf(2,num_photo) / transf(3,num_photo) - 
c     &           transf(4,1)     
c         By0 = - transf(1,num_photo) / transf(3,num_photo) -  
c     &           transf(4,2)  
c         Write(55,*) Bx0,  By0
c      EndDo  

c-----lecture des coordonnees des 4 coins
      Call ReadCorners(NbrePhotos,FichDIC,r_coins,
     &                 idisto,DistoParam)
      Call TransfCorners(NbrePhotos,transf,r_coins,fact_echelle)

c-----Calcul des deformations a partir des 4 coins (hyp: grandes def)
c      Write(6,*) 'Strains computation from 1g2e corners...'
c      Call ComputeDefHPP(NbrePhotos,r_coins,Fxx,Fxy,Fyx,Fyy,Exx,Exy,
c     &                Eyy,Evol,Evide)

c-----Calcul des tailles initialles des segements
      num_photo = 1
c        -- lecture des pos des grains
         Call ReadGrains(FichDIC(num_photo),idisto,DistoParam,
     &        ngrains,r_grains)
c   
      n_seg = ngrains / 2
      Write(6,*) "Number of sgments : ",n_seg
      k1 = 0
      k2 = 0
      Do i = 1,n_seg
         k1 = k2 + 1
         k2 = k1 + 1
         x12 = r_grains(1,k2) - r_grains(1,k1)
         y12 = r_grains(2,k2) - r_grains(2,k1)
         seg0(i) = sqrt(x12*x12+y12*y12)
         segx_0(i) = x12/seg0(i)
         segy_0(i) = y12/seg0(i)
      EndDo

      ncpt = 0
      Do num_photo = 1,NbrePhotos
c        -- lecture des pos des grains
         Call ReadGrains(FichDIC(num_photo),idisto,DistoParam,
     &        ngrains,r_grains)
c        -- calcul des tailles courantes de segments
         k1 = 0
         k2 = 0
         Do i = 1,n_seg
            k1 = k2 + 1
            k2 = k1 + 1
            x12 = r_grains(1,k2) - r_grains(1,k1)
            y12 = r_grains(2,k2) - r_grains(2,k1)
	    seg2(i) = sqrt(x12*x12+y12*y12) ! longueur 
            seg1(i) = (seg2(i) - seg0(i)) ! variation de longueur
c            segx_1(i) = seg1(i)*segx_0(i)
c            segy_1(i) = seg1(i)*segy_0(i)
c            Write(43,*) segx_1(i),segy_1(i),
c     &                  atan(segy_0(i)/segx_0(i))*180/3.14159
            Write(43,*) seg0(i),seg2(i)
         EndDo
         Call MeanStdev(n_seg,seg1,rmean1,stdev1)
         Call MeanStdev(n_seg,seg2,rmean2,stdev2)
         Call MeanStdev(n_seg,segx_1,rmeanx1,stdevx1)
         Call MeanStdev(n_seg,segy_1,rmeany1,stdevy1)
	 Do i = 1,ngrains
	    CC(i) = r_grains(5,i)
	    rlon(i) = r_grains(6,i)
	 EndDo
	 Call MeanStDev(ngrains,CC,rmeanCC,stdevCC)
	 Call MeanStDev(ngrains,rlon,rmeanrl,stdevrl)
	 ncpt = ncpt + 1
         Write(39,'(i4,15(1x,e12.5))') ncpt,Exx(ncpt),
     &                                Eyy(ncpt),Exy(ncpt),
     &                                rmean1,stdev1,
     &                                rmean2,stdev2,	 
     &                                rmeanCC,stdevCC,
     &                                rmeanrl,stdevrl,
     &                                rmeanx1,stdevx1,
     &                                rmeany1,stdevy1
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

c      Write(6,*) 'Writing file ',FileName
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

      Subroutine UnDisto(DistoParam,xd,yd,xu,yu,idisto)
      
      Implicit Real*8 (a-h,o-z)
      Real*8 DistoParam(*)
      Real*8 xd,yd,xu,yu
c      Real*8 Xc,Yc ! centre of the distorsion
c      Real*8 K1,K2,K3 ! coeff for radial disto
c      Real*8 P1,P2,P3 ! coeef for ortho-radial disto
      Parameter (niterMax=100,nconv=3)
      Parameter (tol=1.d-3)
      Real*8 fx,fy
      Real*8 dxd_dxu,dyd_dyu

c     -- no undisto
      If (idisto.GT.1) Then ! on ne corrige rien
         xu = xd
         yu = yd
         Return
      EndIf
c     -- fake_undisto // fake_undisto = 1 dans tracker
      If (idisto.EQ.1) Then
         Call Disto(DistoParam,xu,yu,xd,yd) ! on l'utilise à l'envers
	 Return
      EndIf
c     -- real_undisto // fake_undisto = 0 dans tracker 
      xu = xd  ! initialisation
      yu = yd  ! initialisation
            	 
      kconv = 0
      nbiter = 0
      
      Do i = 1,niterMax
      	 nbiter = nbiter + 1
	 Call diff_Disto(DistoParam,xu,yu,dxd_dxu,dyd_dyu) ! derivee
	 Call Disto(DistoParam,fx,fy,xu,yu) ! nouvelle estimation
	 If (abs(dxd_dxu).GT.1.d-15) Then
	    deltax = (xd-fx) / dxd_dxu
	 Else
	    deltax = 0.d0
	 EndIf
	 If (abs(dyd_dyu).GT.1.d-15) Then
	    deltay = (yd-fy) / dyd_dyu
	 Else
	    deltay = 0.d0
	 EndIf
	 xu = xu + deltax
	 yu = yu + deltay  
	 If ((deltax.LT.tol).AND.(deltay.LT.tol)) Then
	    kconv = kconv + 1
	 Else
	    kconv = 0
	 EndIf
	 If (kconv.GE.nconv) Return
      EndDo
      
      Return
      End

c=======================================================================

      Subroutine Disto(DistoParam,xd,yd,xu,yu)

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

      dx = xu - Xc
      dy = yu - Yc
      r2 = dx*dx+dy*dy
      r4 = r2 * r2
      r6 = r2 * r2 * r2
      xd = Xc + dx*(1.d0 + K1*r2 + K2*r4 + K3*r6) +
     &    (P1*(r2+2.0*dx*dx) + 2.0*P2*dx*dy)*(1.0+P3*r2)
      yd = Yc + dy*(1.d0 + K1*r2 + K2*r4 + K3*r6) + 
     &    (P1*(r2+2.0*dy*dy) + 2.0*P2*dx*dy)*(1.0+P3*r2)

      Return
      End
c=======================================================================

      Subroutine diff_Disto(DistoParam,xu,yu,dxd_dxu,dyd_dyu)

      Implicit Real*8 (a-h,o-z)
      Real*8 DistoParam(*)
      Real*8 xu,yu,dxd_dxu,dyd_dyu
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

      dx = xu - Xc
      dy = yu - Yc
      r2 = dx*dx+dy*dy
      r4 = r2 * r2
      r6 = r2 * r2 * r2
      dxd_dxu = 1.d0 + K1*r2 + K2*r4 + K3*r6 + 
     &          2.d0*dx*dx*(2.d0*K2*r2 + K1 + 3.d0*K3*r4) +
     &          (6.d0*P1*dx + 2.d0*P2*dy)*(1.d0 + P3*r2)  +
     &          (P1*(r2 + 2.d0*dx*dx) + 2.d0*P2*dx*dy)*P3*2.d0*dx

      dyd_dyu = 1.d0 + K1*r2 + K2*r4 + K3*r6 + 
     &          2.d0*dy*dy*(2.d0*K2*r2 + K1 + 3.d0*K3*r4) +
     &          (6.d0*P1*dy + 2.d0*P2*dx)*(1.d0 + P3*r2)  +
     &          (P1 * (r2 + 2.d0*dy*dy) + 2.d0*P2*dy*dx)*P3*2.d0*dy

      Return
      End

c=======================================================================

      Subroutine TransfCorners(NbrePhotos,transf,r_coins,fact_echelle)

      Implicit Real*8 (a-h,o-z)

      Parameter (NbrePhotosMax = 5000) 
      Parameter (ncoins = 4)
      Real*8 r_coins(2,ncoins,NbrePhotosMax)
      Dimension transf(5,NbrePhotosMax)
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
      Dimension transf(5,NbrePhotosMax)
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

      Write(6,*) 'Reading corners coordinates for each DIC file...'
      Do i = 1,NbrePhotos
         Open(1,file=FichDIC(i),status='old')
         Read(1,*) nligne
         ncn = 0
         Do k = 1,nligne
            Read(1,*) xini,yini,rot,rad,dx,dy
            irad = int(rad)
            If (irad.EQ.1) Then
               ncn = ncn + 1
               xd = xini + dx
               yd = yini + dy
               xu = xd
               yu = yd
               Call UnDisto(DistoParam,xd,yd,xu,yu,idisto)
               r_coins(1,ncn,i) = xu
               r_coins(2,ncn,i) = yu
            EndIf
         EndDo
         Close(1)
      EndDo

      Return
      End
c=======================================================================

      Subroutine ComputeDefHPP(NbrePhotos,r_coins,Fxx,Fxy,Fyx,Fyy,
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
      xmean = 0.d0
      ymean = 0.d0
      Do i = 1,ncoins
         xmean = xmean + r_coins(1,i,num_photo)
         ymean = ymean + r_coins(2,i,num_photo)
      EndDo
      xmean = xmean / dfloat(ncoins)
      ymean = ymean / dfloat(ncoins)
c     -- the vector from the center to i1 gives (cos,sin) < 0
      Do i = 1,ncoins
         xij = r_coins(1,i,num_photo) - xmean
         yij = r_coins(2,i,num_photo) - ymean
         If ( (xij.LT.0.d0).AND.(yij.LT.0.d0) ) i1 = i
         If ( (xij.GT.0.d0).AND.(yij.LT.0.d0) ) i2 = i
         If ( (xij.GT.0.d0).AND.(yij.GT.0.d0) ) i3 = i
         If ( (xij.LT.0.d0).AND.(yij.GT.0.d0) ) i4 = i
      EndDo

c-----Forme intiale du trapeze et determinant, num_photo = 1
      hxx0 = r_coins(1,i2,num_photo) - r_coins(1,i1,num_photo)
      hyx0 = r_coins(2,i2,num_photo) - r_coins(2,i1,num_photo)
      hxy0 = r_coins(1,i4,num_photo) - r_coins(1,i1,num_photo)
      hyy0 = r_coins(2,i4,num_photo) - r_coins(2,i1,num_photo)
      deth = hxx0 * hyy0 - hyx0 * hxy0

c-----Boucle sur toutes les photos
      Do num_photo = 1,NbrePhotos
         hxx = r_coins(1,i2,num_photo) - r_coins(1,i1,num_photo)
         hyx = r_coins(2,i2,num_photo) - r_coins(2,i1,num_photo)
         hxy = r_coins(1,i4,num_photo) - r_coins(1,i1,num_photo)
         hyy = r_coins(2,i4,num_photo) - r_coins(2,i1,num_photo)
         Evide(num_photo) = hxx*hyy-hyx*hxy
c        -- gradient de transformation
         Fxx(num_photo) = (hxx * hyy0 - hxy * hyx0) / deth
         Fxy(num_photo) = (hxy * hxx0 - hxx * hxy0) / deth
         Fyx(num_photo) = (hyx * hyy0 - hyy * hyx0) / deth
         Fyy(num_photo) = (hyy * hxx0 - hyx * hxy0) / deth
c        -- tenseur de def (Green-Lagrange)
c         Exx(num_photo) = 0.5d0 * (Fxx(num_photo) * Fxx(num_photo) +
c     &                             Fyx(num_photo) * Fyx(num_photo) - 1.d0)
c         Exy(num_photo) = 0.5d0 * (Fxx(num_photo) * Fxy(num_photo) +
c     &                             Fyy(num_photo) * Fyx(num_photo))
c         Eyy(num_photo) = 0.5d0 * (Fyy(num_photo) * Fyy(num_photo) +
c     &                             Fxy(num_photo) * Fxy(num_photo) - 1.d0)
c        -- tenseur de def (HPP)
         Exx(num_photo) = 0.5d0*(Fxx(num_photo)+Fxx(num_photo))-1.d0
         Exy(num_photo) = 0.5d0*(Fxy(num_photo)+Fyx(num_photo))
         Eyy(num_photo) = 0.5d0*(Fyy(num_photo)+Fyy(num_photo))-1.d0
c        -- deformation volumique
         Evol(num_photo) = Exx(num_photo) + Eyy(num_photo)
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
               Call UnDisto(DistoParam,xd,yd,xu,yu,idisto)
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

      Open(1,file=FileName,status='old')
      Read(1,*) nligne
      ngrains = 0
      Do i = 1,nligne
         Read(1,*) xini,yini,rotini,radii,dx,dy,drot,z,z,z,z,
     &             zncc1,zncc2
         If (radii.GT.3.d0) Then ! c'est un grain
            ngrains = ngrains + 1
            xd = xini + dx
            yd = yini + dy
            xu = xd
            yu = yd
            Call UnDisto(DistoParam,xd,yd,xu,yu,idisto)
            r_grains(1,ngrains) = xu
            r_grains(2,ngrains) = yu
            r_grains(3,ngrains) = radii
            r_grains(4,ngrains) = rotini + drot
            r_grains(5,ngrains) = zncc1
            r_grains(6,ngrains) = zncc2
         EndIf
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
      Character*24 cdisto
      Logical lgYES

c-----Input DIC files
      Write(6,*) '-------------------------------'
      Write(6,*) ' DIC files to read '
      Write(6,*) ' - - - - - - - - - - - - - - -'
      Write(6,*) ' Ex: '
      Write(6,*) char(8)//'../5-GrainTracking/dic_out_%04d.txt'
      Read(5,'(a50)') NomFichDIC1 
      Write(6,*) 'NUMBERS: First, Last dic file number '//
     &           ' --> they must exist!' 
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
      Write(6,*) 'Total number of DIC files: ', NbrePhotos 

c-----Parameters needed to undistor the coordinates
      Write(6,*) ' - - - - - - - - - - - - - - -'
      fichTMP = 'corrDisto.log'
      Inquire(file=fichTMP,exist=lgYES)
      If (.NOT.lgYES) Then
         Write(6,*) 'please copy ''corrDisto.log'' file'//
     &              ' in this folder'
         Write(6,*) 'hit ctrl-C or Enter to continue without'//
     &              ' correction of the disto'
         Read(5,*)
	 idisto = 2
	 Do i = 1,8
	    DistoParam(i) = 0.d0
	 EndDo   
	 Return !on sort, post-process sans correction disto
      EndIf
c      
      Open(1,file=fichTMP,status='old')
      Read(1,*) cdisto,idisto
      Write(6,*) ' --> fake_undistor = ', idisto    
      Write(6,*) ' - - - - - - - - - - - - - - -'
      Read(1,*) cdisto
      Read(1,*) cdisto ! on saute 3 lignes
      Read(1,*) cdisto
      Read(1,*) cdisto,DistoParam(1)
      Read(1,*) cdisto,DistoParam(2)
      Read(1,*) cdisto,DistoParam(3)
      Read(1,*) cdisto,DistoParam(4)
      Read(1,*) cdisto,DistoParam(5)
      Read(1,*) cdisto,DistoParam(6)
      Read(1,*) cdisto,DistoParam(7)
      Read(1,*) cdisto,DistoParam(8)
      Write(6,*) (DistoParam(i),i=1,8)

      Return
      End
c=======================================================================

      Subroutine MeanStdev(n_seg,seg1,rmean,stdev)
 
      Implicit Real*8 (a-h,o-z)
      
      Dimension seg1(50000)

      rmean = 0.d0
      Do i = 1,n_seg
         rmean = rmean + seg1(i)
      EndDo
      rmean = rmean / dfloat(n_seg)

      stdev = 0.d0
      Do i = 1,n_seg
         stdev = stdev + (seg1(i)-rmean)**2.d0
      EndDo
      stdev = sqrt(stdev/dfloat(n_seg))

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
      Write(6,*) 'ver du 2 octobre 2021, GCombe'
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
