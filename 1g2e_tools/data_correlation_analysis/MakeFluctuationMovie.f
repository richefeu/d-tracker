c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c  This program make a png file of fluctuation map
c
c  It uses BilinearInterpol and pgplot library
c          
c
c  input
c  iDICdeb: num of the firt file of positions of grains
c  iDICfin: num of the last file of positions of grains
c  iINC : increment between files of positions
c
c
c
c
c
c  vincent.richefeu@hmg.inpg.fr
c  gael.combe@hmg.inpg.fr (*)
c
c  Janvier 2012
c=======================================================================
      program TemporalFlucDispCorrelation_f

      Parameter (NbrePhotosMax = 5000) 
      parameter(npam=30000,nlm=4*npam)
      Character*50 FichierCONFini,FichierCONFfin
      Integer PGOPEN
      Real*4 xIni(npam),yIni(npam)
      Real*4 xFin(npam),yFin(npam)
      Real*4 depX(npam),depY(npam),dep(npam)
      Real*4 x1(2,npam),x2(2,npam)
      Real*8 RotIni(npam),RotFin(npam)
      Real*8 errIni(npam),errFin(npam)
      Real*8 coinsIni(2,4),coinsFin(2,4)
      Real*8 XX(2),UU(2)
      Character*50 a50
      Character*50 PATH
      Character*1 a1num
      Character*2 a2num
      Character*3 a3num
      Character*4 a4num
      Character*6 PATHfig
      Character*7 chartau
      Character*50 FichDIC,FichConf
      Character*50 NomFichDIC
      Character*22 FigName
      Dimension FichDIC(NbrePhotosMax)
      Dimension FichConf(NbrePhotosMax)
      Dimension FigName(NbrePhotosMax)

      Write(6,*) '--------------------------'
      Write(6,*) 'Numero des fichiers a lire :'
      Write(6,*) ' NOM = CONF'
      Write(6,*) ' - - - - -'
      Write(6,*) 'NUMERO : Premier, Dernier, increment'
      Read(5,*) iDICdeb,iDICfin,iINC
      NomFichDIC = 'CONF'
      n1 = iCharlen(NomFichDIC)
      PATH = './'
      PATHfig = './FIG/'
      npath = iCharlen(PATH)

      k = 0
      Do i = iDICdeb,iDICfin
         k = k + 1
c        -- files to be readed
         If (i.LT.10) Then
            Write(a1num,'(i1)') i
            FichConf(i) ='CONF000'//a1num
         EndIf
         If ((i.GE.10).AND.(i.LT.100)) Then
            Write(a2num,'(i2)') i
            FichConf(i) ='CONF00'//a2num
         EndIf
         If ((i.GE.100).AND.(i.LT.1000)) Then
            Write(a3num,'(i3)') i
            FichConf(i) ='CONF0'//a3num
         EndIf
         If ((i.GE.1000).AND.(i.LT.10000)) Then
            Write(a4num,'(i4)') i
            FichConf(i) ='CONF'//a4num
         EndIf
         FichDIC(k) = PATH(1:npath)//FichConf(i)
         Write(6,*) FichDIC(k)
c        -- images to be created
         If (k.LT.10) Then
            Write(a1num,'(i1)') k
            FigName(k) = PATHfig//'fig_000'//a1num//'.png/PNG'
         EndIf
         If ((k.GE.10).AND.(k.LT.100)) Then
            Write(a2num,'(i2)') k
            FigName(k) = PATHfig//'fig_00'//a2num//'.png/PNG'
         EndIf
         If ((i.GE.100).AND.(i.LT.1000)) Then
            Write(a3num,'(i3)') k
            FigName(k) = PATHfig//'fig_0'//a3num//'.png/PNG'
         EndIf
         If ((i.GE.1000).AND.(i.LT.10000)) Then
            Write(a4num,'(i4)') k
            FigName(k) = PATHfig//'fig_'//a4num//'.png/PNG'
         EndIf
         Write(6,*) FigName(k)
      EndDo
      nfich = k
      Write(6,*) 'Nombre de couple de fichier = ',nfich 
c
      Write(6,*) 'Facteur d''echelle pour le dessin ?'
      ReaD(5,*) fact

      Do i = 1,nfich-iINC-1
c     ---Premier couple
         FichierCONFini = FichDIC(i)
         n1 = iCharLen(FichierCONFini)
         FichierCONFfin = FichDIC(i+iINC)
         n2 = iCharLen(FichierCONFfin)
         Write(6,*) FichierCONFini(1:n1),' ->',FichierCONFfin(1:n2)
c        --Lecture de la configuration initiale
         diamoy = 0.d0
         Open(1,File=FichierCONFini,status='old',err=991)
         Read(1,*) npa
         Do k = 1,npa
            Read(1,*) xIni(k),yIni(k),diam,rotIni(k)
            diamoy = diamoy + diam
         EndDo
         diamoy = diamoy / float(npa)
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Do j = 1,4
            Read(1,*) coinsIni(1,j),coinsIni(2,j)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Read(1,*) z,z,CisINI
         Close(1)
c        --Lecture de la configuration finale
         Open(1,File=FichierCONFfin,status='old',err=992)
         Read(1,*) z
         Do k = 1,npa
            Read(1,*) xFin(k),yFin(k),diam,rotFin(k)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 2
         Do j = 1,4
            Read(1,*) coinsFin(1,j),coinsFin(2,j)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Read(1,*) z,z,CisFIN
         Close(1)
c        --Calcul de tau
         tau = abs(CisFIN-CisINI) !abs pour cisail. retour
c        --Calcul des fluctuations de deplacements.
         dmin =  1.e10
         dmax = -1.e10
         dmoy = 0.
         Do k = 1,npa
            XX(1) = xIni(k)
            XX(2) = yIni(k)
            Call MeanField(coinsIni,coinsFin,XX,UU)
            depX(k) = (xFin(k)-xIni(k)) - UU(1)
            depX(k) = depX(k) / tau / diamoy
            depY(k) = (yFin(k)-yIni(k)) - UU(2)
            depY(k) = depY(k) / tau / diamoy
            dep(k) = sqrt(depX(k)*depX(k)+depY(k)*depY(k))
            dmin = min(dmin,dep(k))
            dmax = max(dmax,dep(k))
            dmoy = dmoy + dep(k)
         EndDo
         dmoy = dmoy / float(npa)
         dvar = 0.
         Do k = 1,npa
            dvar = dvar + (dep(k)-dmoy)**2.
         EndDo
         dvar = dvar / float(npa)
         dvar = sqrt(dvar)
c         Write(6,*) 'dmin,dmax = ',dmin,dmax
c         Write(6,*) 'dmoy,dvar = ',dmoy,dvar

c--------Dessin de la carte
         xmin = -0.60
         xmax = 0.1
         ymin = -0.1
         ymax = 0.5
         Do k = 1,npa
            x1(1,k) = xIni(k)
            x1(2,k) = yIni(k)
            x2(1,k) = xIni(k)+depX(k)*fact
            x2(2,k) = yIni(k)+depY(k)*fact
         EndDo
         IF (PGOPEN(FigName(i)).LE.0) STOP

         CALL PGQCIR(iC1,iC2)
         NC = MAX(0,iC2-iC1+1)
c         WRITE(6,*) 'Nbre de couleur utilisees: ',NC
c         Write(6,*) IC1,IC2
c         IC1 = 16
c         IC2 = 150
         BRIGHT = 0.5
         CONTRA  = 1.0
c         Write(6,*) 'luminosite = ',BRIGHT
c         Write(6,*) 'Contraste  = ',CONTRA
         CALL PALETT(3, CONTRA, BRIGHT)     


         CALL PGENV(xmin,xmax,ymin,ymax,1,-2)
         CALL PGSLW(1)
         CALL PGTEXT(0.,-0.025,FichierCONFfin)
         CALL PGTEXT(0.,0.,FichierCONFini)
         Write(chartau,'(f7.4)') tau
         CALL PGTEXT(0.3,0.,'Tau = '//chartau)
         CALL PGSLW(3)
         Do k = 1,npa
            CALL PGSAH(1,35.,0.3)
            a = (1.5-0.1)/(dmax-dmin)
            b = 0.1 - a*dmin
            CALL PGSCH(a*dep(k)+b) !taille de la tete de fleche
            a = (ic2-ic1) / (dmoy+dvar)
            b = ic1
            icoul = a*dep(k)+b
            If (icoul.LT.ic1) icoul = IC1
            If (icoul.GT.ic2) icoul = IC2
            Call PGSCI(icoul)
            CALL PGARRO(x1(1,k),x1(2,k),x2(1,k),x2(2,k))
         EndDo
c         CALL PGUNSA
         CALL PGCLOS
      EndDo

      Stop
991   Write(6,*) 'Fichier Inexistant : '
      Write(6,*) FichierCONFini
992   Write(6,*) 'Fichier Inexistant : '
      Write(6,*) FichierCONFfin

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


