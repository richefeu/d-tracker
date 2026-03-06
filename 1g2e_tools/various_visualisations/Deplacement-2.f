c=======================================================================
      program Deplacement_f

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

      Write(6,*) 'Nom du premier fichier ?'
      Read(5,*) FichierCONFini
      Write(6,*) 'Nom du dernier fichier ?'
      Read(5,*) FichierCONFfin

      diamoy = 0.d0
c-----Lecture de la configuration initiale
      Open(1,File=FichierCONFini,status='old')
      Read(1,*) npa
      Do k = 1,npa
         Read(1,*) xIni(k),yIni(k),ray,rotIni(k)
         diamoy = diamoy + ray*2.
      EndDo
      diamoy = diamoy / dfloat(npa)
      Read(1,'(a50)') a50 ! on saute une ligne fichier 1
      Do k = 1,4
         Read(1,*) coinsIni(1,k),coinsIni(2,k)
      EndDo
      Close(1)

c-----Traitement de la configuration finale
      Open(1,File=FichierCONFfin,status='old')
      Read(1,*) z
      Do k = 1,npa
         Read(1,*) xFin(k),yFin(k),diam,rotFin(k)
      EndDo
      Read(1,'(a50)') a50 ! on saute une ligne fichier 2
      Do k = 1,4
         Read(1,*) coinsFin(1,k),coinsFin(2,k)
      EndDo
      Close(1)

      Write(6,*) 'Carte des dep ou des fluct ?'
      Write(6,*) ' 0 = dep. 1 = fluct'
      Read(5,*) irep

c-----deplacement max des 4 coins
      dmax = 0.
      Do k = 1,4
         dxi = coinsFin(1,k) - coinsIni(1,k)
         dyi = coinsFin(2,k) - coinsFin(2,k)
         dd = sqrt(dxi*dxi + dyi*dyi)
         dmax = max(dmax,dd)
      EndDo
      Write(6,*) 'Distance max parcouru par les coins : ',dmax

c-----fluctuation : verification avec les 4 coins. ok, ca marche !
      If (irep.GT.0) Then
         dmax = 0.
         Do k = 1,4
            XX(1) = coinsIni(1,k)
            XX(2) = coinsIni(2,k)
            Call MeanField(coinsIni,coinsFin,XX,UU)
            depX(k) = (coinsFin(1,k)-coinsIni(1,k)) - UU(1)
            depY(k) = (coinsFin(2,k)-coinsIni(2,k)) - UU(2)
            dep(k) = sqrt(depX(k)*depX(k) + depY(k)*depY(k))
            dmax = max(dmax,dep(k))
         EndDo
         If (dmax.GT.1.d-10) Then
            Write(6,*)'pb sur le calcul du champ moyen'
            Write(6,*) 'dmax = ',dmax
            Call Exit(2)
         EndIf
      EndIf
 
c-----Calcul des deplacements.
      dmax = 0.
      dmoy = 0.
      dmin = 1d10
      Do k = 1,npa
         XX(1) = xIni(k)
         XX(2) = yIni(k)
         Call MeanField(coinsIni,coinsFin,XX,UU)
         depX(k) = (xFin(k)-xIni(k)) - dfloat(irep)*UU(1)
         depY(k) = (yFin(k)-yIni(k)) - dfloat(irep)*UU(2)
         dep(k) = sqrt(depX(k)*depX(k) + depY(k)*depY(k))
         dmax = max(dmax,dep(k))
         dmin = min(dmin,dep(k))
         dmoy = dmoy + dep(k)
c         Write(77,*) UU(1),UU(2),depX(k),depY(k)
      EndDo
      dmoy = dmoy / float(npa)
      dvar = 0.
      Do k = 1,npa
         dvar = dvar + (dep(k)-dmoy)**2.
      EndDo
      dvar = dvar / float(npa)
      dvar = sqrt(dvar)

      Write(6,*) 'dmin,dmax = ',dmin,dmax
      Write(6,*) 'dmoy,dvar = ',dmoy,dvar

c     --deplacement des coins
c      Do k = 1,4
c         dx = (coinsFin(1,k)-coinsIni(1,k)) - 
c     &             (coinsIni(1,k)*ExxInc + coinsIni(2,k)*ExyInc)
c         dy = (coinsFin(2,k)-coinsIni(2,k)) - 
c     &             (coinsIni(1,k)*EyxInc + coinsIni(2,k)*EyyInc)
c         Write(6,*) k,dx,dy
c      EndDo
        


c     -- calcul de l'echelle en fonction du diam moy
c     -- deplacement moy  = diam moy
      fact = diamoy / dmoy            
      fact = 1.5*fact
      Write(6,*) 'fact = ',fact
c      Read(5,*) fact



      xmin = 1e7
      xmax = -1e7
      ymin = 1e7
      ymax = -1e7
      Do i = 1,npa
         x1(1,i) = xFin(i)
         x1(2,i) = yFin(i)
         x2(1,i) = xFin(i)+depX(i)*fact
         x2(2,i) = yFin(i)+depY(i)*fact
         xmin = min(xmin,x1(1,i)-diamoy*2.d0)
         xmax = max(xmax,x1(1,i)+diamoy*2.d0)
         ymin = min(ymin,x1(2,i)-diamoy*2.d0)
         ymax = max(ymax,x1(2,i)+diamoy*2.d0)
         Write(19,*) xFin(i),yFin(i),depX(i)*fact,depY(i)*fact
      EndDo


c-----Affichage graphique des deplacements
      IF (PGOPEN('?') .LE. 0) STOP

c      xmin = -31.4352989
c      xmax = -0.282500267
c      ymin = 0.291249275
c      ymax = 22.3056488

      CALL PGQCIR(iC1,iC2)
      NC = MAX(0,iC2-iC1+1)
c      WRITE(6,*) 'Nbre de couleur utilisees: ',NC
c      Write(6,*) IC1,IC2
      IC1 = 16
      IC2 = 150
      BRIGHT = 0.5
      CONTRA  = 1.0
c      Write(6,*) 'luminosite = ',BRIGHT
c      Write(6,*) 'Contraste  = ',CONTRA
      CALL PALETT(3, CONTRA, BRIGHT)     

      CALL PGSAVE
      CALL PGENV(xmin,xmax,ymin,ymax,1,-2)
c      Write(6,*) xmin,xmax,ymin,ymax
      
      CALL PGSLW(3)
      Do i = 1,npa
            CALL PGSAH(1,35.,0.3)
            a = (1.-0.1)/(dmoy-dmin)
            b = 0.25 - a*dmin
            CALL PGSCH(0.2*(a*dep(i)+b)) !taille de la tete de fleche
            a = (ic2-ic1) / (dmoy+dvar)
            b = ic1
            icoul = a*dep(i)+b
            If (icoul.LT.ic1) icoul = ic1
            If (icoul.GT.ic2) icoul = ic2
            Call PGSCI(icoul)
            CALL PGARRO(x1(1,i),x1(2,i),x2(1,i),x2(2,i))
      EndDo
      CALL PGUNSA
c
      CALL PGCLOS

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
