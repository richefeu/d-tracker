c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c  This program compute pearson correlation coefficient between
c  two set of fluctuation of displacement
c
c  It uses BilinearInterpol and pearsn library
c          
c
c  input
c  iDICdeb: num of the firt file of positions of grains
c  iDICfin: num of the last file of positions of grains
c  iINC : increment between files of positions
c
c
c  Create a file called 'fort.77'. See line 224 for details
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
      Dimension SQR_DEP(npam)
      Real*4 xIni(npam),yIni(npam)
      Real*4 xFin(npam),yFin(npam)
      Real*4 depX1(npam),depY1(npam),dep1(npam)
      Real*4 depX2(npam),depY2(npam),dep2(npam)
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
      Character*50 FichDIC,FichConf
      Character*50 NomFichDIC
      Dimension FichDIC(NbrePhotosMax)
      Dimension FichConf(NbrePhotosMax)
      Dimension XX1(npam),YY1(npam),XX2(npam),YY2(npam)
      Dimension xv1(npam),yv1(npam),xv2(npam),yv2(npam)
      Dimension PROB_TAB(npam)

      Write(6,*) '--------------------------'
      Write(6,*) 'Numero des fichiers a lire :'
      Write(6,*) ' NOM = dic_out_'
      Write(6,*) ' - - - - -'
      Write(6,*) 'NUMERO : Premier, Dernier, increment'
      Read(5,*) iDICdeb,iDICfin,iINC
      NomFichDIC = 'CONF'
      n1 = iCharlen(NomFichDIC)
      PATH = '../../../DEPOUILLE/'
      npath = iCharlen(PATH)

      k = 0
      Do i = iDICdeb,iDICfin
         k = k + 1
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
      EndDo
      nfich = k
      Write(6,*) 'Nombre de couple de fichier = ',nfich

      open(77,file='fort')

      Do i = 1,nfich-iINC-1
c     ---Premier couple
         FichierCONFini = FichDIC(i)
         n1 = iCharLen(FichierCONFini)
         FichierCONFfin = FichDIC(i+iINC)
         n2 = iCharLen(FichierCONFfin)
         Write(6,*) FichierCONFini(1:n1),' ->',FichierCONFfin(1:n2)
c        --Lecture de la configuration initiale
         Open(1,File=FichierCONFini,status='old',err=991)
         Read(1,*) npa
         Do k = 1,npa
            Read(1,*) xIni(k),yIni(k),diam,rotIni(k)
            XX1(k) = xIni(k)
            YY1(k) = yIni(k)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Do j = 1,4
            Read(1,*) coinsIni(1,j),coinsIni(2,j)
         EndDo
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
         Close(1)
c        --Calcul des fluctuations de deplacements.
         Do k = 1,npa
            XX(1) = xIni(k)
            XX(2) = yIni(k)
            Call MeanField(coinsIni,coinsFin,XX,UU)
            depX1(k) = (xFin(k)-xIni(k)) - UU(1)
            depY1(k) = (yFin(k)-yIni(k)) - UU(2)
            dep1(k) = sqrt(depX1(k)*depX1(k) + depY1(k)*depY1(k))
c            Write(66,*) xIni(k),yIni(k),depX1(k),depY1(k)
         EndDo
c         stop
c     ---Second couple
         FichierCONFini = FichDIC(i+1)
         n1 = iCharLen(FichierCONFini)
         FichierCONFfin = FichDIC(i+iINC+1)
         n2 = iCharLen(FichierCONFfin)
         Write(6,*) FichierCONFini(1:n1),' ->',FichierCONFfin(1:n2)
c        --Lecture de la configuration initiale
         Open(1,File=FichierCONFini,status='old',err=991)
         Read(1,*) npa
         Do k = 1,npa
            Read(1,*) xIni(k),yIni(k),diam,rotIni(k)
            XX2(k) = xIni(k)
            YY2(k) = yIni(k)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Do j = 1,4
            Read(1,*) coinsIni(1,j),coinsIni(2,j)
         EndDo
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
         Close(1)
c        --Calcul des fluctuations de deplacements.
         Do k = 1,npa
            XX(1) = xIni(k)
            XX(2) = yIni(k)
            Call MeanField(coinsIni,coinsFin,XX,UU)
            depX2(k) = (xFin(k)-xIni(k)) - UU(1)
            depY2(k) = (yFin(k)-yIni(k)) - UU(2)
            dep2(k) = sqrt(depX2(k)*depX2(k) + depY2(k)*depY2(k))
         EndDo
c
c      --Calcul des moindres carrées  
c         Do k = 1,npa
c            SQR_DEP(k) = (depX2(k)-depX1(k))**2.
c         EndDo
c         Call moment(SQR_DEP,npa,rmoyX,z,stdevX,z,z,z)
c         Do k = 1,npa
c            SQR_DEP(k) = (depY2(k)-depY1(k))**2.
c         EndDo
c         Call moment(SQR_DEP,npa,rmoyY,z,stdevY,z,z,z)

c         Do k = 1,npa
c            YY1(k) = depX1(k)
c            YY2(k) = depX2(k)
c         EndDo

c--------projection des depX suivant X ou Y
         Do k = 1,npa
            xv1(k) = XX1(k)
            yv1(k) = depX1(k)
            xv2(k) = XX2(k)
            yv2(k) = depX2(k)
         EndDo
         d = 0
         prod1 = 0
         Call ks2d2s(xv1,yv1,npa,xv2,yv2,npa,d,prod1)
         Do k = 1,npa
            xv1(k) = YY1(k)
            yv1(k) = depX1(k)
            xv2(k) = YY2(k)
            yv2(k) = depX2(k)
         EndDo
         d = 0
         prod2 = 0
         Call ks2d2s(xv1,yv1,npa,xv2,yv2,npa,d,prod2)
c--------projection des depY suivant X ou Y
         Do k = 1,npa
            xv1(k) = XX1(k)
            yv1(k) = depY1(k)
            xv2(k) = XX2(k)
            yv2(k) = depY2(k)
         EndDo
         d = 0
         prod3 = 0
         Call ks2d2s(xv1,yv1,npa,xv2,yv2,npa,d,prod3)
         Do k = 1,npa
            xv1(k) = YY1(k)
            yv1(k) = depY1(k)
            xv2(k) = YY2(k)
            yv2(k) = depY2(k)
         EndDo
         d = 0
         prod4 = 0
         Call ks2d2s(xv1,yv1,npa,xv2,yv2,npa,d,prod4)

         Write(6,*) i,prod1,prod2,prod3,prod4
         Write(77,'(i4,6(1x,e12.5))') i,prod1,prod2,prod3,prod4
      EndDo
      Close(77)

c      inum = (nfich-iINC-1)/2
c      Call four1(PROB_TAB,inum,1)
c      Do i = 1,2*inum
c         Write(55,*) i,PROB_TAB(i)
c      EndDo

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
      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
      Implicit Real*8 (a-h,o-z)
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
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

      SUBROUTINE ks2d2s(x1,y1,n1,x2,y2,n2,d,prob)
      INTEGER n1,n2
      REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
CU    USES pearsn,probks,quadct
      INTEGER j
      REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen,probks
      d1=0.0
      do 11 j=1,n1
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x1(j),y1(j),x2,y2,n2,ga,gb,gc,gd)
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
11    continue
      d2=0.0
      do 12 j=1,n2
        call quadct(x2(j),y2(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x2(j),y2(j),x2,y2,n2,ga,gb,gc,gd)
        d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
12    continue
      d=0.5*(d1+d2)
      sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
      call pearsn(x1,y1,n1,r1,dum,dumm)
      call pearsn(x2,y2,n2,r2,dum,dumm)
      rr=sqrt(1.0-0.5*(r1**2+r2**2))
      prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
c======================================================================
      FUNCTION probks(alam)
      REAL probks,alam,EPS1,EPS2
      PARAMETER (EPS1=0.001, EPS2=1.e-8)
      INTEGER j
      REAL a2,fac,term,termbf
      a2=-2.*alam**2
      fac=2.
      probks=0.
      termbf=0.
      do 11 j=1,100
        term=fac*exp(a2*j**2)
        probks=probks+term
        if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
        fac=-fac
        termbf=abs(term)
11    continue
      probks=1.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
c======================================================================
      SUBROUTINE quadct(x,y,xx,yy,nn,fa,fb,fc,fd)
      INTEGER nn
      REAL fa,fb,fc,fd,x,y,xx(nn),yy(nn)
      INTEGER k,na,nb,nc,nd
      REAL ff
      na=0
      nb=0
      nc=0
      nd=0
      do 11 k=1,nn
        if(yy(k).gt.y)then
          if(xx(k).gt.x)then
            na=na+1
          else
            nb=nb+1
          endif
        else
          if(xx(k).gt.x)then
            nd=nd+1
          else
            nc=nc+1
          endif
        endif
11    continue
      ff=1.0/nn
      fa=ff*na
      fb=ff*nb
      fc=ff*nc
      fd=ff*nd
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
c======================================================================
c======================================================================

