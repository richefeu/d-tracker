c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c  This program create files with fluctuation velocities
c
c  It uses BilinearInterpol 
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
      Character*50 OutputFILE
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
      Dimension val(NbrePhotosMax)
      Dimension ss1(10000000)
      Dimension ss2(10000000)
      Dimension ss3(10000000)

      Write(6,*) '--------------------------'
      Write(6,*) 'Numero des fichiers a lire :'
      Write(6,*) ' NOM = dic_out_'
      Write(6,*) ' - - - - -'
      Write(6,*) 'NUMERO : Premier, Dernier, increment'
      Read(5,*) iDICdeb,iDICfin,iINC
      NomFichDIC = 'CONF'
      n1 = iCharlen(NomFichDIC)
      PATH = '../../DEPOUILLE/'
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
      EndDo
      nfich = k
      Write(6,*) 'Nombre de couple de fichier = ',nfich 
c

      Write(6,*) 'Output file name for velocities ?'
      Read(5,*) OutputFILE
      Open(99,file=OutputFILE)

      nbre = 0
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
            diamoy = diamoy + diam*2.
         EndDo
         diamoy = diamoy / float(npa)
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Do j = 1,4
            Read(1,*) coinsIni(1,j),coinsIni(2,j)
         EndDo
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Read(1,*) z,z,CisINI
         Read(1,'(a50)') a50 ! on saute une ligne fichier 1
         Read(1,*) fact_echelle
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

c        seuil de confiance dans les fluctuations 
c        precision de la mesure : 0.1 pixels
c        les fluctuations inferieures sont supprimes     

         Do k = 1,npa
            XX(1) = xIni(k)
            XX(2) = yIni(k)
            Call MeanField(coinsIni,coinsFin,XX,UU)
            dx = xFin(k)-xIni(k)
            dy = yFin(k)-yIni(k)
            depX(k) = (xFin(k)-xIni(k)) - UU(1)
            depY(k) = (yFin(k)-yIni(k)) - UU(2)
            dep(k) = sqrt(dx*dx+dy*dy)

            depX(k) = depX(k) / tau / diamoy
            depY(k) = depY(k) / tau / diamoy

            dmin = min(dmin,dep(k))
            dmax = max(dmax,dep(k))
            dmoy = dmoy + dep(k)
            fluct = dep(k) * fact_echelle
c            If (fluct.GT.0.1) Then
               Write(99,'(3(e12.5,1x))') depX(k),depY(k),fluct
c            EndIf
         EndDo
         dmoy = dmoy / float(npa)
         dvar = 0.
         Do k = 1,np
            dvar = dvar + (dep(k)-dmoy)**2.
         EndDo
         dvar = dvar / float(npa)
         dvar = sqrt(dvar)
c         Write(6,*) 'dmin,dmax = ',dmin,dmax
c         Write(6,*) 'dmoy,dvar = ',dmoy,dvar
         nbre = nbre + 1
         val(nbre) = tau
      EndDo
      
      rmoy = 0.
      Do i = 1,nbre
         rmoy = rmoy + val(i)
      EndDo
      rmoy = rmoy / float(nbre)
      rvar = 0.     
      Do i = 1,nbre
         rvar = rvar + (val(i)-rmoy)**2.
      EndDo
      rvar = rvar / float(nbre)
      rvar = sqrt(rvar)
      Write(6,*) 'Tau (moy,var) =  ',rmoy,rvar
      Open(23,file='TAU',access='append')
      Write(23,'(i3,5x,e12.5,1x,e12.5)') iInc,rmoy,rvar
      Close(23)

      Close(99)

      Open(99,file=OutputFILE)
      n = 0
1     Read(99,*,end=2) z
          n = n + 1
          GoTo 1
2     Rewind(99)
      Do i = 1,n
         Read(99,*) ss1(i),ss2(i),ss3(i)
      EndDo
      Close(99)
c      Call sort(n,ss1)
c      Call sort(n,ss2)
c      Call sort(n,ss3)

      Open(99,file=OutputFILE)
      Do i = 1,n
         Write(99,'(3(e12.5,1x))') ss1(i),ss2(i),ss3(i)
      EndDo
      Close(99)

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
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
