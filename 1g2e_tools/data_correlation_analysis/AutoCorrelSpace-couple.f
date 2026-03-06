c23456789012345678901234567890123456789012345678901234567890123456789012
c        1         2         3         4         5         6         7
c=======================================================================
c
c Auto Correlation spatiale entre deux configurations
c
c
c - Lecture des deux config 
c - on retire les grains du bord
c - On retranche le champ moyen via les def pour obtenir les fluct
c - on fait comme si on avait des CLP
c - on calcul les coeff d'auto-correlation
c
c
c=======================================================================

      Program AutoCorrelSpace_couple

      Implicit Real*8 (a-h,o-z)

      parameter(npam=8000)
      Real*8 xIni(npam),yIni(npam)
      Real*8 xFin(npam),yFin(npam)
      Real*8 depX(npam),depY(npam)
      Real*8 RotIni(npam),RotFin(npam)
      Real*8 errIni(npam),errFin(npam)
      Real*8 rcoin1(2,4),rcoin2(2,4)
      Real*8 xn(2,4)
      Dimension fluc(npam)
      Parameter(nlm=npam*npam)
      Real*8 dist(nlm),List(2,nlm)
      Parameter(nclassemax = 5000)
      Dimension CC12(nclassemax),N12(nclassemax)
      Dimension CC(nclassemax)
      Real*8 hxx,hyy,hxy,hyx
      Character*50 FichCONF1,FichCONF2
      Character*50 a50
      Logical lg(npam)
      Real*8 XX(2),UU(2)
      Character*50 FichierCONFini,FichierCONFfin
      Character*50 PATH
      Character*1 a1num
      Character*2 a2num
      Character*3 a3num
      Character*4 a4num
      Character*50 FichDIC,FichConf
      Dimension Correl(3,npam)
      Real*4 xr(npam)

c     -- fichier de sortie
      Open(55,file='check-fluct') !fichier pour tracker les cartes de fluct
      Open(65,file='auto-correl') !correlogramme


      PATH = '../'
      npath = iCharlen(PATH)

      Write(6,*) 'Fichier Conf 1 ?'
      Read(5,*) FichCONF1
      Write(6,*) 'Fichier Conf 2 ?'
      Read(5,*) FichCONF2
      
      FichierCONFini = PATH(1:npath)//FichCONF1
      FichierCONFfin = PATH(1:npath)//FichCONF2     
      Write(6,*) '--'
      Write(6,*) FichierCONFini
      Write(6,*) FichierCONFfin
      
c=====Lecture des fichiers
      Open(1,file=FichierCONFini,status='old')
      Open(2,file=FichierCONFfin,status='old')
      Read(1,*) npa1
      Read(2,*) npa2
      npa = npa1
      Write(6,*) '--'
      Write(6,*) 'npa = ',npa
      If (npa1.NE.npa2) Then
         Write(6,*) 'Ouups, Conf files seems to have '//
     &              'diff number of grains'
         Write(6,*) 'We stop here :-('
         Call Exit(2)
      EndIf
c     -- lecture des positions 
      diamoy = 0.d0
      Do k = 1,npa
         Read(1,*) xIni(k),yIni(k),diam
         xr(k) = diam
         diamoy = diamoy + diam*2.d0
      EndDo 
      diamoy = diamoy / float(npa)
      Write(6,*) '--'
      Write(6,'(a32,e10.3)') 'diametre moyen des grains (m)= ',diamoy
      Write(6,'(a21,e10.3)') 'd50 des grains (m) = ',Compute_d50(npa,xr)
      Do k = 1,npa
         Read(2,*) xFin(k),yFin(k)
      EndDo
c     -- lecture des coord des coins
      Read(1,'(a50)') a50 ! on saute une ligne fichier 1
      Read(2,'(a50)') a50 ! on saute une ligne fichier 2
      Do k = 1,4
         Read(1,*) rcoin1(1,k),rcoin1(2,k)
         Read(2,*) rcoin2(1,k),rcoin2(2,k)
      EndDo
c     -- lecture des def
      Read(1,'(a50)') a50 ! on saute une ligne fichier 1
      Read(2,'(a50)') a50 ! on saute une ligne fichier 2
      Read(1,*) Exx1,Eyy1,Exy1,Eyx1,z,z,z,Erot1
      Read(2,*) Exx2,Eyy2,Exy2,Eyx2,z,z,z,Erot2
      Close(1)
      Close(2)

c     -- deplacement max des 4 coins
      dmax = 0.d0
      Do k = 1,4
         dxi = rcoin2(1,k) - rcoin1(1,k)
         dyi = rcoin2(2,k) - rcoin1(2,k)
         dd = sqrt(dxi*dxi + dyi*dyi)
         dmax = max(dmax,dd)
      EndDo
      Write(6,*) '--'
      Write(6,'(a41,e10.3)') 'Distance max parcouru '//
     &                       'par les coins (m): ',dmax

c     -- incr de def entre les deux config
      ExxInc = Exx2 - Exx1
      EyyInc = Eyy2 - Eyy1
      ExyInc = Exy2 - Exy1
      EyxInc = Eyx2 - Eyx1
      dRotGamma = - (Erot2 - Erot1) ! <0 si sens inverse trigo
      Write(6,*) '--'
      Write(6,'(a12,2(e10.3,1x))') 'dExx,dEyy = ',ExxInc,EyyInc
      Write(6,'(a12,2(e10.3,1x))') 'dExy,dEyx = ',ExyInc,EyxInc
      Write(6,'(a12,1(e10.3,1x))') 'dRotGamma = ',dRotGamma
c
c-----calcul des fluctuations
      Do k = 1,npa
         XX(1) = xIni(k)
         XX(2) = yIni(k)
         Call MeanField(rcoin1,rcoin2,XX,UU)
         depX(k) = (xFin(k)-xIni(k)) - UU(1)
         depY(k) = (yFin(k)-yIni(k)) - UU(2)
c         Write(34,*) XX(1),XX(2),depX(k),depY(k) !fluctuations
      EndDo

c=====calcul des correlations spatiales
c     -- variable a correler
      fmin = 1e10
      fmax = -1e10
      Do i = 1,npa
         fluc(i) = depX(i)
c         fluc(i) = sqrt(depX(i)*depX(i)+depY(i)*depY(i))
c         fluc(i) = atan(tan(RotFin(i)-RotIni(i))) - dRotGamma
         fmin = min(fmin,fluc(i))
         fmax = max(fmax,fluc(i))
      EndDo
      Write(6,'(a21,2(e10.3,1x))') 'fluc_min,fluc_max = ',fmin,fmax
c     -- on ecrete la variable a correler
      Write(6,*) 'fluc_min, fluc_max ?'
      Read(5,*) fmin,fmax

c     -- valeur min et max des distances inte-grains (fichier 1)
      dmin = 1.d10
      dmax = 0.d0
      nval = 0
      factvit = 1.! 200.d0 utile pour un sur-echantillonnage
      Do i = 1,npa
         Do j = 1,npa
            dx = xIni(j) - xIni(i)
            dy = yIni(j) - yIni(i)
            dij2 = sqrt(dx*dx + dy*dy)/diamoy*factvit
            nval = nval + 1
            dist(nval) = (dij2)
            dmin = min(dmin,(dij2))
            dmax = max(dmax,(dij2))
            List(1,nval) = i
            List(2,nval) = j
         EndDo
      EndDo
      Write(6,*) '--'
      Write(6,*) 'Distance intergrain normalisee (diametre moyen):'
      Write(6,'(a12,2(e10.3,1x))') 'dmin,dmax = ',dmin,dmax
      Write(6,*) '--'
      Write(6,'(a27,i8)') 'nbre de couple de grains = ',nval
c     -- nbre de classes  = nbre de diam sur la diagonale
      nclasse = int(dmax)+1
      Write(6,*) '--'
      Write(6,*) 'Nbre de classes pour les correlations = ',nclasse

c     -- moyenne des fluctuations
      rmoy = 0.d0
      Do i = 1,npa
         rmoy = rmoy + fluc(i)
      EndDo
      rmoy = rmoy / float(nval)
      Write(6,*) '--'
      Write(6,'(a27,e10.3)') 'Moyenne des fluctuations : ',rmoy

c     -- variance (denominateur)
      ncc1 = 0
      CC1 = 0.d0
      Do i = 1,npa
         If (fluc(i).GT.fmin.AND.fluc(i).LT.fmax) Then
            ncc1 = ncc1 + 1
            CC1 = CC1 + (fluc(i)-rmoy)**2.d0
            Write(55,*) xini(i),yini(i),depX(i),depY(i)
         EndIf
      EndDo
      CC1 = CC1 / float(ncc1)
      Write(6,*) 'ncc1, CC1 = ',ncc1,CC1

c     -- covariance (numerateur)
      Do il = 1,nclassemax
         CC12(il) = 0.d0
         N12(il) = 0
      EndDo
      Do il = 1,nval
         i = List(1,il)
         j = List(2,il)
         k = int(dist(il))
         If (fluc(i).GT.fmin.AND.fluc(i).LT.fmax) Then
            CC12(k) = CC12(k) + (fluc(i)-rmoy)*(fluc(j)-rmoy)
            N12(k) = N12(k) + 1
         EndIf
      EndDo
      Write(65,*) 0,1.,1.
      Do k = 1,nclasse
         If (N12(k).GT.0) Then
            CC12(k) = CC12(k)/N12(k)
            CC(k) = CC12(k)/CC1
         Else
            CC(k) = 0.d0
         EndIf
         Write(65,*) k/factvit,CC(k),N12(k)
      EndDo

c      Write(67,*) 0,0,1.
c      kmin = 1
c      kmax = nclasse
c      Write(6,*) 'diamoy*factvit = ',factvit*diamoy
c      lx = 20
c      ncpt = 0
c      Do k = kmin+lx,kmax-lx,lx
c         m = k - lx
c         n = k + lx
c         rmoy = 0.d0
c         Do l = m,n
c            rmoy = rmoy + CC(l)
c         EndDo
c         rmoy = rmoy / dfloat(2*lx+1)
c         Write(67,*) k,dfloat(k)/factvit,rmoy
c         ncpt = ncpt + 1
c         Correl(1,ncpt) = rmoy
c         Correl(3,ncpt) = dfloat(k)/factvit
c      EndDo





c=====calcul des correlations spatiales
c     -- variable a correler
c      Do i = 1,npa
c         fluc(i) = depY(i)
c      EndDo

c     -- moyenne des fluctuations
c      rmoy = 0.d0
c      Do i = 1,npa
c         rmoy = rmoy + fluc(i)
c      EndDo
c      rmoy = rmoy / float(nval)
c      Write(6,*) 'Moyenne des fluctuations :'
c      Write(6,*) 'rmoy =  ',rmoy

c     -- variance (denominateur)
c      ncc1 = 0
c      CC1 = 0.d0
c      Do i = 1,npa
c         ncc1 = ncc1 + 1
c         CC1 = CC1 + (fluc(i)-rmoy)**2.d0
c      EndDo
c      CC1 = CC1 / float(ncc1)
c      Write(6,*) 'ncc1, CC1 = ',ncc1,CC1

c     -- covariance (numerateur)
c      Do il = 1,nclassemax
c         CC12(il) = 0.d0
c         N12(il) = 0
c      EndDo
c      Do il = 1,nval
c         i = List(1,il)
c         j = List(2,il)
c         k = int(dist(il))
c         CC12(k) = CC12(k) + (fluc(i)-rmoy)*(fluc(j)-rmoy)
c         N12(k) = N12(k) + 1
c      EndDo
c      Write(65,*) 0,1.
c      Do k = 1,nclasse
c         If (N12(k).GT.0) Then
c            CC12(k) = CC12(k)/N12(k)
c            CC(k) = CC12(k)/CC1
c         Else
c            CC(k) = 0.d0
c         EndIf
c         Write(65,*) k/factvit,CC(k)
c      EndDo

c      Write(68,*) 0,0,1.
c      kmin = 1
c      kmax = nclasse
c      Write(6,*) 'diamoy*factvit = ',factvit*diamoy
c      lx = 20
c      ncpt = 0
c      Do k = kmin+lx,kmax-lx,lx
c         m = k - lx
c         n = k + lx
c         rmoy = 0.d0
c         Do l = m,n
c            rmoy = rmoy + CC(l)
c         EndDo
c         rmoy = rmoy / dfloat(2*lx+1)
c         Write(68,*) k,dfloat(k)/factvit,rmoy
c         ncpt = ncpt + 1
c         Correl(2,ncpt) = rmoy
c      EndDo


c=====moyenne
c      Write(69,*) 0,1,1,1
c      Do i = 1,ncpt
c         Write(69,*) Correl(3,i),Correl(1,i),Correl(2,i),
c     &               (Correl(1,i)+Correl(2,i))*0.5
c      EndDo

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

      Real*8 Function Compute_d50(npa,xr)

      Integer npa
      Dimension xr(npa)
      Dimension cumul(10000)

c     -- Sort
      Call sort(npa,xr)
c     -- cumulative aera
      cumulval = 0.0
      Do i = 1,npa
         cumulval = cumulval + xr(i)*xr(i)*4.d0*3.14159
         cumul(i) = cumulval 
      EndDo
      Do i = 1,npa
         cumul(i) = cumul(i) / cumul(npa)
      EndDo
c     -- which size for 50% or the grading
      distm = 1e10
      Do i = 1,npa
         dist = abs(cumul(i)-0.5)
         distm = min(distm,dist)
         If (distm.EQ.dist) i0 = i
      EndDo

      Compute_d50 = (xr(i0)*2.d0)

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
        if(jstack.gt.NSTACK) stop
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
c=======================================================================




