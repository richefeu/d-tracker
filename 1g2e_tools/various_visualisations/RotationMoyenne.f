      Program RotationMoyenne

      Implicit None

      Integer npam,npa
      Parameter (npam=50000)
      Real*4 r(4,npam)
      Integer i,ii,ifich,k    
      Real*4 z,Gxy,rotmoy,rotvar,rot
      Character*1 a1
      Character*50 CONFfile1,CONffile2

      CONFfile1 = '../CONF%04d'
      Write(6,*) 'File name ?'
      Write(6,*) 'ex: ',CONFfile1
      Read(5,'(a)') CONFfile1
      Write(6,*) CONFfile1

      Write(6,*) 'Nombre de fichier CONF a lire ?'
      Read(5,*) ifich
      
      Open(2,file='RotMoyVar')

      Do ii = 1,ifich
         CONFfile2 = CONFfile1
         Call NumberingString(CONFfile2,ii)
c     ---Lecture du Fichier    
         Open(1,file=CONFfile2,status='old',err=999)
         Read(1,*) npa
         k = 0
         Do i = 1,npa
            Read(1,*) z,z,z,rot !r(4,i) = inclinaison angulaire
               k = k + 1
               r(4,k) = rot
         EndDo
         npa = k
         Do i = 1,6
            Read(1,'(a1)') a1 ! on saute des lignes
         EndDo
         Read(1,*) z,z,Gxy
         Close(1)
c     ---Calcul de la rotation moyenne
         rotmoy = 0.0
         Do i = 1,npa
            rotmoy = rotmoy + r(4,i)
         EndDo
         rotmoy = rotmoy / float(npa)
         rotvar = 0.0   
         Do i = 1,npa
            rotvar = rotvar + (r(4,i)-rotmoy)**2.
         EndDo
         rotvar = rotvar / float(npa)
         rotvar = sqrt(rotvar)
c     ---Ecriture des resultats dans un fichier
         Write(2,*) Gxy,rotmoy,rotvar         
      EndDo

      Stop

999   Write(6,*) 'Fichier ',CONFfile2
      Write(6,*) ' inexistant !'
      Stop

      End
