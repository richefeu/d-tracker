c23456789012345678901234567890123456789012345678901234567890123456789012
c=======================================================================
      Program Fabric4Points

      Implicit Real*4 (a-h,o-z)

      Parameter(nmax=4*20000)
      Dimension r_ncoins(4,nmax)
      Dimension r_nptfix(4,nmax)
      Dimension r_ngrains(4,nmax)
      Dimension r_newgrain(4,nmax)
	  Dimension angle(nmax)
      Dimension pattern(2,nmax)
      Character*50 FichierInput

c-----Titre
      Call Title

c-----lecture du fichier initial des positions
      Write(6,*) 'Coordinates files to read ?'
      Read(5,*) FichierInput
      Open(1,file=FichierInput,status='old')
      Read(1,*) npt
      Write(6,*) ' '
      Write(6,*) 'Total number of points = ',npt
      ncoins = 0
      nptfix = 0
      ngrain = 0
      Do i = 1,npt
         Read(1,*) x,y,rot,rad
         If (rad.EQ.1) Then
            ncoins = ncoins + 1
            r_ncoins(1,ncoins) = x
            r_ncoins(2,ncoins) = y
            r_ncoins(3,ncoins) = rot
            r_ncoins(4,ncoins) = rad
         EndIf
         If ((rad.GE.2).AND.(rad.LT.3)) Then
            nptfix = nptfix + 1
            r_nptfix(1,nptfix) = x
            r_nptfix(2,nptfix) = y
            r_nptfix(3,nptfix) = rot
            r_nptfix(4,nptfix) = rad
         EndIf
         If (rad.GT.3) Then
            ngrain = ngrain + 1
            r_ngrains(1,ngrain) = x
            r_ngrains(2,ngrain) = y
            r_ngrains(3,ngrain) = rot
            r_ngrains(4,ngrain) = rad
         EndIf
      EndDo
      Close(1)
      Write(6,*) 'Numbers of points for corners = ',ncoins
      Write(6,*) 'Numbers of fixed points = ',nptfix
      Write(6,*) 'Numbers of GRAINS =',ngrain 

c-----Stat sur les rayons des points
      rmin = 1e10
      rmax = -1e10
      rmoy = 0.
      rvar = 0.
      Do i = 1,ngrain
         rad = r_ngrains(4,i)
         rmin = min(rmin,rad)
         rmax = max(rmax,rad)
         rmoy = rmoy + rad
      EndDo
      rmoy = rmoy / float(ngrain)
      Do i = 1,ngrain
         rad = r_ngrains(4,ngrain)
         rvar = rvar + (rmoy - rad)**2.
      EndDo
      rvar = sqrt(rvar / float(ngrain))
      Write(6,*) '--> Stats on points radii'
      Write(6,*) 'min = ',rmin
      Write(6,*) 'max = ',rmax
      Write(6,*) 'mean = ',rmoy
      Write(6,*) 'stdev = ',rvar
      If ((rvar.GT.rmoy).OR.(rvar.GT.rmin)) Then
         Write(6,*) 'It seems that your file contains grains'
         Write(6,*) 'of different sizes...'
         Write(6,*) 'You will probably want to change that...'
         Write(6,*) 'To do that, change the file to process'
         Write(6,*) 'Strike any key to continue anyway ...'
         Read(*,*) ! Pause
      EndIf

c-----Creation des 4 points autour du point suppose origine
c     --rayon du pattern en pixel
      Write(6,*) ' '
      Write(6,*) 'Length of the segments ? '
      Read(5,*) decal
	  
c     -- orientation aleatoire des croix
	  igr = 123
	  Pi = 4.*atan(1.)
	  smin = 0.
	  smax = Pi*0.5
	  ds = smax-smin
	  Do i = 1,ngrain
		 angle(i) = smin + ds*rand2(igr) 
		 Write(67,*) angle(i)
      EndDo	  	 
      xdecal = decal*0.5   
      ydecal = decal*0.5 
      Write(6,*) 'Give a radius for new the fictive points '
      Read(5,*) iradpatt 
c     -- text on the radius
      total_length = 2*iradpatt + decal
      If (total_length.GT.2*rmin) Then
         Write(6,*) 'Be carefull, the points associated'
         Write(6,*) 'to the segments will be out of the'
         Write(6,*) 'physical grains !'
         Write(6,*) 'Strike any key to continue anyway ...'
         Read(*,*) ! Pause
      EndIf
      
      nnewgrain = 0
      Do i = 1,ngrain   
		  cs = int(xdecal*cos(angle(i)))
		  sn = int(ydecal*sin(angle(i)))  
		  Write(6,*) cs,sn      
         nnewgrain = nnewgrain + 1
         r_newgrain(1,nnewgrain) = r_ngrains(1,i) + cs
         r_newgrain(2,nnewgrain) = r_ngrains(2,i) + sn
         r_newgrain(3,nnewgrain) = r_ngrains(3,i)
         r_newgrain(4,nnewgrain) = r_ngrains(4,i)
         nnewgrain = nnewgrain + 1
         r_newgrain(1,nnewgrain) = r_ngrains(1,i) - cs
         r_newgrain(2,nnewgrain) = r_ngrains(2,i) - sn
         r_newgrain(3,nnewgrain) = r_ngrains(3,i)
         r_newgrain(4,nnewgrain) = r_ngrains(4,i)
         nnewgrain = nnewgrain + 1
         r_newgrain(1,nnewgrain) = r_ngrains(1,i) - sn
         r_newgrain(2,nnewgrain) = r_ngrains(2,i) + cs
         r_newgrain(3,nnewgrain) = r_ngrains(3,i)
         r_newgrain(4,nnewgrain) = r_ngrains(4,i)
         nnewgrain = nnewgrain + 1
         r_newgrain(1,nnewgrain) = r_ngrains(1,i) + sn
         r_newgrain(2,nnewgrain) = r_ngrains(2,i) - cs
         r_newgrain(3,nnewgrain) = r_ngrains(3,i)
         r_newgrain(4,nnewgrain) = r_ngrains(4,i)
      EndDo

c-----Enregistrement du nouveau fichier de position des points
      ntotal = nnewgrain + ncoins + nptfix
      Write(6,*) 'Nbre total de pt = ',ntotal
      Write(6,*) 'Creation of the new coordinates files for Tracker'
      Write(6,*) '--> New4PosGrain.data'
      Open(1,file='New4PosGrain.data')
      Call Ecriture1Nbre(ntotal,1)
      Do i = 1,ncoins
         x   = r_ncoins(1,i)
         y   = r_ncoins(2,i)
         rot = r_ncoins(3,i)
         rad = r_ncoins(4,i)
         Call Ecriture4Nbre(x,y,rot,rad,1)
      EndDo
      Do i = 1,nptfix
         x   = r_nptfix(1,i)
         y   = r_nptfix(2,i)
         rot = r_nptfix(3,i)
         rad = r_nptfix(4,i)
         Call Ecriture4Nbre(x,y,rot,rad,1)
      EndDo
      Do i = 1,nnewgrain
         x   = r_newgrain(1,i)
         y   = r_newgrain(2,i)
         rot = r_newgrain(3,i)
         rad = iradpatt
         Call Ecriture4Nbre(x,y,rot,rad,1)
      EndDo
      Close(1)

c-----creation du fichier de pattern centres sur r_grain 
c                                (et non pas r_newgrain)
c

      Open(1,file='Pattern.data')
c     -- nbre total de pattern
      npattern = ntotal
      Write(6,*) 'Nbre total de pattern = ',npattern
      Call Ecriture1Nbre(npattern,1)
c     -- pattern associes aux 4 coins
      Write(6,*) 'pattern associes aux 4 coins'
      Do i = 1,ncoins
         dx = 0.
         dy = 0.
         Call MakeCircularPattern(iradpatt*2,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
      EndDo
c     -- pattern associes aux pt fixes
      Write(6,*) 'pattern associes aux pt fixes'
      Do i = 1,nptfix
         dx = 0.
         dy = 0.
         Call MakeCircularPattern(iradpatt,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
      EndDo
c     -- pattern associes aux 4 pts des grains. 
      Write(6,*) 'pattern associes aux 4 pts des grains'
      Do i = 1,nnewgrain,4
c        -- pt 1
         dx = xdecal
         dy = 0.
         Call MakeCircularPattern(iradpatt,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
c        -- pt 2
         dx = - xdecal
         dy = 0.
         Call MakeCircularPattern(iradpatt,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
c        -- pt 3
         dx = 0.
         dy = ydecal
         Call MakeCircularPattern(iradpatt,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
c        -- pt 4
         dx = 0.
         dy = - ydecal
         Call MakeCircularPattern(iradpatt,dx,dy,npatt,pattern)
         Call Ecriture1Nbre(npatt,1) ! nbre de pt du pattern
         Do k = 1,npatt
            x = pattern(1,k)
            y = pattern(2,k)
            Call Ecriture2Nbre(x,y,1)
         EndDo
      EndDo

      End
c=======================================================================
      Subroutine MakeCircularPattern(iradii,dx,dy,npatt,pattern)

      Parameter(nmax=4*20000)
      Dimension pattern(2,nmax)

c-----creation du pattern
      npatt = 0
      Do ix = -iradii,iradii
         Do iy = -iradii,iradii
            dist = sqrt(float(ix*ix+iy*iy))
            If (dist.LE.iradii) Then
               npatt = npatt + 1
               pattern(1,npatt) = ix
               pattern(2,npatt) = iy
            EndIf
         EndDo
      EndDo 
c-----décalage du pattern
      Do i = 1,npatt
         pattern(1,i) = pattern(1,i) + dx
         pattern(2,i) = pattern(2,i) + dy
      EndDo

      Return
      End
      
c=======================================================================
      Subroutine Ecriture1Nbre(ix,nfich)

      If (ix.LT.10) Then
c         Write(6,20) ix    
         Write(nfich,20) ix
      Else
         If (ix.LT.100) Then
c            Write(6,21) ix
            Write(nfich,21) ix
         Else
            If (ix.LT.1000) Then
c               Write(6,22) ix
               Write(nfich,22) ix
            Else
c               If (ix.LT.10000) Write(6,23) ix
               If (ix.LT.10000) Write(nfich,23) ix
            EndIf
         EndIf
       EndIf

20     Format(i1)
21     Format(i2)
22     Format(i3)
23     Format(i4)

       Return
       End

c=======================================================================
      Subroutine Ecriture2Nbre(x,y,nfich)

      If (abs(x).LT.10) Then
         If (x.GE.0) Then
            Write(nfich,30) int(x),int(y)
         Else
            Write(nfich,40) int(x),int(y)
         EndIf
      Else
         If (x.LT.100) Then
            If (x.GE.0) Then
               Write(nfich,31) int(x),int(y)
            Else
               Write(nfich,41) int(x),int(y)
            EndIf
         Else
            If (x.LT.1000) Then
               Write(nfich,32) int(x),int(y)
            Else
               If (x.LT.10000) Write(nfich,33) int(x),int(y)
            EndIf
         EndIf
       EndIf

30     Format(i1,4x,i4)
40     Format(i2,3x,i4)
31     Format(i2,3x,i4)
41     Format(i3,2x,i4)
32     Format(i3,2x,i4)
33     Format(i4,1x,i4)

       Return
       End

c=======================================================================
      Subroutine Ecriture4Nbre(x,y,rot,rad,nfich)

      If (x.LT.10) Then
         Write(nfich,10) int(x),int(y),rot,rad
      Else
         If (x.LT.100) Then
            Write(nfich,11) int(x),int(y),rot,rad
         Else
            If (x.LT.1000) Then
               Write(nfich,12) int(x),int(y),rot,rad
            Else
               If (x.LT.10000) Then
                  Write(nfich,13) int(x),int(y),rot,rad
               Else
                  Write(nfich,14) int(x),int(y),rot,rad
               EndIf
            EndIf
         EndIf
       EndIf

10     Format(i1,5x,i5,2(1x,e12.5))
11     Format(i2,4x,i5,2(1x,e12.5))
12     Format(i3,3x,i5,2(1x,e12.5))
13     Format(i4,2x,i5,2(1x,e12.5))
14     Format(i5,1x,i5,2(1x,e12.5))

       Return
       End
c=======================================================================

      Subroutine Title

      Write(6,*) '-----------------------------------------------'
      Write(6,*) '         __    _______  _______  _______       '
      Write(6,*) '        /  \  (  ____ \/ ___   )(  ____ \      '
      Write(6,*) '        \/) ) | (    \/\/   )  || (    \/      '
      Write(6,*) '          | | | |          /   )| (__          ' 
      Write(6,*) '          | | | | ____   _/   / |  __)         '
      Write(6,*) '          | | | | \_  ) /   _/  | (            '
      Write(6,*) '        __) (_| (___) |(   (__/\| (____/\      '
      Write(6,*) '        \____/(_______)\_______/(_______/      '
      Write(6,*) '                                               '
      Write(6,*) '==============================================='
      Write(6,*) '  _                          _                 '
      Write(6,*) ' |_ _. |_  ._ o  _   |_|_   |_) _  o ._ _|_  _ '
      Write(6,*) ' | (_| |_) |  | (_     |    |  (_) | | | |_ _> '
      Write(6,*) '                                               '
      Write(6,*) '==============================================='
      Write(6,*) 'ver 1.0 --> 2015-02-13'
      Write(6,*) ' '
      Write(6,*) '-----------------------------------------------'

      Return
      End
c======================================================================

      Real*4 Function rand2(isem)

      Implicit Double Precision (a-h,o-z)
      Parameter(iA = 843314861, iB = 453816693, iM = 1073741824)

      aux = 0.5 / float(iM)

      isem = isem*iA + iB
      If (isem < 0) Then
        isem = (isem + iM) + iM
      EndIf

      x = isem*aux

      rand2 = x
c     Write(6,*) x,aux

      Return
      End
c=======================================================================








