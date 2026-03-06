      Program FollowZNCCvalues_f

      Parameter (NbrePhotosMax = 5000) 
      Character*50 FichDIC
      Dimension FichDIC(NbrePhotosMax)
      Character*50 NomFichDIC1,NomFichDIC2
      Dimension zncc(100000)
      Logical lgYES
      
      
c-----Input DIC files
      Write(6,*) '-------------------------------'
      Write(6,*) ' DICfiles to be readed '
      Write(6,*) ' - - - - - - - - - - - - - - -'
      Write(6,*) ' Ex: ./dic_out_%05d.txt or '//
     &           '../DIC/dic_out_%2d.txt or ...'
      Read(5,'(a50)') NomFichDIC1 
      Write(6,*) NomFichDIC1
      Write(6,*) 'NUMBERS : First, last dic file number'
      Read(5,*) iDICdeb,iDICfin
      Write(6,*) iDICdeb,iDICfin

      k = 0
      Do i = iDICdeb,iDICfin
         k = k + 1
         NomFichDIC2 = NomFichDIC1
         Call NumberingString(NomFichDIC2,i) 
         FichDIC(k) = NomFichDIC2
         lgYES = .FALSE.
         Inquire(file=FichDIC(k),exist=lgYES)
         If (.NOT.lgYES) Then
            Write(6,*) 'The file ',FichDIC(k),'does not exist'
            Write(6,*) 'Program STOPPED'
            Call Exit(1)
         EndIf
      EndDo
 
      NbrePhotos = k
      Write(6,*) 'Total number of DICfiles : ',NbrePhotos 
  
c-----Read DIC files, compute average, stdev, min,max, ZNCC
      Open(66,file='zncc.txt')
      Do num = 1,NbrePhotos
         zncc_min = 1e10
	 zncc_max = -1e10
	 zncc_average = 0.0
	 zncc_stdev = 0.0
         Open(1,file=FichDIC(num),status='old')
	 Read(1,*) npt
	 Do i = 1,npt
	    Read(1,*) z,z,z,z,z,z,z,z,z,z,z,z,zncc(i) 
	    zncc_min = min(zncc_min,zncc(i))
	    zncc_max = max(zncc_max,zncc(i))
	    zncc_average = zncc_average + zncc(i)
	 EndDo
	 Close(1)
	 zncc_average = zncc_average / float(npt)
	 Do i = 1,npt
	    zncc_stdev = zncc_stdev + (zncc(i) - zncc_average)**2.0
	 EndDo
	 zncc_stdev = zncc_stdev / float(npt)
	 zncc_stdev = sqrt(zncc_stdev)
	 nval9 = 0
	 nval8 = 0
	 Do i = 1,npt
	    If (zncc(i).LT.0.9) nval9 = nval9 + 1  
	    If (zncc(i).LT.0.8) nval8 = nval8 + 1
	 EndDo       
	 Write(66,*) num,zncc_min,zncc_max,zncc_average,zncc_stdev,
     &	             nval9,nval8,FichDIC(num)
1        Format(i4,4(1x,e12.5),2(1x,i4),1x,a50)
      EndDo		               
      
      End
  
      
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c
c  Subroutine NumberingString(string,i)
c
c     ex1: string = "qwerty_%5d.txt", i = 1245
c         -> string = "qwerty_1245.txt"
c     ex2: string = "qwerty_%05d.txt", i = 1245
c         -> string = "qwerty_01245.txt"  
c
c  
c  Function CharLen(string)
c
c     ex1: string = "qwerty"
c          -> CharLen = 6
c
c
c
c
c
c
c=======================================================================
      Subroutine NumberingString(String,Num)

      Character*50 String
      Character*50 Debut,Fin
      Integer Num
      Character*1 a1
      Character*10 a9
      Integer CharLen

      a9 = "          " ! 10 espaces 
      Do kk = 1,50
         Debut(kk:kk) = " "
         Fin(kk:kk) = " "
      EndDo
      k = CharLen(String)
c-----look for % format
      i = 0
1     i = i + 1
      a1 = String(i:i)
      If (a1.NE."%") GoTo 1
c-----deux cas possible : %xd ou %0xd
      j1 = i + 2
      j2 = i + 3
      a1 = String(j1:j1)
      If (a1.EQ."d") Then ! numerotation sans 0
         Debut = String(1:i-1)
         l1 = CharLen(Debut)
         Fin = String(j1+1:k)
         l3 = CharLen(Fin)
         If (Num.LT.10) 
     &      Write(a9(1:1),'(i1)') Num
         If ((Num.GE.10).AND.(Num.LT.100)) 
     &      Write(a9(1:2),'(i2)') Num
         If ((Num.GE.100).AND.(Num.LT.1000)) 
     &      Write(a9(1:3),'(i3)') Num
         If ((Num.GE.1000).AND.(Num.LT.10000)) 
     &      Write(a9(1:4),'(i4)') Num
         If ((Num.GE.10000).AND.(Num.LT.100000)) 
     &      Write(a9(1:5),'(i5)') Num
         If ((Num.GE.100000).AND.(Num.LT.1000000)) 
     &      Write(a9(1:6),'(i6)') Num
         If ((Num.GE.1000000).AND.(Num.LT.10000000)) 
     &      Write(a9(1:7),'(i7)') Num
         If ((Num.GE.10000000).AND.(Num.LT.100000000)) 
     &      Write(a9(1:8),'(i8)') Num
         If ((Num.GE.100000000).AND.(Num.LT.1000000000)) 
     &      Write(a9(1:9),'(i9)') Num
         l2 = CharLen(a9)
         String = Debut(1:l1)//a9(1:l2)//Fin(1:l3)
         GoTo 9
      EndIf
      a1 = String(j2:j2)
      If (a1.EQ."d") Then ! numerotation avec 0
         Debut = String(1:i-1)
         l1 = CharLen(Debut)
         Fin = String(j2+1:k)
         l3 = CharLen(Fin)
         k = j2 - 1 ! nombre de digits
         Read(String(k:k),*) j
         If (j.EQ.1) Write(a9(1:j),'(i1.1)') Num
         If (j.EQ.2) Write(a9(1:j),'(i2.2)') Num
         If (j.EQ.3) Write(a9(1:j),'(i3.3)') Num
         If (j.EQ.4) Write(a9(1:j),'(i4.4)') Num
         If (j.EQ.5) Write(a9(1:j),'(i5.5)') Num
         If (j.EQ.6) Write(a9(1:j),'(i6.6)') Num
         If (j.EQ.7) Write(a9(1:j),'(i7.7)') Num
         If (j.EQ.8) Write(a9(1:j),'(i8.8)') Num
         If (j.EQ.9) Write(a9(1:j),'(i9.9)') Num
         String = Debut(1:l1)//a9(1:j)//Fin(1:l3)
         GoTo 9
      EndIf
      Write(6,*) "Big PB : %xd ou %0xd non trouve"
      Stop

9     Continue
      Return
      End                  

c======================================================================
      Integer Function CharLen(text)
      
      Integer i
      Character*(*) text
      Character*1 cc

      i = 0
1     i = i + 1
      cc = text(i:i)
      If (cc.NE.'') GoTo 1
      
      CharLen = i - 1

      Return
      End

c=======================================================================
