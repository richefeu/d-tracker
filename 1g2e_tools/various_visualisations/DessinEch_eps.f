c23456789012345678901234567890123456789012345678901234567890123456789012
c=======================================================================
      Program DessinPS

      Implicit Real*4 (a-h,o,z)
      parameter(npam=100000,nvois=4*npam)
      Real*4 x(2,npam),xr(npam),hxx,hxy,hyx,hyy
      Real*4 ior(nvois),iex(nvois),fn(nvois)
      Integer npoint,npa,ncont,ntotal,nl,nlt
      Real*4 px(2),py(2)
      Real*4 s(2,npam),r(2,npam)

      Character*20 Fichier

      Write(6,*) 'Fichier ?'
      Read(5,*) Fichier

      Open(1,file=Fichier,status='old')
      Read(1,*) npa
      Write(6,*)  'npa = ',npa
      Do i = 1,npa
         Read(1,*) x(1,i),x(2,i),xr(i)
      EndDo

      Close(1)
      xmin = 0.
      xmax = 0.
      ymin = 0.
      ymax = 0.
      Do i = 1,npa
         xmin = min(xmin,x(1,i)-xr(i))
         xmax = max(xmax,x(1,i)+xr(i))
         ymin = min(ymin,x(2,i)-xr(i))
         ymax = max(ymax,x(2,i)+xr(i))
      EndDo

c-----Constants
      xdil=0.75
      sca0=2200.
      id = 1

c-----Scale positions 
      xsca=xdil*sca0/(xmax-xmin)
      ysca=xdil*sca0/(ymax-ymin)
      sca=min(xsca,ysca)
      xoff=-sca*xmin
      yoff=-sca*ymin

      Call psinit(id)

c-----draw circles
      write(id,*) '0 setlinewidth'
      g = 0.65
      Do i = 1,npa
        rxi = x(1,i)*sca + xoff
         ryi = x(2,i)*sca + yoff
         rdi = xr(i)*sca
	 write(id,198) rxi,ryi,rdi,g,' C'
      EndDo


c-----Close PS file
10    call psclose(id)

196   format(4(f6.1,1x),e12.5,3(1x,f8.2),a4)
198   format(3(f9.4,1x),f9.4,a3)
 98   format(a10)
 
      End

c=======================================================================
      subroutine psclose(id)

      write(id,*) 'showpage'
      write(id,*) '%%Trailer '
      write(id,*) '%%BoundingBox:  80 280 450 650'
      write(id,*) "%%Pages: 1"
      write(id,*) "%%EOF"
      close(id)

      return
      end

c=======================================================================
c     SR psinit
c     Initializes the PS file for the image

      subroutine  psinit(id)

      open(id,file='image.ps')
      write(id,96) '%!PS-Adobe-2.0'
      write(id,95) '%%Creator: SR    '
      write(id,95) '%%Title: image.ps'
      write(id,95) '%%Pages: (atend) '
      write(id,98) '%%BoundingBox:  80 280 450 650'
      write(id,95) '%%EndComments    '
      write(id,*) '1 setlinecap'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure cprint: texte centre' 
      write(id,*) '% (texte) x y cprint  :'
      write(id,*) '/cprint{moveto dup stringwidth pop 2 div neg 0'
      write(id,*) '        rmoveto show}def'
      write(id,*) '/Times-Roman findfont 15 scalefont setfont'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure C: trace un cercle de centre xy,' 
      write(id,*) '% de rayon r et niveau de gris g  :'
      write(id,*) '% x y r g C '
      write(id,*) '/C{newpath 4 1 roll 0 360 arc gsave setgray fill'
      write(id,*) '   grestore stroke}def'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Ca: i.e. C plus un rayon' 
      write(id,*) '% trace selon l''angle t'
      write(id,*) '% t x y r g Ca'
      write(id,*) '/Ca{4 1 roll 3 copy 7 -1 roll C'
      write(id,*) '   /r exch def newpath moveto dup cos r mul' 
      write(id,*) '   exch sin r mul' 
      write(id,*) '   rlineto stroke} def'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Lic: ligne intercentre de largeur donnee'
      write(id,*) '% largeur t, du point 0 au point 1'
      write(id,*) '% x0 y0 x1 y1 t Lic'
      write(id,*) '/Lic{gsave setrgbcolor 0.1 mul setlinewidth newpath m 
     &oveto lineto' 
      write(id,*) '    stroke grestore}def'   
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Lis: ligne de largeur fixe'
      write(id,*) '% du point 0 au point 1'
      write(id,*) '% x0 y0 x1 y1 Lis'
      write(id,*) '/Lis{gsave newpath moveto lineto' 
      write(id,*) '    stroke grestore}def'   
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Ve: Vecteur'
      write(id,*) '% module argument x0 y0 Ve'
      write(id,*) '% x0 y0 depart du vecteur'
      write(id,*) '/fref{2 setlinewidth' 
      write(id,*) '      newpath 0 0 moveto 86 0 lineto stroke'
      write(id,*) '      newpath 86 0 moveto 80 -6 lineto 100 0 lineto'
      write(id,*) '      80 6 lineto closepath fill}def'
      write(id,*) '/Ve{gsave translate rotate vesca mul '
      write(id,*) '    dup scale fref grestore}def'
      write(id,*) '/vesca 5 def'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Cg: Contact glissant'
      write(id,*) '% ang x0 y0 Cg'
      write(id,*) '% x0 y0 point de contact'
      write(id,*) '% ang: angle de la normale au contact'
      write(id,*) '%  Le symbole est une bande blanche de largeur'
      write(id,*) '%  2*cgwid et d''espacement 2*cgspa'
      write(id,*) '/Cg{gsave translate rotate 0.5 dup scale'
      write(id,*) ' newpath cgspa neg dup cgwid neg moveto cgwid lineto'
      write(id,*) ' cgspa dup cgwid lineto cgwid neg lineto closepath'
      write(id,*) ' 1 setgray fill 0 setgray 0 setlinewidth newpath'
      write(id,*) ' cgspa neg dup cgwid neg moveto cgwid lineto'
      write(id,*) ' cgspa dup cgwid moveto cgwid neg lineto stroke'
      write(id,*) ' grestore} def'
      write(id,*) '/cgwid{50}def'
      write(id,*) '/cgspa{10}def'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(id,*) '% Procedure Cng: Contact non-glissant'
      write(id,*) '% ang x0 y0 Cng'
      write(id,*) '% x0 y0 point de contact'
      write(id,*) '% ang: angle de la normale au contact'
      write(id,*) '%  Le symbole est un disque noir de rayon cngsiz'
      write(id,*) '/Cng{gsave translate rotate newpath 0 0 cngsiz'
      write(id,*) ' 0 360 arc fill grestore} def'
      write(id,*) '/cngsiz{8}def'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c	write(id,97) '(Pas = ',pas,') 300 250 cprint'
      write(id,*) '100 300 translate 0.2 dup scale' 
c	write(id,*) ' 650 070 translate 0.35 dup scale 90 rotate'
      write(id,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c        write(id,*) ''

98    format(a30)
97    format(a7,f7.0,a16)
96    format(a14)
95    format(a17)

      return
      end
c==============================================================
