      program tdep
c
c-----TAD Diff Exp Plot
c-----Written by: Stephen W. Schaeffer
c-----The Pennsylvania State University
c-----Dept. of Biology
c-----Version: 1.0
c-----Date: 25 January 2022
c
c-----This fortran program generates a postscript file
c-----that plots TADs for the third chromosome of 
c-----Drosophila pseudoobscura based on data from
c-----Liao, Y., X. Zhang, M. Chakraborty and J. J. Emerson, 2021
c-----Topologically associating domains and their role in the
c-----evolution of genome structure and function in Drosophila.
c-----Genome Research 31: 397-410.
c
      integer nx,ny,nx1,nx2,ix,iy
c
c-----Open the output TAD_Diff_Expression_Plot_Raw postscript file
c-----Adobe postscript commands are from Adobe Systems Incorporated
c-----1985 Postscript Language: Tutorial and Cookbook.
c-----Addison-Wesley Publishing Company, Inc.
c
c-----This plot will show TADs on five vertical lines of 5 Megabase
c-----segments.  This is a raw plot that has added details that
c-----are integrated in Adobe Illustrator.
c
      open(unit=1,file='TAD_Diff_Expression_Plot_Raw.ps')
      write(1,'("%!")')
      write(1,'("1600 setlinewidth")') 
      write(1,'("0.0006 0.00015 scale")')
c
c-----Plot five vertical lines
c-----ix - x coordinate
c-----iy - y coordinate
c
      ix=201600
      iy=10
      do i=1,5
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,iy
      write(1,'("0 5000000 rlineto")') 
      write(1,'("stroke")')
c
c-----Increment the x coordinate
c
      ix=ix+200000
      end do
c
c-----Map the TADs across the chromosome moving from the bottom to the
c-----top of the page.  Each TAD is drawn as a triangle to the left of
c-----the line.
c
c-----nf   - factor to subtract from TAD in multiples of 5000000 
c-----ny   - maximum y coordinate
c
      nf=0
      ny=5000000
      ix=200000
      iy=10 
c
c-----The "TADs.txt" file has the list of 167 TADs for chromosome 3 
c
      open(unit=2,file='TADs.txt')
c
c-----nx1  - beginning coordinate of TAD
c-----nx2  - end coordinate of TAD
c-----ndeg - number of differentially expressed genes within the TAD
c-----ng   - total number of genes within the TADs
c
50    read(2,*,end=100) nx1,nx2,ndeg,ng
c
c-----fde - fraction of differentially expressed genes
c
      fde=float(ndeg)/float(ng)
c
c-----Check that nx1 and nx2 are less than the maximum for each
c-----vertical column.
c
c-----Column 1 goes from       0  to  5000000 bottom to top
c-----Column 2 goes from 5000000  to 10000000 bottom to top
c-----Column 3 goes from 10000000 to 15000000 bottom to top
c-----Column 4 goes from 15000000 to 20000000 bottom to top
c-----Column 5 goes from 20000000 to 25000000 bottom to top
c
      if(nx1.le.ny.and.nx2.le.ny) then
c
c-----If the fraction of differentially expressed genes fde is 0
c-----draw an unfilled triangle.
c
      if(ndeg.eq.0) then
      write(1,'("1 1 1 1 setcmykcolor")')
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,nx1+iy-nf
      write(1,'(2i11," rlineto")') -1*(nx2-nx1)/2,(nx2-nx1)/2
      write(1,'(2i11," rlineto")') (nx2-nx1)/2,(nx2-nx1)/2
      write(1,'("stroke")')
      else
c
c-----If the fraction of differentially expressed genes fde is >0
c-----draw a filled triangle with a brown color from light brown
c-----to dark brown based on the fraction of fde, low to high
c
      if(fde.gt.0.0.and.fde.le.0.1) then
      write(1,'("0.25 0.40 0.65 0.00  setcmykcolor")') 
      elseif(fde.gt.0.1.and.fde.le.0.2) then
      write(1,'("0.30 0.50 0.75 0.10  setcmykcolor")') 
      elseif(fde.gt.0.2.and.fde.le.0.3) then
      write(1,'("0.35 0.60 0.80 0.25  setcmykcolor")') 
      elseif(fde.gt.0.3.and.fde.le.0.4) then
      write(1,'("0.40 0.65 0.90 0.35  setcmykcolor")') 
      elseif(fde.gt.0.4.and.fde.le.0.5) then
      write(1,'("0.40 0.70 1.00 0.50  setcmykcolor")') 
      end if
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,nx1+iy-nf
      write(1,'(2i11," rlineto")') -1*(nx2-nx1)/2,(nx2-nx1)/2
      write(1,'(2i11," rlineto")') (nx2-nx1)/2,(nx2-nx1)/2
      write(1,'("fill")')
      end if
c
c-----If the TAD extends past the maximum, plot the last triangle
c-----and move to the next vertical column by re-initializing
c-----the x (ix) and y (iy) coordinates.  Update the base
c-----coordinate nf for the next vertical column.  Similar coloring
c-----schemes are used as above.
c
      elseif (nx1.le.ny.and.nx2.gt.ny) then
      if(ndeg.eq.0) then
      write(1,'("1 1 1 1 setcmykcolor")')
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,nx1+iy-nf
      write(1,'(2i11," rlineto")') -1*(nx2-nx1)/2,(nx2-nx1)/2
      write(1,'(2i11," rlineto")') (nx2-nx1)/2,(nx2-nx1)/2
      write(1,'("stroke")')
      else
      if(fde.gt.0.0.and.fde.le.0.1) then
      write(1,'("0.25 0.40 0.65 0.00  setcmykcolor")') 
      elseif(fde.gt.0.1.and.fde.le.0.2) then
      write(1,'("0.30 0.50 0.75 0.10  setcmykcolor")') 
      elseif(fde.gt.0.2.and.fde.le.0.3) then
      write(1,'("0.35 0.60 0.80 0.25  setcmykcolor")') 
      elseif(fde.gt.0.3.and.fde.le.0.4) then
      write(1,'("0.40 0.65 0.90 0.35  setcmykcolor")') 
      elseif(fde.gt.0.4.and.fde.le.0.5) then
      write(1,'("0.40 0.70 1.00 0.50  setcmykcolor")') 
      end if
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,nx1+iy-nf
      write(1,'(2i11," rlineto")') -1*(nx2-nx1)/2,(nx2-nx1)/2
      write(1,'(2i11," rlineto")') (nx2-nx1)/2,(nx2-nx1)/2
      write(1,'("fill")')
      end if
      ix=ix+200000
      iy=10
      nf=nf+5000000
      ny=ny+5000000
      endif
      goto 50     
100	 close(unit=2)
c
c-----Map Breakpoints found in the "BPs.txt" file
c-----A similar approach was used to map the breakpoints
c-----except that the breakpoint is drawn as a tick mark to 
c-----to the left of the vertical line.
c
      write(1,'("1 1 1 1 setcmykcolor")')
      nf=0
      ny=5000000
      ix=201600
      iy=10 
      open(unit=2,file='BPs.txt')
150   read(2,*,end=200) nx1,nx2
c
c-----nx - the midpoint of the breakpoint
c
      nx=(nx1+nx2)/2
      if(nx.gt.0.and.nx.le.5000000) then
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') ix,nx+iy-nf
      write(1,'("16000 0 rlineto")') 
      write(1,'("stroke")')
      elseif (nx.gt.5000000.and.nx.le.10000000) then
      nf=5000000
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') 2*ix,nx+iy-nf
      write(1,'("16000 0 rlineto")')
      write(1,'("stroke")')
      elseif (nx.gt.10000000.and.nx.le.15000000) then
      nf=10000000
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') 3*ix,nx+iy-nf
      write(1,'("16000 0 rlineto")')
      write(1,'("stroke")')
      elseif (nx.gt.15000000.and.nx.le.20000000) then
      nf=15000000
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') 4*ix,nx+iy-nf
      write(1,'("16000 0 rlineto")') 
      write(1,'("stroke")')
      elseif (nx.gt.20000000.and.nx.le.25000000) then
      nf=20000000
      write(1,'(" newpath")')
      write(1,'(2i11," moveto")') 5*ix,nx+iy-nf
      write(1,'("16000 0 rlineto")') 
      write(1,'("stroke")')
      endif
      goto 150     
200	 close(unit=2)
      write(1,'(" showpage")')	
      stop
      end
