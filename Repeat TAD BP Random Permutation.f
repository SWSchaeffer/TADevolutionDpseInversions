      program rrp
c
c-----Repeat_Random_Permutation
c
c-----Performs a random permutation test
c
c-----Written by: Stephen W. Schaeffer
c-----The Pennsylvania State University
c-----Dept. of Biology
c-----Version: 1.0
c-----Date: 4 February 2022
c
      integer ir(3000),it(3000),ib(3000),i,j,k,l,nran,itmp
      real otab(2,2,2),irs(2),its(2),ibs(2),as,etab(2,2,2),chisq
      real rtab(2,2,2),ntab(2,2,2),xtab(2,2,2),nchi,xchi
      character cl(2)
c
c-----Set the random number seed - seconds past midnight
c-----nsec  - seed passed to the random number generator
c-----iseed - record the initial seed
c
 
      nsec=-1*int(secnds(0.0)*10)
      iseed=nsec
c
c-----cl(i) - classification as either absent (-) or present (+)
c
      cl(1)='-'
      cl(2)='+'
c
c-----Input the data from the 'Repeat_TAD_BP.txt' file 
c-----column 1 intergenic region          - [0=without repeat; 1=with repeat]
c-----column 2 TAD boundary designation   - [0=no; 1=yes]
c-----column 3 Breakpoint                 - [0=no; 1=yes]
c
c-----ir(i) - repeat of the ith intergenic region
c-----it(i) - TAD designation of the ith region
c-----ib(i) - Breakpoint of the ith region
c
c-----otab(i,j,k) - observed count of the ith, jth, kth class of region
c
      nd=0
      open(unit=1,file='Repeat_TAD_BP.txt')
50    nd=nd+1
      read(1,*,end=100) ir(nd),it(nd),ib(nd)
      goto 50
100   close(unit=1)
      nd=nd-1
      otab=0.0
      do i=1,nd
      j=ir(i)+1
      k=it(i)+1
      l=ib(i)+1
      otab(j,k,l)=otab(j,k,l)+1.0
      end do
c
c-----Estimate (-) and (+) for Repeats
c
c-----irs(i) - count of the (-) and (+) classification of the intergenic region
c
      irs=0.
      do i=1,2
      do j=1,2
      do k=1,2
      irs(i)=irs(i)+otab(i,j,k)
      end do
      end do
      end do
c
c-----Estimate (-) and (+) for TADs
c
c-----its(i) - count of the (-) and (+) TAD classification of the intergenic region
c

      its=0.
      do j=1,2
      do i=1,2
      do k=1,2
      its(j)=its(j)+otab(i,j,k)
      end do
      end do
      end do
c
c-----Estimate (-) and (+) for BPs
c
c-----ibs(i) - count of the (-) and (+) Breakpoint classification of the intergenic region
c
      ibs=0.
      do k=1,2
      do i=1,2
      do j=1,2
      ibs(k)=ibs(k)+otab(i,j,k)
      end do
      end do
      end do
c
c-----Calculate the total number
c
c-----as - total number of intergenic regions
c
      as=0.0
      do i=1,2
      as=as+irs(i)
      end do
c
c-----Estimate the Expected values
c
c-----etab(i,j,k) - expected number of the ith, jth, kth class of region
c
      do i=1,2
      do j=1,2
      do k=1,2
      etab(i,j,k)=(irs(i)*its(j)*ibs(k))/(as*as)
      end do
      end do
      end do
c
c-----Estimate Observed Chi Square
c
c-----chisq - observed Chi Square value
c
      chisq=0.0
      do i=1,2
      do j=1,2
      do k=1,2
      chisq=chisq+((otab(i,j,k)-etab(i,j,k))**2.)/etab(i,j,k)
      end do
      end do
      end do
c
c-----Randomly Permute the Data
c
c-----isim - number of permutations
c-----ntab(i,j,k) - minimum of the ith, jth, kth class of region from the permutations
c-----xtab(i,j,k) - maximum of the ith, jth, kth class of region from the permutations
c-----nchi        - minimum Chi Square value from the permutations
c-----xchi        - maximum Chi Square value from the permutations
c
	 ntab=as
	 nchi=10000.
      xtab=0.
      xchi=0.
c
c-----Open an output file 'Repeat_TAD_BP_Results.csv' to collect the random permutation data
c-----from the isim random permutations
c
      open(unit=2,file='Repeat_TAD_BP_Results.csv')
	 write(2,'(\"(---),(-+-),(--+),(-++),(+--),(++-),(+-+),(+++)")')
      write(2,'(",ChiSq")')
	 do isim=1,1000000
c
c-----Permute the repeat vector
c
      do i=nd,1,-1
	 nran=int(ran1(nsec)*float(i))+1
	 itmp=ir(nran)
	 ir(nran)=ir(i)
	 ir(i)=itmp
	 end do
c
c-----Permute the TAD vector
c
	 do i=nd,1,-1
	 nran=int(ran1(nsec)*float(i))+1
	 itmp=it(nran)
	 it(nran)=it(i)
	 it(i)=itmp
	 end do
c
c-----Permute the Breakpoint vector
c
	 do i=nd,1,-1
	 nran=int(ran1(nsec)*float(i))+1
	 itmp=ib(nran)
	 ib(nran)=ib(i)
	 ib(i)=itmp
	 end do
c
c-----Permute the repeat vector
c
	 rtab=0
c
c-----rtab(i,j,k) - count of the ith, jth, kth class of region from permutation vectors
c
      do i=1,nd
      j=ir(i)+1
      k=it(i)+1
      l=ib(i)+1
      rtab(j,k,l)=rtab(j,k,l)+1.0
      end do
c
c-----rchisq - Chi Square value from the permuted data set
c
      rchisq=0.0
      do i=1,2
      do j=1,2
      do k=1,2
      rchisq=rchisq+((rtab(i,j,k)-etab(i,j,k))**2.)/etab(i,j,k)
      end do
      end do
      end do
c
c-----Test the random permuted data for a new minimum or maximum
c
      do i=1,2
      do j=1,2
      do k=1,2
      if(rtab(i,j,k).lt.ntab(i,j,k)) ntab(i,j,k)=rtab(i,j,k)
      if(rtab(i,j,k).gt.xtab(i,j,k)) xtab(i,j,k)=rtab(i,j,k)
      end do
      end do
      end do
c
c-----Test the random permuted Chi Square values for a new minimum or maximum
c
      if(rchisq.lt.nchi) nchi=rchisq
      if(rchisq.gt.xchi) xchi=rchisq      
	 write(2,'(\f8.2,7(",",f8.2))') (((rtab(i,j,k),j=1,2),k=1,2),i=1,2)
      write(2,'(",",f8.2)') rchisq
      end do
	 close (unit=2)
c
c-----Output the results to the 'Repeat_TAB_BP_Permutation_Results.txt' file
c
      open(unit=1,file='Repeat_TAB_BP_Permutation_Results.txt')
      write(1,'("Repeat, TAD, and Breakpoint data")')
      write(1,'("Drosophila pseudoobscura Chromosome 3")')
      write(1,'("Intergenic Regions = ",i4)') nd
      write(1,'("Random seed= ",i15)') iseed
      write(1,'(1x)')
      write(1,'(1x)')
      do i=1,2
      write(1,'("Repeat(",a1,")")') cl(i)
      write(1,'(6x,2(8x,"TAD(",a1,")",4x))') (cl(j),j=1,2)
      do k=1,2
      write(1,'(\"BP(",a1)') cl(k)
      write(1,'(")   Obs",2(3x,f8.2,7x))') (otab(i,j,k),j=1,2)
      write(1,'(8x,"Exp",2x,2(1x,f8.2,9x))') (etab(i,j,k),j=1,2)
      write(1,'(1x)')
      write(1,'("Permute Min",2x,2(1x,f8.2,9x))') (ntab(i,j,k),j=1,2)
      write(1,'("Permute Max",2x,2(1x,f8.2,9x))') (xtab(i,j,k),j=1,2)
      write(1,'(1x)')
      end do
      write(1,'(1x)')
      end do
      write(1,'("Frequencies")')
      write(1,'(6x,"Repeat",7x,"TAD",8x,"BP")')
      do i=1,2
      write(1,'("(",a1,")",3f10.0)') cl(i),irs(i),its(i),ibs(i)
      end do
      write(1,'("Total",f8.0,2f10.0)') as,as,as
      write(1,'(1x)')
      write(1,'(1x)')
      write(1,'("        Observed Chi-Square = ",f10.2)') chisq
      write(1,'(1x)')      
      write(1,'("Permute Minimum  Chi-Square = ",f10.2)') nchi      
      write(1,'("Permute Maximum  Chi-Square = ",f10.2)') xchi      
      stop
      end
c
c-----Random number generator
c	
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
