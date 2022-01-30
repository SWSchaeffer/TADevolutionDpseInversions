      program bptad
c
c-----BP TAD Power
c-----Written by: Stephen W. Schaeffer
c-----The Pennsylvania State University
c-----Dept. of Biology
c-----Version: 1.0
c-----Date: 16 September 2021
c
c-----nts  - Number of total nucleotide sites
c-----isim - Number of simulations 
c
      parameter(nts=23510042,isim=1000000)
      integer seq(nts),nsec,nbeg,nend,ir,nout,nhist(0:99)
      integer nsite(0:2)
c
c-----nbp  - Number of breakpoints
c
      write(*,'("How many breakpoints on the chromosome?")')
      read(*,*) nbp
c
c-----Choose a random number seed
c-----nsec   - random seed
c
      open(unit=2,file='BP_in_TADs_Bootstrap_Summary_Data.txt')
      write(*,'("1. New seed")')
      write(*,'("2. Old seed")')
      read(*,*) nseed
      select case(nseed)
      case(1)
c
c-----Initialize the random number seed - number of seconds past midnight
c
      nsec=-1*int(secnds(0.0)*10)
c
c-----Output the seed
c      
      write(2,'("seed = ",i15)') nsec
      case(2)
c
c-----User chosen random number seed - if positive value multiply by -1
c
      write(*,'("Input a seed.")')
      read(*,*) nsec
      if(nsec.gt.0) nsec=-1*nsec
c
c------Output the seed
c      
      write(2,'("seed = ",i15)') nsec
      end select
c
c-----Initialize the sequence array to 0
c-----seq(i) - array of chromosome nucleotides
c
      seq=0
c
c-----Map TAD boundary nucleotides from boundary.tsv
c-----Change seq(i) to 1
c-----nbeg   - beginning nucleotide of gene i in the boundary.tsv file
c-----nend   - end       nucleotide of gene i in the boundary.tsv file
      open(1,file='boundary.tsv')
50    read(1,*,end=100) nbeg,nend
      do i=nbeg,nend
      seq(i)=1
      end do
      goto 50
100   close(1)
c
c-----Map gene nucleotides from Gene_Coor.tsv
c-----Change seq(i) value to 2 if not already 1 i.e., a TAD boundary
c-----If a gene overlaps with the defined TAD boundary, then the site is called as a gene
c-----nbeg   - beginning nucleotide of gene i in the boundary.tsv file
c-----nend   - end       nucleotide of gene i in the boundary.tsv file
c
      open(1,file='Gene_Coor.tsv')
150   read(1,*,end=200) nbeg,nend
      do i=nbeg,nend
      seq(i)=2
      end do
      goto 150
200   close(1)
c
c-----Summary of seq(i) values
c-----seq(i)=0 - Nucleotide in noncoding regions
c-----seq(i)=1 - Nucleotide in a TAD boundary region
c-----seq(i)=2 - Nucleotide in a gene region
c-----nsite(0) - Total number of nucleotides in noncoding regions
c-----nsite(1) - Total number of nucleotides in TAD boundary regions
c-----nsite(2) - Total number of nucleotides in gene regions
c
      nsite=0
      do i=1,nts
      nsite(seq(i))=nsite(seq(i))+1
      end do
      write(2,'("Nucleotide Site Summary Data")')
      write(2,'(\"Noncoding Sites   = ",i8)') nsite(0)
      write(2,'(4x,"% = ",f5.1)') (float(nsite(0))*100.)/float(nts)
      write(2,'(\"TAD Boundary Sites= ",i8)') nsite(1)
      write(2,'(4x,"% = ",f5.1)') (float(nsite(1))*100.)/float(nts)
      write(2,'(\"Gene Sites        = ",i8)') nsite(2)
      write(2,'(4x,"% = ",f5.1)') (float(nsite(2))*100.)/float(nts)
      write(2,'("Sum Sites         = ",i8)') nsite(0)+nsite(1)+nsite(2)
      write(2,'(1x)')
c
c-----Run the bootstrap test
c-----Boundary_test.csv has the number of breakpoints that occur in TAD boundaries for each simulation
c
      open(1,file='Boundary_test.csv')
      write(1,'("Simulation,TAD_Boundary_Hit_Number")')
c
c-----Do isim random draws of nbp breakpoints
c
      nhist=0
      do nsim=1,isim
c
c-----nout - Number of breakpoints within TAD boundaries
c
      nout=0
c
c-----Choose nbp breakpoints
c
      do i=1,nbp
c
c-----Choose a random breakpoint
c-----If the random nucleotide site is within a gene, then resample
c-----A priori knowledge is that breakpoints do not disrupt gene regions 
c-----ir        - random breakpoint nucleotide site
c-----nhist(i) - Number of simulations with i breakpoints within TAD boundaries
c
250   ir=int(ran1(nsec)*float(nts))+1
      if(seq(ir).eq.2) goto 250
      nout=nout+seq(ir)
      end do
      write(1,'(i10,",",i3)') nsim,nout
      nhist(nout)=nhist(nout)+1
      end do
      close(1)
      write(2,'("Summary of Simulation Data for ",i2," BPs")') nbp
      write(2,'(\"TAD Hits",4x,"Sim Numbers",11x,"Prob")')
      write(2,'(5x,"Cumul Prob",2x,"1-(Cumul Prob)")')
      ist=0
      afr=0.
      do i=0,nbp
      ist=ist+nhist(i)
      afr=afr+(float(nhist(i))/float(isim))
      write(2,'(\6x,i2,7x,i8)') i,nhist(i)
      write(2,'(\4x,f11.9)') float(nhist(i))/float(isim)
      write(2,'(\4x,f11.9)') afr
      write(2,'(4x,f11.9)') 1.-afr
      end do
      write(2,'("Total",i18,4x,f11.9)') ist,afr
      stop
      end
c
c-----Random number generator
c	
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     +NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
