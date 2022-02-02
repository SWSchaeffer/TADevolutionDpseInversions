      program esfcb
c
c-----Extract Sequence from Contigs Batch
c-----Written by: Stephen W. Schaeffer
c-----The Pennsylvania State University
c-----Dept. of Biology
c-----Version: 1.0
c-----Date: 4 May 2021
c
      integer ity,nctg,ip
      integer nbeg,nend
      character*1 seqi(35000000),seqo(35000000),ori
      character*4 ctg
      character*8 cbeg,cend
      character*30 fileo
      character*40 filein
c
c-----ity - type identification
c
      write(*,'("1. Contig")')
      write(*,'("2. Scaffold")')
      write(*,'("Extract a Contig or Scaffold? [1 or 2]")')
      read(*,*) ity
      write(*,'("Choose the contig number")')
      read(*,*) nctg
      write(*,'("How many bases per line in *.FAS file")')
      read(*,*) nc 
      if(nctg.lt.10) then
      write(ctg,'(i1)') nctg
      elseif(nctg.lt.100) then
      write(ctg,'(i2)') nctg
      elseif(nctg.lt.1000) then
      write(ctg,'(i3)') nctg
      elseif(nctg.lt.10000) then
      write(ctg,'(i4)') nctg
      end if
      if(ity.eq.1) then
      filein='contig_'//trim(ctg)//'.fas'
      elseif(ity.eq.2) then
      filein='scaffold_'//trim(ctg)//'.fas'
      end if
      ilen=len(trim(filein))
      open(unit=1,file=trim(filein))
      read(1,'(1x)')
      select case (nc)
      case(50)
      nb=-49
      ne=0
50    nb=nb+50
      ne=ne+50
      print *,'Input bases ',nb,' - ',ne
      read(1,'(50a1)',end=100) (seqi(j),j=nb,ne)
      goto 50
100   close(unit=1)
      ne=ne-50
      case(60)
      nb=-59
      ne=0
150   nb=nb+60
      ne=ne+60
      read(1,'(60a1)',end=200) (seqi(j),j=nb,ne)
      goto 150
200   close(unit=1)
      ne=ne-60
      case(80)
      nb=-79
      ne=0
250   nb=nb+80
      ne=ne+80
      read(1,'(80a1)',end=300) (seqi(j),j=nb,ne)
      goto 250
300   close(unit=1)
      ne=ne-80
      end select
      write(*,'(1x)')
      write(*,'("The contig or scaffold has ",i10," nucleotides.")') ne
      open(unit=2,file='extract_list.txt')
350   read(2,*,end=400) fileo,nbeg,nend,ori
      ip=scan(fileo,' ')
      nbeg=nbeg+1
      nend=nend-1
      write(cbeg,'(i8.8)') nbeg
      write(cend,'(i8.8)') nend
      if(ori.eq.'+') then
      open(unit=1,file=fileo//'_'//cbeg//'_'//cend//'.fas')
      write(1,'(">",a30,"_",a8,"_",a8)') fileo,cbeg,cend
      write(1,'(60a1)') (seqi(j),j=nbeg,nend)
      close(unit=1)
      else
      open(unit=1,file=fileo(:ip-1)//'.fas')
      write(1,'(">",a30,"_",a8,"_",a8)') fileo,cbeg,cend
      do i=nbeg,nend
      seqo(i)=seqi(nend-i+1)
      if(seqo(i).eq.'A') then
      seqo(i)='T'
      elseif(seqo(i).eq.'C') then
      seqo(i)='G'
      elseif(seqo(i).eq.'G') then
      seqo(i)='C'
      elseif(seqo(i).eq.'T') then
      seqo(i)='A'
      end if
      end do
      write(1,'(60a1)') (seqo(j),j=nbeg,nend)
      close(unit=1)
      end if
      goto 350
400   close(unit=2)
      stop
      end
