c      subroutine Mazpartition
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      character*12 fname1,fname2
      character*20 filename1(20)
      include 'partdata.h'
      include 'memalloc.h'
      logical filflgdisp(20)
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
      nsource=0
c     Get the element data file and domain partition number'
      filename(19) = 'elem0.elm'
c
c     Check if the file name and the partition number are correct !
c
c      write(*,*)'The input partition filname and number of subdomains:'
c      write(*,*) filename(1), nparts
c
c     get the lmddm control file
c
      open(1,file='partition.dat',form='formatted',status='unknown')
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) numblk,numtyp,nmdof,numtypl,kdgofl,keymt,lgio
      read(1,*)
      read(1,*) t0,tmax,dt
      close(1)
      nparts = numblk
      maxlm=0
      idisp1 = 0
      idisp2 = 0
      inquire(file='disp2',exist=filflgdisp(1))
      if (filflgdisp(1)) idisp1 = 1
      inquire(file='disp3',exist=filflgdisp(2))
      if (filflgdisp(2)) idisp2 = 1
c
c
      open(2,file='part',form='formatted',status='unknown')
      write(2,*) numblk
      close(2)
c
      nsource=0
c
c     open data file from GID preprocessing
c
C...... OPEN ELEM0 FILE
        OPEN (30,FILE='elem0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN COOR0 FILE
        OPEN (31,FILE='coor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ID0 FILE
        OPEN (32,FILE='id0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP0 FILE (Boundary condition file)
        OPEN (33,FILE='disp0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP1 FILE (Initial value file for displacement)
        OPEN (34,FILE='disp1',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C     open data file for other initial value problem
C
C...... OPNE DISP2 FILE (Initial value file for velocity)
      inquire(file='disp2',exist=filflgdisp(1))
      if (filflgdisp(1)) then
      open(35,file='disp2',form='unformatted',status='unknown')
      endif
C...... OPNE DISP3 FILE (Initial value file for acceleration)
      inquire(file='disp3',exist=filflgdisp(2))
      if (filflgdisp(2)) then
      open(36,file='disp3',form='unformatted',status='unknown')
      endif
C
c
c     write partition file for domain decomposition algorithms
c
C...... OPEN ELEM0 FILE
      open(21,file='melem0',form='unformatted',status='unknown')
C...... OPEN COOR0 FILE
      open(22,file='mcoor0',form='unformatted',status='unknown')
C...... OPEN ID0 FILE
      open(23,file='mid0',form='unformatted',status='unknown')
C...... OPEN DISP0 FILE (Boundary condition file)
      open(24,file='mdisp0',form='unformatted',status='unknown')
C...... OPEN DISP1 FILE (Initial value file for displacement)
      open(25,file='mdisp1',form='unformatted',status='unknown')
C
      if (filflgdisp(1)) then
      open(26,file='mdisp2',form='unformatted',status='unknown')
      endif
      if (filflgdisp(2)) then
      open(27,file='mdisp3',form='unformatted',status='unknown')
      endif
C...... OPEN ELEMENT TYPE FILE FOR EACH SUBDOMAIN (etype0)
      open(29,file='metype0',form='unformatted',status='unknown')
c...... open subnodes file
      open(40,file='subnodes0',form='unformatted',status='unknown')
c...... open lmddm file
      open(41,file='mlmddm',form='unformatted',status='unknown')
c
      READ(33) KNODE,KDGOF
      REWIND(33)
      READ(30) NUM0,NNODE0
      REWIND(30)
      READ (31) NUMNOD,NCOOR
      REWIND(31)
      num = num0
      nnode = nnode0
      KVAR=KNODE*KDGOF
      kelem=num*nnode
      knode3 = knode*3
      kelem0 = kelem
      kelemsub = kelem/numblk*2
c
      kelemi=kelem/numblk+15000
      knodei=knode/numblk*3
c
        DO 1000 iblk=1,NUMBLK
c
c
c     write lmddm file for each subdomain
c
      write(41) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
     & idisp1,idisp2
      write(41) lgio,t0,tmax,dt
c
      knb1=kdgof*knode*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      kna4=ncoor*knode*1
      kna1=kdgof*knode*1
      kna2=kdgof*knode*1
      kna3=kdgof*knode*1
      kna5=kdgof*knode*1
      knb7=kelem0*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb2=knode*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb6=knode*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb3=num0*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb5=numnod*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      knb4=numnod*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb9=kelemi*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      kna8=ncoor*knodei*1
      knb10=kdgof*knodei*1
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      kna9=kdgof*knodei*1
      knb12=numtyp*1
      if (knb12/2*2 .lt. knb12) knb12=knb12+1
      knb13=numtyp*1
      if (knb13/2*2 .lt. knb13) knb13=knb13+1
      knb14=numtyp*1
      if (knb14/2*2 .lt. knb14) knb14=knb14+1
      knb15=numtyp*1
      if (knb15/2*2 .lt. knb15) knb15=knb15+1
      knb11=numtyp*1
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      knb8=kelemsub*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      kna7=9150000*1
      kna6=360000*1
      knb16=kdgof*knode*1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
      kna8=kna8+kna7
      kna9=kna9+kna8
      cknb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      knb5=knb5+knb4
      knb6=knb6+knb5
      knb7=knb7+knb6
      knb8=knb8+knb7
      knb9=knb9+knb8
      knb10=knb10+knb9
      knb11=knb11+knb10
      knb12=knb12+knb11
      knb13=knb13+knb12
      knb14=knb14+knb13
      knb15=knb15+knb14
      knb16=knb16+knb15
      call azpartition(knode,kdgof,kcoor,time,dt,iblk,
     *nsource,num,nnode,numnod,ncoor,nparts,numblk,kelem0,
     *kelemi,knodei,filflgdisp,numtyp,lgio,ndual,knode3,filename1,
     *num0,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),aa(kna7),aa(kna8),
     *ia(knb0),ia(knb1),ia(knb2),ia(knb3),ia(knb4),
     *ia(knb5),ia(knb6),ia(knb7),ia(knb8),ia(knb9),
     *ia(knb10),ia(knb11),ia(knb12),ia(knb13),ia(knb14),
     *filename,kelemsub,ia(knb15))
c
1000    CONTINUE
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
      close(36)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(40)
      close(41)
 
      end
      subroutine azpartition(knode,kdgof,kcoor,time,dt,iblk,
     *nsource,num,nnode,numnod,ncoor,nparts,numblk,kelem0,
     *kelemi,knodei,filflgdisp,numtyp,lgio,ndual,knode3,filename1,
     *num0,u0,u1,u2,coor,bfu,emate,emateall,
     *coori,ui,nodvar,inode,idelem,idnode1,idnode,jnode,
     *node,nodeall,nodei,nodvari,idet,numa,nnea,mmta,
     *nmta,
     *filename,kelemsub,nodvar1)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION  NODVAR(KDGOF,KNODE),COOR(NCOOR,KNODE),
     &  U0(KDGOF,KNODE),U1(KDGOF,KNODE),U2(KDGOF,KNODE),
     &  BFU(KDGOF,KNODE),
     &  NODE(kelem0),inode(knode),jnode(knode),
     &  IDELEM(NUM0),IDNODE(NUMNOD),IDNODE1(numnod),
     &  nodei(kelemi),coori(ncoor,knodei),nodvari(kdgof,knodei),
     &  ui(kdgof,knodei),numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),
     &  nodeall(kelemsub),emate(9150000),emateall(360000),
     &  nodvar1(KDGOF,KNODE)
       character*12 fname1,fname2
       character*20 filename1(20)
       character*3 uname(20)
       logical filflgdisp(20)
       include 'partdata.h'
       dimension disid(10)
c      parameter (fname1='elem0.elmc
       iprint  = 0
       kcoor=ncoor
       if(iblk.eq.1) then
c
       if(nparts.ne.numblk) then
       write(*,*)'Error when get partition information of the dimain'
       stop 1111
       end if
       READ (30) NUM,NNODE,
     *           ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
cc      WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
cc      WRITE(*,*) 'NODE ='
cc      WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
C.......OPEN AND READ COOR file
        READ (31) NUMNOD,NCOOR,((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
c        WRITE(*,*) 'NUMNOD,NCOOR=',NUMNOD,NCOOR
C.......OPEN AND READ ID FILE
        READ (32) NSKP,NSKP,((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C.......OPEN AND READ DISP0 FILE
        READ (33) NUMNOD,NODDOF,((BFU(I,J),I=1,NODDOF),J=1,NUMNOD)
C.......OPEN AND READ DISP1 FILE
        READ (34) NUMNOD,NODDOF,((U0(I,J),I=1,NODDOF),J=1,NUMNOD)
C.......OPEN AND READ DISP2 FILE
      if (filflgdisp(1)) then
        READ (35) NUMNOD,NODDOF,((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
        endif
C.......OPEN AND READ DISP3 FILE
      if (filflgdisp(2)) then
        READ (36) NUMNOD,NODDOF,((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
        endif
C.......Start Data partition.....................................
      NNE = NNODE
      nne = nne-1
      nelemtype=0
      if(nne.eq.3) nelemtype=1
      if((nne.eq.4).and.(ncoor.eq.2)) nelemtype=4
      if((nne.eq.4).and.(ncoor.eq.3)) nelemtype=2
      if((nne.eq.8).and.(ncoor.eq.3)) nelemtype=3
      if((nne.eq.10).and.(ncoor.eq.3)) nelemtype=2
      if(nelemtype.eq.1) filename1(20) = 'Triangle'
      if(nelemtype.eq.2) filename1(20) = 'Tetrahedra'
      if(nelemtype.eq.3) filename1(20) = 'Hexahedra'
      if(nelemtype.eq.4) filename1(20) = 'Quadrilateral'
 
C...... OPEN ELEMENT FILE FOR DATA PARTITIONING
        OPEN (10,FILE='elem0.elm',FORM='FORMATTED',STATUS='UNKNOWN')
      if((nne.eq.10).and.(ncoor.eq.3)) then
        WRITE(10,*) NUM,nelemtype
      do i=1,num
      write(10,8001) (NODE((I-1)*NNODE+J),J=1,4)
      end do
      else
      WRITE(10,*) NUM,nelemtype
      do i=1,num
      write(10,8001) (NODE((I-1)*NNODE+J),J=1,nne)
      end do
      end if
      close(10)
8001    format(8i10)
c
c       get the partition data for each subdomain
c
        call partdmain(3,fname1,fname2,nparts,idelem,idnode,ne,np)
c
c       Chech the partition data
c
        if(ne.ne.num) then
      write(*,*) 'fatal error when partition the whole domain'
      write(*,*) 'Different element number!'
      stop 1111
      end if
        if(np.ne.numnod) then
      write(*,*) 'fatal error when partition the whole domain'
      write(*,*) 'Different node number!'
      stop 1111
      end if
 
c.......Check the partition using GID postprocessing function
c
      do i=1,numnod
      idnode1(i) = idnode(i)
      end do
c
       if(iprint.eq.1) then

c...... OPEN ELEMENT FILE FOR DATA PARTITIONING
        open(10,file='partition.flavia.msh',status='unknown',
     +        form='formatted')
        if((nne.eq.10).and.(ncoor.eq.3)) then
        write(10,*) ' Mesh "Whole" Dimension',ncoor,
     + ' Elemtype ',filename1(20),' Nnode',4
        else
        write(10,*) ' Mesh "Whole" Dimension',ncoor,
     +' Elemtype ',filename1(20),' Nnode',nnode-1
        end if
        write(10,*) 'Coordinates'
        do i=1,numnod
        write(10,1000) i,(coor(j,i),j=1,ncoor)
        end do
        write(10,*) 'End coordinates'
        write(10,*) 'Elements'
      if((nne.eq.10).and.(ncoor.eq.3)) then
        do i=1,num
        write(10,1100) i,(NODE((I-1)*NNODE+J),J=1,4),idelem(i)
        end do
      else
        do i=1,num
        write(10,1100) i,(NODE((I-1)*NNODE+J),J=1,nnode-1),idelem(i)
        end do
      end if
        write(10,*)'End elements'
        close(10)
 
        open(10,file='partition.flavia.res',status='unknown',
     +        form='formatted')
        uname(1) = 'id1'
        uname(2) = 'id2'
        uname(3) = 'id3'
        uname(4) = 'id4'
        uname(5) = 'id5'
        uname(6) = 'id6'
        uname(7) = 'id7'
        uname(8) = 'id8'
        uname(9) = 'id9'
        uname(10) = 'id10'
      write(10,*)'GID Post Results File 1.0'
        do j = 1,noddof
        write(10,*)
        write(10,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
        write(10,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(10,*)'Values'
        do i=1,numnod
      disid(j) = nodvar(j,i)
        write(10,1201) i,disid(j)
        end do
        write(10,*)'end Values'
        end do
c
      uname(1) = 'bd1'
        uname(2) = 'bd2'
        uname(3) = 'bd3'
        uname(4) = 'bd4'
        uname(5) = 'bd5'
        uname(6) = 'bd6'
        uname(7) = 'bd7'
        uname(8) = 'bd8'
        uname(9) = 'bd9'
        uname(10) = 'bd0'
        do j = 1,noddof
        write(10,*)
        write(10,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
        write(10,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(10,*)'Values'
        do i=1,numnod
        write(10,1201) i,bfu(j,i)
        end do
        write(10,*)'end Values'
        end do
c
      uname(1) = 'u01'
        uname(2) = 'u02'
        uname(3) = 'u03'
        uname(4) = 'u04'
        uname(5) = 'u05'
        uname(6) = 'u06'
        uname(7) = 'u07'
        uname(8) = 'u08'
        uname(9) = 'u09'
        uname(10) = 'u10'
        do j = 1,noddof
        write(10,*)
        write(10,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
        write(10,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(10,*)'Values'
        do i=1,numnod
        write(10,1201) i,u0(j,i)
        end do
        write(10,*)'end Values'
        end do
 
      if (filflgdisp(1)) then
        uname(1) = 'u11'
        uname(2) = 'u12'
        uname(3) = 'u13'
        uname(4) = 'u14'
        uname(5) = 'u15'
        uname(6) = 'u16'
        uname(7) = 'u17'
        uname(8) = 'u18'
        uname(9) = 'u19'
        uname(10) = 'u10'
        do j = 1,noddof
        write(10,*)
        write(10,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
        write(10,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(10,*)'Values'
        do i=1,numnod
        write(10,1201) i,u1(j,i)
        end do
        write(10,*)'end Values'
        end do
      end if
        if (filflgdisp(2)) then
        uname(1) = 'u21'
        uname(2) = 'u22'
        uname(3) = 'u23'
        uname(4) = 'u24'
        uname(5) = 'u25'
        uname(6) = 'u26'
        uname(7) = 'u27'
        uname(8) = 'u28'
        uname(9) = 'u29'
        uname(10) = 'u20'
        do j = 1,noddof
        write(10,*)
        write(10,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
        write(10,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(10,*)'Values'
        do i=1,numnod
        write(10,1201) i,u2(j,i)
        end do
        write(10,*)'end Values'
        end do
      end if
c
        write(10,*)
      write(10,*)'Result "Dualpoint" "Load Analysis" 1 Scalar OnNodes'
        write(10,*)'ComponentNames "Dual"'
        write(10,*)'Values'
        do i=1,numnod
        dip = idnode(i)
        write(10,1200) i,dip
        end do
        write(10,*)'end Values'
      close(10)
c
        end if
c
c
c
        end if
c
       npartition = iblk
       print *,'Begin ',iblk,'  subdomain data partition...'
c      nchar=len(filename(1))
c
       if(iprint.eq.1) then
c
       n1000=(npartition/1000)
       n100=((npartition-1000*n1000)/100)
       n10=((npartition-1000*n1000-100*n100)/10)
       n1=npartition-1000*n1000-100*n100-10*n10
       if((npartition.ge.1000).and.(npartition.le.9999)) nchar=4
       if((npartition.ge.100).and.(npartition.le.999)) nchar=3
       if((npartition.ge.10).and.(npartition.le.99)) nchar=2
       if((npartition.ge.1).and.(npartition.le.9)) nchar=1
       if(nchar.gt.1) goto 2001
       open(10,file='part' //char(48+n1)// '.flavia.msh',
     +        status='unknown',
     +        form='formatted')
       open(11,file='part' //char(48+n1)// '.flavia.res',
     +        status='unknown',
     +        form='formatted')
       goto 2005
2001   continue
       if(nchar.gt.2) goto 2002
       open(10,file='part' //char(48+n10)//char(48+n1)//
     +  '.flavia.msh',
     +  status='unknown',form='formatted')
       open(11,file='part' //char(48+n10)//char(48+n1)//
     +  '.flavia.res',
     +  status='unknown',form='formatted')
       goto 2005
2002   continue
       if(nchar.gt.3) goto 2003
       open(10,file='part' //char(48+n100)//
     + char(48+n10)//char(48+n1)// '.flavia.msh',
     + status='unknown',form='formatted')
       open(11,file='part' //char(48+n100)//
     + char(48+n10)//char(48+n1)// '.flavia.res',
     + status='unknown',form='formatted')
       goto 2005
2003   continue
       if(nchar.gt.4) goto 2004
       open(10,file='part' //char(48+n1000)//char(48+n100)//
     + char(48+n10)//char(48+n1)// '.flavia.msh',
     + status='unknown',form='formatted')
       open(11,file='part' //char(48+n1000)//char(48+n100)//
     + char(48+n10)//char(48+n1)// '.flavia.res',
     + status='unknown',form='formatted')
       goto 2005
2004   continue
       write(*,*) 'fatal error, partition number more that 9999'
       stop 1111
2005   continue
c
c
c
       end if
c
c      get the node number of the iblk subdomain
c
       nodesiblk=0
       do i = 1,knode
       idnode(i) = 0
       end do
       do i = 1,knode
       if(idnode1(i).eq.(iblk-1)) then
       nodesiblk = nodesiblk + 1
       idnode(nodesiblk) = i
       end if
       end do
c
c
       do i = 1,knode
	 do j = 1,kdgof
	 nodvar1(j,i) = nodvar(j,i)
	 end do
	 end do
c
c
c       neqall=0
c       do  j = 1,knode
c         do  i = 1,kdgof
c         if (nodvar1(i,j).ge.1) then
c         neqall =  neqall + 1
c         nodvar1(i,j) = neqall
c         end if
c       end do
c       end do
        NEQ=0
        DO 20 J=1,KNODE
          DO 18 I=1,KDGOF
            IF (NODVAR1(I,J).LT.1) GOTO 18
            NEQ = NEQ + 1
            NODVAR1(I,J) = NEQ
18        CONTINUE
20      CONTINUE

c
c       write(*,*) 'nodvar,neqall =====',neqall
c       do j = 1,knode
c       write(*,6) j,(nodvar1(i,j),i=1,kdgof)
c       end do
        do  j = 1,knode
          do 28 i=1,kdgof
            if (nodvar1(i,j).ge.-1) goto 28
            n = -nodvar1(i,j)-1
27          continue
            if (nodvar1(i,n).lt.-1) then
              n=-nodvar1(i,n)-1
              goto 27
            end if
            nodvar1(i,j) = nodvar1(i,n)
28        continue
        end do
c
c        DO 30 J=1,KNODE
c          DO 28 I=1,KDGOF
c            IF (NODVAR1(I,J).GE.-1) GOTO 28
c            N = -NODVAR1(I,J)-1
c27          CONTINUE
c            IF (NODVAR1(I,N).LT.-1) THEN
c              N=-NODVAR1(I,N)-1
c              GOTO 27
c            ENDIF
c            NODVAR1(I,J) = NODVAR1(I,N)
c28        CONTINUE
c30      CONTINUE
c
c
c       write(*,*) 'nodvar after constrain conditions=====',neqall
c       do j = 1,knode
c       write(*,6) j,(nodvar1(i,j),i=1,kdgof)
c       end do
c
c     
        do i = 1,knode
        if(idnode1(i).ne.(iblk-1)) then
        do j = 1,kdgof
        if(nodvar1(j,i).gt.0) then 
        nodvar1(j,i) = -100-nodvar1(j,i)
        else if(nodvar1(j,i).eq.0) then
        nodvar1(j,i) = 0
        else
        nodvar1(j,i) = -1
        end if
        end do
        end if
        end do
c
c      write out which nodes belongs to iblk subdomain
c
       write(40) knode,kdgof
       write(40) nodesiblk,(idnode(i),i=1,1)
c
       do i=1,knode
       idnode(i)=-100
       inode(i)=0
       jnode(i)=0
       end do
       num_iblk=0
       do i=1,num
       noyes = 0
       do j = 1,nnode-1
       nodi = NODE((I-1)*NNODE+J)
       if(idnode1(nodi).eq.(iblk-1)) noyes = 1
       end do
       if(noyes.eq.1) then
       num_iblk=num_iblk+1
       do j=1,nnode
       NODEI((num_iblk-1)*NNODE+J)=NODE((I-1)*NNODE+J)
       end do
       do j=1,nnode-1
       idnode(NODE((I-1)*NNODE+J))=iblk-1
       end do
       end if
       end do
       nod_iblk=0
       do i=1,knode
       if(idnode(i).eq.(iblk-1)) then
       nod_iblk=nod_iblk+1
       inode(nod_iblk)=i
       end if
       end do
 
c       do i=1,nod_iblk
c       jnode(inode(i)) =i
c       end do
c
c    I modified here for NEU project
c
c       get the local nodvar dimension
c
        nodsub = 0
        do i = 1,knode
        if(idnode1(i).eq.(iblk-1)) nodsub = nodsub + 1
        end do
        if(nodsub.ne.nodesiblk) print *,'partition error in azpartition'
c
c       reorder of inode dimension to ascent ordering
c
        do i = 1,nod_iblk-1
        do j = i,nod_iblk
        if(inode(j).lt.inode(i)) then
        ntt  = inode(i)
        inode(i) = inode(j)
        inode(j) = ntt
        end if
        end do
        end do
c
       do i=1,nod_iblk
       jnode(inode(i)) =i
       end do
c
        N_update = 0
        do i = 1,nod_iblk
        nodi = inode(i)
        N_update = N_update + 1
        do j = 1,kdgof
        nodvar1(j,N_update) = nodvar1(j,nodi)
        end do
        end do
c
c      write out index file to subnodes file
c
       write(40) nod_iblk,(inode(i),i = 1,1)
       write(40) knode,(jnode(i),i=1,1)
       write(40) knode,(idnode1(i),i=1,1)
       write(40) knode,kdgof,((nodvar1(j,i),j=1,kdgof),i=1,nod_iblk)
       if(node_iblk.gt.knodei) print *,'Warning, knodei < nod_iblk'
c
c      deal with the coor id and disp0 disp1 disp2 disp3 file for each subdomain
c
c
c     get coor0 file for each subdomain
c
      do i=1,nod_iblk
      do j=1,ncoor
      coori(j,i)=coor(j,inode(i))
      end do
      end do
c
      write(22) nod_iblk,ncoor
      write(22) ((COORI(I,J),I=1,KCOOR),J=1,nod_iblk)
c
c     write out the partition data for each subdomain
c
c
c
       if(iprint.eq.1) then
c
c
       if((nne.eq.10).and.(ncoor.eq.3)) then
       write(10,*) 'Mesh "Part" Dimension',ncoor,' Elemtype ',
     + filename1(20),' Nnode',4
       else
       write(10,*) 'Mesh "Part" Dimension',ncoor,' Elemtype ',
     + filename1(20),' Nnode',nnode-1
       end if
       write(10,*) 'Coordinates'
       do i=1,nod_iblk
       write(10,1000) i,(coori(j,i),j=1,ncoor)
       end do
       write(10,*) 'End coordinates'
       write(10,*) 'Elements'
       do i=1,num_iblk
       if((nne.eq.10).and.(ncoor.eq.3)) then
       write(10,1100) i,(jnode(NODEI((I-1)*NNODE+J)),J=1,4),
     + NODEI((I-1)*NNODE+11)
       else
       write(10,1100) i,(jnode(NODEI((I-1)*NNODE+J)),J=1,nnode-1),
     + NODEI(I*NNODE)
       end if
       end do
       write(10,*)'End elements'
       close(10)
c
       end if
c
c
c     get id0 file for each subdomain
c
      do i=1,nod_iblk
      do j=1,kdgof
      nodvari(j,i)=nodvar(j,inode(i))
      end do
      end do
c
      write(23) nod_iblk,kdgof
      write(23) ((nodvari(I,J),I=1,kdgof),J=1,nod_iblk)

      if(iprint.eq.1) then
 
        uname(1) = 'id1'
        uname(2) = 'id2'
        uname(3) = 'id3'
        uname(4) = 'id4'
        uname(5) = 'id5'
        uname(6) = 'id6'
        uname(7) = 'id7'
        uname(8) = 'id8'
        uname(9) = 'id9'
        uname(10) = 'id10'
      write(11,*)'GID Post Results File 1.0'
        do j = 1,kdgof
        write(11,*)
        write(11,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
c        write(11,*)'Result "ID" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(11,*)'Values'
        do i=1,nod_iblk
      disid(j) = nodvari(j,i)
        write(11,1201) i,disid(j)
        end do
        write(11,*)'end Values'
        end do
c
c
        end if
c
c
c
c
c     get disp0 file for each subdomain
c
      do i=1,nod_iblk
      do j=1,kdgof
      ui(j,i)=bfu(j,inode(i))
      end do
      end do
c
      write(24) nod_iblk,kdgof
      write(24) ((ui(I,J),I=1,kdgof),J=1,nod_iblk)
c
       if(iprint.eq.1) then
      uname(1) = 'bd1'
        uname(2) = 'bd2'
        uname(3) = 'bd3'
        uname(4) = 'bd4'
        uname(5) = 'bd5'
        uname(6) = 'bd6'
        uname(7) = 'bd7'
        uname(8) = 'bd8'
        uname(9) = 'bd9'
        uname(10) = 'bd0'
        do j =1,kdgof
        write(11,*)
        write(11,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
c        write(11,*)'Result "BD_disp" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(11,*)'Values'
        do i=1,nod_iblk
        write(11,1201) i,ui(j,i)
        end do
        write(11,*)'end Values'
        end do

        end if
c
c     get disp1 file for each subdomain
c
      do i=1,nod_iblk
      do j=1,kdgof
      ui(j,i)=u0(j,inode(i))
      end do
      end do
c
      write(25) nod_iblk,kdgof
      write(25) ((ui(I,J),I=1,kdgof),J=1,nod_iblk)
c
       if(iprint.eq.1) then
      uname(1) = 'u01'
        uname(2) = 'u02'
        uname(3) = 'u03'
        uname(4) = 'u04'
        uname(5) = 'u05'
        uname(6) = 'u06'
        uname(7) = 'u07'
        uname(8) = 'u08'
        uname(9) = 'u09'
        uname(10) = 'u00'
        do j = 1,kdgof
        write(11,*)
        write(11,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
c        write(11,*)'Result "Init_u0" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(11,*)'Values'
        do i=1,nod_iblk
        write(11,1201) i,ui(j,i)
        end do
        write(11,*)'end Values'
        end do
        end if
c
c     get disp2 file for each subdomain
c
      if (filflgdisp(1)) then
      do i=1,nod_iblk
      do j=1,kdgof
      ui(j,i)=u1(j,inode(i))
      end do
      end do
c
      write(26) nod_iblk,kdgof
      write(26) ((ui(I,J),I=1,kdgof),J=1,nod_iblk)
      end if
c
      if(iprint.eq.1) then
      if (filflgdisp(1)) then
        uname(1) = 'u11'
        uname(2) = 'u12'
        uname(3) = 'u13'
        uname(4) = 'u14'
        uname(5) = 'u15'
        uname(6) = 'u16'
        uname(7) = 'u17'
        uname(8) = 'u18'
        uname(9) = 'u19'
        uname(10) = 'u10'
        do j = 1,kdgof
        write(11,*)
        write(11,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
c        write(11,*)'Result "Init_u1" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames  ', '"'//uname(j)//'"'
        write(11,*)'Values'
        do i=1,nod_iblk
        write(11,1201) i,ui(j,i)
        end do
        write(11,*)'end Values'
        end do
      end if
      end if
c
c
c     get disp3 file for each subdomain
c
      if (filflgdisp(2)) then
      do i=1,nod_iblk
      do j=1,kdgof
      ui(j,i)=u2(j,inode(i))
      end do
      end do
      write(27) nod_iblk,kdgof
      write(27) ((ui(I,J),I=1,kdgof),J=1,nod_iblk)
      end if
c
       if(iprint.eq.1) then
        if (filflgdisp(2)) then
        uname(1) = 'u21'
        uname(2) = 'u22'
        uname(3) = 'u23'
        uname(4) = 'u24'
        uname(5) = 'u25'
        uname(6) = 'u26'
        uname(7) = 'u27'
        uname(8) = 'u28'
        uname(9) = 'u29'
        uname(10) = 'u20'
        do j = 1,kdgof
        write(11,*)
        write(11,*)'Result "'//uname(j)//'" "Load Analysis" 1'//
     &' Scalar OnNodes'
c        write(11,*)'Result "Init_u2" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames  ', '"'//uname(i)//'"'
        write(11,*)'Values'
        do i=1,nod_iblk
        write(11,1201) i,ui(j,i)
        end do
        write(11,*)'end Values'
        end do
      end if
      end if
c
        do i=1,nod_iblk
        if(idnode1(inode(i)).eq.(iblk-1)) then
        ui(1,i)=iblk
        else
        ui(1,i) = 0.d0
        end if
        end do
c
       if(iprint.eq.1) then
        write(11,*)
      write(11,*)'Result "Dualpoint" "Load Analysis" 1 Scalar OnNodes'
        write(11,*)'ComponentNames "Dual"'
        write(11,*)'Values'
        do i=1,nod_iblk
        write(11,1200) i,ui(1,i)
        end do
        write(11,*)'end Values'
      end if
      close(11)
c
c     deal with element and mate information for each subdomain and multiplier
c
      rewind(30)
      nodall=0
      matall=0
      do itype=1,numtyp
      numa(itype)=0
      nnea(itype)=0
      mmta(itype)=0
      nmta(itype)=0
      READ (30) NUM,NNODE,
     *           ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
cc      WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
cc      WRITE(*,*) 'NODE ='
cc      WRITE(*,*) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
      READ (30) MMATE,NMATE,((EMATE((I-1)*NMATE+J),J=1,NMATE),
     *I=1,MMATE)
c      WRITE(*,*) 'MMATE =',MMATE,' NMATE =',NMATE
c      WRITE (*,*) 'EMATE ='
c      WRITE (*,*) ((EMATE((I-1)*NMATE+J),J=1,NMATE),
c     *I=1,MMATE)
      nne=nnode-1
      nblk=0
      do i=1,num
      noyes=1
      do j=1,nne
      if(idnode(NODE((I-1)*NNODE+J)).ne.(iblk-1)) noyes=0
      end do
      if(noyes.eq.1) then
      nblk=nblk+1
      do j=1,nnode-1
      nodei((nblk-1)*nnode+j)=NODE((I-1)*NNODE+J)
      end do
      nodei(nblk*nnode)=nblk
      imate=node(i*nnode)
      mmateall=nblk
      do j=1,nmate
      emateall(matall+(mmateall-1)*nmate+j)=emate((imate-1)*nmate+j)
      end do
      end if
      end do
c
      do i=1,nblk
      do j=1,nnode-1
      NODEI((i-1)*NNODE+J)=jnode(NODEI((i-1)*NNODE+J))
      end do
      end do
c
      do i=1,nblk
      do j=1,nnode
      NODEALL(nodall+(i-1)*NNODE+J)=NODEI((i-1)*NNODE+J)
      end do
      end do
      nodall=nodall+nblk*nnode
      numa(itype)=nblk
      nnea(itype)=nnode
c
c      do i=1,MMATE
c      do j=1,NMATE
c      EMATEall(matall+(I-1)*NMATE+J)=EMATE((I-1)*NMATE+J)
c      end do
c      end do
      matall = matall + mmateall*nmate
      print *,'In partition::: nblk,mmateall = ',nblk,mmateall
      mmta(itype) = mmateall
      nmta(itype) = nmate
      end do
C     end do of itype
      iidet = 1
C...  Write etype0 file
      write(29) numtyp,nodall,matall,iblk
      write(29) (iidet,numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
     &(inode(i),i=1,nod_iblk)
C ... Write elem0 file
      write(21) (nodeall(i),i=1,nodall)
      write(21) (emateall(i),i=1,matall)
c
c     read the initial data for next subdomain
c
c      rewind(30)
c      READ (30) NUM,NNODE,
c     *           ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
cc      WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
cc      WRITE(*,*) 'NODE ='
cc      WRITE(*,*) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
c      READ (30) MMATE,NMATE,((EMATE((I-1)*NMATE+J),J=1,NMATE),
c     *I=1,MMATE)
c      WRITE(*,*) 'MMATE =',MMATE,' NMATE =',NMATE
c      WRITE (*,*) 'EMATE ='
c      WRITE (*,*) ((EMATE((I-1)*NMATE+J),J=1,NMATE),
c     *I=1,MMATE)
c
1000  format(i10,3e15.5)
c1100  format(6i10)
1100  format(10i10)
1200  format(i10,e15.5)
1201   format(i10,10e15.5)
      return
      end
 
 
