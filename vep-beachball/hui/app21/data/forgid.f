      program trangid
      implicit real*8 (a-h,o-z)
c      character*12 fname,filename(20)
c      character*12 fname1,fname2
c      character*20 filename1(20)
c       common /f/ filename1
      logical filflgdisp(20)
      DIMENSION  NODVAR(3,1000000),COOR(3,1000000),
     &  U0(3,1000000),U1(3,1000000),U2(3,1000000),
     &  BFU(3,1000000),NODE(10000000)
c     $ ,emate(10000)
        character*3 uname(20)
       dimension disid(10)
      
C...... OPEN ELEM0 FILE
        OPEN (30,FILE='elem0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN COOR0 FILE
        OPEN (31,FILE='coor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ID0 FILE
        OPEN (32,FILE='id0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP0 FILE (Boundary condition file)
        OPEN (33,FILE='disp0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP1 FILE (Initial value file for displacement)
C        OPEN (34,FILE='disp1',FORM='UNFORMATTED',STATUS='UNKNOWN')
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
c 
       iprint  = 1
       kcoor=ncoor
       READ (30) NUM,NNODE,
     *           ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
      WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
cc      WRITE(*,*) 'NODE ='
cc      WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
C.......OPEN AND READ COOR file
        READ (31) NUMNOD,NCOOR,((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
c        WRITE(*,*) 'NUMNOD,NCOOR=',NUMNOD,NCOOR
C.......OPEN AND READ ID FILE
        READ (32) KNODE,KDGOF,((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C.......OPEN AND READ DISP0 FILE
        READ (33) NUMNOD,NODDOF,((BFU(I,J),I=1,NODDOF),J=1,NUMNOD)
C.......OPEN AND READ DISP1 FILE
C        READ (34) NUMNOD,NODDOF,((U0(I,J),I=1,NODDOF),J=1,NUMNOD)
C.......OPEN AND READ DISP2 FILE
      if (filflgdisp(1)) then
        READ (35) NUMNOD,NODDOF,((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
        endif
C.......OPEN AND READ DISP3 FILE
      if (filflgdisp(2)) then
        READ (36) NUMNOD,NODDOF,((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
        endif
        print *,numnod,num
        KNODE=NUMNOD
        KDGOF=NCOOR
        NODDOF=KDGOF
        NNE = NNODE
        nne = nne-1
c      if(nelemtype.eq.1) filename1(20) = 'Triangle'
c      if(nelemtype.eq.2) filename1(20) = 'Tetrahedra'
c      if(nelemtype.eq.3) filename1(20) = 'Hexahedra'
c      if(nelemtype.eq.4) filename1(20) = 'Quadrilateral'
 
C...... OPEN ELEMENT FILE FOR DATA PARTITIONING
       if(iprint.eq.1) then
c...... OPEN ELEMENT FILE FOR DATA PARTITIONING
        open(10,file='partition.flavia.msh',status='unknown',
     +        form='formatted')
        write(10,*)'Mesh "Whole" Dimension 3 Elemtype Hexahedra Nnode',8
c        write(10,*)'Mesh "Whole" Dimension 3 Elemtype Prism Nnode',6
        write(10,*)'Coordinates'
        do i=1,numnod
        write(10,1000) i,(coor(j,i),j=1,ncoor)
        end do
        write(10,*) 'End coordinates'
        write(10,*) 'Elements'
        do i=1,num
        write(10,1100) i,(NODE((I-1)*NNODE+j),j=1,9)
c        write(10,1100) i,NODE((I-1)*NNODE+1),NODE((I-1)*NNODE+2),
c     +                 NODE((I-1)*NNODE+3),
c     +                 NODE((I-1)*NNODE+4),NODE((I-1)*NNODE+5),
c     +                 NODE((I-1)*NNODE+6),
c     +                 NODE((I-1)*NNODE+7)
        end do
        write(10,*)'End elements'
        close(10)
c        GOTO 9999
c        NODVAR(1,100)=1
c        NODVAR(2,100)=1
c        NODVAR(3,100)=1
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
        write(10,*)
        close(10)
c
9999    CONTINUE        
        end if
1000    format(i10,3e15.5)
1100    format(10i10)
1200    format(i10,e15.5)
1201    format(i10,10e15.5)
        stop
        end
c
c
