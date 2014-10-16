      subroutine Mmbft(iblk,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)

      if (iblk.ne.0) return

        nsource=0
c ....  open lmddm file
c       open (21,file='lmddm',form='unformatted',status='old')
      ISTATUS = 1
      call openf(21,1,ISTATUS)
c       read (21) numblk,numtyp,nskip,nskip,numtypl,nskp,nskp
      call initrwf(21,iui,iur)
      iui = iui+1
      numblk=ipool(iui)
      iui = iui+1
      numtyp=ipool(iui)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      numtypl=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      num0=ipool(iui)
      iui = iui+1
      knode0=ipool(iui)
      call endrwf(21,iui,iur)
c       close(21)
cc        write(*,*) 'numblk,numtyp,nskip,numtypl =='
cc        write(*,*)  numblk,numtyp,nskip,numtypl


c        OPEN (34,file='metype0',form='unformatted',status='unknown')


C...... OPNE DISP0 FILE
c       OPEN (21,FILE='sdisp0',FORM='UNFORMATTED',STATUS='OLD')
C...... OPEN NV FILE
c       OPEN (23,FILE='nv',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN BFD FILE
c       OPEN (24,FILE='bfd',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN COOR0.BAK file
c       OPEN (25,FILE='scoor0',FORM='UNFORMATTED',STATUS='OLD')
C...... OPEN UV.DDA file
c       OPEN (26,FILE='uv.dda',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN UNOD FILE
c       OPEN (27,FILE='unod',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN COOR0 FILE
c       OPEN (28,FILE='coor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
c        DO 1000 iblk=1,NUMBLK
      ISTATUS = 4
      call openf(21,4,ISTATUS)
c       READ(21) KNODE,KDGOF
      call initrwf(21,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      call endrwf(21,iui,iur)
        KVAR=KNODE*KDGOF
c        WRITE(*,*) 'KNODE,KDGOF,KVAR ='
c        WRITE(*,*) KNODE,KDGOF,KVAR
      ISTATUS = 2
      call openf(25,2,ISTATUS)
c        READ(25) NSKIP,KCOOR
      call initrwf(25,iui,iur)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      kcoor=ipool(iui)
      call endrwf(25,iui,iur)

      knb1=numtyp*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=numtyp*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numtyp*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb4=numtyp*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb5=numtyp*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      knb6=kdgof*knode*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb7=knode0
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb8=num0
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      knb9=num0*8
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      knb10=num0*8
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      knb11=knode0
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      kna7=kcoor*knode*1
      kna1=kdgof*knode*1
      kna2=kdgof*knode*1
      kna5=kdgof*knode*1
      kna3=kdgof*knode*1
      kna4=kdgof*knode*1
      kna6=kdgof*knode*1
      kna8=kdgof*knode*1
      kna9=knode0*6
      kna10=knode0*6
      kna11=knode0*6
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
      kna10=kna10+kna9
      kna11=kna11+kna10
      knb0=1
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
      call bftmas(knode,kdgof,kcoor,nsource,
     *iblk,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),aa(kna7),ia(knb0),
     *ia(knb1),ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *filename,numtyp,knode0,ia(knb6),ia(knb7),istop,kend,
     *num0,ia(knb8),ia(knb9),numblk,aa(kna8),aa(kna9),aa(kna10),
     *ia(knb10))
 
1000    CONTINUE

c        close(34)
c       CLOSE(21)
c       CLOSE(23)
c       CLOSE(24)
c       CLOSE(25)
c       CLOSE(26)
c       CLOSE(27)
c       CLOSE(28)
 
      end
      subroutine bftmas(knode,kdgof,kcoor,nsource,
     *iblk,u0,u1,v0,v1,w0,w1,coor,
     *bfu,numa,nnea,mmta,nmta,idet,nodvar,
     *filename,numtyp,knode0,inode,ielem,istop,kend,
     *num0,ilabel,itemp,numblk,stri,str,atemp,idnode)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
      logical filflg
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION  NODVAR(KDGOF,KNODE),COOR(KCOOR,KNODE),R(3),
     &  U0(KDGOF,KNODE),U1(KDGOF,KNODE),W0(KDGOF,KNODE),
     &  V0(KDGOF,KNODE),V1(KDGOF,KNODE),W1(KDGOF,KNODE),
     &  BFU(KDGOF,KNODE)
      dimension NUMA(NUMTYP),NNEA(NUMTYP),MMTA(NUMTYP),NMTA(NUMTYP)
      dimension idet(numtyp),inode(knode0),ielem(num0)
      dimension ilabel(num0*8),itemp(num0*8),atemp(knode0*6)
      dimension stri(6,knode0),str(6,knode0),idnode(knode0)

6       FORMAT (1X,10I5)
7       FORMAT (1X,6E12.5)

      do i=1,num0*8
      ilabel(i)=0
      enddo

c      read(126) (idnode(i),i=1,knode0)
      NUMUNIT = 50
      ISTATUS = 50
      call openf(50,40,ISTATUS)
      call initrwf(50,iui,iur)
      do i = 1,knode0
      iui = iui+1
      idnode(i) = ipool(iui)
      enddo
      call endrwf(50,iui,iur)


      ISTATUS = 35
      call openf(35,6,ISTATUS)


      kfail=0
      do ibk = 0, numblk-1
C...... READ ETYPE0 FILE
C      READ (34) NUMTYPi,NODALL,MATALL
C      read(34) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtypi),
C     &nknode,(inode(i),i=1,nknode)
C      read(34) nelem,(ielem(i),i=1,nelem)

      call initrwf(35,iui,iur)
      iui = iui+1
      numtypi = ipool(iui)
      iui = iui+1
      nodall = ipool(iui)
      iui = iui+1
      matall = ipool(iui)
      iui = iui+1
      nknode = ipool(iui)
      iui = iui+1
      nelem = ipool(iui)
      do i=1,numtypi
      iui = iui+1
      idet(i) = ipool(iui)
      iui = iui+1
      numa(i) = ipool(iui)
      iui = iui+1
      nnea(i) = ipool(iui)
      iui = iui+1
      mmta(i) = ipool(iui)
      iui = iui+1
      nmta(i) = ipool(iui)
      enddo
      do i=1,nknode
      iui = iui+1
      inode(i) = ipool(iui)
      enddo
      do i=1,nelem
      iui = iui+1
      ielem(i) = ipool(iui)
      enddo
      call endrwf(35,iui,iur)


       if (ibk.eq.0) then


c          read(22) nknode,kdgofs
        ISTATUS = 31
        call openf(26,31,ISTATUS)
        call initrwf(26,iui,iur)
        do J=1,6
        do I=1,nknode
        iur = iur+1
        stri(j,i)=rpool(iur)
        enddo
        enddo
        call endrwf(26,iui,iur)



        ISTATUS = 37
        call openf(47,37,ISTATUS)
        call initrwf(47,iui,iur)
        do i = 1, nelem*8
        iui = iui+1
        itemp(i) = ipool(iui)
        end do
        call endrwf(47,iui,iur)

       else


c       receive stress file

        iur0=1
        call recvar(nsource,ibk,atemp(iur0),6*nknode)
        iur=0
        do j=1,6
        do i=1,nknode
        iur = iur+1
        stri(j,i)=atemp(iur)
        enddo
        enddo


c
c       receive label file
c
        iui0=1
        call recvai(nsource,ibk,itemp(iui0),nelem*8)

       endif

c
c
c       combine stress file
c
c
        do i = 1,nknode
          ii = inode(i)
          if (idnode(ii).eq.ibk) then
            do j = 1,6
            str(j,ii) = stri(j,i)
            end do
          endif
        end do


       do i=1,nelem
       do j=1,8
        if (itemp((i-1)*8+j).ge.2) then
         kfail=1
         do k=1,8
         ilabel((ielem(i)-1)*8+k)=1
         enddo
        endif
        if (itemp((i-1)*8+j).eq.1) then
         ilabel((ielem(i)-1)*8+j)=1
        endif
       enddo
       enddo

      enddo


C...... OPEN TIME File
      if (istop.lt.0) then
      OPEN(1,FILE='time',FORM='UNFORMATTED')
      read(1) TMAX,DT,TIME,IT
      CLOSE(1)
      tmax=1.0d0
      istop=0
      else
      ISTATUS = 11
      call openf(51,11,ISTATUS)
      call initrwf(51,iui,iur)
      iur = iur+1
      tmax=rpool(iur)
      iur = iur+1
      dt=rpool(iur)
      iur = iur+1
      time=rpool(iur)
      iui = iui+1
      it=ipool(iui)
      call endrwf(51,iui,iur)
      endif

      if (kend.eq.1) then
      IT = IT+1
      TIME = TIME+DT
      endif

cc      if (kend.eq.1.or.kfail.eq.1) then
c
c     renew unod file
c      
cc      ISTATUS = 33
cc      call openf(43,33,ISTATUS)
cc      call initrwf(43,iui,iur)
cc      do j = 1,kdgof
cc      do i = 1,knode
cc      iur = iur+1
cc      rpool(iur) = 0.d0
cc      end do
cc      end do     
cc      do j = 1,kdgof
cc      do i = 1,knode
cc      iur = iur+1
cc      rpool(iur) = 0.d0
cc      end do
cc      end do
cc      call endrwf(43,iui,iur)
cc      endif

      inquire(file='stop',exist=filflg)
      if (filflg.or.it.ge.50000) then
      istop=1
      tmax=-1.0d0
      OPEN(1,FILE='time',FORM='UNFORMATTED')
      WRITE(1) TMAX,DT,TIME,IT
      CLOSE(1)
      endif
      write(*,*) 'time,dt,it',time,dt,it


      ISTATUS = 35
      call openf(35,6,ISTATUS)

      do ibk = 0, numblk-1

C...... READ ETYPE0 FILE
C      READ (34) NUMTYPi,NODALL,MATALL
C      read(34) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtypi),
C     &nknode,(inode(i),i=1,nknode)
C      read(34) nelem,(ielem(i),i=1,nelem)

      call initrwf(35,iui,iur)
      iui = iui+1
      numtypi = ipool(iui)
      iui = iui+1
      nodall = ipool(iui)
      iui = iui+1
      matall = ipool(iui)
      iui = iui+1
      nknode = ipool(iui)
      iui = iui+1
      nelem = ipool(iui)
      do i=1,numtypi
      iui = iui+1
      idet(i) = ipool(iui)
      iui = iui+1
      numa(i) = ipool(iui)
      iui = iui+1
      nnea(i) = ipool(iui)
      iui = iui+1
      mmta(i) = ipool(iui)
      iui = iui+1
      nmta(i) = ipool(iui)
      enddo
      do i=1,nknode
      iui = iui+1
      inode(i) = ipool(iui)
      enddo
      do i=1,nelem
      iui = iui+1
      ielem(i) = ipool(iui)
      enddo
      call endrwf(35,iui,iur)
 
c
c
c       divide stress file
c
c
        do i = 1,nknode
          ii = inode(i)
          do j = 1,6
          stri(j,i) = str(j,ii)
          end do
        end do


       do i=1,nelem
       do j=1,8
       itemp((i-1)*8+j)=ilabel((ielem(i)-1)*8+j)
       enddo
       enddo

       if (ibk.eq.0) then


c          write(22) nknode,kdgofs
        ISTATUS = 31
        call openf(26,31,ISTATUS)
        call initrwf(26,iui,iur)
        do J=1,6
        do I=1,nknode
        iur = iur+1
        rpool(iur)=stri(j,i)
        enddo
        enddo
        call endrwf(26,iui,iur)



        ISTATUS = 37
        call openf(47,37,ISTATUS)
        call initrwf(47,iui,iur)
        do i = 1, nelem*8
        iui = iui+1
        ipool(iui)=itemp(i)
        end do
        call endrwf(47,iui,iur)

        ISTATUS = 11
        call openf(51,11,ISTATUS)
        call initrwf(51,iui,iur)
        iur = iur+1
        rpool(iur)=tmax
        iur = iur+1
        rpool(iur)=dt
        iur = iur+1
        rpool(iur)=time
        iui = iui+1
        ipool(iui)=it
        call endrwf(51,iui,iur)

       else


c       send stress file

        iur=0
        do j=1,6
        do i=1,nknode
        iur = iur+1
        atemp(iur)=stri(j,i)
        enddo
        enddo
        iur0=1
        call sendar(ibk,nsource,atemp(iur0),6*nknode)


c
c      send label file
c
       iui0=1
       call sendai(ibk,nsource,itemp(iui0),nelem*8)

c
c      send time file
c

       call sendr(ibk,nsource,tmax)
       call sendr(ibk,nsource,dt)
       call sendr(ibk,nsource,time)
       call sendint(ibk,nsource,kfail)
       call sendint(ibk,nsource,it)

       endif

      enddo

C...... READ NV FILE
      ISTATUS = 12
      call openf(23,12,ISTATUS)
c       READ(23) NSKP
      call initrwf(23,iui,iur)
      iui = iui+1
      nskp=ipool(iui)
      call endrwf(23,iui,iur)
c       READ(23) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(23,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      nodvar(I,J)=ipool(iui)
      enddo
      enddo
      call endrwf(23,iui,iur)
CC        WRITE(*,*) 'NV =============='
CC        WRITE(*,6) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C.......READ COOR0.BAK FILE
c       READ(25) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
      call initrwf(25,iui,iur)
      do J=1,KNODE
      do I=1,KCOOR
      iur = iur+1
      coor(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(25,iui,iur)
CC        WRITE(*,*) 'COOR,KCOOR,KNODE ===',KCOOR,KNODE
CC        WRITE(*,7) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
C.......READ UV.DDA FILE
cc      ISTATUS = 14
cc      call openf(26,14,ISTATUS)
c       READ(26) ((U0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((U1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W1(I,J),J=1,KNODE),I=1,KDGOF)
cc      call initrwf(26,iui,iur)
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      u0(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      u1(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      v0(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      v1(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      w0(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      w1(I,J)=rpool(iur)
cc      enddo
cc      enddo
cc      call endrwf(26,iui,iur)
C.......WRITE UNOD FILE
cc      ISTATUS = 18
cc      call openf(27,18,ISTATUS)
c       WRITE(27) ((U0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((U1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W1(I,J),J=1,KNODE),I=1,KDGOF)
cc      call initrwf(27,iui,iur)
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=u0(I,J)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=u1(I,J)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=v0(I,J)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=v1(I,J)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=w0(I,J)
cc      enddo
cc      enddo
cc      do I=1,KDGOF
cc      do J=1,KNODE
cc      iur = iur+1
cc      rpool(iur)=w1(I,J)
cc      enddo
cc      enddo
cc      call endrwf(27,iui,iur)
C...... TO RENEW COOR0
c        DO J=1,KNODE
c        DO I=1,KDGOF
c        IF (KDGOF.GT.1) COOR(I,J)=COOR(I,J)+U0(I,J)
c        ENDDO
c        ENDDO
C.......WRITE COOR0 FILE
cc      ISTATUS = 0
cc      call openf(28,22,ISTATUS)
c       WRITE(28) KNODE,KCOOR
cc      call initrwf(28,iui,iur)
cc      iui = iui+1
cc      ipool(iui)=knode
cc      iui = iui+1
cc      ipool(iui)=kcoor
cc      call endrwf(28,iui,iur)
c       WRITE(28) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
cc      call initrwf(28,iui,iur)
cc      do J=1,KNODE
cc      do I=1,KCOOR
cc      iur = iur+1
cc      rpool(iur)=coor(I,J)
cc      enddo
cc      enddo
cc      call endrwf(28,iui,iur)
CC        WRITE(*,*) 'COOR,KCOOR,KNODE ===',KCOOR,KNODE
CC        WRITE(*,7) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
C...... READ DISP0 FILE
c       READ(21) ((BFU(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(21,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iur = iur+1
      bfu(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(21,iui,iur)
c
c     amend by wanghui
c
      do J=1,KNODE
      do I=1,KDGOF
      bfu(I,J)=bfu(i,j)*(dt/3.15e7)
      enddo
      enddo
      goto 700
        DO 400 N=1,KNODE
          DO J=1,KCOOR
            R(J) = COOR(J,N)
          ENDDO
        DO 500 J=1,KDGOF
      ID = NODVAR(J,N)
      BFU(J,N) = 0.0d0
      IF (ID.LT.0) BFU(J,N) = BOUND(R,TIME,J,dt,it)
C     IF (ID.GT.0) BFU(J,N) = FORCE(R,TIME,J)
c        IF (ID.LT.0)
c     *  BFU(J,N) = BOUNDb(R,TIME,J,DT,N,NBLK)
500     CONTINUE
400     CONTINUE
700     CONTINUE
cc        WRITE(*,*) ' BFU ==========='
cc        WRITE(*,7) ((BFU(J,N),J=1,KDGOF),N=1,KNODE)
C.....  WRITE BFD file
      ISTATUS = 13
      call openf(24,13,ISTATUS)
c       WRITE(24) ((BFU(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(24,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iur = iur+1
      rpool(iur)=bfu(I,J)
      enddo
      enddo
      call endrwf(24,iui,iur)
      RETURN
      END
