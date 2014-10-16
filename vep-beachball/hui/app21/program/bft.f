      subroutine Mbft(iblk,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)

      if (iblk.eq.0) return

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
      call endrwf(21,iui,iur)
c       close(21)
cc        write(*,*) 'numblk,numtyp,nskip,numtypl =='
cc        write(*,*)  numblk,numtyp,nskip,numtypl
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

C...... READ ETYPE0 FILE
      ISTATUS = 6
      call openf(30,6,ISTATUS)
c       READ(30) NUMTYP,NUMNOD,NUMMAT
      call initrwf(30,iui,iur)
      iui = iui+1
      numtyp=ipool(iui)
      iui = iui+1
      numnod=ipool(iui)
      iui = iui+1
      nummat=ipool(iui)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      nskip=ipool(iui)
      call endrwf(30,iui,iur)


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
      knb7=knode
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb8=numnod
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      kna7=kcoor*knode*1
      kna1=kdgof*knode*1
      kna2=kdgof*knode*1
      kna5=kdgof*knode*1
      kna3=kdgof*knode*1
      kna4=kdgof*knode*1
      kna6=kdgof*knode*1
      kna8=kdgof*knode*1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
      kna8=kna8+kna7
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      knb5=knb5+knb4
      knb6=knb6+knb5
      knb7=knb7+knb6
      knb8=knb8+knb7
      call bft(knode,kdgof,kcoor,nsource,
     *iblk,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),aa(kna7),ia(knb0),
     *ia(knb1),ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *filename,numtyp,numnod,ia(knb6),ia(knb7),istop,kend)
 
1000    CONTINUE
c       CLOSE(21)
c       CLOSE(23)
c       CLOSE(24)
c       CLOSE(25)
c       CLOSE(26)
c       CLOSE(27)
c       CLOSE(28)
 
      end
      subroutine bft(knode,kdgof,kcoor,nsource,
     *iblk,u0,u1,v0,v1,w0,w1,coor,
     *bfu,numa,nnea,mmta,nmta,idet,nodvar,
     *filename,numtyp,numnod,inode,ielem,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION  NODVAR(KDGOF,KNODE),COOR(KCOOR,KNODE),R(3),
     &  U0(KDGOF,KNODE),U1(KDGOF,KNODE),W0(KDGOF,KNODE),
     &  V0(KDGOF,KNODE),V1(KDGOF,KNODE),W1(KDGOF,KNODE),
     &  BFU(KDGOF,KNODE)
      dimension NUMA(NUMTYP),NNEA(NUMTYP),MMTA(NUMTYP),NMTA(NUMTYP)
      dimension idet(numtyp),inode(knode),ielem(numnod)

6       FORMAT (1X,10I5)
7       FORMAT (1X,6E12.5)

c ....  read etype0 file
c       read (30) (idet(i),NUMA(I),NNEA(I),
c    &             MMTA(I),NMTA(I),i=1,numtyp)
      call initrwf(30,iui,iur)
      do i=1,numtyp
      iui = iui+1
      idet(i)=ipool(iui)
      iui = iui+1
      numa(I)=ipool(iui)
      iui = iui+1
      nnea(I)=ipool(iui)
      iui = iui+1
      mmta(I)=ipool(iui)
      iui = iui+1
      nmta(I)=ipool(iui)
      enddo
      iui0=iui+1
      iur0=iur+1
      do i=1,knode
      iui = iui+1
      inode(i) = ipool(iui)
      enddo
      iui0=iui+1
      iur0=iur+1
      do i=1,numa(1)
      iui = iui+1
      ielem(i) = ipool(iui)
      enddo
      call endrwf(30,iui,iur)


c
c     send stress file
c

      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call sendar(nsource,iblk,rpool(iur0),knode*6)
      call endrwf(26,iui,iur)


c
c     send label file
c
      ISTATUS = 37
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call sendai(nsource,iblk,ipool(iui0),numa(1)*8)
      call endrwf(47,iui,iur)


c
c     receive stress file
c

      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),knode*6)
      call endrwf(26,iui,iur)



c
c     receive label file
c
      ISTATUS = 37
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),numa(1)*8)
      call endrwf(47,iui,iur)

c
c     receive time file
c

      ISTATUS = 11
      call openf(51,11,ISTATUS)
      call initrwf(51,iui,iur)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      tmax=rpool(iur)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      dt=rpool(iur)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      time=rpool(iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kfail=ipool(iui)
      call recvint(iblk,nsource,ipool(iui))
      it=ipool(iui)
      call endrwf(51,iui,iur)
      if (tmax.lt.0) istop=1

c
c     renew unod file
c
cc      if (kend.eq.1.or.kfail.eq.1) then
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
