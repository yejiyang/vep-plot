      subroutine Mstart(iblk)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
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
c       read (21) lgio,t0,tmax,dt
      call initrwf(21,iui,iur)
      iui = iui+1
      lgio=ipool(iui)
      iur = iur+1
      t0=rpool(iur)
      iur = iur+1
      tmax=rpool(iur)
      iur = iur+1
      dt=rpool(iur)
      call endrwf(21,iui,iur)
c       close(21)
c        write(*,*) 'numblk,numtyp,nskip,numtypl =='
c        write(*,*)  numblk,numtyp,nskip,numtypl
        time=t0
c        write(*,*) 'tmax,dt,time,it ===',tmax,dt,time,it
        it=0
c...... open time file
c       open(21,file='time',form='unformatted',status='unknown')
      ISTATUS = 0
      call openf(21,11,ISTATUS)
c       write(21) tmax,dt,time,it
      call initrwf(21,iui,iur)
      iur = iur+1
      rpool(iur)=tmax
      iur = iur+1
      rpool(iur)=dt
      iur = iur+1
      rpool(iur)=time
      iui = iui+1
      ipool(iui)=it
      call endrwf(21,iui,iur)
c       close(21)
C...... OPNE DISP0 FILE
c       OPEN (21,FILE='sdisp0',FORM='UNFORMATTED',STATUS='OLD')
C...... OPEN ID0 FILE
c       OPEN (22,FILE='sid0',FORM='UNFORMATTED',STATUS='OLD')
C...... OPEN NV FILE
c       OPEN (23,FILE='nv',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN BFD FILE
c       OPEN (24,FILE='bfd',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN UV.DDA file
c       OPEN (26,FILE='uv.dda',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN UVIT.DDA FILE
c       OPEN (27,FILE='uvit.dda',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN COOR00 FILE
c       OPEN (28,FILE='scoor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ORDER.NOD FILE
c       OPEN (29,FILE='order.nod',FORM='UNFORMATTED',STATUS='UNKNOWN')
c...... open end0 file
c       open(10,file='end0',form='formatted',status='unknown')
C...... OPEN COOR0 FILE
c       OPEN (11,FILE='scoor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN unod FILE
c       OPEN (12,FILE='unod',FORM='UNFORMATTED',STATUS='UNKNOWN')
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
      call openf(28,2,ISTATUS)
c        READ(28) NSKIP,KCOOR
      call initrwf(28,iui,iur)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      kcoor=ipool(iui)
      call endrwf(28,iui,iur)
c        WRITE(*,*) 'NSKIP,KCOOR ======',NSKIP,KCOOR
         kend=1
      ISTATUS = 0
      call openf(10,17,ISTATUS)
c        write(10,*) kend
      call initrwf(10,iui,iur)
      iui = iui+1
      ipool(iui)=kend
      call endrwf(10,iui,iur)
 
      knb1=kdgof*knode*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      kna7=kcoor*knode*1
      kna1=kdgof*knode*1
      kna2=kdgof*knode*1
      kna5=kdgof*knode*1
      kna3=kdgof*knode*1
      kna4=kdgof*knode*1
      kna6=kdgof*knode*1
      kna8=kdgof*knode*1
      knb2=knode*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=knode*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
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
      call start(knode,kdgof,kcoor,time,dt,iblk,
     *nsource,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),aa(kna7),ia(knb0),
     *ia(knb1),ia(knb2),
     *filename)
 
1000    CONTINUE
c       CLOSE(21)
c       CLOSE(22)
c       CLOSE(23)
c       CLOSE(24)
c       CLOSE(26)
c       CLOSE(27)
c       CLOSE(28)
c       CLOSE(29)
c       close(10)
c       close(11)
c       close(12)
 
      end
      subroutine start(knode,kdgof,kcoor,time,dt,iblk,
     *nsource,u0,u1,v0,v1,w0,w1,coor,
     *bfu,nodvar,inodvar,jnodvar,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION  NODVAR(KDGOF,KNODE),COOR(KCOOR,KNODE),R(3),
     &  U0(KDGOF,KNODE),U1(KDGOF,KNODE),W0(KDGOF,KNODE),
     &  V0(KDGOF,KNODE),V1(KDGOF,KNODE),W1(KDGOF,KNODE),
     &  BFU(KDGOF,KNODE),INODVAR(KNODE),
     &  JNODVAR(KNODE)
 
6       FORMAT (1X,10I5)
7       FORMAT (1X,6E12.5)
        DO J=1,KNODE
          JNODVAR(J)=0
        ENDDO
      goto 21
C...... READ ORDER.NOD FILE
cc      ISTATUS = 0
cc      call openf(29,16,ISTATUS)
c       READ(29,ERR=21) (INODVAR(I),I=1,KNODE)
cc      call initrwf(29,iui,iur)
cc      do I=1,KNODE
cc      iui = iui+1
cc      inodvar(I)=ipool(iui)
cc      enddo
cc      call endrwf(29,iui,iur)
        DO J=1,KNODE
          INODVAR(J) = J
          JNODVAR(INODVAR(J))=JNODVAR(INODVAR(J))+1
        ENDDO
        DO J=1,KNODE
          IF (JNODVAR(J).EQ.0.OR.JNODVAR(J).GT.1) THEN
          WRITE(*,*)
          WRITE(*,*) 'IT IS ERROR: J,JNODVAR ==='
          WRITE(*,6)  J,JNODVAR(J)
          STOP 0000
          ENDIF
        ENDDO
        GOTO 22
21      CONTINUE
        DO I=1,KNODE
          INODVAR(I)=I
        ENDDO
22      CONTINUE
CC        WRITE(*,*) 'INODVAR ======='
CC        WRITE(*,6) (INODVAR(I),I=1,KNODE)
C...... READ ID0 FILE
      ISTATUS = 3
      call openf(22,3,ISTATUS)
c       READ (22) NSKP,NSKP
      call initrwf(22,iui,iur)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      call endrwf(22,iui,iur)
c       READ (22) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(22,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      nodvar(I,J)=ipool(iui)
      enddo
      enddo
      call endrwf(22,iui,iur)
c       WRITE (*,*) 'ID ='
c       WRITE (*,6) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
        NEQ=0
        DO 20 J=1,KNODE
          J1=INODVAR(J)
          DO 18 I=1,KDGOF
            IF (NODVAR(I,J1).LT.1) GOTO 18
            NEQ = NEQ + 1
            NODVAR(I,J1) = NEQ
18        CONTINUE
20      CONTINUE
        DO 30 J=1,KNODE
          J1=INODVAR(J)
          DO 28 I=1,KDGOF
            IF (NODVAR(I,J1).GE.-1) GOTO 28
            N = -NODVAR(I,J1)-1
27          CONTINUE
            IF (NODVAR(I,N).LT.-1) THEN
              N=-NODVAR(I,N)-1
              GOTO 27
            ENDIF
            NODVAR(I,J1) = NODVAR(I,N)
28        CONTINUE
30      CONTINUE
CC        WRITE(*,*) 'NODVAR SECOND ============'
CC        DO J=1,KNODE
CC        WRITE(*,6) J,(NODVAR(I,J),I=1,KDGOF)
CC        ENDDO
C.......WRITE NV FILE
      ISTATUS = 0
      call openf(23,12,ISTATUS)
c       WRITE(23) NEQ
      call initrwf(23,iui,iur)
      iui = iui+1
      ipool(iui)=neq
      call endrwf(23,iui,iur)
c       WRITE(23) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(23,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      ipool(iui)=nodvar(I,J)
      enddo
      enddo
      call endrwf(23,iui,iur)
CC        WRITE(*,*) 'NV =============='
CC        WRITE(*,6) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C.......READ COOR0.BAK FILE
c       READ(28) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
      call initrwf(28,iui,iur)
      do J=1,KNODE
      do I=1,KCOOR
      iur = iur+1
      coor(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(28,iui,iur)
CC        WRITE(*,*) 'COOR,KCOOR,KNODE ===',KCOOR,KNODE
CC        WRITE(*,7) ((COOR(I,J),I=1,KCOOR),J=1,KNODE)
c        WRITE(*,*) 'T0 ======',T0
        DO 300 N=1,KNODE
        DO 100 J=1,KCOOR
100     R(J) = COOR(J,N)
        DO 200 J=1,KDGOF
        T0=TIME
        T1=TIME-DT
        U0(J,N) = BOUNDOLD(R,T0,J,DT,N,iblk)
        U1(J,N) = BOUNDOLD(R,T1,J,DT,N,iblk)
        V0(J,N) = BOUNDV(R,T0,J,DT,N,iblk)
        V1(J,N) = BOUNDV(R,T1,J,DT,N,iblk)
        W0(J,N) = BOUNDA(R,T0,J,DT,N,iblk)
        W1(J,N) = BOUNDA(R,T1,J,DT,N,iblk)
200     CONTINUE
300     CONTINUE
c        WRITE(*,*) ' U0 ============'
c        WRITE(*,7) ((U0(J,N),J=1,KDGOF),N=1,KNODE)
C        WRITE(*,*) ' U1 ============'
C        WRITE(*,7) ((U1(J,N),J=1,KDGOF),N=1,KNODE)
C.......WRITE UV.DDA FILE
      ISTATUS = 0
      call openf(26,14,ISTATUS)
c       WRITE(26) ((U0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((U1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W1(I,J),J=1,KNODE),I=1,KDGOF)
      call initrwf(26,iui,iur)
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=u0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=u1(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=v0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=v1(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=w0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=w1(I,J)
      enddo
      enddo
      call endrwf(26,iui,iur)
C.......WRITE UV.DDA FILE
      ISTATUS = 0
      call openf(12,18,ISTATUS)
c       WRITE(12) ((U0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((U1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((V1(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W0(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((W1(I,J),J=1,KNODE),I=1,KDGOF)
      call initrwf(12,iui,iur)
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=u0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=u1(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=v0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=v1(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=w0(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=w1(I,J)
      enddo
      enddo
      call endrwf(12,iui,iur)
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
C        WRITE(*,*) ' BFU ==========='
C        WRITE(*,7) ((BFU(J,N),J=1,KDGOF),N=1,KNODE)
C.......WRITE BFD file
      ISTATUS = 0
      call openf(24,13,ISTATUS)
c       WRITE(24) ((BFU(I,J),J=1,KNODE),I=1,KDGOF)
      call initrwf(24,iui,iur)
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=bfu(I,J)
      enddo
      enddo
      call endrwf(24,iui,iur)
        DO I=1,KNODE
          DO J=1,KDGOF
            INV=NODVAR(J,I)
            IF (INV.GT.0) BFU(J,I)=0.0d0
          ENDDO
        ENDDO
C.......WRITE UVIT.DDA FILE
      ISTATUS = 0
      call openf(27,15,ISTATUS)
c       WRITE(27) ((BFU(I,J),J=1,KNODE),I=1,KDGOF),
c    *  ((BFU(I,J),J=1,KNODE),I=1,KDGOF)
      call initrwf(27,iui,iur)
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=bfu(I,J)
      enddo
      enddo
      do I=1,KDGOF
      do J=1,KNODE
      iur = iur+1
      rpool(iur)=bfu(I,J)
      enddo
      enddo
      call endrwf(27,iui,iur)
c...... write coor0 file
      ISTATUS = 2
      call openf(11,2,ISTATUS)
c       write(11) knode,kcoor
      call initrwf(11,iui,iur)
      iui = iui+1
      ipool(iui)=knode
      iui = iui+1
      ipool(iui)=kcoor
      call endrwf(11,iui,iur)
c       write(11) ((coor(i,j),i=1,kcoor),j=1,knode)
      call initrwf(11,iui,iur)
      do j=1,knode
      do i=1,kcoor
      iur = iur+1
      rpool(iur)=coor(i,j)
      enddo
      enddo
      call endrwf(11,iui,iur)
      RETURN
      END
