      subroutine Mazrecvpart(iblk,istop)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      character*12 fname1,fname2
      include 'partdata.h'
      logical filflgdisp(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
C
      if(iblk.eq.0) return
      nsource=0
      maxlm=0
      idisp1 = 0
      idisp2 = 0
c      open(1,file='lmddm',form='unformatted',status='unknown')
c      read(1) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
c     & idisp1,idisp2
c      read(1) lgio,t0,tmax,dt
c      close(1)
c      nparts = numblk
C
c
c     open data file from master processor
c
C...... OPEN LMDDM FILE
c     open(1,file='lmddm',form='unformatted',status='unknown')
 
C...... OPEN COOR0 FILE
c     OPEN (31,FILE='scoor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPEN ID0 FILE
c     OPEN (32,FILE='sid0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPNE DISP0 FILE (Boundary condition file)
c     OPEN (33,FILE='sdisp0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPNE DISP1 FILE (Initial value file for displacement)
c     OPEN (34,FILE='sdisp1',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPEN ELEMENT TYPE FILE FOR EACH SUBDOMAIN (etype0)
c     OPEN (35,file='setype0',form='unformatted',status='unknown')
C
C
C...... OPEN ELEM0 FILE
c       OPEN (30,FILE='selem0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C
C...... OPEN SSUBNODES0 file for node index of each subdomain
c       OPEN (39,file='ssubnodes0',form='unformatted',status='unknown')
C
C
C     open data file for other initial value problem
C
C...... OPNE DISP2 FILE (Initial value file for velocity)
      if (idisp1.eq.1) then
c     open(37,file='sdisp2',form='unformatted',status='unknown')
      END IF
C...... OPNE DISP3 FILE (Initial value file for acceleration)
      if (idisp2.eq.1) then
c     open(38,file='sdisp3',form='unformatted',status='unknown')
      END IF
C
C
c      DO 1000 IBLK = 0,numblk-1
C     READ LMDDM FILE
      ISTATUS = 0
      call openf(1,1,ISTATUS)
c     read(1) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
c    & idisp1,idisp2
      call initrwf(1,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      numblk=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      numtyp=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      maxlm=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      nmdof=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      numtypl=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgofl=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      keymt=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      idisp1=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      idisp2=ipool(iui)
      call endrwf(1,iui,iur)
c     read(1) lgio,t0,tmax,dt
      call initrwf(1,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      lgio=ipool(iui)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      t0=rpool(iur)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      tmax=rpool(iur)
      iur = iur+1
      call recvr(iblk,nsource,rpool(iur))
      dt=rpool(iur)
      call endrwf(1,iui,iur)
C
C     READ COOR0 FILE
      NUMUNIT = 31
      ISTATUS = 0
      call openf(31,2,ISTATUS)
c     READ (31) NUMNOD,NCOOR
      call initrwf(31,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      numnod=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      ncoor=ipool(iui)
      call endrwf(31,iui,iur)
      CALL READCOOR(AA(1),NUMNOD,NCOOR,NUMUNIT,IBLK)
C
C     READ ID0 FILE
      NUMUNIT = 32
      ISTATUS = 0
      call openf(32,3,ISTATUS)
c     READ (32) KNODE,KDGOF
      call initrwf(32,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knode=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgof=ipool(iui)
      call endrwf(32,iui,iur)
      CALL READID(IA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP0 FILE
      NUMUNIT = 33
      ISTATUS = 0
      call openf(33,4,ISTATUS)
c     READ (33) KNODE,KDGOF
      call initrwf(33,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knode=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgof=ipool(iui)
      call endrwf(33,iui,iur)
      CALL READDISP0(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP1 FILE
      NUMUNIT = 34
      ISTATUS = 0
      call openf(34,5,ISTATUS)
c     READ (34) KNODE,KDGOF
      call initrwf(34,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knode=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgof=ipool(iui)
      call endrwf(34,iui,iur)
      CALL READDISP1(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP2 FILE
      if (idisp1.eq.1) then
      NUMUNIT = 37
      ISTATUS = 0
      call openf(37,9,ISTATUS)
c     READ (37) KNODE,KDGOF
      call initrwf(37,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knode=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgof=ipool(iui)
      call endrwf(37,iui,iur)
      CALL READDISP2(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
      endif
C
C     READ DISP3 FILE
      if (idisp2.eq.1) then
      NUMUNIT = 38
      ISTATUS = 0
      call openf(38,10,ISTATUS)
c     READ (38) KNODE,KDGOF
      call initrwf(38,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knode=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgof=ipool(iui)
      call endrwf(38,iui,iur)
      CALL READDISP3(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
      endif
C
C     READ ETYPE0 FILE
      NUMUNIT = 35
      ISTATUS = 0
      call openf(35,6,ISTATUS)
c     READ (35) NUMTYP,NODALL,MATALL
      call initrwf(35,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      numtyp=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      nodall=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      matall=ipool(iui)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = knode
      call endrwf(35,iui,iur)
      K1 = 1
      K2 = K1 + KNODE + 10
      K3 = K2 + NUMTYP + 10
      K4 = K3 + NUMTYP + 10
      K5 = K4 + NUMTYP + 10
      K6 = K5 + NUMTYP + 10
      K7 = K6 + NUMTYP + 10
      K8 = K7 + NODALL/6+10
      CALL READEYTPE(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),
     & IA(K6),NUMTYP,KNODE,KDGOF,NUMUNIT,IBLK,nodall,ia(k7),num)
C
C     READ ELEM0 FILE
      NUMUNIT = 30
      K8 = K7 + NODALL
      KVAR=KNODE*KDGOF
        kelem=NODALL
        CALL READELEM(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),
     &IA(K6),IA(K7),AA(1),NUMTYP,KELEM,MATALL,NUMUNIT,KNODE,IBLK)
C
C     READ SSUBNODES0 FILE
      NUMUNIT = 39
      ISTATUS = 0
      call openf(39,8,ISTATUS)
c     READ(39) KNODEALL,KDGOFALL
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knodeall=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgofall=ipool(iui)
      call endrwf(39,iui,iur)
      K1 = 1
      K2 = K1 + KNODEALL + 10
      K3 = K2 + KNODEALL + 10
      K4 = K3 + KNODEALL + 10
      K5 = K4 + KNODEALL + 10
      K6 = K6 + KNODEALL*KDGOFALL + 10
 
      CALL READSSUBNODES(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),
     & KNODEALL,KDGOFALL,IBLK)
C

C
C     principle stress at nodes at previous time step
C

      ISTATUS = 0
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui = iui+1
      ipool(iui)=knode
      iui = iui+1
      ipool(iui)=6
      call endrwf(26,iui,iur)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,6
      do I=1,KNODE
      iur = iur+1
      rpool(iur)=0.d0
      enddo
      enddo
      if (istop.lt.0) then
      call recvar(iblk,nsource,rpool(iur0),1*6*KNODE)
      endif
      call endrwf(26,iui,iur)


C
C     strw6.tmp for stress at previous time step
C
      NUMUNIT = 42
      ISTATUS = 0
      call openf(42,32,ISTATUS)
      call initrwf(42,iui,iur)
      iur0=iur+1
      do i = 1,6*NODALL
      iur = iur+1
      rpool(iur) = 0.d0
      enddo
      if (istop.lt.0) then
      call recvar(iblk,nsource,rpool(iur0),num*6*8)
      endif
      call endrwf(42,iui,iur)

C
C     unod for displacement and velocity at previous time step
C
      NUMUNIT = 43
      ISTATUS = 0
      call openf(43,33,ISTATUS)
      call initrwf(43,iui,iur)
      iur0=iur+1
      do i = 1,1*knode*kdgof
      iur = iur+1
      rpool(iur) = 0.d0
      enddo
      if (istop.lt.0) then
      call recvar(iblk,nsource,rpool(iur0),1*knode*kdgof)
      endif
      iur0=iur+1
      do i = 1,1*knode*kdgof
      iur = iur+1
      rpool(iur) = 0.d0
      enddo
      if (istop.lt.0) then
      call recvar(iblk,nsource,rpool(iur0),1*knode*kdgof)
      endif
      call endrwf(43,iui,iur)

C
C     label for failure switch
C
      NUMUNIT = 47
      ISTATUS = 0
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      iui0=iui+1
      do i = 1,NODALL
      iui = iui+1
      ipool(iui) = 0
      enddo
      call endrwf(47,iui,iur)

C
C     iteration steps
C
      NUMUNIT = 48
      ISTATUS = 0
      call openf(48,38,ISTATUS)
      call initrwf(48,iui,iur)
      iui = iui+1
      ipool(iui) = 0
      if (istop.lt.0) then
      istop=0
      endif
      call endrwf(48,iui,iur)

C
C     energy release accumulation
C
      NUMUNIT = 49
      ISTATUS = 0
      call openf(49,39,ISTATUS)
      call initrwf(49,iui,iur)
      do i = 1,2*NODALL
      iur = iur+1
      rpool(iur) = 0.d0
      enddo
      call endrwf(49,iui,iur)

      call azrecvpart(iblk,
     *filename)
C
1000  CONTINUE
c     close(1)
c     close(30)
c     close(31)
c     close(32)
c     close(33)
c     close(34)
c     close(35)
c     close(36)
c     close(37)
c     close(38)
c     close(39)
 
      end
      subroutine azrecvpart(iblk,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
         include 'memalloc.h'
c         WRITE(*,*) 'RECV PARTITION DATA OK..........'
       RETURN
       END
 
      SUBROUTINE  READCOOR(COOR,NUMNOD,NCOOR,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      DIMENSION COOR(NCOOR,NUMNOD)
C
C
      nsource = 0
C.......OPEN AND READ COOR file
c       READ (31) ((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
      call initrwf(31,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*NUMNOD*NCOOR)
      do J=1,NUMNOD
      do I=1,NCOOR
      iur = iur+1
      coor(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(31,iui,iur)
c        WRITE(*,*) 'NUMNOD,NCOOR=',NUMNOD,NCOOR
c        WRITE(*,*)
c        DO I = 1, NUMNOD
c        WRITE(*,1100) I,(COOR(J,I),J=1,NCOOR)
c        END DO
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  READID(NODVAR,KNODE,KDGOF,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION NODVAR(KDGOF,KNODE)
C.......OPEN AND READ ID FILE
        NUMNOD = KNODE
       nsource = 0
c      READ (32) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(32,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*KNODE*KDGOF)
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      nodvar(I,J)=ipool(iui)
      enddo
      enddo
      call endrwf(32,iui,iur)
C        WRITE(*,*) 'KNODE,KDGOF=',NUMNOD,KDGOF
C        WRITE(*,*)
C        DO I = 1, KNODE
C        WRITE(*,1000) I,(NODVAR(J,I),J=1,KDGOF)
C        END DO
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
 
       SUBROUTINE  READDISP0(BFU,KNODE,KDGOF,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION BFU(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
c       READ (33) ((BFU(I,J),I=1,NODDOF),J=1,NUMNOD)
      call initrwf(33,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*NUMNOD*NODDOF)
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      bfu(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(33,iui,iur)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(BFU(J,I),J=1,NODDOF)
C        END DO
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  READDISP1(U0,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION U0(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
c       READ (34) ((U0(I,J),I=1,NODDOF),J=1,NUMNOD)
      call initrwf(34,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*NUMNOD*NODDOF)
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      u0(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(34,iui,iur)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U0(J,I),J=1,NODDOF)
C        END DO
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  READDISP2(U1,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION U1(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
c       READ (37) ((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
      call initrwf(37,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*NUMNOD*NODDOF)
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      u1(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(37,iui,iur)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U1(J,I),J=1,NODDOF)
C        END DO
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  READDISP3(U2,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION U2(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
c       READ (38) ((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
      call initrwf(38,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*NUMNOD*NODDOF)
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      u2(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(38,iui,iur)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U2(J,I),J=1,NODDOF)
C        END DO
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
      SUBROUTINE READEYTPE(INODE,NUMA,NNEA,MMTA,NMTA,
     & IDET,NUMTYP,KNODE,KDGOF,NUMUNIT,IBLK,nodall,ielem,num)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION  inode(knode),
     &  numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp)
       dimension ielem(nodall/6)
C
       nsource = 0
C...  READ etype0 file
c     read(35) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
c    &(inode(i),i=1,knode)
      call initrwf(35,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),5*numtyp)
      do i=1,numtyp
      iui = iui+1
      idet(i)=ipool(iui)
      iui = iui+1
      numa(i)=ipool(iui)
      iui = iui+1
      nnea(i)=ipool(iui)
      iui = iui+1
      mmta(i)=ipool(iui)
      iui = iui+1
      nmta(i)=ipool(iui)
      enddo
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*knode)
      do i=1,knode
      iui = iui+1
      inode(i)=ipool(iui)
      enddo
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*numa(1))
      do i=1,numa(1)
      iui = iui+1
      ielem(i)=ipool(iui)
      enddo
      call endrwf(35,iui,iur)
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        num=numa(1)
        return
        end
 
       SUBROUTINE READELEM(INODE,NUMA,NNEA,MMTA,NMTA,
     & IDET,NODEALL,EMATEALL,NUMTYP,KELEM,MATALL,NUMUNIT,KNODE,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
       DIMENSION  inode(knode),
     &  numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),nodeall(kelem),emateall(matall)
 
       nsource = 0
       nodall = kelem
C ... READ elem0 file
      ISTATUS = 0
      call openf(30,7,ISTATUS)
c     read(30) (nodeall(i),i=1,nodall)
      call initrwf(30,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*nodall)
      do i=1,nodall
      iui = iui+1
      nodeall(i)=ipool(iui)
      enddo
      call endrwf(30,iui,iur)
c     read(30) (emateall(i),i=1,matall)
      call initrwf(30,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(iblk,nsource,rpool(iur0),1*matall)
      do i=1,matall
      iur = iur+1
      emateall(i)=rpool(iur)
      enddo
      call endrwf(30,iui,iur)
C
      return
      nadd = 0
      madd = 0
      do itype = 1,numtyp
      mmate = mmta(itype)
      nmate = nmta(itype)
      num = numa(itype)
      nnode = nnea(itype)
      WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
      WRITE(*,*) 'NODE ='
      DO I = 1,NUM
      WRITE(*,1000) I,(NODEALL((I-1)*NNODE+J+nadd),J=1,NNODE)
      END DO
      NADD = NADD + NUM * NNODE
      WRITE(*,*) ' MATE FOR EACH ELEMENT TYPE'
      WRITE(*,*) 'MMATE ====',MMATE,'NMATE====',NMATE
      WRITE(*,*) 'MATE ======'
      DO I = 1,MMATE
      WRITE(*,1100) I,(EMATEALL((I-1)*NMATE+J+MADD),J = 1,NMATE)
      END DO
      MADD = MADD + MMATE*NMATE
      END DO
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
       return
       end
 
 
      SUBROUTINE READSSUBNODES(idnode,inodeall,jnode,idnode1,nodvarall,
     & KNODEALL,KDGOFALL,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      DIMENSION idnode(knodeall),inodeall(knodeall),jnode(knodeall),
     &  idnode1(knodeall),nodvarall(kdgofall,knodeall)
C
C     READ IDL FILE
       nsource = 0
C........read index of subdomain nodes file
c        read(39) nodesiblk,(idnode(i),i=1,nodesiblk)
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      nodesiblk=ipool(iui)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*1)
      do i=1,1
      iui = iui+1
      idnode(i)=ipool(iui)
      enddo
      call endrwf(39,iui,iur)
c        read(39) nod_iblk,(inodeall(i),i = 1,nod_iblk)
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      nod_iblk=ipool(iui)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*1)
      do i=1,1
      iui = iui+1
      inodeall(i)=ipool(iui)
      enddo
      call endrwf(39,iui,iur)
c        read(39) knodeall,(jnode(i),i=1,knodeall)
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knodeall=ipool(iui)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*1)
      do i=1,1
      iui = iui+1
      jnode(i)=ipool(iui)
      enddo
      call endrwf(39,iui,iur)
c        read(39) knodeall,(idnode1(i),i=1,knodeall)
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knodeall=ipool(iui)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*1)
      do i=1,1
      iui = iui+1
      idnode1(i)=ipool(iui)
      enddo
      call endrwf(39,iui,iur)
c        read(39) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c    &                                ,i=1,nod_iblk)
      call initrwf(39,iui,iur)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      knodeall=ipool(iui)
      iui = iui+1
      call recvint(iblk,nsource,ipool(iui))
      kdgofall=ipool(iui)
      iui0=iui+1
      iur0=iur+1
      call recvai(iblk,nsource,ipool(iui0),1*nod_iblk*kdgofall)
      do i=1,nod_iblk
      do j=1,kdgofall
      iui = iui+1
      nodvarall(j,i)=ipool(iui)
      enddo
      enddo
      call endrwf(39,iui,iur)
C
C         WRITE(*,*) ' subnodes '
c         write(*,*) knodeall,kdgofall
c         write(*,*) nodesiblk,(idnode(i),i=1,nodesiblk)
c         write(*,*) nod_iblk,(inodeall(i),i = 1,nod_iblk)
c         write(*,*) knodeall,(jnode(i),i=1,knodeall)
c         write(*,*) knodeall,(idnode1(i),i=1,knodeall)
c         write(*,*) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c     &                                ,i=1,knodeall)
C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
