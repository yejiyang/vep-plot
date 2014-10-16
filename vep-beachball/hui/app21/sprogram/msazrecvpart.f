      subroutine Mmsazrecvpart(istop)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      character*12 fname1,fname2
      include 'partdata.h'
      logical filflgdisp(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
C
c      if(iblk.gt.0) return
      nsource=0
      maxlm=0
      idisp1 = 0
      idisp2 = 0
      open(1,file='mlmddm',form='unformatted',status='unknown')
      read(1) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
     & idisp1,idisp2
      read(1) lgio,t0,tmax,dt,num0,knode0
      close(1)
c
      ISTATUS = 0
      call openf(1,1,ISTATUS)
      call initrwf(1,iui,iur)
      iui = iui+1
      ipool(iui) = numblk
      iui = iui+1
      ipool(iui) = numtyp
      iui = iui+1
      ipool(iui) = maxlm
      iui = iui+1
      ipool(iui) = nmdof
      iui = iui+1
      ipool(iui) = numtypl
      iui = iui+1
      ipool(iui) = kdgofl
      iui = iui+1
      ipool(iui) = keymt
      iui = iui+1
      ipool(iui) = idisp1
      iui = iui+1
      ipool(iui) = idisp2
      iui = iui+1
      ipool(iui) = lgio
      iur = iur+1
      rpool(iur) = t0
      iur = iur+1
      rpool(iur) = tmax
      iur = iur+1
      rpool(iur) = dt
      iui = iui+1
      ipool(iui) = num0
      iui = iui+1
      ipool(iui) = knode0
      call endrwf(1,iui,iur)
C
c
c     open data file from master processor
c
C...... OPEN LMDDM FILE
      open(40,file='mlmddm',form='unformatted',status='unknown')
 
C...... OPEN COOR0 FILE
      OPEN (31,FILE='mcoor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPEN ID0 FILE
      OPEN (32,FILE='mid0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPNE DISP0 FILE (Boundary condition file)
      OPEN (33,FILE='mdisp0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPNE DISP1 FILE (Initial value file for displacement)
      OPEN (34,FILE='mdisp1',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C...... OPEN ELEMENT TYPE FILE FOR EACH SUBDOMAIN (etype0)
      OPEN (35,file='metype0',form='unformatted',status='unknown')
C
C
C...... OPEN ELEM0 FILE
        OPEN (30,FILE='melem0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C
C...... OPEN SSUBNODES0 file for node index of each subdomain
        OPEN (39,file='subnodes0',form='unformatted',status='unknown')
        open(123,file='u1v1',form='unformatted',status='unknown')
        open(124,file='fstr1',form='unformatted',status='unknown')
        open(126,file='nmbnodes0',form='unformatted',status='unknown')
C
C
C     open data file for other initial value problem
C
C...... OPNE DISP2 FILE (Initial value file for velocity)
      if (idisp1.eq.1) then
      open(37,file='mdisp2',form='unformatted',status='unknown')
      END IF
C...... OPNE DISP3 FILE (Initial value file for acceleration)
      if (idisp2.eq.1) then
      open(38,file='mdisp3',form='unformatted',status='unknown')
      END IF
C
C
      DO 1000 IBLK = 0,0
C     READ LMDDM FILE
      read(40) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
     & idisp1,idisp2
      read(40) lgio,t0,tmax,dt,nskip,nskip
c
c
      ISTATUS = 1
      call openf(1,1,ISTATUS)
      call initrwf(1,iui,iur)
      iui = iui+1
      ipool(iui) = numblk
      iui = iui+1
      ipool(iui) = numtyp
      iui = iui+1
      ipool(iui) = maxlm
      iui = iui+1
      ipool(iui) = nmdof
      iui = iui+1
      ipool(iui) = numtypl
      iui = iui+1
      ipool(iui) = kdgofl
      iui = iui+1
      ipool(iui) = keymt
      iui = iui+1
      ipool(iui) = idisp1
      iui = iui+1
      ipool(iui) = idisp2
      iui = iui+1
      ipool(iui) = lgio
      iur = iur+1
      rpool(iur) = t0
      iur = iur+1
      rpool(iur) = tmax
      iur = iur+1
      rpool(iur) = dt
      call endrwf(1,iui,iur)
C
      nparts = numblk
c
C
C     READ COOR0 FILE
      NUMUNIT = 31
      READ (31) NUMNOD,NCOOR
      ISTATUS = 0
      call openf(31,2,ISTATUS)
      call initrwf(31,iui,iur)
      iui = iui+1
      ipool(iui) = numnod
      iui = iui+1
      ipool(iui) = ncoor
      call endrwf(31,iui,iur)
      CALL MREADCOOR(AA(1),NUMNOD,NCOOR,NUMUNIT,IBLK)
C
C     READ ID0 FILE
      NUMUNIT = 32
      READ (32) KNODE,KDGOF
      ISTATUS = 0
      call openf(32,3,ISTATUS)
      call initrwf(32,iui,iur)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = kdgof
      call endrwf(32,iui,iur)
      CALL MREADID(IA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP0 FILE
      NUMUNIT = 33
      READ (33) KNODE,KDGOF
      ISTATUS = 0
      call openf(33,4,ISTATUS)
      call initrwf(33,iui,iur)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = kdgof
      call endrwf(33,iui,iur)
      CALL MREADDISP0(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP1 FILE
C
      NUMUNIT = 34
      READ (34) KNODE,KDGOF
      ISTATUS = 0
      call openf(34,5,ISTATUS)
      call initrwf(34,iui,iur)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = kdgof
      call endrwf(34,iui,iur)
      CALL MREADDISP1(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
C
C     READ DISP2 FILE
C
      if (idisp1.eq.1) then
      NUMUNIT = 37
      READ (37) KNODE,KDGOF
      ISTATUS = 0
      call openf(37,9,ISTATUS)
      call initrwf(37,iui,iur)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = kdgof
      call endrwf(37,iui,iur)
      CALL MREADDISP2(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
      endif
C
C     READ DISP3 FILE
C
      if (idisp2.eq.1) then
      NUMUNIT = 38
      READ (38) KNODE,KDGOF
      ISTATUS = 0
      call openf(38,10,ISTATUS)
      call initrwf(38,iui,iur)
      iui = iui+1
      ipool(iui) = knode
      iui = iui+1
      ipool(iui) = kdgof
      call endrwf(38,iui,iur)
      CALL MREADDISP3(AA(1),KNODE,KDGOF,NUMUNIT,IBLK)
      endif
C
C     READ ETYPE0 FILE
      NUMUNIT = 35
      ISTATUS = 0
      call openf(35,6,ISTATUS)

      do imas = 0, numblk-1
      READ (35) NUMTYP,NODALL,MATALL
      call initrwf(35,iui,iur)
      iui = iui+1
      ipool(iui) = numtyp
      iui = iui+1
      ipool(iui) = nodall
      iui = iui+1
      ipool(iui) = matall
      call endrwf(35,iui,iur)
      K1 = 1
      K2 = K1 + KNODE0 + 10
      K3 = K2 + NUMTYP + 10
      K4 = K3 + NUMTYP + 10
      K5 = K4 + NUMTYP + 10
      K6 = K5 + NUMTYP + 10
      K7 = K6 + NUMTYP + 10
      K8 = K7 + NUM0+10
      CALL MREADEYTPE(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),IA(K6),
     & NUMTYP,KNODE,knode0,KDGOF,NUMUNIT,IBLK,NUM0,ia(k7),nummas)
      if (imas.eq.0) then
      num=nummas
      nmtpbak=numtyp
      ndllbak=nodall
      mtllbak=matall
      endif
      enddo
      numtyp=nmtpbak
      nodall=ndllbak
      matall=mtllbak

C
C     READ ELEM0 FILE
      NUMUNIT = 30
      K8 = K7 + NODALL
      KVAR=KNODE*KDGOF
        kelem=NODALL
        CALL MREADELEM(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),
     &IA(K6),IA(K7),AA(1),NUMTYP,KELEM,MATALL,NUMUNIT,KNODE,IBLK)
C
C     READ SSUBNODES0 FILE
      NUMUNIT = 39
      READ(39) KNODEALL,KDGOFALL
      ISTATUS = 0
      call openf(39,8,ISTATUS)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui) = knodeall
      iui = iui+1
      ipool(iui) = kdgofall
      call endrwf(39,iui,iur)
      K1 = 1
      K2 = K1 + KNODEALL + 10
      K3 = K2 + KNODEALL + 10
      K4 = K3 + KNODEALL + 10
      K5 = K4 + KNODEALL + 10
      K6 = K6 + KNODEALL*KDGOFALL + 10
 
      CALL MREADSSUBNODES(IA(K1),IA(K2),IA(K3),IA(K4),IA(K5),
     & KNODEALL,KDGOFALL,IBLK)

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
      read(124) (rpool(iur0-1+i),i=1,6*KNODE)
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
      read(124) (rpool(iur0-1+i),i=1,num*6*8)
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
      read(123) (rpool(iur0-1+i),i=1,knode*kdgof)
      endif
      iur0=iur+1
      do i = 1,1*knode*kdgof
      iur = iur+1
      rpool(iur) = 0.d0
      enddo
      if (istop.lt.0) then
      read(123) (rpool(iur0-1+i),i=1,knode*kdgof)
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

C
C     idnode used in nzrecvstr.f
C
      NUMUNIT = 50
      ISTATUS = 0
      call openf(50,40,ISTATUS)
      call initrwf(50,iui,iur)
      iui0=iui+1
      do i = 1,knode0
      iui = iui+1
      ipool(iui) = 0
      enddo
      read(126) (ipool(iui0-1+i),i=1,knode0)
      call endrwf(50,iui,iur)

C
      call msazrecvpart(iblk,
     *filename)
C
1000  CONTINUE
      close(40)
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
      close(37)
      close(38)
      close(39)
      close(123)
      close(124)
      close(126)
 
      end

      subroutine msazrecvpart(iblk,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
         include 'memalloc.h'
c         WRITE(*,*) 'RECV PARTITION DATA OK..........'
       RETURN
       END
 
      SUBROUTINE  MREADCOOR(COOR,NUMNOD,NCOOR,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      DIMENSION COOR(NCOOR,NUMNOD)
C
C
      nsource = 0
C.......OPEN AND READ COOR file
        READ (31) ((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
c        WRITE(*,*) 'NUMNOD,NCOOR=',NUMNOD,NCOOR
c        WRITE(*,*)
c        DO I = 1, NUMNOD
c        WRITE(*,1100) I,(COOR(J,I),J=1,NCOOR)
c        END DO
      call initrwf(31,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,NUMNOD
      do I=1,NCOOR
      iur = iur+1
      rpool(iur) = coor(I,J)
      enddo
      enddo
      call endrwf(31,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  MREADID(NODVAR,KNODE,KDGOF,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION NODVAR(KDGOF,KNODE)
C.......OPEN AND READ ID FILE
        NUMNOD = KNODE
       nsource = 0
       READ (32) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C        WRITE(*,*) 'KNODE,KDGOF=',NUMNOD,KDGOF
C        WRITE(*,*)
C        DO I = 1, KNODE
C        WRITE(*,1000) I,(NODVAR(J,I),J=1,KDGOF)
C        END DO
C
      call initrwf(32,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      ipool(iui) = nodvar(I,J)
      enddo
      enddo
      call endrwf(32,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
 
       SUBROUTINE  MREADDISP0(BFU,KNODE,KDGOF,NUMUNIT,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION BFU(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
        READ (33) ((BFU(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(BFU(J,I),J=1,NODDOF)
C        END DO
C
      call initrwf(33,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      rpool(iur) = bfu(I,J)
      enddo
      enddo
      call endrwf(33,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  MREADDISP1(U0,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION U0(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
        READ (34) ((U0(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U0(J,I),J=1,NODDOF)
C        END DO
C
      call initrwf(34,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      rpool(iur) = u0(I,J)
      enddo
      enddo
      call endrwf(34,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  MREADDISP2(U1,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION U1(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
        READ (37) ((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U1(J,I),J=1,NODDOF)
C        END DO
C
      call initrwf(37,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      rpool(iur) = u1(I,J)
      enddo
      enddo
      call endrwf(37,iui,iur)
 
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
       SUBROUTINE  MREADDISP3(U2,KNODE,KDGOF,NUMUNIT,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION U2(KDGOF,KNODE)
       NUMNOD = KNODE
       NODDOF = KDGOF
       nsource = 0
C.......OPEN AND READ DISP0 FILE
        READ (38) ((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(U2(J,I),J=1,NODDOF)
C        END DO
C
      call initrwf(38,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,NUMNOD
      do I=1,NODDOF
      iur = iur+1
      rpool(iur) = u2(I,J)
      enddo
      enddo
      call endrwf(38,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        return
        end
 
      SUBROUTINE MREADEYTPE(INODE,NUMA,NNEA,MMTA,NMTA,
     & IDET,NUMTYP,KNODE,knode0,KDGOF,NUMUNIT,IBLK,num0,ielem,num)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION  inode(knode0),
     &  numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp)
       dimension ielem(num0)
C
       nsource = 0
C...  READ etype0 file
      read(35) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
     &nknode,(inode(i),i=1,nknode)
      read(35) nskip,(ielem(i),i=1,numa(1))
C
      call initrwf(35,iui,iur)
      iui0=iui+1
      iur0=iur+1
      iui = iui+1
      ipool(iui) = nknode
      iui = iui+1
      ipool(iui) = numa(1)
      do i=1,numtyp
      iui = iui+1
      ipool(iui) = idet(i)
      iui = iui+1
      ipool(iui) = numa(i)
      iui = iui+1
      ipool(iui) = nnea(i)
      iui = iui+1
      ipool(iui) = mmta(i)
      iui = iui+1
      ipool(iui) = nmta(i)
      enddo
      iui0=iui+1
      iur0=iur+1
      do i=1,nknode
      iui = iui+1
      ipool(iui) = inode(i)
      enddo
      iui0=iui+1
      iur0=iur+1
      do i=1,numa(1)
      iui = iui+1
      ipool(iui) = ielem(i)
      enddo
      call endrwf(35,iui,iur)
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
        num=numa(1)
        return
        end
 
       SUBROUTINE MREADELEM(INODE,NUMA,NNEA,MMTA,NMTA,
     & IDET,NODEALL,EMATEALL,NUMTYP,KELEM,MATALL,NUMUNIT,KNODE,IBLK)
       implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION  inode(knode),
     &  numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),nodeall(kelem),emateall(matall)
 
       nsource = 0
       nodall = kelem
C ... READ elem0 file
      read(30) (nodeall(i),i=1,nodall)
      read(30) (emateall(i),i=1,matall)
C
      ISTATUS = 0
      call openf(30,7,ISTATUS)
      call initrwf(30,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do i=1,nodall
      iui = iui+1
      ipool(iui) = nodeall(i)
      enddo
      call endrwf(30,iui,iur)
      call initrwf(30,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do i=1,matall
      iur = iur+1
      rpool(iur) = emateall(i)
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
 
 
      SUBROUTINE MREADSSUBNODES(idnode,inodeall,jnode,idnode1,nodvarall,
     & KNODEALL,KDGOFALL,IBLK)
      implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      DIMENSION idnode(knodeall),inodeall(knodeall),jnode(knodeall),
     &  idnode1(knodeall),nodvarall(kdgofall,knodeall)
C
C     READ IDL FILE
       nsource = 0
C........read index of subdomain nodes file
         read(39) nodesiblk,(idnode(i),i=1,1)
         read(39) nod_iblk,(inodeall(i),i = 1,1)
         read(39) knodeall,(jnode(i),i=1,1)
         read(39) knodeall,(idnode1(i),i=1,1)
         read(39) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
     &                                ,i=1,nod_iblk)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui)= nodesiblk
      iui0=iui+1
      iur0=iur+1
      do i=1,1
      iui = iui+1
      ipool(iui) = idnode(i)
      enddo
      call endrwf(39,iui,iur)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui) = nod_iblk
      iui0=iui+1
      iur0=iur+1
      do i=1,1
      iui = iui+1
      ipool(iui) = inodeall(i)
      enddo
      call endrwf(39,iui,iur)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui) = knodeall
      iui0=iui+1
      iur0=iur+1
      do i=1,1
      iui = iui+1
      ipool(iui) = jnode(i)
      enddo
      call endrwf(39,iui,iur)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui) = knodeall
      iui0=iui+1
      iur0=iur+1
      do i=1,1
      iui = iui+1
      ipool(iui) = idnode1(i)
      enddo
      call endrwf(39,iui,iur)
c        read(39) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c    &                                ,i=1,knodeall)
      call initrwf(39,iui,iur)
      iui = iui+1
      ipool(iui) = knodeall
      iui = iui+1
      ipool(iui) = kdgofall
      iui0=iui+1
      iur0=iur+1
      do i=1,nod_iblk
      do j=1,kdgofall
      iui = iui+1
      ipool(iui) = nodvarall(j,i)
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
 
