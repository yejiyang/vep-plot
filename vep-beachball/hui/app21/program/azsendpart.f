      subroutine Mazsendpart(istop)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      character*12 fname1,fname2
      include 'memalloc.h'
      logical filflgdisp(20)
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
      nsource=0
      maxlm=0
      idisp1 = 0
      idisp2 = 0
      open(1,file='mlmddm',form='unformatted',status='unknown')
      read(1) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
     & idisp1,idisp2
      read(1) lgio,t0,tmax,dt,nskip,nskip
      close(1)

      nparts = numblk
c
c     open data file from master processor
c
C...... OPEN COOR0 FILE
        OPEN (30,FILE='mcoor0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ID0 FILE
        OPEN (31,FILE='mid0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP0 FILE (Boundary condition file)
        OPEN (32,FILE='mdisp0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP1 FILE (Initial value file for displacement)
        OPEN (33,FILE='mdisp1',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ELEMENT TYPE FILE FOR EACH SUBDOMAIN (etype0)
        OPEN (34,file='metype0',form='unformatted',status='unknown')
C...... OPEN ELEM0 FILE
        OPEN (35,FILE='melem0',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN subnodes0 file for contact boundary of each subdomain
        OPEN (36,file='subnodes0',form='unformatted',status='unknown')
        open(123,file='u1v1',form='unformatted',status='unknown')
        open(124,file='fstr1',form='unformatted',status='unknown')
C
C
C
C.......WRITE OR SEND THE INITIAL DATA FOR EACH SUBDOMAIN
C...... OPEN LMDDM file for contact boundary of each subdomain
C        OPEN (40,file=' ',form='unformatted',status='unknown')
C...... OPEN COOR0 FILE
C        OPEN (41,FILE=' ',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ID0 FILE
C        OPEN (42,FILE=' ',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP0 FILE (Boundary condition file)
C        OPEN (43,FILE=' ',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPNE DISP1 FILE (Initial value file for displacement)
C        OPEN (44,FILE=' ',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ELEMENT TYPE FILE FOR EACH SUBDOMAIN (etype0)
C        OPEN (45,file=' ',form='unformatted',status='unknown')
C...... OPEN ELEM0 FILE
C        OPEN (46,FILE=' ',FORM='UNFORMATTED',STATUS='UNKNOWN')
C...... OPEN ssubnodes0 file for contact boundary of each subdomain
C        OPEN (47,file=' ',form='unformatted',status='unknown')
C
C
C     open data file for other initial value problem
C
C...... OPNE DISP2 FILE (Initial value file for velocity)
      if (idisp1.eq.1) then
      open(37,file='mdisp2',form='unformatted',status='unknown')
C      open(48,file=' ',form='unformatted',status='unknown')
      endif
C...... OPNE DISP3 FILE (Initial value file for acceleration)
      if (idisp2.eq.1) then
      open(38,file='mdisp3',form='unformatted',status='unknown')
C      open(49,file=' ',form='unformatted',status='unknown')
      endif
C
      DO 1000 iblk=0,NUMBLK-1
      READ (30) NUMNOD,NCOOR
      READ (31) KNODE,KDGOF
      READ (34) NUMTYP,NODALL,MATALL
      READ (36) KNODEALL,KDGOFALL
C
      KVAR=KNODE*KDGOF
      KVARALL=KNODEALL*KDGOFALL
      kelem=NODALL
      KVARALL = MAX(kelem*6,KVARALL)
      knodeall0 = knodeall/numblk*3
c      WRITE(*,*) 'KNODE,KDGOF,KVAR =',iblk
c      WRITE(*,'(1X,4I7)') KNODE,KDGOF,KVAR
C
C
      knb1=kdgof*knode*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      kna4=ncoor*knode*1
      kna1=kdgof*knode*1
      kna2=kdgof*knode*1
      kna3=kdgof*knode*1
      kna5=kdgof*knode*1
      knb8=knode*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      knb4=numtyp*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb5=numtyp*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      knb6=numtyp*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb7=numtyp*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb3=numtyp*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb2=kelem*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      kna6=matall*1
      knb9=1*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      knb10=1*1
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      knb11=1*1
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      knb12=1*1
      if (knb12/2*2 .lt. knb12) knb12=knb12+1
      knb13=kdgofall*knodeall0*1
      if (knb13/2*2 .lt. knb13) knb13=knb13+1
      knb14=kvarall*1
      if (knb14/2*2 .lt. knb14) knb14=knb14+1
      knb15=kelem/6
      if (knb15/2*2 .lt. knb15) knb15=knb15+1
      kna7=kvarall*1
	
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      kna7=kna7+kna6
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
      knb12=knb12+knb11
      knb13=knb13+knb12
      knb14=knb14+knb13
      knb15=knb15+knb14
      call azsendpart(knodeall0,knode,kdgof,kcoor,iblk,nsource,num,
     *nnode,numnod,ncoor,nparts,numblk,kelem,filflgdisp,numtyp,
     *matall,nodall,ncoor2,idisp1,idisp2,knodeall,kdgofall,kvarall,
     *maxlm,nmdof,numtypl,kdgofl,keymt,lgio,t0,tmax,
     *dt,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),ia(knb0),ia(knb1),
     *ia(knb2),ia(knb3),ia(knb4),ia(knb5),ia(knb6),
     *ia(knb7),ia(knb8),ia(knb9),ia(knb10),ia(knb11),
     *ia(knb12),ia(knb13),
     *filename,ia(knb14),istop)
 
1000    CONTINUE
      close(30)
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
      close(36)
      close(37)
      close(38)
      close(123)
      close(124)
      end
      subroutine azsendpart(knodeall0,knode,kdgof,kcoor,iblk,nsource,
     *num,
     *nnode,numnod,ncoor,nparts,numblk,kelem,filflgdisp,numtyp,
     *matall,nodall,ncoor2,idisp1,idisp2,knodeall,kdgofall,kvarall,
     *maxlm,nmdof,numtypl,kdgofl,keymt,lgio,t0,tmax,
     *dt,u0,u1,u2,coor,bfu,emateall,rtemp,
     *nodvar,nodeall,idet,numa,nnea,mmta,nmta,inode,
     *idnode,inodeall,jnode,idnode1,nodvarall,itemp,
     *filename,ielem,istop)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
       DIMENSION  NODVAR(KDGOF,KNODE),COOR(NCOOR,KNODE),
     &  U0(KDGOF,KNODE),U1(KDGOF,KNODE),U2(KDGOF,KNODE),
     &  BFU(KDGOF,KNODE),inode(knode),
     &  numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),
     &  nodeall(kelem),emateall(matall),
     &  idnode(1),inodeall(1),jnode(1),
     &  idnode1(1),nodvarall(kdgofall,knodeall0),
     &  itemp(kvarall),rtemp(kvarall)
       dimension ielem(kelem/6)
       logical filflgdisp(20)
c 
       kcoor=ncoor
       NODDOF = KDGOF
       if(nparts.ne.numblk) then
       write(*,*)'Error when get partition information of the dimain'
       stop 1111
       end if
C
       if(numnod.ne.knode) then
       write(*,*)'Error when get partition information of the dimain'
       stop 1111
       end if
C
C...  READ etype0 file
      read(34) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
     &nskip,(inode(i),i=1,knode)
      read(34) nskip,(ielem(i),i=1,numa(1))
C
C ... READ elem0 file
      read(35) (nodeall(i),i=1,nodall)
      read(35) (emateall(i),i=1,matall)
C
      goto 1212
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
1212  continue
C
C
C.......OPEN AND READ COOR file
        READ (30) ((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
c        WRITE(*,*) 'NUMNOD,NCOOR=',NUMNOD,NCOOR
c        WRITE(*,*)
c        DO I = 1, NUMNOD
c        WRITE(*,1100) I,(COOR(J,I),J=1,NCOOR)
c        END DO
C.......OPEN AND READ ID FILE
      READ (31) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C        WRITE(*,*) 'KNODE,KDGOF=',NUMNOD,KDGOF
C        WRITE(*,*)
C        DO I = 1, KNODE
C        WRITE(*,1000) I,(NODVAR(J,I),J=1,KDGOF)
C        END DO
C
C.......OPEN AND READ DISP0 FILE
        READ (32) NSKP,NSKP
      READ (32) ((BFU(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1000) I,(BFU(J,I),J=1,NODDOF)
C        END DO
C
C.......OPEN AND READ DISP1 FILE
        READ (33) NUMNOD,NODDOF
      READ (33) ((U0(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1100) I,(U0(J,I),J=1,NODDOF)
C        END DO
C
 
C........read index of subdomain nodes file
         read(36) nodesiblk,(idnode(i),i=1,1)
         read(36) nod_iblk,(inodeall(i),i = 1,1)
         read(36) knodeall,(jnode(i),i=1,1)
         read(36) knodeall,(idnode1(i),i=1,1)
         read(36) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
     &                                ,i=1,nod_iblk)
       NMUL = MAXLM/NODDOF
c         if(iblk.eq.8) then
c         WRITE(*,*) ' subnodes '
c         write(*,*) knodeall,kdgofall
c         write(*,*) nodesiblk,(idnode(i),i=1,1)
c         write(*,*) nod_iblk,(inodeall(i),i = 1,1)
c         write(*,*) knodeall,(jnode(i),i=1,1)
c         write(*,*) knodeall,(idnode1(i),i=1,1)
c         write(*,*) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c     &                                ,i=1,200)
cc     &                                ,i=1,nod_iblk)
c         open(56,file='aa',status='unknown',form='formatted')
c         write(56,*) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c     &                                ,i=1,nod_iblk)
c         close(56)
c         print *,'=============nod_iblk,knodeall0==',nod_iblk,knodeall0
c         end if
C
C.......OPEN AND READ DISP2 FILE
      if (idisp1.eq.1) then
        READ (37) NUMNOD,NODDOF
      READ (37)((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1100) I,(U1(J,I),J=1,NODDOF)
C        END DO
C
        endif
C.......OPEN AND READ DISP3 FILE
      if (idisp2.eq.1) then
        READ (38) NUMNOD,NODDOF
      READ (38) ((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
C
C        WRITE(*,*) 'NUMNOD,KDGOF=',NUMNOD,NODDOF
C        WRITE(*,*)
C        DO I = 1, NUMNOD
C        WRITE(*,1100) I,(U2(J,I),J=1,NODDOF)
C        END DO
C
      endif

C.......OPEN AND READ u1v1,fstr1 FILE
        if (istop.lt.0) then
          if (iblk.eq.0) then
          read(124) (rtemp(i),i=1,NUMNOD*6)
          read(124) (rtemp(i),i=1,numa(1)*6*8)
          read(123) (rtemp(i),i=1,NODDOF*NUMNOD)
          read(123) (rtemp(i),i=1,NODDOF*NUMNOD)
          endif
        endif

C
C       PLEASE UNCOMMAND THE LINE BELOW IF YOU JUST RUN THIS PROGRAM SEQUENTIALLY
C     GOTO 9988
C
C       WRITTE OUT ALL THE DATA FOR EACH SUBDOMAIN FROM THE MASTER NODE
C
        if(iblk.eq.0) return
C
C     WRITE LMDDM FILE FOR EACH SUBDOMAIN
C
C      WRITE(40)numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt,
C     & idisp1,idisp2
C      WRITE(40) lgio,t0,tmax,dt
      call sendint(iblk,nsource,numblk)
      call sendint(iblk,nsource,numtyp)
      call sendint(iblk,nsource,maxlm)
      call sendint(iblk,nsource,nmdof)
      call sendint(iblk,nsource,numtypl)
      call sendint(iblk,nsource,kdgofl)
      call sendint(iblk,nsource,keymt)
      call sendint(iblk,nsource,idisp1)
      call sendint(iblk,nsource,idisp2)
      call sendint(iblk,nsource,lgio)
      call sendr(iblk,nsource,t0)
      call sendr(iblk,nsource,tmax)
      call sendr(iblk,nsource,dt)
C
C
C       WRITE OUT THE COOR0 FILE
C
C        WRITE(41) NUMNOD,NCOOR
C     WRITE(41) ((COOR(I,J),I=1,NCOOR),J=1,NUMNOD)
        call sendint(iblk,nsource,numnod)
        call sendint(iblk,nsource,ncoor)
        iur = 0
        iur0=iur+1
        do J=1,NUMNOD
        do I=1,NCOOR
        iur = iur+1
        rtemp(iur)=coor(I,J)
        enddo
        enddo
        call sendar(iblk,nsource,rtemp(iur0),NUMNOD*NCOOR)
C
C       WRITE OUT THE ID0 FILE
C
c        WRITE(42) KNODE,KDGOF
c     WRITE(42) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
c
        call sendint(iblk,nsource,knode)
        call sendint(iblk,nsource,kdgof)
        iui = 0
        iui0=iui+1
        do J=1,KNODE
        do I=1,KDGOF
        iui = iui+1
        itemp(iui)=nodvar(I,J)
        enddo
        enddo
        call sendai(iblk,nsource,itemp(iui0),KNODE*KDGOF)
 
 
C
C       WRITE OUT THE BFU FILE
C
c        WRITE(43) KNODE,KDGOF
c     WRITE(43) ((BFU(I,J),I=1,KDGOF),J=1,KNODE)
        call sendint(iblk,nsource,knode)
        call sendint(iblk,nsource,kdgof)
        iur = 0
        iur0=iur+1
        do J=1,KNODE
        do I=1,KDGOF
        iur = iur+1
        rtemp(iur)=bfu(I,J)
        enddo
        enddo
        call sendar(iblk,nsource,rtemp(iur0),KNODE*KDGOF)
C
C       WRITE OUT THE DISP0 FILE
C
c        WRITE(44) KNODE,KDGOF
c        WRITE(44) ((U0(I,J),I=1,KDGOF),J=1,KNODE)
        call sendint(iblk,nsource,knode)
        call sendint(iblk,nsource,kdgof)
        iur = 0
        iur0=iur+1
        do J=1,KNODE
        do I=1,KDGOF
        iur = iur+1
        rtemp(iur)=u0(I,J)
        enddo
        enddo
        call sendar(iblk,nsource,rtemp(iur0),KNODE*KDGOF)
C
C
C       WRITE OUT THE ETYPE0 FILE
C
c        WRITE(45) NUMTYP,NODALL,MATALL
c        WRITE(45) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
c     & (inode(i),i=1,knode)
        call sendint(iblk,nsource,numtyp)
        call sendint(iblk,nsource,nodall)
        call sendint(iblk,nsource,matall)
        iui = 0
        iur = 0
        iui0=iui+1
        iur0=iur+1
        do i=1,numtyp
        iui = iui+1
        itemp(iui)=idet(i)
        iui = iui+1
        itemp(iui)=numa(i)
        iui = iui+1
        itemp(iui)=nnea(i)
        iui = iui+1
        itemp(iui)=mmta(i)
        iui = iui+1
        itemp(iui)=nmta(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),5*numtyp)
        iui0=iui+1
        iur0=iur+1
        do i=1,knode
        iui = iui+1
        itemp(iui)=inode(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*knode)
        iui0=iui+1
        iur0=iur+1
        do i=1,numa(1)
        iui = iui+1
        itemp(iui)=ielem(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),numa(1))
C
C       WRITE OUT THE ELEM0 FILE
C
c        WRITE(46) (nodeall(i),i=1,nodall)
c        WRITE(46) (emateall(i),i=1,matall)
       iui = 0
       iur = 0
       iui0=iui+1
       iur0=iur+1
       do i=1,nodall
       iui = iui+1
       itemp(iui)=nodeall(i)
       enddo
       call sendai(iblk,nsource,itemp(iui0),1*nodall)
       iui = 0
       iur = 0
       iui0=iui+1
       iur0=iur+1
       do i=1,matall
       iur = iur+1
       rtemp(iur)=emateall(i)
       enddo
       call sendar(iblk,nsource,rtemp(iur0),1*matall)
C........write out index of subdomain nodes file
c        write(47) knodeall,kdgofall
c        write(47) nodesiblk,(idnode(i),i=1,nodesiblk)
c        write(47) nod_iblk,(inodeall(i),i = 1,nod_iblk)
c        write(47) knodeall,(jnode(i),i=1,knodeall)
c        write(47) knodeall,(idnode1(i),i=1,knodeall)
c        write(47) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c     &                                ,i=1,knodeall)
        call sendint(iblk,nsource,knodeall)
        call sendint(iblk,nsource,kdgofall)
c       write(47) nodesiblk,(idnode(i),i=1,nodesiblk)
        iui = 0
        iur = 0
        iui = iui+1
        itemp(iui)=nodesiblk
        call sendint(iblk,nsource,itemp(iui))
        iui0=iui+1
        iur0=iur+1
        do i=1,1
        iui = iui+1
        itemp(iui)=idnode(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*1)
c       write(47) nod_iblk,(inodeall(i),i = 1,nod_iblk)
        iui = 0
        iur = 0
        iui = iui+1
        itemp(iui)=nod_iblk
        call sendint(iblk,nsource,itemp(iui))
        iui0=iui+1
        iur0=iur+1
        do i=1,1
        iui = iui+1
        itemp(iui)=inodeall(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*1)
c       write(47) knodeall,(jnode(i),i=1,knodeall)
        iui = 0
        iur = 0
        iui = iui+1
        itemp(iui)=knodeall
        call sendint(iblk,nsource,itemp(iui))
        iui0=iui+1
        iur0=iur+1
        do i=1,1
        iui = iui+1
        itemp(iui)=jnode(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*1)
c       write(47) knodeall,(idnode1(i),i=1,knodeall)
        iui = 0
        iur = 0
        iui = iui+1
        itemp(iui)=knodeall
        call sendint(iblk,nsource,itemp(iui))
        iui0=iui+1
        iur0=iur+1
        do i=1,1
        iui = iui+1
        itemp(iui)=idnode1(i)
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*1)
c       write(47) knodeall,kdgofall,((nodvarall(j,i),j=1,kdgofall)
c    &                                ,i=1,knodeall)
        iui = 0
        iur = 0
        iui = iui+1
        itemp(iui)=knodeall
        call sendint(iblk,nsource,itemp(iui))
        iui = iui+1
        itemp(iui)=kdgofall
        call sendint(iblk,nsource,itemp(iui))
        iui0=iui+1
        iur0=iur+1
        do i=1,nod_iblk
        do j=1,kdgofall
        iui = iui+1
        itemp(iui)=nodvarall(j,i)
        enddo
        enddo
        call sendai(iblk,nsource,itemp(iui0),1*nod_iblk*kdgofall)
C
C       WRITE OUT DISP2 FOR INITIAL DU/DT
C
      if (idisp1.eq.1) then
c        WRITE(48)NUMNOD,NODDOF
c     WRITE(48)((U1(I,J),I=1,NODDOF),J=1,NUMNOD)
        call sendint(iblk,nsource,numnod)
        call sendint(iblk,nsource,noddof)
        iur = 0
        iur0=iur+1
        do J=1,KNODE
        do I=1,KDGOF
        iur = iur+1
        rtemp(iur)=u1(I,J)
        enddo
        enddo
        call sendar(iblk,nsource,rtemp(iur0),KNODE*KDGOF)
        endif
C
C
C       WRITE OUT DISP3 FOR INITIAL (DU/DT)2
C
      if (idisp2.eq.1) then
c        WRITE(49)NUMNOD,NODDOF
c     WRITE(49)((U2(I,J),I=1,NODDOF),J=1,NUMNOD)
        call sendint(iblk,nsource,numnod)
        call sendint(iblk,nsource,noddof)
        iur = 0
        iur0=iur+1
        do J=1,KNODE
        do I=1,KDGOF
        iur = iur+1
        rtemp(iur)=u2(I,J)
        enddo
        enddo
        call sendar(iblk,nsource,rtemp(iur0),KNODE*KDGOF)
        endif

        if (istop.lt.0) then
          read(124) (rtemp(i),i=1,6*NUMNOD)
          call sendar(iblk,nsource,rtemp(1),NUMNOD*6)
          read(124) (rtemp(i),i=1,numa(1)*6*8)
          call sendar(iblk,nsource,rtemp(1),numa(1)*6*8)
          read(123) (rtemp(i),i=1,NODDOF*NUMNOD)
          call sendar(iblk,nsource,rtemp(1),NODDOF*NUMNOD)
          read(123) (rtemp(i),i=1,NODDOF*NUMNOD)
          call sendar(iblk,nsource,rtemp(1),NODDOF*NUMNOD)
        endif

C
1000    format(10i6)
1100    format(i6,8e16.5)
1200    format(6i6)
      return
      end
 
 
