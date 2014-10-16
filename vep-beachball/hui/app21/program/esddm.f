      subroutine Mesddm(iblk,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
6       format(2x,10i5)
7       format(2x,6e12.5)
        kelem0=400000
        nsource=0
c...... open lmddm file
c       open(21,file='lmddm',form='unformatted',status='unknown')
      ISTATUS = 1
      call openf(21,1,ISTATUS)
c       read(21) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt
      call initrwf(21,iui,iur)
      iui = iui+1
      numblk=ipool(iui)
      iui = iui+1
      numtyp=ipool(iui)
      iui = iui+1
      maxlm=ipool(iui)
      iui = iui+1
      nmdof=ipool(iui)
      iui = iui+1
      numtypl=ipool(iui)
      iui = iui+1
      kdgofl=ipool(iui)
      iui = iui+1
      keymt=ipool(iui)
      call endrwf(21,iui,iur)
c       read(21) lgio,t0,tmax,dt
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
c        write(*,*) 'numblk,numtyp =======',numblk,numtyp
c...... open time file
c       open(21,file='time',form='unformatted',status='unknown')
      ISTATUS = 11
      call openf(21,11,ISTATUS)
c       read(21) tmax,dt,time,it
      call initrwf(21,iui,iur)
      iur = iur+1
      tmax=rpool(iur)
      iur = iur+1
      dt=rpool(iur)
      iur = iur+1
      time=rpool(iur)
      iui = iui+1
      it=ipool(iui)
      call endrwf(21,iui,iur)
c       close(21)
c        write(*,*) 'tmax,dt,time,it =',tmax,dt,time,it
c       hzhang modified here for there is no need of the files: ids0 and disps0
c...... open disp0 file
c        open(21,file=' ',form='unformatted',status='unknown')
c...... open coor0 file
c       open(22,file='scoor0',form='unformatted',status='unknown')
c       hzhang modified here for there is no need of the files: ids0 and disps0
c...... open nv file
c        open(23,file=' ',form='unformatted',status='unknown')
c...... open etype0 file
c       open(24,file='setype0',form='unformatted',status='unknown')
c...... open elem0 file
c       open(25,file='selem0',form='unformatted',status='unknown')
c...... open stress file
        if(iblk.eq.0) then
1        open(26,file='munods0',form='unformatted',status='unknown')
        else
c       open(26,file='munods0',form='unformatted',status='unknown')
        end if
cc...... write out the lmddm file for the recvstr file usage
c        open(28,file=' ',form='unformatted',status='unknown')
c...... open unod file
c     open(11,file='unod',form='unformatted',status='old')
        kdgof = 6
c        keymt = kdgof
c
c        DO 1000 iblk=0,NUMBLK-1
c       hzhang modified here for there is no need of the files: ids0 and disps0
c        if(iblk .eq. 1) then
c        write(28) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt
c        write(28) lgio,t0,tmax,dt
c        end if
      ISTATUS = 2
      call openf(22,2,ISTATUS)
c       READ(22) knode,ncoor
      call initrwf(22,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      ncoor=ipool(iui)
      call endrwf(22,iui,iur)
        KVAR=KNODE*KDGOF
c        READ(22) NSKP,NCOOR
c        WRITE(*,*) 'KNODE,KDGOF,KVAR,NCOOR ===='
c        WRITE(*,6)  KNODE,KDGOF,KVAR,NCOOR
C.....  READ ETYPE00 FILE
      ISTATUS = 6
      call openf(24,6,ISTATUS)
c       read(24) numtyp,numnod,nummat
      call initrwf(24,iui,iur)
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
      call endrwf(24,iui,iur)

        kelem = numnod*kdgof

      knb5=kdgof*knode*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      kna1=kdgof*knode*1
      kna3=ncoor*knode*1
      kna6=kdgof*knode*1
      kna7=kdgof*knode*1
      kna8=knode*1
      kna9=knode*1
      kna10=knode*1
      kna11=kvar*1
      kna2=kvar*1
      kna4=150000*1
      kna12=kelem0
      kna13=knode*3
      kna14=knode*3
      kna15=kelem*1
      kna16=numnod*2
      kna17=knode*6
      knb6=kelem*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb7=numtyp*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb8=numnod*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      knb9=numnod*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      kna5=nummat*1
      knb1=numtyp*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=numtyp*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numtyp*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb4=numtyp*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
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
      kna12=kna12+kna11
      kna13=kna13+kna12
      kna14=kna14+kna13
      kna15=kna15+kna14
      kna16=kna16+kna15
      kna17=kna17+kna16
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
      call esddm(knode,kdgof,kvar,ncoor,numtyp,numel,
     *neq,kelem,time,dt,it,numnod,nummat,nsource,
     *iblk,aa(kna0),aa(kna1),aa(kna2),aa(kna3),
     *aa(kna4),aa(kna5),aa(kna6),aa(kna7),aa(kna8),
     *aa(kna9),aa(kna10),aa(kna11),ia(knb0),ia(knb1),
     *ia(knb2),ia(knb3),ia(knb4),ia(knb5),ia(knb6),
     *ia(knb7),
     *filename,aa(kna12),aa(kna13),aa(kna14),aa(kna15),
     *ia(knb8),istop,aa(kna16),kend)
 
1000    CONTINUE
c       CLOSE(21)
c       CLOSE(22)
c        CLOSE(23)
c       CLOSE(24)
c       CLOSE(25)
        if(iblk.eq.0) then
2        CLOSE(26)
        else
c       CLOSE(26)
        end if
c        close(28)
      end
      subroutine esddm(knode,kdgof,kvar,ncoor,numtyp,numel,
     *neq,kelem,time,dt,it,numnod,nummat,nsource,
     *iblk,u,f,coor,emate,emata,estr,ezmain,
     *eu,ev,ew,emass,sml,numa,nnea,mmta,
     *nmta,nodvar,node,idet,nodea,
     *filename,eu1,ev1,ess,enrg,ilabel,istop,smain,kend)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION NODVAR(KDGOF,KNODE),U(KDGOF,KNODE),COOR(NCOOR,KNODE),
     *estr(kdgof,knode),ezmain(kdgof,knode),eu(knode),ev(knode),
     *ew(knode),Emass(kvar),
     *  F(KVAR),EMATE(150000),SML(KELEM),NODE(KELEM),IDET(NUMTYP),
     *  nodea(numnod),emata(nummat),numa(numtyp),nnea(numtyp),
     *  mmta(numtyp),nmta(numtyp)
        dimension eu1(3,knode),ev1(3,knode),ess(numnod*6)
        dimension enrg(numnod*2),ilabel(numnod),smain(6,knode)

6       FORMAT(2X,10I5)
7       FORMAT(2X,6E12.5)
c       hzhang modified here for there is no need of the files: ids0 and disps0
cC...... READ DISPS0 FILE
c        READ(21) ((U(J,I),J=1,KDGOF),I=1,KNODE)
        do i = 1,knode
      do j = 1,kdgof
      u(j,i) = 0.d0
      end do
      end do
c       hzhang modified here for there is no need of the files: ids0 and disps0
cC...... READ ID0 FILE
c        READ(23) nskp,nskp
c        READ(23) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
        do i = 1,knode
      do j = 1,kdgof
      nodvar(j,i) = 1
      end do
      end do
C...... COMPUTE NODVAR
        NEQ = 0
        DO 100 J=1,KNODE
        DO 100 I=1,KDGOF
          IF (NODVAR(I,J).LE.0) GOTO 100
          NEQ = NEQ+1
          NODVAR(I,J) = NEQ
100     CONTINUE
        DO 200 J=1,KNODE
        DO 200 I=1,KDGOF
          IF (NODVAR(I,J).GT.-1) GOTO 200
          N=-NODVAR(I,J)-1
127        CONTINUE
          IF (NODVAR(I,N).LT.-1) THEN
          N=-NODVAR(I,N)-1
          GOTO 127
          ENDIF
          NODVAR(I,J) = NODVAR(I,N)
200     CONTINUE
CC        WRITE(*,*) 'KNODE,KDGOF,NODVAR ==='
CC        WRITE (*,6) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
C...... READ COOR0 FILE
c       READ(22) ((COOR(I,J),I=1,NCOOR),J=1,KNODE)
      call initrwf(22,iui,iur)
      do J=1,KNODE
      do I=1,NCOOR
      iur = iur+1
      coor(I,J)=rpool(iur)
      enddo
      enddo
      call endrwf(22,iui,iur)
      ISTATUS = 18
      call openf(11,18,ISTATUS)
c     read(11) (eu(i),i=1,knode),
c    &  (ev(i),i=1,knode),
c    &  (ew(i),i=1,knode)
      call initrwf(11,iui,iur)
      do i=1,knode
      iur = iur+1
      eu(i)=rpool(iur)
      enddo
      do i=1,knode
      iur = iur+1
      ev(i)=rpool(iur)
      enddo
      do i=1,knode
      iur = iur+1
      ew(i)=rpool(iur)
      enddo
      call endrwf(11,iui,iur)

C
C     principle stess at nodes at previous step
C

      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,6
      do I=1,KNODE
      iur = iur+1
      smain(j,i)=rpool(iur)
      enddo
      enddo
      call endrwf(26,iui,iur)

C
C     str6.tmp for Stress at previous time step
C
      ISTATUS = 32
      call openf(42,32,ISTATUS)
      call initrwf(42,iui,iur)
      do i = 1,NUMNOD*6
      iur = iur+1
      ess(i) = rpool(iur)
      end do
      call endrwf(42,iui,iur)

C
C     unod for displacement and velocity at previous time step
C
      NUMUNIT = 43
      ISTATUS = 33
      call openf(43,33,ISTATUS)
      call initrwf(43,iui,iur)
      do j = 1,3
      do i = 1,knode
      iur = iur+1
      eu1(j,i) = rpool(iur)
      end do
      end do
c      
      do j = 1,3
      do i = 1,knode
      iur = iur+1
      ev1(j,i) = rpool(iur)
      end do
      end do
c      
      call endrwf(43,iui,iur)

C
C     label for failure switch
C
      ISTATUS = 37
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      do i = 1,NUMNOD
      iui = iui+1
      ilabel(i) = ipool(iui)
      end do
      call endrwf(47,iui,iur)

C
C     energy release accumulation
C
      ISTATUS = 39
      call openf(49,39,ISTATUS)
      call initrwf(49,iui,iur)
      do i = 1,NUMNOD*2
      iur = iur+1
      enrg(i) = rpool(iur)
      end do
      call endrwf(49,iui,iur)


        DO I=1,NEQ
        Emass(i) = 0.0d0
        ENDDO
C...... READ ETYPE0 FILE
C        READ(24) nskp,nskp
C        READ(24) (IDET(I),I=1,NUMTYP)
c       read(24) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp)
      call initrwf(24,iui,iur)
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
      call endrwf(24,iui,iur)
c        WRITE(*,*) 'IDET ============='
c        WRITE(*,6) (IDET(I),I=1,NUMTYP)
C...... READ ELEM0 FILE
      ISTATUS = 7
      call openf(25,7,ISTATUS)
c       read(25) (nodea(i),i=1,numnod)
      call initrwf(25,iui,iur)
      do i=1,numnod
      iui = iui+1
      nodea(i)=ipool(iui)
      enddo
      call endrwf(25,iui,iur)
        if(nummat.gt.0) then
c       read(25) (emata(i),i=1,nummat)
      call initrwf(25,iui,iur)
      do i=1,nummat
      iur = iur+1
      emata(i)=rpool(iur)
      enddo
      call endrwf(25,iui,iur)
        endif
 
        NUMEL=0
 
        nodadd=0
        nmtadd=0

c    ......... NUMTYP in esddm is different from that in eddm

        DO 300 ITYP=1,1
        IF (IDET(ITYP).LE.0) GOTO 300
C...... READ ELEM FILE
C        READ(25) NUM,NNODE
C        READ(25) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        num=numa(ityp)
        nnode=nnea(ityp)
        do i=1,num
           do j=1,nnode
              node((i-1)*nnode+j)=nodea(nodadd+(i-1)*nnode+j)
           enddo
        enddo
        nodadd=nodadd+num*nnode
 
CC        WRITE(*,*) 'NUM,NNODE,NODE ===',NUM,NNODE
CC        WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        NNE = NNODE
      nne = nne-1
        IF (NNE.EQ.NNODE) THEN
          IF (NCOOR.EQ.2.AND.NNE.EQ.2) GOTO 300
          IF (NCOOR.EQ.3.AND.(NNE.EQ.3.OR.NNE.EQ.4)) GOTO 300
        ENDIF
        IF (NNE.LT.NNODE) THEN
          IF (NCOOR.EQ.2.AND.NNE.EQ.2) THEN
C            READ(25) nmat,mmat
C            READ(25) (xskp,i=1,nmat*mmat)
            nmat=nmta(ityp)
            mmat=mmta(ityp)
            do i=1,nmat*mmat
               xskp=emata(nmtadd+i)
            enddo
C            nmtadd=nmtadd+nmat*mmat
            GOTO 300
          ENDIF
          IF (NCOOR.EQ.3.AND.(NNE.EQ.3.OR.NNE.EQ.4)) THEN
C            READ(25) nmat,mmat
C            READ(25) (xskp,i=1,nmat*mmat)
            nmat=nmta(ityp)
            mmat=mmta(ityp)
            do i=1,nmat*mmat
               xskp=emata(nmtadd+i)
            enddo
C            nmtadd=nmtadd+nmat*mmat
            GOTO 300
          ENDIF
        ENDIF
        K=0
        DO 115 J=1,NNE
        JNOD = NODE(J)
        DO 115 L=1,KDGOF
        IF (NODVAR(L,JNOD).NE.0) K=K+1
115     CONTINUE
C        WRITE(*,*) 'K =',K
      kk=k*k
      k0=1
      k1=k0+k*k
      k2=k1+k
      k3=k2+k
      k4=k3+k*k
      k5=k4+k*k
      k6=k5+k
      k7=k6+k
        if (k7.ge.400000) then
        write(*,*) 'it is error: k7 gt kelem0'
        stop 0000 
        endif     

        CALL EsTSUB(KNODE,KDGOF,kvar,IT,KELEM,K,KK,NNODE,NNE,
     *  NUMEL,ITYP,NCOOR,NUM,TIME,DT,NODVAR,COOR,NODE,EMATE,
     *  nmtadd,emata,mmta,nmta,
     *sml(k0),sml(k1),sml(k2),sml(k3),sml(k4),sml(k5),sml(k6),
     &eu,ev,ew,estr,ezmain,Emass,
     *U,eu1,ev1,ess,enrg,ilabel,smain)
        NUMEL = NUMEL + NUM
300     CONTINUE
        DO II=1,NEQ
          F(II)=0.0d0
        ENDDO
c        WRITE(*,*) 'EMASS ====='
c        WRITE(*,7) (EMASS(I),I=1,NEQ)
        DO 400 I=1,KNODE
        DO 400 J=1,KDGOF
          IJ=NODVAR(J,I)
          IF (IJ.LE.0) GOTO 400
          F(IJ) = F(IJ)+U(J,I)/EMASS(IJ)
400     CONTINUE
        DO 500 I=1,KNODE
          DO 600 J=1,KDGOF
          IJ=NODVAR(J,I)
          estr(J,I) = 0.0d0
          IF (IJ.LE.0) GOTO 600
          estr(J,I) = F(IJ)
600       CONTINUE
500     CONTINUE
c        WRITE(*,*) 'estr(J,I),KNODE,KDGOF ===',KNODE,KDGOF
c        DO I=1,KNODE
c          WRITE(*,'(2X,I5,6E12.5)') I,(estr(J,I),J=1,KDGOF)
c        ENDDO
C...... WRITE STRESS FILE
        if(iblk.eq.0) then
c1        WRITE(26) knode,kdgof
c2        WRITE(26) ((estr(J,I),I=1,KNODE),J=1,KDGOF)
      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui = iui+1
      ipool(iui)=knode
      iui = iui+1
      ipool(iui)=kdgof
      call endrwf(26,iui,iur)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,KDGOF
      do I=1,KNODE
      iur = iur+1
      rpool(iur)=estr(J,I)
      enddo
      enddo
      call endrwf(26,iui,iur)
             if (istop.ge.1.and.kend.eq.1) then
             WRITE(26) (ess(i),i=1,numa(1)*6*8)
             endif
        else
      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui = iui+1
      ipool(iui)=knode
      if (kend.eq.1) then
        call sendint(nsource,iblk,ipool(iui))
      endif
      iui = iui+1
      ipool(iui)=kdgof
      if (kend.eq.1) then
        call sendint(nsource,iblk,ipool(iui))
      endif
      call endrwf(26,iui,iur)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,KDGOF
      do I=1,KNODE
      iur = iur+1
      rpool(iur)=estr(J,I)
      enddo
      enddo
      if (kend.eq.1) then
        call sendar(nsource,iblk,rpool(iur0),1*KDGOF*KNODE)
      endif
      call endrwf(26,iui,iur)
             if (istop.ge.1.and.kend.eq.1) then
             call sendar(nsource,iblk,ess(1),1*numa(1)*6*8)
             endif
        end if

      if (kend.eq.1) then

c
c     renew ess file
c
      ISTATUS = 32
      call openf(42,32,ISTATUS)
      call initrwf(42,iui,iur)
      do i = 1,NUMNOD*6
      iur = iur+1
      rpool(iur) = ess(i)
      end do
      call endrwf(42,iui,iur)

c
c     renew label file
c
      ISTATUS = 37
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      do i = 1,NUMNOD
      iui = iui+1
      ipool(iui) = 0
      end do
      call endrwf(47,iui,iur)

c
c     renew energy release accumulation file
c
      ISTATUS = 39
      call openf(49,39,ISTATUS)
      call initrwf(49,iui,iur)
      do i = 1,NUMNOD*2
      iur = iur+1
      rpool(iur) = enrg(i)
      end do
      call endrwf(49,iui,iur)

      endif

        RETURN
      END
 
 
        SUBROUTINE EsTSUB(KNODE,KDGOF,kvar,IT,KELEM,K,KK,NNODE,NNE,
     *  NUMEL,ITYP,NCOOR,NUM,TIME,DT,NODVAR,COOR,NODE,EMATE,
     *  nmadd,emata,mmta,nmta,
     *es,em,ef,Estifn,Estifv,Emassn,Emassv,
     *eu,ev,ew,estr,ezmain,Emass,
     *U,eu1,ev1,ess,enrg,ilabel,smain)
        IMPLICIT REAL*8 (A-H,O-Z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION NODVAR(KDGOF,KNODE),COOR(NCOOR,KNODE),NODE(KELEM),
     *  U(KDGOF,KNODE),EMATE(*),emata(*),mmta(*),nmta(*),
     *es(k,k),em(k),ef(k),eu(knode),
     *ev(knode),ew(knode),estr(kdgof,knode),
     *ezmain(kdgof,knode),Estifn(k,k),Estifv(kk),Emassn(k),Emassv(k),
     *Emass(kvar),
     *  R(3000),PRMT(3000),COEF(3000),LM(3000)
        dimension eu1(3,knode),ev1(3,knode),ess(KELEM),fstr(6,8)
        dimension enrg(KELEM/3+1),erg(2,8),ilabel(kelem/6+1)
        dimension label(8),smain(6,knode)
17      FORMAT (1X,10I5)
18      FORMAT (1X,6E12.5)
8       FORMAT (1X,6E12.5)
c        READ (25) MMATE,NMATE
c        READ (25) ((EMATE((I-1)*NMATE+J),J=1,NMATE),I=1,MMATE)
      MMATE=MMTA(ITYP)
      NMATE=NMTA(ITYP)
      DO I=1,MMATE
         DO J=1,NMATE
            EMATE((I-1)*NMATE+J)=EMATA(NMADD+(I-1)*NMATE+J)
         ENDDO
      ENDDO
      NMADD=NMADD+MMATE*NMATE
CC        WRITE(*,*) 'MMATE =',MMATE,' NMATE =',NMATE
CC        WRITE (*,*) 'EMATE ='
CC        WRITE (*,8) ((EMATE((I-1)*NMATE+J),J=1,NMATE),
CC     *  I=1,MMATE)
      DO 1000 NE=1,NUM
      NR=0
      DO 130 J=1,NNE
      JNOD = NODE((NE-1)*NNODE+J)
      IF (JNOD.LT.0) JNOD = -JNOD
      coef(j+0*nne)=eu(jnod)
      coef(j+1*nne)=ev(jnod)
      coef(j+2*nne)=ew(jnod)
      coef(j+3*nne)=(smain(1,jnod)+smain(2,jnod)+smain(3,jnod))/3
      coef(j+4*nne)=0.d0
      coef(j+5*nne)=0.d0
      coef(j+6*nne)=0.d0
      coef(j+7*nne)=0.d0
      coef(j+8*nne)=0.d0
      DO 120 I=1,NCOOR
      NR=NR+1
120   R(NR) = COOR(I,JNOD)
130   CONTINUE
      IMATE = NODE(NNODE*NE)
      DO 140 J=1,NMATE
140   PRMT(J) = EMATE((IMATE-1)*NMATE+J)
      PRMT(NMATE+1)=TIME
      PRMT(NMATE+2)=DT
      prmt(nmate+3)=0.d0+imate
      nnee=ne 

      goto 1
1     continue
      do j=1,8
       do i=1,6
       fstr(i,j)=ess((ne-1)*6*8+(j-1)*6+i)
       enddo
       label(j)=ilabel((ne-1)*8+j)
       do i=1,2
       erg(i,j)=enrg((ne-1)*2*8+(j-1)*2+i)
       enddo
      enddo
      call ssc8g2(r,coef,prmt,es,em,ef,nnee,fstr,label,erg)
      do j=1,8
       do i=1,6
       ess((ne-1)*6*8+(j-1)*6+i)=fstr(i,j)
       enddo
       do i=1,2
       enrg((ne-1)*2*8+(j-1)*2+i)=erg(i,j)
       enddo
      enddo

      goto 2
2     continue
C       WRITE(*,*) 'ES EM EF ='
C       DO 555 I=1,K
C555    WRITE(*,18) (ES(I,J),J=1,K)
C       WRITE(*,18) (EM(I),I=1,K)
C       WRITE(*,18) (EF(I),I=1,K)
 
CC      IF (IT.GT.0) THEN
      do 201 i=1,k
      do 201 j=1,k
      Estifn(i,j)=0.0d0
201   continue
      do 202 i=1,k
      Estifn(i,i)=Estifn(i,i)
      do 202 j=1,k
      Estifn(i,j)=Estifn(i,j)+es(i,j)
202   continue
      do 203 i=1,k
      Emassn(i)=0.0d0
203   continue
      do 204 i=1,k
      Emassn(i)=Emassn(i)+em(i)
204   continue
 
      L=0
      M=0
      I=0
      DO 700 INOD=1,NNE
      NODI=NODE((NE-1)*NNODE+INOD)
      DO 600 IDGF=1,KDGOF
      INV=NODVAR(IDGF,NODI)
      IF (INV.EQ.0) GOTO 600
      I=I+1
      IF (INV.LT.0) GOTO 305
      L=L+1
      LM(L)=INV
      U(IDGF,NODI)=U(IDGF,NODI)
     *+ef(i)
      Emassv(l)=Emassn(i)
305     J=0
      DO 500 JNOD=1,NNE
      NODJ=NODE((NE-1)*NNODE+JNOD)
      DO 400 JDGF=1,KDGOF
      JNV=NODVAR(JDGF,NODJ)
      IF (JNV.EQ.0) GOTO 400
      J=J+1
 
      IF (JNV.LT.0) GOTO 400
      IF (INV.LT.0) GOTO 310
      M=M+1
      Estifv(m)=Estifn(i,j)
310     CONTINUE
 
 
      IF (INV.LT.0)
     *  U(JDGF,NODJ)=U(JDGF,NODJ)-ESTIFN(J,I)*U(IDGF,NODI)
400     CONTINUE
500     CONTINUE
600     CONTINUE
700     CONTINUE
C       WRITE (*,*) 'U ='
C       WRITE (*,18) ((U(J,I),J=1,KDGOF),I=1,KNODE)
 
      LRD=M
      NER=NUMEL+NE
      DO 800 I=1,L
      J=LM(I)
      Emass(j) = Emass(j) + Emassv(i)
800     CONTINUE
1000    CONTINUE
        RETURN
        END
