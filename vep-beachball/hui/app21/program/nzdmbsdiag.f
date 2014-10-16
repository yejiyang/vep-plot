      subroutine Mnzdmbsdiag(iblk)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
6     format (1x,10i10)
7     format (1x,6e12.5)
 
1001  format(1x,7i7)
      nsource=0
      kelem=500000
      maxnd=100
c.....open lmddm file
c     open(21,file='lmddm',form='unformatted',status='old')
      ISTATUS = 1
      call openf(21,1,ISTATUS)
c     read(21) numblk,numtyp,nskp,nskp,nskp,nskp,nskp
      call initrwf(21,iui,iur)
      iui = iui+1
      numblk=ipool(iui)
      iui = iui+1
      numtyp=ipool(iui)
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
      call endrwf(21,iui,iur)
c     read(21) lgio
      call initrwf(21,iui,iur)
      iui = iui+1
      lgio=ipool(iui)
      call endrwf(21,iui,iur)
c     close(21)
c ....open etype0 file
c     open (30,file='setype0',form='unformatted',status='unknown')
c ....open sys file
c     open (31,file='sys',form='unformatted',status='unknown')
c.....open nodvar file
c     open (32,file='nv',form='unformatted',status='unknown')
c.....open elem0 file
c     open (33,file='selem0',form='unformatted',status='unknown')
c.....open disp0 file
c     open (34,file='sdisp0',form='unformatted',status='unknown')
c.....open diag file
c     open (35,file='jdiag',form='unformatted',status='unknown')
c.....open na file
c     open (36,file='na',form='unformatted',status='unknown')
c
c      do 2300 iblk=1,numblk
C....   READ ETYPE0 FILE
      ISTATUS = 6
      call openf(30,6,ISTATUS)
c       read(30) numtyp,numnod,nummat
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
C....   READ NV FILE
      ISTATUS = 12
      call openf(32,12,ISTATUS)
c       read(32) neq
      call initrwf(32,iui,iur)
      iui = iui+1
      neq=ipool(iui)
      call endrwf(32,iui,iur)
        neq1=neq+1
        maxna=neq*maxnd*3
c        maxna = (maxia/3)*2
C....   READ DISP0 FILE
      ISTATUS = 4
      call openf(34,4,ISTATUS)
c       read(34) knode,kdgof
      call initrwf(34,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      call endrwf(34,iui,iur)
c        read(34) ((xskp,j=1,kdgof),i=1,knode)
        kvar=knode*kdgof
c       write(*,*) 'knode,kdgof,kvar,neq,maxnd,maxna ==='
c       write(*,6)  knode,kdgof,kvar,neq,maxnd,maxna
 
      knb2=kdgof*knode*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=kelem*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      kna2=150000*1
      knb4=maxnd*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb6=neq1*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb7=neq1*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      knb8=neq1*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      knb1=numtyp*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb5=numnod*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      kna1=nummat*1
      knb13=maxna*1
      if (knb13/2*2 .lt. knb13) knb13=knb13+1
      knb9=numtyp*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      knb10=numtyp*1
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      knb11=numtyp*1
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      knb12=numtyp*1
      if (knb12/2*2 .lt. knb12) knb12=knb12+1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
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
      call nzdmbsdiag(knode,kdgof,kelem,numtyp,iblk,numblk,
     *kvar,numnod,nummat,maxnd,maxna,maxt,neq,neq1,
     *nsource,lgio,aa(kna0),aa(kna1),ia(knb0),
     *ia(knb1),ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *ia(knb6),ia(knb7),ia(knb8),ia(knb9),ia(knb10),
     *ia(knb11),ia(knb12),
     *filename)
 
2300  continue
 
c     close(21)
c     close(22)
c     close(25)
c     close(28)
c     close(29)
c     close(26)
c     close(30)
      end
      subroutine nzdmbsdiag(knode,kdgof,kelem,numtyp,iblk,numblk,
     *kvar,numnod,nummat,maxnd,maxna,maxt,neq,neq1,
     *nsource,lgio,emata,emate,idet,nodvar,node,lm,
     *nodea,jdiag,jdiagr,jdiagaz,numa,nnea,nmata,nmatpa,
     *na,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension nodvar(kdgof,knode),node(kelem),emate(150000),
     *lm(maxnd),jdiag(neq1),jdiagr(neq1),jdiagaz(neq1),idet(numtyp)
      dimension nodea(numnod),emata(nummat),na(maxna)
      dimension numa(numtyp),nnea(numtyp),nmata(numtyp),
     *          nmatpa(numtyp)
      CHARACTER*1 MATERIAL
 
6       format (1x,10i5)
7       format (1x,6e12.5)
          DO I=1,NEQ1
            JDIAG(I)=0
            JDIAGR(I)=0
            JDIAGAZ(I)=0
          ENDDO
          DO I=1,MAXNA
            NA(I)=0
          ENDDO
c
c...... READ NV FILE
c       READ(32) ((NODVAR(I,J),I=1,KDGOF),J=1,KNODE)
      call initrwf(32,iui,iur)
      do J=1,KNODE
      do I=1,KDGOF
      iui = iui+1
      nodvar(I,J)=ipool(iui)
      enddo
      enddo
      call endrwf(32,iui,iur)
c        if (lgio.eq.0) material='y'
c        if (lgio.eq.1) material='n'
        material='y'
C.......READ ETYPE0 FILE
c       READ(30) (IDET(I),numa(i),nnea(i),
c    &   nmata(i),nmatpa(i),I=1,NUMTYP)
      call initrwf(30,iui,iur)
      do I=1,NUMTYP
      iui = iui+1
      idet(I)=ipool(iui)
      iui = iui+1
      numa(i)=ipool(iui)
      iui = iui+1
      nnea(i)=ipool(iui)
      iui = iui+1
      nmata(i)=ipool(iui)
      iui = iui+1
      nmatpa(i)=ipool(iui)
      enddo
      call endrwf(30,iui,iur)
C.......READ ELEM0 FILE
      ISTATUS = 7
      call openf(33,7,ISTATUS)
c       READ(33) (NODEA(I),I=1,NUMNOD)
      call initrwf(33,iui,iur)
      do I=1,NUMNOD
      iui = iui+1
      nodea(I)=ipool(iui)
      enddo
      call endrwf(33,iui,iur)
        IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
c       READ(33) (EMATA(I),I=1,NUMMAT)
      call initrwf(33,iui,iur)
      do I=1,NUMMAT
      iur = iur+1
      emata(I)=rpool(iur)
      enddo
      call endrwf(33,iui,iur)
        ENDIF
 
        NADD=0
        NADM=0
        NUMEL=0
        NDMAX=0
        nlast = 0
        nnai = 0
        DO 2000 ITYP=1,NUMTYP
        IF (IDET(ITYP).EQ.0) GOTO 2000
C.......READ ELEM0 FILE
c        READ(25) NUM,NNODE
c        READ(25) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        NUM=NUMA(ITYP)
        NNODE=NNEA(ITYP)
        DO I=1,NUM
           DO J=1,NNODE
              NODE((I-1)*NNODE+J)=NODEA(NADD+(I-1)*NNODE+J)
           ENDDO
        ENDDO
        NADD=NADD+NUM*NNODE
c        WRITE(*,*) 'NUM =',NUM,' NNODE =',NNODE
c        WRITE(*,*) 'NODE ='
c        WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        NNE = NNODE
        IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
c        READ (25) MMATE,NMATE
c        READ (25) ((EMATE((I-1)*NMATE+J),J=1,NMATE),I=1,MMATE)
        MMATE=NMATA(ITYP)
        NMATE=NMATPA(ITYP)
        DO I=1,MMATE
           DO J=1,NMATE
c              write(*,*) 'nadm,kk == ',nadm,(I-1)*NMATE+J
              EMATE((I-1)*NMATE+J)=EMATA(NADM+(I-1)*NMATE+J)
           ENDDO
        ENDDO
        NADM=NADM+NMATE*MMATE
 
        NNE = NNE-1
        ENDIF
c        WRITE(*,*) 'NODE ='
c        WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
 
        DO 3000 NE=1,NUM
        L=0
        DO 3600 INOD=1,NNE
        NODI=NODE((NE-1)*NNODE+INOD)
        DO 3700 IDGF=1,KDGOF
        INV=NODVAR(IDGF,NODI)
        IF (INV.LE.0) GOTO 3700
        L=L+1
        LM(L)=INV
3700    CONTINUE
3600    CONTINUE
        NUMEL=NUMEL+1
C        WRITE (*,*) 'L,LM =',L
C        WRITE (*,'(1X,15I5)') (LM(I),I=1,L)
        NDMAX=MAX(NDMAX,L)
c----------------------------------------------------------
c
c       setup a link table which get the connection message of all unknowns
c
        do i = 1,l
        nunsi = lm(i)
c
        if(jdiagr(nunsi).eq.0) then
        jdiag(nunsi) = nlast + 1
        do j = 1,l
        nunsj = lm(j)
        nlast = nlast + 1
        if(nlast.ge.maxna) print *,'no enough memory for NA'
        na(nlast) = nunsj
        nlast = nlast + 1
        if(nlast.ge.maxna) print *,'no enough memory for NA'
        na(nlast) = nlast + 1
        end do
        jdiagr(nunsi) = l
        end if
c
        if(jdiagr(nunsi).gt.0) then
         nbtemp = jdiag(nunsi)
         do j = 1,l
          nunsj = lm(j)
            noyes = 0
            nbtemp = jdiag(nunsi)
            do k = 1,jdiagr(nunsi)
             kun = na(nbtemp)
             nbtemp0 = nbtemp
             nbtemp = na(nbtemp+1)
             if(kun.eq.nunsj) noyes = 1
            end do
c
            if(noyes.eq.0) then
             jdiagr(nunsi) = jdiagr(nunsi) + 1
             nlast = nlast + 1
             if(nlast.ge.maxna) print *,'no enough memory for NA'
             na(nbtemp0+1) = nlast
             na(nlast) = nunsj
             nlast = nlast + 1
             if(nlast.ge.maxna) print *,'no enough memory for NA'
             na(nlast) = nlast + 1
            end if
          end do
         end if
        end do
c
c----------------------------------------------------------
3000    CONTINUE
2000    CONTINUE
c------------------------------------------------------------
c
c       get the who Na memory structur and add again!
c
        jdiagaz(1) = 1
        do i = 1,neq
        jdiagaz(i+1) = jdiagaz(i) + jdiagr(i)
        end do
        maxa = jdiagaz(neq+1) -1
c
c        print *,(jdiag(i),i=1,neq)
c        print *,(jdiagr(i),i=1,neq)
c        print *,(jdiagaz(i),i=1,neq1)
c        do i = 1,neq
c        np = jdiag(i)
c        do j = 1,jdiagr(i)
c        lm(j) = na(np)
c        np = na(np+1)
c        end do
c        print *,(lm(j),j=1,jdiagr(i))
c        end do
c        print *,(na(i),i=1,maxa*2+2)
c
c       Then get the who matrix information again
c       Get Na dimension contains
c
          do i = 1,maxa
            na(i)=0
          end do
c
          do i=1,neq1
            jdiag(i)=0
          end do
c
c          do i=1,neq1
c            jdiagr(i)=0
c          end do
c
        NADD=0
        NADM=0
        NUMEL=0
        NDMAX=0
        DO 5000 ITYP=1,NUMTYP
        IF (IDET(ITYP).EQ.0) GOTO 5000
        NUM=NUMA(ITYP)
        NNODE=NNEA(ITYP)
        DO I=1,NUM
           DO J=1,NNODE
              NODE((I-1)*NNODE+J)=NODEA(NADD+(I-1)*NNODE+J)
           ENDDO
        ENDDO
        NADD=NADD+NUM*NNODE
        NNE = NNODE
        IF (MATERIAL.EQ.'Y' .OR. MATERIAL.EQ.'y') THEN
        MMATE=NMATA(ITYP)
        NMATE=NMATPA(ITYP)
        DO I=1,MMATE
           DO J=1,NMATE
              EMATE((I-1)*NMATE+J)=EMATA(NADM+(I-1)*NMATE+J)
           ENDDO
        ENDDO
        NADM=NADM+NMATE*MMATE
        NNE = NNE-1
        ENDIF
cc      WRITE(*,*) 'NODE ='
cc      WRITE(*,6) ((NODE((I-1)*NNODE+J),J=1,NNODE),I=1,NUM)
        DO 6000 NE=1,NUM
        L=0
        DO 6600 INOD=1,NNE
        NODI=NODE((NE-1)*NNODE+INOD)
        DO 6700 IDGF=1,KDGOF
        INV=NODVAR(IDGF,NODI)
        IF (INV.LE.0) GOTO 6700
        L=L+1
        LM(L)=INV
6700    CONTINUE
6600    CONTINUE
        NUMEL=NUMEL+1
C       WRITE (*,*) 'L,LM =',L
C       WRITE (*,'(1X,15I5)') (LM(I),I=1,L)
c
        call nzdmbsadd(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
c
6000    CONTINUE
5000    CONTINUE
c
c       order Na dimension
c
        call nzdmbsorder(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
c
c       get position of diag elements, error check
c
        call dmbsgetdiag(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
c
c
c.......WRITE SYS FILE
      ISTATUS = 0
      call openf(31,19,ISTATUS)
c       WRITE(31) NUMEL,NEQ,MAXA
      call initrwf(31,iui,iur)
      iui = iui+1
      ipool(iui)=numel
      iui = iui+1
      ipool(iui)=neq
      iui = iui+1
      ipool(iui)=maxa
      call endrwf(31,iui,iur)
c.......WRITE DIAG FILE
      ISTATUS = 0
      call openf(35,20,ISTATUS)
c       WRITE(35) (JDIAG(I),I=1,NEQ)
      call initrwf(35,iui,iur)
      do I=1,NEQ
      iui = iui+1
      ipool(iui)=jdiag(I)
      enddo
      call endrwf(35,iui,iur)
c       WRITE(35) (JDIAGAZ(I),I=1,NEQ1)
      call initrwf(35,iui,iur)
      do I=1,NEQ1
      iui = iui+1
      ipool(iui)=jdiagaz(I)
      enddo
      call endrwf(35,iui,iur)
c
c       write out Na file
c
      ISTATUS = 0
      call openf(36,21,ISTATUS)
c       write(36)(na(i),i=1,maxa)
      call initrwf(36,iui,iur)
      do i=1,maxa
      iui = iui+1
      ipool(iui)=na(i)
      enddo
      call endrwf(36,iui,iur)
        RETURN
        END
 
 
        subroutine nzdmbsadd(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
        implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        dimension na(maxa),jdiag(neq1),jdiagr(neq1),lm(l),jdiagaz(neq1)
        if (neq1.le.1) return
        do 300 i=1,l
        ii = lm(i)
        do 280 j=1,l
        jj = lm(j)
        noyes = 0
        nl = jdiagaz(ii)
        nu = jdiagaz(ii+1)-1
        do k = nl,nu
        if(na(k).eq.jj) noyes = 1
        end do
c
c       noyes = 0 add a new one
c
        if(noyes.eq.0) then
        jdiag(ii) = jdiag(ii) + 1
        nk = jdiagaz(ii) + jdiag(ii) -1
        if((nk.lt.jdiagaz(ii)).or.(nk.gt.(jdiagaz(ii+1)-1))) then
        print *,' Na dimension boundary error'
        end if
        na(nk) = jj
        end if
280     continue
300     continue
        return
        end
 
        subroutine nzdmbsorder(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
        implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension na(maxa),jdiag(neq1),jdiagr(neq1),lm(100),jdiagaz(neq1)
        if (neq1.le.1) return
        neq = neq1 - 1
        do i = 1,neq
        nl = jdiagaz(i)
        nu = jdiagaz(i+1)-1
        k = 0
        do j = nl,nu
        k = k + 1
        lm(k) = na(j)
        end do
c
        do l =1, k-1
        do n = l+1,k
        if(lm(n).lt.lm(l)) then
        if(lm(n).eq.lm(l)) print *,'Na contains error'
        ntemp = lm(l)
        lm(l) = lm(n)
        lm(n) = ntemp
        end if
        end do
        end do
c
        k = 0
        do j = nl,nu
        k = k + 1
        na(j) = lm(k)
        end do
c
        end do
c
c       Na dimension reorder ok, acent ordering
c
        return
        end
 
        subroutine dmbsgetdiag(neq1,maxa,na,jdiag,jdiagr,jdiagaz,l,lm)
        implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension na(maxa),jdiag(neq1),jdiagr(neq1),lm(100),jdiagaz(neq1)
        if (neq1.le.1) return
        neq = neq1 - 1
        do i = 1,neq
        nl = jdiagaz(i)
        nu = jdiagaz(i+1)-1
        noyes = 0
        do j = nl,nu
        if(na(j).eq.i) then
        jdiag(i) = j
        noyes = noyes + 1
        end if
        end do
c        if (noyes.ne.1) print *, ' Get diag position error'
        end do
        return
        end
 
