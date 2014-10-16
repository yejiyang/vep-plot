      subroutine Meddm(iblk)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
16     format (1x,10i5)
17     format (1x,6e12.5)
c      if(iblk.ne.1) return
c      call timer(2,1)
        nsource=0
      ksml=30000
      kelem=500000
      maxnd=300
c.....open time file
c     open(21,file='time',form='unformatted',status='old')
      ISTATUS = 11
      call openf(21,11,ISTATUS)
c     read(21) tmax,dt,time,it
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
c     close(21)
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
c     open (21,file='sys',form='unformatted',status='unknown')
c.....open nodvar file
c     open (22,file='nv',form='unformatted',status='unknown')
c.....open coor0 file
c     open (23,file='scoor0',form='unformatted',status='unknown')
c.....open bfd file
c     open (24,file='bfd',form='unformatted',status='unknown')
c.....open elem0 file
c     open (25,file='selem0',form='unformatted',status='unknown')
c.....open jraf file
c     open (27,file='jraf',form='unformatted',status='unknown')
c.....open disp0 file
c     open (28,file='sdisp0',form='unformatted',status='unknown')
c.....open diag file
c     open (29,file='jdiag',form='unformatted',status='unknown')
c.....open na file
c     open (31,file='na',form='unformatted',status='unknown')
c.....open ex_mcs file
c     open (10,file='ex_mcs',form='unformatted',status='unknown')
c.... open unod file
      IF(IT.EQ.0) THEN
      ELSE
      ENDIF
c      do 2300 iblk=1,numblk
        nt1=21
        nt2=22
      ISTATUS = 19
      call openf(21,19,ISTATUS)
c       READ(21) NSKP,NEQ,MAXA
      call initrwf(21,iui,iur)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      neq=ipool(iui)
      iui = iui+1
      maxa=ipool(iui)
      call endrwf(21,iui,iur)
      ISTATUS = 12
      call openf(22,12,ISTATUS)
c       read(22) nskp
      call initrwf(22,iui,iur)
      iui = iui+1
      nskp=ipool(iui)
      call endrwf(22,iui,iur)
        NEQ1=NEQ+1
        MAXA1=MAXA+1
        MAXNA=1
        MAXB=1
      ISTATUS = 2
      call openf(23,2,ISTATUS)
c       read(23) nskip,ncoor
      call initrwf(23,iui,iur)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      ncoor=ipool(iui)
      call endrwf(23,iui,iur)
      ISTATUS = 4
      call openf(28,4,ISTATUS)
c       read(28) knode,kdgof
      call initrwf(28,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      call endrwf(28,iui,iur)
c       read(28) ((skip,j=1,kdgof),i=1,knode)
      call initrwf(28,iui,iur)
      do i=1,knode
      do j=1,kdgof
      iur = iur+1
      skip=rpool(iur)
      enddo
      enddo
      call endrwf(28,iui,iur)
        kvar=knode*kdgof
c        write(*,*) 'knode,kdgof,ncoor,kvar =='
c        write(*,16) knode,kdgof,ncoor,kvar
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
 
      knb6=kdgof*knode*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      kna5=kdgof*knode*1
      kna7=ncoor*knode*1
      kna6=neq1*1
      knb7=kelem*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      kna12=ksml*1
      knb8=maxnd*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      kna11=maxa1*1
      kna1=maxb*1
      knb11=neq1*1
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      knb5=numtyp*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      kna8=neq1*1
      knb13=maxa1*1
      if (knb13/2*2 .lt. knb13) knb13=knb13+1
      kna9=neq1*1
      kna10=neq1*1
      kna2=kvar*1
      kna3=neq1*1
      knb9=kvar*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      knb10=numnod*1
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      kna4=nummat*1
      knb1=numtyp*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=numtyp*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numtyp*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb4=numtyp*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb12=neq1*1
      if (knb12/2*2 .lt. knb12) knb12=knb12+1
      kna13=150000*1
      kna14=numnod*6
      kna15=knode*kdgof
      kna16=knode*kdgof
      kna17=knode*6
      knb14=numnod
      if (knb14/2*2 .lt. knb14) knb14=knb14+1
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
      knb10=knb10+knb9
      knb11=knb11+knb10
      knb12=knb12+knb11
      knb13=knb13+knb12
      knb14=knb14+knb13
      call eddm(knode,kdgof,ncoor,kelem,numtyp,maxa,
     *iblk,numblk,ksml,it,time,dt,kkk,maxnd,
     *maxna,maxt,neq,neq1,maxb,maxa1,kvar,nsource,
     *lgio,numnod,nummat,aa(kna0),aa(kna1),aa(kna2),
     *aa(kna3),aa(kna4),aa(kna5),aa(kna6),aa(kna7),
     *aa(kna8),aa(kna9),aa(kna10),aa(kna11),aa(kna12),ia(knb0),
     *ia(knb1),ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *ia(knb6),ia(knb7),ia(knb8),ia(knb9),ia(knb10),
     *ia(knb11),ia(knb12),
     *filename,aa(kna13),aa(kna14),aa(kna15),
     *ia(knb13),aa(kna16))
 
2300  continue
c      call timer(2,2)
 
c     close(21)
c     close(22)
c     close(23)
c     close(24)
c     close(25)
c     close(27)
c     close(28)
c     close(29)
c     close(10)
c     close(30)
c     close(31)
      end
      subroutine eddm(knode,kdgof,ncoor,kelem,numtyp,maxa,
     *iblk,numblk,ksml,it,time,dt,kkk,maxnd,
     *maxna,maxt,neq,neq1,maxb,maxa1,kvar,nsource,
     *lgio,numnod,nummat,b,emadd,emnew,emata,u,
     *f,coor,uu,exem,exec,a,sml,emate,numa,
     *nnea,mmta,nmta,idet,nodvar,node,lm,iemdex,
     *nodea,jdiag,jdiagaz,na,
     *filename,ess,eu1,ev1,ilabel,smain)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension nodvar(kdgof,knode),u(kdgof,knode),coor(ncoor,knode),
     *
     *f(neq1),node(kelem),sml(ksml),emate(nummat),lm(maxnd),
     *a(maxa1),b(maxb),jdiag(neq1),idet(numtyp),uu(neq1),na(maxa1),
     *exem(neq1),exec(neq1),emadd(kvar),emnew(neq1),iemdex(kvar),
     *NODEA(NUMNOD),EMATA(NUMMAT),NUMA(NUMTYP),NNEA(NUMTYP),
     *MMTA(NUMTYP),NMTA(NUMTYP),jdiagaz(neq1)
      dimension eu1(kdgof,knode),ev1(kdgof,knode),ess(6*NUMNOD)
      dimension ilabel(NUMNOD),smain(6,knode)
      CHARACTER*1 MATERIAL
 
6       format (1x,10i5)
7       format (1x,6e12.5)
 
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
      call endrwf(30,iui,iur)
c       write(*,*) 'idet ===================='
c       write(*,6) (idet(i),i=1,numtyp)
C.....  READ ELEM0 FILE
      ISTATUS = 7
      call openf(25,7,ISTATUS)
c       READ(25) (NODEA(I),I=1,NUMNOD)
      call initrwf(25,iui,iur)
      do I=1,NUMNOD
      iui = iui+1
      nodea(I)=ipool(iui)
      enddo
      call endrwf(25,iui,iur)
        IF(nummat.gt.0) THEN
c       READ(25) (EMATA(I),I=1,NUMMAT)
      call initrwf(25,iui,iur)
      do I=1,NUMMAT
      iur = iur+1
      emata(I)=rpool(iur)
      enddo
      call endrwf(25,iui,iur)
        ENDIF
c ....  read nv file
c       read (22) ((nodvar(i,j),i=1,kdgof),j=1,knode)
      call initrwf(22,iui,iur)
      do j=1,knode
      do i=1,kdgof
      iui = iui+1
      nodvar(i,j)=ipool(iui)
      enddo
      enddo
      call endrwf(22,iui,iur)
c       write (*,*) 'nodvar ========================'
c       write (*,6) ((nodvar(i,j),i=1,kdgof),j=1,knode)
c ....  read coor file
c       read (23) ((coor(i,j),i=1,ncoor),j=1,knode)
      call initrwf(23,iui,iur)
      do j=1,knode
      do i=1,ncoor
      iur = iur+1
      coor(i,j)=rpool(iur)
      enddo
      enddo
      call endrwf(23,iui,iur)
c       write(*,*) 'coor ==========================',iblk
c       write(*,7) ((coor(j,i),j=1,ncoor),i=1,knode)
c ....  read bfd file
      ISTATUS = 13
      call openf(24,13,ISTATUS)
c       read (24) ((u(j,i),j=1,kdgof),i=1,knode)
      call initrwf(24,iui,iur)
      do i=1,knode
      do j=1,kdgof
      iur = iur+1
      u(j,i)=rpool(iur)
      enddo
      enddo
      call endrwf(24,iui,iur)
c       write (*,*) 'bf ===========================',iblk
c       write(*,7) ((u(j,i),j=1,kdgof),i=1,knode)


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
C     strw6.tmp for stress at previous time step
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
      ISTATUS = 33
      call openf(43,33,ISTATUS)
      call initrwf(43,iui,iur)
      do j = 1,kdgof
      do i = 1,knode
      iur = iur+1
      eu1(j,i) = rpool(iur)
      end do
      end do
c      
      do j = 1,kdgof
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

c.....  read diag file
      ISTATUS = 20
      call openf(29,20,ISTATUS)
c       read(29) (jdiag(i),i=1,neq)
      call initrwf(29,iui,iur)
      do i=1,neq
      iui = iui+1
      jdiag(i)=ipool(iui)
      enddo
      call endrwf(29,iui,iur)
c       read(29) (jdiagaz(i),i=1,neq1)
      call initrwf(29,iui,iur)
      do i=1,neq1
      iui = iui+1
      jdiagaz(i)=ipool(iui)
      enddo
      call endrwf(29,iui,iur)
c       write (*,*) 'jdiag jdiagaz ================'
c        write(*,*) (jdiag(i),i=1,neq)
c        write(*,*) (jdiagaz(i),i=1,neq1)
c.....  read na file
      ISTATUS = 21
      call openf(31,21,ISTATUS)
c       read(31) (na(i),i=1,maxa)
      call initrwf(31,iui,iur)
      do i=1,maxa
      iui = iui+1
      na(i)=ipool(iui)
      enddo
      call endrwf(31,iui,iur)
c       write (*,*) 'na============================'
c        do i = 1,neq
c        nl = jdiagaz(i)
c        nu = jdiagaz(i+1)-1
c        write(*,*) (na(j),j=nl,nu)
c        end do
 
        do 11 i=1,neq
11      continue
        do i=1,maxa
        a(i)=0.0d0
        enddo
        do i=1,maxb
        b(i)=0.0d0
        enddo
        do i=1,neq
          exem(i)=0.0d0
          exec(i)=0.0d0
        enddo
        numel=0
        NODADD=0
        NMADD=0
        do 200 ityp=1,numtyp
        if (idet(ityp).eq.0) goto 200
c ....  read elem file
c        read (25) num,nnode
c        read (25) ((node((i-1)*nnode+j),j=1,nnode),i=1,num)
        NUM=NUMA(ITYP)
        NNODE=NNEA(ITYP)
        DO I=1,NUM
           DO J=1,NNODE
              NODE((I-1)*NNODE+J)=NODEA(NODADD+(I-1)*NNODE+J)
           ENDDO
        ENDDO
        NODADD=NODADD+NUM*NNODE
c        if(iblk.eq.1) then
c        do i = 1,num
c        write(*,8888) (node((i-1)*nnode+j),j=1,nnode)
c        end do
c        end if
c8888    format(10I10)
c        return
c        stop
        knum=num*nnode
        if (knum.gt.kelem) then
        write(*,*) 'it is error: knum gt kelem'
        stop 0000
        endif
c       write(*,*) 'node ='
c       write(*,6) ((node((i-1)*nnode+j),j=1,nnode),i=1,num)
        nne = nnode
      nne = nne-1
        k=0
        do 115 j=1,nne
        jnod = node(j)
        do 115 l=1,kdgof
        if (nodvar(l,jnod).ne.0) k=k+1
115     continue
c        write(*,*) 'k =',k
        kk=k*k
        k0=1
        k1=k0+k*k
        k2=k1+k
        k3=k2+k
        k4=k3+k
        k5=k4+k*k
        k6=k5+k*k
        if (k6.ge.ksml) then
        write(*,*) 'it is error: k5 gt ksml'
        stop 0000
        endif

        call etsub(iblk,knode,kdgof,it,kelem,k,kk,nnode,nne,
     *  neq1,numel,ityp,ncoor,num,time,dt,nodvar,coor,node,emate,
     *  a,b,JDIAG,jdiagaz,maxa,maxb,neq,numtyp,na,exem,exec,
     *  emata,mmta,nmta,nmadd,
     *sml(k0),sml(k1),sml(k2),sml(k3),sml(k4),sml(k5),
     &
     *u,nummat,ess,eu1,ev1,ilabel,smain)
        numel = numel + num
200     continue
        do i=1,neq
        f(i)=0.0d0
        enddo
        do 220 i=1,knode
        do 210 j=1,kdgof
        ij=nodvar(j,i)
        if (ij.le.0) goto 210
        if (ij.gt.neq) neq = ij
        f(ij)=f(ij)+u(j,i)
210     continue
220     continue
c       write(*,*) 'a(ii) ======='
c       write(*,7) (a(jdiag(ii)),ii=1,neq)
c       write(*,*) 'b(ii) ======='
c       write(*,7) (b(jdiag(ii)),ii=1,neq)
c       write(*,*) '******** call redu ***********'
c...... read emassadd file
ccccc        read(15,err=225) (emadd(i),i=1,kvar)
        do i=1,kvar
          emadd(i)=0.0d0
        enddo
        goto 226
225     continue
        do i=1,kvar
          emadd(i)=0.0d0
        enddo
226     continue
c        call reduin(a,f,jdiag,neq,maxa,3)
      ISTATUS = 0
      call openf(27,24,ISTATUS)
c       write(27) neq,maxa
      call initrwf(27,iui,iur)
      iui = iui+1
      ipool(iui)=neq
      iui = iui+1
      ipool(iui)=maxa
      call endrwf(27,iui,iur)
c       write(27) (jdiag(ii),ii=1,neq)
      call initrwf(27,iui,iur)
      do ii=1,neq
      iui = iui+1
      ipool(iui)=jdiag(ii)
      enddo
      call endrwf(27,iui,iur)
c       write(27) (jdiagaz(ii),ii=1,neq1)
      call initrwf(27,iui,iur)
      do ii=1,neq1
      iui = iui+1
      ipool(iui)=jdiagaz(ii)
      enddo
      call endrwf(27,iui,iur)
c       write(27) (a(ii),ii=1,maxa)
      call initrwf(27,iui,iur)
      do ii=1,maxa
      iur = iur+1
      rpool(iur)=a(ii)
      enddo
      call endrwf(27,iui,iur)
c       write(27) (f(ii),ii=1,neq)
      call initrwf(27,iui,iur)
      do ii=1,neq
      iur = iur+1
      rpool(iur)=f(ii)
      enddo
      call endrwf(27,iui,iur)
cc      write(*,*) 'a =========='
cc      write(*,7) (a(jdiag(ii)),ii=1,neq)
cc      write(*,*) 'f =========='
cc      write(*,7) (f(ii),ii=1,neq)

c
c     renew label file
c
      ISTATUS = 37
      call openf(47,37,ISTATUS)
      call initrwf(47,iui,iur)
      do i = 1,NUMNOD
      iui = iui+1
      ipool(iui) = ilabel(i)
      end do
      call endrwf(47,iui,iur)

        return
        end
 
        subroutine etsub(iblk,knode,kdgof,it,kelem,k,kk,nnode,nne,
     *  neq1,numel,ityp,ncoor,num,time,dt,nodvar,coor,node,emate,
     *  a,b,jdiag,jdiagaz,maxa,maxb,neq,numtyp,na,exem,exec,
     *  emata,mmta,nmta,nmadd,
     *es,em,ec,ef,Estifn,Estifv,
     *  u,nummat,ess,eu1,ev1,ilabel,smain)
        implicit real*8 (a-h,o-z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        dimension nodvar(kdgof,knode),coor(ncoor,knode),
     *  u(kdgof,knode),emate(nummat),node(kelem),
     *  exem(neq1),exec(neq1),emata(nummat),mmta(numtyp),nmta(numtyp),
     *es(k,k),em(k),ec(k),ef(k),Estifn(k,k),Estifv(kk),
     *
     *  a(maxa),b(maxb),jdiag(neq1),jdiagaz(neq1),na(maxa)
        dimension r(3000),prmt(3000),coef(3000),lm(3000),
     *  emassm(3000),emassc(3000)
        dimension eu1(kdgof,knode),ev1(kdgof,knode),ess(6*kelem)
        dimension fstr(6,8),ilabel(kelem)
        dimension label(8),smain(6,knode)

7       format (1x,10i5)
8       format (1x,6e12.5)
c        print *,'nummat=======k========',nummat,k 
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

      DO 1000 NE=1,NUM
      NR=0
      DO 130 J=1,NNE
      JNOD = NODE((NE-1)*NNODE+J)
      IF (JNOD.LT.0) JNOD = -JNOD
      do 110 l=1,kdgof
      coef(j+(l-1)*nne)=eu1(l,jnod)
110   continue
      coef(j+(0+1*kdgof)*nne)=(smain(1,jnod)
     &+smain(2,jnod)+smain(3,jnod))/3
      coef(j+(1+1*kdgof)*nne)=0.d0
      coef(j+(2+1*kdgof)*nne)=0.d0
      coef(j+(3+1*kdgof)*nne)=0.d0
      coef(j+(4+1*kdgof)*nne)=0.d0
      coef(j+(5+1*kdgof)*nne)=0.d0
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
 
      if(ityp.eq.1) then
       do j=1,8
       do i=1,6
       fstr(i,j)=ess((ne-1)*6*8+(j-1)*6+i)
       enddo
       label(j)=ilabel((ne-1)*8+j)
       enddo
       call suc8g2(r,coef,prmt,es,em,ec,ef,nnee,
     *fstr,label)
       do j=1,8
       ilabel((ne-1)*8+j)=label(j)
       enddo
      endif

      if(ityp.eq.2) call hugq4g2(r,coef,prmt,es,em,ec,ef,nnee)


      do 201 i=1,k
      do 201 j=1,k
      Estifn(i,j)=0.0d0
201   continue
      do 202 i=1,k
      Estifn(i,i)=Estifn(i,i)
      do 202 j=1,k
      Estifn(i,j)=Estifn(i,j)+Es(i,j)
202   continue
 
      l=0
      m=0
      i=0
      do 700 inod=1,nne
      nodi=node((ne-1)*nnode+inod)
      do 600 idgf=1,kdgof
      inv=nodvar(idgf,nodi)
      if (inv.eq.0) goto 600
      i=i+1
      if (inv.lt.0) goto 305
      l=l+1
      lm(l)=inv
      u(idgf,nodi)=u(idgf,nodi)
     *+Ef(i)
c      emassm(l)=em(i)
c      emassc(l)=ec(i)
305   j=0
      do 500 jnod=1,nne
      nodj=node((ne-1)*nnode+jnod)
      do 400 jdgf=1,kdgof
      jnv=nodvar(jdgf,nodj)
      if (jnv.eq.0) goto 400
      j=j+1
 
      if (jnv.lt.0) goto 400
      if (inv.lt.0) goto 310
      m=m+1
      Estifv(m)=Estifn(i,j)
310   continue

      u(jdgf,nodj)=u(jdgf,nodj)
     *+es(j,i)*eu1(idgf,nodi)
c 
        if (inv.lt.0)
     *  U(JDGF,NODJ)=U(JDGF,NODJ)-ESTIFN(J,I)*U(IDGF,NODI)
400   continue
500   continue
600   continue
700   continue
      if(ityp.eq.3) print *,'aaaaaaaaaaaaaaaaaaa' 
      call nzaddmbs(a,na,jdiag,jdiagaz,l,lm,estifv,neq1,maxa)
      if(ityp.eq.3) print *,'eeeeene ityp numeee',ne,ityp,num
      do i=1,l
        j=lm(i)
c        exem(j)=exem(j)+emassm(i)
c        exec(j)=exec(j)+emassc(i)
      enddo
1000  continue
        end
 
        SUBROUTINE NZADDMBS(A,NA,JDIAG,JDIAGAZ,ND,LM,ESTIF,NEQ1,MAXA)
        IMPLICIT real*8 (A-H,O-Z)
        include 'memalloc.h'
        common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION A(MAXA),JDIAG(NEQ1),JDIAGAZ(NEQ1),LM(ND),ESTIF(ND,ND)
        DIMENSION NA(MAXA)
        IF (NEQ1.LE.1) RETURN
        DO 300 I=1,ND
        II = LM(I)
        DO 280 J=1,ND
        JJ = LM(J)
c
c       add here
c
        nl = jdiagaz(ii)
        nu = jdiagaz(ii+1)-1
        noyes = 0
        do k = nl,nu
        if(na(k).eq.jj) then
        noyes = 1
        nup = k
        end if
        end do
c        if(noyes.eq.0) print *,' Na pointer dimension error'
        a(nup) = a(nup) + estif(i,j)
280     CONTINUE
300     CONTINUE
        RETURN
        END
 
 
        SUBROUTINE ADDMBS(A,JDIAG,JDIAGAZ,ND,LM,ESTIF,NEQ1,MAXA)
        IMPLICIT real*8 (A-H,O-Z)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
        DIMENSION A(MAXA),JDIAG(NEQ1),JDIAGAZ(NEQ1),LM(ND),ESTIF(ND,ND)
        IF (NEQ1.LE.1) RETURN
        DO 300 I=1,ND
        II = LM(I)
        DO 280 J=1,ND
        JJ = LM(J)
        k = jdiag(ii) - ii + jj
        if((k.lt.jdiagaz(ii)).or.(k.gt.(jdiagaz(ii+1)-1))) then
        print *,' Na dimension boundary error'
        end if
        A(K) = A(K) + ESTIF(I,J)
c
c        IF (II.LT.JJ) THEN
c        K = JDIAG(JJ) - JJ + II
c        A(K) = A(K) + ESTIF(I,J)
c        ENDIF
c        IF (II.GE.JJ) THEN
c        K = JDIAG(II) - II + JJ
c        A(K) = A(K) + ESTIF(I,J)
c        ENDIF
c
280     CONTINUE
300     CONTINUE
        RETURN
        END
