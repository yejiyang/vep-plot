      subroutine Mdmbs2dmsr(iblk)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
6     format (1x,10i5)
7     format (1x,6e12.5)
 
1001  format(1x,7i7)
        nsource=0
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
c ....open jraf file
c     open (21,file='jraf',form='unformatted',status='unknown')
c.....open na file
c     open (22,file='na',form='unformatted',status='unknown')
c.....open dmsr file
c     open (23,file='dmsr',form='unformatted',status='unknown')
c
c      do 2300 iblk=1,numblk
      ISTATUS = 24
      call openf(21,24,ISTATUS)
c       read(21) neq,maxa
      call initrwf(21,iui,iur)
      iui = iui+1
      neq=ipool(iui)
      iui = iui+1
      maxa=ipool(iui)
      call endrwf(21,iui,iur)
        neq1=neq+1
        maxa1=maxa+neq1
 
      knb1=neq1*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=neq1*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=maxa1*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      kna3=maxa1*1
      kna1=neq1*1
      kna2=neq1*1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      call dmbs2dmsr(iblk,numblk,maxa,maxa1,neq,neq1,nsource,
     *aa(kna0),aa(kna1),aa(kna2),ia(knb0),ia(knb1),
     *ia(knb2),
     *filename)
 
2300  continue
 
c     close(21)
c     close(22)
c     close(23)
c     close(24)
c     close(25)
c     close(26)
c     close(27)
      end
      subroutine dmbs2dmsr(iblk,numblk,maxa,maxa1,neq,neq1,nsource,
     *diag,f,a,jdiag,jdiagaz,na,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension jdiag(neq1),jdiagaz(neq1),na(maxa1),a(maxa1),diag(neq1),
     *          f(neq1)
c
6       format (1x,10i5)
7       format (1x,6e12.5)
c          DO I=1,NEQ1
c            JDIAG(I)=0
c          END DO
c          DO I=1,NEQ1
c            JDIAGAZ(I)=0
c          END DO
c          DO I=1,MAXA
c            NA(I)=0
c          END DO
C ...   TO READ JRAF FILE
c       READ (21) (JDIAG(I),I=1,NEQ)
      call initrwf(21,iui,iur)
      do I=1,NEQ
      iui = iui+1
      jdiag(I)=ipool(iui)
      enddo
      call endrwf(21,iui,iur)
c       READ (21) (JDIAGAZ(I),I=1,NEQ1)
      call initrwf(21,iui,iur)
      do I=1,NEQ1
      iui = iui+1
      jdiagaz(I)=ipool(iui)
      enddo
      call endrwf(21,iui,iur)
c       READ (21) (A(I),I=1,MAXA)
      call initrwf(21,iui,iur)
      do I=1,MAXA
      iur = iur+1
      a(I)=rpool(iur)
      enddo
      call endrwf(21,iui,iur)
c       READ (21) (F(I),I=1,NEQ)
      call initrwf(21,iui,iur)
      do I=1,NEQ
      iur = iur+1
      f(I)=rpool(iur)
      enddo
      call endrwf(21,iui,iur)
c        WRITE(*,*) 'EMASS ======='
c        WRITE(*,7) (A(I),I=1,MAXA)
c        WRITE(*,*) 'F ========'
c        WRITE(*,7) (F(I),I=1,NEQ)
c
c
c...... READ Na FILE
      ISTATUS = 21
      call openf(22,21,ISTATUS)
c       READ(22) (NA(i),i=1,MAXA)
      call initrwf(22,iui,iur)
      do i=1,MAXA
      iui = iui+1
      na(i)=ipool(iui)
      enddo
      call endrwf(22,iui,iur)
c
c        WRITE(*,*) 'Einform ======='
c        WRITE(*,7) (NA(I),I=1,MAXA)
c
c        DMSR----distributed modified sparse row
c        change the distrubuted gross stiff matrix to DMSR format
c
        do i = 1,neq
        diag(i) = a(jdiag(i))
        end do
c
        nonzero = 0
        ka = 0
        do i = 1,neq
        nl = jdiagaz(i)
        nu = jdiagaz(i+1)-1
        jdiagaz(i) = nonzero + 1
        do j = nl,nu
        ka = ka + 1
        if(na(j).eq.i) then
        if(dabs(diag(i)-a(ka)).gt.1.0e-35) print *,'DMSR error'
        end if
        if((na(j).gt.0).and.(na(j).ne.i)) then
        nonzero = nonzero + 1
        a(nonzero) = a(ka)
        na(nonzero) = na(ka)
        end if
        end do
        end do
        jdiagaz(neq+1)  = nonzero + 1
 
c
        do i = nonzero,1,-1
        na(i+neq1) = na(i)
        a(i+neq1) = a(i)
        end do
c
        do i = 1,neq1
        na(i) = jdiagaz(i) + neq1
        end do
c
        do i = 1,neq
        a(i) = diag(i)
        end do
        a(neq1) = 0.0d0
        maxa = na(neq1)
c
c
c       write back DCSR format gross matrix to dcsr
c
      ISTATUS = 0
      call openf(23,26,ISTATUS)
c       write(23) neq,maxa
      call initrwf(23,iui,iur)
      iui = iui+1
      ipool(iui)=neq
      iui = iui+1
      ipool(iui)=maxa
      call endrwf(23,iui,iur)
c       write(23) (JDIAG(I),I=1,NEQ)
      call initrwf(23,iui,iur)
      do I=1,NEQ
      iui = iui+1
      ipool(iui)=jdiag(I)
      enddo
      call endrwf(23,iui,iur)
c       write(23) (JDIAGAZ(I),I=1,NEQ1)
      call initrwf(23,iui,iur)
      do I=1,NEQ1
      iui = iui+1
      ipool(iui)=jdiagaz(I)
      enddo
      call endrwf(23,iui,iur)
c       write(23) (NA(I),I=1,MAXA)
      call initrwf(23,iui,iur)
      do I=1,MAXA
      iui = iui+1
      ipool(iui)=na(I)
      enddo
      call endrwf(23,iui,iur)
c       write(23) (A(I),I=1,MAXA)
      call initrwf(23,iui,iur)
      do I=1,MAXA
      iur = iur+1
      rpool(iur)=a(I)
      enddo
      call endrwf(23,iui,iur)
c       write(23) (F(I),I=1,NEQ)
      call initrwf(23,iui,iur)
      do I=1,NEQ
      iur = iur+1
      rpool(iur)=f(I)
      enddo
      call endrwf(23,iui,iur)
c
c       print out the dcsr matrix for error check
c
c        write(*,*) 'neq, new maxa=====',neq,maxa
c        write(*,*) 'new jdiag file ====='
c        write(*,*) (JDIAG(I),I=1,NEQ)
c        write(*,*) 'new jdiagaz file ====='
c        write(*,*) (JDIAGAZ(I),I=1,NEQ1)
c        write(*,*) 'new na file ====='
c        write(*,*) (na(I),I=1,maxa)
c        write(*,*) 'new a file ====='
c        write(*,*) (a(I),I=1,maxa)
c        write(*,*) 'new f file ====='
c        write(*,*) (f(I),I=1,NEQ)
c
c
        return
        end
