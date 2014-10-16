      subroutine Muddm(iblk,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
      character*8 xyzr(6),string(6)
        nsource=0
c ....  open lmddm file
c       open (21,file='lmddm',form='unformatted',status='old')
      ISTATUS = 1
      call openf(21,1,ISTATUS)
c       read(21) numblk,numtyp,nskp,kdgof,nskp,nskp,nskp
      call initrwf(21,iui,iur)
      iui = iui+1
      numblk=ipool(iui)
      iui = iui+1
      numtyp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      iui = iui+1
      nskp=ipool(iui)
      call endrwf(21,iui,iur)
c       close(21)
c ....  open time file
c       open(21,file='time',form='unformatted',status='old')
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
        t=time
c ....  open sys file
c       open(21,file='sys',form='unformatted',status='old')
c ....  open bfd file
c       open(22,file='bfd',form='unformatted',status='old')
c ....  open nv file
c       open(23,file='nv',form='unformatted',status='old')
c ....  open u file
c       open(24,file='u',form='unformatted',status='old')
c ....  open displace file
        if(iblk.eq.0) then
       open(25,file='displace',form='formatted',status='unknown')
        endif
c       open(25,file='displace',form='formatted',status='unknown')
c ....  open disp0 file
c       open(27,file='sdisp0',form='unformatted',status='old')
c.....  open dispname file
cc        open(28,file=' ',form='formatted',status='unknown')
c...... open coor0 file
c       open(30,file='scoor0',form='formatted',status='unknown')
c...... open munod0 file
        if(iblk.eq.0) then
2        open(20,file='munod0',form='unformatted',status='unknown')
        else
c       open(20,file='munod0',form='unformatted',status='unknown')
        end if
c
c       initial value for the format of output solution file
c
        do i=1,kdgof
          xyzr(i)='u'
        enddo
        k=0
        do i=1,kdgof
        k=k+1
        string(k) = xyzr(i)
        end do
c        write (25,1500) (string(i),i=1,k)
1500    format (1x,4hnode,2x,6(4x,a8))
c
c
c     open(11,file='unod',form='unformatted',status='unknown')
c
c        do 2300 iblk=0,numblk-1
c.....  read sys file
      ISTATUS = 19
      call openf(21,19,ISTATUS)
c       read(21) numel,neq
      call initrwf(21,iui,iur)
      iui = iui+1
      numel=ipool(iui)
      iui = iui+1
      neq=ipool(iui)
      call endrwf(21,iui,iur)
c        write(*,*) 'numel,neq =',numel,neq
c
c       read u file (lsol file)
c
      ISTATUS = 27
      call openf(24,27,ISTATUS)
c       read(24) keq
      call initrwf(24,iui,iur)
      iui = iui+1
      keq=ipool(iui)
      call endrwf(24,iui,iur)
c.....  read disp0 file
      ISTATUS = 4
      call openf(27,4,ISTATUS)
c       read(27) knode,kdgof
      call initrwf(27,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      call endrwf(27,iui,iur)
c       read(27) ((xskp,j=1,kdgof),i=1,knode)
      call initrwf(27,iui,iur)
      do i=1,knode
      do j=1,kdgof
      iur = iur+1
      xskp=rpool(iur)
      enddo
      enddo
      call endrwf(27,iui,iur)
        kvar=knode*kdgof
c
c       check if kvar equal to keq,
c
c        if(keq.ne.kvar) then
        if(keq.ne.neq) then
        write(*,*) 'There must be errors exist between nzazsolv
     &  and u file, check the data communication part'
        return
        end if
c.....  read coor0 file
      ISTATUS = 2
      call openf(30,2,ISTATUS)
c       read(30) nskip,kcoor
      call initrwf(30,iui,iur)
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      kcoor=ipool(iui)
      call endrwf(30,iui,iur)
c        write(*,*) 'knode,kdgof,kvar ===',knode,kdgof,kvar
c
      knb1=kdgof*knode*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      kna3=kdgof*knode*1
      kna2=kvar*1
      kna1=kcoor*knode*1
      kna0=1
      kna4=kdgof*knode*1
      kna5=kdgof*knode*1
      kna6=kdgof*knode*1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
      knb0=1
      knb1=knb1+knb0
      call uddm(knode,kdgof,kvar,neq,kcoor,t,
     *dt,numblk,it,nsource,iblk,iter,tmax,xyzr,
     *keq,aa(kna0),aa(kna1),aa(kna2),
     *ia(knb0),
     *filename,aa(kna3),aa(kna4),aa(kna5),istop,kend)
c
2300    continue
c     close(11)
c       close(21)
c       close(22)
c       close(23)
c       close(24)
        if(iblk.eq.0) then
        close(25)
        endif
c       close(25)
c       close(27)
c       close(30)
        if(iblk.eq.0) then
5        close(20)
        else
        end if
      end
      subroutine uddm(knode,kdgof,kvar,neq,kcoor,t,
     *dt,numblk,it,nsource,iblk,iter,tmax,xyzr,
     *keq,coor,uvar,eu,nodvar,
     *filename,eu1,ev1,ue,istop,kend)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension nodvar(kdgof,knode),rtemp(10),itemp(10),
     *  eu(kdgof,knode),eu1(kdgof,knode),ev1(kdgof,knode),
     *  uvar(kvar),r(3),coor(kcoor,knode),ue(kdgof,knode)
        character*8 xyzr(6),string(6)
 
6       format (1x,10i5)
7       format (1x,6e12.5)
c.......read bfd file
      ISTATUS = 13
      call openf(22,13,ISTATUS)
c       read (22) ((eu(i,j),i=1,kdgof),j=1,knode)
      call initrwf(22,iui,iur)
      do j=1,knode
      do i=1,kdgof
      iur = iur+1
      eu(i,j)=rpool(iur)
      enddo
      enddo
      call endrwf(22,iui,iur)
cc      write (*,*) 'u ='
cc      write(*,7) ((eu(i,j),i=1,kdgof),j=1,knode)
c...... read coor file
c       read(30) ((coor(i,j),i=1,kcoor),j=1,knode)
      call initrwf(30,iui,iur)
      do j=1,knode
      do i=1,kcoor
      iur = iur+1
      coor(i,j)=rpool(iur)
      enddo
      enddo
      call endrwf(30,iui,iur)
c...... read nv file
      ISTATUS = 12
      call openf(23,12,ISTATUS)
c       read (23) nskp
      call initrwf(23,iui,iur)
      iui = iui+1
      nskp=ipool(iui)
      call endrwf(23,iui,iur)
c       read (23) ((nodvar(i,j),i=1,kdgof),j=1,knode)
      call initrwf(23,iui,iur)
      do j=1,knode
      do i=1,kdgof
      iui = iui+1
      nodvar(i,j)=ipool(iui)
      enddo
      enddo
      call endrwf(23,iui,iur)
cc      write (*,*) 'nodvar ='
cc      write (*,6) ((nodvar(i,j),i=1,kdgof),j=1,knode)
c
c       I have modified here for convenient, old code left here
c
c...... read u file
c        read (24) (uvar(i),i=1,neq)
c       write (*,*) 'uvar ='
c       write (*,7) (uvar(i),i=1,neq)
c        do 200 inod=1,knode
c        do 100 idfg=1,kdgof
c        n=nodvar(idfg,inod)
c        if (n.le.0) goto 100
c        eu(idfg,inod)=uvar(n)
c100     continue
c200     continue
c
c...... read u file
c       read (24) (uvar(i),i=1,keq)
      call initrwf(24,iui,iur)
      do i=1,keq
      iur = iur+1
      uvar(i)=rpool(iur)
      enddo
      call endrwf(24,iui,iur)
c        write (*,*) 'iblk==== ',iblk,'uvar ='
c        write (*,7) (uvar(i),i=1,keq)
        k = 0
        do 200 inod=1,knode
        do 100 idfg=1,kdgof
        n=nodvar(idfg,inod)
        if (n.le.0) goto 100
        k = k + 1
        eu(idfg,inod)=uvar(k)
100     continue
200     continue

C
C     unod for du and ddu at previous time step
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
      do j = 1,kdgof
      do i = 1,knode
      iur = iur+1
      ev1(j,i) = rpool(iur)
      end do
      end do      
      call endrwf(43,iui,iur)

      eaa=0.d0
      eab=0.d0
      ebb=0.d0
      do j = 1,kdgof
      do i = 1,knode
      ue(j,i)=eu(j,i)-eu1(j,i)
      eaa=eaa+ue(j,i)*ue(j,i)
      eab=eab+ue(j,i)*ev1(j,i)
      ebb=ebb+ev1(j,i)*ev1(j,i)
      end do
      end do      

      if (iblk.eq.0) then
       do ibk = 1,numblk-1
        call recvr(nsource,ibk,rtemp(1))
        eaa=eaa+rtemp(1)
        call recvr(nsource,ibk,rtemp(1))
        eab=eab+rtemp(1)
        call recvr(nsource,ibk,rtemp(1))
        ebb=ebb+rtemp(1)
       enddo
       ecc=1.0d0
cc       if (ebb.gt.1.0d-11) then
       if (kend.ne.1) then
c        read numit
         NUMUNIT = 48
         ISTATUS = 38
         call openf(48,38,ISTATUS)
         call initrwf(48,iui,iur)
         iui = iui+1
         numit = ipool(iui) 
         call endrwf(48,iui,iur)
         numit=numit+1
         rab=dsqrt(eaa)*dsqrt(ebb)
         if (eaa.gt.ebb) ecc=dsqrt(ebb/eaa)
         if (eab.gt.0.8*rab) ecc=ecc*2
         if (eab.lt.0.3*rab) ecc=ecc/2
         if (eab.lt.-0.4*rab) ecc=ecc/2
       else
         numit=1
       endif
       if (ecc.gt.1.0d0) ecc=1.0d0
       if (ecc.lt.1.d-3) ecc=1.0d0
       do ibk = 1,numblk-1
        call sendr(ibk,nsource,ecc)
       enddo
      else
       call sendr(nsource,iblk,eaa)
       call sendr(nsource,iblk,eab)
       call sendr(nsource,iblk,ebb)
       call recvr(iblk,nsource,rtemp(1))
       ecc=rtemp(1)
      endif
      err=0.d0
      do j = 1,kdgof
      do i = 1,knode
      err=err+ue(j,i)*ue(j,i)
      ue(j,i)=ue(j,i)*ecc
      eu(j,i)=eu1(j,i)+ue(j,i)
      end do
      end do
      if (iblk.eq.0) then
       do ibk = 1,numblk-1
        call recvr(nsource,ibk,rtemp(1))
        err=err+rtemp(1)
       enddo
       write(*,*) 'err ecc  numit =',err,numit,ecc
       write(*,*) 'aa  bb  ab rab =',eaa,ebb,eab,rab
       NMI=4
       ERRN=1.d-12
       if (err.lt.ERRN*knode*numblk.or.numit.ge.NMI) then
        kend=1
        numit=0
       else
        kend=0
       endif
c      write numit
       NUMUNIT = 48
       ISTATUS = 38
       call openf(48,38,ISTATUS)
       call initrwf(48,iui,iur)
       iui = iui+1
       ipool(iui) = numit
       call endrwf(48,iui,iur)
       do ibk = 1,numblk-1
        call sendint(ibk,nsource,kend)
       enddo
      else
       call sendr(nsource,iblk,err)
       call recvint(iblk,nsource,itemp(1))
       kend=itemp(1)
      endif


      ISTATUS = 18
      call openf(11,18,ISTATUS)
c     write(11)  ((eu(j,i),i=1,knode),j=1,kdgof)
      call initrwf(11,iui,iur)
      do j=1,kdgof
      do i=1,knode
      iur = iur+1
      rpool(iur)=eu(j,i)
      enddo
      enddo
      call endrwf(11,iui,iur)

c
c     renew unod file
c      
      ISTATUS = 33
      call openf(43,33,ISTATUS)
      call initrwf(43,iui,iur)
      do j = 1,kdgof
      do i = 1,knode
      iur = iur+1
      rpool(iur) = eu(j,i)
      end do
      end do     
      do j = 1,kdgof
      do i = 1,knode
      iur = iur+1
      rpool(iur) = ue(j,i)
      end do
      end do
      call endrwf(43,iui,iur)

c...... write munod0 file
        if(iblk.eq.0) then
      ISTATUS = 0
      call openf(20,30,ISTATUS)
c       write(20) ((eu(j,i),i=1,knode),j=1,kdgof)
      call initrwf(20,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do j=1,kdgof
      do i=1,knode
      iur = iur+1
      rpool(iur)=eu(j,i)
      enddo
      enddo
      call endrwf(20,iui,iur)
           if (istop.ge.1.and.kend.eq.1) then
           write(20) ((eu1(j,i),i=1,knode),j=1,kdgof)
           write(20) ((ev1(j,i),i=1,knode),j=1,kdgof)
           endif
        else
      ISTATUS = 0
      call openf(20,30,ISTATUS)
c       write(20) ((eu(j,i),i=1,knode),j=1,kdgof)
      call initrwf(20,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do j=1,kdgof
      do i=1,knode
      iur = iur+1
      rpool(iur)=eu(j,i)
      enddo
      enddo
      if (kend.eq.1) then
        call sendar(nsource,iblk,rpool(iur0),1*kdgof*knode)
      endif
      call endrwf(20,iui,iur)

           if (istop.ge.1.and.kend.eq.1) then
           ISTATUS = 30
           call openf(20,30,ISTATUS)
           call initrwf(20,iui,iur)
           iui0=iui+1
           iur0=iur+1
           do j=1,kdgof
           do i=1,knode
           iur = iur+1
           rpool(iur)=eu1(j,i)
           enddo
           enddo
           call sendar(nsource,iblk,rpool(iur0),1*kdgof*knode)
           call endrwf(20,iui,iur)
           ISTATUS = 30
           call openf(20,30,ISTATUS)
           call initrwf(20,iui,iur)
           iui0=iui+1
           iur0=iur+1
           do j=1,kdgof
           do i=1,knode
           iur = iur+1
           rpool(iur)=ev1(j,i)
           enddo
           enddo
           call sendar(nsource,iblk,rpool(iur0),1*kdgof*knode)
           call endrwf(20,iui,iur)
           endif


        end if

        
1600    format (1x,i4,6e12.5)
        return
        end
