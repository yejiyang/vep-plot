      subroutine Mnzrecvdip(istop)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
 
        nsource=0
c...... open lmddm file
c        open(21,file='mlmddm',form='unformatted')
c        read(21) numblk,numtyp
c        close(21)

      ISTATUS = 1
      call openf(21,1,ISTATUS)
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
      nskip=ipool(iui)
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
      iui = iui+1
      nskip=ipool(iui)
      iui = iui+1
      knode=ipool(iui)
      call endrwf(21,iui,iur)

c        write(*,*) 'numblk ===',numblk
c...... open munod0 file
2        open(22,file='munod0',form='unformatted',status='unknown')
c       open(22,file='munod0',form='unformatted',status='unknown')
c...... open munod file
        open(23,file='munod',form='unformatted',status='unknown')
             open(123,file='u1v1',form='unformatted',status='unknown')
       open(133,file='disp-vector1',form='formatted',status='unknown',
     *position='append')
c       open(134,file='disp-vector2',form='formatted',status='unknown',
c     *position='append')
c       open(135,file='disp-vector3',form='formatted',status='unknown',
c     *position='append')
c...... open pdisp0 file
c        open(24,file='mdisp0',form='unformatted',status='unknown')
c...... open petype0 file
c        open(25,file='metype0',form='unformatted',status='unknown')
c...... open disp0 file
c        open(26,file='disp0',form='unformatted',status='unknown')
 
        kdgof=3

      ISTATUS = 35
      call openf(35,6,ISTATUS)

        do 1000 iblk=0,numblk-1
c          read(24) knodei,kdgof
c          print *,'In recvdip :: knodei == ',knodei
c          read(24) (xskp,i=1,knodei*kdgof)
 
      kna1=kdgof*knode*1
      knb2=numtyp*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      knb3=numtyp*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      knb4=numtyp*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb5=numtyp*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      knb1=numtyp*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb6=knode*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      kna2=kdgof*knode*1
      kna3=kdgof*knode*1
      kna4=kdgof*knode*1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      knb3=knb3+knb2
      knb4=knb4+knb3
      knb5=knb5+knb4
      knb6=knb6+knb5
      call nzrecvdip(kdgof,knode,iblk,numblk,numtyp,
     *nsource,aa(kna0),aa(kna1),ia(knb0),ia(knb1),
     *ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *filename,istop,aa(kna2),aa(kna3))
 
1000    continue
7        close(22)
        close(23)
             close(123)
             close(133)
c             close(134)
c             close(135)
c        close(24)
c      close(25)
c        close(26)

      end

      subroutine nzrecvdip(kdgof,knode,iblk,numblk,numtyp,
     *nsource,u,ui,idet,numa,nnea,mmta,nmta,
     *inode,
     *filename,istop,u1i,v1i)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
        dimension u(kdgof,knode),numa(numtyp),nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),inode(knode),ui(kdgof,knode)
        dimension u1i(kdgof,knode),v1i(kdgof,knode)
        dimension nleft(26),nright(26)
c
c       read etype0 file
c
c
c        READ (25) NUMTYP,NODALL,MATALL
c      READ (25) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
c     &nskip,(inode(i),i=1,knodei)
c      read(25) nliq,(nskip,i=1,nliq)

      call initrwf(35,iui,iur)
      iui = iui+1
      numtyp = ipool(iui)
      iui = iui+1
      nodall = ipool(iui)
      iui = iui+1
      matall = ipool(iui)
      iui = iui+1
      knodei = ipool(iui)
      iui = iui+1
      nelem = ipool(iui)
      do i=1,numtyp
      iui = iui+1
      idet(i) = ipool(iui)
      iui = iui+1
      numa(i) = ipool(iui)
      iui = iui+1
      nnea(i) = ipool(iui)
      iui = iui+1
      mmta(i) = ipool(iui)
      iui = iui+1
      nmta(i) = ipool(iui)
      enddo
      do i=1,knodei
      iui = iui+1
      inode(i) = ipool(iui)
      enddo
      do i=1,nelem
      iui = iui+1
      nskip = ipool(iui)
      enddo
      call endrwf(35,iui,iur)


c
c       read munod0 file
c
        if(iblk.eq.0) then
c        read(22) ((ui(j,i),i=1,knodei),j=1,kdgof)
      ISTATUS = 30
      call openf(40,30,ISTATUS)
      call initrwf(40,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do j=1,kdgof
      do i=1,knodei
      iur = iur+1
      ui(j,i)=rpool(iur)
      enddo
      enddo
      call endrwf(40,iui,iur)
           if (istop.ge.1) then
           read(22) ((u1i(j,i),i=1,knodei),j=1,kdgof)
           read(22) ((v1i(j,i),i=1,knodei),j=1,kdgof)
           endif
        else
c       read(22) ((ui(j,i),i=1,knodei),j=1,kdgof)
      ISTATUS = 0
      call openf(44,34,ISTATUS)
      call initrwf(44,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(nsource,iblk,rpool(iur0),1*kdgof*knodei)
      do j=1,kdgof
      do i=1,knodei
      iur = iur+1
      ui(j,i)=rpool(iur)
      enddo
      enddo
      call endrwf(44,iui,iur)
             if (istop.ge.1) then
             ISTATUS = 34
             call openf(44,34,ISTATUS)
             call initrwf(44,iui,iur)
             iui0=iui+1
             iur0=iur+1
             call recvar(nsource,iblk,rpool(iur0),1*kdgof*knodei)
             do j=1,kdgof
             do i=1,knodei
             iur = iur+1
             u1i(j,i)=rpool(iur)
             enddo
             enddo
             call endrwf(44,iui,iur)
             ISTATUS = 34
             call openf(44,34,ISTATUS)
             call initrwf(44,iui,iur)
             iui0=iui+1
             iur0=iur+1
             call recvar(nsource,iblk,rpool(iur0),1*kdgof*knodei)
             do j=1,kdgof
             do i=1,knodei
             iur = iur+1
             v1i(j,i)=rpool(iur)
             enddo
             enddo
             call endrwf(44,iui,iur)
             endif
        end if
c
c       write munod file
c
        do i = 1,knodei
      ii = inode(i)
      do j = 1,kdgof
      u(j,ii) = ui(j,i)
      end do
      end do
c
             if (istop.ge.1) then
             write(123) ((u1i(j,i),i=1,knodei),j=1,kdgof)
             write(123) ((v1i(j,i),i=1,knodei),j=1,kdgof)
             endif


        if(iblk.eq.(numblk-1)) then
        write(23) ((u(j,i),i=1,knode),j=1,kdgof)

c        nx=25
c        ny=24
c        nz=6
c        do i=1,25
c        nleft(i)=(nx+1)*(ny+1)*nz+(nx+1)*(i-1)+(nx+1)/2
c        nright(i)=(nx+1)*(ny+1)*nz+(nx+1)*(i-1)+(nx+1)/2+1
c        enddo
        np2=1978

        nleft(1)  =1422+np2*5
        nleft(2)  =1311+np2*5

        nright(1) =1420+np2*5
        nright(2) =1310+np2*5
        
        write(133,*) '---------------------------------'
        write(133,'(40e11.3)') (u(1,nleft(i))-u(1,nright(i)),
     *u(2,nleft(i))-u(2,nright(i)),i=1,2)

c        do i=1,19
c        du=u0(1,nleft(i))-u0(1,nright(i))
c        dv=u0(2,nleft(i))-u0(2,nright(i))
c        duv=sqrt(du*du+dv*dv)
c        write(*,1100) i,du,dv,duv
c        enddo
        
        nleft(1) = 236    +np2*5
        nleft(2) = 261    +np2*5
        
        nright(1) = 182   +np2*5
        nright(2) = 208   +np2*5
               
        write(*,*)
        write(133,'(40e11.3)') (u(1,nleft(i))-u(1,nright(i)),
     *u(2,nleft(i))-u(2,nright(i)),i=1,2)

        end if
        return
        end
