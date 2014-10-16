      subroutine Mnzrecvstr(istop)
      implicit real*8 (a-h,o-z)
      character*12 fname,filename(20)
      include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
      common /aa/ aa(maxaa)
      common /ia/ ia(maxia)
        nsource = 0
        iblk = 1
c...... open lmddm file
c        open(21,file='mlmddm',form='unformatted')
c        read(21) numblk,numtyp,maxlm,nmdof,numtypl,kdgofl,keymt
c        read(21) lgio,t0,tmax,dt
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
c.....open time file
      ISTATUS = 11
      call openf(21,11,ISTATUS)
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

c...... open munods0 file
1        open(22,file='munods0',form='unformatted',status='unknown')
c       open(22,file='munods0',form='unformatted',status='unknown')
c...... open munods file
        open(23,file='munod1',form='unformatted',status='unknown')
        open(124,file='fstr1',form='unformatted',status='unknown')
c...... open block number for each node file
c      open(126,file='nmbnodes0',form='unformatted',status='unknown')
      open(127,file='curve1',form='formatted',status='unknown',
     *position='append')
c      open(128,file='curve2',form='formatted',status='unknown',
c     *position='append')
c      open(129,file='curve3',form='formatted',status='unknown',
c     *position='append')
c      open(130,file='curve4',form='formatted',status='unknown',
c     *position='append')
c      open(131,file='curve5',form='formatted',status='unknown',
c     *position='append')
c      open(132,file='curve6',form='formatted',status='unknown',
c     *position='append')
c...... open pdisp0 file
c        open(24,file='mdisp0',form='unformatted',status='unknown')
c...... open petype0 file
c        open(25,file='metype0',form='unformatted',status='unknown')
c...... open disp0 file
c        open(26,file='disp0',form='unformatted',status='unknown')
c      read(26) knode,nskip

      ISTATUS = 35
      call openf(35,6,ISTATUS)

        do 1000 iblk=0,numblk-1
          if(iblk.eq.0) then
c          read(22) knodei,kdgofs
      ISTATUS = 31
      call openf(26,31,ISTATUS)
      call initrwf(26,iui,iur)
      iui = iui+1
      knodei=ipool(iui)
      iui = iui+1
      kdgofs=ipool(iui)
      call endrwf(26,iui,iur)
          else
c         read(22) knodei,kdgofs
      ISTATUS = 0
      call openf(45,35,ISTATUS)
      call initrwf(45,iui,iur)
      iui = iui+1
      call recvint(nsource,iblk,ipool(iui))
      knodei=ipool(iui)
      iui = iui+1
      call recvint(nsource,iblk,ipool(iui))
      kdgofs=ipool(iui)
      call endrwf(45,iui,iur)
          end if
c          read(24) knodei,kdgof
c          read(24) (xskp,i=1,knodei*kdgof)
 
      kna1=kdgofs*knode*1
      kna2=kdgofs*knodei*1
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
      knb6=knodei*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb7=knode*1
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
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
      call nzrecvstr(kdgofs,knode,iblk,numblk,numtyp,knodei,
     *nsource,aa(kna0),aa(kna1),ia(knb0),ia(knb1),
     *ia(knb2),ia(knb3),ia(knb4),ia(knb5),
     *filename,istop,ia(knb6),dt,time,it)
 
1000    continue
3        close(22)
        close(23)
        close(124)
c        close(126)
        close(127)
c        close(128)
c        close(129)
c        close(130)
c        close(131)
c        close(132)
c        close(24)
c        close(25)
c        close(26)
      end

      subroutine nzrecvstr(kdgofs,knode,iblk,numblk,numtyp,knodei,
     *nsource,us,usi,idet,numa,nnea,mmta,nmta,
     *inode,
     *filename,istop,idnode,dt,time,it)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpoolm),ipool(maxrpoolm)
        dimension us(kdgofs,knode),usi(kdgofs,knodei),numa(numtyp),
     &  nnea(numtyp),mmta(numtyp),
     &  nmta(numtyp),idet(numtyp),inode(knodei)
        dimension idnode(knode)
        dimension nleft(26),nright(26)
c
c       read etype0 file
c
c
c        READ (25) NUMTYP,NODALL,MATALL
c      READ (25) (idet(i),numa(i),nnea(i),mmta(i),nmta(i),i=1,numtyp),
c     & nskip,(inode(i),i=1,knodei)
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

c      read(126) (idnode(i),i=1,knode)
      NUMUNIT = 50
      ISTATUS = 50
      call openf(50,40,ISTATUS)
      call initrwf(50,iui,iur)
      do i = 1,knode
      iui = iui+1
      idnode(i) = ipool(iui)
      enddo
      call endrwf(50,iui,iur)


      num=numa(1)

c
c       read munods0 file
c
        if(iblk.eq.0) then
c        read(22) ((usi(j,i),i=1,knodei),j=1,kdgofs)
      call initrwf(26,iui,iur)
      iui0=iui+1
      iur0=iur+1
      do J=1,KDGOFS
      do I=1,KNODEI
      iur = iur+1
      usi(j,i)=rpool(iur)
      enddo
      enddo
             if (istop.ge.1) then
             write(*,*)'hellox, fstr'
             write(124) (rpool(iur0+i-1),i=1,kdgofs*knodei)
             endif
      call endrwf(26,iui,iur)
             if (istop.ge.1) then
             ISTATUS = 0
             call openf(45,35,ISTATUS)
             call initrwf(45,iui,iur)
             iui0=iui+1
             iur0=iur+1
             read(22) (rpool(iur0+i-1),i=1,num*6*8)
             write(124) (rpool(iur0+i-1),i=1,num*6*8)
             call endrwf(45,iui,iur)
             endif
        else
c       read(22) ((usi(j,i),i=1,knodei),j=1,kdgofs)
      ISTATUS = 35
      call openf(45,35,ISTATUS)
      call initrwf(45,iui,iur)
      iui0=iui+1
      iur0=iur+1
      call recvar(nsource,iblk,rpool(iur0),1*kdgofs*knodei)
      do j=1,kdgofs
      do i=1,knodei
      iur = iur+1
      usi(j,i)=rpool(iur)
      enddo
      enddo
             if (istop.ge.1) then
             write(124) (rpool(iur0+i-1),i=1,kdgofs*knodei)
             endif
      call endrwf(45,iui,iur)
             if (istop.ge.1) then
             ISTATUS = 35
             call openf(45,35,ISTATUS)
             call initrwf(45,iui,iur)
             iui0=iui+1
             iur0=iur+1
             call recvar(nsource,iblk,rpool(iur0),num*6*8)
             write(124) (rpool(iur0+i-1),i=1,num*6*8)
             call endrwf(45,iui,iur)
             endif
        end if
 
c
c
c       write munods file
c
c
      do i = 1,knodei
        ii = inode(i)
        if (idnode(ii).eq.iblk) then
          do j = 1,kdgofs
          us(j,ii) = usi(j,i)
          end do
        endif
      end do
c
        if(iblk.eq.(numblk-1)) then
          write(23) ((us(j,i),i=1,knode),j=1,kdgofs)

c        nx=25
c        ny=24
c        nz=6
c        do i=1,25
c        nleft(i)=(nx+1)*(ny+1)*nz+(nx+1)*(i-1)+(nx+1)/2
c        nright(i)=(nx+1)*(ny+1)*nz+(nx+1)*(i-1)+(nx+1)/2+1
c        enddo

        nleft(1) =1422+1978*5
        nleft(2) =1311+1978*5

        nright(1) =1420+1978*5
        nright(2) =1310+1978*5

        write(127,*) 'it=',it
        write(127,'(30e11.4)') ((us(4,nleft(i))+us(4,nright(i)))/2,
     *i=1,2)
c        write(128,'(10e11.4)') ((us(4,nleft(i))+us(4,nright(i)))/2,
c     *i=9,16)
c        write(129,'(10e11.4)') ((us(4,nleft(i))+us(4,nright(i)))/2,
c     *i=17,25)
c        write(130,'(10e11.4)') ((us(4,1992+i)+us(4,2039+i))/2,i=1,8)
c        write(131,'(10e11.4)') ((us(4,2000+i)+us(4,2047+i))/2,i=1,8)
c        write(132,'(10e11.4)') ((us(4,2008+i)+us(4,2055+i))/2,i=1,7)
        end if
        return
        end
