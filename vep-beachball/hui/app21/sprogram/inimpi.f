c******************************************************************************c
c FILE: inimpi.f                                                               c
c Author: Huai Zhang                                                           c
c DESCRIPTION:                                                                 c
c  Called by ddm.f code to obtain message passing interface(MPI) environment   c
c  and initialize the user level communication protocol, then return the rank  c
c  of the distributed parallel machine node to mail program.                   c 
c******************************************************************************c
      subroutine inimpi(myrank,numblk)

      implicit real*8 (a-h,o-z)
      common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
      include 'mpif.h'
c
c     initial message common block
c
      master = 0
      nsource = 0
      do i=1,2
      msidm(i)=0
      enddo
      do i=1,514
      do j = 1,2
      msids(j,i) = 0
      end do
      end do
      iprint = 0 
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numnodes, ierr )
      numcpus = numnodes
      myrank = myid
      if (myid .eq. master) then
      write(*,*) '**********  Starting MPI Master Task  ************'
c
c     get the lmddm control file
c
      open(1,file='partition.dat',form='formatted',status='unknown')
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) numblk,numtyp,nmdof,numtypl,kdgofl,keymt,lgio
      read(1,*)
      read(1,*) t0,tmax,dt
      close(1)
c      write(*,*) numcpus,' CPUs compute ',numblk,' subdomains.'//
c     +           ' Node ',myid,' is assigned as master node.'//
c     +           ' and the other nodes are all assigned as slave nodes.'
c
c     send out numblk parameter 
c
      do iblk = 1,numblk
      call sendint(iblk,nsource,numblk)
      end do
c
      end if  ! master node
c
c     slave node task 
c
      if(myid.gt.master) then
c
c     receive numblk from master node
c
      call recvint(myid,nsource,numblk)
c
c     data check for node balancing 
c
      ndomains = numblk
      if((myid).gt.numblk) then  
      write(*,*) 'There are more CPUs than blocknumbers of subdomains'
      write(*,*) 'Fatal error ............'
      call endjob(ierr)
      end if
      if(myid.gt.512) then
      write(*,*) ' CPUs number greater then 512, stop......'
      call endjob(ierr)
      end if
      if(numblk.gt.512) then
      write(*,*) ' subdomain number greater then 512, stop......'
      call endjob(ierr)
      end if 
c
c     Current version must has condition as following:
c
      if (numcpus.ne.(numblk)) then
      write(*,*)' There are less CPUs to compute all the subdomains...'
      call endjob(ierr)
      end if
c
      end if ! slave node
c  
      return
      end

      subroutine findworkers(numworkers)
      implicit real*8 (a-h,o-z)
      integer  ierr,numworkers
      include 'mpif.h'
c
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numcpu, ierr )
      numworkers = numcpu-1
      if (numworkers.le.0) then
      print *,'Number_workers less than Minworkers,fatal errer!!'
      call endjob(ierr)
      end if
      return
      end

      subroutine findrunk(numrunk)
      implicit real*8 (a-h,o-z)
      integer  ierr,numrunk
      include 'mpif.h'
c
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrunk, ierr )
      numrunk = myrunk
      if (numrunk.le.0) then
      print *,'Number_runk less than Minworkers,fatal errer!!'
      call endjob(ierr)
      end if
      return
      end

      subroutine endjob(ierr)
      implicit real*8 (a-h,o-z)
      integer  ierr,numrunk
      include 'mpif.h'
c
      call MPI_FINALIZE(ierr)
      print *,'The process has been successfully finalized !!'
      return
      end


        subroutine mendit(kend)
        implicit real*8 (a-h,o-z)
        logical filflg
 
        open(10,file='lmddm',form='unformatted')
        read(10) numblk
        close(10)
        inquire(file='endit',exist=filflg)
        kend=1
        if (filflg) kend=0
        nsource=0
        do iblk=1,numblk
        call sendint(iblk,nsource,kend)
        enddo
        return
        end
 
        subroutine ssendit(kend,iblk)
        implicit real*8 (a-h,o-z)
 
        nsource=0
        call recvint(iblk,nsource,kend)
        write(*,*) 'kend iblk ==',kend,iblk
        return
        end
 
        subroutine mend(kend)
        implicit real*8 (a-h,o-z)
        logical filflg
 
        kend=1
        inquire(file='end',exist=filflg)
        if (filflg) kend=0
        write(*,*) 'kend ==',kend
        open(10,file='lmddm',form='unformatted',status='old')
        read(10) numblk
        close(10)
        nsource=0
        do iblk=1,numblk
        call sendint(iblk,nsource,kend)
        enddo
 
        return
        end
 
        subroutine ssend(kend,iblk)
        implicit real*8 (a-h,o-z)
c 
        nsource=0
        call recvint(iblk,nsource,kend)
c 
        return
        end
