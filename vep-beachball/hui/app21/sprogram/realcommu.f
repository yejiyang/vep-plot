c******************************************************************************c
c FILE: realcommu.f                                                            c
c Author: Huai Zhang                                                           c
c DESCRIPTION:                                                                 c
c  Library for lagrange multiplier domain discomposition method                c
c  Called by communication  subroutines. For all the integer communication of  c
c  single integer, integer array and two dimension interarray. A user level    c
c  communication protocol is used here to guarantee the communication safty    c
c  among all the nodes.                                                        c
c******************************************************************************c
c      the arguments in the commom block for mseeage passing  are as following:c
c      ns -- nsource    the source of the message (the rank of the sender)     c
c      nd -- ndest      the destination of the message (the rank of the sender)c 
c      nw -- nworkers   the number of slavers (the total number of slave nodes)c
c      nb -- nblocks    the total number of subdomains                         c
c      mstag -- messagetag the message tag                                     c
c      msb -- the begin message tag for each communication                     c
c      mse -- the end message tag for each communication                       c
c      mst -- the times of message communication                               c
c      msid -- the message id                                                  c
c******************************************************************************c
       subroutine sendr (iblk,nsource,rsend)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c      
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when mpisend...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when mpisend...'
       call endjob(ierr)
       end if
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(nd.eq.0) then
       mtimes = msidm(1)
       else 
       mtimes = msids(1,iblk)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start send one real message................'
       write(*,*) 'This real message is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This real is',rsend 
       end if
c
       call MPI_SEND(rsend,1,MPI_DOUBLE_PRECISION,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr) 
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'double precision one variable success!!!'
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(nd.eq.0) then
       msidm(1) = mtimes
       else
       msids(1,iblk) = mtimes
       endif
c       call timer(4,2)
       return
       end

       subroutine sendar (iblk,nsource,array_r,nal)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  array_r(nal)
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when mpisend...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when mpisend...'
       call endjob(ierr)
       end if
       if((nal.lt.0).or.(nal.gt.100000000)) then
       write(*,*) 'Fatal error, array longth error when mpisend...'
       call endjob(ierr)
       end if
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(nd.eq.0) then
       mtimes = msidm(1)
       else 
       mtimes = msids(1,iblk)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start send one double precision message...........'
       write(*,*) 'This real message is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This longth of this real array is',nal
       write(*,9000) (array_r(i),i=1,nal) 
       end if
c
       rnal = nal    
       call MPI_SEND(rnal,1,MPI_DOUBLE_PRECISION,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes
c
       call MPI_SEND(array_r,nal,MPI_DOUBLE_PRECISION,nrecver,
     &                messagetag
     &               ,MPI_COMM_WORLD,ierr) 
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'double precision one dimension array success!!!'
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(nd.eq.0) then
       msidm(1) = mtimes
       else
       msids(1,iblk) = mtimes
       endif
9000   format(6e16.5) 
c       call timer(4,2)      
       return
       end



c      subroutine send_2d_array_r4 (ndest,nsource,array_r4,n1,n2)
       subroutine sendar_2d (iblk,nsource,array_r,nrow,ncol)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  array_r(ncol*nrow)
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when mpisend...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when mpisend...'
       call endjob(ierr)
       end if
       if((nrow.lt.0).or.(nrow.gt.100000000)) then
       write(*,*) 'Fatal error, row longth error when mpisend...'
       call endjob(ierr)
       end if
       if((ncol.lt.0).or.(ncol.gt.100000000)) then
       write(*,*) 'Fatal error, colum longth error when mpisend...'
       call endjob(ierr)
       end if
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(nd.eq.0) then
       mtimes = msidm(1)
       else 
       mtimes = msids(1,iblk)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start send two dinension real array message....'
       write(*,*) 'This 2d double precision from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This real 2d array col and row are :',ncol,nrow
       write(*,9000) ((array_r((j-1)*ncol+i),j=1,ncol),i = 1,nrow) 
       end if
c
       rnrow = nrow
       call MPI_SEND(rnrow,1,MPI_DOUBLE_PRECISION,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c      
       rncol = ncol 
       call MPI_SEND(rncol,1,MPI_DOUBLE_PRECISION,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       call MPI_SEND(array_r,ncol*nrow,MPI_DOUBLE_PRECISION,nrecver,
     &                messagetag
     &               ,MPI_COMM_WORLD,ierr) 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) '2d real ',nrow,ncol,'array success!!!'
       write(*,*) 'and the block number is: ',nblock
       write(*,*) 'the send message tag is :', messagetag
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(nd.eq.0) then
       msidm(1) = mtimes
       else
       msids(1,iblk) = mtimes
       endif
9000   format(6e15.5)
c       call timer(4,2)
       return
       end
       

       subroutine recvr (iblk,nsource,rrecv)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when recving...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when recving...'
       call endjob(ierr)
       end if
c
c      cleaning buffer
c       
       rrecv = 0.d0
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(ns.eq.0) then
       mtimes = msidm(2)
       else 
       mtimes = msids(2,ns)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start receive one real message................'
       write(*,*) 'This real message is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
        call MPI_RECV(rrecv,1,MPI_DOUBLE_PRECISION,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,IERR)
 
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'double precision one variable success!!!'
       write(*,*) 'This real is',rrecv
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(ns.eq.0) then
       msidm(2) = mtimes
       else
       msids(2,ns) = mtimes
       endif
c       call timer(4,2)
       return
       end


       subroutine recvar (iblk,nsource,array_r,nal)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  array_r(nal)
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when recving...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when recving...'
       call endjob(ierr)
       end if
       if((nal.lt.0).or.(nal.gt.100000000)) then
       write(*,*) 'Fatal error, array longth error when recving...'
       call endjob(ierr)
       end if
c
c      cleaning buffer
c       
       do i =1,nal
       array_r(i) = 0.d0
       end do   
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(ns.eq.0) then
       mtimes = msidm(2)
       else 
       mtimes = msids(2,ns)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start receive one double precision message.........'
       write(*,*) 'This real message is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
       call MPI_RECV(rnals,1,MPI_DOUBLE_PRECISION,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       nals = int(rnals)   
       if(nals.ne.nal) then
       write(*,*) 'Different longth of double precision array received'
       write(*,*) 'Fatal error'
       call endjob(ierr)
       end if  
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes
c
        call MPI_RECV(array_r,nal,MPI_DOUBLE_PRECISION,nsender,
     &               messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       if( iprint .eq. 1 ) then
       write(*,*) 'receive from',nsource, ' to ', ndest
       write(*,*) 'double precision one dimension array success!!!'
       write(*,*) 'This longth of this real array is',nal
       write(*,9000) (array_r(i),i=1,nal) 
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(ns.eq.0) then
       msidm(2) = mtimes
       else
       msids(2,ns) = mtimes
       endif
9000   format(6e16.5)       
c       call timer(4,2)
       return
       end

c      subroutine recv_2d_array_r4 (ndest,nsource,array_r4,n1,n2)
       subroutine recvar_2d (iblk,nsource,array_r,nrow,ncol)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  array_r(ncol*nrow)
       include 'mpif.h'
       dimension nstat(MPI_STATUS_SIZE)
c
c      check parameters 
c
c       call timer(4,1)
       if(iblk.lt.0) then 
       write(*,*) 'Fatal error, destination id error when recving...'
       call endjob(ierr)
       end if
       if(nsource.lt.0) then
       write(*,*) 'Fatal error, source id error when recving...'
       call endjob(ierr)
       end if
       if((nrow.lt.0).or.(nrow.gt.100000000)) then
       write(*,*) 'Fatal error, row longth error when recving...'
       call endjob(ierr)
       end if
       if((ncol.lt.0).or.(ncol.gt.100000000)) then
       write(*,*) 'Fatal error, colum longth error when recving...'
       call endjob(ierr)
       end if
c
c      cleaning buffer 
c
       do i = 1,ncol*nrow
       array_r(i) = 0.d0
       end do      
c
c      get message tag first
c
       ns = nsource
       nd = iblk
       ndest = iblk
       if(ns.eq.0) then
       mtimes = msidm(2)
       else 
       mtimes = msids(2,ns)
       endif
       nworkers = numcpus - 1
       nsender = nsource
       nrecver = iblk 
       if(nsource.gt.nworkers) nsender = mod(nsource,nworkers)
       if(iblk.gt.nworkers)    nrecver = mod(iblk,nworkers)
c
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       if( iprint .eq. 1 ) then
       write(*,*) 'Start send two dinension real array message....'
       write(*,*) 'This 2d double precision from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
       call MPI_RECV(rnrows,1,MPI_DOUBLE_PRECISION,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       nrows = int(rnrows)
       if(nrows.ne.nrow) then
       write(*,*) 'Different row number of integer array received'
       write(*,*) 'Fatal error'
       call endjob(ierr)
       end if  
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c       
       call MPI_RECV(rncols,1,MPI_DOUBLE_PRECISION,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       ncols = int(rncols)   
       if(ncols.ne.ncol) then
       write(*,*) 'Different clolum number of integer array received'
       write(*,*) 'Fatal error'
       call endjob(ierr)
       end if  
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       call MPI_RECV(array_r,ncol*nrow,MPI_DOUBLE_PRECISION,nsender,
     &               messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) '2d real ',nrow,ncol,'array success!!!'
       write(*,*) 'and the block number is: ',nblock
       write(*,*) 'the send message tag is :', messagetag
       write(*,*) 'This real 2d array col and row are :',ncol,nrow
       write(*,9000) ((array_r((j-1)*ncol+i),j=1,ncol),i = 1,nrow) 
       end if
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       nd = iblk
       if(ns.eq.0) then
       msidm(2) = mtimes
       else
       msids(2,ns) = mtimes
       endif
9000   format(6e15.5)
c       call timer(4,2)
       return
       end
