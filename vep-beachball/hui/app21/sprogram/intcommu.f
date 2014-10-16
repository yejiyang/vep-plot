c******************************************************************************c
c FILE: intcommu.f                                                             c
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
       subroutine sendint (iblk,nsource,intsend)

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
       write(*,*) 'Start send one integer message................'
       write(*,*) 'This integer is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This integer is',intsend 
       end if
c
       call MPI_SEND(intsend,1,MPI_INTEGER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'integer one variable success!!!'
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



       subroutine sendai(iblk,nsource,iarray_int,nal)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  iarray_int(nal)
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
       write(*,*) 'Start send one integer array message................'
       write(*,*) 'This integer is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This integer array longth is :',nal
       write(*,9000) (iarray_int(i),i = 1,nal) 
       end if
c
       call MPI_SEND(nal,1,MPI_INTEGER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
       call MPI_SEND(iarray_int,nal,MPI_INTEGER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'integer ',nal,'long array success!!!'
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
9000   format(6i10)
c       call timer(4,2)
       return
       end

c      subroutine send_2d_array_int4 (ndest,nsource,iarray_int4,n1,n2)
       subroutine sendai_2d (iblk,nsource,iarray_int,nrow,ncol)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  iarray_int(ncol*nrow)
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
       write(*,*) 'Start send two dinension integer array message....'
       write(*,*) 'This 2d integer array is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This integer 2d array col and row are :',ncol,nrow
       write(*,9000) ((iarray_int((j-1)*ncol+i),j=1,ncol),i = 1,nrow) 
       end if
c
       call MPI_SEND(nrow,1,MPI_INTEGER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c       
       call MPI_SEND(ncol,1,MPI_INTEGER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
c
       call MPI_SEND(iarray_int,nrow*ncol,MPI_INTEGER,nrecver,
     &               messagetag,
     &               MPI_COMM_WORLD,ierr)
       
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) '2d integer ',nrow,ncol,'array success!!!'
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
9000   format(6i10)
c       call timer(4,2)
       return
       end
       
       subroutine recvint (iblk,nsource,intrecv)

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
       intrecv = 0
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
       write(*,*) 'Start receive one integer message................'
       write(*,*) 'This integer is sent from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
       call MPI_RECV(intrecv,1,MPI_INTEGER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'receiv message from',nsource,' to ', ndest
       write(*,*) 'integer one variable success!!!'
       write(*,*) 'This integer is',intrecv
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

       subroutine recvai(iblk,nsource,iarray_int,nal)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  iarray_int(nal)
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
       iarray_int(i) = 0
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
       write(*,*) 'Start receive integer array message............'
       write(*,*) 'This integer array is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
       call MPI_RECV(nals,1,MPI_INTEGER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
       if(nals.ne.nal) then
       write(*,*) 'Different longth of integer array received'
       write(*,*) 'Fatal error'
       call endjob(ierr)
       end if  
c
       mtimes = mtimes + 1
       if(mtimes.gt.997) mtimes = mod(mtimes,997)
       messagetag = nsource*1000000+ndest*1000+mtimes 
        call MPI_RECV(iarray_int,nals,MPI_INTEGER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
       if( iprint .eq. 1 ) then
       write(*,*) 'receive from',nsource, ' to ', ndest
       write(*,*) 'integer ',nal,'long array success!!!'
       write(*,*) 'and the block number is: ',iblk
       write(*,*) 'the send message tag is :', messagetag
       write(*,*) 'This integer array longth is :',nal
       write(*,9000) (iarray_int(i),i = 1,nal) 
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
9000   format(6i10)
c       call timer(4,2)
       return
       end

       subroutine recvai_2d (iblk,nsource,iarray_int,nrow,ncol)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       dimension  iarray_int(ncol*nrow)
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
       iarray_int(i) = 0
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
       write(*,*) 'Start receive two dinension integer array message...'
       write(*,*) 'This 2d integer array is from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
       call MPI_RECV(nrows,1,MPI_INTEGER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
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
       call MPI_RECV(ncols,1,MPI_INTEGER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c
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
       call MPI_RECV(iarray_int,nrow*ncol,MPI_INTEGER,nsender,
     &               messagetag,
     &               MPI_COMM_WORLD,nstat,ierr)
c       
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) '2d integer ',nrow,ncol,'array success!!!'
       write(*,*) 'and the block number is: ',iblk
       write(*,*) 'the send message tag is :', messagetag
       write(*,*) 'This integer 2d array col and row are :',ncol,nrow
       write(*,9000) ((iarray_int((j-1)*ncol+i),j=1,ncol),i = 1,nrow) 
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
9000   format(6i10)
c       call timer(4,2)
       return
       end
       
       
