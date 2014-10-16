c******************************************************************************c
c FILE: charcommu.f                                                            c
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
       subroutine sendchr (iblk,nsource,chr)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       character*1 chr
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
       write(*,*) 'This character from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       write(*,*) 'This character is',chr 
       end if
c
       call MPI_SEND(chr,1,MPI_CHARACTER,nrecver,messagetag,
     &               MPI_COMM_WORLD,ierr)
c      
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'character one variable success!!!'
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

       subroutine recvchr (iblk,nsource,chrrecv)

       implicit real*8 (a-h,o-z)
       common/messg/numcpus,myid,ndomains,msids(2,0:514),msidm(2),iprint
       character*1 chrrecv
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
       chrrecv = ' '
c
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
       write(*,*) 'This character from ',nsource,' to ',iblk
       write(*,*) 'The source computer node is',nsender
       write(*,*) 'The destination computer node is',nrecver
       write(*,*) 'Message tag is:',messagetag
       end if
c
        call MPI_RECV(chrrecv,1,MPI_CHARACTER,nsender,messagetag,
     &               MPI_COMM_WORLD,nstat,IERR)
c
       if( iprint .eq. 1 ) then
       write(*,*) 'sent from',nsource, ' to ', ndest
       write(*,*) 'character one variable success!!!'
       write(*,*) 'This character is',chrrecv 
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

