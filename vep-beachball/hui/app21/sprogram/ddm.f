C**************************************************************************C
C              & & &     & & & &    & & & &    & & &      & & & &          C
C              &    &    &          &          &    &    &                 C
C              & & &     & & & &    & & & &    & & &     &     & &         C
C              &         &          &          &         &      &          C
C              &         &          & & & &    &           && & &          C
C                                                                          C
C     1026 - 2002 (C) COPYRIGHT INISTITUTE OF MATHEMATICS,                 C
C                               CHINESE ACADEMY OF SCIENCES.               C
C                     All Rights Reserved                                  C
C                     AUTHORS:                                             C
C                       Guoping Liang                                      C
C                       Email:  guoping@fegensoft.com                      C
C                               ling@fegensoft.com                         C
C                       Huai Zhang                                         C
C                       Email:  huaizhang@hotmail.com                      C
C                               hzhang@mail.amss.ac.cn                     C
C                       Shaopeng Liu                                       C
C                       Email:                                             C
C                       Date :  Oct. 2002                                  C
C**************************************************************************C
C -------------------------------------------------------------------------C
C Ddm.f - Parallelized Fortran Version                                     C
C FILE: Ddm.f                                                              C
C Ddm.f: main program, based on Lagrange Multiplier Domain Decomposition   C
C        Method. (LMDDM)                                                   C
C OTHER FILES:                                                             C
C                                                                          C
C fepglib.f: lib subroutines for fepg                                      C
C incore.f:  memmory file administrate subroutines                         C
C ddmm.f:    master program for LMDDM algorithm                            C
C ddms.f:    slave  program for LMDDM algorithm                            C
C partition.f: graph partition and data partition based on Metis packages  C
C sendpart.f: send all the initial data from master node to slave node     C
C recvpart.f: slave receive data for master node                           C
C order.f: reordering the node number of unstructured meshes               C
C start.f: initial node id and initial values and boundary values          C
C ediag.f: caculate the one-dimension compressed gross stiff matrix index  C
C bft.f: boundary values and initial values update for each time step      C
C srecv.f: slave node receive some data from master node during computing  C
C e*.f: element caculate subroutines for all types of meshes               C
C sinsub.f/ninsub.f: subroutines for caculate A(-1)                        C
C mlmm.f: to caculate D metrix, B(t)A(-1)B lamada = B(t)A(-1)f             C
C u*.f: postprocessing subroutines for each subdomain                      C
C intcommu.f: MPI communication subroutines for each Integer data          C
C realcommu.f: MPI communication subroutines for each real data            C
C                                                                          C
C Notice:                                                                  C
C  In this parallelized version, the grid is decomposed by the master      C
C process and then mapped and distributed to workers processes.            C
C -------------------------------------------------------------------------C
C Explanation of constants and variables                                   C
C   MAXWORKER                      =  maximum number of workers tasks      C
C   MINWORKER                      =  minimum number of workers tasks      C
C   NWORKERS                       =  number of workers processes          C
C   NDEST,NSOURCE                  =  to - from for message send-receive   C
C -------------------------------------------------------------------------C
      program ddm
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
c     Initializing the parallel computation encirenment and find out
c     how many nodes this programm can use,if the number of the
c     processes is less than the minium number if processes assigned,
c     then quit........
      master = 0
      maxnodeid = 512
      minnodeid = 2 
      write(*,*) 'initializing the parallel computational environment..'
      call inimpi(myrank,numblk)
      call initpool
c      print *,'myrank=====',myrank
c      print *, '*********** My node rank is:',myrank,'***************'
c      print *, '*********** Starting Slave Nodes Programing *********'
      call slavep(myrank)
c      call timer(9,2)
c      write(*,*)'.......... Slave Program Success !!! ..............'
c      call timer(9,1)
c 
c     time and load imbalance report and analysis 
c
c      call timer(6,1)
c      call timer(7,1)
c      call timer(8,1)
c      call timer(9,2)
c      if (myrank.eq.master) then
c      open(88,file='report',form='formatted', status='unknown')
c      end if
c      ioport = 88
c      isreport = 1
c      call  timereport(myrank,myrank,numblk,ioport,isreport)
c      if (myrank.eq.master) then
c      close(88)
c      end if 
c
      call endjob(ierr)
c
c
      stop
      end
