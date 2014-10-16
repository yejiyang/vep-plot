      subroutine Mazsolv(iblk)
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
c ....open dmsr file
c     open (22,file='dmsr',form='unformatted',status='unknown')
c...... open subnodes file
c     open (23,file='ssubnodes0',form='unformatted',status='unknown')
c...... open nv file
c     open (24,file='nv',form='unformatted',status='unknown')
c...... open lsol file
c     open (25,file='u',form='unformatted',status='unknown')
c
c      do 2300 iblk=1,numblk
      ISTATUS = 26
      call openf(22,26,ISTATUS)
c       read(22) neq,maxa
      call initrwf(22,iui,iur)
      iui = iui+1
      neq=ipool(iui)
      iui = iui+1
      maxa=ipool(iui)
      call endrwf(22,iui,iur)
        neq1=neq+1
      ISTATUS = 8
      call openf(23,8,ISTATUS)
c       read(23) knode,kdgof
      call initrwf(23,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      call endrwf(23,iui,iur)
        neqmax = knode*kdgof

      neqmax1=neqmax+1
      maxa1=maxa+1
 
      knb1=neq1*1
      if (knb1/2*2 .lt. knb1) knb1=knb1+1
      knb2=neq1*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      kna4=neq1*1
      knb15=neqmax1*1
      if (knb15/2*2 .lt. knb15) knb15=knb15+1
      knb3=maxa1*1
      if (knb3/2*2 .lt. knb3) knb3=knb3+1
      kna1=maxa1*1
      kna2=neqmax1*1
      kna3=neqmax1*1
      knb4=neqmax1*1
      if (knb4/2*2 .lt. knb4) knb4=knb4+1
      knb5=neqmax1*1
      if (knb5/2*2 .lt. knb5) knb5=knb5+1
      knb6=neqmax1*1
      if (knb6/2*2 .lt. knb6) knb6=knb6+1
      knb14=neqmax1*1
      if (knb14/2*2 .lt. knb14) knb14=knb14+1
      knb8=neqmax1*1
      if (knb8/2*2 .lt. knb8) knb8=knb8+1
      knb13=neqmax1*1
      if (knb13/2*2 .lt. knb13) knb13=knb13+1
      if (knb16/2*2 .lt. knb16) knb16=knb16+1
      knb9=1*1
      if (knb9/2*2 .lt. knb9) knb9=knb9+1
      knb10=1*1
      if (knb10/2*2 .lt. knb10) knb10=knb10+1
      knb11=1*1
      if (knb11/2*2 .lt. knb11) knb11=knb11+1
      knb12=1*1
      if (knb12/2*2 .lt. knb12) knb12=knb12+1
      knb7=kdgof*knode*1/numblk*3
      knodei0 = knode/numblk*3
      if (knb7/2*2 .lt. knb7) knb7=knb7+1
      kna5=kdgof*knode*1
      kna6=neqmax1*1
      kna0=1
      kna1=kna1+kna0
      kna2=kna2+kna1
      kna3=kna3+kna2
      kna4=kna4+kna3
      kna5=kna5+kna4
      kna6=kna6+kna5
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
      knb15=knb15+knb14
      call azsolv(knodei0,iblk,numblk,maxa,neq,neq1,nsource,
     *neqmax,knode,kdgof,neqmax1,maxa1,aa(kna0),aa(kna1),aa(kna2),
     *aa(kna3),aa(kna4),aa(kna5),ia(knb0),ia(knb1),
     *ia(knb2),ia(knb3),ia(knb4),ia(knb5),ia(knb6),
     *ia(knb7),ia(knb8),ia(knb9),ia(knb10),ia(knb11),
     *ia(knb12),ia(knb13),ia(knb14),
     *filename)
 
2300  continue
 
c     close(21)
c     close(22)
c     close(23)
c     close(24)
c     close(25)
      end
      subroutine azsolv(knodei0,iblk,numblk,maxa,neq,neq1,nsource,
     *neqmax,knode,kdgof,neqmax1,maxa1,val,x,b,diag,uu,
     *x1,jdiag,jdiagaz,nbindex,nupdate,ndata_org,nexternal,nodvar,
     *jupdate,idnode,idnode1,inode,jnode,nupdate_index,nextern_index,
     *mapgl,
     *filename)
      implicit real*8 (a-h,o-z)
      character*12 filename(20)
        include 'memalloc.h'
      common /pool/ rpool(maxrpools),ipool(maxrpools)
      dimension jdiag(neq1),jdiagaz(neq1),diag(neq1),mapgl(neqmax1),
     &     nbindex(maxa1),val(maxa1),x(neqmax1),b(neqmax1),
     &     nupdate(neqmax1),
     &     ndata_org(neqmax1),nexternal(neqmax1),
     &     nextern_index(neqmax1),jupdate(neqmax1),
     &     nupdate_index(neqmax1),
     &     idnode(1),idnode1(1),inode(1),jnode(1),
     &     nodvar(kdgof,knodei0),uu(kdgof,knode),x1(neqmax1)
c
c
6       format (1x,10i5)
7       format (1x,6e12.5)
c
c
c       read DMSR format gross matrix and index information
c
c       read(22) (JDIAG(I),I=1,NEQ)
      call initrwf(22,iui,iur)
      do I=1,NEQ
      iui = iui+1
      jdiag(I)=ipool(iui)
      enddo
      call endrwf(22,iui,iur)
c       read(22) (JDIAGAZ(I),I=1,NEQ1)
      call initrwf(22,iui,iur)
      do I=1,NEQ1
      iui = iui+1
      jdiagaz(I)=ipool(iui)
      enddo
      call endrwf(22,iui,iur)
c       read(22) (nbindex(I),I=1,MAXA)
      call initrwf(22,iui,iur)
      do I=1,MAXA
      iui = iui+1
      nbindex(I)=ipool(iui)
      enddo
      call endrwf(22,iui,iur)
c       read(22) (val(I),I=1,MAXA)
      call initrwf(22,iui,iur)
      do I=1,MAXA
      iur = iur+1
      val(I)=rpool(iur)
      enddo
      call endrwf(22,iui,iur)
c       read(22) (b(I),I=1,NEQ)
      call initrwf(22,iui,iur)
      do I=1,NEQ
      iur = iur+1
      b(I)=rpool(iur)
      enddo
      call endrwf(22,iui,iur)
c      if(iblk.eq.0) then 
c
c       print out the dmsr matrix for error check
c
c        write(*,*) 'neq, new maxa=====',neq,maxa
c        write(*,*) 'new jdiag file ====='
c        write(*,*) (JDIAG(I),I=1,NEQ)
c        write(*,*) 'new jdiagaz file ====='
c        write(*,*) (JDIAGAZ(I),I=1,NEQ1)
c        write(*,*) 'new na file ====='
c        write(*,*) (nbindex(I),I=1,maxa)
c        write(*,*) 'new a file ====='
c        write(*,*) (val(I),I=1,maxa)
c        write(*,*) 'new f file      ====='
c        write(*,*) (b(I),I=1,NEQ)
c        do i = 1,neq
c        nl = nbindex(i)
c        nu = nbindex(i+1)-1
c        write(*,*) i,(nbindex(j),j=nl,nu)
c        write(*,*) val(i),(val(j),j=nl,nu)
c        end do
c        end if
c        return
c        do i = 1,maxa
c        if(nbindex(i).lt.0) then
c        print *,'Negative row number, fatal error!!!',iblk
c        return
c        end if
c        end do
c        return
c
c
c       read subnodes file to get the local function number index
c
c       read(23) nodesiblk,(idnode(i),i=1,nodesiblk)
      call initrwf(23,iui,iur)
      iui = iui+1
      nodesiblk=ipool(iui)
      do i=1,1
      iui = iui+1
      idnode(i)=ipool(iui)
      enddo
      call endrwf(23,iui,iur)
c       read(23) nod_iblk,(inode(i),i = 1,nod_iblk)
      call initrwf(23,iui,iur)
      iui = iui+1
      nod_iblk=ipool(iui)
      do i=1,1
      iui = iui+1
      inode(i)=ipool(iui)
      enddo
      call endrwf(23,iui,iur)
c       read(23) knode,(jnode(i),i=1,knode)
      call initrwf(23,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      do i=1,1
      iui = iui+1
      jnode(i)=ipool(iui)
      enddo
      call endrwf(23,iui,iur)
c       read(23) knode,(idnode1(i),i=1,knode)
      call initrwf(23,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      do i=1,1
      iui = iui+1
      idnode1(i)=ipool(iui)
      enddo
      call endrwf(23,iui,iur)
c       read(23) knode,kdgof,((nodvar(j,i),j=1,kdgof),i=1,knode)
      call initrwf(23,iui,iur)
      iui = iui+1
      knode=ipool(iui)
      iui = iui+1
      kdgof=ipool(iui)
      do i=1,nod_iblk
      do j=1,kdgof
      iui = iui+1
      nodvar(j,i)=ipool(iui)
      enddo
      enddo
      call endrwf(23,iui,iur)
c      if(iblk.eq.8)then
c
c       print out the subnodes file for error check
c
c        write(*,*) 'nodesiblk,nod_iblk,knode,kdgof====='
c        write(*,*) nodesiblk,nod_iblk,knode,kdgof
c        write(*,*) 'idnode file ====='
c        write(*,*) (idnode(i),i=1,1)
c        write(*,*) 'inode file ====='
c        write(*,*) (inode(i),i=1,1)
c        write(*,*) 'jnode file ====='
c        write(*,*) (jnode(i),i=1,1)
c        write(*,*) 'idnode1 file ====='
c        write(*,*) (idnode1(i),i=1,1)
c        write(*,*) 'nodvar file ====='
c        write(*,*) ((nodvar(j,i),j=1,kdgof),i=1,200)
c        write(*,*) ((nodvar(j,i),j=1,kdgof),i=1,nod_iblk)
c         open(56,file='aaa',status='unknown',form='formatted')
c        write(56,*) knode,kdgof,((nodvar(j,i),j=1,kdgof)
c     &                                ,i=1,nod_iblk)
c         close(56)
c        print *,'=========knode0,nod_iblk===',knode0,nod_iblk
c       end if
c        return
c
c
c     I modified here for new NEU project
c
        N_update = 0
        do i = 1,nod_iblk
        do j = 1,kdgof
        if(nodvar(j,i).gt.0) then
        N_update = N_update + 1
        nupdate(N_update) = nodvar(j,i)-1
        end if
        end do
        end do
c        print *,N_update,nod_iblk*kdgof
c
c       delete useless na and val and regenerate new gross matrix
c
c
c
        neq2 = 0
        do i = 1,nod_iblk
        do j = 1,kdgof
        if((nodvar(j,i).gt.0).or.(nodvar(j,i).lt.-100)) then
        neq2 = neq2 + 1
        jdiag(neq2) = neq2
        if(nodvar(j,i).lt.0) then
        jdiag(neq2) = -100
        end if
        end if
        end do
        end do
c        print *,"neq,neq2======",neq,neq2,iblk
        if(neq.ne.neq2) print *,'error, different number functions',iblk
c
c       setup an index map between local and global function number
 
        do i = 1,neq1
        mapgl(i)=0
        end do
        neq2 = 0
        do i = 1,nod_iblk
        do j = 1,kdgof
        if(nodvar(j,i).gt.0) then
        neq2 = neq2 + 1
        mapgl(neq2) = nodvar(j,i)
        end if
        if(nodvar(j,i).lt.-100) then
        neq2 = neq2 + 1
        mapgl(neq2) = -nodvar(j,i)-100
        end if
        end do
        end do
c 
c        print *,"N_update,neq2======",N_update,neq2,iblk
c        if(N_update.ne.neq2) print *,'error function number '
c        return
c
c       Now, delete extra boundary function matrix and information
c
c        do i = 1,maxa
c        if(nbindex(i).lt.0) then
c        print *,'Negative row number, fatal error!!!',iblk,1
c        return
c        end if
c        end do
c        return

        do i = 1,neq1
        jdiagaz(i) = nbindex(i)
        end do
c
        do i = 1,neq1
        diag(i) = val(i)
        end do
c
        nnew = 0
        k = 0
        nbindex(1) = N_update + 1
        do i = 1,neq
        if(jdiag(i).gt.0) then
        nl = jdiagaz(i)
        nu = jdiagaz(i+1)
        nadd = 0
        do j = nl,nu-1
c        if(jdiag(nbindex(j)).gt.0) then
        nnew = nnew + 1
        nadd = nadd + 1
        nbindex(nnew+N_update+1) = mapgl(nbindex(j))-1
        val(nnew+N_update+1) = val(j)
c        end if
        end do
        k = k + 1
        nbindex(k+1) = nbindex(k)+nadd
        val(k) = diag(i)
        end if
        end do
        val(k+1) = 0.0d0
c
        k = 0
        do i = 1,neq
        if(jdiag(i).gt.0) then
        k = k + 1
        b(k) = b(i)
        end if
        end do
        if(k.ne.N_update) print *,'error when delete index and val'
        maxa = nnew+N_update+1
c
c       call aztec solver and got the solution.
c
        neqmax0 = neqmax-1

       maxa0 = maxa - 1
        knode0 = knode - 1
c
c         data check!!!
c
        do i = 1,maxa
        if(nbindex(i).lt.0) then
        print *,'Negative row number, fatal error!!!',iblk
        return
        end if
        end do
c        return
        
        

        call callpsolv(neqmax0,knode0,maxa0,N_update,nupdate,
     &                     nbindex,val,b,x,ndata_org,nexternal,
     &                     nupdate_index,nextern_index,neqall,x1)
 
c
c       write out the solution and return
c
        do i = 1,neq
        b(i) = 0.0d0
        end do
        k = 0
        do i = 1,neq
        if(jdiag(i).gt.0) then
        k = k + 1
        b(i) = x(k)
        end if
        end do
c
c        write(25) (b(i),i=1,neq)
c        write(25) (x(i),i=1,N_update)
c
        do i = 1,neqmax
        x(i) = 0.0d0
        end do
        nsource = 0
        if(iblk.ge.1) then
        call sendint(nsource,iblk,neq)
        call sendai(nsource,iblk,mapgl,neq)
        call sendar(nsource,iblk,b,neq)
        end if
        if(iblk.eq.0) then
        neqm = neq
        do k = 1,neq
        x(mapgl(k)) = x(mapgl(k)) + b(k)
        end do
        do i = 1,numblk-1
        call recvint(nsource,i,neq)
        call recvai(nsource,i,mapgl,neq)
        call recvar(nsource,i,b,neq)
        do k = 1,neq
        x(mapgl(k)) = x(mapgl(k)) + b(k)
        end do
        end do
c
        end if
c
c
c        return
c
c
c       server gather all the solutions from all slave nodes ,and send
c       them back including boundary data.
c
        nsource = 0
        if(iblk.ge.1) then
        do k = 1,neq
        b(k) = 0.0d0
        end do
        call sendint(nsource,iblk,neq)
        call sendai(nsource,iblk,mapgl,neq)
        call recvar(iblk,nsource,b,neq)
        end if
        if(iblk.eq.0) then
        do i = 1,numblk-1
        do k = 1,neq
        b(k) = 0.0d0
        end do
        call recvint(nsource,i,neq)
        call recvai(nsource,i,mapgl,neq)
        do k = 1,neq
        b(k) = x(mapgl(k)) 
        end do
        call sendar(i,nsource,b,neq)
        end do
        end if
c
c       restore master node values include boundary 
c
        if(iblk.eq.0) then
        do i = 1,neq1
        mapgl(i)=0
        end do
        neq2 = 0
        do i = 1,nod_iblk
        do j = 1,kdgof
        if(nodvar(j,i).gt.0) then
        neq2 = neq2 + 1
        mapgl(neq2) = nodvar(j,i)
        end if
        if(nodvar(j,i).lt.-100) then
        neq2 = neq2 + 1
        mapgl(neq2) = -nodvar(j,i)-100
        end if
        end do
        end do
        neq = neq2
        do k = 1,neq
        b(k) = x(mapgl(k)) 
        end do
        end if
c
c       write out lsol file
c
      ISTATUS = 0
      call openf(25,27,ISTATUS)
c       write(25) neq
      call initrwf(25,iui,iur)
      iui = iui+1
      ipool(iui)=neq
      call endrwf(25,iui,iur)
c       write(25) (x(k),k = 1,neq)
      call initrwf(25,iui,iur)
      do k=1,neq
      iur = iur+1
      rpool(iur)=b(k)
      enddo
      call endrwf(25,iui,iur)
c
1000    format(i5,10f25.16)
        return
        end
 
 
 
       subroutine callpsolv(neqmax,knode,maxa,N_update,nupdate,
     &                     nbindx,val,b,x,ndata_org,nexternal,
     &                     nupdate_index,nextern_index,neqall,x1)
      implicit real*8 (a-h,o-z)
c
       include "az_aztecf.h"
       include "mpif.h"
c
       integer nproc_config(0:AZ_PROC_SIZE), noptions(0:AZ_OPTIONS_SIZE)
       double precision params(0:AZ_PARAMS_SIZE)
       double precision status(0:AZ_STATUS_SIZE)
       integer N_update,ierror
       integer i
       integer n, neqmax
c
c             See Aztec User's Guide for the variables that follow:
c
       integer nbindx(0:maxa)
       double  precision val(0:maxa)
       double precision b(0:neqmax),x(0:neqmax),x1(0:neqmax)
       integer ndata_org(0:neqmax)
       integer nupdate(0:neqmax), nexternal(0:neqmax)
       integer nupdate_index(0:neqmax), nextern_index(0:neqmax)
c
c      get number of processors and the name of this processor
c
c       call timer(5,1)
       do i = 0,AZ_PROC_SIZE
       nproc_config(i)=0
       end do
       do i = 0,neqmax
       x(i)=0.0d0
       end do
       do i = 0,neqmax
       x1(i)=0.0d0
       end do
       do i = 0,neqmax
       ndata_org(i)=0
       end do
       do i = 0,neqmax
       nexternal(i)=0
       end do
       do i = 0,neqmax
       nupdate_index(i)=0
       end do
       do i = 0,neqmax
       nextern_index(i)=0
       end do
c
       call AZ_set_proc_config(nproc_config, MPI_COMM_WORLD)
c
c       nrow = N_update
c       nrow = neqall*neqall
       nrow = neqall
c       call AZ_read_update(N_update,nupdate,nproc_config,nrow,1,0)
c
c      convert matrix to a local distributed matrix */
c
c       print *,'nproc_config =====',nproc_config
c       print *, 'N_update=======',N_update
c       print *, 'nupdate=========',(nupdate(i),i=0,N_update-1)
c       print *, 'val============',(val(i),i=0,maxa)
c       print *, 'nbindx==========',(nbindx(i),i=0,maxa)
c       print *, 'new f file ====='
c       print *, (b(I),I=0,N_update-1)
c       print *, 'nexternal==========',(nexternal(i),i=0,neqmax)
c       print *, 'nupdate_index==========',(nupdate_index(i),i=0,neqmax)
c       print *, 'nextern_index==========',(nextern_index(i),i=0,neqmax)
c       print *, 'ndata_orgndex==========',(ndata_org(i),i=0,neqmax)
       call AZ_transform(nproc_config,nexternal,nbindx,val,nupdate,
     $                   nupdate_index,nextern_index,ndata_org,
     $                   N_update,0,0,0,0,AZ_MSR_MATRIX)
c
c      Set rhs (delta function at grid center) and initialize guess
c
      call AZ_reorder_vec(b,ndata_org,nupdate_index,0)
       do 350 i = 0, N_update-1
          x(nupdate_index(i)) = 0.0d0
350    continue
 
c      initialize AZTEC options
c
       call AZ_defaults(noptions, params)
c       params(AZ_tol) = 1.0e-8
       params(AZ_tol) = 1.0e-06
c       noptions(AZ_solver) = AZ_cg
c       noptions(AZ_solver) = AZ_gmres
       noptions(AZ_solver) = AZ_bicgstab
c       noptions(AZ_solver) = AZ_cg
       noptions(AZ_precond) = AZ_sym_GS
       noptions(AZ_poly_ord) =5 
       noptions(AZ_scaling) = AZ_sym_diag
c      noptions(AZ_omega) = 2
c      noptions(AZ_omega) = 15
c       noptions(AZ_precond) = AZ_dom_decomp
c       noptions(AZ_precond) = 0
c       noptions(AZ_overlap) = AZ_diag
c       noptions(AZ_overlap) = 10
       noptions(AZ_subdomain_solve) = AZ_ilut
c       noptions(AZ_subdomain_solve) = AZ_ilu
       noptions(AZ_max_iter) = 15000
       params(AZ_ilut_fill) = 3
c      params(AZ_ilut_fill) = 15
       params(AZ_drop) = 1.0e-3
c       params(AZ_drop) = 0.0d0
c
 
c
c      solve the system of equations using b  as the right hand side
c
       call AZ_solve(x,b, noptions, params, 0,nbindx,0,0,
     $               0,val, ndata_org, status, nproc_config)
c
       call AZ_invorder_vec(x,ndata_org,nupdate_index,0,x1)
c        print *,(x1(i),i=0,N_update-1)
        do i = 0,N_update-1
        x(i) = x1(i)
        end do
c       call timer(5,2)
       return
       end
 
