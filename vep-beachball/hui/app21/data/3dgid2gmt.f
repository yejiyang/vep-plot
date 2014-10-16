c     this program convert the gid file to GMT file

      parameter (n=500000,nf=90000)
      implicit real*8 (a-h,o-z)      
      dimension xt(n),yt(n),zt(n)
      dimension line(2,n)
      dimension itri(9,n),mat(n)

c.....2D mesh
      open(121,file='data',form='formatted',status='old')
      read(121,*) nnode
      do i=1,nnode
      read(121,*) ii,xx,yy
      xt(i)=xx/1.d3
      yt(i)=yy/1.d3
      enddo
      read(121,*) ntri
      do i=1,ntri
      read(121,*) ii,itri(1,i),itri(2,i),itri(3,i),itri(4,i),
     *  mat(i)
      enddo
      close(121)
      write(*,*) '  nnode = ', nnode, ';      ntri =',ntri

c.....mesh boundaries
      open(11,file='boundaries',status='old')
      mbb=1000
c         write(*,*)'surrounding sides',mbb
      k=0
      do i=1,mbb
      read(11,*,end=37)ia,ib,ic,idx,idy,idm
        write(*,*)i,ic,idx,idy,idm
        do j=1,ic
        k=k+1
        read(11,*)ia,line(1,k),line(2,k)
        enddo
      enddo
37      close(11)
      write(*,*) '  Number of boundaries = ', k        
        
c.....rewrite mesh

      z0=0.d0
      step=60.d0
      open(20,file='3dgrid.txt',status='unknown',form='formatted')

      do i=1,ntri
      write(20,'(A,i1)') '> -Z',1
      do j=1,4
      write(20,'(3f15.5)') xt(itri(j,i)),yt(itri(j,i)),z0
      enddo
      enddo
      
      zstep=10.d0
      nlayers=6
      do ii=1,nlayers
      do i=1,k
      if ((xt(line(1,i)).le.1).and.(xt(line(2,i)).le.1)
     *  .or.(yt(line(1,i)).le.1).and.(yt(line(2,i)).le.1)) then
      if (ii.le.4) then
      write(20,'(A,i1)') '> -Z',1
      else
      write(20,'(A,i1)') '> -Z',2
      endif
      write(20,'(3f15.5)') xt(line(1,i)),yt(line(1,i)),
     *   z0+(ii-1)*zstep
      write(20,'(3f15.5)') xt(line(2,i)),yt(line(2,i)),
     *   z0+(ii-1)*zstep
      write(20,'(3f15.5)') xt(line(2,i)),yt(line(2,i)),z0+ii*zstep 
      write(20,'(3f15.5)') xt(line(1,i)),yt(line(1,i)),z0+ii*zstep
      endif 
      enddo
      enddo

      do i=1,ntri
      write(20,'(A,i1)') '> -Z',2
      do j=1,4
      write(20,'(3f15.5)') xt(itri(j,i)),yt(itri(j,i)),z0+step
      enddo
      enddo
      
      
      
      close(20)

1000  format(10f14.4)    
1100  format(A,10f14.4)
      end
      
      
