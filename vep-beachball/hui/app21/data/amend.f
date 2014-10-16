c       to make 3D mesh, repeat the 2d mesh in the 3rd dimension
        parameter (n=500000,nf=90000)
        implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      dimension nod(5,n),npoint(2,1000)

      knode=1221
      nqudr=1152
      open(1,file='Salton.flavia.msh',status='old')
      x0=0.d0
      y0=0.d0
      read(1,*)
      read(1,*)
      do i=1,knode
      read(1,*) ii,x(i),y(i)
      enddo
      read(1,*)
      read(1,*)      
      do i=1,nqudr
      read(1,*) ii,nod(1,i),nod(2,i),nod(3,i),nod(4,i),nod(5,i)
      enddo
      write(*,*) 'nqudr =',nqudr
      close(1)
      open(1,file='duplicates',status='old')
      read(1,*) nump
      do i=1,nump
      read(1,*) npoint(1,i),npoint(2,i)
      enddo
      close(1)
      
      do i=1,nqudr
      do j=1,4
      do k=1,nump
c      write(*,*) nod(j,i),npoint(2,k)
      if (nod(j,i).eq.npoint(2,k)) then 
      write(*,*) i,j,nod(j,i),npoint(1,k)
      nod(j,i)=npoint(1,k)
      endif
      enddo
      enddo
      enddo

1111  continue
      do i=1,nqudr
      x1=x(nod(1,i))
      y1=y(nod(1,i))
      x2=x(nod(2,i))
      y2=y(nod(2,i))
      x3=x(nod(3,i))
      y3=y(nod(3,i))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      if (abs(s1).lt.1.0e-6) then
      write(*,*)i,nod(1,i),nod(2,i),nod(3,i)
      do j=i,nqudr-1
      do k=1,5
      nod(k,j)=nod(k,j+1)
      enddo
      enddo
      nqudr=nqudr-1
      goto 1111
      endif      
      enddo
      write(*,*) 'nqudr =',nqudr

      open(30,file='amend.flavia.msh',status='unknown',
     +        form='formatted')
c      write(30,*)
      write(30,*)'Mesh "q4" Dimension 2 Elemtype Quadrilateral Nnode 4'
c      write(30,*)'Mesh "q4" Dimension 2 Elemtype Triangle Nnode 3'
      write(30,*)'Coordinates'
      do i = 1,knode
      write(30,1001) i,x(i),y(i)
      end do
      write(30,*)'End coordinates'
      write(30,*) 'Elements'
      do i=1,nqudr
      write(30,1101)i,(nod(j,i),j=1,5)
      end do
      write(30,*)'End elements'
      close(30)

      open(30,file='data1',status='unknown',
     +        form='formatted')
      write(30,*)knode
      do i = 1,knode
      write(30,1001) i,x(i),y(i)
      end do
      write(30,*) nqudr
      do i=1,nqudr
      write(30,1101)i,(nod(j,i),j=1,5)
      end do
      close(30)

1001  format(i10,3es15.6)
1101  format(15i10)

      end 
