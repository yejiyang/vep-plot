      implicit real*8 (a-h,o-z)
      character*22 fname,filename(20)
      common /aa/ ia(20000000)
      open(2,file='coor0',form='unformatted',status='old')
      read(2) np,nd
      close(2)
      WRITE(*,*) 'np,nd =',np,nd
      open(3,file='elem0',form='unformatted',status='old')
      read(3) num,nnode
      close(3)
      WRITE(*,*) 'num,nnode =',num,nnode
      kna1=np*20*2
      knb1=np*nd*2
      knb2=num*nnode*1
      if (knb2/2*2 .lt. knb2) knb2=knb2+1
      kna0=1
      kna1=kna1+kna0
      knb0=1
      knb1=knb1+knb0
      knb2=knb2+knb1
      call gidpostf(np,nd,num,nnode,ia(kna0),ia(knb0),
     *ia(knb1),
     *filename)
      end
      subroutine gidpostf(np,nd,num,nnode,u,coor,node,
     *filename)
      implicit real*8 (a-h,o-z)
      character*22 fname,filename(20)
      logical filflg
      character row*60,res*6,c*1,mate
      dimension ndof(10),ns(10),u(np,20),coor(np,nd),node(num*nnode)
1000  format(1x,i10,6e11.4)
1100  format(1x,i10,28i10)
 
      mate='y'
      inquire(file='material',exist=filflg)
      if (filflg) then
      open(1,file='material',form='formatted',status='old')
      read(1,'(a1)') mate
      close(1)
      endif
 
      open(1,file='tempgid',form='formatted',status='old')
      numarg = 0
      numarg = numarg+1
      call getarg(numarg,fname)
      write(*,*) 'fname = ',fname
      open(3,file=fname,form='formatted',status='unknown')
      read(1,*) nf,(ndof(i),i=1,nf),(ns(i),i=1,nf)
      res(1:6) = 'munod '
      write(3,*) 'GID Post Results File 1.0'
      write(3,*)
 
      do n=1,nf
      WRITE(UNIT=C,FMT='(I1)') n-1
      if (n.gt.1) res(6:6) = c
      write(*,*) 'res ==== ',res
      open(2,file=res,form='unformatted',status='old')
      read(2) ((u(i,j),i=1,np),j=1,ndof(n))
      close(2)
      do k=1,ns(n)
       do m=1,2
        read(1,'(a60)') row
        do i=60,1,-1
        if (row(i:i).ne.' ') goto 1
        enddo
1       continue
        write(3,*) (row(j:j),j=1,i)
c        write(3,*) row
       enddo
       write(3,*) 'Values'
       do i=1,np
        if (ns(n).eq.1) then
        write(3,1000) i,(u(i,j),j=1,ndof(n))
        else
        write(3,1000) i,u(i,k)
        endif
       enddo
       write(3,*) 'end Values'
       write(3,*)
      enddo
      enddo
      close(3)
 
      numarg = numarg+1
      call getarg(numarg,fname)
      write(*,*) 'fname = ',fname
      open(3,file=fname,form='formatted',status='unknown')
      inquire(file='coor1',exist=filflg)
      if (filflg) then
      open(2,file='coor1',form='unformatted',status='old')
      else
      open(2,file='coor0',form='unformatted',status='old')
      endif
      read(2) np,nd,((coor(i,j),j=1,nd),i=1,np)
      close(2)
 
c      open(2,file=' ',form='formatted',status='old')
c      read(2,*) c
c      read(2,*) nelem
c      close(2)
 
      open(2,file='elem0',form='unformatted',status='old')
      n = 0
3     continue
      n = n+1
       read(1,'(a60)') row
        do i=60,1,-1
         if (row(i:i).ne.' ') goto 2
        enddo
2       continue
      if (row(1:1).ne.'#') then
       write(3,*) (row(j:j),j=1,i)
c       write(3,*) row
       write(3,*) 'Coordinates'
       if (n.eq.1) then
       do i=1,np
        write(3,1000) i,(coor(i,j),j=1,nd)
       enddo
       endif
       write(3,*) 'End coordinates'
       read(2) ne,nnode,(node(i),i=1,ne*nnode)
       write(3,*) 'Elements'
       do i=1,ne
        write(3,1100) i,(node((i-1)*nnode+j),j=1,nnode-1)
       enddo
       write(3,*) 'End elements'
       if (mate.eq.'y' .or. mate.eq.'Y') read(2) ne
       goto 3
      endif
      close(2)
      close(3)
 
      end
