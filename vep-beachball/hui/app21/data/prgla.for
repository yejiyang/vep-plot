      implicit real*8 (a-h,l,o-z)
      parameter (n=10000)
      dimension point(2,n),cx(n),cy(n),itri(4,n),mat(n)
 
      open(10,file='mesh.cor',status='old')
      read(10,*) num
      do i=1,num
      read(10,*) point(1,i),point(2,i)
      enddo
      close (10)
      open(10,file='mesh.elm',status='old')
      read(10,*) ntri
      do i=1,ntri
      read(10,*) itri(1,i),itri(2,i),itri(3,i),itri(4,i),mat(i)
      enddo


c....just for the case
      open(30,file='data',status='unknown',
     +        form='formatted')
      write(30,*)num
      do i = 1,num
      write(30,1001) i,point(1,i),point(2,i)
      end do
      write(30,*) ntri
      do i=1,ntri
      write(30,1101)i,(itri(j,i),j=1,4),mat(i)
      end do
      close(30)
      goto 1111



c
c      calculate the barycenter of the fault system
c
      xcent=0.d0
      ycent=0.d0
      do i=1,num
      xcent=xcent+point(1,i)
      ycent=ycent+point(2,i)
      enddo
      xcent=xcent/(num*1.0)
      ycent=ycent/(num*1.0)
      
      write(*,1200) num,xcent,ycent
      
c
c      projection translation to the plane
c
      pi=4*datan(1.d0)/180.d0
      R=6378.e3
      a=xcent*pi
      b=ycent*pi
      write(*,*) 'a=',a/pi,'  b=',b/pi
      do i=1,num
      xx=point(1,i)*pi
      yy=point(2,i)*pi
      call trans(R,xx,yy,a,b,c,d)
      cx(i)=c
      cy(i)=d
      enddo
c
c      NW40
c
      do i=1,num
      alpha=0.d0*pi
      a11=cos(alpha)
      a12=sin(alpha)
      a21=-sin(alpha)
      a22=cos(alpha)
      point(1,i)=a11*cx(i)+a12*cy(i)
      point(2,i)=a21*cx(i)+a22*cy(i)      
      enddo

      xmax=point(1,1)
      xmin=point(1,1)
      ymax=point(2,1)
      ymin=point(2,1)
      do i=1,num
      if (xmax.lt.point(1,i))      xmax=point(1,i)
      if (xmin.gt.point(1,i))      xmin=point(1,i)
      if (ymax.lt.point(2,i))      ymax=point(2,i)
      if (ymin.gt.point(2,i))      ymin=point(2,i)
      enddo
      write(*,1100)xmax,xmin,ymax,ymin

      open(30,file='data',status='unknown',
     +        form='formatted')
      write(30,*)num
      do i = 1,num
      write(30,1001) i,point(1,i),point(2,i)
      end do
      write(30,*) ntri
      do i=1,ntri
      write(30,1101)i,(itri(j,i),j=1,4),mat(i)
      end do
      close(30)

1111  continue

1100  format(10f15.5)      
1200  format(i9,10f15.5)      
1001  format(i10,3es15.5)
1101  format(15i10)
      end

        subroutine trans(R,x,y,a,b,c,d)
        implicit real*8 (a-h,l,o-z)
        dimension p(3,3),q(3,3),s(3,3),xx(3),xxx(3),xy(3)
        p(1,1)=cos(a)
        p(1,2)=sin(a)
        p(1,3)=0.
        p(2,1)=-sin(a)
        p(2,2)=cos(a)
        p(2,3)=0.
        p(3,1)=0.
        p(3,2)=0.
        p(3,3)=1.
        q(1,1)=cos(b)
        q(1,2)=0
        q(1,3)=sin(b)
        q(2,1)=0.
        q(2,2)=1.
        q(2,3)=0.
        q(3,1)=-sin(b)
        q(3,2)=0.
        q(3,3)=cos(b)
        do i=1,3
        do j=1,3
        s(i,j)=0.
        enddo
        enddo
        xx(1)=R*cos(x)*cos(y)
        xx(2)=R*sin(x)*cos(y)
        xx(3)=R*sin(y)
        do i=1,3
        do j=1,3
        xxx(i)=0.
        enddo
        enddo
        do i=1,3
        do j=1,3
        xxx(i)=xxx(i)+p(i,j)*xx(j)
        enddo
        enddo
        do i=1,3
        xy(i)=0.
        do j=1,3
        xy(i)=xy(i)+q(i,j)*xxx(j)
        enddo
        enddo
        c=xy(2)
        d=xy(3)
        return
        end
