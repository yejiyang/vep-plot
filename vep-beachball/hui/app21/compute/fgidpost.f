c      this program deal with the results data
c      creat the selected data files, the data is on surface, 
c      upper crust and lower crust ...

      program gidpost
      implicit real*8 (a-h,l,o-z)
       
      dimension nqudr(5,1000000),coor(2,100000)
      dimension u0(3,1000000),u1(6,1000000)
      dimension coor3(3,2000000),nqudr3(9,2000000)
       
      open(10,file='../data/data',status='old')
      read(10,*) np2
      do i=1,np2
      read(10,*) ii,coor(1,i),coor(2,i)
      enddo
      read(10,*) nqu2
      do i=1,nqu2
      read(10,*) ii,(nqudr(j,i),j=1,5)
      enddo      
      close(10)
      open(10,file='coor0',form='unformatted',status='old')
      read(10) knode,knodof,((coor3(j,i),j=1,knodof),i=1,knode)
      close(10)
      open(10,file='elem0',form='unformatted',status='old')
      read(10) nelem,melem,((nqudr3(j,i),j=1,melem),i=1,nelem)
      close(10)

      nlayer=6
      np3=np2*nlayer
      nqu3=nqu2*(nlayer-1)

c     munod
      kdgof=3
      knode=np3
      write(*,*)knode
      open(10,file='munod',status='old',form='unformatted')
      read(10) ((u0(j,i),i=1,knode),j=1,kdgof)
      close(10)
      kdgof=6
      knode=np3
      write(*,*)knode
      open(10,file='munod1',status='old',form='unformatted')
      read(10) ((u1(j,i),i=1,knode),j=1,kdgof)
      close(10)
c
c     write out results file
c
      open(40,file='res.flavia.msh',status='unknown',
     +        form='formatted')
      write(40,*)  '   Mesh "Whole" Dimension 3  Elemtype',
     +' Hexahedra Nnode 8'
      write(40,*)  ' Coordinates'
      do i=1,knode
      write(40,1000) i,(coor3(j,i),j=1,3)
      enddo
      write(40,*) ' End coordinates'
      write(40,*) ' Elements'
      do i=1,nelem
      write(40,1200) i,(nqudr3(j,i),j=1,melem)
      enddo
      write(40,*) ' End elements'
      close(40)

      open(21,file='res.flavia.res',form='formatted',
     *  status='unknown')
      write(21,*)'GID Post Results File 1.0'
      write(21,*)
      write(21,*) 'Result "disp" "Load Analysis" ',1,
     *  ' Vector OnNodes'
      write(21,*) 'ComponentNames "u" "v" "w"'
      write(21,*) 'Values'
      do i=1,knode
      write(21,1100) i,u0(1,i),u0(2,i),u0(3,i)
      end do
      write(21,*) 'end Values'
      write(21,*)

c     stress matrix
      write(21,*) 'Result "stress" "Load Analysis" ',1,
     *  ' Matrix OnNodes'
      write(21,*) 'ComponentNames "dxx" "dyy"',
     +' "dzz" "dyz" "dxz" "dxy"'
      write(21,*) 'Values'
      do i=1,knode
      write(21,1100) i,(u1(j,i),j=1,6)
      enddo
      write(21,*) 'end Values'
      close(21)
      
      open(11,file='res.txt',status='unknown',form='formatted')
      do i=1,np2
      k=np2*(nlayer-1)+i
      write(11,1000)i,(coor(j,i),j=1,2),
     *  (u0(j,k),j=1,3),(u1(j,k),j=1,6)
      enddo
      close(11)

c      sample the node of element to grid 1 degree multi 1 degree 
c      velocity and energy
      write(*,*)
      write(*,*)'... Sample the grid ...'
      call process2(nlayer,np2,nqu2,coor,nqudr,u0,u1)


1000  format(i9,2f15.5,20es15.5)
1100  format(i9,20es15.5)
1200  format(10i9)
      
      end

       
c******************************************************************      
c
c      sample the node of element to grid 0.1 degree multi 0.1 degree 
c      velocity and energy
c
      SUBROUTINE process2(nlayer,np2,nqu2,coor,nqudr,u0,u1)
      implicit real*8 (a-h,l,o-z)

      dimension nqudr(5,1000000),coor(2,100000)
      dimension u0(3,1000000),u1(6,1000000)
      dimension coor3(3,2000000)
      dimension xxd(4),yyd(4),zzd(4)
      
      pi=4*datan(1.d0)/180.d0
      R=6378.e3

      xcent=-119.0911376567164
      ycent=35.628061910447770
      x0=250000.d0
      y0=300000.d0
      xmin=-98800.43893
      ymin=-491461.11950

      step=0.1d0
      lon=-125.0d0-step
      lat=32.0d0-step
      a=xcent*pi
      b=ycent*pi
      write(*,*) 'a=',a/pi,'  b=',b/pi
      
      alpha=37.5d0*pi
      a11=cos(alpha)
      a12=sin(alpha)
      a21=-sin(alpha)
      a22=cos(alpha)

      open(12,file='grid.txt',status='unknown',form='formatted')

      do while ( lon .le. -115.0d0 ) 
      lon=lon+step
      lat=32.0d0-step
      do while ( lat .le. 40.0d0 ) 
      flag=-1.0
      lat=lat+step
c     projection
      xx=lon*pi
      yy=lat*pi
      call trans(R,xx,yy,a,b,c,d)
      cx=c
      cy=d
      xxx0=a11*cx+a12*cy-xmin+x0
      yyy0=a21*cx+a22*cy-ymin+y0
      x3=xxx0
      y3=yyy0      
      

      do k=1,nqu2
      x1=coor(1,nqudr(1,k))
      y1=coor(2,nqudr(1,k))
      x2=coor(1,nqudr(2,k))
      y2=coor(2,nqudr(2,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
      if (s1 .lt. 0.d0) goto 100
      
      x1=coor(1,nqudr(2,k))
      y1=coor(2,nqudr(2,k))
      x2=coor(1,nqudr(3,k))
      y2=coor(2,nqudr(3,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
      if (s1 .lt. 0.d0) goto 100

      x1=coor(1,nqudr(3,k))
      y1=coor(2,nqudr(3,k))
      x2=coor(1,nqudr(4,k))
      y2=coor(2,nqudr(4,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
      if (s1 .lt. 0.d0) goto 100

      x1=coor(1,nqudr(4,k))
      y1=coor(2,nqudr(4,k))
      x2=coor(1,nqudr(1,k))
      y2=coor(2,nqudr(1,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
      if (s1 .lt. 0.d0) goto 100
      
      flag=1.0
      
c      write(*,*) k

       xxd(1) = coor(1,nqudr(1,k))
       xxd(2) = coor(1,nqudr(2,k))
       xxd(3) = coor(1,nqudr(3,k))
       xxd(4) = coor(1,nqudr(4,k))
      
       yyd(1) = coor(2,nqudr(1,k))
       yyd(2) = coor(2,nqudr(2,k))
       yyd(3) = coor(2,nqudr(3,k))
       yyd(4) = coor(2,nqudr(4,k))
      
       zzd(1) = u0(1,nqudr(1,k))
       zzd(2) = u0(1,nqudr(2,k))
       zzd(3) = u0(1,nqudr(3,k))
       zzd(4) = u0(1,nqudr(4,k))
      
       xxd0 = xxx0
       yyd0 = yyy0
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)       
       u01= top*1000.d0

       zzd(1) = u0(2,nqudr(1,k))
       zzd(2) = u0(2,nqudr(2,k))
       zzd(3) = u0(2,nqudr(3,k))
       zzd(4) = u0(2,nqudr(4,k))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u02= top*1000.d0
      
c      u01=(u0(1,nqudr(1,k))+u0(1,nqudr(2,k))
c     +  +u0(1,nqudr(3,k))+u0(1,nqudr(4,k)))*1000.d0/4.0
c      u02=(u0(2,nqudr(1,k))+u0(2,nqudr(2,k))
c     +  +u0(2,nqudr(3,k))+u0(2,nqudr(4,k)))*1000.d0/4.0
      length=sqrt(u01*u01+u02*u02)
      if (abs(u01) .lt. 1e-1) then
       ang=90.0d0
      else ang=atan(u02/u01)/pi
      endif
      uu01=length*cos(ang*pi+alpha)
      uu02=length*sin(ang*pi+alpha)

      write(12,1000) 1,lon,lat,uu01,uu02
      exit
100   continue      
      enddo   
      
      enddo
      enddo
      close(12)

1000  format(i9,10f13.5,20es15.5)
      return
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






c******************************************************************      

      SUBROUTINE CJCBJ(A,N,EPS,V)
      implicit real*8 (a-h,l,o-z)
      DIMENSION A(N,N),V(N,N)
      DOUBLE PRECISION A,V,FF,FM,CN,SN,OMEGA,X,Y
c      INTEGER P,Q
      DO 20 I=1,N
        V(I,I)=1.0
        DO 10 J=1,N
          IF (I.NE.J) V(I,J)=0.0
10      CONTINUE
20    CONTINUE
          FF=0.0
      DO 500 I=2,N
      DO 500 J=1,I-1
500   FF=FF+A(I,J)*A(I,J)
      FF=SQRT(2.0*FF)
205   FF=FF/(1.0*N)
25    DO 30 I=2,N
      DO 30 J=1,I-1
        IF (ABS(A(I,J)).GE.FF) THEN
          k=I
          m=J
          GOTO 600
        END IF
30    CONTINUE
      IF (FF.GE.EPS) GOTO 205
      RETURN
600   X=-A(k,m)
      Y=(A(m,m)-A(k,k))/2.0
      OMEGA=X/SQRT(X*X+Y*Y)
      IF (Y.LT.0.0) OMEGA=-OMEGA
      SN=1.0+SQRT(1.0-OMEGA*OMEGA)
      SN=OMEGA/SQRT(2.0*SN)
      CN=SQRT(1.0-SN*SN)
      FM=A(k,k)
      A(k,k)=FM*CN*CN+A(m,m)*SN*SN+A(k,m)*OMEGA
      A(m,m)=FM*SN*SN+A(m,m)*CN*CN-A(k,m)*OMEGA
      A(k,m)=0.0
      A(m,k)=0.0
      DO 60 J=1,N
        IF ((J.NE.k).AND.(J.NE.m)) THEN
          FM=A(k,J)
          A(k,J)=FM*CN+A(m,J)*SN
          A(m,J)=-FM*SN+A(m,J)*CN
        END IF
60    CONTINUE
      
      DO 70 I=1,N
        IF ((I.NE.k).AND.(I.NE.m)) THEN
          FM=A(I,k)
          A(I,k)=FM*CN+A(I,m)*SN
          A(I,m)=-FM*SN+A(I,m)*CN
        END IF
70    CONTINUE
      DO 80 I=1,N
        FM=V(I,k)
        V(I,k)=FM*CN+V(I,m)*SN
        V(I,m)=-FM*SN+V(I,m)*CN
80    CONTINUE
      GOTO 25
      END
       


        subroutine locate_qudr(xx,yy,x0,y0,r,s)
        implicit real*8 (a-h,o-z)
        dimension rnod(2,4),p(2),rc(2,2),cr(2,2),r0(2)
        dimension xx(4),yy(4)                    
        do  i=1,4
         rnod(1,i)=xx(i)
         rnod(2,i)=yy(i)
        end do
        r0(1)=x0
        r0(2)=y0
        call locate_d2q4(4,2,rnod,r0,p,rc,cr)
        r=p(1)
        s=p(2)
c----------------------------
        return
        end

        subroutine locate_d2q4(nnode,ndm,rnod,r0,p,rc,cr)
c	dimension rnod(2,4),r(2),p(2),rc(2,2),cr(2,2),r0(2),dr(2)
        implicit real*8 (a-h,o-z)
	dimension rnod(ndm,nnode),r0(ndm),p(ndm),rc(ndm,ndm),
     &            cr(ndm,ndm),r(3),dr(3)
c	real*8 err,errmax,d
c .......................................................... c
c .... F(p) = r(p) - r0
c .... F(p+dp) = F(p) + F'(p)*dp = r(p)-r0 + {r/p}*dp = 0
c .... dp = - {p/r}*(r(p)-r0)
c .... p+dp = p - {p/r}*(r(p)-r0)
c .... p+dp = p - cr*( r(p) - r0 )
c .......................................................... c
	d=0.0d0
	do 40 i=1,ndm
	do 20 j=1,nnode
	d = d + (rnod(i,j)-r0(i))**2
20	continue
40	continue
c        write(*,*) 'd = ',d
	errmax = d*1.d-6/nnode
	itmax = 10
	it = 0
	do 100 i=1,ndm
	p(i)=0.0
100	continue
1	continue
	call locate_elemt_qudr(nnode,ndm,ndm,p,r,rnod,rc,cr,det)
	do 200 i=1,ndm
	dr(i) = r(i) - r0(i)
200	continue
	err = 0.0d0
	do 300 i=1,ndm
	err = err+dr(i)**2
300	continue
c        write(*,*) 'err,dr =',err,(dr(i),i=1,ndm)
	if (err.lt.errmax .or. it.ge.itmax) goto 2
	do 400 i=1,ndm
	do 401 j=1,ndm
	p(i) = p(i) - cr(i,j)*dr(j)
401	continue
400	continue
	it = it+1
	goto 1
2	continue
c        write(*,*) 'p = ',p
	return
	end
      
      subroutine locate_elemt_qudr(nnode,nrefc,ncoor,refc,
     * coor,coorr,rc,cr,det)
c      implicit real*4 (a-h,o-z)
      implicit real*8 (a-h,o-z)
c     implicit integer*2 (i-n)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     *            coorr(ncoor,nnode),coor(ncoor)
      call locate_telem_qudr(refc,coor,coorr,rc)
      n=nrefc
      m=n*2
      det = 1.0
      do 10 i=1,n
      do 11 j=1,n
      if (i.le.ncoor) a(i,j) = rc(i,j)
      if (i.gt.ncoor) a(i,j)=1.0
      a(i,n+j)=0.0
      if (i.eq.j) a(i,n+i) = 1.0
11    continue
10    continue
c     write(*,*) 'a ='
c     do 21 i=1,n
c21   write(*,8) (a(i,j),j=1,m)
      do 400 i=1,n
      amax = 0.0
      l = 0
      do 50 j=i,n
      c = abs(a(j,i))
      if (c.le.amax) goto 50
      amax = c
      l = j
50    continue
c      write(*,*) amax
c      pause
      do 60 k=1,m
      c = a(l,k)
      a(l,k) = a(i,k)
      a(i,k) = c
60    continue
      c = a(i,i)
      det = c*det
      do 100 k=i+1,m
      a(i,k) = a(i,k)/c
100   continue
      do 300 j=1,n
      if (i.eq.j) goto 300
      do 200 k=i+1,m
      a(j,k) = a(j,k)-a(i,k)*a(j,i)
200   continue
c     write(*,*) 'i =',i,'  j =',j,'  a ='
c     do 11 ii=1,n
c11   write(*,8) (a(ii,jj),jj=1,m)
300   continue
400   continue
      do 500 i=1,nrefc
      do 501 j=1,ncoor
      cr(i,j) = a(i,n+j)
501   continue
500   continue
c     write(*,*) 'a ='
c     do 22 i=1,n
c22   write(*,8) (a(i,j),j=1,m)
c     write(*,*) 'rc ='
c     do 24 i=1,ncoor
c24   write(*,8) (rc(i,j),j=1,nrefc)
c     write(*,*) 'cr ='
c     do 23 i=1,nrefc
c23   write(*,8) (cr(i,j),j=1,ncoor)
c     write(*,*) 'det =',det
      if (det.lt.0.0) det=-det
c     write(*,*) 'det =',det
8     format(1x,6f12.3)
      end

      subroutine locate_telem_qudr(refc,coor,coorr,rc)
      implicit real*8 (a-h,o-z)
      dimension x(4),y(4),refc(2),coor(2),coorr(2,4),rc(2,2)
      p=refc(1)
      q=refc(2)
      do 1 n=1,4
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
1     continue
      coor(1)=+(+(+1.-p)*(+1.-q)/4.)*x(1)+(+(+1.+p)*(+1.-
     /q)/4.)*x(2)+(+(+1.+p)*(+1.+q)/4.)*x(3)+(+(+1.-p)*(+1.+
     /q)/4.)*x(4)
      coor(2)=+(+(+1.-p)*(+1.-q)/4.)*y(1)+(+(+1.+p)*(+1.-
     /q)/4.)*y(2)+(+(+1.+p)*(+1.+q)/4.)*y(3)+(+(+1.-p)*(+1.+
     /q)/4.)*y(4)
      rc(1,1)=+(+(-1.)*(+1.-q)/4.)*x(1)+(+(+1.)*(+1.-q)/
     /4.)*x(2)+(+(+1.)*(+1.+q)/4.)*x(3)+(+(-1.)*(+1.+q)/4.)*x(4)
      rc(2,1)=+(+(-1.)*(+1.-q)/4.)*y(1)+(+(+1.)*(+1.-q)/
     /4.)*y(2)+(+(+1.)*(+1.+q)/4.)*y(3)+(+(-1.)*(+1.+q)/4.)*y(4)
      rc(1,2)=+(+(+1.-p)*(-1.)/4.)*x(1)+(+(+1.+p)*(-1.)/
     /4.)*x(2)+(+(+1.+p)*(+1.)/4.)*x(3)+(+(+1.-p)*(+1.)/4.)*x(4)
      rc(2,2)=+(+(+1.-p)*(-1.)/4.)*y(1)+(+(+1.+p)*(-1.)/
     /4.)*y(2)+(+(+1.+p)*(+1.)/4.)*y(3)+(+(+1.-p)*(+1.)/4.)*y(4)
      end

      subroutine inter_quadri(r,s,zz,top)
      implicit real*8 (a-h,o-z)
      dimension zz(4),x(24)
      rx=r
      ry=s
      do n=1,4
      x(n)=zz(n)
      end do
           top=+(+(+1.-rx)/2.*(+1.-ry)/2.)*x(1)
     *         +(+(+1.+rx)/2.*(+1.-ry)/2.)*x(2)
     *         +(+(+1.+rx)/2.*(+1.+ry)/2.)*x(3)
     *         +(+(+1.-rx)/2.*(+1.+ry)/2.)*x(4)
      return
      end
