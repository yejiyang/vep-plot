c      this program deal with the results data-geography cordination
c      creat the selected data files, the data is on surface, 
c      upper crust and lower crust ...
c       including both Hui and Yang's method to plot beach ball

      program gidpost
      implicit real*8 (a-h,l,o-z)
       
      dimension nqudr(5,1000000),coor(2,100000)
      dimension u0(3,1000000),u1(6,1000000),u2(6,1000000)
      dimension coor3(3,2000000),nqudr3(9,2000000)
      dimension nleft(100),nright(100)
      dimension label(1000000)
  
        open(11,file='partition.dat',status='old')
        do i=1,17
        read(11,*)
        enddo
        read(11,*)t0,tmax,dt
        close (11)
        dt=dt/3.15e7
        write(*,*)'dt=',dt

c      open(10,file='../data4_check/data1',status='old')
      open(10,file='data1',status='old')
      read(10,*) np2
      do i=1,np2
      read(10,*)
      enddo
      read(10,*) nqu2
      close(10)
      open(10,file='amend.flavia.msh',status='old')
      read(10,*)
      read(10,*)
      do i=1,np2
      read(10,*) ii,coor(1,i),coor(2,i)
      enddo
      read(10,*)
      read(10,*)
      do i=1,nqu2
      read(10,*) ii,(nqudr(j,i),j=1,5)
      enddo
      close(10)
      do i=1,np2
      label(i)=-1
      enddo
      do i=1,nqu2
      do j=1,4
      label(nqudr(j,i))=1
      enddo
      enddo
      
      open(10,file='coor0',form='unformatted',status='old')
      read(10) knode,kdgof,((coor3(j,i),j=1,kdgof),i=1,knode)
      close(10)
      open(10,file='elem0',form='unformatted',status='old')
      read(10) nelem,melem,((nqudr3(j,i),j=1,melem),i=1,nelem)
      close(10)

      nlayer=knode/np2
      np3=knode
      nqu3=nelem

      write(*,*)'      knode       nelem       nlayer'
      write(*,*)knode,nelem,nlayer

      open(10,file='munod',status='old',form='unformatted')
      read(10) ((u0(j,i),i=1,knode),j=1,kdgof)
      close(10)
      kdgof=6
      open(10,file='munod1',status='old',form='unformatted')
      read(10) ((u1(j,i),i=1,knode),j=1,kdgof)
      close(10)

c
c
c
       do i=1,knode
       do j=1,3
       u0(j,i)=u0(j,i)/dt
       enddo
       do j=1,6
       u1(j,i)=u1(j,i)/dt
       enddo
       enddo

        
c        write(*,*)np2,nlayer
c        write(*,*)i-1,nleft(i-1),u0(1,nleft(i-1)),u0(1,nleft(i-1))
c        write(*,*)i-1,nright(i-1),u0(1,nright(i-1)),u0(2,nright(i-1))

c
c      Calculate strain rates on nodes according to predicted 
c      velocity
c
      write(*,*)
c      write(*,*)'... Recalculate strain rates ...'
c      call process1(nlayer,np2,nqu2,coor,nqudr,u0,u2)

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
     +' "dzz" "j2" "p" "e"'
      write(21,*) 'Values'
      do i=1,knode
      write(21,1100) i,(u1(j,i),j=1,6)
      enddo
      write(21,*) 'end Values'

c     strain matrix
      write(21,*) 'Result "strain" "Load Analysis" ',1,
     *  ' Matrix OnNodes'
      write(21,*) 'ComponentNames "dxx" "dyy"',
     +' "dxy" "non" "non" "non"'
      write(21,*) 'Values'
      do i=1,knode
      write(21,1100) i,(u2(j,i),j=1,6)
      enddo
      write(21,*) 'end Values'
      close(21)
      
c      open(11,file='res.txt',status='unknown',form='formatted')
c      do i=1,np2
c      k=np2*(nlayer-1)+i
c      if (label(i).gt.0) write(11,1000)i,(coor(j,i),j=1,2),
c     *  (u0(j,k),j=1,3),(u1(j,k),j=1,6)
c      enddo
c      close(11)


c      sample the node of element to grid 1 degree multi 1 degree 
c      velocity and energy
      write(*,*)
      write(*,*)'... Sample the grid ...'
      call process2(nlayer,np2,nqu2,coor,nqudr,u0,u1)



1000  format(i9,2f15.5,20es15.5)
1100  format(i9,20es15.5)
1200  format(10i9)
1300  format(2i9,10es13.5)
1400  format(i9,2f13.3,f13.1)      
      end

       
c******************************************************************      
c
c      Calculate strain rates on nodes according to predicted 
c      velocity
c
      SUBROUTINE process1(nlayer,np2,nqu2,coor,nqudr,u0,u1)
      implicit real*8 (a-h,l,o-z)

      dimension nqudr(5,1000000),coor(2,100000)
      dimension u0(3,1000000),u1(6,1000000)
      dimension index(1000000),ipoint(8)
      dimension iline(2,1000000)

      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          MMAX, NB, NMAX
      PARAMETER        (MMAX=16,NB=64,NMAX=8)
      INTEGER          LDA, LWORK
      PARAMETER        (LDA=MMAX,LWORK=3*NMAX+NB*(MMAX+NMAX))
*     .. Local Scalars ..
      DOUBLE PRECISION RCOND, RNORM
      INTEGER          I, INFO, J, M, N, RANK
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), B(MMAX), S(NMAX), WORK(LWORK)
c*     .. External Functions ..
c      DOUBLE PRECISION DNRM2
c      EXTERNAL         DNRM2
c*     .. External Subroutines ..
c      EXTERNAL         DGELSS
*      
c    
c     find all lines
c
      k=0
      do i=1,nqu2
      k=k+1
      iline(1,k)=nqudr(1,i)
      iline(2,k)=nqudr(2,i)
      k=k+1
      iline(1,k)=nqudr(2,i)
      iline(2,k)=nqudr(3,i)
      k=k+1
      iline(1,k)=nqudr(3,i)
      iline(2,k)=nqudr(4,i)
      k=k+1
      iline(1,k)=nqudr(4,i)
      iline(2,k)=nqudr(1,i)
      enddo
      nline=k
      write(*,*) 'number of all lines is: ',nline
      write(*,*)
      do i=1,nline
      index(i)=0 
      if (iline(1,i).ge.iline(2,i)) then
      ittt=iline(1,i)
      iline(1,i)=iline(2,i)
      iline(2,i)=ittt
      endif
      enddo
c      do i=1,nline
c      write(*,*)i,iline(1,i),iline(2,i),index(i)
c      enddo
c    
c     find all duplicate lines
c                  
      kk=0
      do i=1,nline
      if (index(i).eq.0) then
      do j=i+1,k
      idist=(iline(1,i)-iline(1,j))*(iline(1,i)-iline(1,j))+
     *  (iline(2,i)-iline(2,j))*(iline(2,i)-iline(2,j)) 
      if (idist.eq.0) then 
      index(j)=1
      kk=kk+1
      endif
      enddo
      endif
      enddo
      nk=nline-kk
      write(*,*) 'number of all lines without duplicate is: ',nk
      write(*,*)
c      do i=1,nline
c      write(*,*)i,iline(1,i),iline(2,i),index(i)
c      enddo
c      write(*,*)

c    
c     calculate the strains on each node
c
      open(40,file='strains.txt',status='unknown',form='formatted')
      do i=1,np2
c     find all point around selected point
      k=0
      do j=1,nline
      if (index(j).eq.0) then
      if (iline(1,j).eq.i) then
      k=k+1
      ipoint(k)=iline(2,j)
      endif
      if (iline(2,j).eq.i) then
      k=k+1
      ipoint(k)=iline(1,j)
      endif
      endif
      enddo
      
c      write(*,'(10i8)') i,k,(ipoint(ii),ii=1,k)
c      write(*,*)

c     data prepare for solve overdetermined matrix
      M=k*2
      N=3
      
      do j=1,k
      dx=coor(1,i)-coor(1,ipoint(j))
      dy=coor(2,i)-coor(2,ipoint(j))
      
      du=u0(1,i+np2*(nlayer-1))-u0(1,ipoint(j)+np2*(nlayer-1))
      dv=u0(2,i+np2*(nlayer-1))-u0(2,ipoint(j)+np2*(nlayer-1))
c      write(*,'(i6,4es12.4)') i,coor(1,i),coor(2,i),u0(1,i),u0(2,i)
c      write(*,'(i6,4es12.4)') ipoint(j),coor(1,ipoint(j)),
c     *   coor(2,ipoint(j)),u0(1,ipoint(j)),u0(2,ipoint(j))
      
      a(j*2-1,1)=dx
      a(j*2-1,2)=0.d0
      a(j*2-1,3)=0.5*dy
      a(j*2,1)=0.d0
      a(j*2,2)=dy
      a(j*2,3)=0.5*dx
      b(j*2-1)=du
      b(j*2)=dv
c      write(*,'(4es12.4)') (a(j*2-1,jj),jj=1,3),b(j*2-1)
c      write(*,'(4es12.4)') (a(j*2,jj),jj=1,3),b(j*2)
c      write(*,*)
      enddo
c      do j=1,k*2
c      write(*,'(4es12.4)') (a(j,jj),jj=1,3),b(j)
c      enddo
      
      RCOND = 0.01D0
c      CALL DGELSS(M,N,1,A,LDA,B,M,S,RCOND,RANK,WORK,LWORK,INFO)
      IF (INFO.EQ.0) THEN
      WRITE (40,'(i6,4es14.4)') i,(B(j),j=1,N)
      u1(1,i+np2*(nlayer-1))=b(1)
      u1(2,i+np2*(nlayer-1))=b(2)
      u1(3,i+np2*(nlayer-1))=b(3)
         ELSE
            WRITE (40,*) 'The SVD algorithm failed to converge'
         END IF
c*
c*           Print solution
c*
c         IF (INFO.EQ.0) THEN
c*
c*           Print solution
c*
c            WRITE (NOUT,*) 'Least squares solution'
c            WRITE (NOUT,99999) (B(j),j=1,N)
c*
c*           Print the effective rank of A
c*
c            WRITE (NOUT,*)
c            WRITE (NOUT,*) 'Tolerance used to estimate the rank of A'
c            WRITE (NOUT,99998) RCOND
c            WRITE (NOUT,*) 'Estimated rank of A'
c            WRITE (NOUT,99997) RANK
c*
c*           Print singular values of A
c*
c            WRITE (NOUT,*)
c            WRITE (NOUT,*) 'Singular values of A'
c            WRITE (NOUT,99999) (S(j),j=1,N)
c*
c*           Compute and print estimate of the square root of the
c*           residual sum of squares
c*
c            IF (RANK.EQ.N) THEN
c               RNORM = DNRM2(M-N,B(N+1),1)
c               WRITE (NOUT,*)
c               WRITE (NOUT,*)
c     +           'Square root of the residual sum of squares'
c               WRITE (NOUT,99998) RNORM
c            END IF
c         ELSE
c            WRITE (NOUT,*) 'The SVD algorithm failed to converge'
c         END IF
      
      
      enddo
      close(40)
      
      
      
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
      dimension xxd(4),yyd(4),zzd(4)   
c       add by jiyang
      dimension s(6),bb(6)
      
      PARAMETER(NP=3)
      REAL dd(NP),vv(NP,NP),aa(NP,NP),ress(np),resd(np)
      
      pi=4.d0*datan(1.d0)/180.d0
      R=6378.0d3
cc      
cc      xcent=-116.00551
cc      ycent=33.3025296000000
cc      a=xcent*pi
cc      b=ycent*pi
cc      x0=100.d3
cc      y0=100.d3
cc      xmin=-39695.88935
cc      ymin=-113105.03175
cc      alpha=35.d0*pi
cc      a11=cos(alpha)
cc      a12=sin(alpha)
cc      a21=-sin(alpha)
cc      a22=cos(alpha)
      
      xx0=0.0d0
      yy0=0.0d0
      xx1=3.0d5
      yy1=4.0d5     
      step=1.0d3

      
      open(12,file='grid.txt',status='unknown',form='formatted')
c       add by jiyang
      open(122,file='beachball.txt',status='unknown',form='formatted')
      
      do lon=xx0,xx1,step
      do lat=yy0,yy1,step
c      write(*,*)lon,lat

ccc    projection
cc      xx=lon*pi
cc      yy=lat*pi
ccc      xx=120.70d0*pi
ccc      yy=35.80d0*pi
cc      call trans(R,xx,yy,a,b,c,d)
cc      cx=c
cc      cy=d
cc
cc      xxx0=a11*cx+a12*cy-xmin+x0
cc      yyy0=a21*cx+a22*cy-ymin+y0
cc      
cc      x3=xxx0
cc      y3=yyy0
      x3=lon
      y3=lat
c      write(*,'(2f10.4,2es12.4)') xx/pi,yy/pi,x3,y3
      uu01=-999
      uu02=-999
      
      do k=1,nqu2
      
      x1=coor(1,nqudr(1,k))
      y1=coor(2,nqudr(1,k))
      x2=coor(1,nqudr(2,k))
      y2=coor(2,nqudr(2,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
c      write(*,'(i,2es12.4)') 1,s1
      if (s1 .lt. 0.d0) goto 100
      
      x1=coor(1,nqudr(2,k))
      y1=coor(2,nqudr(2,k))
      x2=coor(1,nqudr(3,k))
      y2=coor(2,nqudr(3,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
c      write(*,'(i,2es12.4)') 2,s1
      if (s1 .lt. 0.d0) goto 100

      x1=coor(1,nqudr(3,k))
      y1=coor(2,nqudr(3,k))
      x2=coor(1,nqudr(4,k))
      y2=coor(2,nqudr(4,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
c      write(*,'(i,2es12.4)') 3,s1
      if (s1 .lt. 0.d0) goto 100

      x1=coor(1,nqudr(4,k))
      y1=coor(2,nqudr(4,k))
      x2=coor(1,nqudr(1,k))
      y2=coor(2,nqudr(1,k))
      s1=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)      
c      write(*,'(i,2es12.4)') 4,s1
      if (s1 .lt. 0.d0) goto 100
      
c      flag=1.0
      
c      write(*,*) k
c
c.....surfae
c
c
c..... u0
c

       xxd(1) = coor(1,nqudr(1,k))
       xxd(2) = coor(1,nqudr(2,k))
       xxd(3) = coor(1,nqudr(3,k))
       xxd(4) = coor(1,nqudr(4,k))
      
       yyd(1) = coor(2,nqudr(1,k))
       yyd(2) = coor(2,nqudr(2,k))
       yyd(3) = coor(2,nqudr(3,k))
       yyd(4) = coor(2,nqudr(4,k))
      
       zzd(1) = u0(1,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u0(1,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u0(1,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u0(1,nqudr(4,k)+np2*(nlayer-1))
      
       xxd0 = x3
       yyd0 = y3
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)       
       u01= top*1000.d0

       zzd(1) = u0(2,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u0(2,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u0(2,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u0(2,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u02= top*1000.d0
c       write(*,*)zzd,top,u02
c       pause

       zzd(1) = u0(3,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u0(3,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u0(3,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u0(3,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u03= top*1000.d0
       
       length=sqrt(u01*u01+u02*u02)
       if (abs(u01) .lt. 1e-1) then
        ang=90.0d0
       else 
        ang=atan(u02/u01)/pi
       endif
       if (ang .lt. 0.d0) ang=180.0+ang
       uu01=length*cos(ang*pi+alpha)
       uu02=length*sin(ang*pi+alpha)
       uu03=u03

c
c..... u1
c
       zzd(1) = u1(1,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(1,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(1,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(1,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u11= top
      
       zzd(1) = u1(2,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(2,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(2,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(2,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u12= top

       zzd(1) = u1(3,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(3,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(3,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(3,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u13= top

       zzd(1) = u1(4,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(4,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(4,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(4,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u14= top

       zzd(1) = u1(5,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(5,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(5,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(5,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u15= top

       zzd(1) = u1(6,nqudr(1,k)+np2*(nlayer-1))
       zzd(2) = u1(6,nqudr(2,k)+np2*(nlayer-1))
       zzd(3) = u1(6,nqudr(3,k)+np2*(nlayer-1))
       zzd(4) = u1(6,nqudr(4,k)+np2*(nlayer-1))
       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
       call inter_quadri(rr,ss,zzd,top)
       u16= top
    
c       write(*,*)lat,u02

c       write(12,1000) 1,lon,lat,uu01,uu02,u03,u14,u15,u16
c       
cc
cc.....moho
cc
cc
cc..... u0
cc
c
c       xxd(1) = coor(1,nqudr(1,k))
c       xxd(2) = coor(1,nqudr(2,k))
c       xxd(3) = coor(1,nqudr(3,k))
c       xxd(4) = coor(1,nqudr(4,k))
c      
c       yyd(1) = coor(2,nqudr(1,k))
c       yyd(2) = coor(2,nqudr(2,k))
c       yyd(3) = coor(2,nqudr(3,k))
c       yyd(4) = coor(2,nqudr(4,k))
c      
c       zzd(1) = u0(1,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u0(1,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u0(1,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u0(1,nqudr(4,k)+np2*(5-1))
c      
c       xxd0 = x3
c       yyd0 = y3
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)       
c       u01= top*1000.d0
c
c       zzd(1) = u0(2,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u0(2,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u0(2,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u0(2,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u02= top*1000.d0
cc       write(*,*)zzd,top,u02
cc       pause
c
c       zzd(1) = u0(3,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u0(3,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u0(3,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u0(3,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u03= top*1000.d0
c
c       length=sqrt(u01*u01+u02*u02)
c       if (abs(u01) .lt. 1e-1) then
c        ang=90.0d0
c       else 
c        ang=atan(u02/u01)/pi
c       endif
c       if (ang .lt. 0.d0) ang=180.0+ang
c       uu01=length*cos(ang*pi+alpha)
c       uu02=length*sin(ang*pi+alpha)
c       uu03=u03
c
cc
cc..... u1
cc
c       zzd(1) = u1(1,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(1,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(1,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(1,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u11= top
c      
c       zzd(1) = u1(2,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(2,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(2,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(2,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u12= top
c
c       zzd(1) = u1(3,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(3,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(3,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(3,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u13= top
c
c       zzd(1) = u1(4,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(4,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(4,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(4,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u14= top
c
c       zzd(1) = u1(5,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(5,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(5,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(5,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u15= top
c
c       zzd(1) = u1(6,nqudr(1,k)+np2*(5-1))
c       zzd(2) = u1(6,nqudr(2,k)+np2*(5-1))
c       zzd(3) = u1(6,nqudr(3,k)+np2*(5-1))
c       zzd(4) = u1(6,nqudr(4,k)+np2*(5-1))
c       call locate_qudr(xxd,yyd,xxd0,yyd0,rr,ss)
c       call inter_quadri(rr,ss,zzd,top)
c       u16= top
c    
cc       write(*,*)lat,u02
c
c       write(12,1000) 2,lon,lat,uu01,uu02,u03,u14,u15,u16       

c-----------------------------------------------------------------
c.....Principal stress
c-----------------------------------------------------------------
      ttt_s=(u11+u12+u13)/3.d0
      uu11=u11-ttt_s
      uu12=u12-ttt_s
      uu13=u13-ttt_s
      
      Nn=3
      do ii=1,Nn
      do jj=1,Nn
      aa(ii,jj)=0.
      enddo
      enddo
      aa(1,1)=uu11
      aa(1,2)=u16
      aa(1,3)=u15
      aa(2,1)=aa(1,2)
      aa(2,2)=uu12
      aa(2,3)=u14
      aa(3,1)=aa(1,3)
      aa(3,2)=aa(2,3)
      aa(3,3)=uu13

      call jacobi(aa,NP,NP,dd,vv,nrot)
      call eigsrt(dd,vv,NP,NP)

      do ii=1,3
c      dip(ii)=acos(v(3,ii))/pi
      resd(ii)=abs(90-abs(acos(vv(3,ii))/pi))
      if (abs(vv(2,ii)) .lt. 1.0e-3) then
      ress(ii)=90.
      else
      ress(ii)=atan(vv(1,ii)/vv(2,ii))/pi
      endif
c      ress(ii)=ress(ii)+alpha/pi
      ress(ii)=ress(ii)
      
c      write(*,*) ii,ress(ii),resd(ii),acos(0.5d0)
      enddo

      ttt_s=sqrt(u11*u11+u12*u12+2*u16*u16)
      ttt_s=log(abs(ttt_s))/log(10.d0)

      write(12,1000) 1,lon,lat,u01,u02,u03,u11,u12,u13,u14,u15,u16,
     +  dd(1),ress(1),resd(1),dd(2),ress(2),resd(2),
     +  dd(3),ress(3),resd(3),ttt_s

c       add by jiyang
      s(1)=u11
      s(2)=u12
      s(3)=u13
      s(4)=u14
      s(5)=u15
      s(6)=u16
      idx=0
      call beachpall(s,bb,idx)
      write(122,1001) lon,lat, 5., (bb(j),j=1,6), idx,lon, lat
c       write(2,9)x(k),y(k),5.,(s(kk,j),j=1,6),idx(kk),x(k),y(k)
       
      exit
100   continue

      enddo   
      
      enddo
      enddo
      close(12)
c       add by jiyang
      close(122)

c1000  format(i9,5f13.5,20es15.5)
1000  format(i9,2f13.0,3f13.5,20es15.5)
c       add by jiyang
1001  format(2f13.0,f5.0,6f8.2,i3,2f13.0)
      return
      end

c*********************add by jiyang*************start******************
c  subroutine from Yang's beach ball
        subroutine beachpall(s,bb,idx)
c        parameter (n=1000000)
        implicit real*8 (a-h,o-z)
        dimension s(6),bb(6)
        dimension t(6),ps(3),dir(3,2)

          ace=0.12

         do j=1,6
         t(j)=s(j)*1.e-6
        !print*,t(j)
         enddo

        call spstr(t,ps,dir)

        av=(ps(1)+ps(2)+ps(3))/3.
        !print*,av
        tal=(ps(1)-ps(3))/2
        ta=tal
        !print*,tal

        If(tal.lt.1.e-6) tal=1.e-6
        !print*,tal

         do j=1,3
         t(j)=(t(j)-av)/tal
         enddo

         do j=4,6
         t(j)=t(j)/tal
         enddo

        tx=tal*ace+20.

        if(tx.gt.28)tx=28.

        idx=int(tx)
        frc=exp((tx-idx)*log(10.))

        if(frc.gt.10.)then
        write(*,*)'i,tal,ta,tx,idx,frc'
        write(*,*)i,tal,ta,tx,idx,frc
        frc=10.
        endif

         do j=1,6
         t(j)=t(j)*frc
         enddo

         bb(1)=t(3)
         bb(2)=t(2)
         bb(3)=t(1)
         bb(4)=-t(4)
         bb(5)=t(5)
         bb(6)=-t(6)
c         write(2,9)t(3),t(2),t(1),-t(4),t(5),-t(6),idx,av,ta
c          write(2,9)x(k),y(k),5.,(s(kk,j),j=1,6),idx(kk),x(k),y(k)
         return
         END


c     6-stress-conponent --> principal stress and its direction
      subroutine spstr(s,sa,sb)
      implicit real*8 (a-h,o-z)
      dimension s(6),sa(3),sb(3,2)
        pi=atan(1.)*4/180
      do i=1,6
        !print*,s(i)
      s(i)=s(i)*1.e-6
        !print*,s(i)
      enddo
      a=-s(1)-s(2)-s(3)
      b=s(1)*s(2)+s(2)*s(3)+s(3)*s(1)-s(4)*s(4)-s(5)*s(5)-s(6)*s(6)
      c=s(1)*s(4)*s(4)+s(2)*s(5)*s(5)+s(3)*s(6)*s(6)-s(1)*s(2)*s(3)
     /-2*s(4)*s(5)*s(6)
      call cub(a,b,c,x1,x2,x3)
      sa(1)=x1
      sa(2)=x2
      sa(3)=x3

      do 10 i=1,3
      a1=s(1)-sa(i)
      a2=s(6)
      a3=s(5)
      a4=s(6)
      a5=s(2)-sa(i)
      a6=s(4)
c       a1*X+a2*Y+a3*Z=0
c       a4*X+a5*Y+a6*Z=0
      dt=a1*a5-a2*a4
      if(abs(dt).gt.1.e-6)then
      ax=(-a3*a5+a2*a6)/dt
      ay=(a3*a4-a1*a6)/dt
      az=1.
      else
      az=0.
      ax=-a2
      zy=-a1
      endif
        if(ax.lt.0.0)then
           ay=-ay
           az=-az
        endif
      ad=sqrt(ax*ax+ay*ay+az*az)
      if(ad.lt.1e-5)ad=1.e-5
      sb(i,2)=asin(az/ad)/pi
      sb(i,1)=0.
      if(abs(sb(i,2)-90).ge.0.1)then
      add=sqrt(ax*ax+ay*ay)
      if(add.lt.1e-5)add=1.e-5
      shad=ay/add
c      write(*,*)'shad= ',shad
      sb(i,1)=asin(shad)/pi
      endif
10    continue
        do i=1,6
        s(i)=s(i)*1.e6
        enddo
        do i=1,3
        sa(i)=sa(i)*1.e6
        enddo
      return
      end

c       cubic equation xxx+a*xx+b*x+c=0 solving
      subroutine cub(a,b,c,x1,x2,x3)
      implicit real*8 (a-h,o-z)
        pi=atan(1.)*4
      p=-a*a/3+b
      q=a*a*a*2/27-a*b/3+c
      dq=p*p*p/27+q*q/4
        if(p.ge.0.0)p=-0.000001
        r=2*sqrt(-p/3)
        pp=-q/2/sqrt(-p*p*p/27)
        if(pp.gt.1.)pp=1.
        if(pp.lt.-1.)pp=-1.
c        write(*,*)'pp= ',pp
        af=acos(pp)
        x1=r*cos(af/3)-a/3
        x2=-r*cos(af/3+pi/3)-a/3
        x3=-r*cos(af/3-pi/3)-a/3
      return
      end
       
c*********************add by jiyang*******end**************************


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

      
      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))

14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)

                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END      

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

c
c       reverse interplation subroutines for different elements
c
        subroutine locate_tri(xx,yy,x0,y0,r,s)
        implicit real*8 (a-h,o-z)
        dimension rnod(2,3),p(2),rc(2,2),cr(2,2),r0(2)
        dimension xx(3),yy(3)                    
        do 1 i=1,3
        rnod(1,i)=xx(i)
        rnod(2,i)=yy(i)
1       continue
        r0(1)=x0
        r0(2)=y0
        call locate_tri2d(3,2,rnod,r0,p,rc,cr)
        r=p(1)
        s=p(2)
        end
        
        subroutine locate_tri2d(nnode,ndm,rnod,r0,p,rc,cr)
        implicit real*8 (a-h,o-z)
	    dimension rnod(ndm,nnode),r0(ndm),p(ndm),rc(ndm,ndm),
     &            cr(ndm,ndm),r(3),dr(3)
c        real*8 err,errmax,d
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
20      continue
40      continue
c        write(*,*) 'd = ',d
        errmax = d*1.d-6/nnode
        itmax = 10
        it = 0
        do 100 i=1,ndm
        p(i)=0.0
100     continue
1       continue
        call locate_elemt_tri2d(nnode,ndm,ndm,p,r,rnod,rc,cr,det)
        do 200 i=1,ndm
        dr(i) = r(i) - r0(i)
200     continue
        err = 0.0d0
        do 300 i=1,ndm
        err = err+dr(i)**2
300     continue
c        write(*,*) 'err,dr =',err,(dr(i),i=1,ndm)
        if (err.lt.errmax .or. it.ge.itmax) goto 2
        do  i=1,ndm
        do  j=1,ndm
        p(i) = p(i) - cr(i,j)*dr(j)
        end do
        end do 
        it = it+1
        goto 1
2       continue
c        write(*,*) 'p = ',p
        return
        end

      subroutine locate_elemt_tri2d(nnode,nrefc,ncoor,refc,coor,
     *coorr,rc,cr,det)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     *            coorr(ncoor,nnode),coor(ncoor)
      call locate_telem_tri2d(refc,coor,coorr,rc)
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
         
      subroutine locate_telem_tri2d(refc,coor,coorr,rc)
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),refc(2),coor(2),coorr(2,3),rc(2,2)
      s1=refc(1)
      s2=refc(2)
      do 1 n=1,3
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
1     continue
      coor(1)=+(+s1)*x(1)+(+s2)*x(2)+(+1.-s1-s2)*x(3)
      coor(2)=+(+s1)*y(1)+(+s2)*y(2)+(+1.-s1-s2)*y(3)
      rc(1,1)=+(+1.)*x(1)+(-1.)*x(3)
      rc(2,1)=+(+1.)*y(1)+(-1.)*y(3)
      rc(1,2)=+(+1.)*x(2)+(-1.)*x(3)
      rc(2,2)=+(+1.)*y(2)+(-1.)*y(3)
      end
c
c
c
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
c       dimension rnod(2,4),r(2),p(2),rc(2,2),cr(2,2),r0(2),dr(2)
        implicit real*8 (a-h,o-z)
        dimension rnod(ndm,nnode),r0(ndm),p(ndm),rc(ndm,ndm),
     &            cr(ndm,ndm),r(3),dr(3)
c	    real*8 err,errmax,d
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
20	    continue
40	    continue
c        write(*,*) 'd = ',d
	    errmax = d*1.d-6/nnode
	    itmax = 10
	    it = 0
	    do 100 i=1,ndm
	    p(i)=0.0
100	    continue
1	    continue
	    call locate_elemt_qudr(nnode,ndm,ndm,p,r,rnod,rc,cr,det)
	    do 200 i=1,ndm
	    dr(i) = r(i) - r0(i)
200	    continue
	    err = 0.0d0
	    do 300 i=1,ndm
	    err = err+dr(i)**2
300	    continue
c        write(*,*) 'err,dr =',err,(dr(i),i=1,ndm)
	    if (err.lt.errmax .or. it.ge.itmax) goto 2
	    do 400 i=1,ndm
	    do 401 j=1,ndm
	    p(i) = p(i) - cr(i,j)*dr(j)
401	    continue
400	    continue
	    it = it+1
	    goto 1
2	    continue
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

c
c
c
c
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
c
      subroutine inter_triangle(r,s,zz,top)
      implicit real*8 (a-h,o-z)
      dimension zz(4),x(24)
      rx=r
      ry=s
      do n=1,3
      x(n)=zz(n)
      end do
           top=+(+rx)*x(1)
     *         +(+ry)*x(2)
     *         +(+(+1.-rx-ry))*x(3)
      return
      end
      


