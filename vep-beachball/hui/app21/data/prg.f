      program d3xyz
      open(10,file='coor0',form='unformatted',status='unknown')
      call coor0
      close(10)
      open(10,file='id0',form='unformatted',status='unknown')
      call id0
      close(10)
      open(10,file='disp0',form='unformatted',status='unknown')
      call disp0
      close(10)
      open(10,file='tli0',form='unformatted',status='unknown')
      call tli0
      close(10)
      open(10,file='elem0',form='unformatted',status='unknown')
      call w6m0
      call mate0
      call lq4g20
      call bf0
      close(10)
      end
      subroutine coor0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /coor / x(100000),y(100000),z(100000)
      integer n
      open(121,file='data',form='formatted',status='old')
      read(121,*) nnode
      do i=1,nnode
      read(121,*) ii,xt(i),yt(i)
      enddo

c      do i=1,nnode
c      if(yt(i).lt.100) yt(i)=-200000.d0
c      if(yt(i).gt.1649900) yt(i)=1850000.d0
c      enddo

      read(121,*) ntri
      do i=1,ntri
      read(121,*) ii,itri(1,i),itri(2,i),itri(3,i),itri(4,i),
     *mat(i)
      enddo
      close(121)
      open(122,file='prg.dat',form='formatted',status='old')
      read(122,*) stif1
      read(122,*) stif2
      do i=1,6
      read(122,*) pee(i),pvv(i),yitaa(i)
      enddo
      close(122)
      nz=7
      dz=60.d3/(nz-1.d0)
      n=0
      do  2 j=1,nz
      do  1 i=1,nnode
      n=n+1
      x(n)=xt(i)
      y(n)=yt(i)
      z(n)=0.d0+(j-1)*dz
   1  continue
   2  continue
c      do i=1,nnode*5
c      if (y(i).ge.490000) y(i)=y(i)+200000
c      if (y(i).ge.470000) y(i)=y(i)+100000
c      if (y(i).ge.30000) y(i)=y(i)+100000
c      if (y(i).ge.10000) y(i)=y(i)+200000
c      enddo
      n=nnode*nz
      nn=nnode*nz
      mtij=n
      mtik=  3
      write(10) mtij,mtik,
     *(x(mtii),y(mtii),z(mtii),mtii=1,n)
      return
      end
 
      subroutine id0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /id / idu(100000),idv(100000),idw(100000)
      common /coor / x(100000),y(100000),z(100000)
      integer n,idu,idv,idw
      do  5 i=1,nn
      n=i
c      if (n.le.0) goto  5
      idu(n)=1
      idv(n)=1
      idw(n)=1
   5  continue
c      do  6 i=1,nnode
c      n=i
c      if (n.le.0) goto  6
c      idu(n)=1
c      idv(n)=1
c      idw(n)=-1
c   6  continue

      open(11,file='boundaries',status='old')
      mbb=1000
c         write(*,*)'surrounding sides',mbb
      do i=1,mbb
      read(11,*,end=37)ia,ib,ic,idx,idy,idm
c        write(*,*)i,ic,idx,idy,idm
        do j=1,ic
        read(11,*)ia,k1,k2
         if(idx.eq.-1) idu(k1)=-1
         if(idy.eq.-1) idv(k1)=-1
         if(idm.eq.-1) idw(k1)=-1
         if(idx.eq.-1) idu(k2)=-1
         if(idy.eq.-1) idv(k2)=-1
         if(idm.eq.-1) idw(k2)=-1
c         write(*,*)'j ',k1,idu(k1),idv(k1),idw(k1)
        enddo
      enddo
37      close(11)

      do iz=1,nz
      do i=1,nnode
      k=i+(iz-1)*nnode
c      idw(k)=1
      idu(k)=idu(i)
      idv(k)=idv(i)
      idw(k)=idw(i)
      enddo
      enddo
      do iz=1,nz
      do i=1,nnode
      k=i+(iz-1)*nnode
      if(iz.eq.1)idw(k)=-1
      enddo
      enddo

      n=nn
      mtij=n
      mtik=  3
      write(10) mtij,mtik,
     *(idu(mtii),idv(mtii),idw(mtii),mtii=1,n)
      return
      end
 
      subroutine disp0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /disp / u(100000),v(100000),w(100000)
      integer n
      common /coor/x(100000),y(100000),z(100000)

      dimension xtt(100000),ytt(100000),r(3)
      dimension idu(100000),idv(100000),idw(100000)
      common/bdline/nop,Lin,p(1000,4),Ln(1000,2)
            
      do 22 i=1,nn
      n=i
c      if (n.le.0) goto 22
      u(n)=0.0
      v(n)=0.0
      w(n)=0.0
  22  continue

      OPEN (11,file='id0',form='unformatted',status='unknown')
      read(11) mtij,mtik,
     *(idu(i),idv(i),idw(i),i=1,mtij)
      close(11)

        open(11,file='mesh.cor',status='old')
        read(11,*)np
        do i=1,np
        read(11,*)xtt(i),ytt(i)
        enddo
        close (11)

        open(11,file='t2.dat',status='old')
        read(11,*)nop
        do i=1,nop
        read(11,*)m,p(i,1),p(i,2)
        p(i,3)=0.
        p(i,4)=0.
        enddo
        read(11,*)num
        do i=1,num
        read(11,*)
        enddo
        read(11,*)Lin
        do i=1,Lin
        read(11,*)m,Ln(i,1),Ln(i,2)
        enddo
        close (11)
        
        open(11,file='gps.dat',form='formatted',status='old')
        do i=1,1000
        read(11,*,end=9988)k,p(k,3),p(k,4)
        enddo
        close (11)
9988    continue

c        write(*,*)xtt(1),ytt(1),p(1,1),p(1,2),p(1,3),p(1,4) 

        do i=1,m0
         u(i)=0.0
         v(i)=0.0
         w(i)=0.0
        enddo       
        open(31,file='disp.dat',form='formatted',status='unknown')
        do iz=1,nz
        do i=1,nnode
         k=i+(iz-1)*nnode
         r(1)=xtt(i)
         r(2)=ytt(i)
         if (idu(k).eq.-1) u(k)=bound(r,t,1)
         if (idv(k).eq.-1) v(k)=bound(r,t,2)
         w(k)=0.0
         if(iz.eq.nz) 
     *    write(31,'(i6,3es12.2)')i,u(k),v(k),sqrt(u(k)*u(k)+v(k)*v(k))
        enddo
        enddo
        close(31)

      n=nn
      mtij=n
      mtik=  3
      write(10) mtij,mtik,
     *(u(mtii),v(mtii),w(mtii),mtii=1,n)
      return
      end
 
      subroutine tli0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /tli / dxx(100000),dyy(100000),dzz(100000),
     *dyz(100000),dxz(100000),dxy(100000)
      integer n
      do 23 i=1,nn
      n=i
      if (n.le.0) goto 23
      dxx(n)=0.0
      dyy(n)=0.0
      dzz(n)=0.0
      dyz(n)=0.0
      dxz(n)=0.0
      dxy(n)=0.0
  23  continue
      mtij=n
      mtik=  6
      write(10) mtij,mtik,
     *(dxx(mtii),dyy(mtii),dzz(mtii),dyz(mtii),
     *dxz(mtii),dxy(mtii),mtii=1,n)
      return
      end
 
      subroutine w6m0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /w6m / nod1(100000),nod2(100000),nod3(100000),
     *nod4(100000),nod5(100000),nod6(100000),nod7(100000),
     *nod8(100000),nod9(100000),nod10(100000),nod11(100000),
     *nod12(100000),nod13(100000),nod14(100000),nod15(100000),
     *nod16(100000),nod17(100000),nod18(100000),nod19(100000),
     *nod20(100000),nod21(100000),nod22(100000),nod23(100000),
     *nod24(100000),nod25(100000),nod26(100000),nod27(100000),
     *mate(100000)
      integer n,nod1,nod2,nod3,nod4,nod5,
     *nod6,nod7,nod8,nod9,nod10,nod11,nod12,nod13,nod14,nod15,nod16,
     *nod17,nod18,nod19,nod20,nod21,nod22,nod23,nod24,nod25,nod26,
     *nod27,mate
c      do i=1,ntri
c      if (mat(i).eq.2) mat(i)=1
c      if (mat(i).eq.3) mat(i)=2
c      enddo
      do 25 j=1,4
      do 24 i=1,ntri
      n=(j-1)*ntri+i
c      if (n.le.0) goto 24
      nod1(n)=itri(1,i)+nnode*(j-1)
      nod2(n)=itri(2,i)+nnode*(j-1)
      nod3(n)=itri(3,i)+nnode*(j-1)
      nod4(n)=itri(4,i)+nnode*(j-1)
      nod5(n)=itri(1,i)+nnode*j
      nod6(n)=itri(2,i)+nnode*j
      nod7(n)=itri(3,i)+nnode*j
      nod8(n)=itri(4,i)+nnode*j
      mate(n)=mat(i)+3
  24  continue
  25  continue
      do 27 j=5,6
      do 26 i=1,ntri
      n=(j-1)*ntri+i
c      if (n.le.0) goto 26
      nod1(n)=itri(1,i)+nnode*(j-1)
      nod2(n)=itri(2,i)+nnode*(j-1)
      nod3(n)=itri(3,i)+nnode*(j-1)
      nod4(n)=itri(4,i)+nnode*(j-1)
      nod5(n)=itri(1,i)+nnode*j
      nod6(n)=itri(2,i)+nnode*j
      nod7(n)=itri(3,i)+nnode*j
      nod8(n)=itri(4,i)+nnode*j
      mate(n)=mat(i)
  26  continue
  27  continue
      n=ntri*(nz-1)
      mtij=n
      mtik=  9
      write(10) mtij,mtik,
     *(nod1(mtii),nod2(mtii),nod3(mtii),nod4(mtii),
     *nod5(mtii),nod6(mtii),nod7(mtii),nod8(mtii),
     *mate(mtii),mtii=1,n)
      return
      end
 
      subroutine mate0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /mate / pe(100000),pv(100000),fx(100000),
     *fy(100000),fz(100000),rou(100000),alpha(100000),
     *yita(100000)
      integer n
      do 28 i=1,7
      n=i
      if (n.le.0) goto 28
      pe(n)=pee(i)
      pv(n)=pvv(i)
      fx(n)=0.
      fy(n)=0.
      fz(n)=0.
      rou(n)=0.
      alpha(n)=0.
      yita(n)=yitaa(i)
  28  continue
      do 29 i=1,5
      n=i
      if (n.le.0) goto 29
      pe(n)=pee(i)
      pv(n)=pvv(i)
      fx(n)=0.
      fy(n)=0.
      fz(n)=0.
      rou(n)=10.e6
      alpha(n)=0.
      yita(n)=yitaa(i)
  29  continue
      do 30 i=6,15
      n=i
      if (n.le.0) goto 30
      pe(n)=pee(i)
      pv(n)=pvv(i)
      fx(n)=0.
      fy(n)=0.
      fz(n)=0.
      rou(n)=0.
      alpha(n)=0.
      yita(n)=yitaa(i)
  30  continue
      mtij=n
      mtik=  8
      write(10) mtij,mtik,
     *(pe(mtii),pv(mtii),fx(mtii),fy(mtii),
     *fz(mtii),rou(mtii),alpha(mtii),yita(mtii)
     *,mtii=1,n)
      return
      end
 
      subroutine lq4g20
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /lq4g2 / nod1(100000),nod2(100000),nod3(100000),
     *nod4(100000),mate(100000)
      integer n,nod1,nod2,nod3,nod4,mate
      do 32 j=1,1
      do 31 i=1,1
      n=(j-1)*nz+i
      if (n.le.0) goto 31
      nod1(n)=itri(1,i)
      nod2(n)=itri(2,i)
      nod3(n)=itri(3,i)
      nod4(n)=itri(4,i)
      mate(n)=1
  31  continue
  32  continue
      mtij=n
      mtik=  5
      write(10) mtij,mtik,
     *(nod1(mtii),nod2(mtii),nod3(mtii),nod4(mtii),
     *mate(mtii),mtii=1,n)
      return
      end
 
      subroutine bf0
      implicit real*8 (a-h,o-z)
      common /nxy/ nz,nn,nnode,ntri,itri(9,100000),mat(100000),nnode1
      common /xy/ xt(100000),yt(100000)      
      common /mat/ pee(100000),pvv(100000),yitaa(100000),stiff(100000)
      common /xy1/ press,faz(1000),fax(1000),fay(1000)
      common /bf / fx(100000),fy(100000),fz(100000)
      integer n
      n=1
      fx(n)=0.0
      fy(n)=0.0
      fz(n)=0.0
      mtij=n
      mtik=  3
      write(10) mtij,mtik,
     *(fx(mtii),fy(mtii),fz(mtii),mtii=1,n)
      return
      end

      real*8 function bound(r,t,j)
      implicit real*8 (a-h,o-z)
       common/bdline/nop,Lin,p(1000,4),Ln(1000,2)

c      implicit  (a-h,o-z)
      dimension r(2)
      bound=0.0

      do i=1,Lin
      k1=Ln(i,1)
      k2=Ln(i,2)
      call proc(p(k1,1),p(k1,2),p(k1,3),p(k1,4),
     # p(k2,1),p(k2,2),p(k2,3),p(k2,4),r(1),r(2),u3,v3,L)
      if(L.gt.0)then
      if(j.eq.1)bound=u3
      if(j.eq.2)bound=v3
      return
      endif
      enddo
c     write(*,*) 'bound =',bound
      return
      end
 
        subroutine proc(x1,y1,u1,v1,x2,y2,u2,v2,x,y,u3,v3,L)
        implicit real*8 (a-h,o-z)
        u3=0.
        v3=0.
        L=-1
        a1=x2-x1
        a2=y2-y1
        b1=x-x1
        b2=y-y1
        c1=x-x2
        c2=y-y2
        d=sqrt(a1*a1+a2*a2)
        s=abs(a1*b2-a2*b1)
        h=2*s/d
        f=b1*c1+b2*c2

        if(h.lt.0.05.and.f.le.0)then
        L=1
        d1=sqrt(b1*b1+b2*b2)
        r=d1/d
        u3=(1-r)*u1+r*u2
        v3=(1-r)*v1+r*v2
        return
        endif
        return
        end
