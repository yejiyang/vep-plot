      subroutine hulq4g2(coorr,coefr,prmt,estif,emass,edamp,eload,num)
      implicit real*8 (a-h,o-z)
      dimension estif(12,12),elump(12),emass(12),
     &edamp(12),eload(12)
      dimension prmt(3000),coef(3),coefr(4,3),coorr(2,4),coor(2)
      common /rhulq4g2/ru(4,12),rv(4,12),rw(4,12),
     & cu(4,3),cv(4,3),cw(4,3)
      common /vhulq4g2/rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      common /dhulq4g2/ refc(2,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(12,3)
      if (num.eq.1) call hulq4g2i
      do 10 i=1,nvar
      emass(i)=0.0
      edamp(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do 999 igaus=1,ngaus
      call hulq4g2t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      call ehulq4g2(refc(1,igaus),coef,coorr,coefr,coefd)
      iu=(igaus-1)*3+1
      iv=(igaus-1)*3+1
      iw=(igaus-1)*3+1
      if (num.gt.1) goto 2
      call hulq4g21(refc(1,igaus),ru(1,iu),rctr,crtr)
      call hulq4g22(refc(1,igaus),rv(1,iv),rctr,crtr)
      call hulq4g23(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,4,ru(1,iu),cu,crtr,1,3,3)
      call shapn(nrefc,ncoor,4,rv(1,iv),cv,crtr,1,3,3)
      call shapn(nrefc,ncoor,4,rw(1,iw),cw,crtr,1,3,3)
      call shapc(nrefc,ncoor,3,coefd,coefc,crtr,2,5,5)
      gx=coef(1)
      gy=coef(2)
      gz=coef(3)
      weigh=det*gaus(igaus)
      fx=prmt(1)
      fy=prmt(2)
      fz=prmt(3)
      imate=prmt(6)+0.001
      do 202 i=1,4
      iv=kvord(i,3)
      do 201 j=1,4
      jv=kvord(j,3)
      stif=+cw(i,1)*cw(j,1)*0.0
      estif(jv,iv)=estif(jv,iv)+stif*weigh
201    continue
202    continue
      stif= 0.
      elump(1)=stif*weigh
      stif= 0.
      elump(4)=stif*weigh
      stif= 0.
      elump(7)=stif*weigh
      stif= 0.
      elump(10)=stif*weigh
      stif= 0.
      elump(2)=stif*weigh
      stif= 0.
      elump(5)=stif*weigh
      stif= 0.
      elump(8)=stif*weigh
      stif= 0.
      elump(11)=stif*weigh
      stif= 0.
      elump(3)=stif*weigh
      stif= 0.
      elump(6)=stif*weigh
      stif= 0.
      elump(9)=stif*weigh
      stif= 0.
      elump(12)=stif*weigh
      do 301 i=1,nnode
      if (nvard(1).lt.i) goto 301
      iv = kvord(i,1)
      emass(iv)=emass(iv)+elump(iv)*cu(i,1)
      if (nvard(2).lt.i) goto 301
      iv = kvord(i,2)
      emass(iv)=emass(iv)+elump(iv)*cv(i,1)
      if (nvard(3).lt.i) goto 301
      iv = kvord(i,3)
      emass(iv)=emass(iv)+elump(iv)*cw(i,1)
301   continue
      stif=0.
      elump(1)=stif*weigh
      stif=0.
      elump(4)=stif*weigh
      stif=0.
      elump(7)=stif*weigh
      stif=0.
      elump(10)=stif*weigh
      stif=0.
      elump(2)=stif*weigh
      stif=0.
      elump(5)=stif*weigh
      stif=0.
      elump(8)=stif*weigh
      stif=0.
      elump(11)=stif*weigh
      stif=0.
      elump(3)=stif*weigh
      stif=0.
      elump(6)=stif*weigh
      stif=0.
      elump(9)=stif*weigh
      stif=0.
      elump(12)=stif*weigh
      do 401 i=1,nnode
      if (nvard(1).lt.i) goto 401
      iv = kvord(i,1)
      edamp(iv)=edamp(iv)+elump(iv)*cu(i,1)
      if (nvard(2).lt.i) goto 401
      iv = kvord(i,2)
      edamp(iv)=edamp(iv)+elump(iv)*cv(i,1)
      if (nvard(3).lt.i) goto 401
      iv = kvord(i,3)
      edamp(iv)=edamp(iv)+elump(iv)*cw(i,1)
401   continue
      do 501 i=1,4
      iv=kvord(i,1)
      stif=+cu(i,1)*fx
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,4
      iv=kvord(i,2)
      stif=+cv(i,1)*fy
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,4
      iv=kvord(i,3)
      stif=+cw(i,1)*fz
      eload(iv)=eload(iv)+stif*weigh
503   continue
999   continue
      return
      end

      subroutine hulq4g2i
      implicit real*8 (a-h,o-z)
      common /dhulq4g2/ refc(2,4),gaus(4),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(12,3)
      ngaus=  4
      ndisp=  3
      nrefc=  2
      ncoor=  2
      nvar = 12
      nnode=  4
      kdord(1)=1
      nvard(1)=4
      kvord(1,1)=1
      kvord(2,1)=4
      kvord(3,1)=7
      kvord(4,1)=10
      kdord(2)=1
      nvard(2)=4
      kvord(1,2)=2
      kvord(2,2)=5
      kvord(3,2)=8
      kvord(4,2)=11
      kdord(3)=1
      nvard(3)=4
      kvord(1,3)=3
      kvord(2,3)=6
      kvord(3,3)=9
      kvord(4,3)=12
      refc(1,1)=5.77350e-001
      refc(2,1)=5.77350e-001
      gaus(1)=1.00000e+000
      refc(1,2)=5.77350e-001
      refc(2,2)=-5.77350e-001
      gaus(2)=1.00000e+000
      refc(1,3)=-5.77350e-001
      refc(2,3)=5.77350e-001
      gaus(3)=1.00000e+000
      refc(1,4)=-5.77350e-001
      refc(2,4)=-5.77350e-001
      gaus(4)=1.00000e+000
      end


      subroutine hulq4g2t(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call thulq4g2(refc,coor,coorr,coefr,rc)
      n=nrefc
      m=n*2
      det = 1.0
      do 10 i=1,n
      do 10 j=1,n
      if (i.le.ncoor) a(i,j) = rc(i,j)
      if (i.gt.ncoor) a(i,j)=1.0
      a(i,n+j)=0.0
      if (i.eq.j) a(i,n+i) = 1.0
10    continue
c     write(*,*) 'a ='
c     do 21 i=1,n
c21   write(*,8) (a(i,j),j=1,m)
      do 400 i=1,n
      amax = 0.0
      l = 0
      do 50 j=i,n
      c = a(j,i)
      if (c.lt.0.0) c = -c
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
100   a(i,k) = a(i,k)/c
      do 300 j=1,n
      if (i.eq.j) goto 300
      do 200 k=i+1,m
200   a(j,k) = a(j,k)-a(i,k)*a(j,i)
c     write(*,*) 'i =',i,'  j =',j,'  a ='
c     do 11 ii=1,n
c11   write(*,8) (a(ii,jj),jj=1,m)
300   continue
400   continue
      do 500 i=1,nrefc
      do 500 j=1,ncoor
500   cr(i,j) = a(i,n+j)
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

      subroutine hulq4g21(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(4,3),rctr(2,2),crtr(2,2)
      external fhulq4g21
      rx=refc(1)
      ry=refc(2)
      call dshap(fhulq4g21,refc,shpr,2,4,1)
      return
      end

      real*8 function fhulq4g21(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchulq4g2/ xa(4),ya(4),gxa(4),gya(4),
     &gza(4)
      common /vhulq4g2/ rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      rx=refc(1)
      ry=refc(2)
      goto (1,2,3,4) n
1     fhulq4g21=+(+1.-rx)/2.*(+1.-ry)/2. 
      goto 1000
2     fhulq4g21=+(+1.+rx)/2.*(+1.-ry)/2. 
      goto 1000
3     fhulq4g21=+(+1.+rx)/2.*(+1.+ry)/2. 
      goto 1000
4     fhulq4g21=+(+1.-rx)/2.*(+1.+ry)/2. 
      goto 1000
1000  return
      end

      subroutine hulq4g22(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(4,3),rctr(2,2),crtr(2,2)
      external fhulq4g22
      rx=refc(1)
      ry=refc(2)
      call dshap(fhulq4g22,refc,shpr,2,4,1)
      return
      end

      real*8 function fhulq4g22(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchulq4g2/ xa(4),ya(4),gxa(4),gya(4),
     &gza(4)
      common /vhulq4g2/ rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      rx=refc(1)
      ry=refc(2)
      goto (1,2,3,4) n
1     fhulq4g22=+(+1.-rx)/2.*(+1.-ry)/2. 
      goto 1000
2     fhulq4g22=+(+1.+rx)/2.*(+1.-ry)/2. 
      goto 1000
3     fhulq4g22=+(+1.+rx)/2.*(+1.+ry)/2. 
      goto 1000
4     fhulq4g22=+(+1.-rx)/2.*(+1.+ry)/2. 
      goto 1000
1000  return
      end

      subroutine hulq4g23(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(2),shpr(4,3),rctr(2,2),crtr(2,2)
      external fhulq4g23
      rx=refc(1)
      ry=refc(2)
      call dshap(fhulq4g23,refc,shpr,2,4,1)
      return
      end

      real*8 function fhulq4g23(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchulq4g2/ xa(4),ya(4),gxa(4),gya(4),
     &gza(4)
      common /vhulq4g2/ rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      dimension refc(2)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      rx=refc(1)
      ry=refc(2)
      goto (1,2,3,4) n
1     fhulq4g23=+(+1.-rx)/2.*(+1.-ry)/2. 
      goto 1000
2     fhulq4g23=+(+1.+rx)/2.*(+1.-ry)/2. 
      goto 1000
3     fhulq4g23=+(+1.+rx)/2.*(+1.+ry)/2. 
      goto 1000
4     fhulq4g23=+(+1.-rx)/2.*(+1.+ry)/2. 
      goto 1000
1000  return
      end

      subroutine thulq4g2(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(2),coor(2),coorr(2,4),coefr(4,3),rc(2,2)
      common /cchulq4g2/ x(4),y(4),gx(4),gy(4),gz(4)
      external fthulq4g2
      do 100 n=1,4
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
100   continue
      do 200 n=1,4
      gx(n)=coefr(n,1)
      gy(n)=coefr(n,2)
      gz(n)=coefr(n,3)
200   continue
      rx=refc(1)
      ry=refc(2)
      call dcoor(fthulq4g2,refc,coor,rc,2,2,1)
      return
      end

      real*8 function fthulq4g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(2)
      common /cchulq4g2/ x(4),y(4),gx(4),gy(4),gz(4)
      common /vhulq4g2/ rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      rx=refc(1)
      ry=refc(2)
      goto (1,2) n
1     fthulq4g2=+(+(+1.-rx)/2.*(+1.-ry)/2.)*x(1)+(+(+1.+rx)/
     & 2.*(+1.-ry)/2.)*x(2)+(+(+1.+rx)/2.*(+1.+ry)/2.)*x(3)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.)*x(4)
      goto 1000
2     fthulq4g2=+(+(+1.-rx)/2.*(+1.-ry)/2.)*y(1)+(+(+1.+rx)/
     & 2.*(+1.-ry)/2.)*y(2)+(+(+1.+rx)/2.*(+1.+ry)/2.)*y(3)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.)*y(4)
      goto 1000
1000  return
      end

      subroutine ehulq4g2(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(2),coef(3),coorr(2,4),coefr(4,3),coefd(3,2)
      external fehulq4g2
      rx=refc(1)
      ry=refc(2)
      call dcoef(fehulq4g2,refc,coef,coefd,2,3,2)
      return
      end

      real*8 function fehulq4g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(2)
      common /cchulq4g2/ xa(4),ya(4),gx(4),gy(4),gz(4)
      common /vhulq4g2/ rctr(2,2),crtr(2,2),coefd(3,5),coefc(3,5)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      rx=refc(1)
      ry=refc(2)
      goto (1,2,3) n
1     fehulq4g2=+(+(+1.-rx)/2.*(+1.-ry)/2.)*gx(1)+(+(+1.+rx)/
     & 2.*(+1.-ry)/2.)*gx(2)+(+(+1.+rx)/2.*(+1.+ry)/2.)*gx(3)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.)*gx(4)
      goto 1000
2     fehulq4g2=+(+(+1.-rx)/2.*(+1.-ry)/2.)*gy(1)+(+(+1.+rx)/
     & 2.*(+1.-ry)/2.)*gy(2)+(+(+1.+rx)/2.*(+1.+ry)/2.)*gy(3)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.)*gy(4)
      goto 1000
3     fehulq4g2=+(+(+1.-rx)/2.*(+1.-ry)/2.)*gz(1)+(+(+1.+rx)/
     & 2.*(+1.-ry)/2.)*gz(2)+(+(+1.+rx)/2.*(+1.+ry)/2.)*gz(3)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.)*gz(4)
      goto 1000
1000  return
      end


