      subroutine huw6m(coorr,coefr,prmt,estif,emass,edamp,eload,num
     &,fstr6,label6)
      implicit real*8 (a-h,o-z)
      dimension estif(18,18),elump(18),emass(18),fstr6(6,6),
     &edamp(18),eload(18),label6(6)
      dimension prmt(3000),coef(9),coefr(6,9),
     & efuna(18),efunb(18),efunc(18),efund(18),
     & efune(18),efunf(18),coorr(3,6),coor(3)
      common /rhuw6m/ru(6,24),rv(6,24),rw(6,24),
     & cu(6,4),cv(6,4),cw(6,4)
      common /vhuw6m/rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      common /dhuw6m/ refc(3,6),gaus(6),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(18,3)
      dimension d(6,6),f(6),stress(7),strain(6)
      dimension dp(6,6),p(4),dv(6),ddv(6)
      external prager
      pe = prmt(1)
      pv = prmt(2)
      fu = prmt(3)
      fv = prmt(4)
      fw = prmt(5)
      rou = prmt(6)
      alpha = prmt(7)
      yita = prmt(8)
      yita1 = prmt(8)
      time = prmt(9)
      dt = prmt(10)
      imate = prmt(11)+0.01
      if (num.eq.1) call huw6mi
      do 10 i=1,nvar
      emass(i)=0.0d0
      edamp(i)=0.0d0
      eload(i)=0.0d0
      do 10 j=1,nvar
      estif(i,j)=0.0d0
10    continue
      do 999 igaus=1,ngaus
      call huw6mt(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call ehuw6m(refc(1,igaus),coef,coorr,coefr,coefd)
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call huw6m1(refc(1,igaus),ru(1,iu),rctr,crtr)
      call huw6m2(refc(1,igaus),rv(1,iv),rctr,crtr)
      call huw6m3(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,6,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rw(1,iw),cw,crtr,1,4,4)
      call shapc(nrefc,ncoor,9,coefd,coefc,crtr,2,9,9)
      un=coef(1)
      vn=coef(2)
      wn=coef(3)
      fxxn=coef(4)
      fyyn=coef(5)
      fzzn=coef(6)
      fyzn=coef(7)
      fxzn=coef(8)
      fxyn=coef(9)
      weigh=det*gaus(igaus)
      do 100 i=1,18
      efuna(i) = 0.0d0
      efunb(i) = 0.0d0
      efunc(i) = 0.0d0
      efund(i) = 0.0d0
      efune(i) = 0.0d0
      efunf(i) = 0.0d0
100   continue
      do 101 i=1,6
      iv=kvord(i,1)
      stif=+cu(i,2) 
      efuna(iv)=efuna(iv)+stif
101   continue
      do 102 i=1,6
      iv=kvord(i,2)
      stif=+cv(i,3) 
      efunb(iv)=efunb(iv)+stif
102   continue
      do 103 i=1,6
      iv=kvord(i,3)
      stif=+cw(i,4) 
      efunc(iv)=efunc(iv)+stif
103   continue
      do 104 i=1,6
      iv=kvord(i,2)
      stif=+cv(i,4) 
      efund(iv)=efund(iv)+stif
104   continue
      do 105 i=1,6
      iv=kvord(i,3)
      stif=+cw(i,3) 
      efund(iv)=efund(iv)+stif
105   continue
      do 106 i=1,6
      iv=kvord(i,1)
      stif=+cu(i,4) 
      efune(iv)=efune(iv)+stif
106   continue
      do 107 i=1,6
      iv=kvord(i,3)
      stif=+cw(i,2) 
      efune(iv)=efune(iv)+stif
107   continue
      do 108 i=1,6
      iv=kvord(i,1)
      stif=+cu(i,3) 
      efunf(iv)=efunf(iv)+stif
108   continue
      do 109 i=1,6
      iv=kvord(i,2)
      stif=+cv(i,2) 
      efunf(iv)=efunf(iv)+stif
109   continue
      fxx=fstr6(1,igaus)
      fyy=fstr6(2,igaus)
      fzz=fstr6(3,igaus)
      fyz=fstr6(4,igaus)
      fxz=fstr6(5,igaus)
      fxy=fstr6(6,igaus)
      stress(1)=fxx
      stress(2)=fyy
      stress(3)=fzz
      stress(4)=fyz
      stress(5)=fxz
      stress(6)=fxy
      stress(7)=0.d0
      strain(1)=+coefc(1,1)
      strain(2)=+coefc(2,2)
      strain(3)=+coefc(3,3)
      strain(4)=+coefc(2,3)+coefc(3,2)
      strain(5)=+coefc(3,1)+coefc(1,3)
      strain(6)=+coefc(1,2)+coefc(2,1)
      yita=yita1
      call getdfu(pe,pv,yita,dt,stress,strain,d,f)
      do i=1,6
      ddv(i)=f(i)-stress(i)
      enddo
      if (imate.ge.2.and.imate.le.4) then
      iimate=1
      p(1)=0.d0
      p(2)=10.0d6
      p(3)=prmt(2)
      p(4)=0.0d0
c      if (label6(igaus).eq.1) p(2)=1.0d4
      else
      iimate=2
      p(1)=0.d0
      p(2)=50.0d6
      p(3)=prmt(2)
      p(4)=0.0d0
      endif
      call getstr(p,stress,strain,dv,ddv,d,dp,1,prag,prager)
      if (prag.lt.-p(2)*0.0001) goto 1100
      label6(igaus)=1
      do 1000 i=1,6
      do 1000 j=1,6
      d(i,j)=d(i,j)-dp(i,j)
1000   continue
1100  continue
      do 202 iv=1,18
      do 201 jv=1,18
      stif=+efuna(iv)*efuna(jv)*d(1,1)
     & +efuna(iv)*efunb(jv)*d(1,2)
     & +efuna(iv)*efunc(jv)*d(1,3)
     & +efuna(iv)*efund(jv)*d(1,4)
     & +efuna(iv)*efune(jv)*d(1,5)
     & +efuna(iv)*efunf(jv)*d(1,6)
     & +efunb(iv)*efuna(jv)*d(2,1)
     & +efunb(iv)*efunb(jv)*d(2,2)
     & +efunb(iv)*efunc(jv)*d(2,3)
     & +efunb(iv)*efund(jv)*d(2,4)
     & +efunb(iv)*efune(jv)*d(2,5)
     & +efunb(iv)*efunf(jv)*d(2,6)
     & +efunc(iv)*efuna(jv)*d(3,1)
     & +efunc(iv)*efunb(jv)*d(3,2)
     & +efunc(iv)*efunc(jv)*d(3,3)
     & +efunc(iv)*efund(jv)*d(3,4)
     & +efunc(iv)*efune(jv)*d(3,5)
     & +efunc(iv)*efunf(jv)*d(3,6)
     & +efund(iv)*efuna(jv)*d(4,1)
     & +efund(iv)*efunb(jv)*d(4,2)
     & +efund(iv)*efunc(jv)*d(4,3)
     & +efund(iv)*efund(jv)*d(4,4)
     & +efund(iv)*efune(jv)*d(4,5)
     & +efund(iv)*efunf(jv)*d(4,6)
     & +efune(iv)*efuna(jv)*d(5,1)
     & +efune(iv)*efunb(jv)*d(5,2)
     & +efune(iv)*efunc(jv)*d(5,3)
     & +efune(iv)*efund(jv)*d(5,4)
     & +efune(iv)*efune(jv)*d(5,5)
     & +efune(iv)*efunf(jv)*d(5,6)
     & +efunf(iv)*efuna(jv)*d(6,1)
     & +efunf(iv)*efunb(jv)*d(6,2)
     & +efunf(iv)*efunc(jv)*d(6,3)
     & +efunf(iv)*efund(jv)*d(6,4)
     & +efunf(iv)*efune(jv)*d(6,5)
     & +efunf(iv)*efunf(jv)*d(6,6)
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
      stif=0.d0
      elump(1)=stif*weigh
      stif=0.d0
      elump(4)=stif*weigh
      stif=0.d0
      elump(7)=stif*weigh
      stif=0.d0
      elump(10)=stif*weigh
      stif=0.d0
      elump(13)=stif*weigh
      stif=0.d0
      elump(16)=stif*weigh
      stif=0.d0
      elump(2)=stif*weigh
      stif=0.d0
      elump(5)=stif*weigh
      stif=0.d0
      elump(8)=stif*weigh
      stif=0.d0
      elump(11)=stif*weigh
      stif=0.d0
      elump(14)=stif*weigh
      stif=0.d0
      elump(17)=stif*weigh
      stif=0.d0
      elump(3)=stif*weigh
      stif=0.d0
      elump(6)=stif*weigh
      stif=0.d0
      elump(9)=stif*weigh
      stif=0.d0
      elump(12)=stif*weigh
      stif=0.d0
      elump(15)=stif*weigh
      stif=0.d0
      elump(18)=stif*weigh
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
      stif=0.d0*alpha
      elump(1)=stif*weigh
      stif=0.d0*alpha
      elump(4)=stif*weigh
      stif=0.d0*alpha
      elump(7)=stif*weigh
      stif=0.d0*alpha
      elump(10)=stif*weigh
      stif=0.d0*alpha
      elump(13)=stif*weigh
      stif=0.d0*alpha
      elump(16)=stif*weigh
      stif=0.d0*alpha
      elump(2)=stif*weigh
      stif=0.d0*alpha
      elump(5)=stif*weigh
      stif=0.d0*alpha
      elump(8)=stif*weigh
      stif=0.d0*alpha
      elump(11)=stif*weigh
      stif=0.d0*alpha
      elump(14)=stif*weigh
      stif=0.d0*alpha
      elump(17)=stif*weigh
      stif=0.d0*alpha
      elump(3)=stif*weigh
      stif=0.d0*alpha
      elump(6)=stif*weigh
      stif=0.d0*alpha
      elump(9)=stif*weigh
      stif=0.d0*alpha
      elump(12)=stif*weigh
      stif=0.d0*alpha
      elump(15)=stif*weigh
      stif=0.d0*alpha
      elump(18)=stif*weigh
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
      do 501 i=1,6
      iv=kvord(i,1)
      stif=+cu(i,1)*fu*0.0d0
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,6
      iv=kvord(i,2)
      stif=+cv(i,1)*fv*0.0d0
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,6
      iv=kvord(i,3)
      stif=-cw(i,1)*fw*0.0d0
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 iv=1,18
      stif=-efuna(iv)*(stress(1)+dv(1))
     & -efunb(iv)*(stress(2)+dv(2))
     & -efunc(iv)*(stress(3)+dv(3))
     & -efund(iv)*(stress(4)+dv(4))
     & -efune(iv)*(stress(5)+dv(5))
     & -efunf(iv)*(stress(6)+dv(6))
      eload(iv)=eload(iv)+stif*weigh
504   continue
999   continue
      return
      end

      subroutine huw6mi
      implicit real*8 (a-h,o-z)
      common /dhuw6m/ refc(3,6),gaus(6),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(18,3)
      ngaus=  6
      ndisp=  3
      nrefc=  3
      ncoor=  3
      nvar = 18
      nnode=  6
      kdord(1)=1
      nvard(1)=6
      kvord(1,1)=1
      kvord(2,1)=4
      kvord(3,1)=7
      kvord(4,1)=10
      kvord(5,1)=13
      kvord(6,1)=16
      kdord(2)=1
      nvard(2)=6
      kvord(1,2)=2
      kvord(2,2)=5
      kvord(3,2)=8
      kvord(4,2)=11
      kvord(5,2)=14
      kvord(6,2)=17
      kdord(3)=1
      nvard(3)=6
      kvord(1,3)=3
      kvord(2,3)=6
      kvord(3,3)=9
      kvord(4,3)=12
      kvord(5,3)=15
      kvord(6,3)=18
      refc(1,1)=0.000000d+00
      refc(2,1)=0.000000d+00
      refc(3,1)=-1.000000d+00
      gaus(1)=1.666666d-01
      refc(1,2)=0.000000d+00
      refc(2,2)=1.000000d+00
      refc(3,2)=-1.000000d+00
      gaus(2)=1.666666d-01
      refc(1,3)=1.000000d+00
      refc(2,3)=0.000000d+00
      refc(3,3)=-1.000000d+00
      gaus(3)=1.666666d-01
      refc(1,4)=0.000000d+00
      refc(2,4)=0.000000d+00
      refc(3,4)=1.000000d+00
      gaus(4)=1.666666d-01
      refc(1,5)=0.000000d+00
      refc(2,5)=1.000000d+00
      refc(3,5)=1.000000d+00
      gaus(5)=1.666666d-01
      refc(1,6)=1.000000d+00
      refc(2,6)=0.000000d+00
      refc(3,6)=1.000000d+00
      gaus(6)=1.666666d-01
      end


      subroutine huw6mt(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call thuw6m(refc,coor,coorr,coefr,rc)
      n=nrefc
      m=n*2
      det = 1.0d0
      do 10 i=1,n
      do 10 j=1,n
      if (i.le.ncoor) a(i,j) = rc(i,j)
      if (i.gt.ncoor) a(i,j)=1.0d0
      a(i,n+j)=0.0d0
      if (i.eq.j) a(i,n+i) = 1.0d0
10    continue
c     write(*,*) 'a ='
c     do 21 i=1,n
c21   write(*,8) (a(i,j),j=1,m)
      do 400 i=1,n
      amax = 0.0d0
      l = 0
      do 50 j=i,n
      c = a(j,i)
      if (c.lt.0.0d0) c = -c
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
      if (det.lt.0.0d0) det=-det
c     write(*,*) 'det =',det
8     format(1x,6f12.3)
      end

      subroutine huw6m1(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhuw6m1
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhuw6m1,refc,shpr,3,6,1)
      return
      end

      real*8 function fhuw6m1(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchuw6m/ xa(6),ya(6),za(6),una(6),
     &vna(6),wna(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhuw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhuw6m1=+rx*(+1.-rz)/2. 
      goto 1000
2     fhuw6m1=+ry*(+1.-rz)/2. 
      goto 1000
3     fhuw6m1=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhuw6m1=+rx*(+1.+rz)/2. 
      goto 1000
5     fhuw6m1=+ry*(+1.+rz)/2. 
      goto 1000
6     fhuw6m1=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine huw6m2(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhuw6m2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhuw6m2,refc,shpr,3,6,1)
      return
      end

      real*8 function fhuw6m2(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchuw6m/ xa(6),ya(6),za(6),una(6),
     &vna(6),wna(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhuw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhuw6m2=+rx*(+1.-rz)/2. 
      goto 1000
2     fhuw6m2=+ry*(+1.-rz)/2. 
      goto 1000
3     fhuw6m2=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhuw6m2=+rx*(+1.+rz)/2. 
      goto 1000
5     fhuw6m2=+ry*(+1.+rz)/2. 
      goto 1000
6     fhuw6m2=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine huw6m3(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhuw6m3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhuw6m3,refc,shpr,3,6,1)
      return
      end

      real*8 function fhuw6m3(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchuw6m/ xa(6),ya(6),za(6),una(6),
     &vna(6),wna(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhuw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhuw6m3=+rx*(+1.-rz)/2. 
      goto 1000
2     fhuw6m3=+ry*(+1.-rz)/2. 
      goto 1000
3     fhuw6m3=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhuw6m3=+rx*(+1.+rz)/2. 
      goto 1000
5     fhuw6m3=+ry*(+1.+rz)/2. 
      goto 1000
6     fhuw6m3=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine thuw6m(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,6),coefr(6,9),rc(3,3)
      common /cchuw6m/ x(6),y(6),z(6),un(6),vn(6),wn(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      external fthuw6m
      do 100 n=1,6
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,6
      un(n)=coefr(n,1)
      vn(n)=coefr(n,2)
      wn(n)=coefr(n,3)
      fxx(n)=coefr(n,4)
      fyy(n)=coefr(n,5)
      fzz(n)=coefr(n,6)
      fyz(n)=coefr(n,7)
      fxz(n)=coefr(n,8)
      fxy(n)=coefr(n,9)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(fthuw6m,refc,coor,rc,3,3,1)
      return
      end

      real*8 function fthuw6m(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /cchuw6m/ x(6),y(6),z(6),un(6),vn(6),wn(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      common /vhuw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     fthuw6m=+(+rx*(+1.-rz)/2.)*x(1)+(+ry*(+1.-rz)/2.)*x(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*x(3)+(+rx*(+1.+rz)/2.)*x(4)
     & +(+ry*(+1.+rz)/2.)*x(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*x(6)
      goto 1000
2     fthuw6m=+(+rx*(+1.-rz)/2.)*y(1)+(+ry*(+1.-rz)/2.)*y(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*y(3)+(+rx*(+1.+rz)/2.)*y(4)
     & +(+ry*(+1.+rz)/2.)*y(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*y(6)
      goto 1000
3     fthuw6m=+(+rx*(+1.-rz)/2.)*z(1)+(+ry*(+1.-rz)/2.)*z(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*z(3)+(+rx*(+1.+rz)/2.)*z(4)
     & +(+ry*(+1.+rz)/2.)*z(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*z(6)
      goto 1000
1000  return
      end

      subroutine ehuw6m(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(9),coorr(3,6),coefr(6,9),coefd(9,3)
      external fehuw6m
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fehuw6m,refc,coef,coefd,3,9,2)
      return
      end

      real*8 function fehuw6m(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /cchuw6m/ xa(6),ya(6),za(6),un(6),vn(6),wn(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      common /vhuw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9) n
1     fehuw6m=+(+rx*(+1.-rz)/2.)*un(1)+(+ry*(+1.-rz)/2.)*un(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*un(3)+(+rx*(+1.+rz)/2.)*un(4)
     & +(+ry*(+1.+rz)/2.)*un(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*un(6)
      goto 1000
2     fehuw6m=+(+rx*(+1.-rz)/2.)*vn(1)+(+ry*(+1.-rz)/2.)*vn(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*vn(3)+(+rx*(+1.+rz)/2.)*vn(4)
     & +(+ry*(+1.+rz)/2.)*vn(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*vn(6)
      goto 1000
3     fehuw6m=+(+rx*(+1.-rz)/2.)*wn(1)+(+ry*(+1.-rz)/2.)*wn(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*wn(3)+(+rx*(+1.+rz)/2.)*wn(4)
     & +(+ry*(+1.+rz)/2.)*wn(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*wn(6)
      goto 1000
4     fehuw6m=+(+rx*(+1.-rz)/2.)*fxx(1)+(+ry*(+1.-rz)/2.)*fxx(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxx(3)+(+rx*(+1.+rz)/2.)*fxx(4)
     & +(+ry*(+1.+rz)/2.)*fxx(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxx(6)
      goto 1000
5     fehuw6m=+(+rx*(+1.-rz)/2.)*fyy(1)+(+ry*(+1.-rz)/2.)*fyy(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fyy(3)+(+rx*(+1.+rz)/2.)*fyy(4)
     & +(+ry*(+1.+rz)/2.)*fyy(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fyy(6)
      goto 1000
6     fehuw6m=+(+rx*(+1.-rz)/2.)*fzz(1)+(+ry*(+1.-rz)/2.)*fzz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fzz(3)+(+rx*(+1.+rz)/2.)*fzz(4)
     & +(+ry*(+1.+rz)/2.)*fzz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fzz(6)
      goto 1000
7     fehuw6m=+(+rx*(+1.-rz)/2.)*fyz(1)+(+ry*(+1.-rz)/2.)*fyz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fyz(3)+(+rx*(+1.+rz)/2.)*fyz(4)
     & +(+ry*(+1.+rz)/2.)*fyz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fyz(6)
      goto 1000
8     fehuw6m=+(+rx*(+1.-rz)/2.)*fxz(1)+(+ry*(+1.-rz)/2.)*fxz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxz(3)+(+rx*(+1.+rz)/2.)*fxz(4)
     & +(+ry*(+1.+rz)/2.)*fxz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxz(6)
      goto 1000
9     fehuw6m=+(+rx*(+1.-rz)/2.)*fxy(1)+(+ry*(+1.-rz)/2.)*fxy(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxy(3)+(+rx*(+1.+rz)/2.)*fxy(4)
     & +(+ry*(+1.+rz)/2.)*fxy(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxy(6)
      goto 1000
1000  return
      end


