      subroutine suc27g3(coorr,coefr,prmt,estif,emass,edamp,eload,num
     &,fstr6,label6)
      implicit real*8 (a-h,o-z)
      dimension estif(81,81),elump(81),emass(81),
     &eload(81)
      dimension prmt(*),coef(3),coefr(27,3),
     & efuna(81),efunb(81),efunc(81),efund(81),
     & efune(81),efunf(81),coorr(3,27),coor(3)
      common /rsuc27g3/ru(27,32),rv(27,32),rw(27,32),
     & cu(27,4),cv(27,4),cw(27,4)
      common /vsuc27g3/rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      common /dsuc27g3/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(81,3)
      dimension d(6,6),dp(6,6),p(4),dv(6)
      dimension fstr6(6,8),label6(8),f(6)
      dimension stress(7),strain(6),ddv(6)
      external prager
      pe = prmt(1)
      pv = prmt(2)
      fu=prmt(3)
      fv=prmt(4)
      fw=prmt(5)
      rou=prmt(6)
      alpha=prmt(7)
      yita=prmt(8)
      time=prmt(9)
      dt=prmt(10)
      imate=prmt(11)+0.01
      if (num.eq.1) call suc27g3i
      do 10 i=1,nvar
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      ssss=0.0
      tttt=0.0
      do 999 igaus=1,ngaus
      call suc27g3t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call esuc27g3(refc(1,igaus),coef,coorr,coefr,coefd)
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call suc27g31(refc(1,igaus),ru(1,iu),rctr,crtr)
      call suc27g32(refc(1,igaus),rv(1,iv),rctr,crtr)
      call suc27g33(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,27,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rw(1,iw),cw,crtr,1,4,4)
      call shapc(nrefc,ncoor,3,coefd,coefc,crtr,2,9,9)
      un=coef(1)
      vn=coef(2)
      wn=coef(3)
      weigh=det*gaus(igaus)
      do 100 i=1,81
      efuna(i) = 0.0
      efunb(i) = 0.0
      efunc(i) = 0.0
      efund(i) = 0.0
      efune(i) = 0.0
      efunf(i) = 0.0
100   continue
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
      strain(1)=coefc(1,1)
      strain(2)=coefc(2,2)
      strain(3)=coefc(3,3)
      strain(4)=coefc(2,3)+coefc(3,2)
      strain(5)=coefc(1,3)+coefc(3,1)
      strain(6)=coefc(1,2)+coefc(2,1)
      call getdfu(pe,pv,yita,dt,stress,strain,d,f)
      do i=1,6
      ddv(i)=f(i)-stress(i)
      enddo
      if (imate.eq.2) then
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
      do 1000 i=1,6
      do 1000 j=1,6
      d(i,j)=d(i,j)-dp(i,j)
1000   continue
1100  continue
      do 101 i=1,27
      iv=kvord(i,1)
      stif=+cu(i,2) 
      efuna(iv)=efuna(iv)+stif
101   continue
      ssss=ssss+(cu(1,2)+cv(1,3)+cw(1,4))*weigh
      tttt=tttt+weigh
      if (igaus.eq.ngaus) then
      write(*,*) 'num',num,ssss
      endif
      do 102 i=1,27
      iv=kvord(i,2)
      stif=+cv(i,3) 
      efunb(iv)=efunb(iv)+stif
102   continue
      do 103 i=1,27
      iv=kvord(i,3)
      stif=+cw(i,4) 
      efunc(iv)=efunc(iv)+stif
103   continue
      do 104 i=1,27
      iv=kvord(i,2)
      stif=+cv(i,4) 
      efund(iv)=efund(iv)+stif
104   continue
      do 105 i=1,27
      iv=kvord(i,3)
      stif=+cw(i,3) 
      efund(iv)=efund(iv)+stif
105   continue
      do 106 i=1,27
      iv=kvord(i,1)
      stif=+cu(i,4) 
      efune(iv)=efune(iv)+stif
106   continue
      do 107 i=1,27
      iv=kvord(i,3)
      stif=+cw(i,2) 
      efune(iv)=efune(iv)+stif
107   continue
      do 108 i=1,27
      iv=kvord(i,1)
      stif=+cu(i,3) 
      efunf(iv)=efunf(iv)+stif
108   continue
      do 109 i=1,27
      iv=kvord(i,2)
      stif=+cv(i,2) 
      efunf(iv)=efunf(iv)+stif
109   continue
      do 202 iv=1,81
      do 201 jv=1,81
      stif=+efuna(iv)*efuna(jv)*(d(1,1))
     & +efuna(iv)*efunb(jv)*(d(1,2))
     & +efuna(iv)*efunc(jv)*(d(1,3))
     & +efuna(iv)*efund(jv)*(d(1,4))
     & +efuna(iv)*efune(jv)*(d(1,5))
     & +efuna(iv)*efunf(jv)*(d(1,6))
     & +efunb(iv)*efuna(jv)*(d(2,1))
     & +efunb(iv)*efunb(jv)*(d(2,2))
     & +efunb(iv)*efunc(jv)*(d(2,3))
     & +efunb(iv)*efund(jv)*(d(2,4))
     & +efunb(iv)*efune(jv)*(d(2,5))
     & +efunb(iv)*efunf(jv)*(d(2,6))
     & +efunc(iv)*efuna(jv)*(d(3,1))
     & +efunc(iv)*efunb(jv)*(d(3,2))
     & +efunc(iv)*efunc(jv)*(d(3,3))
     & +efunc(iv)*efund(jv)*(d(3,4))
     & +efunc(iv)*efune(jv)*(d(3,5))
     & +efunc(iv)*efunf(jv)*(d(3,6))
     & +efund(iv)*efuna(jv)*(d(4,1))
     & +efund(iv)*efunb(jv)*(d(4,2))
     & +efund(iv)*efunc(jv)*(d(4,3))
     & +efund(iv)*efund(jv)*(d(4,4))
     & +efund(iv)*efune(jv)*(d(4,5))
     & +efund(iv)*efunf(jv)*(d(4,6))
     & +efune(iv)*efuna(jv)*(d(5,1))
     & +efune(iv)*efunb(jv)*(d(5,2))
     & +efune(iv)*efunc(jv)*(d(5,3))
     & +efune(iv)*efund(jv)*(d(5,4))
     & +efune(iv)*efune(jv)*(d(5,5))
     & +efune(iv)*efunf(jv)*(d(5,6))
     & +efunf(iv)*efuna(jv)*(d(6,1))
     & +efunf(iv)*efunb(jv)*(d(6,2))
     & +efunf(iv)*efunc(jv)*(d(6,3))
     & +efunf(iv)*efund(jv)*(d(6,4))
     & +efunf(iv)*efune(jv)*(d(6,5))
     & +efunf(iv)*efunf(jv)*(d(6,6))
      estif(iv,jv)=estif(iv,jv)+stif*weigh
201    continue
202    continue
      do 501 i=1,27
      iv=kvord(i,1)
      stif=+cu(i,1)*fu*0.d0
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,27
      iv=kvord(i,2)
      stif=+cv(i,1)*fv*0.d0
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,27
      iv=kvord(i,3)
      stif=+cw(i,1)*fw*0.d0
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 iv=1,81
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

      subroutine suc27g3i
      implicit real*8 (a-h,o-z)
      common /dsuc27g3/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(81,3)
      ngaus=  8
      ndisp=  3
      nrefc=  3
      ncoor=  3
      nvar = 81
      nnode= 27
      kdord(1)=1
      nvard(1)=27
      kvord(1,1)=1
      kvord(2,1)=25
      kvord(3,1)=4
      kvord(4,1)=34
      kvord(5,1)=61
      kvord(6,1)=28
      kvord(7,1)=10
      kvord(8,1)=31
      kvord(9,1)=7
      kvord(10,1)=37
      kvord(11,1)=64
      kvord(12,1)=40
      kvord(13,1)=73
      kvord(14,1)=79
      kvord(15,1)=67
      kvord(16,1)=46
      kvord(17,1)=70
      kvord(18,1)=43
      kvord(19,1)=13
      kvord(20,1)=49
      kvord(21,1)=16
      kvord(22,1)=58
      kvord(23,1)=76
      kvord(24,1)=52
      kvord(25,1)=22
      kvord(26,1)=55
      kvord(27,1)=19
      kdord(2)=1
      nvard(2)=27
      kvord(1,2)=2
      kvord(2,2)=26
      kvord(3,2)=5
      kvord(4,2)=35
      kvord(5,2)=62
      kvord(6,2)=29
      kvord(7,2)=11
      kvord(8,2)=32
      kvord(9,2)=8
      kvord(10,2)=38
      kvord(11,2)=65
      kvord(12,2)=41
      kvord(13,2)=74
      kvord(14,2)=80
      kvord(15,2)=68
      kvord(16,2)=47
      kvord(17,2)=71
      kvord(18,2)=44
      kvord(19,2)=14
      kvord(20,2)=50
      kvord(21,2)=17
      kvord(22,2)=59
      kvord(23,2)=77
      kvord(24,2)=53
      kvord(25,2)=23
      kvord(26,2)=56
      kvord(27,2)=20
      kdord(3)=1
      nvard(3)=27
      kvord(1,3)=3
      kvord(2,3)=27
      kvord(3,3)=6
      kvord(4,3)=36
      kvord(5,3)=63
      kvord(6,3)=30
      kvord(7,3)=12
      kvord(8,3)=33
      kvord(9,3)=9
      kvord(10,3)=39
      kvord(11,3)=66
      kvord(12,3)=42
      kvord(13,3)=75
      kvord(14,3)=81
      kvord(15,3)=69
      kvord(16,3)=48
      kvord(17,3)=72
      kvord(18,3)=45
      kvord(19,3)=15
      kvord(20,3)=51
      kvord(21,3)=18
      kvord(22,3)=60
      kvord(23,3)=78
      kvord(24,3)=54
      kvord(25,3)=24
      kvord(26,3)=57
      kvord(27,3)=21
      refc(1,1)=5.77350e-001
      refc(2,1)=5.77350e-001
      refc(3,1)=5.77350e-001
      gaus(1)=1.00000e+000
      refc(1,2)=5.77350e-001
      refc(2,2)=5.77350e-001
      refc(3,2)=-5.77350e-001
      gaus(2)=1.00000e+000
      refc(1,3)=5.77350e-001
      refc(2,3)=-5.77350e-001
      refc(3,3)=5.77350e-001
      gaus(3)=1.00000e+000
      refc(1,4)=5.77350e-001
      refc(2,4)=-5.77350e-001
      refc(3,4)=-5.77350e-001
      gaus(4)=1.00000e+000
      refc(1,5)=-5.77350e-001
      refc(2,5)=5.77350e-001
      refc(3,5)=5.77350e-001
      gaus(5)=1.00000e+000
      refc(1,6)=-5.77350e-001
      refc(2,6)=5.77350e-001
      refc(3,6)=-5.77350e-001
      gaus(6)=1.00000e+000
      refc(1,7)=-5.77350e-001
      refc(2,7)=-5.77350e-001
      refc(3,7)=5.77350e-001
      gaus(7)=1.00000e+000
      refc(1,8)=-5.77350e-001
      refc(2,8)=-5.77350e-001
      refc(3,8)=-5.77350e-001
      gaus(8)=1.00000e+000
      end


      subroutine suc27g3t(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tsuc27g3(refc,coor,coorr,coefr,rc)
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

      subroutine suc27g31(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fsuc27g31
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc27g31,refc,shpr,3,27,1)
      return
      end

      real*8 function fsuc27g31(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc27g3/ xa(27),ya(27),za(27),una(27),
     &vna(27),wna(27)
      common /vsuc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     & 16,17,18,19,20,21,22,23,24,25,26,27) n
1     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fsuc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fsuc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fsuc27g31=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fsuc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fsuc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fsuc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fsuc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fsuc27g31=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fsuc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fsuc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fsuc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fsuc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fsuc27g31=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fsuc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fsuc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fsuc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fsuc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine suc27g32(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fsuc27g32
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc27g32,refc,shpr,3,27,1)
      return
      end

      real*8 function fsuc27g32(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc27g3/ xa(27),ya(27),za(27),una(27),
     &vna(27),wna(27)
      common /vsuc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     & 16,17,18,19,20,21,22,23,24,25,26,27) n
1     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fsuc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fsuc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fsuc27g32=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fsuc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fsuc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fsuc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fsuc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fsuc27g32=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fsuc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fsuc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fsuc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fsuc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fsuc27g32=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fsuc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fsuc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fsuc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fsuc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine suc27g33(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fsuc27g33
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc27g33,refc,shpr,3,27,1)
      return
      end

      real*8 function fsuc27g33(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc27g3/ xa(27),ya(27),za(27),una(27),
     &vna(27),wna(27)
      common /vsuc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     & 16,17,18,19,20,21,22,23,24,25,26,27) n
1     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fsuc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fsuc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fsuc27g33=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fsuc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fsuc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fsuc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fsuc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fsuc27g33=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fsuc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fsuc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fsuc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fsuc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fsuc27g33=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fsuc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fsuc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fsuc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fsuc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine tsuc27g3(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,27),coefr(27,3),rc(3,3)
      common /ccsuc27g3/ x(27),y(27),z(27),un(27),vn(27),wn(27)
      external ftsuc27g3
      do 100 n=1,27
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,27
      un(n)=coefr(n,1)
      vn(n)=coefr(n,2)
      wn(n)=coefr(n,3)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftsuc27g3,refc,coor,rc,3,3,1)
      return
      end

      real*8 function ftsuc27g3(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccsuc27g3/ x(27),y(27),z(27),un(27),vn(27),wn(27)
      common /vsuc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftsuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*x(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*x(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*x(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*x(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*x(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*x(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*x(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*x(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*x(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*x(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*x(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*x(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*x(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*x(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*x(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*x(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*x(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*x(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*x(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*x(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*x(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*x(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*x(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*x(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*x(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*x(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*x(7)
      goto 1000
2     ftsuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*y(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*y(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*y(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*y(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*y(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*y(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*y(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*y(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*y(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*y(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*y(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*y(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*y(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*y(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*y(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*y(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*y(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*y(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*y(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*y(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*y(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*y(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*y(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*y(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*y(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*y(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*y(7)
      goto 1000
3     ftsuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*z(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*z(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*z(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*z(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*z(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*z(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*z(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*z(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*z(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*z(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*z(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*z(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*z(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*z(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*z(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*z(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*z(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*z(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*z(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*z(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*z(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*z(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*z(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*z(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*z(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*z(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*z(7)
      goto 1000
1000  return
      end

      subroutine esuc27g3(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(3),coorr(3,27),coefr(27,3),coefd(3,3)
      external fesuc27g3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fesuc27g3,refc,coef,coefd,3,3,2)
      return
      end

      real*8 function fesuc27g3(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccsuc27g3/ xa(27),ya(27),za(27),un(27),vn(27),wn(27)
      common /vsuc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     fesuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*un(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*un(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*un(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*un(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*un(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*un(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*un(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*un(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*un(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*un(13)+(+(+1.-rx**
     & 2)*ry*(+ry-1.)/2.*(+1.-rz**2))*un(22)+(+rx*(+1.+rx)/2.*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*un(14)+(+rx*(+rx-1.)/2.*(+
     & 1.-ry**2)*(+1.-rz**2))*un(25)+(+(+1.-rx**2)*(+1.-ry**
     & 2)*(+1.-rz**2))*un(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+
     & 1.-rz**2))*un(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-
     & rz**2))*un(16)+(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**
     & 2))*un(24)+(+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*un(15)
     & +(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*un(5)
     & +(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*un(17)+(+
     & rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*un(6)+(+
     & rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*un(20)+(+(+
     & 1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2.)*un(26)+(+rx*(+1.+
     & rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*un(18)+(+rx*(+rx-1.)/
     & 2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*un(8)+(+(+1.-rx**2)*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*un(19)+(+rx*(+1.+rx)/2.*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*un(7)
      goto 1000
2     fesuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*vn(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*vn(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*vn(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*vn(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*vn(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*vn(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*vn(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*vn(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*vn(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*vn(13)+(+(+1.-rx**
     & 2)*ry*(+ry-1.)/2.*(+1.-rz**2))*vn(22)+(+rx*(+1.+rx)/2.*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*vn(14)+(+rx*(+rx-1.)/2.*(+
     & 1.-ry**2)*(+1.-rz**2))*vn(25)+(+(+1.-rx**2)*(+1.-ry**
     & 2)*(+1.-rz**2))*vn(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+
     & 1.-rz**2))*vn(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-
     & rz**2))*vn(16)+(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**
     & 2))*vn(24)+(+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*vn(15)
     & +(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*vn(5)
     & +(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*vn(17)+(+
     & rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*vn(6)+(+
     & rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*vn(20)+(+(+
     & 1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2.)*vn(26)+(+rx*(+1.+
     & rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*vn(18)+(+rx*(+rx-1.)/
     & 2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*vn(8)+(+(+1.-rx**2)*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*vn(19)+(+rx*(+1.+rx)/2.*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*vn(7)
      goto 1000
3     fesuc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*wn(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*wn(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*wn(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*wn(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*wn(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*wn(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*wn(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*wn(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*wn(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*wn(13)+(+(+1.-rx**
     & 2)*ry*(+ry-1.)/2.*(+1.-rz**2))*wn(22)+(+rx*(+1.+rx)/2.*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*wn(14)+(+rx*(+rx-1.)/2.*(+
     & 1.-ry**2)*(+1.-rz**2))*wn(25)+(+(+1.-rx**2)*(+1.-ry**
     & 2)*(+1.-rz**2))*wn(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+
     & 1.-rz**2))*wn(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-
     & rz**2))*wn(16)+(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**
     & 2))*wn(24)+(+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*wn(15)
     & +(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*wn(5)
     & +(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*wn(17)+(+
     & rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*wn(6)+(+
     & rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*wn(20)+(+(+
     & 1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2.)*wn(26)+(+rx*(+1.+
     & rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*wn(18)+(+rx*(+rx-1.)/
     & 2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*wn(8)+(+(+1.-rx**2)*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*wn(19)+(+rx*(+1.+rx)/2.*
     & ry*(+1.+ry)/2.*rz*(+1.+rz)/2.)*wn(7)
      goto 1000
1000  return
      end

