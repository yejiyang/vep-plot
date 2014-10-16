      subroutine suc8g2(coorr,coefr,prmt,estif,emass,edamp,eload,num
     &,fstr6,label6)
      implicit real*8 (a-h,o-z)
      dimension estif(24,24),elump(24),emass(24),
     &eload(24)
      dimension prmt(*),coef(4),coefr(8,4),
     & efuna(24),efunb(24),efunc(24),efund(24),
     & efune(24),efunf(24),coorr(3,8),coor(3)
      common /rsuc8g2/ru(8,32),rv(8,32),rw(8,32),
     & cu(8,4),cv(8,4),cw(8,4)
      common /vsuc8g2/rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      common /dsuc8g2/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(24,3)
      dimension d(6,6),dp(6,6),p(7),dv(6)
      dimension fstr6(6,8),label6(8),f(6)
      dimension stress(7),strain(6),ddv(6)
      dimension averu(8),averv(8),averw(8)
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
      if (num.eq.1) call suc8g2i
      do 10 i=1,nvar
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do i=1,8
      averu(i)=0.d0
      averv(i)=0.d0
      averw(i)=0.d0
      enddo
      avere=0.d0
      tttt=0.d0
      do igaus=1,ngaus
      call suc8g2t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      call esuc8g2(refc(1,igaus),coef,coorr,coefr,coefd)
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 20
      call suc8g21(refc(1,igaus),ru(1,iu),rctr,crtr)
      call suc8g22(refc(1,igaus),rv(1,iv),rctr,crtr)
      call suc8g23(refc(1,igaus),rw(1,iw),rctr,crtr)
20     continue
      call shapn(nrefc,ncoor,8,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rw(1,iw),cw,crtr,1,4,4)
      call shapc(nrefc,ncoor,4,coefd,coefc,crtr,2,9,9)
      weigh=det*gaus(igaus)
      avere=avere+(coefc(1,1)+coefc(2,2)+coefc(3,3))*weigh
      tttt=tttt+weigh
      do i=1,8
      averu(i)=averu(i)+cu(i,2)*weigh
      averv(i)=averv(i)+cv(i,3)*weigh
      averw(i)=averw(i)+cw(i,4)*weigh
      enddo
      enddo
      avere=avere/tttt
      do i=1,8
      averu(i)=averu(i)/tttt
      averv(i)=averv(i)/tttt
      averw(i)=averw(i)/tttt
      enddo
      do 999 igaus=1,ngaus
      call suc8g2t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call esuc8g2(refc(1,igaus),coef,coorr,coefr,coefd)
      iu=(igaus-1)*4+1
      iv=(igaus-1)*4+1
      iw=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call suc8g21(refc(1,igaus),ru(1,iu),rctr,crtr)
      call suc8g22(refc(1,igaus),rv(1,iv),rctr,crtr)
      call suc8g23(refc(1,igaus),rw(1,iw),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,8,ru(1,iu),cu,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rv(1,iv),cv,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rw(1,iw),cw,crtr,1,4,4)
      call shapc(nrefc,ncoor,4,coefd,coefc,crtr,2,9,9)
      un=coef(1)
      vn=coef(2)
      wn=coef(3)
      fxyz=coef(4)
      weigh=det*gaus(igaus)
      do 100 i=1,24
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
      stress(7)=fxyz
      strain(1)=coefc(1,1)-coefc(1,1)/3-coefc(2,2)
     & /3-coefc(3,3)/3+avere/3
      strain(2)=coefc(2,2)-coefc(1,1)/3-coefc(2,2)
     & /3-coefc(3,3)/3+avere/3
      strain(3)=coefc(3,3)-coefc(1,1)/3-coefc(2,2)
     & /3-coefc(3,3)/3+avere/3
      strain(4)=coefc(2,3)+coefc(3,2)
      strain(5)=coefc(1,3)+coefc(3,1)
      strain(6)=coefc(1,2)+coefc(2,1)
      call getdfu(pe,pv,yita,dt,stress,strain,d,f)
      do i=1,6
      ddv(i)=f(i)-stress(i)
      enddo
c.....amend by wanghui
      iimate=2
      p(1)=0.4d0
      p(2)=50.0d6
      p(3)=prmt(2)
      p(4)=0.0d0
      if (imate.eq.1.or.imate.eq.4) then
      iimate=1
      p(1)=0.d0
      p(2)=10.0d6
      p(3)=prmt(2)
      p(4)=0.0d0
      endif

c      if (imate.eq.1.or.imate.eq.7) then
c      iimate=1
c      p(1)=0.d0
c      p(2)=10.0d6
c      p(3)=prmt(2)
c      p(4)=0.0d0
cc      if (time.ge.30000*3.1536d7) then
cc      p(2)=10.0d6-1.d2*(time/3.1536d7-30000)
cc      endif
cc      if (time.ge.40000*3.1536d7) then
cc      p(2)=9.d6
cc      endif
cc      if (imate.eq.2) p(2)=30.0d6
cc      if (label6(igaus).eq.1) p(2)=0.5d6
c      else
c      iimate=2
c      p(1)=0.4d0
c      p(2)=50.0d6
c      p(3)=prmt(2)
c      p(4)=0.0d0
c      endif
      call getstr(p,stress,strain,dv,ddv,d,dp,1,prag,prager)
      if (prag.lt.-p(2)*1.d-5) goto 1100
c      if (label6(igaus).eq.0) then
c      label6(igaus)=2
c      if (iimate.eq.2) label6(igaus)=1
c      endif
      do 1000 i=1,6
      do 1000 j=1,6
      d(i,j)=d(i,j)-dp(i,j)
1000   continue
1100  continue
      do 101 i=1,8
      iv=kvord(i,1)
      stif=+cu(i,2)
     & -cu(i,2)/3+averu(i)/3
      efuna(iv)=efuna(iv)+stif
101   continue
      do 102 i=1,8
      iv=kvord(i,2)
      stif=-cv(i,3)/3+averv(i)/3
      efuna(iv)=efuna(iv)+stif
102   continue
      do 103 i=1,8
      iv=kvord(i,3)
      stif=-cw(i,4)/3+averw(i)/3
      efuna(iv)=efuna(iv)+stif
103   continue
      do 104 i=1,8
      iv=kvord(i,1)
      stif=-cu(i,2)/3+averu(i)/3
      efunb(iv)=efunb(iv)+stif
104   continue
      do 105 i=1,8
      iv=kvord(i,2)
      stif=+cv(i,3)
     & -cv(i,3)/3+averv(i)/3
      efunb(iv)=efunb(iv)+stif
105   continue
      do 106 i=1,8
      iv=kvord(i,3)
      stif=-cw(i,4)/3+averw(i)/3
      efunb(iv)=efunb(iv)+stif
106   continue
      do 107 i=1,8
      iv=kvord(i,1)
      stif=-cu(i,2)/3+averu(i)/3
      efunc(iv)=efunc(iv)+stif
107   continue
      do 108 i=1,8
      iv=kvord(i,2)
      stif=-cv(i,3)/3+averv(i)/3
      efunc(iv)=efunc(iv)+stif
108   continue
      do 109 i=1,8
      iv=kvord(i,3)
      stif=+cw(i,4)
     & -cw(i,4)/3+averw(i)/3
      efunc(iv)=efunc(iv)+stif
109   continue
      do 110 i=1,8
      iv=kvord(i,2)
      stif=+cv(i,4) 
      efund(iv)=efund(iv)+stif
110   continue
      do 111 i=1,8
      iv=kvord(i,3)
      stif=+cw(i,3) 
      efund(iv)=efund(iv)+stif
111   continue
      do 112 i=1,8
      iv=kvord(i,1)
      stif=+cu(i,4) 
      efune(iv)=efune(iv)+stif
112   continue
      do 113 i=1,8
      iv=kvord(i,3)
      stif=+cw(i,2) 
      efune(iv)=efune(iv)+stif
113   continue
      do 114 i=1,8
      iv=kvord(i,1)
      stif=+cu(i,3) 
      efunf(iv)=efunf(iv)+stif
114   continue
      do 115 i=1,8
      iv=kvord(i,2)
      stif=+cv(i,2) 
      efunf(iv)=efunf(iv)+stif
115   continue
      do 202 iv=1,24
      do 201 jv=1,24
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
      do 501 i=1,8
      iv=kvord(i,1)
      stif=+cu(i,1)*fu*0.d0
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,8
      iv=kvord(i,2)
      stif=+cv(i,1)*fv*0.d0
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,8
      iv=kvord(i,3)
      stif=+cw(i,1)*fw*0.d0
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 iv=1,24
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

      subroutine suc8g2i
      implicit real*8 (a-h,o-z)
      common /dsuc8g2/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(3),kdord(3),kvord(24,3)
      ngaus=  8
      ndisp=  3
      nrefc=  3
      ncoor=  3
      nvar = 24
      nnode=  8
      kdord(1)=1
      nvard(1)=8
      kvord(1,1)=1
      kvord(2,1)=4
      kvord(3,1)=7
      kvord(4,1)=10
      kvord(5,1)=13
      kvord(6,1)=16
      kvord(7,1)=19
      kvord(8,1)=22
      kdord(2)=1
      nvard(2)=8
      kvord(1,2)=2
      kvord(2,2)=5
      kvord(3,2)=8
      kvord(4,2)=11
      kvord(5,2)=14
      kvord(6,2)=17
      kvord(7,2)=20
      kvord(8,2)=23
      kdord(3)=1
      nvard(3)=8
      kvord(1,3)=3
      kvord(2,3)=6
      kvord(3,3)=9
      kvord(4,3)=12
      kvord(5,3)=15
      kvord(6,3)=18
      kvord(7,3)=21
      kvord(8,3)=24
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


      subroutine suc8g2t(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tsuc8g2(refc,coor,coorr,coefr,rc)
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

      subroutine suc8g21(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fsuc8g21
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc8g21,refc,shpr,3,8,1)
      return
      end

      real*8 function fsuc8g21(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc8g2/ xa(8),ya(8),za(8),una(8),
     &vna(8),wna(8),fxyza(8)
      common /vsuc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fsuc8g21=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fsuc8g21=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fsuc8g21=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fsuc8g21=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fsuc8g21=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fsuc8g21=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fsuc8g21=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fsuc8g21=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine suc8g22(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fsuc8g22
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc8g22,refc,shpr,3,8,1)
      return
      end

      real*8 function fsuc8g22(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc8g2/ xa(8),ya(8),za(8),una(8),
     &vna(8),wna(8),fxyza(8)
      common /vsuc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fsuc8g22=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fsuc8g22=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fsuc8g22=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fsuc8g22=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fsuc8g22=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fsuc8g22=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fsuc8g22=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fsuc8g22=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine suc8g23(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fsuc8g23
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fsuc8g23,refc,shpr,3,8,1)
      return
      end

      real*8 function fsuc8g23(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccsuc8g2/ xa(8),ya(8),za(8),una(8),
     &vna(8),wna(8),fxyza(8)
      common /vsuc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fsuc8g23=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fsuc8g23=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fsuc8g23=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fsuc8g23=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fsuc8g23=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fsuc8g23=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fsuc8g23=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fsuc8g23=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine tsuc8g2(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,8),coefr(8,4),rc(3,3)
      common /ccsuc8g2/ x(8),y(8),z(8),un(8),vn(8),wn(8),
     &fxyz(8)
      external ftsuc8g2
      do 100 n=1,8
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,8
      un(n)=coefr(n,1)
      vn(n)=coefr(n,2)
      wn(n)=coefr(n,3)
      fxyz(n)=coefr(n,4)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftsuc8g2,refc,coor,rc,3,3,1)
      return
      end

      real*8 function ftsuc8g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccsuc8g2/ x(8),y(8),z(8),un(8),vn(8),wn(8),
     &fxyz(8)
      common /vsuc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftsuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*x(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*x(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*x(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*x(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*x(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*x(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*x(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*x(8)
      goto 1000
2     ftsuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*y(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*y(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*y(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*y(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*y(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*y(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*y(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*y(8)
      goto 1000
3     ftsuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*z(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*z(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*z(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*z(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*z(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*z(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*z(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*z(8)
      goto 1000
1000  return
      end

      subroutine esuc8g2(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(4),coorr(3,8),coefr(8,4),coefd(4,3)
      external fesuc8g2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fesuc8g2,refc,coef,coefd,3,4,2)
      return
      end

      real*8 function fesuc8g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccsuc8g2/ xa(8),ya(8),za(8),un(8),vn(8),wn(8),
     &fxyz(8)
      common /vsuc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fesuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*un(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*un(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*un(3)+(+(+1.-rx)/2.*
     & (+1.+ry)/2.*(+1.-rz)/2.)*un(4)+(+(+1.-rx)/2.*(+1.-ry)/
     & 2.*(+1.+rz)/2.)*un(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+
     & rz)/2.)*un(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*un(7)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*un(8)
      goto 1000
2     fesuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*vn(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*vn(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*vn(3)+(+(+1.-rx)/2.*
     & (+1.+ry)/2.*(+1.-rz)/2.)*vn(4)+(+(+1.-rx)/2.*(+1.-ry)/
     & 2.*(+1.+rz)/2.)*vn(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+
     & rz)/2.)*vn(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*vn(7)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*vn(8)
      goto 1000
3     fesuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*wn(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*wn(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*wn(3)+(+(+1.-rx)/2.*
     & (+1.+ry)/2.*(+1.-rz)/2.)*wn(4)+(+(+1.-rx)/2.*(+1.-ry)/
     & 2.*(+1.+rz)/2.)*wn(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+
     & rz)/2.)*wn(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*wn(7)
     & +(+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*wn(8)
      goto 1000
4     fesuc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*fxyz(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*fxyz(2)+(+(+
     & 1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*fxyz(3)+(+(+1.-rx)/
     & 2.*(+1.+ry)/2.*(+1.-rz)/2.)*fxyz(4)+(+(+1.-rx)/2.*(+
     & 1.-ry)/2.*(+1.+rz)/2.)*fxyz(5)+(+(+1.+rx)/2.*(+1.-ry)/
     & 2.*(+1.+rz)/2.)*fxyz(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+
     & 1.+rz)/2.)*fxyz(7)+(+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/
     & 2.)*fxyz(8)
      goto 1000
1000  return
      end

