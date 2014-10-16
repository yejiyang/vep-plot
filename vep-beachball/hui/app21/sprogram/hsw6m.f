      subroutine hsw6m(coorr,coefr,prmt,estif,emass,eload,num,
     &fstr6,label6,erg)
      implicit real*8 (a-h,o-z)
      dimension estif(36,36),elump(36),emass(36),
     &eload(36)
      dimension prmt(3000),coef(9),coefr(6,9),
     & coorr(3,6),coor(3)
      common /rhsw6m/rsa(6,24),rsb(6,24),rsc(6,24),
     & rsd(6,24),rse(6,24),rsf(6,24),
     & csa(6,4),csb(6,4),csc(6,4),csd(6,4),
     & cse(6,4),csf(6,4)
      common /vhsw6m/rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      common /dhsw6m/ refc(3,6),gaus(6),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(36,6)
      dimension d(6,6),f(6),stress(7),strain(6),zmain(6),fstr6(6,6),
     &erg(2,6),ddv(6),p(4),label6(6),dp(6,6),dv(6),de(6)
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
      if (num.eq.1) call hsw6mi
      do 10 i=1,nvar
      emass(i)=0.0d0
      eload(i)=0.0d0
      do 10 j=1,nvar
      estif(i,j)=0.0d0
10    continue
      do 999 igaus=1,ngaus
      call hsw6mt(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call ehsw6m(refc(1,igaus),coef,coorr,coefr,coefd)
      isa=(igaus-1)*4+1
      isb=(igaus-1)*4+1
      isc=(igaus-1)*4+1
      isd=(igaus-1)*4+1
      ise=(igaus-1)*4+1
      isf=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call hsw6m1(refc(1,igaus),rsa(1,isa),rctr,crtr)
      call hsw6m2(refc(1,igaus),rsb(1,isb),rctr,crtr)
      call hsw6m3(refc(1,igaus),rsc(1,isc),rctr,crtr)
      call hsw6m4(refc(1,igaus),rsd(1,isd),rctr,crtr)
      call hsw6m5(refc(1,igaus),rse(1,ise),rctr,crtr)
      call hsw6m6(refc(1,igaus),rsf(1,isf),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,6,rsa(1,isa),csa,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rsb(1,isb),csb,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rsc(1,isc),csc,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rsd(1,isd),csd,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rse(1,ise),cse,crtr,1,4,4)
      call shapn(nrefc,ncoor,6,rsf(1,isf),csf,crtr,1,4,4)
      call shapc(nrefc,ncoor,9,coefd,coefc,crtr,2,9,9)
      u=coef(1)
      v=coef(2)
      w=coef(3)
      fxx=coef(4)
      fyy=coef(5)
      fzz=coef(6)
      fyz=coef(7)
      fxz=coef(8)
      fxy=coef(9)
      weigh=det*gaus(igaus)
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
      call getdfs(pe,pv,yita,dt,stress,strain,d,f)
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
      fsa=stress(1)+dv(1)
      fsb=stress(2)+dv(2)
      fsc=stress(3)+dv(3)
      fsd=stress(4)+dv(4)
      fse=stress(5)+dv(5)
      fsf=stress(6)+dv(6)
      fstr6(1,igaus)=fsa
      fstr6(2,igaus)=fsb
      fstr6(3,igaus)=fsc
      fstr6(4,igaus)=fsd
      fstr6(5,igaus)=fse
      fstr6(6,igaus)=fsf
      call getde(pe,pv,dv,de)
      eng=0.0d0
      do iq=1,6
      eng=eng+(stress(iq)+dv(iq)/2)*(strain(iq)-de(iq))
      enddo
      if (iimate.eq.1) then
      erg(1,igaus)=erg(1,igaus)+eng
      else
      erg(2,igaus)=erg(2,igaus)+eng
      endif
      stress(1)=fsa
      stress(2)=fsb
      stress(3)=fsc
      stress(4)=fsd
      stress(5)=fse
      stress(6)=fsf
      kdgof=6
      call mstress6(kdgof,stress,zmain)
      if (imate.ge.2.and.imate.le.4) then
      fa=(zmain(6)+zmain(4))/2
      qie=dabs(zmain(6)-zmain(4))/2
      comp=qie+fa*0.d0
      else
      fa=(zmain(6)+zmain(4))/2
      qie=dabs(zmain(6)-zmain(4))/2
      comp=qie+fa*0.d0
      endif
      sp1=zmain(4)
      sp2=zmain(5)
      sp3=zmain(6)
      fse=(sp1**2+sp2**2+sp3**2+(sp1*sp2+sp1*sp3+sp2*sp3)*2*pv)/pe/2
c      fsc=erg(1,igaus)
      fsa=sp1
      fsb=sp2
      fsc=sp3
      fsf=erg(2,igaus)
      fsd=qie
      do 202 i=1,6
      iv=kvord(i,1)
      do 201 j=1,6
      jv=kvord(j,1)
      stif=+csa(i,1)*csa(j,1)*0.0d0
      estif(jv,iv)=estif(jv,iv)+stif*weigh
201    continue
202    continue
      stif=1.
      elump(1)=stif*weigh
      stif=1.
      elump(7)=stif*weigh
      stif=1.
      elump(13)=stif*weigh
      stif=1.
      elump(19)=stif*weigh
      stif=1.
      elump(25)=stif*weigh
      stif=1.
      elump(31)=stif*weigh
      stif=1.
      elump(2)=stif*weigh
      stif=1.
      elump(8)=stif*weigh
      stif=1.
      elump(14)=stif*weigh
      stif=1.
      elump(20)=stif*weigh
      stif=1.
      elump(26)=stif*weigh
      stif=1.
      elump(32)=stif*weigh
      stif=1.
      elump(3)=stif*weigh
      stif=1.
      elump(9)=stif*weigh
      stif=1.
      elump(15)=stif*weigh
      stif=1.
      elump(21)=stif*weigh
      stif=1.
      elump(27)=stif*weigh
      stif=1.
      elump(33)=stif*weigh
      stif=1.
      elump(4)=stif*weigh
      stif=1.
      elump(10)=stif*weigh
      stif=1.
      elump(16)=stif*weigh
      stif=1.
      elump(22)=stif*weigh
      stif=1.
      Elump(28)=stif*weigh
      stif=1.
      elump(34)=stif*weigh
      stif=1.
      elump(5)=stif*weigh
      stif=1.
      elump(11)=stif*weigh
      stif=1.
      elump(17)=stif*weigh
      stif=1.
      elump(23)=stif*weigh
      stif=1.
      elump(29)=stif*weigh
      stif=1.
      elump(35)=stif*weigh
      stif=1.
      elump(6)=stif*weigh
      stif=1.
      elump(12)=stif*weigh
      stif=1.
      elump(18)=stif*weigh
      stif=1.
      elump(24)=stif*weigh
      stif=1.
      elump(30)=stif*weigh
      stif=1.
      elump(36)=stif*weigh
      do 301 i=1,nnode
      if (nvard(1).lt.i) goto 301
      iv = kvord(i,1)
      emass(iv)=emass(iv)+elump(iv)*csa(i,1)
      if (nvard(2).lt.i) goto 301
      iv = kvord(i,2)
      emass(iv)=emass(iv)+elump(iv)*csb(i,1)
      if (nvard(3).lt.i) goto 301
      iv = kvord(i,3)
      emass(iv)=emass(iv)+elump(iv)*csc(i,1)
      if (nvard(4).lt.i) goto 301
      iv = kvord(i,4)
      emass(iv)=emass(iv)+elump(iv)*csd(i,1)
      if (nvard(5).lt.i) goto 301
      iv = kvord(i,5)
      emass(iv)=emass(iv)+elump(iv)*cse(i,1)
      if (nvard(6).lt.i) goto 301
      iv = kvord(i,6)
      emass(iv)=emass(iv)+elump(iv)*csf(i,1)
301   continue
      do 501 i=1,6
      iv=kvord(i,1)
      stif=+csa(i,1)*fsa
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,6
      iv=kvord(i,2)
      stif=+csb(i,1)*fsb
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,6
      iv=kvord(i,3)
      stif=+csc(i,1)*fsc
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 i=1,6
      iv=kvord(i,4)
      stif=+csd(i,1)*fsd
      eload(iv)=eload(iv)+stif*weigh
504   continue
      do 505 i=1,6
      iv=kvord(i,5)
      stif=+cse(i,1)*fse
      eload(iv)=eload(iv)+stif*weigh
505   continue
      do 506 i=1,6
      iv=kvord(i,6)
      stif=+csf(i,1)*fsf
      eload(iv)=eload(iv)+stif*weigh
506   continue
999   continue
      return
      end

      subroutine hsw6mi
      implicit real*8 (a-h,o-z)
      common /dhsw6m/ refc(3,6),gaus(6),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(36,6)
      ngaus=  6
      ndisp=  6
      nrefc=  3
      ncoor=  3
      nvar = 36
      nnode=  6
      kdord(1)=1
      nvard(1)=6
      kvord(1,1)=1
      kvord(2,1)=7
      kvord(3,1)=13
      kvord(4,1)=19
      kvord(5,1)=25
      kvord(6,1)=31
      kdord(2)=1
      nvard(2)=6
      kvord(1,2)=2
      kvord(2,2)=8
      kvord(3,2)=14
      kvord(4,2)=20
      kvord(5,2)=26
      kvord(6,2)=32
      kdord(3)=1
      nvard(3)=6
      kvord(1,3)=3
      kvord(2,3)=9
      kvord(3,3)=15
      kvord(4,3)=21
      kvord(5,3)=27
      kvord(6,3)=33
      kdord(4)=1
      nvard(4)=6
      kvord(1,4)=4
      kvord(2,4)=10
      kvord(3,4)=16
      kvord(4,4)=22
      kvord(5,4)=28
      kvord(6,4)=34
      kdord(5)=1
      nvard(5)=6
      kvord(1,5)=5
      kvord(2,5)=11
      kvord(3,5)=17
      kvord(4,5)=23
      kvord(5,5)=29
      kvord(6,5)=35
      kdord(6)=1
      nvard(6)=6
      kvord(1,6)=6
      kvord(2,6)=12
      kvord(3,6)=18
      kvord(4,6)=24
      kvord(5,6)=30
      kvord(6,6)=36
      refc(1,1)=0.000000e+00
      refc(2,1)=0.000000e+00
      refc(3,1)=-1.000000e+00
      gaus(1)=1.666666e-01
      refc(1,2)=0.000000e+00
      refc(2,2)=1.000000e+00
      refc(3,2)=-1.000000e+00
      gaus(2)=1.666666e-01
      refc(1,3)=1.000000e+00
      refc(2,3)=0.000000e+00
      refc(3,3)=-1.000000e+00
      gaus(3)=1.666666e-01
      refc(1,4)=0.000000e+00
      refc(2,4)=0.000000e+00
      refc(3,4)=1.000000e+00
      gaus(4)=1.666666e-01
      refc(1,5)=0.000000e+00
      refc(2,5)=1.000000e+00
      refc(3,5)=1.000000e+00
      gaus(5)=1.666666e-01
      refc(1,6)=1.000000e+00
      refc(2,6)=0.000000e+00
      refc(3,6)=1.000000e+00
      gaus(6)=1.666666e-01
      end


      subroutine hsw6mt(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call thsw6m(refc,coor,coorr,coefr,rc)
      n=nrefc
      m=n*2
      det = 1.0
      do 10 i=1,n
      do 10 j=1,n
      if (i.le.ncoor) a(i,j) = rc(i,j)
      if (i.gt.ncoor) a(i,j)=1.0
      a(i,n+j)=0.0d0
      if (i.eq.j) a(i,n+i) = 1.0
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

      subroutine hsw6m1(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m1
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m1,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m1(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m1=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m1=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m1=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m1=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m1=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m1=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine hsw6m2(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m2,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m2(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m2=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m2=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m2=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m2=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m2=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m2=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine hsw6m3(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m3,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m3(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m3=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m3=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m3=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m3=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m3=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m3=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine hsw6m4(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m4
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m4,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m4(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m4=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m4=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m4=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m4=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m4=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m4=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine hsw6m5(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m5
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m5,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m5(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m5=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m5=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m5=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m5=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m5=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m5=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine hsw6m6(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(6,4),rctr(3,3),crtr(3,3)
      external fhsw6m6
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fhsw6m6,refc,shpr,3,6,1)
      return
      end

      real*8 function fhsw6m6(refc,n)
      implicit real*8 (a-h,o-z)
      common /cchsw6m/ xa(6),ya(6),za(6),ua(6),
     &va(6),wa(6),fxxa(6),fyya(6),fzza(6),fyza(6),
     &fxza(6),fxya(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6) n
1     fhsw6m6=+rx*(+1.-rz)/2. 
      goto 1000
2     fhsw6m6=+ry*(+1.-rz)/2. 
      goto 1000
3     fhsw6m6=+(+1.-rx-ry)*(+1.-rz)/2. 
      goto 1000
4     fhsw6m6=+rx*(+1.+rz)/2. 
      goto 1000
5     fhsw6m6=+ry*(+1.+rz)/2. 
      goto 1000
6     fhsw6m6=+(+1.-rx-ry)*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine thsw6m(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,6),coefr(6,9),rc(3,3)
      common /cchsw6m/ x(6),y(6),z(6),u(6),v(6),w(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      external fthsw6m
      do 100 n=1,6
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,6
      u(n)=coefr(n,1)
      v(n)=coefr(n,2)
      w(n)=coefr(n,3)
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
      call dcoor(fthsw6m,refc,coor,rc,3,3,1)
      return
      end

      real*8 function fthsw6m(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /cchsw6m/ x(6),y(6),z(6),u(6),v(6),w(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     fthsw6m=+(+rx*(+1.-rz)/2.)*x(1)+(+ry*(+1.-rz)/2.)*x(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*x(3)+(+rx*(+1.+rz)/2.)*x(4)
     & +(+ry*(+1.+rz)/2.)*x(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*x(6)
      goto 1000
2     fthsw6m=+(+rx*(+1.-rz)/2.)*y(1)+(+ry*(+1.-rz)/2.)*y(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*y(3)+(+rx*(+1.+rz)/2.)*y(4)
     & +(+ry*(+1.+rz)/2.)*y(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*y(6)
      goto 1000
3     fthsw6m=+(+rx*(+1.-rz)/2.)*z(1)+(+ry*(+1.-rz)/2.)*z(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*z(3)+(+rx*(+1.+rz)/2.)*z(4)
     & +(+ry*(+1.+rz)/2.)*z(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*z(6)
      goto 1000
1000  return
      end

      subroutine ehsw6m(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(9),coorr(3,6),coefr(6,9),coefd(9,3)
      external fehsw6m
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fehsw6m,refc,coef,coefd,3,9,2)
      return
      end

      real*8 function fehsw6m(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /cchsw6m/ xa(6),ya(6),za(6),u(6),v(6),w(6),
     &fxx(6),fyy(6),fzz(6),fyz(6),fxz(6),fxy(6)
      common /vhsw6m/ rctr(3,3),crtr(3,3),coefd(9,9),coefc(9,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8,9) n
1     fehsw6m=+(+rx*(+1.-rz)/2.)*u(1)+(+ry*(+1.-rz)/2.)*u(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*u(3)+(+rx*(+1.+rz)/2.)*u(4)
     & +(+ry*(+1.+rz)/2.)*u(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*u(6)
      goto 1000
2     fehsw6m=+(+rx*(+1.-rz)/2.)*v(1)+(+ry*(+1.-rz)/2.)*v(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*v(3)+(+rx*(+1.+rz)/2.)*v(4)
     & +(+ry*(+1.+rz)/2.)*v(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*v(6)
      goto 1000
3     fehsw6m=+(+rx*(+1.-rz)/2.)*w(1)+(+ry*(+1.-rz)/2.)*w(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*w(3)+(+rx*(+1.+rz)/2.)*w(4)
     & +(+ry*(+1.+rz)/2.)*w(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*w(6)
      goto 1000
4     fehsw6m=+(+rx*(+1.-rz)/2.)*fxx(1)+(+ry*(+1.-rz)/2.)*fxx(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxx(3)+(+rx*(+1.+rz)/2.)*fxx(4)
     & +(+ry*(+1.+rz)/2.)*fxx(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxx(6)
      goto 1000
5     fehsw6m=+(+rx*(+1.-rz)/2.)*fyy(1)+(+ry*(+1.-rz)/2.)*fyy(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fyy(3)+(+rx*(+1.+rz)/2.)*fyy(4)
     & +(+ry*(+1.+rz)/2.)*fyy(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fyy(6)
      goto 1000
6     fehsw6m=+(+rx*(+1.-rz)/2.)*fzz(1)+(+ry*(+1.-rz)/2.)*fzz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fzz(3)+(+rx*(+1.+rz)/2.)*fzz(4)
     & +(+ry*(+1.+rz)/2.)*fzz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fzz(6)
      goto 1000
7     fehsw6m=+(+rx*(+1.-rz)/2.)*fyz(1)+(+ry*(+1.-rz)/2.)*fyz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fyz(3)+(+rx*(+1.+rz)/2.)*fyz(4)
     & +(+ry*(+1.+rz)/2.)*fyz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fyz(6)
      goto 1000
8     fehsw6m=+(+rx*(+1.-rz)/2.)*fxz(1)+(+ry*(+1.-rz)/2.)*fxz(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxz(3)+(+rx*(+1.+rz)/2.)*fxz(4)
     & +(+ry*(+1.+rz)/2.)*fxz(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxz(6)
      goto 1000
9     fehsw6m=+(+rx*(+1.-rz)/2.)*fxy(1)+(+ry*(+1.-rz)/2.)*fxy(2)
     & +(+(+1.-rx-ry)*(+1.-rz)/2.)*fxy(3)+(+rx*(+1.+rz)/2.)*fxy(4)
     & +(+ry*(+1.+rz)/2.)*fxy(5)+(+(+1.-rx-ry)*(+1.+rz)/2.)*fxy(6)
      goto 1000
1000  return
      end


