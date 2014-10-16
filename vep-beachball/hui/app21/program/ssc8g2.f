      subroutine ssc8g2(coorr,coefr,prmt,estif,emass,eload,num,
     &fstr6,label6,erg)
      implicit real*8 (a-h,o-z)
      dimension estif(48,48),elump(48),emass(48),
     &eload(48)
      dimension prmt(*),coef(4),coefr(8,4),
     & coorr(3,8),coor(3)
      common /rssc8g2/rsa(8,32),rsb(8,32),rsc(8,32),
     & rsd(8,32),rse(8,32),rsf(8,32),
     & csa(8,4),csb(8,4),csc(8,4),csd(8,4),
     & cse(8,4),csf(8,4)
      common /vssc8g2/rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      common /dssc8g2/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(48,6)
      dimension d(6,6),f(6),stress(7),strain(6),zmain(6),fstr6(6,8)
      dimension erg(2,8),ddv(6),p(7),label6(8),dp(6,6),dv(6),de(6)
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
      if (num.eq.1) call ssc8g2i
      do 10 i=1,nvar
      emass(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      avere=0.d0
      tttt=0.d0
      do igaus=1,ngaus
      call ssc8g2t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      call essc8g2(refc(1,igaus),coef,coorr,coefr,coefd)
      isa=(igaus-1)*4+1
      isb=(igaus-1)*4+1
      isc=(igaus-1)*4+1
      isd=(igaus-1)*4+1
      ise=(igaus-1)*4+1
      isf=(igaus-1)*4+1
      if (num.gt.1) goto 20
      call ssc8g21(refc(1,igaus),rsa(1,isa),rctr,crtr)
      call ssc8g22(refc(1,igaus),rsb(1,isb),rctr,crtr)
      call ssc8g23(refc(1,igaus),rsc(1,isc),rctr,crtr)
      call ssc8g24(refc(1,igaus),rsd(1,isd),rctr,crtr)
      call ssc8g25(refc(1,igaus),rse(1,ise),rctr,crtr)
      call ssc8g26(refc(1,igaus),rsf(1,isf),rctr,crtr)
20     continue
      call shapn(nrefc,ncoor,8,rsa(1,isa),csa,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsb(1,isb),csb,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsc(1,isc),csc,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsd(1,isd),csd,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rse(1,ise),cse,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsf(1,isf),csf,crtr,1,4,4)
      call shapc(nrefc,ncoor,4,coefd,coefc,crtr,2,9,9)
      weigh=det*gaus(igaus)
      avere=avere+(coefc(1,1)+coefc(2,2)+coefc(3,3))*weigh
      tttt=tttt+weigh
      enddo
      avere=avere/tttt
      do 999 igaus=1,ngaus
      call ssc8g2t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call essc8g2(refc(1,igaus),coef,coorr,coefr,coefd)
      isa=(igaus-1)*4+1
      isb=(igaus-1)*4+1
      isc=(igaus-1)*4+1
      isd=(igaus-1)*4+1
      ise=(igaus-1)*4+1
      isf=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call ssc8g21(refc(1,igaus),rsa(1,isa),rctr,crtr)
      call ssc8g22(refc(1,igaus),rsb(1,isb),rctr,crtr)
      call ssc8g23(refc(1,igaus),rsc(1,isc),rctr,crtr)
      call ssc8g24(refc(1,igaus),rsd(1,isd),rctr,crtr)
      call ssc8g25(refc(1,igaus),rse(1,ise),rctr,crtr)
      call ssc8g26(refc(1,igaus),rsf(1,isf),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,8,rsa(1,isa),csa,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsb(1,isb),csb,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsc(1,isc),csc,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsd(1,isd),csd,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rse(1,ise),cse,crtr,1,4,4)
      call shapn(nrefc,ncoor,8,rsf(1,isf),csf,crtr,1,4,4)
      call shapc(nrefc,ncoor,4,coefd,coefc,crtr,2,9,9)
      u=coef(1)
      v=coef(2)
      w=coef(3)
      fxyz=coef(4)
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
      call getdfs(pe,pv,yita,dt,stress,strain,d,f)
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
      if (imate.eq.2.or.imate.eq.5) then
      iimate=1
      p(1)=0.d0
      p(2)=30.0d6
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
c      do iq=1,6
c      eng=eng+(stress(iq)+dv(iq)/2)*(strain(iq)-de(iq))
c      enddo
      eng=(strain(1)-de(1))**2+(strain(2)-de(2))**2
     &+(strain(3)-de(3))**2+2*(strain(4)-de(4))**2
     &+2*(strain(5)-de(5))**2+2*(strain(6)-de(6))**2
      eng=dsqrt(eng*2/3)

      if (iimate.eq.1) then
      erg(1,igaus)=erg(1,igaus)*0.d0+eng
      else
      erg(2,igaus)=erg(2,igaus)*0.d0+eng
      endif
c      stress(1)=fsa
c      stress(2)=fsb
c      stress(3)=fsc
c      stress(4)=fsd
c      stress(5)=fse
c      stress(6)=fsf
c      kdgof=6
c      call mstress6(kdgof,stress,zmain)
c      if (iimate.eq.1) then
c      fa=(zmain(6)+zmain(4))/2
c      qie=dabs(zmain(6)-zmain(4))/2
c      comp=qie+fa*0.d0
c      else
c      fa=(zmain(6)+zmain(4))/2
c      qie=dabs(zmain(6)-zmain(4))/2
c      comp=qie+fa*0.4d0
c      endif
c      sp1=zmain(4)
c      sp2=zmain(5)
c      sp3=zmain(6)
c      fse=(sp1**2+sp2**2+sp3**2+(sp1*sp2+sp1*sp3+sp2*sp3)*2*pv)/pe/2
c      fsa=sp1
c      fsb=sp2
c      fsc=sp3
      strmean=(fsa+fsb+fsc)/3
      sj=((fsa-strmean)**2+(fsb-strmean)**2+(fsc-strmean)**2)/2
     &+fsd**2+fse**2+fsf**2
      fsf=erg(2,igaus)
      fsd=dsqrt(sj)
      fse=fxyz
      do 202 i=1,8
      iv=kvord(i,1)
      do 201 j=1,8
      jv=kvord(j,1)
      stif=+csa(i,1)*csa(j,1)*0.0
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
      elump(37)=stif*weigh
      stif=1.
      elump(43)=stif*weigh
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
      elump(38)=stif*weigh
      stif=1.
      elump(44)=stif*weigh
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
      elump(39)=stif*weigh
      stif=1.
      elump(45)=stif*weigh
      stif=1.
      elump(4)=stif*weigh
      stif=1.
      elump(10)=stif*weigh
      stif=1.
      elump(16)=stif*weigh
      stif=1.
      elump(22)=stif*weigh
      stif=1.
      elump(28)=stif*weigh
      stif=1.
      elump(34)=stif*weigh
      stif=1.
      elump(40)=stif*weigh
      stif=1.
      elump(46)=stif*weigh
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
      elump(41)=stif*weigh
      stif=1.
      elump(47)=stif*weigh
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
      stif=1.
      elump(42)=stif*weigh
      stif=1.
      elump(48)=stif*weigh
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
      do 501 i=1,8
      iv=kvord(i,1)
      stif=+csa(i,1)*fsa
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,8
      iv=kvord(i,2)
      stif=+csb(i,1)*fsb
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,8
      iv=kvord(i,3)
      stif=+csc(i,1)*fsc
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 i=1,8
      iv=kvord(i,4)
      stif=+csd(i,1)*fsd
      eload(iv)=eload(iv)+stif*weigh
504   continue
      do 505 i=1,8
      iv=kvord(i,5)
      stif=+cse(i,1)*fse
      eload(iv)=eload(iv)+stif*weigh
505   continue
      do 506 i=1,8
      iv=kvord(i,6)
      stif=+csf(i,1)*fsf
      eload(iv)=eload(iv)+stif*weigh
506   continue
999   continue
      return
      end

      subroutine ssc8g2i
      implicit real*8 (a-h,o-z)
      common /dssc8g2/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(48,6)
      ngaus=  8
      ndisp=  6
      nrefc=  3
      ncoor=  3
      nvar = 48
      nnode=  8
      kdord(1)=1
      nvard(1)=8
      kvord(1,1)=1
      kvord(2,1)=7
      kvord(3,1)=13
      kvord(4,1)=19
      kvord(5,1)=25
      kvord(6,1)=31
      kvord(7,1)=37
      kvord(8,1)=43
      kdord(2)=1
      nvard(2)=8
      kvord(1,2)=2
      kvord(2,2)=8
      kvord(3,2)=14
      kvord(4,2)=20
      kvord(5,2)=26
      kvord(6,2)=32
      kvord(7,2)=38
      kvord(8,2)=44
      kdord(3)=1
      nvard(3)=8
      kvord(1,3)=3
      kvord(2,3)=9
      kvord(3,3)=15
      kvord(4,3)=21
      kvord(5,3)=27
      kvord(6,3)=33
      kvord(7,3)=39
      kvord(8,3)=45
      kdord(4)=1
      nvard(4)=8
      kvord(1,4)=4
      kvord(2,4)=10
      kvord(3,4)=16
      kvord(4,4)=22
      kvord(5,4)=28
      kvord(6,4)=34
      kvord(7,4)=40
      kvord(8,4)=46
      kdord(5)=1
      nvard(5)=8
      kvord(1,5)=5
      kvord(2,5)=11
      kvord(3,5)=17
      kvord(4,5)=23
      kvord(5,5)=29
      kvord(6,5)=35
      kvord(7,5)=41
      kvord(8,5)=47
      kdord(6)=1
      nvard(6)=8
      kvord(1,6)=6
      kvord(2,6)=12
      kvord(3,6)=18
      kvord(4,6)=24
      kvord(5,6)=30
      kvord(6,6)=36
      kvord(7,6)=42
      kvord(8,6)=48
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


      subroutine ssc8g2t(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tssc8g2(refc,coor,coorr,coefr,rc)
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

      subroutine ssc8g21(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g21
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g21,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g21(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g21=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g21=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g21=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g21=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g21=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g21=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g21=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g21=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine ssc8g22(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g22
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g22,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g22(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g22=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g22=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g22=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g22=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g22=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g22=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g22=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g22=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine ssc8g23(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g23
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g23,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g23(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g23=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g23=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g23=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g23=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g23=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g23=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g23=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g23=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine ssc8g24(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g24
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g24,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g24(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g24=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g24=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g24=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g24=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g24=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g24=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g24=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g24=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine ssc8g25(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g25
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g25,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g25(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g25=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g25=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g25=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g25=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g25=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g25=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g25=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g25=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine ssc8g26(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(8,4),rctr(3,3),crtr(3,3)
      external fssc8g26
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc8g26,refc,shpr,3,8,1)
      return
      end

      real*8 function fssc8g26(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc8g2/ xa(8),ya(8),za(8),ua(8),
     &va(8),wa(8),fxyza(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      dimension refc(3)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4,5,6,7,8) n
1     fssc8g26=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
2     fssc8g26=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2. 
      goto 1000
3     fssc8g26=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
4     fssc8g26=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.-rz)/2. 
      goto 1000
5     fssc8g26=+(+1.-rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
6     fssc8g26=+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/2. 
      goto 1000
7     fssc8g26=+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
8     fssc8g26=+(+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2. 
      goto 1000
1000  return
      end

      subroutine tssc8g2(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,8),coefr(8,4),rc(3,3)
      common /ccssc8g2/ x(8),y(8),z(8),u(8),v(8),w(8),
     &fxyz(8)
      external ftssc8g2
      do 100 n=1,8
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,8
      u(n)=coefr(n,1)
      v(n)=coefr(n,2)
      w(n)=coefr(n,3)
      fxyz(n)=coefr(n,4)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftssc8g2,refc,coor,rc,3,3,1)
      return
      end

      real*8 function ftssc8g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccssc8g2/ x(8),y(8),z(8),u(8),v(8),w(8),
     &fxyz(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftssc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*x(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*x(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*x(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*x(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*x(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*x(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*x(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*x(8)
      goto 1000
2     ftssc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*y(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*y(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*y(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*y(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*y(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*y(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*y(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*y(8)
      goto 1000
3     ftssc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*z(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*z(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*z(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*z(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*z(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*z(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*z(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*z(8)
      goto 1000
1000  return
      end

      subroutine essc8g2(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(4),coorr(3,8),coefr(8,4),coefd(4,3)
      external fessc8g2
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fessc8g2,refc,coef,coefd,3,4,2)
      return
      end

      real*8 function fessc8g2(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccssc8g2/ xa(8),ya(8),za(8),u(8),v(8),w(8),
     &fxyz(8)
      common /vssc8g2/ rctr(3,3),crtr(3,3),coefd(4,9),coefc(4,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3,4) n
1     fessc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*u(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*u(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*u(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*u(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*u(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*u(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*u(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*u(8)
      goto 1000
2     fessc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*v(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*v(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*v(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*v(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*v(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*v(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*v(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*v(8)
      goto 1000
3     fessc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*w(1)
     & +(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*w(2)+(+(+1.+
     & rx)/2.*(+1.+ry)/2.*(+1.-rz)/2.)*w(3)+(+(+1.-rx)/2.*(+
     & 1.+ry)/2.*(+1.-rz)/2.)*w(4)+(+(+1.-rx)/2.*(+1.-ry)/2.*
     & (+1.+rz)/2.)*w(5)+(+(+1.+rx)/2.*(+1.-ry)/2.*(+1.+rz)/
     & 2.)*w(6)+(+(+1.+rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*w(7)+(+
     & (+1.-rx)/2.*(+1.+ry)/2.*(+1.+rz)/2.)*w(8)
      goto 1000
4     fessc8g2=+(+(+1.-rx)/2.*(+1.-ry)/2.*(+1.-rz)/2.)*fxyz(1)
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

