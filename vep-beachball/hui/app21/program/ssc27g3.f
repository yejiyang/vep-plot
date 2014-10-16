      subroutine ssc27g3(coorr,coefr,prmt,estif,emass,eload,num,
     &fstr6,label6,erg)
      implicit real*8 (a-h,o-z)
      dimension estif(162,162),elump(162),emass(162),
     &eload(162)
      dimension prmt(*),coef(3),coefr(27,3),
     & coorr(3,27),coor(3)
      common /rssc27g3/rsa(27,32),rsb(27,32),rsc(27,32),
     & rsd(27,32),rse(27,32),rsf(27,32),
     & csa(27,4),csb(27,4),csc(27,4),csd(27,4),
     & cse(27,4),csf(27,4)
      common /vssc27g3/rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      common /dssc27g3/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(162,6)
      dimension d(6,6),f(6),stress(7),strain(6),zmain(6),fstr6(6,8)
      dimension erg(2,8),ddv(6),p(4),label6(8),dp(6,6),dv(6),de(6)
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
      if (num.eq.1) call ssc27g3i
      do 10 i=1,nvar
      emass(i)=0.0
      eload(i)=0.0
      do 10 j=1,nvar
      estif(i,j)=0.0
10    continue
      do 999 igaus=1,ngaus
      call ssc27g3t(nnode,nrefc,ncoor,refc(1,igaus),coor,coorr,
     & rctr,crtr,det,coefr)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1,igaus)
      ry=refc(2,igaus)
      rz=refc(3,igaus)
      call essc27g3(refc(1,igaus),coef,coorr,coefr,coefd)
      isa=(igaus-1)*4+1
      isb=(igaus-1)*4+1
      isc=(igaus-1)*4+1
      isd=(igaus-1)*4+1
      ise=(igaus-1)*4+1
      isf=(igaus-1)*4+1
      if (num.gt.1) goto 2
      call ssc27g31(refc(1,igaus),rsa(1,isa),rctr,crtr)
      call ssc27g32(refc(1,igaus),rsb(1,isb),rctr,crtr)
      call ssc27g33(refc(1,igaus),rsc(1,isc),rctr,crtr)
      call ssc27g34(refc(1,igaus),rsd(1,isd),rctr,crtr)
      call ssc27g35(refc(1,igaus),rse(1,ise),rctr,crtr)
      call ssc27g36(refc(1,igaus),rsf(1,isf),rctr,crtr)
2     continue
      call shapn(nrefc,ncoor,27,rsa(1,isa),csa,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rsb(1,isb),csb,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rsc(1,isc),csc,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rsd(1,isd),csd,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rse(1,ise),cse,crtr,1,4,4)
      call shapn(nrefc,ncoor,27,rsf(1,isf),csf,crtr,1,4,4)
      call shapc(nrefc,ncoor,3,coefd,coefc,crtr,2,9,9)
      u=coef(1)
      v=coef(2)
      w=coef(3)
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
      strain(1)=coefc(1,1)
      strain(2)=coefc(2,2)
      strain(3)=coefc(3,3)
      strain(4)=coefc(2,3)+coefc(3,2)
      strain(5)=coefc(1,3)+coefc(3,1)
      strain(6)=coefc(1,2)+coefc(2,1)
      call getdfs(pe,pv,yita,dt,stress,strain,d,f)
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
      if (iimate.eq.1) then
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
      fsa=sp1
      fsb=sp2
      fsc=sp3
      fsf=erg(2,igaus)
      fsd=qie
      do 202 i=1,27
      iv=kvord(i,1)
      do 201 j=1,27
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
      elump(49)=stif*weigh
      stif=1.
      elump(55)=stif*weigh
      stif=1.
      elump(61)=stif*weigh
      stif=1.
      elump(67)=stif*weigh
      stif=1.
      elump(73)=stif*weigh
      stif=1.
      elump(79)=stif*weigh
      stif=1.
      elump(85)=stif*weigh
      stif=1.
      elump(91)=stif*weigh
      stif=1.
      elump(97)=stif*weigh
      stif=1.
      elump(103)=stif*weigh
      stif=1.
      elump(109)=stif*weigh
      stif=1.
      elump(115)=stif*weigh
      stif=1.
      elump(121)=stif*weigh
      stif=1.
      elump(127)=stif*weigh
      stif=1.
      elump(133)=stif*weigh
      stif=1.
      elump(139)=stif*weigh
      stif=1.
      elump(145)=stif*weigh
      stif=1.
      elump(151)=stif*weigh
      stif=1.
      elump(157)=stif*weigh
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
      elump(50)=stif*weigh
      stif=1.
      elump(56)=stif*weigh
      stif=1.
      elump(62)=stif*weigh
      stif=1.
      elump(68)=stif*weigh
      stif=1.
      elump(74)=stif*weigh
      stif=1.
      elump(80)=stif*weigh
      stif=1.
      elump(86)=stif*weigh
      stif=1.
      elump(92)=stif*weigh
      stif=1.
      elump(98)=stif*weigh
      stif=1.
      elump(104)=stif*weigh
      stif=1.
      elump(110)=stif*weigh
      stif=1.
      elump(116)=stif*weigh
      stif=1.
      elump(122)=stif*weigh
      stif=1.
      elump(128)=stif*weigh
      stif=1.
      elump(134)=stif*weigh
      stif=1.
      elump(140)=stif*weigh
      stif=1.
      elump(146)=stif*weigh
      stif=1.
      elump(152)=stif*weigh
      stif=1.
      elump(158)=stif*weigh
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
      elump(51)=stif*weigh
      stif=1.
      elump(57)=stif*weigh
      stif=1.
      elump(63)=stif*weigh
      stif=1.
      elump(69)=stif*weigh
      stif=1.
      elump(75)=stif*weigh
      stif=1.
      elump(81)=stif*weigh
      stif=1.
      elump(87)=stif*weigh
      stif=1.
      elump(93)=stif*weigh
      stif=1.
      elump(99)=stif*weigh
      stif=1.
      elump(105)=stif*weigh
      stif=1.
      elump(111)=stif*weigh
      stif=1.
      elump(117)=stif*weigh
      stif=1.
      elump(123)=stif*weigh
      stif=1.
      elump(129)=stif*weigh
      stif=1.
      elump(135)=stif*weigh
      stif=1.
      elump(141)=stif*weigh
      stif=1.
      elump(147)=stif*weigh
      stif=1.
      elump(153)=stif*weigh
      stif=1.
      elump(159)=stif*weigh
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
      elump(52)=stif*weigh
      stif=1.
      elump(58)=stif*weigh
      stif=1.
      elump(64)=stif*weigh
      stif=1.
      elump(70)=stif*weigh
      stif=1.
      elump(76)=stif*weigh
      stif=1.
      elump(82)=stif*weigh
      stif=1.
      elump(88)=stif*weigh
      stif=1.
      elump(94)=stif*weigh
      stif=1.
      elump(100)=stif*weigh
      stif=1.
      elump(106)=stif*weigh
      stif=1.
      elump(112)=stif*weigh
      stif=1.
      elump(118)=stif*weigh
      stif=1.
      elump(124)=stif*weigh
      stif=1.
      elump(130)=stif*weigh
      stif=1.
      elump(136)=stif*weigh
      stif=1.
      elump(142)=stif*weigh
      stif=1.
      elump(148)=stif*weigh
      stif=1.
      elump(154)=stif*weigh
      stif=1.
      elump(160)=stif*weigh
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
      elump(53)=stif*weigh
      stif=1.
      elump(59)=stif*weigh
      stif=1.
      elump(65)=stif*weigh
      stif=1.
      elump(71)=stif*weigh
      stif=1.
      elump(77)=stif*weigh
      stif=1.
      elump(83)=stif*weigh
      stif=1.
      elump(89)=stif*weigh
      stif=1.
      elump(95)=stif*weigh
      stif=1.
      elump(101)=stif*weigh
      stif=1.
      elump(107)=stif*weigh
      stif=1.
      elump(113)=stif*weigh
      stif=1.
      elump(119)=stif*weigh
      stif=1.
      elump(125)=stif*weigh
      stif=1.
      elump(131)=stif*weigh
      stif=1.
      elump(137)=stif*weigh
      stif=1.
      elump(143)=stif*weigh
      stif=1.
      elump(149)=stif*weigh
      stif=1.
      elump(155)=stif*weigh
      stif=1.
      elump(161)=stif*weigh
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
      stif=1.
      elump(54)=stif*weigh
      stif=1.
      elump(60)=stif*weigh
      stif=1.
      elump(66)=stif*weigh
      stif=1.
      elump(72)=stif*weigh
      stif=1.
      elump(78)=stif*weigh
      stif=1.
      elump(84)=stif*weigh
      stif=1.
      elump(90)=stif*weigh
      stif=1.
      elump(96)=stif*weigh
      stif=1.
      elump(102)=stif*weigh
      stif=1.
      elump(108)=stif*weigh
      stif=1.
      elump(114)=stif*weigh
      stif=1.
      elump(120)=stif*weigh
      stif=1.
      elump(126)=stif*weigh
      stif=1.
      elump(132)=stif*weigh
      stif=1.
      elump(138)=stif*weigh
      stif=1.
      elump(144)=stif*weigh
      stif=1.
      elump(150)=stif*weigh
      stif=1.
      elump(156)=stif*weigh
      stif=1.
      elump(162)=stif*weigh
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
      do 501 i=1,27
      iv=kvord(i,1)
      stif=+csa(i,1)*fsa
      eload(iv)=eload(iv)+stif*weigh
501   continue
      do 502 i=1,27
      iv=kvord(i,2)
      stif=+csb(i,1)*fsb
      eload(iv)=eload(iv)+stif*weigh
502   continue
      do 503 i=1,27
      iv=kvord(i,3)
      stif=+csc(i,1)*fsc
      eload(iv)=eload(iv)+stif*weigh
503   continue
      do 504 i=1,27
      iv=kvord(i,4)
      stif=+csd(i,1)*fsd
      eload(iv)=eload(iv)+stif*weigh
504   continue
      do 505 i=1,27
      iv=kvord(i,5)
      stif=+cse(i,1)*fse
      eload(iv)=eload(iv)+stif*weigh
505   continue
      do 506 i=1,27
      iv=kvord(i,6)
      stif=+csf(i,1)*fsf
      eload(iv)=eload(iv)+stif*weigh
506   continue
999   continue
      return
      end

      subroutine ssc27g3i
      implicit real*8 (a-h,o-z)
      common /dssc27g3/ refc(3,8),gaus(8),
     & nnode,ngaus,ndisp,nrefc,ncoor,nvar,
     & nvard(6),kdord(6),kvord(162,6)
      ngaus=  8
      ndisp=  6
      nrefc=  3
      ncoor=  3
      nvar =162
      nnode= 27
      kdord(1)=1
      nvard(1)=27
      kvord(1,1)=1
      kvord(2,1)=49
      kvord(3,1)=7
      kvord(4,1)=67
      kvord(5,1)=121
      kvord(6,1)=55
      kvord(7,1)=19
      kvord(8,1)=61
      kvord(9,1)=13
      kvord(10,1)=73
      kvord(11,1)=127
      kvord(12,1)=79
      kvord(13,1)=145
      kvord(14,1)=157
      kvord(15,1)=133
      kvord(16,1)=91
      kvord(17,1)=139
      kvord(18,1)=85
      kvord(19,1)=25
      kvord(20,1)=97
      kvord(21,1)=31
      kvord(22,1)=115
      kvord(23,1)=151
      kvord(24,1)=103
      kvord(25,1)=43
      kvord(26,1)=109
      kvord(27,1)=37
      kdord(2)=1
      nvard(2)=27
      kvord(1,2)=2
      kvord(2,2)=50
      kvord(3,2)=8
      kvord(4,2)=68
      kvord(5,2)=122
      kvord(6,2)=56
      kvord(7,2)=20
      kvord(8,2)=62
      kvord(9,2)=14
      kvord(10,2)=74
      kvord(11,2)=128
      kvord(12,2)=80
      kvord(13,2)=146
      kvord(14,2)=158
      kvord(15,2)=134
      kvord(16,2)=92
      kvord(17,2)=140
      kvord(18,2)=86
      kvord(19,2)=26
      kvord(20,2)=98
      kvord(21,2)=32
      kvord(22,2)=116
      kvord(23,2)=152
      kvord(24,2)=104
      kvord(25,2)=44
      kvord(26,2)=110
      kvord(27,2)=38
      kdord(3)=1
      nvard(3)=27
      kvord(1,3)=3
      kvord(2,3)=51
      kvord(3,3)=9
      kvord(4,3)=69
      kvord(5,3)=123
      kvord(6,3)=57
      kvord(7,3)=21
      kvord(8,3)=63
      kvord(9,3)=15
      kvord(10,3)=75
      kvord(11,3)=129
      kvord(12,3)=81
      kvord(13,3)=147
      kvord(14,3)=159
      kvord(15,3)=135
      kvord(16,3)=93
      kvord(17,3)=141
      kvord(18,3)=87
      kvord(19,3)=27
      kvord(20,3)=99
      kvord(21,3)=33
      kvord(22,3)=117
      kvord(23,3)=153
      kvord(24,3)=105
      kvord(25,3)=45
      kvord(26,3)=111
      kvord(27,3)=39
      kdord(4)=1
      nvard(4)=27
      kvord(1,4)=4
      kvord(2,4)=52
      kvord(3,4)=10
      kvord(4,4)=70
      kvord(5,4)=124
      kvord(6,4)=58
      kvord(7,4)=22
      kvord(8,4)=64
      kvord(9,4)=16
      kvord(10,4)=76
      kvord(11,4)=130
      kvord(12,4)=82
      kvord(13,4)=148
      kvord(14,4)=160
      kvord(15,4)=136
      kvord(16,4)=94
      kvord(17,4)=142
      kvord(18,4)=88
      kvord(19,4)=28
      kvord(20,4)=100
      kvord(21,4)=34
      kvord(22,4)=118
      kvord(23,4)=154
      kvord(24,4)=106
      kvord(25,4)=46
      kvord(26,4)=112
      kvord(27,4)=40
      kdord(5)=1
      nvard(5)=27
      kvord(1,5)=5
      kvord(2,5)=53
      kvord(3,5)=11
      kvord(4,5)=71
      kvord(5,5)=125
      kvord(6,5)=59
      kvord(7,5)=23
      kvord(8,5)=65
      kvord(9,5)=17
      kvord(10,5)=77
      kvord(11,5)=131
      kvord(12,5)=83
      kvord(13,5)=149
      kvord(14,5)=161
      kvord(15,5)=137
      kvord(16,5)=95
      kvord(17,5)=143
      kvord(18,5)=89
      kvord(19,5)=29
      kvord(20,5)=101
      kvord(21,5)=35
      kvord(22,5)=119
      kvord(23,5)=155
      kvord(24,5)=107
      kvord(25,5)=47
      kvord(26,5)=113
      kvord(27,5)=41
      kdord(6)=1
      nvard(6)=27
      kvord(1,6)=6
      kvord(2,6)=54
      kvord(3,6)=12
      kvord(4,6)=72
      kvord(5,6)=126
      kvord(6,6)=60
      kvord(7,6)=24
      kvord(8,6)=66
      kvord(9,6)=18
      kvord(10,6)=78
      kvord(11,6)=132
      kvord(12,6)=84
      kvord(13,6)=150
      kvord(14,6)=162
      kvord(15,6)=138
      kvord(16,6)=96
      kvord(17,6)=144
      kvord(18,6)=90
      kvord(19,6)=30
      kvord(20,6)=102
      kvord(21,6)=36
      kvord(22,6)=120
      kvord(23,6)=156
      kvord(24,6)=108
      kvord(25,6)=48
      kvord(26,6)=114
      kvord(27,6)=42
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


      subroutine ssc27g3t(nnode,nrefc,ncoor,refc,coor,coorr,rc,cr,det,
     &               coefr)
      implicit real*8 (a-h,o-z)
      dimension refc(nrefc),rc(ncoor,nrefc),cr(nrefc,ncoor),a(5,10),
     &          coorr(ncoor,nnode),coor(ncoor),coefr(nnode,*)
      call tssc27g3(refc,coor,coorr,coefr,rc)
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

      subroutine ssc27g31(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g31
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g31,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g31(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g31=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g31=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g31=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g31=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g31=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g31=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g31=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g31=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g31=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g31=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g31=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine ssc27g32(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g32
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g32,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g32(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g32=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g32=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g32=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g32=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g32=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g32=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g32=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g32=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g32=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g32=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g32=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine ssc27g33(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g33
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g33,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g33(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g33=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g33=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g33=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g33=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g33=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g33=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g33=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g33=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g33=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g33=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g33=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine ssc27g34(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g34
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g34,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g34(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g34=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g34=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g34=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g34=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g34=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g34=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g34=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g34=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g34=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g34=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g34=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g34=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g34=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g34=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g34=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g34=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g34=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g34=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g34=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g34=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g34=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g34=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g34=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g34=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g34=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g34=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g34=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine ssc27g35(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g35
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g35,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g35(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g35=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g35=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g35=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g35=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g35=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g35=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g35=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g35=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g35=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g35=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g35=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g35=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g35=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g35=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g35=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g35=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g35=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g35=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g35=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g35=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g35=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g35=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g35=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g35=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g35=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g35=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g35=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine ssc27g36(refc,shpr,rctr,crtr)
      implicit real*8 (a-h,o-z)
      dimension refc(3),shpr(27,4),rctr(3,3),crtr(3,3)
      external fssc27g36
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dshap(fssc27g36,refc,shpr,3,27,1)
      return
      end

      real*8 function fssc27g36(refc,n)
      implicit real*8 (a-h,o-z)
      common /ccssc27g3/ xa(27),ya(27),za(27),ua(27),
     &va(27),wa(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
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
1     fssc27g36=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
2     fssc27g36=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/2. 
      goto 1000
3     fssc27g36=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
4     fssc27g36=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
5     fssc27g36=+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
6     fssc27g36=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2. 
      goto 1000
7     fssc27g36=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
8     fssc27g36=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2. 
      goto 1000
9     fssc27g36=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/
     *2. 
      goto 1000
10     fssc27g36=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
11     fssc27g36=+(+1.-rx**2)*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
12     fssc27g36=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*(+1.-rz**2) 
      goto 1000
13     fssc27g36=+rx*(+rx-1.)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
14     fssc27g36=+(+1.-rx**2)*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
15     fssc27g36=+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**2) 
      goto 1000
16     fssc27g36=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
17     fssc27g36=+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
18     fssc27g36=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2) 
      goto 1000
19     fssc27g36=+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
20     fssc27g36=+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2. 
      goto 1000
21     fssc27g36=+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
22     fssc27g36=+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
23     fssc27g36=+(+1.-rx**2)*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
24     fssc27g36=+rx*(+1.+rx)/2.*(+1.-ry**2)*rz*(+1.+rz)/2. 
      goto 1000
25     fssc27g36=+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
26     fssc27g36=+(+1.-rx**2)*ry*(+1.+ry)/2.*rz*(+1.+rz)/2. 
      goto 1000
27     fssc27g36=+rx*(+1.+rx)/2.*ry*(+1.+ry)/2.*rz*(+1.+rz)/
     *2. 
      goto 1000
1000  return
      end

      subroutine tssc27g3(refc,coor,coorr,coefr,rc)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coor(3),coorr(3,27),coefr(27,3),rc(3,3)
      common /ccssc27g3/ x(27),y(27),z(27),u(27),v(27),w(27)
      external ftssc27g3
      do 100 n=1,27
      x(n)=coorr(1,n)
      y(n)=coorr(2,n)
      z(n)=coorr(3,n)
100   continue
      do 200 n=1,27
      u(n)=coefr(n,1)
      v(n)=coefr(n,2)
      w(n)=coefr(n,3)
200   continue
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoor(ftssc27g3,refc,coor,rc,3,3,1)
      return
      end

      real*8 function ftssc27g3(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccssc27g3/ x(27),y(27),z(27),u(27),v(27),w(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     ftssc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
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
2     ftssc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
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
3     ftssc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
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

      subroutine essc27g3(refc,coef,coorr,coefr,coefd)
      implicit real*8 (a-h,o-z)
      dimension refc(3),coef(3),coorr(3,27),coefr(27,3),coefd(3,3)
      external fessc27g3
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      call dcoef(fessc27g3,refc,coef,coefd,3,3,2)
      return
      end

      real*8 function fessc27g3(refc,n)
      implicit real*8 (a-h,o-z)
      dimension refc(3)
      common /ccssc27g3/ xa(27),ya(27),za(27),u(27),v(27),w(27)
      common /vssc27g3/ rctr(3,3),crtr(3,3),coefd(3,9),coefc(3,9)
      common /coord/ coor(3),coora(27,3)
      x=coor(1)
      y=coor(2)
      z=coor(3)
      rx=refc(1)
      ry=refc(2)
      rz=refc(3)
      goto (1,2,3) n
1     fessc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*u(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*u(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*u(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*u(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*u(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*u(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*u(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*u(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*u(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*u(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*u(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*u(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*u(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*u(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*u(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*u(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*u(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*u(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*u(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*u(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*u(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*u(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*u(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*u(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*u(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*u(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*u(7)
      goto 1000
2     fessc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*v(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*v(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*v(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*v(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*v(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*v(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*v(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*v(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*v(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*v(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*v(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*v(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*v(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*v(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*v(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*v(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*v(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*v(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*v(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*v(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*v(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*v(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*v(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*v(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*v(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*v(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*v(7)
      goto 1000
3     fessc27g3=+(+rx*(+rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+rz-
     & 1.)/2.)*w(1)+(+(+1.-rx**2)*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*w(9)+(+rx*(+1.+rx)/2.*ry*(+ry-1.)/2.*rz*(+rz-1.)/
     & 2.)*w(2)+(+rx*(+rx-1.)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*w(12)
     & +(+(+1.-rx**2)*(+1.-ry**2)*rz*(+rz-1.)/2.)*w(21)+(+rx*
     & (+1.+rx)/2.*(+1.-ry**2)*rz*(+rz-1.)/2.)*w(10)+(+rx*(+
     & rx-1.)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*w(4)+(+(+1.-
     & rx**2)*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*w(11)+(+rx*(+1.+
     & rx)/2.*ry*(+1.+ry)/2.*rz*(+rz-1.)/2.)*w(3)+(+rx*(+rx-
     & 1.)/2.*ry*(+ry-1.)/2.*(+1.-rz**2))*w(13)+(+(+1.-rx**2)*
     & ry*(+ry-1.)/2.*(+1.-rz**2))*w(22)+(+rx*(+1.+rx)/2.*ry*
     & (+ry-1.)/2.*(+1.-rz**2))*w(14)+(+rx*(+rx-1.)/2.*(+1.-
     & ry**2)*(+1.-rz**2))*w(25)+(+(+1.-rx**2)*(+1.-ry**2)*(+
     & 1.-rz**2))*w(27)+(+rx*(+1.+rx)/2.*(+1.-ry**2)*(+1.-rz**
     & 2))*w(23)+(+rx*(+rx-1.)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*w(16)
     & +(+(+1.-rx**2)*ry*(+1.+ry)/2.*(+1.-rz**2))*w(24)+(+rx*
     & (+1.+rx)/2.*ry*(+1.+ry)/2.*(+1.-rz**2))*w(15)+(+rx*(+
     & rx-1.)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*w(5)+(+(+1.-
     & rx**2)*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*w(17)+(+rx*(+1.+
     & rx)/2.*ry*(+ry-1.)/2.*rz*(+1.+rz)/2.)*w(6)+(+rx*(+rx-
     & 1.)/2.*(+1.-ry**2)*rz*(+1.+rz)/2.)*w(20)+(+(+1.-rx**2)*
     & (+1.-ry**2)*rz*(+1.+rz)/2.)*w(26)+(+rx*(+1.+rx)/2.*(+
     & 1.-ry**2)*rz*(+1.+rz)/2.)*w(18)+(+rx*(+rx-1.)/2.*ry*(+
     & 1.+ry)/2.*rz*(+1.+rz)/2.)*w(8)+(+(+1.-rx**2)*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*w(19)+(+rx*(+1.+rx)/2.*ry*(+1.+
     & ry)/2.*rz*(+1.+rz)/2.)*w(7)
      goto 1000
1000  return
      end

