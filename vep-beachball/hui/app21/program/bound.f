       function bound(r,t,j,dt,it)
      implicit real*8 (a-h,o-z)
      dimension r(3)
      bound= 0.d0
c      pi=3.1415926/180
      if (j.eq.2.and.r(2).lt.100) then
      bound=dt*0.050d0/3.1536d7 
      endif
      return
      end
 
       function bound1(r,t,j)
      implicit real*8 (a-h,o-z)
      dimension r(3)
      bound1=0.d0
      return
      end
 
       function bound2(r,t,j)
      implicit real*8 (a-h,o-z)
      dimension r(3)
      bound2=0.d0
      return
      end
