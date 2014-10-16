       function boundold(r,t,j,dt,i,iblk)
       implicit real*8 (a-h,o-z)
       dimension r(3)
       boundold=0.0d0
       return
       end

       function boundv(r,t,j,dt,i,nblk)
       implicit real*8 (a-h,o-z)
       dimension r(3)
       boundv=0.0d0
       return
       end
      
       function bounda(r,t,j,dt,i,nblk)
       implicit real*8 (a-h,o-z)
       dimension r(3)
       bounda=0.0d0
       return
       end
      
       function boundb(r,t,j,dt,i,nblk)
       implicit real*8 (a-h,o-z)
       dimension r(3)
       boundb=0.d0
       return
       end

