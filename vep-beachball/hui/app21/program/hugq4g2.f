      subroutine hugq4g2(r,coefr,prmt,egs,egm,egd,egl,num)
      implicit real*8 (a-h,o-z)
      dimension egs(12,12),egm(12),egd(12),egl(12),
     & els(12,12),elm(12),eld(12),ell(12),nc(3),
     & r(3,4),coefr(4,3),y(2,4),z(3,3),prmt(3000),t(12,12)
      ngvar = 12
      nlvar = 12
      ngdim = 3
      nldim = 2
      nnode = 4
      nc(1) = 1
      nc(2) = 2
      nc(3) = 3
      ncoef = 0
      do 1 i=1,ngdim
      do 1 j=1,nnode
1     coefr(j,i+ncoef) = r(i,j)
      call smit(ngdim,nldim,nnode,r,y,z,nc)
      call hulq4g2(y,coefr,prmt,els,elm,eld,ell,num)
      do 10 i=1,nlvar
      do 10 j=1,ngvar
10    t(i,j) = 0.0
      t(1,1)=z(1,1)
      t(1,2)=z(1,2)
      t(1,3)=z(1,3)
      t(4,4)=z(1,1)
      t(4,5)=z(1,2)
      t(4,6)=z(1,3)
      t(7,7)=z(1,1)
      t(7,8)=z(1,2)
      t(7,9)=z(1,3)
      t(10,10)=z(1,1)
      t(10,11)=z(1,2)
      t(10,12)=z(1,3)
      t(2,1)=z(2,1)
      t(2,2)=z(2,2)
      t(2,3)=z(2,3)
      t(5,4)=z(2,1)
      t(5,5)=z(2,2)
      t(5,6)=z(2,3)
      t(8,7)=z(2,1)
      t(8,8)=z(2,2)
      t(8,9)=z(2,3)
      t(11,10)=z(2,1)
      t(11,11)=z(2,2)
      t(11,12)=z(2,3)
      t(3,1)=z(3,1)
      t(3,2)=z(3,2)
      t(3,3)=z(3,3)
      t(6,4)=z(3,1)
      t(6,5)=z(3,2)
      t(6,6)=z(3,3)
      t(9,7)=z(3,1)
      t(9,8)=z(3,2)
      t(9,9)=z(3,3)
      t(12,10)=z(3,1)
      t(12,11)=z(3,2)
      t(12,12)=z(3,3)
      call tkt(ngvar,nlvar,t,els,egs)
      call tmt(ngvar,nlvar,t,elm,egm)
      call ntmt(ngvar,nlvar,t,eld,egd)
      call tl(ngvar,nlvar,t,ell,egl)
      return
      end
