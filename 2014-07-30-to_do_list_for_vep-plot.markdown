# to do list for vep-plot

## 1. check the peach-ball plot, from hui, at geo01@/home/jiyang/ss/hui-ss-/results/app3/app31/scompute_5w_test


may have some bug here, for determine the ress, the direction of the principal stress
fortran code
```fortran
      do ii=1,3
c      dip(ii)=acos(v(3,ii))/pi
      resd(ii)=abs(90-abs(acos(vv(3,ii))/pi))
      if (abs(vv(2,ii)) .lt. 1.0e-3) then
      ress(ii)=90.
      else
      ress(ii)=atan(vv(1,ii)/vv(2,ii))/pi
      endif
c      ress(ii)=ress(ii)+alpha/pi
      ress(ii)=ress(ii)

c      write(*,*) ii,ress(ii),resd(ii),acos(0.5d0)
      enddo
```

## 2. plot principal stress
	
	+	test at geo01@/home/jiyang/ss/hui-ss-/results/app3/app31/scompute_5w_test
	
	+	check yang wensheng(Caltech) and David&smith, group stress results, consider provide the product of stress or stress rate for the community use.
	
## 3. check the couple effect of the SAF model

## 4. plot the "pretty" slip rate

## 5. check the roy&royden results, plot sth depth results
	

