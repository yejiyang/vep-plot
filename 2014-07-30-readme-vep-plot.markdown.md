# README for vep-plot

## 1. variables record
	
###	a.	for /scompute, p1234.f, 
		
		? 	在你拷给我的 saltonsea app21 -scompute -  p1234.f 的后处理程序里，  u11 - ... u16，等于 应力的六个分量 ？
		
		u11 = u11
		u12 = u22
		u13 = u33
		u14 = u23 
		u15 = u13
		u16 = u12  
		
		 * or
		| u11 | u16 | u15 |
		| u16 | u12 | u14 |
		| u15 | u14 | u13 |
		
***
		 
		? 	dd(1),ress(1),resd(1),dd(2),ress(2),resd(2), dd(3),ress(3),resd(3),ttt_s, 其他的几个具体代表啥？ 这个导出的结果用gmt画的话是用psmeca 的  s -c, a, m, q的哪一个啊？ 


		*		
```fortran
write(12,1000) 1,lon,lat,u01,u02,u03,u11,u12,u13,u14,u15,u16,dd(1),ress(1),resd(1),dd(2),ress(2),resd(2),dd(3),ress(3),resd(3),ttt_s
```

		*	dd是偏应力的三个主轴分量（特征值），ress是应力主轴的方位，resd为倾角，

		*	awk '{if(($1==1)&&(int($2*10/5)==$2*10/5)&&(int($3*10/5)==$3*10/5)) print $2,$3,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20,$2,$3}' $input1 |psmeca -R $projection -Sx0.30 -G0/0/0 -K -O -o-1 -W2 -L -C >>$output
