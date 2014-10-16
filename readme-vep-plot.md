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
		
		* so for 2D case, its, u11, u12, u16, = $7, $8, $12
		| u11 | u16 | 
		| u16 | u12 | 
		
***
		 
		? 	dd(1),ress(1),resd(1),dd(2),ress(2),resd(2), dd(3),ress(3),resd(3),ttt_s, 其他的几个具体代表啥？ 这个导出的结果用gmt画的话是用psmeca 的  s -c, a, m, q的哪一个啊？ 

+	(fortran code)
```fortran
write(12,1000) 1,lon,lat,u01,u02,u03,u11,u12,u13,u14,u15,u16,dd(1),ress(1),resd(1),dd(2),ress(2),resd(2),dd(3),ress(3),resd(3),ttt_s
```

		*	dd是偏应力的三个主轴分量（特征值），ress是应力主轴的方位，resd为倾角，

+	(bash code)
```bash
awk '{if(($1==1)&&(int($2*10/5)==$2*10/5)&&(int($3*10/5)==$3*10/5)) print $2,$3,1,$14-35,$15,0,$17-35,$18,-1,$20-35,$21,$22+20,$2,$3}' $input1 |psmeca -R $projection -Sx0.30 -G0/0/0 -K -O -o-1 -W2 -L -C >>$output
```

```fortran

c-----------------------------------------------------------------
c.....Principal stress
c-----------------------------------------------------------------
      ttt_s=(u11+u12+u13)/3.d0
      uu11=u11-ttt_s
      uu12=u12-ttt_s
      uu13=u13-ttt_s
      
      Nn=3
      do ii=1,Nn
      do jj=1,Nn
      aa(ii,jj)=0.
      enddo
      enddo
      aa(1,1)=uu11
      aa(1,2)=u16
      aa(1,3)=u15
      aa(2,1)=aa(1,2)
      aa(2,2)=uu12
      aa(2,3)=u14
      aa(3,1)=aa(1,3)
      aa(3,2)=aa(2,3)
      aa(3,3)=uu13

      call jacobi(aa,NP,NP,dd,vv,nrot)
      call eigsrt(dd,vv,NP,NP)

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

      ttt_s=sqrt(u11*u11+u12*u12+2*u16*u16)
      ttt_s=log(abs(ttt_s))/log(10.d0)

      write(12,1000) 1,lon,lat,u01,u02,u03,u11,u12,u13,u14,u15,u16,
     +  dd(1),ress(1),resd(1),dd(2),ress(2),resd(2),
     +  dd(3),ress(3),resd(3),ttt_s
       
```

## 2. From Yang

```bash

psmeca cmt1 -R236/250/33/43 -Jm0.40i  -C -O -K -Sd0.15i/1 -G255/0/0 >> cmt1.ps

```

```fortran

			open(1,file='cmt',status='old')
            do i=1,np
            read(1,*)(s(i,j),j=1,6),idx(i),av(i),ta(i)
            enddo
            close (1)

            open(2,file='cmt1')
            do i=1,kp
            k=npk(i)
            kk=k+iiz*n0
    c        write(2,9)x(k),y(k),-z(k),(s(kk,j),j=1,6),idx(kk),x(k),y(k)
            write(2,9)x(k),y(k),5.,(s(kk,j),j=1,6),idx(kk),x(k),y(k)
            enddo
            close (2)

```

```fortran

	OPEN(1,FILE='tmp')
    OPEN(2,FILE='cmt')
    OPEN(3,FILE='stress')
    ace=0.12
	!    DO  i=1,np
    DO  i=1,1      ! for test
        DO j=1,6
            t(j)=s(i,j)*(1.e-6)
            print*, t(j)
        END DO

        CALL spstr(t,ps,dir)

        av=(ps(1)+ps(2)+ps(3))/3.
        !print*, av !
        tal=(ps(1)-ps(3))/2
        ta=tal
        !print*, tal !
        IF(tal < 1.e-6) tal=1.e-6
        !print*, tal !
        DO j=1,3
            !print*,t(j)
            t(j)=(t(j)-av)/tal
            !print*,t(j)
        END DO

        DO j=4,6
            !print*,t(j)
            t(j)=t(j)/tal
            !print*,t(j)
        END DO

        tx=tal*ace+20.
        !print*,tx
        IF(tx > 28)tx=28.

        idx=INT(tx)
        !print*,LOG(10.)
        !print*, (tx-idx)*LOG(10.)
        frc=EXP((tx-idx)*LOG(10.))
        !print*, frc

        IF(frc > 10.)THEN
            WRITE(*,*)'i,tal,ta,tx,idx,frc'
            WRITE(*,*)i,tal,ta,tx,idx,frc
            frc=10.
        END IF

        DO j=1,6
            !print*,t(j)
            t(j)=t(j)*frc
        END DO

        WRITE(1,9)ps(1),ps(2),ps(3),av,ta
        WRITE(2,9)t(3),t(2),t(1),-t(4),t(5),-t(6),idx,av,ta
        WRITE(3,8) (s(i,j),j=1,6)
    END DO
```

## 3. Sun

```matlab

    for i=1:11*11
        sxx=sxx2(i);sxy=sxy2(i);sxz=sxz2(i);
        syx=sxy2(i);syy=syy2(i);syz=syz2(i);
        szx=sxz2(i);szy=syz2(i);szz=szz2(i);
        s=[sxx sxy sxz;syx syy syz;szx szy szz];
        [v,d]=eig(s);
        dd=eig(s);
        Ta(i)=max(dd);
        Pa(i)=min(dd);
        for j=1:3
            [thet(j),phi(j),r(j)]=cart2sph(v(1,j),v(2,j),v(3,j));
            if dd(j)==Pa(i)
                %             Op=v(j,:);
                if thet(j)<=pi/2
                    Op=[pi/2-thet(j),abs(phi(j))];
                else
                    Op=[2.5*pi-thet(j),abs(phi(j))];
                end
            elseif dd(j)==Ta(i)
                %             Ot=v(j,:);
                if thet(j)<=pi/2
                    Ot=[pi/2-thet(j),abs(phi(j))];
                else
                    Ot=[2.5*pi-thet(j),abs(phi(j))];
                end
            else
                Ba(i)=dd(j);
                %             Ob=v(j,:);
                if thet(j)<=pi/2
                    Ob=[pi/2-thet(j),abs(phi(j))];
                else
                    Ob=[2.5*pi-thet(j),abs(phi(j))];
                end
            end

        end
         TPB=[abs(Ta(i)) abs(Pa(i)) abs(Ba(i))];
        absTPB(i)=max(TPB);
        for k=1:15
            %         las=(Pa(i).^2+Ta(i).^2+Ba(i).^2).^0.5./10.^k ;
            las=absTPB(i)./10.^k ;
            if abs(las)<10.0&abs(las)>=1.0
                pp(i)=k;
            end
        end
        ep(i)=Pa(i)./10.^pp(i);
        et(i)=Ta(i)./10.^pp(i);
        eb(i)=Ba(i)./10.^pp(i);
    fprintf(fp,'%f %f %f %f %.2f %.2f %f %.2f %.2f %f %.2f %.2f %d\n',x1(i),y1(i),abs(zz(i)),et(i),Ot(1)*180/pi,Ot(2)*180/pi,eb(i),Ob(1)*180/pi,Ob(2)*180/pi,ep(i),Op(1)*180/pi,Op(2)*180/pi,pp(i));
    end
    
```

```bash

    psbasemap -R119.0/123.0/20.5/25.5 -JM3i -B1/1WSeN -P -K > $psfile

    #psmeca -R -J -O earthmech4-m.txt -Sy0.13i -W -Gblack -Ewhite -C -K >> $psfile
    psmeca -R -J -O earthmech4-m.txt -Sy0.13i -W -Gblack -Ewhite -C  >> $psfile

    #psxy -R -J -O -K -G240 -L -Wthicker box >> $psfile
    #psmeca -R -J -O epicenter.txt -Sy0.13i -W -C -K >> $psfile
    #echo 122.7     20.6    10      0       4       2       10MPa | pstext -R -J -O -K >> $psfile
    #echo 122.5     20.9    12      0       4       2       10km  | pstext -R -J -O -K >> $psfile

    gnome-open $psfile

```