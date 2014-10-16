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