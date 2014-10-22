clear;
clc;
close all;

%% %%%%%%%%%%%%%%% transfer the degree to meter%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% set a local center, using the center of the region as local(0,0) 
%  plot_boundary of saf716
figure(1) % figure initial degree

longlat=[-121.6,-113.2,-113.2,-121.6,-121.6;
          31.4,31.4,38.1,38.1,31.4]';
line(longlat(:,1),longlat(:,2));
xy=longlat; % initial xy

    x0=-117.774; y0=34.736; % local(0,0)
for i=1:length(longlat)
    x=longlat(i,1); y=longlat(i,2);
    xy(i,1)=111263*abs(cos((pi/180)*((y+y0)/2)))*(x-x0);
    xy(i,2)=111263*(y-y0);
end

figure(2)% show the meter
line(xy(:,1),xy(:,2)); hold on


%% %%%%%%% rotate 45 degrees clockwise, 320 anticlockwise
%  
d=-45*pi/180;
for i=1:length(xy)
    xn=xy(i,1)*cos(d)-xy(i,2)*sin(d);
    yn=xy(i,1)*sin(d)+xy(i,2)*cos(d);
    xy(i,1)=xn;
    xy(i,2)=yn;
end
line(xy(:,1),xy(:,2))

%% %%%%%%%%%%%%%%transfer  meter to degree%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
%close all;
xy=[-2.2e5,2.8e5,2.8e5,-2.2e5,-2.2e5;
    -3.5e5,-3.5e5,3e5,3e5,-3.5e5]';
figure(3) %
line(xy(:,1),xy(:,2)); hold on

%% first rotate "back" 45 degrees clockwise, 
d=45*pi/180;
for i=1:length(xy)
    xn=xy(i,1)*cos(d)-xy(i,2)*sin(d);
    yn=xy(i,1)*sin(d)+xy(i,2)*cos(d);
    xy(i,1)=xn;
    xy(i,2)=yn;
end
line(xy(:,1),xy(:,2)); hold on

%% transfer meter to degree
    x0=-117.774; y0=34.736; % local(0,0)
for i=1:length(xy)
    longlat(i,2)=y0+xy(i,2)/111263;
    longlat(i,1)=x0+xy(i,1)/(111263*abs(cos((pi/180)*(y0+longlat(i,2))/2)));
end
figure(4)
line(longlat(:,1),longlat(:,2)); hold on


%% %%%%%%%% output for saf627bddegree.txt
% 

fid_ans_out=fopen('saf627bddegree.txt','w');
fprintf(fid_ans_out,'%s\n','> saf627 boundary degree');

for i=1:length(longlat)
    fprintf(fid_ans_out,'%f%c%f\n',longlat(i,1),' ',longlat(i,2));    
end

fclose(fid_ans_out);






