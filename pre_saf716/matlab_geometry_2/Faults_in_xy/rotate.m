clear;
clc;
close all;


%% load the fault data in deg format
fault_deg = load('nv.deg');
%longlat=[-121.6,-113.2,-113.2,-121.6,-121.6; 31.4,31.4,38.1,38.1,31.4]';
%fault_deg=longlat;      
% show it
%plot(fault_deg(:,1),fault_deg(:,2));
%


%% %%%%%%%%%%%%%%% transfer the degree to meter
%%%%%%%%% set a local center, using the center of the region as local(0,0) 
%  STEP 2
    x0=-117.774; y0=34.736; % local(0,0)
for i=1:length(fault_deg)
    x=fault_deg(i,1); y=fault_deg(i,2);
    xn=111263*abs(cos((pi/180)*((y+y0)/2)))*(x-x0);
    yn=111263*(y-y0);
    fault_deg(i,1)=xn;
    fault_deg(i,2)=yn;
end

% show the faults in the figure
%figure(2)
%hold on
%plot(fault_deg(:,1),fault_deg(:,2),'r*-');
%

%% %%%%%%% rotate 45 degrees clockwise, 320 anticlockwise
%  STEP 3
d=-45;
xn=0;yn=0;
for i=1:length(fault_deg)
    x=fault_deg(i,1); y=fault_deg(i,2);
    xn=x*cosd(d)-y*sind(d);
    yn=x*sind(d)+y*cosd(d);
    fault_deg(i,1)=xn;
    fault_deg(i,2)=yn;
end

% show the faults in the figure
%figure(3)
%plot(fault_deg(:,1),fault_deg(:,2),'*-');
%hold off
%

%% %%%%%%%% output for gmt
%  STEP4

fid_ans_out=fopen('nv.km','w');
fprintf(fid_ans_out,'%s\n','>nv.km');
fprintf(fid_ans_out,'%3.2f  %3.2f\n',fault_deg'/1000);
fclose(fid_ans_out);


