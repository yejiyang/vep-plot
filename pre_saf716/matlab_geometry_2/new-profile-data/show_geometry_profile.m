clear;
clc;
%close all;
%try read from shapefile directly

%% %%%%%%%%%%% read fault data from arcgis10 export txt file%%%%
%  STEP 1
%%%%% data format%%%%
% Polyline
% 0 0   %%% should start from 0
% 0 -118.861502657 34.8259605674 nan nan
% 1 -118.542735491 34.9779232635 nan nan
% 2 -118.029222225 35.2403497367 nan nan
% 3 -117.923713511 35.3601914946 nan nan
% 4 -117.138593616 35.5829191172 nan nan
% 5 -116.373439099 35.5931625201 nan nan
% 1 0
% 0 -117.407065279 34.2401437074 nan nan
% 1 -117.015484257 33.7977211982 nan nan
% 2 -115.972900224 33.0203162787 nan nan
% 3 -115.508269749 32.8374844402 nan nan
% END
%%%%%%%%%%%%%%%%%%
%%%%%% only used when total faults less than 99 %%%%%%%%%

%fidin=fopen('saf1.txt');    % open saf1.txt file
fidin=fopen('saf-profile-line-no-splitline-new.txt');    % open saf1.txt file
fidout = fopen('MKMATLAB.txt','w'); % creat MKMATLAB.txt file
while ~feof(fidin)  % determine whether the end of file
    tline=fgetl(fidin); % read one line from file
    if double(tline(1))>=48 && double(tline(1))<=57  % determine number or str for one line
     
        % get the id of faults
        if length(tline)==3 || length(tline)==4 % the number of fault must less than 99
            %faultid = double(tline(1))-48
            faultid = eval(tline(1:2));
        end
        
        % record the fault segements, add the fault id first
        if length(tline)>=5
            fprintf(fidout,'%f%s\n',faultid,tline); % if number, write to MKMATLAB.txt
        end
        
        continue        
    end
end
fclose(fidout);
saf_fault=importdata('MKMATLAB.txt');  % import the MKMATLAB into the workspace

% counting the segments for each faults
j=1;
n=1;
for i=2:length(saf_fault)
    if int8(saf_fault(i,1))==int8(saf_fault(i-1,1))
        j=j+1;
    end
    if int8(saf_fault(i,1))~=int8(saf_fault(i-1,1)) || i==length(saf_fault)
        k(n)=j;  % save the segments in k
        n=n+1;
        j=1;
    end
end

n=1;
m=k(1); 
for i=1:length(saf_fault)  % write the number into fourth column
    if i<=m
    saf_fault(i,4)=k(n);
    else 
        n=n+1;
        m=m+k(n);
        saf_fault(i,4)=k(n);
    end
    saf_fault(i,1)=int8(saf_fault(i,1));
end
saf_fault;

% show the faults in the figure
figure(1)
j=1;
for i=1:saf_fault(length(saf_fault),1)+1   
    X=(saf_fault(j:j+saf_fault(j,4)-1,2));
    Y=(saf_fault(j:j+saf_fault(j,4)-1,3));
    line(X,Y);
    j=j+saf_fault(j,4);
end
%





