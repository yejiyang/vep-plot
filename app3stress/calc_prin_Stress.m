%% about this script
% Author: Jiyang Ye (yejiyang@gmail.com)
% initial date: Jul 30, 2014
% revised Oct 14, 2014
% this program aim to calculate 2-d principal stress, same with strain, 
% input:  the 2-d stress, sij = [s11 s12; s21 s22], s12 = s21; 
% output: the principal stress, pij = [max_sij 0; 0 min_sij], 
% output: sitaP_max, the direction of max principal stress countclockwise 
%           from x direction.
% caution:  all the unit of the angle are degrees!

%% begin
clear; clc;

%% read the initial 2-d stress state data
ori_str = load('./input/input_stress_for_prin_stress');
%ori_str = xlsread('./input/input_stress_for_prin_stress.csv');

%% call the function [output] = principal_2d(sij) to calculate 
% the principal stress;
size_ori = size(ori_str);
% the output size of function principal_2d equals 4;
size_out = 4;
%pij =zeros(size_ori(1),size_ori(2)+size_out);

pij =zeros(size_ori(1),size_out);
%pij=zeros(10,4);

for i=1:length(ori_str)
%for i=1:10
	sij = zeros(2);
	sij(1,1) = ori_str(i,4);
	sij(2,2) = ori_str(i,5);
	sij(1,2) = ori_str(i,6);
	sij(2,1) = ori_str(i,6);
	%sij
	
	% call the function
	[output] = principal_2d(sij);
	
	% get the original data
	%pij(i,1:size_ori(2)) = ori_str(i,1:size_ori(2)); % layer number	
    % get the principal data
	%pij(i, (size_ori(2)+1):(size_ori(2)+size_out)) = output;
    pij(i, 1:size_out) = output;
	
end


%% output the results

%csvwrite('principal_stress.data', pij); % output the results

%save prin_stress.out pij -ASCII

fileID = fopen('prin_stress.txt','w');
%fprintf(fileID,'%6s %12s\n','x','exp(x)');
% fprintf write in column order!!!, so pij ->pij'
formatSpec='%4.2f   %4.2f   %3.2f   %4.2f\n';
fprintf(fileID,formatSpec,pij');
fclose(fileID);

%% end
