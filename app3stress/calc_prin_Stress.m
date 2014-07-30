%%
% Author: Jiyang Ye (yejiyang@gmail.com)
% initial date: Jul 30, 2014
% this program aim to calculate 2-d principal stress, same with strain, 
% input:  the 2-d stress, sij = [s11 s12; s21 s22], s12 = s21; 
% output: the principal stress, pij = [max_sij 0; 0 min_sij], 
% output: sitaP_max, the direction of max principal stress countclockwise 
%           from x direction.
% caution:  all the unit of the angle are degrees!

%% begin
clear; clc;

%% read the initial 2-d stress state data
init_str = load('./input/input_stress_for_prin_stress');
%init_str = xlsread('./input/input_stress_for_prin_stress.csv');

%% call the function [pij, sitaP_max] = principal_Stress_trans(sij) to calculate 
% the principal stress
for i=1:length(init_str)
%for i=1:10
	sij = zeros(2);
	sij(1,1) = init_str(i,4);
	sij(2,2) = init_str(i,5);
	sij(1,2) = init_str(i,6);
	sij(2,1) = init_str(i,6);
	%sij
	
	% call the function
	[pij , sitaP_max] = principal_Stress_trans(sij);
	
	% storage the data
	pij_final(i,1) = init_str(i,1); % layer number	
	pij_final(i,2) = init_str(i,2); % x or longi
	pij_final(i,3) = init_str(i,3); % y or lat
	pij_final(i,4) = pij(1,1);      % max_principal_stress
	pij_final(i,5) = pij(2,2);      % min_principal_stress
	pij_final(i,6) = sitaP_max;     % countclockwise angle from x to max_principal_stress	
	
end

%pij_final

% output the results
csvwrite('principal_stress.data', pij_final); % output the results

%end
