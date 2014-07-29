function [pij, sitaP_max] = principal_Stress_trans(sij)

%%
% Author: Jiyang Ye (yejiyang@gmail.com)
% initial date: Jul 29, 2014
% this program aim to calculate 2-d principal stress, same with strain, 
% input:  the 2-d stress, sij = [s11 s12; s21 s22], s12 = s21; 
% output: the principal stress, pij = [max_sij 0; 0 min_sij], 
% output: sitaP_max, the direction of max principal stress countclockwise 
%           from x direction.
% caution:  all the unit of the angle are degrees!
% clear; clc; 

%% test data
% [pij, sita]=   principal_Stress_trans(sij)
%
% test cases, compare the results with
% http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/calc_principal_stress.cfm#calc
%     % max 0.611, deg 20.3    
% sij = [0.5 0.3;
%         0.3 -0.2];
% 
%     %max 70, deg 26.5651
% sij = [50 40;
%         40 -10 ];
%     
%         %max 70, deg 63.4 , ok 
% sij = [-10 40;
%         40 50 ];
%     
%     % max 0, deg -45; min -8, deg 45, ok
% sij = [-4 -4 ;
%         -4 -4];
%     
% % max 11, deg 18.4 , ok       
% sij = [10 3;
%         3 2];
% 
%     % max 1490, deg 13.2, ok
% sij = [1406 -350;
%         -350 0];
% 
%     
%         %max 4.5, deg 40, ok    
% sij = [2 3;
%         3 1];
%     
%         
%     %max 61.1, deg 20.3
% sij = [50 30;
%         30 -20]; 
%     
% 
%     %max 1.65, deg 32.9
% sij = [1 1;
%         1 0.1];
% [pij, sita]=   principal_Stress_trans(sij)
    
%% using eigenvalues and eigenvectors
% using eig() or eigs(), can quickly get the values
% but, have some problems about the determination of direction
%[V, D] = eigs(sij)
%principal_angel_a1 = acosd (V(1,1))

%% here, using the equation, we first calculate the principal direction,
% but don't know whether it is max or min direction;
% then, construct the transformation function to check the direction for
%       max principal stress

% refs http://www.continuummechanics.org/cm/principalstress.html
sitaP = atand(2*sij(1,2)/(sij(1,1)-sij(2,2)))/2;

% get the transformation matrix
T = [cosd(sitaP) sind(sitaP); 
    -sind(sitaP) cosd(sitaP)];

% get the initial principal stress, using the transformation matrix
pij_init = T*sij*T';

% reorder the initial principal stress, so pij(1,1) = max_pincipal_stress
% if pij_init(1,1) is not the max, then swith with the pij_init(2,2)
% also swith the principal direction angle
if pij_init(1,1) < pij_init(2,2)
    pij=zeros(2);
    pij(1,1)=pij_init(2,2);
    pij(2,2) = pij_init(1,1);
    
    sitaP_max = sitaP + 90;
% else, nothing change
else 
    pij = pij_init;
    sitaP_max = sitaP;
    
end

% end of function
end


