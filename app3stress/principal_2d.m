function [output] = principal_2d(sij)
%function [p,sita_p,max_shear] = principal_2d(sij)

%% about this script
% Author: Jiyang Ye (yejiyang@gmail.com)
% initial date: Jul 29, 2014
% revised date: Oct 14, 2014
% this program aim to calculate 2-d principal stress, same with strain, 
% input:  the 2-d stress, sij = [s11 s12; s21 s22], s12 = s21; 
% output: including, 4 float numbers
% output: output[1 2]: the principal stress I & II
% output: output[3]: sitaP, the direction of max principal stress counter-clockwise 
%           from x direction.
% output: output[4]: the maxshear stress
% caution:  all the unit of the angle are degrees!
%clear; clc; 

%% test data
% 
% [output]=   principal_2d(sij)
%
% test cases, compare the results with
% http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/calc_principal_stress.cfm#calc
%     % max 0.611, deg 20.3    
%sij = [0.5 0.3;   0.3 -0.2];
% 
%     %max 70, deg 26.5651
%sij = [50 40; 40 -10 ];
%     
%         %max 70, deg 63.4 , ok 
%sij = [-10 40;    40 50 ];
%     
%     % max 0, deg -45; min -8, deg 45, ok
%sij = [-4 -4 ;    -4 -4];
%     
% % max 11, deg 18.4 , ok       
%sij = [10 3;  3 2];
% 
%     % max 1490, deg -13.2, max shear 785, directon, sita_s1, 31.8, sita-s2, -58.2 ok
%sij = [1406 -350; -350 0];
% 
%     
%         %max 4.5, deg 40, ok    
%sij = [2 3;   3 1];
%     
%         
%     %max 61.1, deg 20.3, ok
%sij = [50 30; 30 -20]; 
%     
% 
%     %max 1.65, deg 32.9
%sij = [1 1;   1 0.1];
% 
%


%% using eigenvalues and eigenvectors
% using eig() or eigs(), can quickly get the values
% original_stress
%sij = [12000 8000; 8000 15000];
%
% v: vectors; p: principal stress or strain
[v, p] = eigs(sij); %Eigen value function in MATLAB, using eigs not eig

% Get the principal stress direction from eigenvector.
% Here, we got the I1 (max principal stress) direction, counter-clockwise
% rotate 45 from x xias.
temp = v(2,1) / v(1,1); % Slope of eigenvector is a ratio of its (y/x)-components
%sita_p = atan(temp);
%sita_p = sita_p * (180/pi); %Principal stress orientation in degrees
sita_p  = atand(temp);

% max_shear: max shear stress or strain
max_shear = (p(1,1) - p(2,2))/2;
% Max shear stress orientation in degrees, counter-clockwise rotate 45 from principal stress I
%sita_I  = sita_p + 45; 
% Max shear stress orientation in degrees, clockwise rotate 45 from principal stress II
%sita_II = sita_p - 45;

output = [p(1,1), p(2,2), sita_p, max_shear];



%% to check the results
% here, using the equation, we first calculate the principal direction,
% but don't know whether it is max or min direction;
% then, construct the transformation function to check the direction for
%        max principal stress
% 
%  refs http://www.continuummechanics.org/cm/principalstress.html
%sitaP = atand(2*sij(1,2)/(sij(1,1)-sij(2,2)))/2;
% 
% % get the transformation matrix
%T = [cosd(sitaP) sind(sitaP); 
%     -sind(sitaP) cosd(sitaP)];
% 
% % get the initial principal stress, using the transformation matrix
%pij_init = T*sij*T';
% 
% % reorder the initial principal stress, so pij(1,1) = max_pincipal_stress
% % if pij_init(1,1) is not the max, then swith with the pij_init(2,2)
% % also swith the principal direction angle
%if pij_init(1,1) < pij_init(2,2)
%     pij=zeros(2);
%     pij(1,1)=pij_init(2,2);
%     pij(2,2) = pij_init(1,1);
%     
%     sitaP_max = sitaP + 90;
% % else, nothing change
%else 
%     pij = pij_init;
%     sitaP_max = sitaP;
%     
%end

%% end of function
end


