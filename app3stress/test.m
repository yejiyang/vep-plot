% A1 = [9.9, 9900];
% A2 = [8.8,  7.7 ; ...
%       8800, 7700];
%   
%  formatSpec = 'X is %4.2f meters or %8.3f mm\n';
%  fprintf(formatSpec,A1,A2)
% 
% A3 = [8.8, 7.7 ; 
%     8800, 7700; 
%     1, 1; 
%     2, 2];
% formatSpec1 = '%4.2f %4.2f \n';
% fprintf(formatSpec1,A3)

x = 0:.1:1;
A = [x; exp(x)];

fileID = fopen('exp.txt','w');
fprintf(fileID,'%6s %12s\n','x','exp(x)');
fprintf(fileID,'%6.2f %12.8f\n',A);
fclose(fileID);