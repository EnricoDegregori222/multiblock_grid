% xy = dlmread('v2c.txt');
% 
% c = 200;
% cd = 1.0;
% x = xy(:,1)-c/2;
% y = xy(:,3);
% xd = x * cd/c;
% yd = y * cd/c;
% 
% fileID = fopen('v2c_d.txt','w');
% for i=1:length(xd)
%     fprintf(fileID,'%e %e\n', xd(i), yd(i));
% end
% fclose(fileID);

fileID = fopen('v2c.dat','w');
xy = dlmread('airfoil_lower.txt');
for i=length(xy):-1:1
    fprintf(fileID,'%e %e\n', xy(i,1), xy(i,2));
end
xy = dlmread('airfoil_upper.txt');
for i=1:length(xy)
    fprintf(fileID,'%e %e\n', xy(i,1), xy(i,2));
end
fclose(fileID);