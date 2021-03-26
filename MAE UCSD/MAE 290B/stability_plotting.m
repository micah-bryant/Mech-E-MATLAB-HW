clc; clear all;

v = -3:0.01:3;  % plotting range from -3 to 3
[x,y] = meshgrid(v);  % get 2-D mesh for x and y
%num = (1+x).^2+y.^2; %for theta = 0
%num = (1+(x./2)).^2+(y./2).^2; %for theta = 0.5
num = 1; %for theta = 1

%den = 1; %for theta = 0
%den = (1-(x./2)).^2+(y./2).^2; %for theta = 0.5
den = (1-(x)).^2+(y).^2; %for theta = 1
cond1 = num./den < 1;  % check conditions for these values
cond1 = double(cond1);  % convert to double for plotting
cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
surf(x,y,cond1)
xlim([-3 3])
ylim([-3 3])
xlabel('h * \lambda R');
ylabel('h * \lambda I');
title('Stability plot for \theta = 1')
view(0,90)    % change to top view
