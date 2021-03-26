%% Part A
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
%% Part C
clc; clear all;

%establishing given variables
h = 1; y0 = 1; w = 0.2; tau = 400; theta = 0.25;
tspan = linspace(0,500,(500/h)+1);

%defining lambda as the coeffecient in dydt
lambda = (-1/tau)+1i*w;

%determining exact y value
y_exact = exp(lambda.*tspan);

%Setting up y values to be solved for
y = zeros(1,length(tspan));
y(1) = y0;

%solving for coeffecient in the scheme
num = 1 + (1-theta)*h*lambda;
den = 1- h*theta*lambda;
coeff = num/den;

%solving step by step for the y values
for i = 2:length(tspan)
    y(i) = y(i-1) * coeff;
end

subplot(2,1,1);
plot(tspan(1:(50/h)+1),y(1:(50/h)+1))
hold on;
plot(tspan(1:(50/h)+1),y_exact(1:(50/h)+1))
xlabel('t');
ylabel('y');
title('y vs t with h=1')
legend('y using theta=0.75', 'y exact');

subplot(2,1,2);
plot(tspan((451/h)+1:(500/h)+1),y((451/h)+1:(500/h)+1))
hold on;
plot(tspan((451/h)+1:(500/h)+1),y_exact((451/h)+1:(500/h)+1))
xlabel('t');
ylabel('y');
title('y vs t with h=1')
legend('y using theta=0.75', 'y exact');