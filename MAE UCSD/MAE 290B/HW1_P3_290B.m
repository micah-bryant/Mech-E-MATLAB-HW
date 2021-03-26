%% Part A
clc; clear all;

v = -3:0.01:3;  % plotting range from -3 to 3
[x,y] = meshgrid(v);  % get 2-D mesh for x and y 
realN = (1 + x + ((x.^2-y.^2-(x.*(y.^2)))./2) + (x.^3)./6).^2; %defining real section
imagN = (y+(x.*y)+(((x.^2).*y)./2) - (y.^3)./6).^2; %defining imaginary section
cond1 = realN+imagN < 1;  % check conditions for these values
cond1 = double(cond1);  % convert to double for plotting
cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
surf(x,y,cond1)
xlim([-3 3])
ylim([-3 3])
xlabel('h * \lambda R');
ylabel('h * \lambda I');
title('Stability plot for RK3')
view(0,90)    % change to top view
%% Part B
clc; clear all;

%establishing given variables
h = 1; y0 = 1; w = 0.2; tau = 400;
tspan = linspace(0,500,(500/h)+1);

%defining lambda as the coeffecient in dydt
lambda = (-1/tau)+1i*w;

%determining exact y value
y_exact = exp(lambda.*tspan);

%Setting up y values to be solved for
y = zeros(1,length(tspan));
y(1) = y0;

%setting up RK3 scheme
a21 = 8/15; a31 = 1/4; a32 = 5/12;
b1 = 1/4; b2 = 0; b3 = 3/4;

%solving step by step for the y values
for i = 2:length(tspan)
    f1 = lambda*(y(i-1));
    f2 = lambda*(y(i-1) + a21*h*f1);
    f3 = lambda*(y(i-1) + a31*h*f1 + a32*h*f2);
    y(i) = y(i-1) + h*(b1*f1 + b2*f2 + b3*f3);
end

subplot(2,1,1);
plot(tspan(1:(50/h)+1),y(1:(50/h)+1))
hold on;
plot(tspan(1:(50/h)+1),y_exact(1:(50/h)+1))
xlabel('t');
ylabel('y');
title('y vs t with h=1')
legend('y using RK3', 'y exact');

subplot(2,1,2);
plot(tspan((451/h)+1:(500/h)+1),y((451/h)+1:(500/h)+1))
hold on;
plot(tspan((451/h)+1:(500/h)+1),y_exact((451/h)+1:(500/h)+1))
xlabel('t');
ylabel('y');
title('y vs t with h=1')
legend('y using RK3', 'y exact');
%% Part C
h = [0.01 0.1 0.5 1]; y0 = 1; w = 0.2; tau = 400;
error = zeros(length(h),1);
for j = 1:length(h)
    tspan = linspace(0,500,(500/h(j))+1);
    
    %defining lambda as the coeffecient in dydt
    lambda = (-1/tau)+1i*w;
    
    %determining exact y value
    y_exact = exp(lambda.*tspan);
    
    %Setting up y values to be solved for
    y = zeros(1,length(tspan));
    y(1) = y0;
    
    %setting up RK3 scheme
    a21 = 8/15; a31 = 1/4; a32 = 5/12;
    b1 = 1/4; b2 = 0; b3 = 3/4;
    
    %solving step by step for the y values
    for i = 2:length(tspan)
        f1 = lambda*(y(i-1));
        f2 = lambda*(y(i-1) + a21*h(j)*f1);
        f3 = lambda*(y(i-1) + a31*h(j)*f1 + a32*h(j)*f2);
        y(i) = y(i-1) + h(j)*(b1*f1 + b2*f2 + b3*f3);
    end
    error(j) = y_exact(end) - y(end);
    error(j) = abs(error(j));
end

loglog(h,error); %visually checking to ensure it is linear
disp((log(error(end))-log(error(1)))/(log(h(end))-log(h(1)))); %determining slope of line