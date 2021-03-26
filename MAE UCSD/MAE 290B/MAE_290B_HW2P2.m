%% Part A-B
clc; clear all;
tspan = linspace(0,1,1000);
y0 = [1.0, 1.0, 1.0, -10., 0.9904];

[t,y] = ode23s(@(t,y)odefcn(t,y),tspan,y0);
%% Part C
clc; clear all;
a = 100.; b = 0.9; c = 1000.; d = 10.;
J = [-a 0 a 0 0;
    a 0 -a 0 0;
    -a/10 0 (c+a)/10 -(b-(0.01*c))/100 -c/10;
    -a 0 a 0 0;
    0 0 c/d 0 -c/d];
[~, D] = eig(J);
d = abs(diag(D));
stiff = max(d)/min(d);
%% Part D
clc; clear all;
numPts = 100;
h = 1/numPts;
tspan = linspace(0,1,numPts);

%setting up y variables
y = zeros(5,length(tspan));
y0 = [1.0, 1.0, 1.0, -10., 0.9904];
y(:,1) = y0;

%setting up 
a21 = 1/2; a31 = -1; a32 = 2;
b1 = 1/6; b2 = 2/3; b3 = 1/6;

for i = 1:length(tspan)-1
    f1 = odefcn(tspan(i),y(:,i));
    f2 = odefcn(tspan(i),(y(:,i) + a21*h*f1));
    f3 = odefcn(tspan(i),(y(:,i) + a31*h*f1 + a32*h*f2));
    y(:,i+1) = y(:,i) + h*(b1*f1 + b2*f2 + b3*f3);
end

t = tspan;
y = y';
%% Plotting
subplot(3,2,1);
plot(t,y(:,1));
xlabel("Time");
ylabel("Y1 Values")
title("Y1 Plot")

subplot(3,2,2);
plot(t,y(:,2));
xlabel("Time");
ylabel("Y2 Values")
title("Y2 Plot")

subplot(3,2,3);
plot(t,y(:,3));
xlabel("Time");
ylabel("Y3 Values")
title("Y3 Plot")

subplot(3,2,4);
plot(t,y(:,4));
xlabel("Time");
ylabel("Y4 Values")
title("Y4 Plot")

subplot(3,2,5);
plot(t,y(:,5));
xlabel("Time");
ylabel("Y5 Values")
title("Y5 Plot")

subplot(3,2,6);
plot(t,y(:,1));
hold on;
plot(t,y(:,2));
hold on;
plot(t,y(:,3));
hold on;
plot(t,y(:,4));
hold on;
plot(t,y(:,5));

xlabel("Time");
ylabel("Y Values")
title("Y Plot")
legend(["Y1", "Y2", "Y3", "Y4", "Y5"], 'location', 'northwest')
