%% Part A
clc; clear all;
dt = 0.0001;
dx = 0.01;
x = 0:dx:5; t = 0:dt:10;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
idx = find(x==2);
u0 = zeros(N,1);
u0(1:idx) = exp((-(x(1:idx)-1).^2)/0.18);
u0(1) = 0;
u=u0;
un = u;
usaved = zeros(5,N);
for i = 1:length(t)
    for j = 2:N-1 %inner node solution using 2nd order central FDA
        f1 = odefcn((u.^2)                  ,j,dx);
        f2 = odefcn((u.^2)+(a21*f1)         ,j,dx);
        f3 = odefcn((u.^2)+(a31*f1)+(a32*f2),j,dx);
        
        un(j) = u(j) + dt*(b1*f1 + b2*f2 + b3*f3);
    end

    %solving end node using 2nd order backward FDA
    f1 = odefcn3(u.^2                  ,N,dx);
    f2 = odefcn3(u.^2+(a21*f1)         ,N,dx);
    f3 = odefcn3(u.^2+(a31*f1)+(a32*f2),N,dx);
    un(N) = u(N) + dt*(b1*f1 + b2*f2 + b3*f3);
    %enforcing BC
    %un(1) = 0;
    %saving desired time stamps
    if t(i) == 0.5
        usaved(1,:) = un;
    elseif t(i) == 1.0
        usaved(2,:) = un;
    elseif t(i) == 1.5
        usaved(3,:) = un;
    elseif t(i) == 5.0
        usaved(4,:) = un;
    elseif t(i) == 10.0
        usaved(5,:) = un;
    end
    u=un;
end

plot(x,u0);
hold on;
plot(x,usaved);
xlabel("Distance (x)")
ylabel("Velocity (u)")
title("Second Order Central Difference Schema")
legend("t = 0.0", "t = 0.5", "t = 1.0", "t = 1.5", "t = 5.0", "t = 10.0")
%% Part B
clc; clear all;
dt = 0.001;
dx = 0.01;
x = 0:dx:5; t = 0:dt:10;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
idx = find(x==2);
u0 = zeros(N,1);
u0(1:idx) = exp((-(x(1:idx)-1).^2)/0.18);
u0(1) = 0;
u=u0;
un = u;
usaved = zeros(5,N);
for i = 1:length(t)
    for j = 2:N %inner node solution
        f1 = odefcn2((u.^2)                  ,j,dx);
        f2 = odefcn2((u.^2)+(a21*f1)         ,j,dx);
        f3 = odefcn2((u.^2)+(a31*f1)+(a32*f2),j,dx);
        un(j) = u(j) + dt*(b1*f1 + b2*f2 + b3*f3);
    end
    %solving for end node using linear extrapolation
    un(1) = 0;
    %saving desired time stamps
    if t(i) == 0.5
        usaved(1,:) = un;
    elseif t(i) == 1.0
        usaved(2,:) = un;
    elseif t(i) == 1.5
        usaved(3,:) = un;
    elseif t(i) == 5.0
        usaved(4,:) = un;
    elseif t(i) == 10.0
        usaved(5,:) = un;
    end
    u=un;
end

plot(x,u0);
hold on;
plot(x,usaved);
xlabel("Distance (x)")
ylabel("Velocity (u)")
title("First Order Upwind Schema")
legend("t = 0.0", "t = 0.5", "t = 1.0", "t = 1.5", "t = 5.0", "t = 10.0")