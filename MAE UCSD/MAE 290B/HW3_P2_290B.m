%% Part A
clc; clear all;
kh = linspace(0,2*pi, 1000);
dkh1 = sin(kh);
dkh2 = (2.*sin(kh/2));
plot(kh,dkh1);
hold on;
plot(kh, dkh2);
hold on;
plot(kh, kh);
xlabel('kh');
ylabel('k''h');
title('Plot of k''h vs kh')
legend('dT/dx', 'd^2T/dx^2', 'Actual kh');

%% Part Di using linear extrapolation
clc; clear all;
u = 0.4; a = 0;
dt = 0.00005;
dx = 0.002;
x = 0:dx:1; t = 0:dt:2;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
T0 = sin(pi.*x) .* (cos(4 .* pi .* x) + sin(20 .* pi .* x));
T=T0;
T(1) = 0; T(end) = 0;
Tn = T;
Tsaved = zeros(3,N);
for i = 1:length(t)
    for j = 2:N-1 %inner node solution
        f1 = odefcn(T                  ,j,dx,a,u);
        f2 = odefcn(T+(a21*f1)         ,j,dx,a,u);
        f3 = odefcn(T+(a31*f1)+(a32*f2),j,dx,a,u);
        Tn(j) = T(j) + dt*(b1*f1 + b2*f2 + b3*f3);
    end
    %solving for end node using linear extrapolation
    Tn(N) = 2*Tn(N-1)-Tn(N-2);
    %saving desired time stamps
    if t(i) == 0.5
        Tsaved(1,:) = Tn;
    elseif t(i) == 1.0
        Tsaved(2,:) = Tn;
    elseif t(i) == 2.0
        Tsaved(3,:) = Tn;
    end
    T=Tn;
end

subplot(2,2,1);
plot(x,T0);
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0 seconds")
subplot(2,2,2);
plot(x,Tsaved(1,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.5 seconds")
subplot(2,2,3);
plot(x,Tsaved(2,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 1 seconds")
subplot(2,2,4);
plot(x,Tsaved(3,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 2 seconds")
sgtitle('Convection w/Linear Extrapolation')

%% Part Di using first order backwards
clc; clear all;
u = 0.4; a = 0;
dt = 0.00005;
dx = 0.002;
x = 0:dx:1; t = 0:dt:2;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
T0 = sin(pi.*x) .* (cos(4 .* pi .* x) + sin(20 .* pi .* x));
T=T0;
T(1) = 0; T(end) = 0;
Tn = T;
Tsaved = zeros(3,N);
for i = 1:length(t)
    for j = 2:N %inner node solution
        if j == N %solves for the last node using backward difference
            f1 = odefcn2(T                  ,j,dx,a,u);
            f2 = odefcn2(T+(a21*f1)         ,j,dx,a,u);
            f3 = odefcn2(T+(a31*f1)+(a32*f2),j,dx,a,u);
            Tn(j) = T(j) + dt*(b1*f1 + b2*f2 + b3*f3);
        else %solves for all other nodes using central difference
            f1 = odefcn(T                  ,j,dx,a,u);
            f2 = odefcn(T+(a21*f1)         ,j,dx,a,u);
            f3 = odefcn(T+(a31*f1)+(a32*f2),j,dx,a,u);
            Tn(j) = T(j) + dt*(b1*f1 + b2*f2 + b3*f3);
        end
    end
    %saving desired time stamps
    if t(i) == 0.5
        Tsaved(1,:) = Tn;
    elseif t(i) == 1.0
        Tsaved(2,:) = Tn;
    elseif t(i) == 2.0
        Tsaved(3,:) = Tn;
    end
    T=Tn;
end

subplot(2,2,1);
plot(x,T0);
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0 seconds")
subplot(2,2,2);
plot(x,Tsaved(1,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.5 seconds")
subplot(2,2,3);
plot(x,Tsaved(2,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 1 seconds")
subplot(2,2,4);
plot(x,Tsaved(3,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 2 seconds")
sgtitle('Convection w/Backward Difference')

%% Part Dii
clc; clear all;
u = 0; a = 0.01;
dt = 0.0005;
dx = 0.02;
x = 0:dx:1; t = 0:dt:10;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
T0 = sin(pi.*x) .* (cos(4 .* pi .* x) + sin(20 .* pi .* x));
T0(1) = 0; T0(end) = 0.5;
T=T0;
Tn = T;
Tsaved = zeros(5,N);
for i = 1:length(t)
    for j = 2:N-1 %inner node solution
        f1 = odefcn(T                  ,j,dx,a,u);
        f2 = odefcn(T+(a21*f1)         ,j,dx,a,u);
        f3 = odefcn(T+(a31*f1)+(a32*f2),j,dx,a,u);
        Tn(j) = T(j) + dt*(b1*f1 + b2*f2 + b3*f3);
    end
    %saving desired time stamps
    if t(i) == 0.02
        Tsaved(1,:) = Tn;
    elseif t(i) == 0.1
        Tsaved(2,:) = Tn;
    elseif t(i) == 0.5
        Tsaved(3,:) = Tn;
    elseif t(i) == 1
        Tsaved(4,:) = Tn;
    elseif t(i) == 10
        Tsaved(5,:) = Tn;
    end
    T=Tn;
end

subplot(3,2,1);
plot(x,T0);
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0 seconds")
subplot(3,2,2);
plot(x,Tsaved(1,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.02 seconds")
subplot(3,2,3);
plot(x,Tsaved(2,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.1 seconds")
subplot(3,2,4);
plot(x,Tsaved(3,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.5 seconds")
subplot(3,2,5);
plot(x,Tsaved(4,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 1 seconds")
subplot(3,2,6);
plot(x,Tsaved(5,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 10 seconds")


sgtitle('Diffusion')

%% Exact Solution to Advection Equation
clc; clear all;
u = 0.4; a = 0;
dt = 0.00005;
dx = 0.002;
x = 0:dx:1; t = 0:dt:2;
N = length(x);
b1=1/4; b2=0; b3=3/4;
a21=8/15; a31=1/4; a32=5/12;
T0 = sin(pi.*x) .* (cos(4 .* pi .* x) + sin(20 .* pi .* x));
T=T0;
T(1) = 0; T(end) = 0;
Tn = T;
Tsaved = zeros(3,N);
for i = 1:length(t)
    for j = 2:N %inner node solution
        if x(j)-(u*t(i)) <=0
            Tn(j) = 0;
        else
            idx = find(x>=x(j)-(u*t(i)),1);
            Tn(j) = T0(idx);
        end
    end
    %saving desired time stamps
    if t(i) == 0.5
        Tsaved(1,:) = Tn;
    elseif t(i) == 1.0
        Tsaved(2,:) = Tn;
    elseif t(i) == 2.0
        Tsaved(3,:) = Tn;
    end
    T=Tn;
end

subplot(2,2,1);
plot(x,T0);
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0 seconds")
subplot(2,2,2);
plot(x,Tsaved(1,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 0.5 seconds")
subplot(2,2,3);
plot(x,Tsaved(2,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 1 seconds")
subplot(2,2,4);
plot(x,Tsaved(3,:));
xlabel("Distance (x)")
ylabel("Temperature (T)")
title("Time = 2 seconds")
sgtitle('Exact Solution')