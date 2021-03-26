%% Part B
clc; clear all;
%Defining Known Quantities
dx = 0.0125; dy = 0.0125; dt = 0.0125;
alpha = 0.1; a = 2; omega = 50;

%Defining X,Y, and T nodes
%final T chosen S.T. steady state is achieved
x = 0:dx:1; y = 0:dy:1; t = 0:dt:10;
M = length(x); N = length(y); tau = length(t);

%Creating Tridiagonal Matrices
beta = (alpha*dt)/(2*dx^2);

d1 = (1+(2*beta))*ones(M-2,1);
d2 = -beta*ones(M-3,1);
A = diag(d1) + diag(d2,1) + diag(d2,-1);

%defining Temperature Arrays
T0 = zeros(M,N);
for i = 1:M
    T0(i,:) = 0.01 .* sin(pi .* x(i)).* sin(pi.*y);
end

%enforcing BCs and creating matrix to store time evolution of T(0.55,0.45)
T0(1,:) = 0; T0(M,:) = 0; T0(:,1) = 0; T0(:,N) = 0;
Thalf = T0; T = T0; Tn = T0;
indX = find(x==0.55); indY = find(y==0.45);
tPlot = zeros(tau,1);
tPlot(1) = T0(indX,indY);

for i = 1:tau
    for k = 2:N-1
        b = (beta.*T(2:M-1,k-1)) + ((1-2*beta).*T(2:M-1,k)) + (beta.*T(2:M-1,k+1)) + ((dt/2) .* Q(x(2:M-1),y(k),t(i),a,omega)');
        Thalf(2:M-1,k) = thomas(A,b);
    end
    for j = 2:M-1
        b = (beta*Thalf(j-1,2:N-1)) + ((1-2*beta)*Thalf(j,2:N-1)) + (beta*Thalf(j+1,2:N-1)) + ((dt/2) * Q(x(j),y(2:N-1),t(i),a,omega));
        temp = thomas(A,b');
        Tn(j,2:N-1) = temp';
    end
    tPlot(i) = Tn(indX,indY);
    T = Tn;
end

%% Part D setting T0 = 0
clc; clear all;
%Defining Known Quantities
dx = 0.0125; dy = 0.0125; dt = 0.0125;
alpha = 0.1; a = 2; omega = 50;

%Defining X,Y, and T nodes
%final T chosen S.T. steady state is achieved
x = 0:dx:1; y = 0:dy:1; t = 0:dt:10;
M = length(x); N = length(y); tau = length(t);

%Creating Tridiagonal Matrices
beta = (alpha*dt)/(2*dx^2);

d1 = (1+(2*beta))*ones(M-2,1);
d2 = -beta*ones(M-3,1);
A = diag(d1) + diag(d2,1) + diag(d2,-1);

%defining Temperature Arrays
T0 = zeros(M,N);

%enforcing BCs and creating matrix to store time evolution of T(0.55,0.45)
T0(1,:) = 0; T0(M,:) = 0; T0(:,1) = 0; T0(:,N) = 0;
Thalf = T0; T = T0; Tn = T0;
indX = find(x==0.55); indY = find(y==0.45);
tPlot = zeros(tau,1);
tPlot(1) = T0(indX,indY);

for i = 1:tau
    for k = 2:N-1
        b = (beta.*T(2:M-1,k-1)) + ((1-2*beta).*T(2:M-1,k)) + (beta.*T(2:M-1,k+1)) + ((dt/2) .* Q(x(2:M-1),y(k),t(i),a,omega)');
        Thalf(2:M-1,k) = thomas(A,b);
    end
    for j = 2:M-1
        b = (beta*Thalf(j-1,2:N-1)) + ((1-2*beta)*Thalf(j,2:N-1)) + (beta*Thalf(j+1,2:N-1)) + ((dt/2) * Q(x(j),y(2:N-1),t(i),a,omega));
        temp = thomas(A,b');
        Tn(j,2:N-1) = temp';
    end
    tPlot(i) = Tn(indX,indY);
    T = Tn;
end
%% Plotting time

plot(t,tPlot);
xlabel("Time (T)")
ylabel("Temperature at (0.55,0.45)")
title("Temperature vs Time")
%% Plotting Contour Plot
[X,Y] = meshgrid(x,y);
contourf(X,Y,T)
xlabel("X Position")
ylabel("Y Position")
title("Temperature Contour Plot")
colorbar