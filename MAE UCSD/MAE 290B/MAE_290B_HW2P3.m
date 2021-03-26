clc; clear all;
%known values
N = 200; L = 2; T0 = 5; TN = 4;
h = L/(N-1);

%preallocating for tridiagonal and answers
alpha = zeros(N,1); beta = zeros(N,1); gamma = zeros(N,1); f = zeros(N,1); x = zeros(N,1); 
T = zeros(N,1);

%setting up values for tridiagonal matrix and answer vector
for i = 1:N
    x(i) = h*(i-1);
    a = -(x(i)+3)/(x(i)+1);
    b = (x(i)+3)/((x(i)+1)^2);
    alpha(i) = (1/(h^2)) - (a/(2*h));
    beta(i) = (-2/(h^2)) + b;
    gamma(i) = (1/(h^2)) + (a/(2*h));
    f(i) = (2*(x(i)+1)) + (3*b);
end

%setting up A and f such that A*T = f
A = diag(alpha(3:N-1),-1) + diag(beta(2:N-1)) + diag(gamma(2:N-2),1);
q = f;
f(2) = f(2) - (alpha(2)*T0);
f(N-1) = f(N-1) - (gamma(N-1)*TN);

w = beta(2);
t = T;
T(2) = f(2)/w;

%forward sweep
for i = 3:N-1
    t(i-1) = gamma(i-1)/w;
    w = beta(i) - alpha(i)*t(i-1);
    T(i) = (f(i) - alpha(i)*T(i-1))/w;
end

%inverse sweep
for j=N-2:-1:2
    T(j) = T(j) - t(j)*T(j+1);
end

T(1) = T0; T(N) = TN;
plot(x,T)
xlabel("Distance")
ylabel("Temperature")
title("Thomas Algorithm for T")