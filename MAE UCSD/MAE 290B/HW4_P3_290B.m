clc; clear all;
dx = 0.05; dy = 0.05; k = 10;
x = 0:dx:20; y = 0:dy:35;
M = length(x); N = length(y);
beta = (-1i*dy)/(4*k*dx^2);
phi = zeros(M,N);
phi(:,1) = exp(-((x-5).^2)/4) + exp((-((x-15).^2)/4) + (10i.*x));
%ensuring BCs
phi(1,1) = 0;
phi(M,1) = 0;
%creating A matrix
d1 = (1+(2*beta))*ones(M-2,1);
d2 = -beta*ones(M-3,1);
A = diag(d1) + diag(d2,1) + diag(d2,-1);

for n = 1:N-1
    b = (beta*phi(3:M,n)) + ((1-2*beta)*phi(2:M-1,n)) + (beta*phi(1:M-2,n));
    phi(2:M-1,n+1) = A\b;
    %ensuring BC are enforced
    phi(1,n+1) = 0;
    phi(M,n+1) = 0;
end

magPhi = (abs(phi).^2)';
h = pcolor(x,y,magPhi);
set(h,'EdgeColor','none');
colorbar
title("Paraxial Helmholtz Equation")