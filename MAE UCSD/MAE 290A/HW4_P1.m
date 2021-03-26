clc; clear all;
load('E.mat');
load('A.mat');
arnMat = A;
nMax = 50;
qrNorm = zeros(nMax,1); iNorm = zeros(nMax,1); qhqNorm = zeros(nMax,1);
nVec = linspace(1,nMax,nMax);
for n = 1:nMax
    [Q,H] = arnoldi(arnMat,n);
    Qn = Q(:,1:n); Hn = H(1:n,:);
    Rn = Hn*Qn';
    [Q1,R1] = qr(Rn);
    I = eye(size(Qn'*Qn));
    qrNorm(n) = norm(arnMat-(Qn*Q1*R1)); 
    iNorm(n) = norm((Qn'*Qn)-I);
    qhqNorm(n) = norm(arnMat-(Qn*Hn*Qn'));
end

subplot(3,1,1)
plot(nVec, qrNorm);
title('QR Norm')
subplot(3,1,2)
plot(nVec, iNorm);
title('Identity Norm')
subplot(3,1,3)
plot(nVec, qhqNorm);
title('Hessenberg Norm')