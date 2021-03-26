%%
clc; clear all;
load('E.mat');
lanMat = E;

nMax = 50;
qrNorm = zeros(nMax,1); iNorm = zeros(nMax,1); qhqNorm = zeros(nMax,1);
nVec = linspace(1,nMax,nMax);

for n = 1:nMax
    [Q,H] = lanczos(lanMat,n);
    R = H*Q';
    [Q1,R1] = qr(R);
    I = eye(size(Q'*Q));
    qrNorm(n) = norm(lanMat-(Q*Q1*R1)); 
    iNorm(n) = norm((Q'*Q)-I);
    qhqNorm(n) = norm(lanMat-(Q*H*Q'));
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
%% Part C
nTime = 30;
tic
for i = 1:20
    [Q,H] = lanczos(lanMat,nTime);
end
t = toc;
t = t/20;
%% Part D
A = zeros(1000,1000);
diagonal = sqrt(linspace(1,1000,1000));
for i = 1:1000
    A(i,i) = diagonal(i);
    if i<1000
        A(i+1,i) = 1;
        A(i,i+1) = 1;
    end
    if i<901
        A(i,i+100) = 1;
        A(i+100,i) = 1;
        disp(i)
    end
end

[Q,H] = lanczos(A,50);
[~,D] = eig(H);
value = min(diag(D));