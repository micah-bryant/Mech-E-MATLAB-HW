clc; clear all;
d = linspace(1,100,100);
sDiag = ones(99,1); b = ones(100,1);
A = diag(d) + diag(sDiag,1) + diag(sDiag,-1);

k = cond(A);

rc = zeros(100,100); rs = zeros(100,100); %residual values
xc = zeros(100,100); xs = zeros(100,100); %x vector values

rcNorm = zeros(100,1); rsNorm = zeros(100,1);
xcNorm = zeros(100,1); estimate = zeros(100,1);

for i = 1:100
    [xc(:,i),rc(:,i)] = CG(A,b,i);
    [xs(:,i),rs(:,i)] = steepDescent(A,b,i);
    rcNorm(i) = norm(rc(:,i));
    rsNorm(i) = norm(rs(:,i));
    xcNorm(i) = norm(b-(A*xc(:,i)));
    estimate(i) = 2*((sqrt(k)-1)^i)/((sqrt(k)+1)^i);
end

nVec = linspace(1,100,100);
p1 = semilogy(nVec, rcNorm, 'b');
hold on
p2 = semilogy(nVec, rsNorm, 'r');
hold on
p3 = semilogy(nVec, xcNorm, 'g');
hold on
p4 = semilogy(nVec, estimate, 'm');
legend([p1,p2,p3,p4],{'CG Residual Norm', 'Steep Descent Residual Norm', 'CG Actual Residual Norm', 'Estimate'},'Location', 'southwest');