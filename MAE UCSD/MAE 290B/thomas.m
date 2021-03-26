function T = thomas(A,f)
% A is a tridiag matrix, b is a vector
alpha = diag(A); beta = diag(A,-1); gamma = diag(A,1);
N = length(f);
T = zeros(N,1);
w = alpha(1);
t = T;
T(1) = f(1)/w;

%forward sweep
for i = 2:N
    t(i-1) = gamma(i-1)/w;
    w = alpha(i) - beta(i-1)*t(i-1);
    T(i) = (f(i) - beta(i-1)*T(i-1))/w;
end

%inverse sweep
for j=N-1:-1:1
    T(j) = T(j) - t(j)*T(j+1);
end