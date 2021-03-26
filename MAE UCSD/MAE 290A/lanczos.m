function [Q,T] = lanczos(A,n)
[~,m] = size(A);
b = rand(m,1); beta = zeros(n); Q = zeros(m,n); T = zeros(n,n);
Q(:,1) = b/norm(b);
for i = 1:n
    v = A*Q(:,i);
    T(i,i) = Q(:,i)'*v;
    v = v-(T(i,i)*Q(:,i));
    if i > 1
        v = v-(beta(i-1)*Q(:,i-1));
    end
    if i < n
        beta(i) = norm(v);
        T(i,i+1) = beta(i);
        T(i+1,i) = beta(i);
        Q(:,i+1) = v/beta(i);
    end
end