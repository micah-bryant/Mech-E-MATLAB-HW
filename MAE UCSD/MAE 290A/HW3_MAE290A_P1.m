m = 50; n = 12;
t = linspace(0,1,m); x = zeros(n,5);
A = fliplr(vander(t));
A = A(:,1:n);
tempA = A'*A;
b = cos(4.*t);
b = b';

%cholesky
R = chol(tempA);
y = forward(R',A'*b);
x(:,1) = backward(R,y);

%Modified Gram-Schmidt
[Q,R] = mgs(A);
x(:,2) = R\(Q'*b);

%built in QR decomp
[Q,R] = qr(A);
x(:,3) = R\(Q'*b);

%Some type of decomp idk what this one is
x(:,4) = A\(b);

%SVD
[U,S,V] = svd(A);
x(:,5) = V*(S\(U'*b));