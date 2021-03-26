function [x,r] = CG(A,b,n)
[m,~] = size(A);
x = zeros(m,1); r0 = b; p = b;
for i = 1:n
    a = (r0'*r0)/(p'*A*p); %defines step size for iteration
    if i > 1
        x = x+(a*p);%calculates approx solution
    else
        x = a*p; %calculates step with polynomial
    end
    r1 = r0-(a*A*p); %determines new residual
    beta = (r1'*r1)/(r0'*r0); %determines thin step length
    p = r1+(beta*p); %determines new polynomial
    r0=r1; %sets new residual to be the old residual for next iteration
end
r=r1;