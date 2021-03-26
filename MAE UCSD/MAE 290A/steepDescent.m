function [x,r] = steepDescent(A,b,n)
[m,~] = size(A);
x = zeros(m,1); r0 = b;
for i = 1:n
    a = (r0'*r0)/(r0'*A*r0); %defines step size for iteration
    if i > 1
        x = x+(a*r0);%calculates approx solution
    else
        x = a*r0; %calculates step with polynomial
    end
    r1 = r0-(a*A*r0); %determines new residual
    r0=r1; %sets new residual to be the old residual for next iteration
end
r=r1;