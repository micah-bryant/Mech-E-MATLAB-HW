function dTdx = odefcn(T,j,dx)
%Temperature using 2nd order central difference FDA
%x is arbitrary spatial variable
dTdx = ((T(j-1) - (2*T(j)) + T(j+1)) / (dx^2));