function dudt = odefcn3(u,j,dx)
%velocity using 2nd order backward difference FDA
dudt = (-1/2) * ((u(j-2)-(4*u(j-1))+(3*u(j))) / (2*dx));