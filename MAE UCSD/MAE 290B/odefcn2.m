function dudt = odefcn2(u,j,dx)
%velocity using upwind First Order FDA
dudt = (-1/2) * ((u(j)-u(j-1)) / (dx));