function q = Q(x,y,t,a,omega)
f = 1 - (exp(-a*sqrt(t)) * sin(omega*t) * cos(2*omega*t)); %scalar
q = 2.5 .* sin(8.*pi.*x) .* sin(2.*pi.*y) .* f; %vector