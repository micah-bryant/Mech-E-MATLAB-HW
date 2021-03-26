%% Part 1
%Define Overarching equation
dxdt = @(t,x) (A*x)-(1/2*(H*(x.^2)));

%Define Variables
T = 3; nt = 500; r = 2*linspace(1,15,15); x0 = x0a;

tspan = linspace(0,T,nt);
eta = linspace(0,1,1024);

options = odeset('mass',E);
[t,x] = ode45(dxdt, tspan, x0, options);
[etaG, tG] = meshgrid(eta, tspan);

%set boundary conditions
x(:,1) = zeros(500,1);
x(:,1024) = zeros(500,1);

h = surf(tG, etaG, x);
set(h,'LineStyle','none')
%% Part 2
X = x'; %get orientation correct to match notes

%compute SVD of Correlation Matrix
[U,S,V] = svd(X'*E*X,'econ');
sValues = diag(S);
xAxis = linspace(1,length(sValues),length(sValues));
semilogy(xAxis, sValues)

%% Part 3
%determine basis of POD-ROM
time=zeros(1,15);
errors = zeros(1,15);
for j=1:15 %run for all r values
    rTemp = r(j);
    phiBasis = zeros(1024,rTemp);
    phiCorr = U(:,1:rTemp);
    for i = 1:rTemp
        phiBasis(:,i) = (1/sqrt(sValues(i))) .* (X*phiCorr(:,i));
    end
    Er=phiBasis'*E*phiBasis;
    opt = odeset("Mass",Er);
    %write equations 1.2 from Lecture Notes
    dxrdt = @(tROM,xROM) (phiBasis'*A*phiBasis*xROM)-(1/2*phiBasis'*H*((phiBasis*xROM).^2));
    xr0 = phiBasis'*E*x0;
    
    %Solve Differential Equation and Measure Time Elapsed
    tic
    [tROM,xROM] = ode45(dxrdt, tspan, xr0, opt);
    toc
    time(j)=toc;
    
    %convert xROM to correct orientation for calculations
    xROM = xROM';
    
    %Take Error
    Xrtilde = phiBasis*xROM;
    
    errors(j) = norm(X-Xrtilde)/norm(X);
end
