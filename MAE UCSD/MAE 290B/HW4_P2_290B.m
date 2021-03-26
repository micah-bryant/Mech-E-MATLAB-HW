%% Part B
clc; clear all;
M = 33; N = 33;
phi = ones(M,N); %initial guess is phi(0) = 1
phi(1,:) = 0; phi(M,:) = 0; phi(:,1) = 0; phi(:,N) = 0;
phiK = phi;
x = linspace(0,1,M); y = linspace(0,1,N);
error = zeros(100,1);
for k = 1:100 %completes 100 iterations
    for i = 2:M-1 %cycle along x coordinates
        for j = 2:N-1 %cycle along y coordinates
            phiK(i,j) = (1/4)*(phiK(i-1,j) + phiK(i,j-1) + phi(i+1,j) + phi(i,j+1));
        end
    end
    error(k) = max(max(abs(phi-phiK)));
    phi = phiK;
end
%% Plotting
[X,Y] = meshgrid(x);
surf(X,Y,phi)
xlabel('X Values');
ylabel('Y Values');
zlabel('Phi');
title('Laplace Solution Using Gauss-Seidel')
%% Error Graph
k = linspace(1,100,100);
plot(k,error)
xlabel("Iteration Number")
ylabel("Error")
title('Error using Difference Between each Iteration')