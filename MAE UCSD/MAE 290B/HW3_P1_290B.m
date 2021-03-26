%% Part c
clc; clear all;
%initialize variables
a = 1e-3; N = 81;
x = linspace(0,1,N); dt = 0.001;
a2 = -1/12; a1 = 4/3; a0 = -5/2; %a-2=a2 and a-1=a1 for 2nd derivative central FDA
dx1 = x(2)-x(1);
dx = dx1^2;
Tn = sin(8.*pi.*x); %to store Tn values
T = zeros(3,N); %to save desired T values at t = 1,5,25
Tn1 = zeros(1,N); %to store Tn+1 values

for i = 1:(25/dt)%cycles through time
    for j = 3:N-2 %solves inner nodes 
        num = ((a2*Tn(j-2)) + (a1*Tn(j-1)) + (a0*Tn(j)) + (a1*Tn(j+1)) + (a2*Tn(j+2)));
        Tn1(j) = Tn(j) + (a*dt*(num))/dx;
    end
    Tn1(1) = Tn(1) + (a*dt*(((a2*Tn(N-2)) + (a1*Tn(N-1)) + (a0*Tn(1)) + (a1*Tn(2)) + (a2*Tn(3))))/dx);
    Tn1(2) = Tn(2) + (a*dt*(((a2*Tn(N-1)) + (a1*Tn(1)) + (a0*Tn(2)) + (a1*Tn(3)) + (a2*Tn(4))))/dx);
    Tn1(N-1) = Tn(N-1) + (a*dt*(((a2*Tn(N-3)) + (a1*Tn(N-2)) + (a0*Tn(N-1)) + (a1*Tn(N)) + (a2*Tn(2))))/dx);
    Tn1(N) = Tn1(1);
    Tn = Tn1;
    if i == (1/dt)
        T(1,:) = Tn;
    end
    if i == (5/dt)
        T(2,:) = Tn;
    end
    if i == (25/dt)
        T(3,:) = Tn;
    end
end

plot(x,T)
xlabel("Distance (x)")
ylabel("Temperature (T)")
legend("T = 1", "T = 5", "T = 25")
%% Part D
clc; clear all;
format long;
%initialize variables
a = 1e-3; 
N = flip([40 50 80 100 200 400 500 800 1000]);
dx1 = 1./N;
N = N+1; %add one to include 0 point
dt = 0.000001;
a2 = -1/12; a1 = 4/3; a0 = -5/2; %a-2=a2 and a-1=a1 for 2nd derivative central FDA
dx = dx1.^2;
T = zeros(1,length(N));
for k = 1:length(N)
    disp(N(k))
    x = linspace(0,1,N(k));
    Tn = sin(8.*pi.*x); %to store Tn values
    Tn1 = zeros(1,N(k)); %to store Tn+1 values
    for i = 1:(1/dt)%cycles through time
        
        for j = 3:N(k)-2 %solves inner nodes
            num = ((a2*Tn(j-2)) + (a1*Tn(j-1)) + (a0*Tn(j)) + (a1*Tn(j+1)) + (a2*Tn(j+2)));
            Tn1(j) = Tn(j) + (a*dt*(num))/dx(k);
        end
        
        Tn1(1) = Tn(1) + (a*dt*(((a2*Tn(N(k)-2)) + (a1*Tn(N(k)-1)) + (a0*Tn(1)) + (a1*Tn(2)) + (a2*Tn(3))))/dx(k));
        
        Tn1(2) = Tn(2) + (a*dt*(((a2*Tn(N(k)-1)) + (a1*Tn(1)) + (a0*Tn(2)) + (a1*Tn(3)) + (a2*Tn(4))))/dx(k));
        
        Tn1(N(k)-1) = Tn(N(k)-1) + (a*dt*(((a2*Tn(N(k)-3)) + (a1*Tn(N(k)-2)) + (a0*Tn(N(k)-1)) + (a1*Tn(N(k))) + (a2*Tn(2))))/dx(k));
        
        Tn1(N(k)) = Tn1(1);
        
        Tn = Tn1;
        if i == (1/dt)
            idx = find(x==0.4);
            T(k) = Tn(idx);
        end
    end
end

err = abs(T(2:end)-T(1));
loglog(dx1(2:end),err); %visually checking to ensure it is linear
xlabel("dx")
ylabel("Error")
title("Log Log plot of Error as a Function of dx")
disp((log(err(end))-log(err(1)))/(log(dx1(end))-log(dx1(2)))); %determining slope of line
%% Part E
clc; clear all;
%initialize variables
a = 1e-3;
N = 80;
dx1 = 1/N;
N = N+1; %add one to include 0 point
dt = [0.0001 0.0005 0.001 0.005 0.01  0.05 0.1 0.5];
a2 = -1/12; a1 = 4/3; a0 = -5/2; %a-2=a2 and a-1=a1 for 2nd derivative central FDA
dx = dx1^2;
x = linspace(0,1,N);

T = zeros(1,length(dt));
for k = 1:length(dt)
    Tn = sin(8.*pi.*x); %to store Tn values
    Tn1 = zeros(1,N); %to store Tn+1 values
    for i = 1:(1/dt(k))%cycles through time
        for j = 3:N-2 %solves inner nodes
            num = ((a2*Tn(j-2)) + (a1*Tn(j-1)) + (a0*Tn(j)) + (a1*Tn(j+1)) + (a2*Tn(j+2)));
            Tn1(j) = Tn(j) + (a*dt(k)*(num))/dx;
        end
        Tn1(1) = Tn(1) + (a*dt(k)*(((a2*Tn(N-2)) + (a1*Tn(N-1)) + (a0*Tn(1)) + (a1*Tn(2)) + (a2*Tn(3))))/dx);
        Tn1(2) = Tn(2) + (a*dt(k)*(((a2*Tn(N-1)) + (a1*Tn(1)) + (a0*Tn(2)) + (a1*Tn(3)) + (a2*Tn(4))))/dx);
        Tn1(N-1) = Tn(N-1) + (a*dt(k)*(((a2*Tn(N-3)) + (a1*Tn(N-2)) + (a0*Tn(N-1)) + (a1*Tn(N)) + (a2*Tn(2))))/dx);
        Tn1(N) = Tn1(1);
        Tn = Tn1;
        %disp(Tn(0.5/dx1))
        if i == (1/dt(k))
            T(k) = Tn(0.5/dx1);
        end
    end
end

err = abs(T(2:end)-T(1));
loglog(dt(2:end),err); %visually checking to ensure it is linear
xlabel("dt")
ylabel("Error")
title("Log Log Plot of Error as a Function of dt")
disp((log(err(6))-log(err(1)))/(log(dt(6))-log(dt(2)))); %determining slope of line