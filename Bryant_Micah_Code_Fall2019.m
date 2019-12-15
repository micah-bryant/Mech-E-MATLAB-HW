clc
clear all;
%initialize all variables
syms theta3c theta4c theta5c theta6c
distance = 0;

%preallocates storage for data points
hDisplace = zeros(1,128);%distance table travels
vDisplace = zeros(1,128);%vertical displacement of table
aDisplace = zeros(1,128);%angular offset of table
aInput = zeros(1,128);%input angle

%vector lengths to decide (in m)
r1 = 0.805;

%mathematical relationships between vector lengths
r3 = r1/2;
r3a = r3/2;
r2 = 1.25*r1;
r4 = r2;
%experimentally determined ratios
r4a = r1*0.6295;
r5 = r1*0.6205;
r6 = r1*0.495;
rLeg = r1*0.2561;
cLen = ((rLeg*r6)+(rLeg*rLeg))/((2*rLeg)+r6);

%initial thetas starting from either thetaMax or thetaMin
thetaMin = 36.87*((2*pi)/360);
thetaMax = 100.51*((2*pi)/360);
theta = [180 36.87 112.38 99.49 38.32 -0.37];
theta = theta*((2*pi)/360);
tLeg = theta(6)+(pi/2);


%number of iterations to run Newton's method
maxiter = 100;
count = 1;
while distance < 0.85
    %initial conditions for theta3/4/5/6
    currVal1 = [theta(3);theta(4)];
    currVal2 = [theta(5);theta(6)];
    
    %functions for first vector loop
    f1 = -r1+(r2*cos(theta(2))+(r3*cos(theta3c))-(r4*cos(theta4c)));
    f2 = (r2*sin(theta(2))+(r3*sin(theta3c))-(r4*sin(theta4c)));
    fl1 = [f1;f2];
    
    %initializing jacobian
    v1 = [theta3c; theta4c];
    jacf1 = jacobian(fl1,v1);

    %creates functions to call in while loop
    func1 = inline(fl1,'theta3c','theta4c');
    jfunc1 = inline(jacf1,'theta3c','theta4c');
    error1 = norm(func1(currVal1(1),currVal1(2)),2);

    iter = 0;
    while error1 >= 0.001
        fxo=func1(currVal1(1),currVal1(2));
        fpxo=jfunc1(currVal1(1),currVal1(2));
        newVal1=currVal1-fpxo\fxo;
        fnewVal=func1(newVal1(1),newVal1(2));
        error1=norm((fnewVal),2);
        if iter > maxiter
            fprintf('exceeded maximum iterations loop 1 \n')
            break
        end
        currVal1=newVal1;
        iter=iter+1;
    end

    %extracting theta3/4 from first vector loop
    theta(3) = currVal1(1);
    theta(4) = currVal1(2);
    %functions for second vector loop
    f3 = -r1+(r2*cos(theta(2)))+(r3a*cos(theta(3)))+(r6*cos(theta6c))-(r5*cos(theta5c))-(r4a*cos(theta(4)));
    f4 = (r2*sin(theta(2)))+(r3a*sin(theta(3)))+(r6*sin(theta6c))-(r5*sin(theta5c))-(r4a*sin(theta(4)));
    fl2 = [f3;f4];
    
    %initializing jacobian
    v2 = [theta5c; theta6c];
    jacf2 = jacobian(fl2,v2);

    %initializing functions to call in while loop
    func2 = inline(fl2,'theta5c','theta6c');
    jfunc2 = inline(jacf2,'theta5c','theta6c');
    error2 = norm(func2(currVal2(1),currVal2(2)),2);

    iter = 0;
    while error2 >= 0.001
        fxo=func2(currVal2(1),currVal2(2));
        fpxo=jfunc2(currVal2(1),currVal2(2));
        newVal2=currVal2-fpxo\fxo;
        fnewVal=func2(newVal2(1),newVal2(2));
        error2=norm((fnewVal),2);
        if iter > maxiter
            fprintf('exceeded maximum iterations loop 2 \n')
            break
        end
        currVal2=newVal2;
        iter=iter+1;
    end
    
    %extracting theta5/6 from second vector loop
    theta(5) = currVal2(1);
    theta(6) = currVal2(2);
    tLeg = theta(6)+(pi/2);
    
    if count == 1
        initLocX = r1+(r4a*cos(theta(4))+(r5*cos(theta(5)))+(cLen*cos(tLeg))-(0.5*r6*cos(theta(6))));
        initLocY = (r4a*sin(theta(4))+(r5*sin(theta(5)))+(cLen*sin(tLeg))-(0.5*r6*sin(theta(6))));
    end

    %solve for new distance, add on additional input displacement
    newLocX = r1+(r4a*cos(theta(4))+(r5*cos(theta(5)))+(cLen*cos(tLeg))-(0.5*r6*cos(theta(6))));
    newLocY = (r4a*sin(theta(4))+(r5*sin(theta(5)))+(cLen*sin(tLeg))-(0.5*r6*sin(theta(6))));

    distance = initLocX-newLocX;
    hDisplace(count) = distance;
    vDisplace(count) = initLocY-newLocY;
    aDisplace(count) = (theta(6)/(2*pi))*360;
    aInput(count) = (theta(2)/(2*pi))*360;
    
    disp(count)
    
    %check if input displacement is out of bounds
    theta(2) = theta(2)+(0.008727);%adds on equivalent of 0.5 degree
    count = count+1;
    if theta(2) > thetaMax || theta(2) < thetaMin
        fprintf('Next input angle out of bounds \n')
        break
    end
end
%plot displacements vs input angle (see BE 182 book for example)

distance = 0;
rho = 2689;
t = 0.03;
t2 = t*t;

%preallocates storage for data points
torque = zeros(1,128);%torque of mechanism (may need to do this separately)
timer = zeros(1,128);

%initial thetas starting from either thetaMax or thetaMin
torqueMin = 36.87*((2*pi)/360);
torqueMax = 97.374*((2*pi)/360);
theta = [180 36.87 112.38 99.49 38.32 -0.37];
theta = theta*((2*pi)/360);
tLeg = theta(6)+(pi/2);
primes = [0 1 1 1 1 1];
dPrimes = [0 1 1 1 1 1];

%initialize masses and inertias
mass = [r1 r2 r3 r4 r5 r6 rLeg rLeg];
mass = mass*(rho*t2);
inertia = [(r1*r1) (r2*r2) (r3*r3) (r4*r4) (r5*r5) (r6*r6) (rLeg*rLeg) (rLeg*rLeg)]+t2;
inertia = (mass.*inertia)/12;

%placeholders for Newton's method
newVal1 = [0.;0.];
newVal2 = [0.;0.];

%number of iterations to run Newton's method
maxiter = 100;
countT = 1;

time = 0;
step = (6*pi)/72;
while time < 6*pi
    theta(2) = (((torqueMax-torqueMin)/2)*sin(time))+((torqueMax+torqueMin)/2);%moves it according to constraint
    inVel = ((torqueMax-torqueMin)/2)*cos(time);
    inAcc = -((torqueMax-torqueMin)/2)*sin(time);
    %initial conditions for theta3/4/5/6
    currVal1 = [theta(3);theta(4)];
    currVal2 = [theta(5);theta(6)];
    
    %functions for first vector loop
    f1 = -r1+(r2*cos(theta(2))+(r3*cos(theta3c))-(r4*cos(theta4c)));
    f2 = (r2*sin(theta(2))+(r3*sin(theta3c))-(r4*sin(theta4c)));
    fl1 = [f1;f2];
    
    %initializing jacobian
    v1 = [theta3c; theta4c];
    jacf1 = jacobian(fl1,v1);

    %creates functions to call in while loop
    func1 = inline(fl1,'theta3c','theta4c');
    jfunc1 = inline(jacf1,'theta3c','theta4c');
    error1 = norm(func1(currVal1(1),currVal1(2)),2);

    iter = 0;
    while error1 >= 0.001
        fxo=func1(currVal1(1),currVal1(2));
        fpxo=jfunc1(currVal1(1),currVal1(2));
        newVal1=currVal1-fpxo\fxo;
        fnewVal=func1(newVal1(1),newVal1(2));
        error1=norm((fnewVal),2);
        if iter > maxiter
            fprintf('exceeded maximum iterations loop 1 \n')
            break
        end
        currVal1=newVal1;
        iter=iter+1;
    end

    %extracting theta3/4 from first vector loop
    theta(3) = currVal1(1);
    theta(4) = currVal1(2);
    %functions for second vector loop
    f3 = -r1+(r2*cos(theta(2)))+(r3a*cos(theta(3)))+(r6*cos(theta6c))-(r5*cos(theta5c))-(r4a*cos(theta(4)));
    f4 = (r2*sin(theta(2)))+(r3a*sin(theta(3)))+(r6*sin(theta6c))-(r5*sin(theta5c))-(r4a*sin(theta(4)));
    fl2 = [f3;f4];
    
    %initializing jacobian
    v2 = [theta5c; theta6c];
    jacf2 = jacobian(fl2,v2);

    %initializing functions to call in while loop
    func2 = inline(fl2,'theta5c','theta6c');
    jfunc2 = inline(jacf2,'theta5c','theta6c');
    error2 = norm(func2(currVal2(1),currVal2(2)),2);

    iter = 0;
    while error2 >= 0.001
        fxo=func2(currVal2(1),currVal2(2));
        fpxo=jfunc2(currVal2(1),currVal2(2));
        newVal2=currVal2-fpxo\fxo;
        fnewVal=func2(newVal2(1),newVal2(2));
        error2=norm((fnewVal),2);
        if iter > maxiter
            fprintf('exceeded maximum iterations loop 2 \n')
            break
        end
        currVal2=newVal2;
        iter=iter+1;
    end
    
    %extracting theta5/6 from second vector loop
    theta(5) = currVal2(1);
    theta(6) = currVal2(2);
    tLeg = theta(6)+(pi/2);
    
    if countT == 1
        initLocX = r1+(r4a*cos(theta(4))+(r5*cos(theta(5)))+(cLen*cos(tLeg))-(0.5*r6*cos(theta(6))));
        initLocY = (r4a*sin(theta(4))+(r5*sin(theta(5)))+(cLen*sin(tLeg))-(0.5*r6*sin(theta(6))));
    end

    %solve for new distance, add on additional input displacement
    jMat1 = [(-r3*sin(theta(3))) (r4*sin(theta(4)));
        (r3*cos(theta(3))) (-r4*cos(theta(4)))];
    sol1 = [r2*sin(theta(2));
        -r2*cos(theta(2))];
    prime1 = jMat1\sol1;
    
    primes(3) = prime1(1);
    primes(4) = prime1(2);
    
    dSol1 = [(r2*cos(theta(2)))+(r3*primes(3)*primes(3)*cos(theta(3)))-(r4*primes(4)*primes(4)*cos(theta(4)));
        (r2*sin(theta(2)))+(r3*primes(3)*primes(3)*sin(theta(3)))-(r4*primes(4)*primes(4)*sin(theta(4)))];
    dPrime1 = jMat1\dSol1;
    dPrimes(3) = dPrime1(1);
    dPrimes(4) = dPrime1(2);
    
    jMat2 = [(r5*sin(theta(5))) (-r6*sin(theta(6)));
        (-r5*cos(theta(5))) (r6*cos(theta(6)))];
    sol2 = [(r2*sin(theta(2)))+(r3a*primes(3)*sin(theta(3)))-(r4a*primes(4)*sin(theta(4)));
        -(r2*cos(theta(2)))-(r3a*primes(3)*cos(theta(3)))+(r4a*primes(4)*cos(theta(4)))];
    prime2 = jMat2\sol2;
    primes(5) = prime2(1);
    primes(6) = prime2(2);
    
    dSol2 = [(r2*cos(theta(2)))+(r3a*primes(3)*primes(3)*cos(theta(3)))-(r4a*primes(4)*primes(4)*cos(theta(4)))-(r5*primes(5)*primes(5)*cos(theta(5)))+(r6*primes(6)*primes(6)*cos(theta(6)))+(r3a*dPrimes(3)*sin(theta(3)))-(r4a*dPrimes(4)*sin(theta(4)));
        (r2*sin(theta(2)))+(r3a*primes(3)*primes(3)*sin(theta(3)))-(r4a*primes(4)*primes(4)*sin(theta(4)))-(r5*primes(5)*primes(5)*sin(theta(5)))+(r6*primes(6)*primes(6)*sin(theta(6)))-(r3a*dPrimes(3)*cos(theta(3)))+(r4a*dPrimes(4)*cos(theta(4)))];
    dPrime2 = jMat2\dSol2;
    dPrimes(5) = dPrime2(1);
    dPrimes(6) = dPrime2(2);
    
    
    %positional coefficient calculations
    xPrime = [0
        (-0.5*r2*sin(theta(2)))
        (-r2*sin(theta(2)))-(0.5*r3*primes(3)*sin(theta(3)))
        (-r2*sin(theta(2)))-(r3*primes(3)*sin(theta(3)))+(0.5*r4*primes(4)*sin(theta(4)))
        (-r4a*primes(4)*sin(theta(4)))-(0.5*r5*primes(5)*sin(theta(5)))
        (-r4a*primes(4)*sin(theta(4)))-(r5*primes(5)*sin(theta(5)))-(rLeg*primes(6)*sin(tLeg))+(0.5*r6*primes(6)*sin(theta(6)))
        (-r4a*primes(4)*sin(theta(4)))-(r5*primes(5)*sin(theta(5)))-(rLeg*primes(6)*sin(tLeg))+(r6*primes(6)*sin(theta(6)))+(0.5*rLeg*primes(6)*sin(tLeg))%leg1
        (-r4a*primes(4)*sin(theta(4)))-(r5*primes(5)*sin(theta(5)))-(0.5*rLeg*primes(6)*sin(tLeg))];%leg2
    yPrime = [0
        (0.5*r2*cos(theta(2)))
        (r2*cos(theta(2)))+(0.5*r3*primes(3)*cos(theta(3)))
        (r2*cos(theta(2)))+(r3*primes(3)*cos(theta(3)))-(0.5*r4*primes(4)*cos(theta(4)))
        (r4a*primes(4)*cos(theta(4)))+(0.5*r5*primes(5)*cos(theta(5)))
        (r4a*primes(4)*cos(theta(4)))+(r5*primes(5)*cos(theta(5)))+(rLeg*primes(6)*cos(tLeg))-(0.5*r6*primes(6)*cos(theta(6)))
        (r4a*primes(4)*cos(theta(4)))+(r5*primes(5)*cos(theta(5)))+(rLeg*primes(6)*cos(tLeg))-(r6*primes(6)*cos(theta(6)))-(0.5*rLeg*primes(6)*cos(tLeg))%leg1
        (r4a*primes(4)*cos(theta(4)))+(r5*primes(5)*cos(theta(5)))+(0.5*rLeg*primes(6)*cos(tLeg))];%leg2
    xDPrime = [0
        (-0.5*r2*cos(theta(2)))
        (-r2*cos(theta(2)))-(0.5*r3*primes(3)*primes(3)*cos(theta(3)))-(0.5*r3*dPrimes(3)*sin(theta(3)))
        (-r2*cos(theta(2)))-(r3*primes(3)*primes(3)*cos(theta(3)))-(r3*dPrimes(3)*sin(theta(3)))+(0.5*r4*primes(4)*primes(4)*cos(theta(4)))+(0.5*r4*dPrimes(4)*sin(theta(4)))
        (-r4a*primes(4)*primes(4)*cos(theta(4)))-(r4a*dPrimes(4)*sin(theta(4)))-(0.5*r5*primes(5)*primes(5)*cos(theta(5)))-(0.5*r5*dPrimes(5)*sin(theta(5)))
        (-r4a*primes(4)*primes(4)*cos(theta(4)))-(r4a*dPrimes(4)*sin(theta(4)))-(r5*primes(5)*primes(5)*cos(theta(5)))-(r5*dPrimes(5)*sin(theta(5)))-(rLeg*primes(6)*primes(6)*cos(tLeg))-(rLeg*dPrimes(6)*sin(tLeg))+(0.5*r6*primes(6)*primes(6)*cos(theta(6)))+(0.5*r6*dPrimes(6)*sin(theta(6)))
        (-r4a*primes(4)*primes(4)*cos(theta(4)))-(r4a*dPrimes(4)*sin(theta(4)))-(r5*primes(5)*primes(5)*cos(theta(5)))-(r5*dPrimes(5)*sin(theta(5)))-(rLeg*primes(6)*primes(6)*cos(tLeg))-(rLeg*dPrimes(6)*sin(tLeg))+(r6*primes(6)*primes(6)*cos(theta(6)))+(r6*dPrimes(6)*sin(theta(6)))+(0.5*rLeg*primes(6)*primes(6)*cos(tLeg))+(0.5*rLeg*dPrimes(6)*sin(tLeg))%leg1
        (-r4a*primes(4)*primes(4)*cos(theta(4)))-(r4a*dPrimes(4)*sin(theta(4)))-(r5*primes(5)*primes(5)*cos(theta(5)))-(r5*dPrimes(5)*sin(theta(5)))-(0.5*rLeg*primes(6)*primes(6)*cos(tLeg))-(0.5*rLeg*dPrimes(6)*sin(tLeg))];%leg2
    yDPrime = [0
        (-0.5*r2*sin(theta(2)))
        (-r2*sin(theta(2)))-(0.5*r3*primes(3)*primes(3)*sin(theta(3)))+(0.5*r3*dPrimes(3)*cos(theta(3)))
        (-r2*sin(theta(2)))-(r3*primes(3)*primes(3)*sin(theta(3)))+(r3*dPrimes(3)*cos(theta(3)))+(0.5*r4*primes(4)*primes(4)*sin(theta(4)))-(0.5*r4*dPrimes(4)*cos(theta(4)))
        (-r4a*primes(4)*primes(4)*sin(theta(4)))+(r4a*dPrimes(4)*cos(theta(4)))-(0.5*r5*primes(5)*primes(5)*sin(theta(5)))+(0.5*r5*dPrimes(5)*cos(theta(5)))
        (-r4a*primes(4)*primes(4)*sin(theta(4)))+(r4a*dPrimes(4)*cos(theta(4)))-(r5*primes(5)*primes(5)*sin(theta(5)))+(r5*dPrimes(5)*cos(theta(5)))-(rLeg*primes(6)*primes(6)*sin(tLeg))+(rLeg*dPrimes(6)*cos(tLeg))+(0.5*r6*primes(6)*primes(6)*sin(theta(6)))-(0.5*r6*dPrimes(6)*cos(theta(6)))
        (-r4a*primes(4)*primes(4)*sin(theta(4)))+(r4a*dPrimes(4)*cos(theta(4)))-(r5*primes(5)*primes(5)*sin(theta(5)))+(r5*dPrimes(5)*cos(theta(5)))-(rLeg*primes(6)*primes(6)*sin(tLeg))+(rLeg*dPrimes(6)*cos(tLeg))+(r6*primes(6)*primes(6)*sin(theta(6)))-(r6*dPrimes(6)*cos(theta(6)))+(0.5*rLeg*primes(6)*primes(6)*sin(tLeg))-(0.5*rLeg*dPrimes(6)*cos(tLeg))%leg1
        (-r4a*primes(4)*primes(4)*sin(theta(4)))+(r4a*dPrimes(4)*cos(theta(4)))-(r5*primes(5)*primes(5)*sin(theta(5)))+(r5*dPrimes(5)*cos(theta(5)))-(0.5*rLeg*primes(6)*primes(6)*sin(tLeg))+(0.5*rLeg*dPrimes(6)*cos(tLeg))];%leg2
    
    %dT/dt Calculations
    A1 = [0
        (xPrime(2)^2)+(yPrime(2)^2)
        (xPrime(3)^2)+(yPrime(3)^2)
        (xPrime(4)^2)+(yPrime(4)^2)
        (xPrime(5)^2)+(yPrime(5)^2)
        (xPrime(6)^2)+(yPrime(6)^2)
        (xPrime(7)^2)+(yPrime(7)^2)
        (xPrime(8)^2)+(yPrime(8)^2)];
    A1 = A1.*mass';
    A2 = [0
        (primes(2)^2)
        (primes(3)^2)
        (primes(4)^2)
        (primes(5)^2)
        (primes(6)^2)
        (primes(6)^2)
        (primes(6)^2)];
    A2 = A2.*inertia';
    A = sum(A1+A2);
    
    B1 = [0
        (xPrime(2)*xDPrime(2))+(yPrime(2)*yDPrime(2))
        (xPrime(3)*xDPrime(3))+(yPrime(3)*yDPrime(3))
        (xPrime(4)*xDPrime(4))+(yPrime(4)*yDPrime(4))
        (xPrime(5)*xDPrime(5))+(yPrime(5)*yDPrime(5))
        (xPrime(6)*xDPrime(6))+(yPrime(6)*yDPrime(6))
        (xPrime(7)*xDPrime(7))+(yPrime(7)*yDPrime(7))
        (xPrime(8)*xDPrime(8))+(yPrime(8)*yDPrime(8))];
    B1 = B1.*mass';
    B2 = [0;
        (primes(2)*dPrimes(2));
        (primes(3)*dPrimes(3));
        (primes(4)*dPrimes(4));
        (primes(5)*dPrimes(5));
        (primes(6)*dPrimes(6));
        (primes(6)*dPrimes(6));
        (primes(6)*dPrimes(6))];
    B2 = B2.*inertia';
    B = sum(B1+B2);
    
    dT = (A*inVel*inAcc)+(B*(inVel^3));
    
    %dU/dt Calculations
    dU = 9.80665*inVel*sum(yPrime.*mass');
    
    %Moment Calculation
    moment = (dT+dU)/inVel;
    
    torque(countT) = moment;
    
    newLocX = r1+(r4a*cos(theta(4))+(r5*cos(theta(5)))+(cLen*cos(tLeg))-(0.5*r6*cos(theta(6))));
    newLocY = (r4a*sin(theta(4))+(r5*sin(theta(5)))+(cLen*sin(tLeg))-(0.5*r6*sin(theta(6))));

    distance = newLocX-initLocX;
    timer(countT) = time;
    disp(countT)
    
    countT = countT+1;
    time = time+step;
end
tiledlayout(2,2)

ax1 = nexttile;
h1 = plot(ax1, aInput(1:count-1), hDisplace(1:count-1)*1000);
ylabel(ax1,'Horizontal Displacement of COM (mm)')
xlabel(ax1, 'Input Angle (Degrees)')

ax2 = nexttile;
h2 = plot(ax2, aInput(1:count-1), vDisplace(1:count-1)*1000);
ylabel(ax2,'Vertical Displacement of COM (mm)')
xlabel(ax2, 'Input Angle (Degrees)')

ax3 = nexttile;
h3 = plot(ax3, aInput(1:count-1), aDisplace(1:count-1));
ylabel(ax3,'Rotation of the Table (Degrees)')
xlabel(ax3, 'Input Angle (Degrees)')

ax4 = nexttile;
h4 = plot(ax4, timer(1:countT-1), torque(1:countT-1));
ylabel(ax4,'Input Torque (N-m)')
xlabel(ax4, 'Time')