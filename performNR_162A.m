clear all
syms theta3 theta4
currVal = [3.2;1.8];
f1 = -0.18+(0.06*cos(1.047))-(0.15*cos(theta3))-(0.08*cos(theta4));
f2 = (0.06*sin(1.047))-(0.15*sin(theta3))-(0.08*sin(theta4));
f = [f1; f2];
v = [theta3;theta4];
jacf = jacobian(f,v);
iter = 1;
maxiter = 30;
func = inline(f);
jfunc = inline(jacf);
error = norm(func(currVal(1),currVal(2)),2);
newVal=[0,0];
while error >= 0.0001
    fxo=func(currVal(1),currVal(2));
    fpxo=jfunc(currVal(1),currVal(2));
    newVal=currVal-inv(fpxo)*fxo;
    fnewVal=func(newVal(1),newVal(2));
    error=norm((fnewVal),2);
    if iter > maxiter
        fprintf('exceeded maximum iterations \n')
        break
    end
    currVal=newVal;
    iter=iter+1;
end
disp(iter)
disp(newVal)
