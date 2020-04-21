function [dstatedt] = statOD_dynamics(t,state,u,w)
global mu
%Define Inputs
x = state(1);
xdot = state(2);
y = state(3);
ydot = state(4);

u1=u(1);
u2=u(2);

w1=w(1);
w2=w(2);

%Useful
r = sqrt(x^2 + y^2);

%Outputs/Derivatives
dstatedt(1) = xdot;
dstatedt(2) = (-mu*x)/r^3+u1+w1;
dstatedt(3) = ydot;
dstatedt(4) = (-mu*y)/r^3+u2+w2;

dstatedt = dstatedt';
end

