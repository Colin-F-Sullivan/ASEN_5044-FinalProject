function [dstatedt] = dnonline_dt(t,state)
global mu
%Define Inputs
x = state(1);
xdot = state(2);
y = state(3);
ydot = state(4);

%Useful
r = sqrt(x^2 + y^2);

%Outputs/Derivatives
dstatedt(1) = xdot;
dstatedt(2) = (-mu*x)/r^3;
dstatedt(3) = ydot;
dstatedt(4) = (-mu*y)/r^3;

dstatedt = dstatedt';
end

