%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODE45_Progress1
% Author: Colin Sullivan
% 
% Date Created: 4/12/20
% Date Last Modified: 4/12/20
%
% Inputs: N/A
% Outputs: Plots of Non-Linearized X and Y state vectors for part c.) of
%          Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Housekeeping
clc;
close all;

%Constants
global mu r0
mu = 398600; %[km^3/s^2]
r0 = 6678;

%Time Vectors
P = 2*pi*sqrt(r0^3/mu);
ti = 0;
dT = 10;
tf = P; %Plot over one period of circular orbit

%Setup ODE45
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
state0 = [r0; 0; 0; r0.*sqrt(mu/r0^3)];
[t,state] = ode45('dnonline_dt',ti:dT:tf,state0,opts);

%Plot X and Y as a function of time;
figure;
subplot(1,2,1);
plot(ti:dT:tf,state(:,1),'LineWidth',1.2);
xlabel('Time [s]','FontSize',14);
ylabel('Position [m]','FontSize',14);
hold on;
grid on;
plot(ti:dT:tf,state(:,3),'Linewidth',1.2);
lgd = legend('X Position','Y-Position');
lgd.FontSize = 14;
xlim([0 P]);
title('Non-Linear Orbital Position Vectors v. Time','FontSize',14);
subplot(1,2,2);

%Plot Velocities as a Function of Time
plot(ti:dT:tf,state(:,2),'LineWidth',1.2);
hold on;
xlabel('Time [s]','FontSize',14);
ylabel('Velocity [km/s]','FontSize',14);
grid on;
plot(ti:dT:tf,state(:,4),'LineWidth',1.2);
lgd = legend('X-Velocity','Y-Velocity');
lgd.FontSize = 14;
xlim([0 P]);
title('Non-Linear Orbital Velocity Vectors v. Time','FontSize',14);

%% Do the same thing for each station
for i = 1:5
y_i = Make_Y_i(i,t,state);
%Plot it
figure;
%Range
subplot(1,3,1);
plot(t,y_i(1,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('Relative Range [m]','FontSize',14);
tit = ['Relative Range v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);
%Range Rate
subplot(1,3,2);
plot(t,y_i(2,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('Relative Range Rate [km/s]','FontSize',14);
tit = ['Relative Range Rate v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);
%Elevation Angle
subplot(1,3,3);
plot(t,y_i(3,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('Elevation Angle [rad]','FontSize',14);
tit = ['Elevation Angle v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);

end

