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
lgd = legend('X-Position','Y-Position');
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

save('state_nonlin_nom.mat','state');

%Do off-nominal condition
state0_offnom = [r0; 0; 100; r0.*sqrt(mu/r0^3)];
[t_offnom,state_offnom] = ode45('dnonline_dt',ti:dT:tf,state0_offnom,opts);
%Plot it
figure;
subplot(1,2,1);
plot(t_offnom,state_offnom(:,1)-state(:,1),'LineWidth',1.2);
hold on;
plot(t_offnom,state_offnom(:,3)-state(:,3),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('\delta Position [m]','FontSize',14);
title('Position Deltas','FontSize',14);
lgd = legend('\delta x_1','\delta x_3');
lgd.FontSize = 12;
xlim([0 P]);
subplot(1,2,2);
plot(t_offnom,state_offnom(:,2)-state(:,2),'LineWidth',1.2);
hold on;
plot(t_offnom,state_offnom(:,4)-state(:,4),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('\delta Position [m]','FontSize',14);
title('Velocity Deltas','FontSize',14);
lgd = legend('\delta x_2','\delta x_4');
lgd.FontSize = 12;
xlim([0 P]);

%% Do the same thing for each station
y = [];
for i = 1:12
y_i = Make_Y_i(i,t,state);
%Plot it
% figure;
% %Range
% subplot(1,3,1);
% plot(t,y_i(1,:),'LineWidth',1.2);
% hold on;
% grid on;
% xlabel('Time [s]','FontSize',14);
% ylabel('Relative Range [m]','FontSize',14);
% tit = ['Relative Range v. Time for Station ' num2str(i)];
% title(tit,'FontSize',14);
% xlim([0 max(t)]);
% %Range Rate
% subplot(1,3,2);
% plot(t,y_i(2,:),'LineWidth',1.2);
% hold on;
% grid on;
% xlabel('Time [s]','FontSize',14);
% ylabel('Relative Range Rate [km/s]','FontSize',14);
% tit = ['Relative Range Rate v. Time for Station ' num2str(i)];
% title(tit,'FontSize',14);
% xlim([0 max(t)]);
% %Elevation Angle
% subplot(1,3,3);
% plot(t,y_i(3,:),'LineWidth',1.2);
% hold on;
% grid on;
% xlabel('Time [s]','FontSize',14);
% ylabel('Elevation Angle [rad]','FontSize',14);
% tit = ['Elevation Angle v. Time for Station ' num2str(i)];
% title(tit,'FontSize',14);
% xlim([0 max(t)]);
y = [y; y_i];
end

%Save y nominal condition
save('y_nonlin_nom.mat','y');
%Make Off nom
y_offnom = [];
for i = 1:12
    y_i_offnom = Make_Y_i(i,t_offnom,state_offnom);
    y_offnom = [y_offnom; y_i_offnom];
end

%Plotting just the first one
figure;
%Range
subplot(1,3,1);
plot(t_offnom,y_offnom(1,:)-y(1,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('\delta Relative Range [m]','FontSize',14);
tit = ['\delta Relative Range v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);
%Range Rate
subplot(1,3,2);
plot(t_offnom,y_offnom(2,:)-y(2,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('\delta Relative Range Rate [km/s]','FontSize',14);
tit = ['\delta Relative Range Rate v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);
%Elevation Angle
subplot(1,3,3);
plot(t_offnom,y_offnom(3,:)-y(3,:),'LineWidth',1.2);
hold on;
grid on;
xlabel('Time [s]','FontSize',14);
ylabel('\delta Elevation Angle [rad]','FontSize',14);
tit = ['\delta Elevation Angle v. Time for Station ' num2str(i)];
title(tit,'FontSize',14);
xlim([0 max(t)]);