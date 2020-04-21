%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: LKF_Main
% Author: Colin Sullivan
% 
% Date Created: 4/16/20
% Date Last Modified: 4/16/20
%
% Purpose:  Perform Linear KF on noisy y data generated from the non-linear
%           system of dnonline_dt.m
%
% Inputs:   None
% Outputs:  None yet cause it dont quite work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;

%% Set Up and Create Noisy Data
global mu r0
mu = 398600; %[km^3/s^2]
r0 = 6678;

%Time Vectors
dT = 10;

input = load('orbitdeterm_finalproj_KFdata.mat');
Q = input.Qtrue;
Rtrue = input.Rtrue;
pert = [0; 0; 1; 0];

[state,state_offnom,ynom,y_offnom,t] = ODE45_Progress1(1,pert);

y = y_offnom-ynom;

%Add Noise
for i = 1:length(y)
    %Find which stations are active
    station = FindStations(y,i);
    %For each of these active stations add the noise we need
    for j = 1:length(station)
    y(station(j)*3-2:station(j)*3,i) = y(station(j)*3-2:station(j)*3,i) + ...
                             mvnrnd([0; 0; 0],Rtrue)';
    end
end

%% KF
P_0 = .01*eye(4);
mu_0 = (state_offnom(:,1)-state(:,1));

%Preallocation
P_plus = zeros(length(P_0),length(P_0),length(t));
P_plus(:,:,1) = P_0;
mu_plus = zeros(length(mu_0),length(t));
mu_plus(:,1) = mu_0;
sigma_plus = zeros(length(mu_0),length(t));
sigma_plus(:,1) = 2.*sqrt(diag(abs(P_plus(:,:,1))));

%Kalman Filter Loop
for k = 2:length(t)
    
    %Make F
    [~,~,~,H_k,~,ObservingStations]...
    = ODEulerDTJacobians(state(:,k),dT,t(k));
    
    %Make y_k, R_k, H_k
    y_k = [];
    R_k = [];
    %Find the stations being used 
    for s = 1:length(ObservingStations)
        y_k_temp = y(ObservingStations(s)*3-2:ObservingStations(s)*3,k);
        y_k = [y_k; y_k_temp];
        R_k = blkdiag(R_k,Rtrue);
    end
    
    %Make F_k_mi1, G_k_minus_1
    [F_k_minus_1,G_k_minus_1,Om_k_minus_1,~,~,ObservingStations]...
    = ODEulerDTJacobians(state(:,k),dT,t(k-1));
    
    %Iterate Loop
    P_k_minus = F_k_minus_1*P_plus(:,:,k-1)*F_k_minus_1' + ...
        Om_k_minus_1*Q*Om_k_minus_1';
    K_k = P_k_minus*H_k'*((H_k*P_k_minus*H_k' + R_k)^-1);
    
    mu_minus = F_k_minus_1*mu_plus(:,k-1);% + G*u(:,k-1);
    mu_plus(:,k) = mu_minus + K_k*(y_k - H_k*mu_minus);
    
    P_plus(:,:,k) = (eye(size(K_k*H_k)) - K_k*H_k)*P_k_minus;
    sigma_plus(:,k) = 2.*sqrt(diag(abs(P_plus(:,:,k))));
    
end

figure;
plot(t,state_offnom(1,:)-state(1,:),'LineWidth',1.2);
hold on;
plot(t,mu_plus(1,:),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',12);
ylabel('\delta X-Position [m]','FontSize',12);
lgd = legend('\delta x_1','KF');
xlim([0 2*pi*sqrt(r0^3/mu)]);
title('Estimation of \delta x_1 v. Time','FontSize',14);

figure;
plot(t,state_offnom(2,:)-state(2,:),'LineWidth',1.2);
hold on;
plot(t,mu_plus(2,:),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',12);
ylabel('\delta X-Velocity [km/s]','FontSize',12);
lgd = legend('\delta x_2','KF');
xlim([0 2*pi*sqrt(r0^3/mu)]);
title('Estimation of \delta x_2 v. Time','FontSize',14);

figure;
plot(t,state_offnom(3,:)-state(3,:),'LineWidth',1.2);
hold on;
plot(t,mu_plus(3,:),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',12);
ylabel('\delta Y-Position [m]','FontSize',12);
lgd = legend('\delta x_3','KF');
xlim([0 2*pi*sqrt(r0^3/mu)]);
title('Estimation of \delta x_3 v. Time','FontSize',14);

figure;
plot(t,state_offnom(4,:)-state(4,:),'LineWidth',1.2);
hold on;
plot(t,mu_plus(4,:),'LineWidth',1.2);
grid on;
xlabel('Time [s]','FontSize',12);
ylabel('\delta Y-Velocity [km/s]','FontSize',12);
lgd = legend('\delta x_4','KF');
xlim([0 2*pi*sqrt(r0^3/mu)]);
title('Estimation of \delta x_4 v. Time','FontSize',14);


