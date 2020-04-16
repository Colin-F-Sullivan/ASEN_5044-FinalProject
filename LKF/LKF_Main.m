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
clear all;
close all;

%% Set Up and Create Noisy Data
global mu r0
mu = 398600; %[km^3/s^2]
r0 = 6678;

%Time Vectors
P = 2*pi*sqrt(r0^3/mu);
ti = 0;
dT = 10;
tf = P+5*dT; %Plot over one period of circular orbit + 2 time step


input = load('orbitdeterm_finalproj_KFdata.mat');
Qtrue = input.Qtrue;
Rtrue = input.Rtrue;

[state,state_offnom,y,t] = ODE45_Progress1(1,[0; 0; 0; 0]);
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
P_0 = eye(4);
mu_0 = [0 0 0 0]';
Q = eye(4);
Q(2,2) = Qtrue(1,1);
Q(4,4) = Qtrue(2,2);


%Preallocation
P_plus = zeros(length(P_0),length(P_0),length(t));
P_plus(:,:,1) = P_0;
mu_plus = zeros(length(mu_0),length(t));
mu_plus(:,1) = mu_0;
sigma_plus = zeros(length(mu_0),length(t));
sigma_plus(:,1) = 2.*sqrt(diag(abs(P_plus(:,:,1))));

%Kalman Filter Loop
for k = 2:length(t)
    
    %Make F_{k-1}
    [F,G,Omtilde,H,M,ObservingStations]...
    = ODEulerDTJacobians(state,dT,t(k));
    
    %Make y_k, R_k, H_k
    y_k = [];
    R_k = [];
    for i = length(ObservingStations)
        %y_k
        y_k_temp = y(ObservingStations(i)*3-2:ObservingStations(i)*3,k);
        y_k = [y_k; y_k_temp];
        %R_k
        R = blkdiag(R_k,Rtrue);
    end
   
    
    %Iterate Loop
    P_k_minus = F*P_plus(:,:,k-1)*F' + Q;
    K_k = P_k_minus*H'*((H*P_k_minus*H' + R)^-1);
    
    mu_minus = F*mu_plus(:,k-1);% + G*u(:,k-1);
    mu_plus(:,k) = mu_minus + K_k*(y_k - H*mu_minus);
    
    P_plus(:,:,k) = (eye(size(K_k*H)) - K_k*H)*P_k_minus;
    sigma_plus(:,k) = 2.*sqrt(diag(abs(P_plus(:,:,k))));
    
end