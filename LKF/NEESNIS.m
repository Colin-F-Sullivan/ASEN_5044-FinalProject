%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: NEESNIS
% Author: Colin Sullivan
% 
% Date Created: 4/20/20
% Date Last Modified: 4/20/20
%
% Purpose:  Perform Linear KF on noisy y data generated from the non-linear
%           system of dnonline_dt.m
%
% Inputs:  None
%
% Outputs:  NEES and NIS tests according to the Monte Carlo simulation
%           described in this report. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
addpath(genpath('../../ASEN_5044-FinalProject'));
%Inputs to Monte Carlo Simulation
NumSims = 50;

%Preallocation
NEESsamps = zeros(NumSims,1401);
NISsamps = zeros(NumSims,1401);
f = waitbar(0,'i = 0');
for i = 1:NumSims
    pert = [0; 0; abs(mvnrnd(0,.2)); 0];
    [NEES,NIS,~] = LKF_Main(0,pert);
    NEESsamps(i,:) = NEES;
    NISsamps(i,:) = NIS;
    waitbar(i/NumSims,f,['Running MC Simulation ' num2str(i) ' out of ' num2str(NumSims)]);
end
close(f);

%NEES Test
E_NEESbar = mean(NEESsamps,2);
alphaNEES = 0.1;
Nnx = NumSims*4;
%Intervals
r1x = chi2inv(alphaNEES/2, Nnx )./ NumSims;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ NumSims;
figure;
plot(E_NEESbar,'ro','MarkerSize',6,'LineWidth',2);
hold on;
plot(r1x*ones(size(E_NEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(E_NEESbar)),'r--','LineWidth',2)
ylabel('NEES','FontSize',14);
grid on;
xlabel('Time Step, k','FontSize',14);
title('NEES Estimation Results','FontSize',14);
lgd = legend('NEES', 'r_1 Bound', 'r_2 Bound');

%NIS Test
E_NIS = mean(NEESsamps,2);
alphaNIS = 0.1;
r1y = chi2inv(alphaNIS/2,3*NumSims)./ NumSims;
r2y = chi2inv(1-alphaNIS/2,3*NumSims)./ NumSims;
figure;
plot(E_NIS,'bo','MarkerSize',6,'LineWidth',2);
hold on;
grid on;
plot(r1y*ones(size(E_NIS)),'b--','LineWidth',2);
plot(r2y*ones(size(E_NIS)),'b--','LineWidth',2);
ylabel('NIS statistic','FontSize',14);
xlabel('Time step, k','FontSize',14);
title('NIS Estimation Results','FontSize',14);
legend('NIS', 'r_1 bound', 'r_2 bound');

