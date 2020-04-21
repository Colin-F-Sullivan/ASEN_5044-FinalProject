%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ProgressReport2_EKF_main
% Author: Ben Wise
% 
% Date Created: 4/20/20
% Date Last Modified: 4/20/20
%
% Purpose: Main execution of Progress Report 2, EKF part.
%
% Inputs: N/A
% Outputs: Plots of Non-Linearized X and Y state vectors for part c.) of
%          Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean Up
clear
close all
clc

%% Constants
DeltaT=10;
perturbation=[0;0.075;0;-0.021];
numTimeSteps=1401;

%% Full Nonlinear Simulation of System Dynamics
[~,OffNominalStateVector,t]=ODE45_EKF(perturbation);

clear perturbation

%% Simulate Measurement Data for Testing
%rng(100)
% 
% Q=[.001,0,0,0;...
%     0,0.1,0,.5;...
%     0,0,.001,0;...
%     0,.5,0,1]/100;
Q=[10,.5;.5,10]/2000;
R=[1,.015,.035;.015,2,.025;.035,.025,0.5]/1000;%2000;

Sw=chol(Q,'lower');
Sv=chol(R,'lower');

wk=Sw*randn(2,numTimeSteps);

wk=[zeros(1,1401);wk(1,:);zeros(1,1401);wk(2,:);];
vk=Sv*randn(3,numTimeSteps);

clear Sw Sv

NomStateVectwithErr=OffNominalStateVector+wk;

yk=zeros(6,numTimeSteps);
sID=zeros(2,numTimeSteps);

for k=0:numTimeSteps-1
    [yk(:,k+1),sID(:,k+1)]=ODGenerateYk([],t(k+1),NomStateVectwithErr(:,k+1),vk(:,k+1));
end

clear k

Q=Q/5;
R=R/5;

%% Generate True yk
ytrue_k=zeros(size(yk));

for k=0:numTimeSteps-1
    [ytrue_k(:,k+1),~]=ODGenerateYk([],t(k+1),OffNominalStateVector(:,k+1),[0;0;0]);
end


%% Extended Kalman Filter
%Preallocate for speed
xEKF=zeros(4,numTimeSteps);
PEKF=[xEKF,xEKF,xEKF,xEKF];

%Initialize with initial guess
xEKF(:,1)=OffNominalStateVector(:,1);
PEKF(1:4,1:4)=diag([300,1,300,1]);
NEES=zeros(1,numTimeSteps);
NIS=NEES;

%No control input for now
u_k=[0 0];

for i=1:numTimeSteps-1
    [xEKF(:,i+1),PEKF(:,(i+1)*4-3:(i+1)*4),NEES(i+1),NIS(i+1)] =...
        OD_EKF(xEKF(:,i),...
        PEKF(:,i*4-3:i*4),...
        u_k,...
        Q,...
        t(i),...
        yk(:,i+1),...
        R,...
        sID(:,i+1),...
        DeltaT,...
        OffNominalStateVector(:,i+1),...
        ytrue_k(:,i+1));
end


%% Plot Results
figure()
for j=1:4
subplot(2,1,mod(j+1,2)+1)
plot(OffNominalStateVector(j,:))
hold on
plot(xEKF(j,:))
end

figure()
for j=1:4
subplot(4,1,mod(j+1,2)+1)
plot(OffNominalStateVector(j,:)-xEKF(j,:))
hold on
end
subplot(4,1,3)
plot(NEES)
subplot(4,1,4)
plot(NIS)
