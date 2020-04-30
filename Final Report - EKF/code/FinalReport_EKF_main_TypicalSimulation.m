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

%% Load Data
load('orbitdeterm_finalproj_KFdata.mat')
Q=Qtrue;
R=Rtrue;
clear Qtrue Rtrue temp ydata measLabels i tvec


%% Full Nonlinear Simulation of System Dynamics
[~,OffNominalStateVector,t]=ODE45_EKF(perturbation);

clear perturbation

%% Simulate Measurement Data for Testing
%rng(100)

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
Q=diag(diag(Q))*10;
R=diag(diag(R))*10;

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
PEKF(1:4,1:4)=diag([30,.5,30,.5].^2)+0.001*(ones(4)-eye(4));
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

%Convert P matrix to nice form
TwoSigmaSq=zeros(4,numTimeSteps);

for i=1:10%numTimeSteps
    TwoSigmaSq(:,i)=2*sqrt(abs(diag(PEKF(:,i*4-3:i*4))));    
end

%% Plot Results

ylabels={'X-Position [km]';'X-Velocity [km/s]';'Y-Position [km]';'Y-Velocity [km/s]'};
graphName={"x_1";"x_2";"x_3";"x_4"};
legendTypes={'_{,ode45}';'_{,EKF}';'_{,EKF} + 2\sigma';'_{,EKF} - 2\sigma'};

%Graphs of Nominal State vs EKF w/2 sigma
for j=1:4
figure(j);
plot(t,OffNominalStateVector(j,:),'LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{1}));
hold on;
plot(t,xEKF(j,:),'LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{2}));
plot(t,xEKF(j,:)+TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{3}));
plot(t,xEKF(j,:)-TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{4}));
grid on;
xlabel('Time [s]','FontSize',12);
ylabel(ylabels{j},'FontSize',12);
legend('location','Best') 
xlim([0 5431]);
title(strcat("Estimation of ",graphName{j},' v. Time'),'FontSize',14);
end

%Graphs of Nominal State -EKF w/2 sigma
ylabels={'X-Position Error [km]';'X-Velocity Error[km/s]';'Y-Position Error[km]';'Y-Velocity Error[km/s]'};
legendTypes={'_{,EKF Error}';'_{,EKF Error} + 2\sigma';'_{,EKF Error} - 2\sigma'};

for j=1:4
figure(j+4);
% plot(t,OffNominalStateVector(j,:)-,'LineWidth',1.2,...
%     'DisplayName',strcat(graphName{j},legendTypes{1}));
plot(t,OffNominalStateVector(j,:)-xEKF(j,:),'LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{1}));
hold on;
plot(t,OffNominalStateVector(j,:)-xEKF(j,:)+TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{2}));
plot(t,OffNominalStateVector(j,:)-xEKF(j,:)-TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{3}));
grid on;
xlabel('Time [s]','FontSize',12);
ylabel(ylabels{j},'FontSize',12);
legend('location','Best') 
xlim([0 5431]);
title(strcat("Error in EKF Estimate of ",graphName{j},'(w.r.t. ode45 Solution) v. Time'),'FontSize',14);
end

%% Use in MC Code
if false
%NEES Test
E_NEESbar=NEES;
NumSims=1
alphaNEES = 0.1;
Nnx = NumSims*4;
%Intervals
r1x = chi2inv(alphaNEES/2, Nnx )./ NumSims;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ NumSims;
figure;
plot(E_NEESbar,'ro','MarkerSize',1,'LineWidth',2);
hold on;
plot(r1x*ones(size(E_NEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(E_NEESbar)),'r--','LineWidth',2)
ylabel('NEES','FontSize',14);
grid on;
xlabel('Time Step, k','FontSize',14);
title('NEES Estimation Results','FontSize',14);
lgd = legend('NEES', 'r_1 Bound', 'r_2 Bound');

%NIS Test
E_NIS=NIS;
alphaNIS = 0.1;
r1y = chi2inv(alphaNIS/2,3*NumSims)./ NumSims;
r2y = chi2inv(1-alphaNIS/2,3*NumSims)./ NumSims;
figure;
plot(E_NIS,'bo','MarkerSize',1,'LineWidth',2);
hold on;
grid on;
plot(r1y*ones(size(E_NIS)),'b--','LineWidth',2);
plot(r2y*ones(size(E_NIS)),'b--','LineWidth',2);
ylabel('NIS statistic','FontSize',14);
xlabel('Time step, k','FontSize',14);
title('NIS Estimation Results','FontSize',14);
legend('NIS', 'r_1 bound', 'r_2 bound');
end