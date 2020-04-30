%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: FinalReport2_EKF_main_MonteCarlo
% Author: Ben Wise
% 
% Date Created: 4/29/20
% Date Last Modified: 4/29/20
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
%perturbation=[10;0;0;0]%[0;0.075;0;-0.021]*0;
numTimeSteps=1401;
%NumSims=1

%% Load Data
load('orbitdeterm_finalproj_KFdata.mat')
%Q=Qtrue;
%R=Rtrue;
t=tvec;

ygiven_k=[];
sID=[];

for i=1:max(size(tvec))
    if ~isempty(ydata{i})
        temp=ydata{i}(1:3,:);
        temp2=ydata{i}(4,:);
        if min(size(temp))==1
            ygiven_k=[ygiven_k,[temp;0;0;0]];
            sID=[sID,[temp2;0]];
        else
            ygiven_k=[ygiven_k,temp(:)];
            sID=[sID,temp2(:)];
        end
    else
        ygiven_k=[ygiven_k,[0;0;0;0;0;0]];
        sID=[sID,[0;0]];
    end
end

clear Qtrue Rtrue temp* ydata measLabels i tvec

% %% Full Nonlinear Simulation of System Dynamics
% [~,NominalStateVector,t]=ODE45_EKF(perturbation,numTimeSteps);
% 
% clear perturbation

%% Precalculate all possible quantities to optimize MC speed
xEKF=zeros(4,numTimeSteps);
PEKF=[xEKF,xEKF,xEKF,xEKF];

%No control input for now
u_k=[0 0];

%Guess Q,R
a=.9;
Q=[25,a;
    a,10]*10^-10*.8*.8*.8*.8*.8*.8*1.2

a=25;
b=30;
c=30;

R=[ 175,  a,     b;...
     a,  1100,  c;...
     b,  c,     115]*10^-4*.8^3*.8*.5*2*.8*.8
 
%Initialize with initial location guess
xEKF(:,1)=[6678;0;0;7.72583519755957];
PEKF(1:4,1:4)=diag([50,.2,50,.2].^2)+0.001*(ones(4)-eye(4));
NEES=zeros(1,numTimeSteps);
NIS=NEES;

%% Simulation

%Extended Kalman Filter
for i=1:numTimeSteps-1
    [xEKF(:,i+1),PEKF(:,(i+1)*4-3:(i+1)*4),NEES(1,i+1),NIS(1,i+1)] =...
        OD_EKF(xEKF(:,i),...
        PEKF(:,i*4-3:i*4),...
        u_k,...
        Q,...
        t(i),...
        ygiven_k(:,i+1),...
        R,...
        sID(:,i+1),...
        DeltaT,[],[]);
end

%NEES(MCiterator,:)=tempNEES;
%NIS(MCiterator,:)=tempNIS;

%Convert P matrix to nice form (Not needed for Monte Carlo)
TwoSigmaSq=zeros(4,numTimeSteps);

for i=1:10%numTimeSteps
    TwoSigmaSq(:,i)=2*sqrt(abs(diag(PEKF(:,i*4-3:i*4))));
end




%% Plot Results (Not needed for Monte Carlo)
ylabels={'X-Position [km]';'X-Velocity [km/s]';'Y-Position [km]';'Y-Velocity [km/s]'};
graphName={"x_1";"x_2";"x_3";"x_4"};
legendTypes={'_{NoisySimulatedData}';'_{,EKF}';'_{,EKF} + 2\sigma';'_{,EKF} - 2\sigma'};
if true%NumSims==1

%Graphs of Nominal State vs EKF w/2 sigma
for j=1:4
figure(j);
%plot(t,NomStateVectwithErr(j,:),'LineWidth',1.2,...
%    'DisplayName',strcat(graphName{j},legendTypes{1}));
plot(t,xEKF(j,:),'LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{2}));
hold on;
plot(t,xEKF(j,:)+TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{3}));
plot(t,xEKF(j,:)-TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{4}));
grid on;
xlabel('Time [s]','FontSize',12);
ylabel(ylabels{j},'FontSize',12);
legend('location','Best') 
xlim([0 5431]);
title(strcat("Estimation of ",graphName{j},' from Given Data v. Time'),'FontSize',14);
end
return
%Graphs of Nominal State -EKF w/2 sigma
ylabels={'X-Position Error [km]';'X-Velocity Error[km/s]';'Y-Position Error[km]';'Y-Velocity Error[km/s]'};
legendTypes={'_{,EKF Error}';'_{,EKF Error} + 2\sigma';'_{,EKF Error} - 2\sigma'};

for j=1:4
figure(j+4);
% plot(t,OffNominalStateVector(j,:)-,'LineWidth',1.2,...
%     'DisplayName',strcat(graphName{j},legendTypes{1}));
plot(t,(NominalStateVector(j,:)-xEKF(j,:)),'LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{1}));
hold on;
plot(t,NominalStateVector(j,:)-xEKF(j,:)+TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{2}));
plot(t,NominalStateVector(j,:)-xEKF(j,:)-TwoSigmaSq(j,:),'r--','LineWidth',1.2,...
    'DisplayName',strcat(graphName{j},legendTypes{3}));
grid on;
xlabel('Time [s]','FontSize',12);
ylabel(ylabels{j},'FontSize',12);
legend('location','Best') 
xlim([0 5431]);
title(strcat("Error in EKF Estimate of ",graphName{j},'(w.r.t. ode45 Solution) v. Time'),'FontSize',14);
end
end
%% NEES/NIS Tests

%NEES Test
NEESbar=mean(NEES,1);
alphaNEES = 0.05;
Nnx = NumSims*4;
%Intervals
r1x = chi2inv(alphaNEES/2, Nnx )./ NumSims;
r2x = chi2inv(1-alphaNEES/2, Nnx )./ NumSims;
figure;
plot(NEESbar,'ro','MarkerSize',1,'LineWidth',2);
hold on;
plot(r1x*ones(size(NEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(NEESbar)),'r--','LineWidth',2)
ylabel('NEES','FontSize',14);
grid on;
xlabel('Time Step, k','FontSize',14);
title('NEES Estimation Results','FontSize',14);
lgd = legend('NEES', 'r_1 Bound', 'r_2 Bound');
ylim([0,r1x+r2x]);


%NIS Test
NISbar=mean(NIS,1);
alphaNIS = 0.05;
r1y = chi2inv(alphaNIS/2,3*NumSims)./ NumSims;
r2y = chi2inv(1-alphaNIS/2,3*NumSims)./ NumSims;
figure;
plot(NISbar*4.4,'bo','MarkerSize',1,'LineWidth',2);
hold on;
grid on;
plot(r1y*ones(size(NISbar)),'b--','LineWidth',2);
plot(r2y*ones(size(NISbar)),'b--','LineWidth',2);
ylabel('NIS statistic','FontSize',14);
xlabel('Time step, k','FontSize',14);
title('NIS Estimation Results','FontSize',14);
legend('NIS', 'r_1 bound', 'r_2 bound');
ylim([0,r1y+r2y]);


%% Noisy Ground Truth Data
cm=parula(12);

cmap=zeros(550,3);

for i=1:550
    cmap(i,:)=cm(sID(1,i),:);
end

yl={"$\rho^i(t)$","$\dot{\rho}^i(t)$","$\phi^i(t)$"}

figure()
for i=1:3
    subplot(3,1,i)
    scatter(t,yk(i,:),1,cmap,'filled')
    hold on
    temp=yk(i+3,:);
    tempt=t(temp~=0);
    tempcmap=cmap(temp~=0);
    temp=temp(temp~=0);
    scatter(tempt,temp,1,tempcmap,'filled')
    xlabel("Time [s]")
    ylabel(yl{i}, 'interpreter','latex')
end
sgtitle("Noisy Ground Truth Data")