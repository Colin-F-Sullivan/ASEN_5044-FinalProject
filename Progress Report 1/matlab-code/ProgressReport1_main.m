%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ProgressReport1_main
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Purpose: Main execution of Progress Report 1, parts b,c (Part a is
% written work only.)
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
perturbation=[0;0;100;0];

%% Full Nonlinear Simulation of System Dynamics
[NominalStateVector,OffNominalStateVector,t]=ODE45_Progress1(false, perturbation);

%% Linearized DT Dynamics and Measurement Model
x=zeros(size(NominalStateVector));
x(:,1)=NominalStateVector(:,1);
x(:,1)=x(:,1)+perturbation;

y=zeros(3*12,size(NominalStateVector,2));

%Now calculate for the remaining time steps
for i=2:length(t)-1
    [x(:,i),ytemp,ObservingStation]=ODSimLinearDTDynamics2(...
        NominalStateVector(:,i-1),...
        NominalStateVector(:,i),...
        x(:,i-1),...
        zeros(2,1),zeros(2,1),... %No control input perturbation, therefore unom=uk=0
        t(i),DeltaT);
        
        k=1;
        for j=1:12
            if ismember(j,ObservingStation)
                y(j*3-2:j*3,i)=ytemp(k*3-2:k*3);
                k=k+1;
            end
        end
    
end

%% Compare and Validate Jacobians and DT model against ode45 Sim

%Augment Graphs from ODE45_Progress1
[~,~,~]=ODE45_Progress1(true,perturbation);

%% ode45 Nominal vs LDT Offnominal
figure(1)
subplot(1,2,1);
plot(t,x(1,:),'LineWidth',1.2, "DisplayName", "Linearized DT X-Position")
plot(t,x(3,:),'LineWidth',1.2, "DisplayName", "Linearized DT Y-Position")
subplot(1,2,2);
plot(t,x(2,:),'LineWidth',1.2, "DisplayName", "Linearized DT X-Velocity")
plot(t,x(4,:),'LineWidth',1.2, "DisplayName", "Linearized DT Y-Velocity")

%% ode45 Offnominal vs LDT Offnominal
figure(18)
hold on
subplot(1,2,1);
plot(t,x(1,:),'LineWidth',1.2, "DisplayName", "Linearized DT X-Position")
plot(t,x(3,:),'LineWidth',1.2, "DisplayName", "Linearized DT Y-Position")
subplot(1,2,2);
plot(t,x(2,:),'LineWidth',1.2, "DisplayName", "Linearized DT X-Velocity")
plot(t,x(4,:),'LineWidth',1.2, "DisplayName", "Linearized DT Y-Velocity")

%% ode45 Nominal vs LDT Offnominal Error
%Compare Error
dx=x-NominalStateVector;
figure(24);
subplot(1,2,1);
 plot(t,dx(1,:),'LineWidth',1.2,"DisplayName","$\delta X_{Linearized DT}$");
% xlabel('Time [s]','FontSize',14);
% ylabel('$\delta X$, $\delta Y$ [km]','FontSize',14,'interpreter','latex');
% hold on;
% grid on;
plot(t,dx(3,:),'Linewidth',1.2,"DisplayName","$\delta Y_{Linearized DT}$");
%lgd = legend('','','interpreter','latex');
% lgd.FontSize = 14;
 xlim([0 max(t)-20]);
% title('Linearized DT Position Deltas','FontSize',14);
legend('location','best','interpreter','latex')
subplot(1,2,2);
plot(t,dx(2,:),'LineWidth',1.2,"DisplayName","$\delta \dot{X}_{Linearized DT}$");
% hold on;
% xlabel('Time [s]','FontSize',14);
% ylabel('$\delta \dot{X}$, $\delta \dot{Y}$ [km/s]','FontSize',14,'interpreter','latex');
% grid on;
 plot(t,dx(4,:),'LineWidth',1.2,"DisplayName","$\delta \dot{Y}_{Linearized DT}$");
 legend('location','best','interpreter','latex')
%lgd = legend('','$\delta \dot{Y}_{Linearized DT}$','interpreter','latex');
% lgd.FontSize = 14;
 xlim([0 max(t)-20]);
% title('Linearized DT Velocity Deltas v. Time','FontSize',14);

%% ode45 Offnominal vs LDT Offnominal Error
%Compare Error
dx=x-OffNominalStateVector;
figure(19);
subplot(1,2,1);
plot(t,dx(1,:),'LineWidth',1.2);
xlabel('Time [s]','FontSize',14);
ylabel('$\delta X$, $\delta Y$ [km]','FontSize',14,'interpreter','latex');
hold on;
grid on;
plot(t,dx(3,:),'Linewidth',1.2);
lgd = legend('$\hat{X}=X_{Linearized DT}-X_{ode45}$','$\hat{Y}=Y_{Linearized DT}-Y_{ode45}$','interpreter','latex');
lgd.FontSize = 14;
xlim([0 max(t)-20]);
title('Difference Between Linearized DT and ode45 Position Solution v. Time','FontSize',14);
subplot(1,2,2);
plot(t,dx(2,:),'LineWidth',1.2);
hold on;
xlabel('Time [s]','FontSize',14);
ylabel('$\delta \dot{X}$, $\delta \dot{Y}$ [km/s]','FontSize',14,'interpreter','latex');
grid on;
plot(t,dx(4,:),'LineWidth',1.2);
lgd = legend('$\hat{\dot{X}}=\dot{X}_{Linearized DT}-\dot{X}_{ode45}$','$\hat{\dot{Y}}=\dot{Y}_{Linearized DT}-\dot{Y}_{ode45}$','interpreter','latex');
lgd.FontSize = 14;
xlim([0 max(t)-20]);
title('Difference Between Linearized DT and ode45 Velocity Solution v. Time','FontSize',14);

%% Stations
%Measurement Stations
for i = [1,7,12]
    figure(i+2);
    %Range
    subplot(1,3,1);
    plot(t,y(i*3-2,:),'LineWidth',1.2, "DisplayName","Linear DT Solution");
    %xlim([0 max(t)]);
    %Range Rate
    subplot(1,3,2);
    plot(t,y(i*3-1,:),'LineWidth',1.2, "DisplayName","Linear DT Solution");
    %xlim([0 max(t)]);
    %Elevation Angle
    subplot(1,3,3);
    plot(t,y(i*3,:),'LineWidth',1.2, "DisplayName","Linear DT Solution");
    %xlim([0 max(t)]);
end

