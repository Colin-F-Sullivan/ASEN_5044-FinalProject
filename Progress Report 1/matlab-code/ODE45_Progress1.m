%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODE45_Progress1
% Author: Colin Sullivan
% 
% v2: Modified by Ben Wise
% Changelog:  
% -Turned script into function so this can be called from a main routine.
% -Plotting Controlled by 'plotOn' boolean.
% -Added Clean Up Section.
%
% Date Created: 4/12/20
% Date Last Modified: 4/12/20
%
% Inputs:   plotOn - if true, generates plots, if false, does not plot
% Outputs:  Generates state vector at all times for one period, to be used
%               as nominal trajectory
%           Plots of Non-Linearized X and Y state vectors for part c.) of
%               Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state,state_offnom,t]=ODE45_Progress1(plotOn, perturbation)
%Constants
global mu r0
mu = 398600; %[km^3/s^2]
r0 = 6678;

%Time Vectors
P = 2*pi*sqrt(r0^3/mu);
ti = 0;
dT = 10;
tf = P+5*dT; %Plot over one period of circular orbit + 2 time step

%Setup ODE45
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
state0 = [r0; 0; 0; r0.*sqrt(mu/r0^3)];
[t,state] = ode45('dnonline_dt',ti:dT:tf,state0,opts);

%Do off-nominal condition
state0_offnom = state0+perturbation;
[t_offnom,state_offnom] = ode45('dnonline_dt',ti:dT:tf,state0_offnom,opts);

if plotOn
    
    %% ode45 Offnominal-ode45State
    %Plot it
    figure(24);
    subplot(1,2,1);
    plot(t_offnom,state_offnom(:,1)-state(:,1),'LineWidth',1.2);
    hold on;
    plot(t_offnom,state_offnom(:,3)-state(:,3),'LineWidth',1.2);
    grid on;
    xlabel('Time [s]','FontSize',14);
    ylabel('\delta Position [m]','FontSize',14);
    title('Position Deltas','FontSize',14);
    lgd = legend('$\delta X_{ode45}$','$\delta Y_{ode45}$','interpreter','latex');
    lgd.FontSize = 12;
    %xlim([0 P]);
    subplot(1,2,2);
    plot(t_offnom,state_offnom(:,2)-state(:,2),'LineWidth',1.2);
    hold on;
    plot(t_offnom,state_offnom(:,4)-state(:,4),'LineWidth',1.2);
    grid on;
    xlabel('Time [s]','FontSize',14);
    ylabel('\delta Position [m]','FontSize',14);
    title('Velocity Deltas','FontSize',14);
    lgd = legend("$\delta \dot{X}_{ode45}$","$\delta \dot{Y}_{ode45}$",'interpreter','latex');
    lgd.FontSize = 12;
    %xlim([0 P]);
    
    %% Nominal ode45 vs Offnominal LDT
    %Plot X and Y as a function of time;
    figure(1);
    subplot(1,2,1);
    plot(ti:dT:tf,state(:,1),'LineWidth',1.2);
    xlabel('Time [s]','FontSize',14);
    ylabel('Position [km]','FontSize',14);
    hold on;
    grid on;
    plot(ti:dT:tf,state(:,3),'Linewidth',1.2);
    lgd = legend('ode45 X-Position','ode45 Y-Position');
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
    lgd = legend('ode45 X-Velocity','ode45 Y-Velocity');
    lgd.FontSize = 14;
    xlim([0 P]);
    title('Non-Linear Orbital Velocity Vectors v. Time','FontSize',14);
    
    save('state_nonlin_nom.mat','state');
    
    %% ode45 Offnominal vs LDT Offnominal
    %Plot X and Y as a function of time;
    figure(18);
    subplot(1,2,1);
    plot(ti:dT:tf,state_offnom(:,1),'LineWidth',1.2);
    xlabel('Time [s]','FontSize',14);
    ylabel('Position [km]','FontSize',14);
    hold on;
    grid on;
    plot(ti:dT:tf,state_offnom(:,3),'Linewidth',1.2);
    lgd = legend('ode45 X-Position','ode45 Y-Position');
    lgd.FontSize = 14;
    xlim([0 P]);
    title('Perturbed Non-Linear Orbital Position Vectors v. Time','FontSize',14);
    subplot(1,2,2);
    %Plot Velocities as a Function of Time
    plot(ti:dT:tf,state_offnom(:,2),'LineWidth',1.2);
    hold on;
    xlabel('Time [s]','FontSize',14);
    ylabel('Velocity [km/s]','FontSize',14);
    grid on;
    plot(ti:dT:tf,state_offnom(:,4),'LineWidth',1.2);
    lgd = legend('ode45 X-Velocity','ode45 Y-Velocity');
    lgd.FontSize = 14;
    xlim([0 P]);
    title('Perturbed Non-Linear Orbital Velocity Vectors v. Time','FontSize',14);
    
    save('state_nonlin_nom.mat','state');
    
    
    
    
    %% Stations
    % Do the same thing for each station
    y = [];
    for i = [1,7,12]
        y_i = Make_Y_i(i,t,state_offnom);
        %Plot it
        figure(i+2);
        %Range
        subplot(1,3,1);
        plot(t,y_i(1,:),'LineWidth',1.2, "DisplayName","ode45 Solution");
        legend('location','best')
        hold on;
        grid on;
        xlabel('Time [s]','FontSize',14);
        ylabel('Relative Range [km]','FontSize',14);
        tit = ['Relative Range v. Time for Station ' num2str(i)];
        title(tit,'FontSize',14);
        xlim([0 max(t)]);
        %Range Rate
        subplot(1,3,2);
        plot(t,y_i(2,:),'LineWidth',1.2, "DisplayName","ode45 Solution");
        legend('location','best')
        hold on;
        grid on;
        xlabel('Time [s]','FontSize',14);
        ylabel('Relative Range Rate [km/s]','FontSize',14);
        tit = ['Relative Range Rate v. Time for Station ' num2str(i)];
        title(tit,'FontSize',14);
        xlim([0 max(t)]);
        %Elevation Angle
        subplot(1,3,3);
        plot(t,y_i(3,:),'LineWidth',1.2, "DisplayName","ode45 Solution");
        legend('location','best')
        hold on;
        grid on;
        xlabel('Time [s]','FontSize',14);
        ylabel('Elevation Angle [rad]','FontSize',14);
        tit = ['Elevation Angle v. Time for Station ' num2str(i)];
        title(tit,'FontSize',14);
        xlim([0 max(t)]);
        y = [y; y_i];
    end
    
    y_i=[]
    for i = [1,7,12]
        y_i = Make_Y_i(i,t,state);
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
    figure(25);
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
    
end

% Clean Up
clear dT i lgd mu opts P r0 state0 tf ti tit y_i

%Transpose state, t vector for use in PR1_main.m
state=state';
state_offnom=state_offnom';

t=t';
end