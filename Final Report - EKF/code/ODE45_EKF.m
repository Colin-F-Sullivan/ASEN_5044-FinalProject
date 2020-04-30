%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODE45_EKF
% Author: Colin Sullivan
% 
% v3: Modified by Ben Wise
% Changelog:
% -Changed name to avoid code collisions.
% -Removed plotting functionality to reduce code complexity
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
% Inputs:   perturbation vector
% Outputs:  Generates state vector at all times for one period, to be used
%               as nominal trajectory
%           Plots of Non-Linearized X and Y state vectors for part c.) of
%               Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state,state_offnom,t]=ODE45_EKF(perturbation,numTimeSteps)
%Constants
global mu r0
mu = 398600; %[km^3/s^2]
r0 = 6678;

%Time Vectors
P = 2*pi*sqrt(r0^3/mu);
ti = 0;
dT = 10;
tf = (numTimeSteps-1)*10; %Plot over one period of circular orbit + 2 time step

%Setup ODE45
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
state0 = [r0; 0; 0; r0.*sqrt(mu/r0^3)];
[t,state] = ode45('dnonline_dt',ti:dT:tf,state0,opts);

%Do off-nominal condition
state0_offnom = state0+perturbation;
[~,state_offnom] = ode45('dnonline_dt',ti:dT:tf,state0_offnom,opts);

% Clean Up
clear dT i lgd mu opts P r0 state0 tf ti tit y_i

%Transpose state, t vector for use in PR1_main.m
state=state';
state_offnom=state_offnom';
t=t';

end