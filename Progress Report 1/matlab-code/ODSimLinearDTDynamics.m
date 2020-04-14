%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODSimLinearDTDynamics
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Purpose:  Simulate Linear DT Dynamics for one time step. 3 indicies are
%           used, 0,1,2. k=0 is the previous time step, k=1 is the current
%           time step, k=2 is the time step we are currently wanting to
%           calculate [next time step]. This program is written in this way
%           to reduce number of calls to ODEulerDTJacobians
% Inputs:   xnomk2 - Nominal state vector for the time step we want to calc
%           xnomk1 - Current time step nominal state vector
%           xk1 - Current time step actual state vector
%           unomk1 - Current time step nominal control input vector
%           uk1 - Current time step actual control input vector
%           ynomk0 - Previous time step nominal observation vector
%           t1 -  Current time step time
%           DeltaT - Length of one time step
%           wk1 - Current process noise
%           vk1 - Current step observation noise
%
% Outputs:  xk2 - Actual state vector we want to calculate [for next time step]
%           yk1 - This time step actual observation vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk2,yk1,ObservingStations]...
    =ODSimLinearDTDynamics(xnomk0,xnomk1,xnomk2,xk1,unomk1,uk1,t1,DeltaT)%,wk1,vk1) <-Assume no measurement,process noise

    dxk=xk1-xnomk1;
    duk=uk1-unomk1;
    
    [Ftilde,Gtilde,Omtilde,Htilde,Mtilde,ObservingStations]...
    = ODEulerDTJacobians(xnomk1, DeltaT,t1);

    xk2=xnomk2+Ftilde*dxk+Gtilde*duk; %+Omtilde*wk1;  <-Assume no measurement,process noise
    yk1=ODStateToYk(ObservingStations,t1-DeltaT,xnomk0)...
        +Htilde*dxk+Mtilde*duk; %+vk1;  <-Assume no measurement,process noise
    
end