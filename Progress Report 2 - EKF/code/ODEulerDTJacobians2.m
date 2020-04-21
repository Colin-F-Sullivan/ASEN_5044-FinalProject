%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODEulerDTJacobians
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Inputs:   Current Nominal State Vector - NominalStateVector
%           Length of Time Step in sec - DeltaT
%           Current Time - t
% Outputs:  Eulerized Discrete Time Jacobians - Ftilde [State Transition Jacobian]
%                                               Gtilde [Control Input Jacobian]
%                                               Omtilde [Process Noise Jacobian]
%                                               Htilde [Observation Matrix Jacobian]
%                                               Mtilde [Input Feedthrough Jacobian]
%           List of stations making observations - ValidStations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ftilde,Gtilde,Omtilde,Htilde,Mtilde,ObservingStations]...
    = ODEulerDTJacobians2(StateVector, DeltaT,t,sID)
    
    x1=StateVector(1);
    x3=StateVector(3);
    
    %OD Specific Constants
    mu=398600; %km^3 s^-2
        
    % OD Specific Matricies
    Atilde = [0,1,0,0;...
        (2*mu*x1^2-mu*x3^2)/((x1^2+x3^2)^(5/2)),0,(3*mu*x1*x3)/((x1^2+x3^2)^(5/2)),0;...
        0,0,0,1;...
        (3*mu*x1*x3)/((x1^2+x3^2)^(5/2)),0,(2*mu*x3^2-mu*x1^2)/((x1^2+x3^2)^(5/2)),0];
    
    Btilde = [0,0;1,0;0,0;0,1];
    
    
    [Ctilde,ObservingStations] = OD_CtildeMatr2(t,StateVector,sID);
        
    Dtilde = zeros(max(size(ObservingStations))*3,2);
    
    Gamma= Btilde;
        
    %General Formula for Eulerized DT Jacobians
    
    
    Ftilde = eye(size(Atilde))+ DeltaT*Atilde;
    Gtilde = DeltaT*Btilde;
    Omtilde = DeltaT*Gamma;
    Htilde = Ctilde;
    Mtilde = Dtilde;

end