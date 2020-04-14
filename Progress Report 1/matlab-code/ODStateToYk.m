%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: StateToYk
% Author: Ben Wise
% Modified from Make_Y_i.m (Colin Sullivan)
% 
% Date Created: 4/12/20
% Date Last Modified: 4/12/20
%
% Purpose: Create the vector of yk measurements for requested stations
%          at some time 't' with corresponding state information
%          'StateVector'. Only requested tracking stations are given.
%
% Inputs:  stations     = Vector of Station Numbers for whick yk should be
%                         calculated
%          t            = Time at which to evalute yk
%          StateVector  = State vector X of the spacecraft
%
% Outputs: yk   = Stacked yk_i's as defined in the project document
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yk] = ODStateToYk(stations,t,StateVector)
    
    %StateVector to Easier to Use Form
    x1 = StateVector(1);
    x2 = StateVector(2);
    x3 = StateVector(3);
    x4 = StateVector(4);

    %Get all station locations at time t
    [Xi,Yi,Xidot,Yidot,~]=ODTrackingStations(t);

    %Get Angles between Sat Position and Tracking Location, Check Tracking
    %Criteria
    phi_i=zeros(size(Xi));
    for i=1:12
        phi_i(i)=atan2((x3-Yi(i)),(x1-Xi(i)));
    end
    
    yk=[];
    
    %Generate Yk matrices
    if ~isempty(stations)
        for iter=stations
            %Equations from doc:
            rho_i = sqrt((x1 - Xi(iter)).^2 + (x3 - Yi(iter)).^2);
            rho_dot_i = (((x1 - Xi(iter)).*(x2-Xidot(iter))) + ((x3 - Yi(iter)).*(x4-Yidot(iter))))./rho_i;

            %Stack them
            yk = [yk;[rho_i;rho_dot_i;phi_i(iter)]];
        end
    end
end
