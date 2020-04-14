%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Make_Y_i
% Author: Colin Sullivan
% 
% v2: Modified by Ben Wise
% Changelog: Used ODSatInView Routine to fix some issues with code used to
% test if sat was in view. Output graphs now consistent with expected.
%
% Date Created: 4/12/20
% Date Last Modified: 4/13/20
%
% Purpose: Create the vector of y_i measurements for some given station 'i'
%          over some time vector 't' with corresponding state information
%          'state'. Additionally the script takes into account the viewing
%          angle information and ensure that the tracking station can only
%          record non-zero information if the spacecraft is in its FOV
%
% Inputs:  i     = Station Number of Y_i that should be generated
%          t     = Some vector of times that this should be calculated over
%          state = state vector X of the spacecraft over the time vector t
%
% Outputs: y_i_out   = y_i_out as defined in the project document
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_i_out] = Make_Y_i(i,t,state)
%% Constants
global mu
R_E = 6378;
omega_E = (2*pi)/mu;

%% Make Station Positions
%Define Initial Conditions;
THETA_0 = (i-1)*(pi/6);
%Define positions as a function of time
X_i = R_E.*cos(omega_E.*t + THETA_0);
Xdot_i = (-omega_E*R_E).*sin(omega_E.*t + THETA_0);
Y_i = R_E.*sin(omega_E.*t + THETA_0);
Ydot_i = (omega_E*R_E).*cos(omega_E.*t + THETA_0);
theta_i = atan2(Y_i,X_i);

%% Make y_i
X = state(:,1);
Xdot = state(:,2);
Y = state(:,3);
Ydot = state(:,4);

%Equations from doc:
rho_i = sqrt((X - X_i).^2 + (Y - Y_i).^2);
rho_dot_i = (((X - X_i).*(Xdot-Xdot_i)) + ((Y - Y_i).*(Ydot-Ydot_i)))./rho_i;
phi_i = atan2((Y-Y_i),(X-X_i));

%Stack them
y_i_temp = [rho_i';rho_dot_i';phi_i'];

%Check angular requirements and set all those that dont meet to zero
y_i_out = zeros(size(y_i_temp)); %Preallocate
for i = 1:length(y_i_temp)
    %If its within the bound, assign normally
    if ODSatInView(phi_i(i),theta_i(i))
        y_i_out(:,i) = y_i_temp(:,i);
    else
        y_i_out(:,i) = zeros(3,1);
    end
end
end
