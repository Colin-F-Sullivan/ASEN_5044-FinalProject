%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODTrackingStations
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Purpose: Generate Tracking Station Locations, Velocities and Valid
%          Observation Window (in radians) at current time t.
%
% Inputs: Current Time - t
% Outputs: Vectors of Current Tracking Station X,Y Locations - Xi,Yi
%          Vectors of Current Tracking Station X,Y Velocities - Xidot,Yidot
%          Vector of Current Tracking Station Angular Location - thetai
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xi,Yi,Xidot,Yidot,thetai] = ODTrackingStations(t)
    %Data Matricies
    Xi=zeros(12,1);
    Yi=Xi;
    Xidot=Xi;
    Yidot=Xi;
    thetai_init=Xi;
    thetai=Xi;
   
    %Constants
    Re=6378; %km
    wE=2*pi/86400; %rad/s
       
    for i=1:12
        thetai_init(i)=(i-1)*pi/6;
        Xi(i)=Re*cos(wE*t+thetai_init(i));
        Xidot(i)=-Re*wE*sin(wE*t+thetai_init(i));
        Yi(i)=Re*sin(wE*t+thetai_init(i));
        Yidot(i)=Re*wE*cos(wE*t+thetai_init(i));        
        thetai(i)=atan2(Yi(i),Xi(i));
    end
end