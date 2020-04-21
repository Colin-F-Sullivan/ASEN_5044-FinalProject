%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ODSatInView
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Purpose: Test whether satellite is in view of station
% Inputs: phi, theta
% Outputs: inView - boolean, 1 for in view, 0 for not in view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inView = ODSatInView(phi,theta)
    
    %Default is not in view
    inView=0;
    
    %If abs(phi-theta) is less than pi/2 then the satellite is in view of 
    %station, i.e. theta-pi/2<phi<theta+pi/2
    if abs(phi-theta)<=pi/2
        inView=1;
    end
    
    %If the abs(phi-theta)>pi then we're measuring the wrong side of the
    %circle, i.e. try the remaining arc 2*pi - abs(diff)
    if abs(phi-theta)>=pi && 2*pi-abs(phi-theta)<=pi/2
        inView=1;
    end
end
