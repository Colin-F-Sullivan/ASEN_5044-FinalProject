%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: OD_CtildeMatr
% Author: Ben Wise
% 
% Date Created: 4/13/20
% Date Last Modified: 4/13/20
%
% Purpose: Generate the current Ctilde Matrix, given State Vector and
%          current time t. Calculates Tracking Station Locations first,
%          then checks for valid observation angle. Only outputs Ctilde for
%          valid (Observing) Tracking Stations.
%
% Inputs: Current Time - t
%         State Vector - StateVector
% Outputs: Current Ctilde Matrix (scaled as appropriate for number of 
%          available tracking stations that can see Sat) - Ctilde
%          List of Valid (Observing) Stations by Number - ObservingStations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ctilde,ObservingStations] = OD_CtildeMatr(t,StateVector)
    
    %StateVector to Easier to Use Form
    x1 = StateVector(1);
    x2 = StateVector(2);
    x3 = StateVector(3);
    x4 = StateVector(4);

    %Calculate current Tracking Station Locations
    [Xi,Yi,Xidot,Yidot,thetai]=ODTrackingStations(t);
    
    %Get Angles between Sat Position and Tracking Location, Check Tracking
    %Criteria
    phi_i=zeros(size(Xi));
    for i=1:12
        phi_i(i)=atan2((x3-Yi(i)),(x1-Xi(i)));
    end

    %Build C Matrix
    Ctilde=[];
    ObservingStations=[];
    
    for i=1:12
       %But only if tracking angles are valid 
       if ODSatInView(phi_i(i),thetai(i))
           ObservingStations=[ObservingStations, i];
           denom=sqrt((x1-Xi(i))^2+(x3-Yi(i))^2);
           
           Ctemp=[(x1-Xi(i))/denom , 0, (x3-Yi(i))/denom, 0; ...
               (x2-Xidot(i))/denom-(((x1-Xi(i))*(x2-Xidot(i))+(x3-Yi(i))*(x4-Yidot(i)))*(x1-Xi(i)))/(denom^3),...
               (x1-Xi(i))/denom,...
               (x4-Yidot(i))/denom-(((x1-Xi(i))*(x2-Xidot(i))+(x3-Yi(i))*(x4-Yidot(i)))*(x3-Yi(i)))/(denom^3),...
               (x3-Yi(i))/denom;...
               (Yi(i)-x3)/denom^2,0, (x1-Xi(i))/denom^2,0];               
           
           Ctilde=[Ctilde;Ctemp];
       end
    end
   
end