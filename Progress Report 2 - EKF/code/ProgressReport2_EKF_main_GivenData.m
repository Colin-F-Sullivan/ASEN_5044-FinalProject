%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: ProgressReport2_EKF_main
% Author: Ben Wise
% 
% Date Created: 4/20/20
% Date Last Modified: 4/20/20
%
% Purpose: Main execution of Progress Report 2, EKF part.
%
% Inputs: N/A
% Outputs: Plots of Non-Linearized X and Y state vectors for part c.) of
%          Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean Up
clear
close all
clc

%% Load Data
load('orbitdeterm_finalproj_KFdata.mat')

Q=Qtrue*100;
R=Rtrue;
t=tvec;

ytrue_k=[];
sID=[];

for i=1:max(size(tvec))
    if ~isempty(ydata{i})
        temp=ydata{i}(1:3,:);
        temp2=ydata{i}(4,:);
        if min(size(temp))==1
            ytrue_k=[ytrue_k,[temp;0;0;0]];
            sID=[sID,[temp2;0]];
        else
            ytrue_k=[ytrue_k,temp(:)];
            sID=[sID,temp2(:)];
        end
    else
        ytrue_k=[ytrue_k,[0;0;0;0;0;0]];
        sID=[sID,[0;0]];
    end
end

clear Qtrue Rtrue temp ydata measLabels i


%% Constants
DeltaT=10;
%perturbation=[0;0.075;0;-0.021];
numTimeSteps=max(size(tvec));

clear tvec

%% Extended Kalman Filter
%Preallocate for speed
xEKF=zeros(4,numTimeSteps);
PEKF=[xEKF,xEKF,xEKF,xEKF];

%Initialize with initial guess
mu = 398600; %[km^3/s^2]
r0 = 6678;
state0 = [r0; 0; 0; r0.*sqrt(mu/r0^3)];

xEKF(:,1)=state0;
PEKF(1:4,1:4)=diag([300,1,300,1]);
NEES=zeros(1,numTimeSteps);
NIS=NEES;

clear mu r0 state0

%No control input for now
u_k=[0 0];

for i=1:numTimeSteps-1
    [xEKF(:,i+1),PEKF(:,(i+1)*4-3:(i+1)*4),NEES(i+1),NIS(i+1)] =...
        OD_EKF(xEKF(:,i),...
        PEKF(:,i*4-3:i*4),...
        u_k,...
        Q,...
        t(i),...
        ytrue_k(:,i+1),...
        R,...
        sID(:,i+1),...
        DeltaT,...
        [],...
        ytrue_k(:,i+1));
end




%% Plot Results
figure()
for j=1:4
subplot(2,1,mod(j+1,2)+1)
%plot(OffNominalStateVector(j,:))
hold on
plot(xEKF(j,:))
end

figure()
for j=1:4
subplot(4,1,mod(j+1,2)+1)
plot(xEKF(j,:))
hold on
end
subplot(4,1,3)
plot(NEES)
subplot(4,1,4)
plot(NIS)
