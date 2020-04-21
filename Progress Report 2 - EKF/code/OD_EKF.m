%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: OD_EKF
% Author: Ben Wise
%
% Date Created: 4/20/20
% Date Last Modified: 4/20/20
%
% Purpose:
%
% Inputs:   perturbation vector
% Outputs:  Generates state vector at all times for one period, to be used
%               as nominal trajectory
%           Plots of Non-Linearized X and Y state vectors for part c.) of
%               Progress report 1 of ASEN 5044 Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xplus_kplus1,Pplus_kplus1,NEES_kplus1,NIS_kplus1] =...
    OD_EKF(x_k,P_k,...
            u_k,Q_k,t_k,...
            y_kplus1,R_kplus1,...
            sID_kplus1,DeltaT,...
            xtruth_kplus1, ytruth_kplus1)

%Constants
ode45opts = odeset('RelTol',1e-13,'AbsTol',1e-13);

%xplus,Pplus are at time k, we will find results for time k+1
xplus_k=x_k;
Pplus_k=P_k;

%% PS: Prediction Step

%% PS.1: Deterministic Nonlinear DT dynamics evaluation on xplus_k
%Use ODE45 to find xminus_kplus1, using xplus_k, u_k, w_k=0
w=[0 0]'; 

[~,new_state] = ode45(...
    @(t,state) statOD_dynamics(t,state,u_k,w),...
    [t_k t_k+DeltaT],...
    xplus_k,...
    ode45opts);

xminus_kplus1 = new_state(end,:)';

%% PS.2: Approximation of Predicted Covariance via linearization about xplus_k
%Find linearized Ftilde_k, Omegatilde_k

[Ftilde_k,~,Omega_k,~,~,~] ...
    = ODEulerDTJacobians(xplus_k,DeltaT,t_k);

%Predict Covariance for k+1
Pminus_kplus1= Ftilde_k*Pplus_k*Ftilde_k'+Omega_k*Q_k*Omega_k';


%% MS: Measurement Step
%Need Htilde for following steps
[~,~,~,Htilde_kplus1,~,ObservingStations_kplus1] ...
    = ODEulerDTJacobians2(xminus_kplus1,DeltaT,t_k+DeltaT,sID_kplus1);

%% MS.1: Deterministic Nonlinear Sensor Function Evaluation
%Find yminus_kplus1 using x_kplus1, v_kplus1=0 using full NL Eqns
yminus_kplus1=ODStateToYk(ObservingStations_kplus1,t_k+DeltaT,xminus_kplus1);

if y_kplus1(1:6)== [0;0;0;0;0;0]
    y_kplus1=[];
elseif y_kplus1(4:6)== [0;0;0]
    y_kplus1=y_kplus1(1:3);
    yminus_kplus1=yminus_kplus1(1:3);
    Htilde_kplus1=Htilde_kplus1(1:3,:);
end

if max(size(yminus_kplus1))<max(size(y_kplus1))
    %y_kplus1=y_kplus1(4:6);
end

%% MS.2: Nonlinear measurement innovation: data-predicted
%Innovation Vector
ey_kplus1=y_kplus1-yminus_kplus1;

%% MS.3: Approximation of Kalman Gain using measurement linearization
%Linear Approx Kalman Gain
if max(size(y_kplus1))>4
    R_kplus1=mdiag(R_kplus1,R_kplus1);
end

if ~max(size(Htilde_kplus1))==0
    S_kplus1=Htilde_kplus1 * Pminus_kplus1 * Htilde_kplus1' + R_kplus1;

    K_kplus1=Pminus_kplus1 * Htilde_kplus1' * ...
        ( S_kplus1)^-1;
end

 
%% MS.4: Update Total State and Approximate Covariance via linearization
%Update Total State Estimate
if exist('K_Kplus1','var')
    xplus_kplus1 = xminus_kplus1 + K_kplus1*ey_kplus1;
else
    xplus_kplus1 = xminus_kplus1;
end

%Update Covariance
if exist('K_kplus1','var')
    KH_kplus1 = K_kplus1 * Htilde_kplus1;
    Pplus_kplus1 = (eye(size(KH_kplus1)) - KH_kplus1)*Pminus_kplus1;
else
    Pplus_kplus1 = Pminus_kplus1;
    NEES_kplus1=NaN;
    NIS_kplus1=NaN;
    return;
end

%% NEES
if exist('xtruth_kplus1','var')
    stateVectErr=xtruth_kplus1-xplus_kplus1;
    NEES_kplus1=stateVectErr'*(Pplus_kplus1)^-1*stateVectErr;
else
    NEES_kplus1=NaN;
end

%% NIS

if ytruth_kplus1(4:6)== [0;0;0]
    ytruth_kplus1=ytruth_kplus1(1:3);
end

if exist('ytruth_kplus1','var')
    measErr=ytruth_kplus1-yminus_kplus1;
    NIS_kplus1=measErr'*(S_kplus1)^-1*measErr;
else
    NIS_kplus1=NaN;
end

end
