% Author: Niladri Das
% Date: June 1, 2019
% Testing NN which are trained using filters
clc; clear; close all;
% Load output data for each training method
load('./data/E_En_U_KF_ref.mat');
% Calculating RMSE for each training method
yEKFrmse = sqrt(squeeze(sum((yEKF-yMeas).*(yEKF-yMeas),1))./size(yMeas,1))';
yEnKFrmse = sqrt(squeeze(sum((yEnKF-yMeas).*(yEnKF-yMeas),1))./size(yMeas,1))';
yUKFrmse = sqrt(squeeze(sum((yUKF-yMeas).*(yUKF-yMeas),1))./size(yMeas,1))';
% yOTrmse = sqrt(squeeze(sum((yOT-yMeas).*(yOT-yMeas),1))./size(yMeas,1))';

%% Find epoc and rep number for each trainig method
EKFminimum = min(min(yEKFrmse));
[EKFrep,EFKepoc]=find(EKFminimum==yEKFrmse);
EnKFminimum = min(min(yEnKFrmse));
[EnKFrep,EnFKepoc]=find(EnKFminimum==yEnKFrmse);
UKFminimum = min(min(yUKFrmse));
[UKFrep,UKFepoc]=find(UKFminimum==yUKFrmse);
% OTminimum = min(min(yOTrmse));
% [OTepoc,OTrep]=find(OTminimum==yOTrmse);
%% We have the repitition number and epoch number to train for each filter
%% Initialize
load('data/trained_NN_complete_data.mat')
x0 = nn2param(NN); % Initial state
nx  = length(x0); % Number of states
ny = Ln(end); % Number of measurements = no. of o/p of NN
Inx  = eye(nx);
Ast = Inx; % State transition matrix is identity

var_meas = 0.005; % variance of measurement noise; sigma^2
var_proc = 0.01; % variance of process noise
var_initState = 1; % initial state covariance

P0 = var_initState*Inx; % Initial state covariance matrix
Q = var_proc*Inx; % Process noise covariance matrix
R = var_meas*eye(ny); % Measurements noise covariance matrix
y1 = y1(1:200);
yMeas = y1 + normrnd(0,sqrt(var_meas),[length(y1),1]); % Synthetic Noisy measurements
kEnd = length(yMeas);
%% EKF Code: just for training
tic ;[yEKF,NN_EKF] = multipleEpochEKFv1(EFKepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,EKFrep));
toc;disp('EKF Training Done.')
yEKF = yEKF(:,end);

%% EnKF Code
tic;[yEnKF,NN_EnKF] = multipleEpochEnKFv1(EnFKepoc,kEnd,nSample,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,EnKFrep));
toc;disp('EnKF Training Done.')
yEnKF = yEnKF(:,end);

%% UKF Code 
tic;[yUKF, NN_UKF] = multipleEpochUKFv1(UKFepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,UKFrep));
toc;disp('UKF Training Done.')
yUKF = yUKF(:,end);

%% OT Code
% TO-DO
% [yOT, NN_OT] = multipleEpochOTv1(OTepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,OTrep));
% disp('UKF Done.')

save('./data/AllNN.mat','NN_EKF','yEKF','NN_EnKF','yEnKF','NN_UKF','yUKF','y1','yMeas');

