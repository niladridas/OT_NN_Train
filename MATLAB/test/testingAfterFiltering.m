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
[EFKepoc,EKFrep]=find(EKFminimum==yEKFrmse);
EnKFminimum = min(min(yEnKFrmse));
[EnFKepoc,EnKFrep]=find(EnKFminimum==yEnKFrmse);
UKFminimum = min(min(yUKFrmse));
[UKFepoc,UKFrep]=find(UKFminimum==yUKFrmse);
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
[yEKF,NN_EKF] = multipleEpochEKFv1(EFKepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,EKFrep));
disp('EKF Training Done.')

%% EnKF Code
[yEnKF,NN_EnKF] = multipleEpochEnKFv1(EnFKepoc,kEnd,nSample,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,EnKFrep));
disp('EnKF Training Done.')

%% UKF Code 
[yUKF, NN_UKF] = multipleEpochUKFv1(UKFepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,UKFrep));
disp('UKF Training Done.')

%% OT Code
% TO-DO
% [yOT, NN_OT] = multipleEpochOTv1(OTepoc,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,OTrep));
% disp('UKF Done.')

%% Performance on Testing data
% Gen testing data
% Out of 583 input data, 200 is used for training, the next 200 is used for
% testing.
y1 = y1(201:400);
kEnd = length(y1);
iP = iP(201:400,:); % Testing input data
% Initializing test output-data set
yEKFtest = zeros(length(y1),1);
yEnKFtest = zeros(length(y1),1);
yUKFtest = zeros(length(y1),1);
%% EKF Code: just for testing
for k = 1:kEnd
    yEKFtest(k) = measModel(NN_EKF,iP(k,:)');
end
disp('EKF Testing Done.')

%% EnKF Code: just for testing
for k = 1:kEnd
    yEnKFtest(k) = measModel(NN_EnKF,iP(k,:)');
end
disp('EnKF Testing Done.')

%% UKF Code: just for testing 
for k = 1:kEnd
    yUKFtest(k) = measModel(NN_UKF,iP(k,:)');
end
disp('UKF Testing Done.')

%% OT Code
% TO-DO

%% Save reference output-data and the NN output-data
yMeas = y1;
save('./data/testresults.mat','y1','yEKFtest','yEnKFtest','yUKFtest');

