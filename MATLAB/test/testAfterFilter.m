
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
% Find epoc and rep number for each trainig method
EKFminimum = min(min(yEKFrmse));[EFKepoc,EKFrep]=find(EKFminimum==yEKFrmse);
EnKFminimum = min(min(yEnKFrmse));[EnFKepoc,EnKFrep]=find(EnKFminimum==yEnKFrmse);
UKFminimum = min(min(yUKFrmse));[UKFepoc,UKFrep]=find(UKFminimum==yUKFrmse);