clc; clear; close all;
% Load output data for each training method
load('./data/E_En_U_KF_ref.mat');
load('./data/OTF_ref.mat');
% Calculating RMSE for each training method
yEKFrmse = sqrt(squeeze(sum((yEKF-yMeas).*(yEKF-yMeas),1))./size(yMeas,1))';
yEnKFrmse = sqrt(squeeze(sum((yEnKF-yMeas).*(yEnKF-yMeas),1))./size(yMeas,1))';
yUKFrmse = sqrt(squeeze(sum((yUKF-yMeas).*(yUKF-yMeas),1))./size(yMeas,1))';
yOTFrmse = sqrt(squeeze(sum((yOTF-yMeas).*(yOTF-yMeas),1))./size(yMeas,1))';
%% Plots 

figure(1); clf; hold on; grid on; box;
plot(mean(yEKFrmse),'bo-'); 
plot(mean(yEnKFrmse),'ko-'); 
plot(mean(yUKFrmse),'go-'); 
plot(mean(yOTFrmse),'mo-'); 
legend('EKF','EnKF','UKF','OTF');
title('Variation of RMSE with epoch');
ylabel('RMSE')
xlabel('Epoch')
set(gcf,'Position',[571   419   849   428]);