clc; clear; close all;
% Load output data for each training method
load('./data/E_En_U_KF_ref.mat');
load('./data/OTF_ref.mat');
% Calculating RMSE for each training method
yEKFrmse = sqrt(squeeze(sum((yEKF-y1).*(yEKF-y1),1))./size(y1,1))';
yEnKFrmse = sqrt(squeeze(sum((yEnKF-y1).*(yEnKF-y1),1))./size(y1,1))';
yUKFrmse = sqrt(squeeze(sum((yUKF-y1).*(yUKF-y1),1))./size(y1,1))';
yOTFrmse = sqrt(squeeze(sum((yOTF-y1).*(yOTF-y1),1))./size(y1,1))';
% Plots 
figure(1); clf; hold on; grid on; box;
plot(mean(yEKFrmse),'b*-'); 
plot(mean(yEnKFrmse),'k*-'); 
plot(mean(yUKFrmse),'g*-'); 
plot(mean(yOTFrmse),'Marker','*','Color',[0.6 0.2 0]); 
legend('EKF','EnKF','UKF','OTF');
title('Variation of RMSE with epoch');
ylabel('RMSE')
xlabel('Epoch')
set(gcf,'Position',[571   419   849   428]);
%%
clear;
load('./data/testresults.mat')
yEKFrmse = sqrt(squeeze(sum((yEKFtest-y2).*(yEKFtest-y2),1))./size(y2,1))';
yEnKFrmse = sqrt(squeeze(sum((yEnKFtest-y2).*(yEnKFtest-y2),1))./size(y2,1))';
yUKFrmse = sqrt(squeeze(sum((yUKFtest-y2).*(yUKFtest-y2),1))./size(y2,1))';
yOTFrmse = sqrt(squeeze(sum((yOTFtest-y2).*(yOTFtest-y2),1))./size(y2,1))';

figure(2); clf; hold on; grid on; box;
plot([yEKFtest],'b--','LineWidth',1);
plot([yEnKFtest],'k-.','LineWidth',1);
plot([yUKFtest],'g','LineWidth',1);
plot([yOTFtest],'Color',[0.6 0.2 0],'LineWidth',1);
plot(y2,'r-','LineWidth',1);
legend('EKF','EnKF','UKF','OTF','True');
title('NN testing for trained and untrained data');
ylabel('y')
xlabel('Time Step')
set(gcf,'position',[182,188,1436,423]);

figure(3); clf; hold on; grid on; box;
bar(categorical({'EKF','EnKF','UKF','OTF'}),[yEKFrmse,yEnKFrmse,yUKFrmse,yOTFrmse],0.3)
title('Combined RMSE for trained and untrained data');
ylabel('RMSE')