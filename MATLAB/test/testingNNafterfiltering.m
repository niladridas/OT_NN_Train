% Author: Niladri Das
% Date: June 1, 2019
clc; clear; close all;
load('data/AllNN.mat');
load('data/trained_NN_complete_data.mat');
%% Performance on Testing data
% Gen testing data
% Out of 583 input data, 200 is used for training, the next 200 is used for
% testing.
y2= y1(201:400);
kEnd = length(y2);
iP = iP(201:400,:); % Testing input data
% Initializing test output-data set
yEKFtest = zeros(length(y2),1);
yEnKFtest = zeros(length(y2),1);
yUKFtest = zeros(length(y2),1);
yOTFtest = zeros(length(y2),1);
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
for k = 1:kEnd
    yOTFtest(k) = measModel(NN_OTF,iP(k,:)');
end
disp('OTF Testing Done.')

%% Save reference output-data and the NN output-data
yMeas = y1;
save('./data/testresults.mat','y2','yEKFtest','yEnKFtest','yUKFtest','yOTFtest');

%% Plotting
figure(1); hold on; box; grid;
plot([yEKF;yEKFtest],'b--','LineWidth',1);
plot([yEnKF;yEnKFtest],'k--','LineWidth',1);
plot([yUKF;yUKFtest],'g--','LineWidth',1);
plot([yUKF;yOTFtest],'c--','LineWidth',1);
plot(y1(1:400),'r','LineWidth',1);
legend('y_{EKF}','y_{EnKF}','y_{UKF}','y_{OT}','y_{real}');
drawnow;