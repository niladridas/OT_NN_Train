% Author: Niladri Das
% Date: June 1, 2019
clc; clear; close all;
load('data/AllNN.mat');
load('data/trained_NN_complete_data.mat');
%% Performance on Testing data
% Gen testing data
% Out of 583 input data, 200 is used for training, the next 200 is used for
% testing.
y2= y1(201:250);
% y2=y1;
kEnd = length(y2);
iP = iP(201:250,:); % Testing input data
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

% % Plotting
h1 = figure(1); hold on; box; grid;
p1 = plot([yEKF;yEKFtest],'LineWidth',2);
p1.Color = [0, 0, 1];
p1.LineStyle = '--';
p2 = plot([yEnKF;yEnKFtest],'LineWidth',2);
p2.Color = 'k';
p2.LineStyle = '-.';
p3 = plot([yUKF;yUKFtest],'LineWidth',2);
p3.Color = [0, 0.5, 0];
p3.LineStyle = '--';
p4 = plot([yOTF;yOTFtest],'LineWidth',2);
p4.Color = [0.6350, 0.0780, 0.1840] 	;
p4.LineStyle = '--';
p5 = plot(y1(1:250),'LineWidth',2);
p5.Color = 'r';
p5.LineStyle = '-';
ax = gca;
ax.FontSize = 20;
ax.FontWeight = 'bold';
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Normalized Output';
ax.Title.String = 'Training and testing performance of NN trained using different filters';
p6 = plot(200*ones(1,10),linspace(0,1,10),'LineWidth',2);
lgd = legend([p1,p2,p3,p4,p5],'y_{EKF}','y_{EnKF}','y_{UKF}','y_{OT}','y_{real}');
lgd.Location = 'northwest';
txt1 = 'Training \rightarrow';
txt2 = '\leftarrow Testing';
text(200,0.95,txt1,'HorizontalAlignment','right','FontSize',20,'FontWeight','bold');
text(200,0.95,txt2,'HorizontalAlignment','left','FontSize',20,'FontWeight','bold');
drawnow;

%%
keyboard;
fig.PaperPositionMode = 'auto';