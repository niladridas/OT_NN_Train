% Author: Vedang Deshpande
% Date: 14th May 2019
% Training NN using non-linear filters
% Note: States here are the parameters of neural network
clc; clear; close all;
%% Initialize
% rng('default');
load('data/trained_NN_complete_data.mat')
% NN_EKF = NNconstruct(ni,Ln,0); % We will train this NN using EKF
% w_true = nn2param(NN); 
% NN_EKF = param2nn(NN_EKF,w_true+1*rand(length(w_true),1));
% NN_UKF = NN_EKF; NN_OTF = NN_EKF; NN_EnKF = NN_EKF;

Ln = [5;5;1]; % no. of nodes in NN
x0 = nn2param(NNconstruct(ni,Ln,rand)); % Initial state
nx  = length(x0); % Number of states
ny = Ln(end); % Number of measurements = no. of o/p of NN
Inx  = eye(nx);
Ast = Inx; % State transition matrix is identity

var_meas = 0.005; % variance of measurement noise; sigma^2
var_proc = 0.01; % variance of process noise
var_initState = 1; % initial state covariance

P0 = var_initState*Inx; % Initial state covariance matrix
% P0 = 2*rand(nx,nx)-1; P0 = (P0*P0.'); P0 = var_initState*P0/max(eig(P0));
Q = var_proc*Inx; % Process noise covariance matrix
R = var_meas*eye(ny); % Measurements noise covariance matrix
y1 = oP(300:500,:);
iP = iP(300:500,:);
yMeas = y1 + normrnd(0,sqrt(var_meas),[length(y1),1]); % Synthetic Noisy measurements
kEnd = length(yMeas);
maxEpoch = 10;
nRepeat = 1;
arrRepeat = 1:20;

% figure(1); clf; hold on; box; grid;
% plot(y1,'r','LineWidth',1) % Clean output
% plot(yMeas,'k+','MarkerSize',3) % Noisy measurements
% figure(2); clf; hold on; box; grid;
% drawnow;

nSample = 2*nx+1;

%% EKF Code
yEKF=zeros(kEnd,maxEpoch,nRepeat);
parfor (iRep = 1:nRepeat,4)
    fprintf('EKF: Rep = %d\n',iRep);
    yEKF(:,:,iRep) = multipleEpochEKF(maxEpoch,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,iRep));
end % repeat

for iRep = 1:nRepeat
    for iEp = 1:maxEpoch
         RMSE_EKF(1,iEp,iRep) = norm(yEKF(:,iEp,iRep) - y1)/sqrt(kEnd);
    end % epoch
end
disp('EKF Done.')

% figure(1); hold on; box; grid;
% plot(yEKF(:,end),'b--','LineWidth',1) 
% figure(2); hold on; box; grid;
% plot(abs(yEKF(:,end)-y1)./abs(y1),'b--','LineWidth',1) 
% drawnow;

%% EnKF Code
yEnKF=zeros(kEnd,maxEpoch,nRepeat);
parfor iRep = 1:nRepeat
    fprintf('EnKF: Rep = %d\n',iRep);
    yEnKF(:,:,iRep) = multipleEpochEnKF(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,iRep));
    % 
end % repeat

for iRep = 1:nRepeat
    for iEp = 1:maxEpoch
         RMSE_EnKF(1,iEp,iRep) = norm(yEnKF(:,iEp,iRep) - y1)/sqrt(kEnd);
    end % epoch
end
disp('EnKF Done.')
% figure(1); hold on; box; grid;
% plot(yEnKF(:,end),'k--','LineWidth',1) 
% figure(2); hold on; box; grid;
% plot(abs(yEnKF(:,end)-y1)./abs(y1),'k--','LineWidth',1) 
% drawnow;

%% UKF Code 
yUKF=zeros(kEnd,maxEpoch,nRepeat);
parfor iRep = 1:nRepeat
    fprintf('UKF: Rep = %d\n',iRep);
    yUKF(:,:,iRep) = multipleEpochUKF(maxEpoch,kEnd,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,iRep));
end % repeat

for iRep = 1:nRepeat
    for iEp = 1:maxEpoch
         RMSE_UKF(1,iEp,iRep) = norm(yUKF(:,iEp,iRep) - y1)/sqrt(kEnd);
    end % epoch
end
disp('UKF Done.')
% figure(1); hold on; box; grid;
% plot(yUKF(:,end),'g--','LineWidth',1) 
% figure(2); hold on; box; grid;
% plot(abs(yUKF(:,end)-y1)./abs(y1),'g--','LineWidth',1) 
% drawnow;

%% OTF Code
tic
yOTF=zeros(kEnd,maxEpoch,nRepeat);
% parfor (iRep = 1:nRepeat,2)
for iRep = 6:10
    clc
    fprintf('OTF: Rep = %d\n',iRep);
    yOTF(:,:,iRep) = multipleEpochOTF(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NNconstruct(ni,Ln,iRep));
end % repeat

for iRep = 1:nRepeat
    for iEp = 1:maxEpoch
         RMSE_OTF(1,iEp,iRep) = norm(yOTF(:,iEp,iRep) - y1)/sqrt(kEnd);
    end % epoch
end
disp('OTF Done.')
toc
save OTF_ref_6_to_10.mat
% figure(1); hold on; box; grid;
% plot(yOTF(:,Ep_OTF),'m--','LineWidth',1) 
% figure(2); hold on; box; grid;
% plot(abs(yOTF(:,Ep_OTF)-y1)./abs(y1),'m--','LineWidth',1) 
% drawnow;

%% Plot
% figure(1); clf; hold on; box; grid;
% plot(y1,'r','LineWidth',1) % Clean output
% plot(yEKF(:,end),'b-','LineWidth',1) % Estimated using EKF
% plot(yEnKF(:,end),'b--','LineWidth',1) % Estimated using EnKF
% plot(yUKF(:,end),'g--','LineWidth',1) % Estimated using UKF
% plot(yOTF(:,end),'m--','LineWidth',1) % Estimated using OTF
% plot(yMeas,'k+','MarkerSize',3) % Noisy measurements
% legend('Clean','EKF','UKF','Noisy')
% set(gcf,'Position',[370,554,929,370])
% % saveas(gcf,'plots/trainingUsingFilters.png')
% 
% figure(2); clf; hold on; box; grid;
% plot(abs(yEKF(:,end)-y1)./abs(y1),'b','LineWidth',1) % Estimated using EKF
% plot(abs(yUKF(:,end)-y1)./abs(y1),'g--','LineWidth',1) % Estimated using UKF
% plot(abs(yOTF(:,end)-y1)./abs(y1),'m--','LineWidth',1) % Estimated using OTF
% legend('EKF','UKF')
% set(gcf,'Position',[368,118,933,347])

figure(3); clf; hold on; grid on;
plot(mean(RMSE_EKF,3),'bo-'); 
plot(mean(RMSE_EnKF,3),'ko-'); 
plot(mean(RMSE_UKF,3),'go-'); 
plot(mean(RMSE_OTF,3),'mo-'); 
legend('EKF','EnKF','UKF','OTF')