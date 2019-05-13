% Author: Vedang Deshpande
% Date: 13th May 2019
% Training NN using non-linear filters
% Note: States here are the parameters of neural network
clc; clear; close all;

load('data/trained_NN_complete_data.mat')
NN_EKF = NNconstruct(ni,Ln); % We will train this NN using EKF
x0 = nn2param(NN_EKF); % Initial state
nx  = length(x0); % Number of states
ny = Ln(end); % Number of measurements = no. of o/p of NN
Inx  = eye(nx);
Ast = Inx; % State transition matrix is identity

var_meas = 0.1; % variance of measurement noise; sigma^2
var_proc = 0.3; % variance of process noise
var_initState = 1; % initial state covariance

P0 = var_initState*Inx; % Initial state covariance matrix
Q = var_proc*Inx; % Process noise covariance matrix
R = var_meas*eye(ny); % Measurements noise covariance matrix
yMeas = y1 + normrnd(0,var_meas,[length(y1),1]); % Synthetic Noisy measurements
kEnd = length(yMeas);

%% EKF Code
x_prev = x0; % estimate of x(k-1)
P_prev = P0; 
yEKF=zeros(kEnd,1);
for k = 1:kEnd
    clc; 
    fprintf('EKF k = %d.\n',k);
    
    % EKF Propagation/Prediction
    x_pr = Ast*x_prev; % x-(k) a priori state estimate
    P_pr = Ast*P_prev*Ast' + Q; % P-(k) a priori state covariance matrix
    
    % EKF Update
    NN_EKF = param2nn(NN_EKF,x_pr); % Update parameters of the NN
    res = yMeas(k,1) - measModel(NN_EKF,iP(k,:)'); % Innovation/ Measurement residual
    H = nnJacobian(NN_EKF,iP(k,:)'); % Jacobian of measurement model w.r.t states (NN Params in this case)
    KK = P_pr*H'/(H*P_pr*H'+R); % Kalman Gain
    x_pst = x_pr + KK*res; % a posteriori state estimate
    P_pst = (Inx-KK*H)*P_pr; % a posteriori state covariance matrix
    
    % Evaluate o/p of NN using a posteriori parameters estimates
    NN_EKF = param2nn(NN_EKF,x_pst); % Update parameters of the NN
    yEKF(k,1) = measModel(NN_EKF,iP(k,:)');
    
    x_prev = x_pst; 
    P_prev = P_pst; 
end
%% Plot
figure(1); clf; hold on; box; grid;
plot(y1,'r','LineWidth',1) % Clean output
plot(yEKF,'b','LineWidth',1) % Estimated using EKF
plot(yMeas,'k+','MarkerSize',3) % Noisy measurements
legend('Clean','EKF','Noisy')
set(gcf,'Position',[ 244,442,1556,457])
% saveas(gcf,'plots/trainingUsingFilters.png')
