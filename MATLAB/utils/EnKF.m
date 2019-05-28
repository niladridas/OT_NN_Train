% Author: Vedang Deshpande
% Date: 27th May 2019
% Inputs: 
% xInitSamples - initial samples of the state
% fstate - function handle for state dynamics
% hmeas - function handle for measurement model
% yMeas - sequence of noisy measurements
% Q - process noise covariance
% R - measuremet noise covariance
% iP - inputs to the system 
% Outputs: 
% xSamples_Post - posterior samples, after kEnd time steps
% W_pr - equal weights for the samples
%
% Propagates init samples for one epoch. One epoch = no. of measurements

function [xSamples_Post,W_pr] = EnKF(xSamples_init,fstate,hmeas,yMeas,Q,R,iP)
    [nx,nSample] = size(xSamples_init); 
    [kEnd,~] = size(yMeas); % no. of time steps
    W_pr = ones(1,nSample)/nSample; % Prior is equally weighted
    xSamples_prev = xSamples_init;
    for k = 1:kEnd % time steps
        % fprintf('  k = %d.\n',k);
        % Propagation/Prediction
        xSamples_pr = zeros(nx,nSample);
        for i = 1:nSample
            xSamples_pr(:,i) = fstate(xSamples_prev(:,i),iP(k,:)');
        end
        xSamples_pr = xSamples_pr +  mvnrnd(zeros(nx,1),Q,nSample)'; % Add process noise

        % Update
        for i = 1:nSample
            ySamples_pr(:,i) = hmeas(xSamples_pr(:,i),iP(k,:)');
        end
        xSamples_pst = enkf_samples(xSamples_pr, yMeas(k,:)',ySamples_pr,R);       
        xSamples_prev = xSamples_pst;
    end % k
    xSamples_Post = xSamples_pst; % Equally weighted samples after one epoch
end % function