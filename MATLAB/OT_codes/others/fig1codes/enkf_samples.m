function a_samples = enkf_samples(X, Y, Y_predict, R, kk_et)
%   enkf_samples:  Returns EnKF updated samples.
% 
%   a_samples = enkf_samples(X, Y, Y_predict, R)
%       Takes X as the prior samples with predicted observations as Y_predict.
%       Y is the actual single observation and R is the variance of the
%       measurement error whose mean is assumed to be zero.
%   
%  Dimension of X is n_states x sample_size
%  Dimension of Y_predict is n_obs x (sample size)
%  Dimension of Y is n_obs x 1
%  Dimension of R is n_obs x n_obs

    Ns = size(X,2);
    nx = size(X,1);
    ny = size(Y,1);

    % Sample mean
    x_mean = mean(X,2);
    % broadcasted observations
    if ny == 1 
        eta_t = normrnd(0,R);
        eta_all = normrnd(0,R,1,Ns);
    else
        eta_t = mvnrnd(zeros(1,ny),R,1)';
        eta_all = mvnrnd(zeros(1,ny),R,Ns)';
    end
    temp_pxyb = zeros(nx,ny,Ns);
    temp_pyy = zeros(ny,ny,Ns);
    temp_pxxb = zeros(nx,nx,Ns);
    for i =1:Ns
        temp_pxyb(:,:,i) = (X(:,i)-x_mean)*(Y_predict(:,i)-Y-eta_t)';
        temp_pyy(:,:,i) = (Y_predict(:,i)-Y-eta_t)*(Y_predict(:,i)-Y-eta_t)';
        temp_pxxb(:,:,i) = (X(:,i)-x_mean)*(X(:,i)-x_mean)';
    end
    P_xyb = (1/(Ns-1))*sum(temp_pxyb,3);
    P_yy = (1/(Ns-1))*sum(temp_pyy,3);
    K = kk_et*P_xyb/P_yy;
    a_samples = X + K*(Y.*ones(ny,Ns)-Y_predict-eta_all);
end