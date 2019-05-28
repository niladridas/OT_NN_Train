function yEKF = multipleEpochEKF(maxEpoch,kEnd,P0,Q,R,yMeas,iP,NN_EKF)
x0 = nn2param(NN_EKF);
x_prev = x0; % estimate of x(k-1)
P_prev = P0; 
yEKF=zeros(kEnd,maxEpoch);
Inx = eye(length(x0));
Ast = Inx;
for iEp = 1:maxEpoch
    for k = 1:kEnd
        % clc; 
        % fprintf('EKF: Epoch = %d, k = %d.\n',iEp,k);

        % EKF Propagation/Prediction
        x_pr = Ast*x_prev; % x-(k) a priori state estimate
        P_pr = Ast*P_prev*Ast' + Q; % P-(k) a priori state covariance matrix

        % EKF Update
        NN_EKF =  param2nn(NN_EKF,x_pr); % Update parameters of the NN
        res = yMeas(k,1) - measModel(NN_EKF,iP(k,:)'); % Innovation/ Measurement residual
        H = nnJacobian(NN_EKF,iP(k,:)'); % Jacobian of measurement model w.r.t states (NN Params in this case)
        KK = P_pr*H'/(H*P_pr*H'+R); % Kalman Gain
        x_pst = x_pr + KK*res; % a posteriori state estimate
        P_pst = (Inx-KK*H)*P_pr; % a posteriori state covariance matrix

        x_prev = x_pst; 
        P_prev = P_pst; 
    end
    % Evaluate o/p of NN using a posteriori parameters estimates at the end
    % of this epoch, for all input sets
    NN_EKF = param2nn(NN_EKF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yEKF(k,iEp) = measModel(NN_EKF,iP(k,:)');
    end
end % epoch
end % function