function yUKF = multipleEpochUKF(maxEpoch,kEnd,P0,Q,R,yMeas,iP,NN_UKF)
% Algorithm 3.1 from the Ref. paper
% Augmented state = [state; process_noise; measurement_noise]; 
x0 = nn2param(NN_UKF); nx = length(x0); ny = size(NN_UKF.B{end},1);
procMean = zeros(nx,1); measMean = zeros(ny,1);
xAug0 = [x0;procMean;measMean];
PAug0 = [P0,                    zeros(nx,nx+ny);
         zeros(nx,nx+ny)', diag([diag(Q);diag(R)])];
xAug_prev = xAug0; 
PAug_prev = PAug0; 
yUKF=zeros(kEnd,maxEpoch);
eta = 1;
for iEp = 1:maxEpoch
    eta = max(0.5*eta,0.1);
    fprintf('UKF: Epoch = %d.\n',iEp);
    for k = 1:kEnd
        % clc;
        
        
        [xAugSP,Wm,Wc] = getSigmaPts(xAug_prev,PAug_prev,(1e-2),0,2); % sigma pts of augemented state
        nPt = size(xAugSP,2); % no. of sigma pts
        xSP = xAugSP(1:nx,:); % sigma points of the actual state
        procSP =  xAugSP(nx+1:nx+nx,:); % sigma points of the process noise
        measSP = xAugSP(nx+nx+1:end,:); % sigma points of the measurement noise
        
%         [xSP,Wm,Wc] = getSigmaPts(x_prev,P_prev,1e-3,0,2);
%         nPt = size(xSP,2); % no. of sigma pts
        
        % UKF Time Update
        xSP_pr = xSP +  procSP; % parameter dynamics with identity state transition matrix
%         xSP_pr = xSP + mvnrnd(procMean,Q,nPt)';
%         xSP_pr = xSP;

        x_pr = xSP_pr*Wm'; % a priori state estimate = weighted sum of prior sigma pts
        tmpP1 = xSP_pr - x_pr; P_pr = zeros(nx,nx);
        for isp = 1:nPt
            P_pr = P_pr + Wc(isp)*(tmpP1(:,isp)*tmpP1(:,isp)');
        end
        
        ySP_pr = zeros(ny,nPt); % prior sigma pts of o/p
        for isp = 1:nPt
            NN_UKF = param2nn(NN_UKF,xSP_pr(:,isp)); % Update parameters of the NN
            ySP_pr(:,isp) = measModel(NN_UKF,iP(k,:)') + measSP(:,isp);
%             ySP_pr(:,isp) = measModel(NN_UKF,iP(k,:)') + mvnrnd(measMean,R)';
%         ySP_pr(:,isp) = measModel(NN_UKF,iP(k,:)');
        end
        y_pr = ySP_pr*Wm'; % a priori o/p estimate
        
        % UKF Measurement Update
        Pyy = zeros(ny,ny); Pxy = zeros(nx,ny);
        tmpP2 = ySP_pr - y_pr;
        for isp = 1:nPt
            Pyy = Pyy +  Wc(isp)*(tmpP2(:,isp)*tmpP2(:,isp)');
            Pxy = Pxy +  Wc(isp)*(tmpP1(:,isp)*tmpP2(:,isp)');
        end

        KK = Pxy/Pyy;
        res = yMeas(k,1) - y_pr;
        x_pst = x_pr + KK*(res);
        P_pst = P_pr - KK*Pyy*KK';
        
        x_prev = x_pst;
        P_prev = P_pst;
         
        xAug_prev = [x_pst;procMean;measMean];
        PAug_prev = [P_pst,                    zeros(nx,nx+ny);
            zeros(nx,nx+ny)', diag([eta*diag(Q);diag(R)])];
    end % k
    % Evaluate o/p of NN using a posteriori parameters estimates
    NN_UKF = param2nn(NN_UKF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yUKF(k,iEp) = measModel(NN_UKF,iP(k,:)');
    end
end % epoch

% MATLAB UKF
% hmeas = @(x,ip)(measModel(param2nn(NN_OTF,x),ip));
% fstate = @(x,ip)(x);
% myUKF = unscentedKalmanFilter(fstate,hmeas,x0); 
% myUKF.MeasurementNoise = R;
% myUKF.ProcessNoise = Q;
% for iEp = 1:maxEpoch
%      fprintf('UKF: Epoch = %d',iEp);
%     for k=1:kEnd
%         predict(myUKF, iP(k,:)');
%         [x_pst, P_pst] = correct(myUKF,yMeas(k,1),iP(k,:)');
%     end
%     NN_UKF = param2nn(NN_UKF,x_pst); % Update parameters of the NN
%     for k = 1:kEnd
%         yUKF(k,iEp) = measModel(NN_UKF,iP(k,:)');
%     end
%     RMSE_UKF(1,iEp) = norm(yUKF(:,iEp) - y1)/sqrt(kEnd);
%     fprintf(' RMSE = %d.\n', RMSE_UKF(1,iEp));
%     if  RMSE_UKF(1,iEp) < tol
%         Ep_UKF = iEp;
%         break;
%         % x_prev_Ep = x_pst;
%     end
% end

end % function