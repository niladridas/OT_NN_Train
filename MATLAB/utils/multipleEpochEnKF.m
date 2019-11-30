function yEnKF = multipleEpochEnKF(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NN_EnKF)
x0 = nn2param(NN_EnKF);
xInitSamples = mvnrnd(x0,P0,nSample)'; % Each column is a sample
hmeas = @(x,ip)(measModel(param2nn(NN_EnKF,x),ip));
fstate = @(x,ip)(x);
xSamples_prev = xInitSamples;
yEnKF=zeros(kEnd,maxEpoch);
eta = 1;
for iEp = 1:maxEpoch
    eta = 0; % max(0.5*eta,0.1);
    % clc;
    fprintf('EnKF: Epoch = %d.\n',iEp);
    [xSamples_pst,W_pst] = EnKF(xSamples_prev,fstate,hmeas,yMeas,eta*Q,R,iP,1);
    x_pst = xSamples_pst*W_pst';
    xSamples_prev = xSamples_pst;
    
    NN_EnKF = param2nn(NN_EnKF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yEnKF(k,iEp) = measModel(NN_EnKF,iP(k,:)');
    end
end % epoch
end % function