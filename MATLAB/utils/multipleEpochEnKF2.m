function [yEnKF, NN_EnKF]  = multipleEpochEnKF2(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NN_EnKF)
x0 = getParams(NN_EnKF);
xInitSamples = mvnrnd(x0,P0,nSample)'; % Each column is a sample
hmeas = @(x,ip)(nntool_eval(NN_EnKF,x,ip));
fstate = @(x,ip)(x);
xSamples_prev = xInitSamples;
yEnKF=zeros(kEnd,maxEpoch);
for iEp = 1:maxEpoch
    % clc;
    fprintf('EnKF: Epoch = %d.\n',iEp);
    [xSamples_pst,W_pst] = EnKF(xSamples_prev,fstate,hmeas,yMeas,Q,R,iP);
    x_pst = xSamples_pst*W_pst';
    xSamples_prev = xSamples_pst;
    
    NN_EnKF = updateParams(NN_EnKF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yEnKF(k,iEp) = NN_EnKF(iP(k,:)');
    end
end % epoch
end % function