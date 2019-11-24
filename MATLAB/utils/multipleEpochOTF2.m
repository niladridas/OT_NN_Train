function [yOTF,NN_OTF] = multipleEpochOTF2(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NN_OTF)
x0 = getParams(NN_OTF);
xInitSamples = mvnrnd(x0,P0,nSample)'; % Each column is a sample
hmeas = @(x,ip)(nntool_eval(NN_OTF,x,ip));
fstate = @(x,ip)(x);
xSamples_prev = xInitSamples;
yOTF=zeros(kEnd,maxEpoch);
for iEp = 1:maxEpoch
    % clc;
    fprintf('  Epoch = %d.\n',iEp);
    [xSamples_pst,W_pst] = OTF(xSamples_prev,fstate,hmeas,yMeas,Q,R,iP);
    x_pst = xSamples_pst*W_pst';
    xSamples_prev = xSamples_pst;
    
    NN_OTF = updateParams(NN_OTF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yOTF(k,iEp) = NN_OTF(iP(k,:)');
    end
end % epoch
end % function