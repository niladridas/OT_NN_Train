function yOTF = multipleEpochOTF(maxEpoch,kEnd,nSample,P0,Q,R,yMeas,iP,NN_OTF)
x0 = nn2param(NN_OTF);
xInitSamples = mvnrnd(x0,P0,nSample)'; % Each column is a sample
hmeas = @(x,ip)(measModel(param2nn(NN_OTF,x),ip));
fstate = @(x,ip)(x);
xSamples_prev = xInitSamples;
yOTF=zeros(kEnd,maxEpoch);
eta = 1;
for iEp = 1:maxEpoch
    eta = max(0.5*eta,0.1);
    % clc;
    fprintf('OTF:  Epoch = %d.\n',iEp);
    [xSamples_pst,W_pst] = OTF(xSamples_prev,fstate,hmeas,yMeas,eta*Q,R,iP,1);
    x_pst = xSamples_pst*W_pst';
    xSamples_prev = xSamples_pst;
    
    NN_OTF = param2nn(NN_OTF,x_pst); % Update parameters of the NN
    for k = 1:kEnd
        yOTF(k,iEp) = measModel(NN_OTF,iP(k,:)');
    end
end % epoch
end % function