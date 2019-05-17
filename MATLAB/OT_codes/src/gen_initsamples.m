function [initsamples] = gen_initsamples(filterparams,simparams,nSamp)
    X = mvnrnd(filterparams.mu,filterparams.Sig,nSamp);
    Rd = ones(nSamp,1)*simparams.rd;
    Thd = ones(nSamp,1)*simparams.thd;
    initsamples = [X(:,1) Rd X(:,2) Thd];
    % The initial samples are themselves generated after propagating the
    % samples for few time points to get the banana shape.
    for i = 1:nSamp
        [~,x1] = ode45(simparams.dyn,[0 simparams.burnperiod],initsamples(i,:)',[]);
        initsamples(i,:) = x1(end,:)';
    end
end