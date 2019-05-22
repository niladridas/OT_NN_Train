function [EnKFfilterdata] = EnKFfilterandprop(initsamples,obs,simparams,filterparams)
    % obs: observations Nt x (dim of observed states)
    % initsamples: nSamp x (dim of system states)
    Nt = simparams.Nt;
    [nSamp,nx] = size(initsamples);
    if size(obs,1)~= Nt
        error("Check Dimension");
    end
    prior_samples = zeros(nSamp,nx,Nt);
    post_samples = zeros(nSamp,nx,Nt);
    MAP_estimate = zeros(Nt,nx);
    run_time = zeros(1,Nt);
    for i = 1:Nt
        %% EnKF based propagation
        for j = 1:nSamp
            [~,x1] = ode45(simparams.dyn,[0 i*simparams.delt],initsamples(j,:)',[]);
            prior_samples(j,:,i) = x1(end,:);
            
        end
        %% EnKF Filter update
        tmp_prior = reshape(prior_samples(:,:,i),nSamp,nx);
        y = obs(i,:)';
        Y_predict = simparams.H*tmp_prior';
        tic; EnkF_samplesX = enkf_samples(tmp_prior', y, Y_predict,simparams.meas_sig); toc;
        run_time(1,i) = toc;
        post_samples(:,:,i) = EnkF_samplesX';
        if filterparams.MAP_Flag == 1
            W_post = zeros(1,nSamp); 
            for l=1:nSamp
                y = obs(i,:)';
                W_post(1,l)  = simparams.likfun(y,EnkF_samplesX(:,l));
            end
            [~,I] = max(W_post);      
            MAP_estimate(i,:) = EnkF_samplesX(:,I)';
        end
        initsamples = post_samples(:,:,i);
    end
    EnKFfilterdata.prior_samples = prior_samples;
    EnKFfilterdata.post_samples = post_samples;
    EnKFfilterdata.run_time = run_time;
    if filterparams.MAP_Flag == 1
        EnKFfilterdata.MAP_estimate = MAP_estimate;  
    end
end