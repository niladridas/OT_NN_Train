function [Otfilterdata] = OTfilterandprop(initsamples,obs,simparams,filterparams)
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
        %% OT Filter update
        tmp_prior = reshape(prior_samples(:,:,i),nSamp,nx);
        tic; 
        W_post = zeros(1,nSamp); 
        for k=1:nSamp
            y = obs(i,:)';
            W_post(1,k)  = simparams.likfun(y,tmp_prior(k,:)');
        end
        W_post = W_post./(sum(W_post));
        W0 = ones(1,nSamp)/nSamp; % Prior is equally weighted
        OT_samplesX = eq_wsamples(prior_samples(:,:,i)',W0,W_post); toc;
        run_time(1,i) = toc;
        post_samples(:,:,i) = OT_samplesX';
        if filterparams.MAP_Flag == 1
            W_post = zeros(1,nSamp); 
            for l=1:nSamp
                y = obs(i,:)';
                W_post(1,l)  = simparams.likfun(y,OT_samplesX(:,l));
            end
            [~,I] = max(W_post);      
            MAP_estimate(i,:) = OT_samplesX(:,I)';
        end
        initsamples = post_samples(:,:,i);
    end
    Otfilterdata.prior_samples = prior_samples;
    Otfilterdata.post_samples = post_samples;
    Otfilterdata.run_time = run_time;
    if filterparams.MAP_Flag == 1
        Otfilterdata.MAP_estimate = MAP_estimate;
    end
end