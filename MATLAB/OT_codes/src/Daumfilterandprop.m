function [Daumfilterdata] = Daumfilterandprop(initsamples,obs,simparams,filterparams)
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
        %% Daum Filter update
        tmp_prior = reshape(prior_samples(:,:,i),nSamp,nx);
        P = cov(tmp_prior);% Prior sample covariance
        mu_0 = mean(tmp_prior)';% Prior sample mean
        R = simparams.meas_sig;
        y = obs(i,:)';
        flowdyn = @(t,x) EDH(t,x,y,simparams.H,mu_0,P,R);
        Daum_samplesX = zeros(4,nSamp);
        tic
        parfor k = 1:nSamp
            [~,xf] = ode23(flowdyn,[0 1],tmp_prior(k,:)');
            Daum_samplesX(:,k) = xf(end,:)';
        end
        toc
        run_time(1,i) = toc;
        post_samples(:,:,i) = Daum_samplesX';
        if filterparams.MAP_Flag == 1
            W_post = zeros(1,nSamp); 
            for l=1:nSamp
                y = obs(i,:)';
                W_post(1,l)  = simparams.likfun(y,Daum_samplesX(:,l));
            end
            [~,I] = max(W_post);      
            MAP_estimate(i,:) = Daum_samplesX(:,I)';
        end
        initsamples = post_samples(:,:,i);
    end
    Daumfilterdata.prior_samples = prior_samples;
    Daumfilterdata.post_samples = post_samples;
    Daumfilterdata.run_time = run_time; 
    if filterparams.MAP_Flag == 1
        Daumfilterdata.MAP_estimate = MAP_estimate;   
    end
end