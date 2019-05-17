function [output] = GSPF(ps,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% ps: a struct that contains model parameters.
% z: a measurement_dim x T matrix
%
% Output:
% output: a struct that contains the filter outputs, including the particle
% estimate, true state, execution time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
T = ps.setup.T; % number of time steps

[vg,output] = initializationFilter(ps);

nParticle = ps.setup.nParticle;

for tt = 1:T
    %% redraw particles componentwise
    if tt~=1
        N = vg.posterior.NumComponents;
        for j = 1:N
            vg.xp(:,:,j) = mvnrnd(vg.posterior.mu(j,:)',vg.posterior.Sigma(:,:,j),ps.setup.nParticle)';
        end
    end
    
    ps.propparams.time_step = tt;
    
    % Propagate the particles one step componentwise
    propparams_no_noise = ps.propparams;
    switch ps.setup.example_name
        case 'Acoustic'
            propparams_no_noise.stateCovarianceSR = zeros(size(ps.propparams.stateCovarianceSR));
        case 'Septier16'
            propparams_no_noise.stateMean = ps.propparams.stateMean;
            propparams_no_noise.stateCovariance = zeros(size(ps.propparams.stateCovariance));
            propparams_no_noise.prop_type = 'gaussian';
        case 'NonStationaryGrowth'
            propparams_no_noise.stateCovarianceSR = zeros(size(ps.propparams.stateCovarianceSR));
        case 'InteractingNonStationaryGrowth'
            propparams_no_noise.stateMean = ps.propparams.stateMean;
            propparams_no_noise.stateCovariance = zeros(size(ps.propparams.stateCovariance));
            propparams_no_noise.prop_type = 'gaussian';
        otherwise
            error('The example name does not match the record');
    end
    %% prediction step
    dynamic_model = ps.propparams.gmmNoiseModel_dyn;
    NumDynComp = dynamic_model.NumComponents;
    if tt == 1
        N = 1;
        priorPredProportion = dynamic_model.ComponentProportion;
        priorPredMean = zeros(size(vg.M,1),NumDynComp);
        priorPredCov = zeros(size(vg.M,1),size(vg.M,1),NumDynComp);
        vg.xp_prop = zeros(size(vg.xp));
        for i = 1:NumDynComp
            vg.xp_prop(:,:,i) = ps.propparams.propagatefcn(vg.xp,propparams_no_noise) + mvnrnd(dynamic_model.mu(i,:)',dynamic_model.Sigma(:,:,i),nParticle)';
            priorPredMean(:,i) = mean(vg.xp_prop(:,:,i),2);
            priorPredCov(:,:,i) = cov(vg.xp_prop(:,:,i)');
            priorPredCov(:,:,i) = (priorPredCov(:,:,i) + priorPredCov(:,:,i))/2;
            priorPredCov(:,:,i) = cov_regularize(priorPredCov(:,:,i));
        end
    else
        N = vg.posterior.NumComponents;
        priorPredProportion = zeros(1,N*NumDynComp);
        priorPredMean = zeros(size(vg.M,1),N*NumDynComp);
        priorPredCov = zeros(size(vg.M,1),size(vg.M,1),N*NumDynComp);
        vg.xp_prop = zeros(size(vg.xp,1),nParticle,N*NumDynComp);
        for i = 1:NumDynComp
            for j= 1:N
                vg.xp_prop(:,:,(i-1)*N+j) = ps.propparams.propagatefcn(vg.xp(:,:,j),propparams_no_noise) + mvnrnd(dynamic_model.mu(i,:)',dynamic_model.Sigma(:,:,i),nParticle)';
                priorPredMean(:,(i-1)*N+j) = mean(vg.xp_prop(:,:,(i-1)*N+j),2);
                priorPredCov(:,:,(i-1)*N+j) = cov(vg.xp_prop(:,:,(i-1)*N+j)');
                priorPredCov(:,:,(i-1)*N+j) = (priorPredCov(:,:,(i-1)*N+j)+priorPredCov(:,:,(i-1)*N+j)')/2;
                priorPredCov(:,:,(i-1)*N+j) = cov_regularize(priorPredCov(:,:,(i-1)*N+j));
                priorPredProportion((i-1)*N+j) = dynamic_model.ComponentProportion(i)*vg.posterior.ComponentProportion(j);
            end
        end
    end
    
    [~,indices] = sort(priorPredProportion,'descend');
    priorPredMean = priorPredMean(:,indices);
    priorPredCov = priorPredCov(:,:,indices);
    vg.xp_prop = vg.xp_prop(:,:,indices);
    priorPredProportion =  priorPredProportion(indices);
    [priorPredMean,priorPredCov, priorPredProportion,idx] = resample_GMM_Comp(priorPredMean,...
        priorPredCov,priorPredProportion,max(N,NumDynComp),0.05);
    vg.prior_predictive = gmdistribution(priorPredMean',priorPredCov,priorPredProportion);
    vg.xp = vg.xp_prop(:,:,idx);
    %% measurement update
    likelihood = ps.likeparams.gmmNoiseModel;
    NumLikeComp = likelihood.NumComponents;
    
    N = vg.prior_predictive.NumComponents;
    
    posteriorMean = zeros(size(vg.M,1),NumLikeComp*N);
    posteriorCov =  zeros(size(vg.M,1),size(vg.M,1),NumLikeComp*N);
    posteriorProportion = zeros(1,NumLikeComp*N);
    logWeights = zeros(NumLikeComp*N,nParticle);
    Weights = zeros(NumLikeComp*N,nParticle);
    for k = 1: NumLikeComp
        for l = 1:N
            
            R_k = likelihood.Sigma(:,:,k);
            zeta_k = likelihood.mu(k,:)';
            
            h_l = ps.likeparams.h_func(vg.xp(:,:,l),ps.likeparams);
            
            logWeights((k-1)*N+l,:) = loggausspdf(h_l,(z(:,tt)-zeta_k),R_k);
            
            Weights((k-1)*N+l,:) = exp(logWeights((k-1)*N+l,:));
            
            posteriorProportion((k-1)*N+l) = vg.prior_predictive.ComponentProportion(l)*...
                likelihood.ComponentProportion(k)* mean(Weights((k-1)*N+l,:));
            
            if any(isinf(Weights((k-1)*N+l,:)))||any(isnan(Weights((k-1)*N+l,:)))||any(~isreal(Weights((k-1)*N+l,:)))
                Weights((k-1)*N+l,:) = ones(1,nParticle)/nParticle;
                logWeights((k-1)*N+l,:) = zeros(1,nParticle);
            else
                Weights((k-1)*N+l,:) = Weights((k-1)*N+l,:)/sum(Weights((k-1)*N+l,:));
            end
            
            posteriorMean(:,(k-1)*N+l) = particle_estimate(logWeights((k-1)*N+l,:),vg.xp(:,:,l));
            posteriorCov(:,:,(k-1)*N+l) = (vg.xp(:,:,l) - ...
                repmat(posteriorMean(:,(k-1)*N+l),1,nParticle))*diag(Weights((k-1)*N+l,:))*...
                (vg.xp(:,:,l) - repmat(posteriorMean(:,(k-1)*N+l),1,nParticle))';
            posteriorCov(:,:,(k-1)*N+l) = (posteriorCov(:,:,(k-1)*N+l)+ posteriorCov(:,:,(k-1)*N+l)')/2;
            posteriorCov(:,:,(k-1)*N+l) = cov_regularize(posteriorCov(:,:,(k-1)*N+l));
            
        end
    end
    
    if isinf(sum(posteriorProportion))||isnan(sum(posteriorProportion))||sum(posteriorProportion)<=0
        posteriorProportion = ones(1,NumLikeComp*N)/(NumLikeComp*N);
    else
        posteriorProportion = posteriorProportion/sum(posteriorProportion);
    end
    
    output.x_est(:,tt) = sum(posteriorMean.*repmat(posteriorProportion,size(posteriorMean,1),1),2);
    proportion_threshold = 0;
    posteriorMean = posteriorMean(:,posteriorProportion>proportion_threshold);
    posteriorCov = posteriorCov(:,:,posteriorProportion>proportion_threshold);
    posteriorProportion = posteriorProportion(posteriorProportion>proportion_threshold);
    posteriorProportion = posteriorProportion/sum(posteriorProportion);
    [~,indices] = sort( posteriorProportion,'descend');
    posteriorMean = posteriorMean(:,indices);
    posteriorCov = posteriorCov(:,:,indices);
    posteriorProportion = posteriorProportion(indices);
    [posteriorMean,posteriorCov,posteriorProportion,~] = resample_GMM_Comp(posteriorMean,...
     posteriorCov,posteriorProportion,min(NumLikeComp,length(posteriorProportion)),0.05);
    vg.posterior = gmdistribution(posteriorMean',posteriorCov,posteriorProportion);
    vg.xp = []; %no need to store the particles from current step
end

output.x = ps.x;
output.execution_time = toc;
alg_name = 'GSPF';
calculateErrors(output,ps,alg_name);
end