function [output] = EKF_UKF_GMM(ps,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The execution of the extended Kalman filter.
%
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

for tt = 1:T
    ps.propparams.time_step = tt;
    
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
        for i = 1:NumDynComp
            if ps.setup.EKF
                [priorPredMean(:,i),priorPredCov(:,:,i)] = ekf_predict1(vg.M,vg.PU,ps.propparams.dpropagatefcn_dx,dynamic_model.Sigma(:,:,i),@propparams_no_noise.propagatefcn,[],propparams_no_noise);
                priorPredMean(:,i) = priorPredMean(:,i) + dynamic_model.mu(i,:)';
            else
                [priorPredMean(:,i),priorPredCov(:,:,i)] = ukf_predict1(vg.M,vg.PU,@propparams_no_noise.propagatefcn,dynamic_model.Sigma(:,:,i),propparams_no_noise);
                priorPredMean(:,i) = priorPredMean(:,i) + dynamic_model.mu(i,:)';
            end
            priorPredCov(:,:,i) = (priorPredCov(:,:,i) + priorPredCov(:,:,i))/2;
            priorPredCov(:,:,i) = cov_regularize(priorPredCov(:,:,i));
        end
    else
        N = vg.posterior.NumComponents;
        priorPredProportion = zeros(1,N*NumDynComp);   %%vg.posterior.ComponentProportion;
        priorPredMean = zeros(size(vg.M,1),N*NumDynComp);
        priorPredCov = zeros(size(vg.M,1),size(vg.M,1),N*NumDynComp);
        for i = 1:NumDynComp
            for j= 1:N
                if ps.setup.EKF
                    [priorPredMean(:,(i-1)*N+j),priorPredCov(:,:,(i-1)*N+j)] = ...
                        ekf_predict1(vg.posterior.mu(j,:)',vg.posterior.Sigma(:,:,j),ps.propparams.dpropagatefcn_dx,dynamic_model.Sigma(:,:,i),@propparams_no_noise.propagatefcn,[],propparams_no_noise);
                    priorPredMean(:,(i-1)*N+j) = priorPredMean(:,(i-1)*N+j) + dynamic_model.mu(i,:)';
                else
                    [priorPredMean(:,(i-1)*N+j),priorPredCov(:,:,(i-1)*N+j)] = ...
                        ukf_predict1(vg.posterior.mu(j,:)',vg.posterior.Sigma(:,:,j),@propparams_no_noise.propagatefcn,dynamic_model.Sigma(:,:,i),propparams_no_noise);
                    priorPredMean(:,(i-1)*N+j) = priorPredMean(:,(i-1)*N+j) + dynamic_model.mu(i,:)';
                end
                priorPredCov(:,:,(i-1)*N+j) = (priorPredCov(:,:,(i-1)*N+j)+priorPredCov(:,:,(i-1)*N+j)')/2;
                priorPredCov(:,:,(i-1)*N+j) = cov_regularize(priorPredCov(:,:,(i-1)*N+j));
                priorPredProportion((i-1)*N+j) = dynamic_model.ComponentProportion(i)*vg.posterior.ComponentProportion(j);
            end
        end
    end
    
    [~,indices] = sort(priorPredProportion,'descend');
    priorPredMean = priorPredMean(:,indices);
    priorPredCov = priorPredCov(:,:,indices);
    priorPredProportion =  priorPredProportion(indices);
        [priorPredMean,priorPredCov, priorPredProportion,~] = resample_GMM_Comp(priorPredMean,...
        priorPredCov,priorPredProportion,max(N,NumDynComp),0.05);
    
    vg.prior_predictive = gmdistribution(priorPredMean',priorPredCov,priorPredProportion);
    %% update step
    likelihood = ps.likeparams.gmmNoiseModel;
    NumLikeComp = likelihood.NumComponents;
    
    N = vg.prior_predictive.NumComponents;
    
    posteriorMean = zeros(size(vg.M,1),NumLikeComp*N);
    posteriorCov =  zeros(size(vg.M,1),size(vg.M,1),NumLikeComp*N);
    posteriorProportion = zeros(1,NumLikeComp*N);
    
    for k = 1: NumLikeComp
        for l = 1:N
            if ps.setup.EKF
                [posteriorMean(:,(k-1)*N+l), posteriorCov(:,:,(k-1)*N+l)] = ...
                    ekf_update1(vg.prior_predictive.mu(l,:)',vg.prior_predictive.Sigma(:,:,l),z(:,tt)-likelihood.mu(k,:)',ps.likeparams.dh_dx_func,likelihood.Sigma(:,:,k),ps.likeparams.h_func,[],ps.likeparams);
            else
                [posteriorMean(:,(k-1)*N+l), posteriorCov(:,:,(k-1)*N+l)] = ...
                    ukf_update1(vg.prior_predictive.mu(l,:)',vg.prior_predictive.Sigma(:,:,l),z(:,tt)-likelihood.mu(k,:)',ps.likeparams.h_func,likelihood.Sigma(:,:,k),ps.likeparams);
            end
            posteriorCov(:,:,(k-1)*N+l) = (posteriorCov(:,:,(k-1)*N+l)+posteriorCov(:,:,(k-1)*N+l)')/2;
            posteriorCov(:,:,(k-1)*N+l) = cov_regularize(posteriorCov(:,:,(k-1)*N+l));
            H = ps.likeparams.dh_dx_func(vg.prior_predictive.mu(l,:)',ps.likeparams);
            h = ps.likeparams.h_func(vg.prior_predictive.mu(l,:)',ps.likeparams);
            cov_kl = H*vg.prior_predictive.Sigma(:,:,l)*H' + likelihood.Sigma(:,:,k);
            cov_kl = (cov_kl + cov_kl')/2;
            posteriorProportion((k-1)*N+l) = vg.prior_predictive.ComponentProportion(l)*...
                likelihood.ComponentProportion(k)*...
                mvnpdf(z(:,tt)',h'+likelihood.mu(k,:),cov_regularize(cov_kl));
        end
    end
    if isinf(sum(posteriorProportion))||isnan(sum(posteriorProportion))||sum(posteriorProportion)<=0
        posteriorProportion = ones(1,NumLikeComp*N)/(NumLikeComp*N);
    else
        posteriorProportion = posteriorProportion/sum(posteriorProportion);
    end
    %% this is resampling of components
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
end

output.x = ps.x;
output.execution_time = toc;
if ps.setup.EKF
    alg_name = 'EKF_GMM';
else
    alg_name = 'UKF_GMM';
end
calculateErrors(output,ps,alg_name);
end