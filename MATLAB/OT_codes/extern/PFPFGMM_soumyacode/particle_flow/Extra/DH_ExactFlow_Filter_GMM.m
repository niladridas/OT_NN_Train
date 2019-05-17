function [output] = DH_ExactFlow_Filter_GMM(ps,z)
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

for tt = 1:T
    %% redraw particles
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
            propparams_no_noise.stateCovarianceSR = zeros(size(ps.propparams.stateCovarianceSR));
        case 'NonStationaryGrowth'
            propparams_no_noise.stateCovarianceSR = zeros(size(ps.propparams.stateCovarianceSR));
        case 'InteractingNonStationaryGrowth'
            propparams_no_noise.stateCovarianceSR = zeros(size(ps.propparams.stateCovarianceSR));
        otherwise
            error('The example name does not match the record');
    end
    
    %% prediction step
    
    if tt == 1
        N = 1;
        [vg.M_prior,vg.PP] = ekf_predict1(vg.M,vg.PU,ps.propparams.dpropagatefcn_dx,ps.propparams.stateCovariance,@propparams_no_noise.propagatefcn,[],propparams_no_noise);
        %[vg.M_prior,vg.PP] = ukf_predict1(vg.M,vg.PU,@propparams_no_noise.propagatefcn,ps.propparams.stateCovariance,propparams_no_noise);
        vg.PP = (vg.PP + vg.PP')/2;
        vg.PP = cov_regularize(vg.PP);
        vg.prior_predictive = gmdistribution(vg.M_prior',vg.PP);
    else
        N = vg.posterior.NumComponents;
        ComponentProportion = vg.posterior.ComponentProportion;
        priorPredMean = zeros(size(vg.M,1),N);
        priorPredCov = zeros(size(vg.M,1),size(vg.M,1),N);
        for j= 1:N
                         [priorPredMean(:,j),priorPredCov(:,:,j)] = ...
        ekf_predict1(vg.posterior.mu(j,:)',vg.posterior.Sigma(:,:,j),ps.propparams.dpropagatefcn_dx,ps.propparams.stateCovariance,@propparams_no_noise.propagatefcn,[],propparams_no_noise);
%             [priorPredMean(:,j),priorPredCov(:,:,j)] = ukf_predict1(vg.posterior.mu(j,:)',vg.posterior.Sigma(:,:,j),...
%                 @propparams_no_noise.propagatefcn,ps.propparams.stateCovariance,propparams_no_noise);
            priorPredCov(:,:,j) = (priorPredCov(:,:,j)+priorPredCov(:,:,j)')/2;
            priorPredCov(:,:,j) = cov_regularize(priorPredCov(:,:,j));
        end
        vg.prior_predictive = gmdistribution(priorPredMean',priorPredCov,ComponentProportion);
    end
    
    %% propagate particle and Step through the lambda
    if tt==1
        N =1;
        vg.xp = ps.propparams.propagatefcn(vg.xp,ps.propparams);
        %         vg.M_prior = mean(vg.xp,2);
        %         vg.PP = cov_regularize(cov(vg.xp'));
        %         vg.prior_predictive = gmdistribution(vg.M_prior',vg.PP);
    else
        N = vg.posterior.NumComponents;
        %         ComponentProportion = vg.posterior.ComponentProportion;
        %         priorPredMean = zeros(size(vg.M,1),N);
        %         priorPredCov = zeros(size(vg.M,1),size(vg.M,1),N);
        for j = 1:N
            vg.xp(:,:,j) = ps.propparams.propagatefcn(vg.xp(:,:,j),ps.propparams);
            %         priorPredMean(:,j) = mean(vg.xp(:,:,j),2);
            %         priorPredCov(:,:,j) = cov_regularize(cov(vg.xp(:,:,j)'));
        end
        %        vg.prior_predictive = gmdistribution(priorPredMean',priorPredCov,ComponentProportion);
    end
    
    likelihood = ps.likeparams.gmmNoiseModel;
    NumLikeComp = likelihood.NumComponents;
    for k = 1: NumLikeComp
        for l = 1:N
            vg.new_particles(:,:,(k-1)*N+l) = vg.xp(:,:,l);
            
            mu_0 = vg.prior_predictive.mu(l,:)';
            P = vg.prior_predictive.Sigma(:,:,l);
            R = likelihood.Sigma(:,:,k);
            zeta = likelihood.mu(k,:)';
            
            lambda_prev = 0;
            for lambda = ps.setup.lambda_range
                step_size = lambda-lambda_prev;
                
                % Calculate the slopes for moving the particles
                
                if ps.setup.LEDH
                    slope = calculateSlopeComponentwise_LEDH(z(:,tt),vg.new_particles(:,:,(k-1)*N+l),mu_0,P,R,zeta,ps,lambda);
                else
                    slope = calculateSlopeComponentwise_EDH(z(:,tt),vg.new_particles(:,:,(k-1)*N+l),mu_0,P,R,zeta,ps,lambda);
                end
                
                vg.new_particles(:,:,(k-1)*N+l) = vg.new_particles(:,:,(k-1)*N+l) + step_size*slope;  % euler update of particles
                
                lambda_prev = lambda;
            end
        end
    end
    %% EKF update and particle flow based estimation, copying paricle means to ekf updated means
    if tt==1
        N =1;
    else
        N = vg.posterior.NumComponents;
    end
    
    posteriorMean = zeros(size(vg.M,1),NumLikeComp*N);
    posteriorCov =  zeros(size(vg.M,1),size(vg.M,1),NumLikeComp*N);
    posteriorProportion = zeros(1,NumLikeComp*N);
    for k = 1: NumLikeComp
        for l = 1:N
                        [~, posteriorCov(:,:,(k-1)*N+l)] = ...
                            ekf_update1(vg.prior_predictive.mu(l,:)',vg.prior_predictive.Sigma(:,:,l),z(:,tt)-likelihood.mu(k,:)',...
                          ps.likeparams.dh_dx_func,likelihood.Sigma(:,:,k),ps.likeparams.h_func,[],ps.likeparams);
%             [~, posteriorCov(:,:,(k-1)*N+l)] = ukf_update1(vg.prior_predictive.mu(l,:)',vg.prior_predictive.Sigma(:,:,l),z(:,tt)-likelihood.mu(k,:)',...
%                 ps.likeparams.h_func,likelihood.Sigma(:,:,k),ps.likeparams);
            posteriorCov(:,:,(k-1)*N+l) = (posteriorCov(:,:,(k-1)*N+l)+posteriorCov(:,:,(k-1)*N+l)')/2;
            posteriorCov(:,:,(k-1)*N+l) = cov_regularize(posteriorCov(:,:,(k-1)*N+l));
            %            posteriorCov(:,:,(k-1)*N+l) =  cov_regularize(cov(vg.new_particles(:,:,(k-1)*N+l)'));
            posteriorMean(:,(k-1)*N+l) = mean(vg.new_particles(:,:,(k-1)*N+l),2);
            
            
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
    
    output.x_est(:,tt) = sum(posteriorMean.*repmat(posteriorProportion,size(posteriorMean,1),1),2);
    proportion_threshold = 0;
    posteriorMean = posteriorMean(:,posteriorProportion>proportion_threshold);
    posteriorCov = posteriorCov(:,:,posteriorProportion>proportion_threshold);
    vg.new_particles = vg.new_particles(:,:,posteriorProportion>proportion_threshold);
    posteriorProportion = posteriorProportion(posteriorProportion>proportion_threshold);
    posteriorProportion = posteriorProportion/sum(posteriorProportion);
    [~,indices] = sort( posteriorProportion,'descend');
    posteriorMean = posteriorMean(:,indices);
    posteriorCov = posteriorCov(:,:,indices);
    vg.new_particles = vg.new_particles(:,:,indices);
    posteriorProportion = posteriorProportion(indices);
    [posteriorMean,posteriorCov,posteriorProportion,idx] = resample_GMM_Comp(posteriorMean,...
        posteriorCov,posteriorProportion,min(NumLikeComp,length(posteriorProportion)),0.05);
    %     [posteriorMean,posteriorCov,posteriorProportion,idx] = prune_RKLD_GMM_Comp(posteriorMean,...
    %     posteriorCov,posteriorProportion,min(NumLikeComp,length(posteriorProportion)));
    vg.posterior = gmdistribution(posteriorMean',posteriorCov,posteriorProportion);
    vg.xp = vg.new_particles(:,:,idx);
end

output.x = ps.x;
output.execution_time = toc;
NumComp = ps.likeparams.gmmNoiseModel.NumComponents;
if NumComp ==1 && ps.setup.LEDH
    alg_name = 'LEDH';
elseif ps.setup.LEDH
    alg_name = 'PF_GMM_LEDH';
elseif NumComp ==1
    alg_name = 'EDH';
else
    alg_name = 'PF_GMM_EDH';
end
calculateErrors(output,ps,alg_name);
end