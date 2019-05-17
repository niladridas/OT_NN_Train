function [output] = PFPF_GMM(ps,z)
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
nParticle = ps.setup.nParticle;

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
    
    %% propagate particle and individual EKF/UKF prediction
    vg.xp_prop_deterministic = ps.propparams.propagatefcn(vg.xp,propparams_no_noise);
    vg.xp_prop = zeros(size(vg.xp,1),nParticle);
    vg.mu = zeros(size(vg.xp,1),nParticle);
    vg.PP = zeros(size(vg.xp,1),size(vg.xp,1),nParticle);
    vg.dynCompMean = zeros(size(vg.xp,1),nParticle);
    vg.dynCompCov = zeros(size(vg.xp,1),size(vg.xp,1),nParticle);
    vg.R = zeros(size(z,1),size(z,1),nParticle);
    vg.zeta = zeros(size(z,1),nParticle);
    
    if tt == 1
        vg.PU = repmat(vg.PU,1,1,nParticle);
    end
    
    dynamic_model = ps.propparams.gmmNoiseModel_dyn;
    likelihood = ps.likeparams.gmmNoiseModel;
    
    for i = 1:nParticle
        dyn_compID = resample(1,dynamic_model.ComponentProportion,'stratified');
        vg.dynCompMean(:,i) = dynamic_model.mu(dyn_compID,:)';
        vg.dynCompCov(:,:,i) = dynamic_model.Sigma(:,:,dyn_compID);
        vg.xp_prop(:,i) = ps.propparams.propagatefcn(vg.xp(:,i),propparams_no_noise) + mvnrnd(vg.dynCompMean(:,i)',vg.dynCompCov(:,:,i))';
        if ps.setup.EKF
            [vg.mu(:,i),vg.PP(:,:,i)] = ekf_predict1(vg.xp(:,i),vg.PU(:,:,i),ps.propparams.dpropagatefcn_dx,vg.dynCompCov(:,:,i),@propparams_no_noise.propagatefcn,[],propparams_no_noise);
            vg.mu(:,i) = vg.mu(:,i) + vg.dynCompMean(:,i);
        else
            [vg.mu(:,i),vg.PP(:,:,i)] = ukf_predict1(vg.xp(:,i),vg.PU(:,:,i),@propparams_no_noise.propagatefcn,vg.dynCompCov(:,:,i),propparams_no_noise);
            vg.mu(:,i) = vg.mu(:,i) + vg.dynCompMean(:,i);
        end
        vg.PP(:,:,i) = (vg.PP(:,:,i)+vg.PP(:,:,i)')/2;
        vg.PP(:,:,i) = cov_regularize(vg.PP(:,:,i));
        like_compID = resample(1,likelihood.ComponentProportion,'stratified');
        vg.R(:,:,i) =  likelihood.Sigma(:,:,like_compID);
        vg.zeta(:,i) =  likelihood.mu(like_compID,:)';
    end
    %% propagate particle and Step through the lambda
    vg.mu_0 = vg.xp_prop_deterministic + vg.dynCompMean;
    vg.new_particles = vg.xp_prop;
    vg.new_particles_deterministic = vg.xp_prop_deterministic + vg.dynCompMean;
    
    lambda_prev = 0;
    for lambda = ps.setup.lambda_range
        step_size = lambda-lambda_prev;
        
        % Calculate the slopes for moving the particles
        [slope_deterministic,slope,logJacobianDet]  = calculateSlope_PFPF_GMM...
            (z(:,tt),vg.new_particles_deterministic,vg.new_particles,vg.mu_0,vg.PP,vg.R,vg.zeta,ps,lambda,step_size);
        vg.logW = vg.logW + logJacobianDet;
        
        vg.new_particles_deterministic = vg.new_particles_deterministic + step_size*slope_deterministic;
        vg.new_particles = vg.new_particles + step_size*slope;  % euler update of particles
        lambda_prev = lambda;
    end
    
    % weight update and estimation and resampling
    logProcessDensity = log(mvnpdf((vg.new_particles-vg.dynCompMean)',vg.xp_prop_deterministic',vg.dynCompCov))';
    logProposalDensity = log(mvnpdf((vg.xp_prop-vg.dynCompMean)',vg.xp_prop_deterministic',vg.dynCompCov))';
    llh = log(mvnpdf(ps.likeparams.h_func(vg.new_particles,ps.likeparams)'+vg.zeta',z(:,tt)',vg.R))';
    vg.logW = vg.logW +  logProcessDensity + llh - logProposalDensity;
    
    wt = exp(vg.logW);
    wt = wt./sum(wt);
    if isnan(wt)
        wt = ones(1,nParticle)/nParticle;
        vg.logW = log(wt);
        output.Neff(tt) = 1;
    else
        output.Neff(tt) = 1/sum(wt.^2);
    end
    
    output.x_est(:,tt) = particle_estimate(vg.logW,vg.new_particles);
    
    %% EKF/UKF update for individual particle copying paricles to ekf updated means
    vg.PU = zeros(size(vg.xp,1),size(vg.xp,1),nParticle);
    for i = 1:nParticle
        if ps.setup.EKF
            [~, vg.PU(:,:,i)] = ...
                ekf_update1(vg.mu(:,i),vg.PP(:,:,i),z(:,tt)-vg.zeta(:,i),ps.likeparams.dh_dx_func,vg.R(:,:,i),ps.likeparams.h_func,[],ps.likeparams);
        else
            [~, vg.PU(:,:,i)] = ...
                ukf_update1(vg.mu(:,i),vg.PP(:,:,i),z(:,tt)-vg.zeta(:,i),ps.likeparams.h_func,vg.R(:,:,i),ps.likeparams);
        end
        vg.PU(:,:,i) = (vg.PU(:,:,i)+vg.PU(:,:,i)')/2;
        vg.PU(:,:,i) = cov_regularize(vg.PU(:,:,i));
    end
    
    [vg.logW,vg.xp,vg.PU,~] = resample_PFPF_GMM(vg.logW,vg.new_particles,vg.PU,ps);
end
output.x = ps.x;
output.execution_time = toc;
alg_name = 'PFPF_GMM';
calculateErrors(output,ps,alg_name);
disp(['avg. ESS: ',num2str(mean(output.Neff))])
end