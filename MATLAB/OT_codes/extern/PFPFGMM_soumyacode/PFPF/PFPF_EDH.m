function [output] = PFPF_EDH(ps,z)
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
    
    %% propagate particle and EKF/UKF prediction
    vg.xp_prop = ps.propparams.propagatefcn(vg.xp,ps.propparams);
    vg.xp_prop_deterministic = ps.propparams.propagatefcn(vg.xp,propparams_no_noise);
    
    if ps.setup.EKF
        [vg.mu,vg.PP] = ekf_predict1(vg.xp_m,vg.PU,ps.propparams.dpropagatefcn_dx,ps.propparams.stateCovariance,@propparams_no_noise.propagatefcn,[],propparams_no_noise);
        vg.mu = vg.mu + ps.propparams.stateMean';
    else
        [vg.mu,vg.PP] = ukf_predict1(vg.xp_m,vg.PU,@propparams_no_noise.propagatefcn,ps.propparams.stateCovariance,propparams_no_noise);
        vg.mu = vg.mu + ps.propparams.stateMean';
    end
    vg.PP = (vg.PP+vg.PP')/2;
    vg.PP = cov_regularize(vg.PP);
    
    %% propagate particle and Step through the lambda
    vg.R = ps.likeparams.totalGmmCovariance;
    vg.zeta = ps.likeparams.totalGmmMean';
    
    vg.mu_0 = ps.propparams.propagatefcn(vg.xp_m,propparams_no_noise);
    vg.new_particles = vg.xp_prop;
    vg.edh_linearization = vg.mu_0;
    
    lambda_prev = 0;
    for lambda = ps.setup.lambda_range
        step_size = lambda-lambda_prev;
        
        % Calculate the slopes for moving the particles
        [slope_deterministic,slope]  = calculateSlope_PFPF_EDH...
            (z(:,tt),vg.edh_linearization,vg.new_particles,vg.mu_0,vg.PP,vg.R,vg.zeta,ps,lambda);
        
        vg.edh_linearization= vg.edh_linearization + step_size*slope_deterministic;
        vg.new_particles = vg.new_particles + step_size*slope;  % euler update of particles
        lambda_prev = lambda;
    end
    
    % weight update and componentwise mean calculation and resampling
    
    logProcessDensity = log_GMM_state_transition(vg.new_particles,vg.xp_prop_deterministic,ps.propparams);
    logProposalDensity = log_GMM_state_transition(vg.xp_prop,vg.xp_prop_deterministic,ps.propparams);
    llh = ps.likeparams.llh(vg.new_particles,z(:,tt),ps.likeparams);
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
    
    vg.xp_m = particle_estimate(vg.logW,vg.new_particles);
    output.x_est(:,tt) = vg.xp_m;
    %% EKF/UKF update
    
    if ps.setup.EKF
        [~, vg.PU] = ...
            ekf_update1(vg.mu,vg.PP,z(:,tt)-vg.zeta,ps.likeparams.dh_dx_func,vg.R,ps.likeparams.h_func,[],ps.likeparams);
    else
        [~, vg.PU] = ...
            ukf_update1(vg.mu,vg.PP,z(:,tt)-vg.zeta,ps.likeparams.h_func,vg.R,ps.likeparams);
    end
    vg.PU = (vg.PU+vg.PU')/2;
    vg.PU = cov_regularize(vg.PU);
    
    [vg.logW,vg.xp,~] = resample_PFPF_EDH(vg.logW,vg.new_particles,ps);
end

output.x = ps.x;
output.execution_time = toc;
alg_name = 'PFPF_EDH';
calculateErrors(output,ps,alg_name);
disp(['avg. ESS: ',num2str(mean(output.Neff))])
end