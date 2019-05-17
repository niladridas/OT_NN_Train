function [output] = LEDH_EDH(ps,z)
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
    %% redraw particles
    if tt~=1
        vg.xp = mvnrnd(vg.xp_m',vg.PU,nParticle)';
    end
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
    vg.xp = ps.propparams.propagatefcn(vg.xp,ps.propparams);
    
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
    
    mu_0 = ps.propparams.propagatefcn(vg.xp_m,propparams_no_noise);
    P = vg.PP;
    R = ps.likeparams.totalGmmCovariance;
    zeta = ps.likeparams.totalGmmMean';
    
    lambda_prev = 0;
    for lambda = ps.setup.lambda_range
        step_size = lambda-lambda_prev;
        
        % Calculate the slopes for moving the particles
        if ps.setup.LEDH
            slope = calculateSlopeComponentwise_LEDH(z(:,tt),vg.xp,mu_0,P,R,zeta,ps,lambda);
        else
            slope = calculateSlopeComponentwise_EDH(z(:,tt),vg.xp,mu_0,P,R,zeta,ps,lambda);
        end
        
        vg.xp = vg.xp + step_size*slope;  % euler update of particles
        lambda_prev = lambda;
    end
    vg.xp_m = mean(vg.xp,2);
    output.x_est(:,tt) = vg.xp_m;
    %% EKF/UKF update
    
    if ps.setup.EKF
        [~, vg.PU] = ...
            ekf_update1(vg.mu,vg.PP,z(:,tt)-zeta,ps.likeparams.dh_dx_func,R,ps.likeparams.h_func,[],ps.likeparams);
    else
        [~, vg.PU] = ...
            ukf_update1(vg.mu,vg.PP,z(:,tt)-zeta,ps.likeparams.h_func,R,ps.likeparams);
    end
    vg.PU = (vg.PU+vg.PU')/2;
    vg.PU = cov_regularize(vg.PU);
end

output.x = ps.x;
output.execution_time = toc;
if ps.setup.LEDH
    alg_name = 'LEDH';
else
    alg_name = 'EDH';
end
calculateErrors(output,ps,alg_name);
end