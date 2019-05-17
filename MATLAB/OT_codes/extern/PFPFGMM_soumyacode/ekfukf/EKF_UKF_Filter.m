function [output] = EKF_UKF_Filter(ps,z)
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
    
    % Propagate one step
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
            error('The example name does not matche the record');
    end
    
    if ps.setup.EKF
        [vg.M_prior,vg.PP] = ekf_predict1(vg.M,vg.PU,ps.propparams.dpropagatefcn_dx,ps.propparams.stateCovariance,@propparams_no_noise.propagatefcn,[],propparams_no_noise);
        vg.M_prior = vg.M_prior + ps.propparams.stateMean';
        [vg.M,vg.PU] = ekf_update1(vg.M_prior,vg.PP,z(:,tt)-ps.likeparams.totalGmmMean',ps.likeparams.dh_dx_func,ps.likeparams.totalGmmCovariance,ps.likeparams.h_func,[],ps.likeparams);
    else
        [vg.M_prior,vg.PP] = ukf_predict1(vg.M,vg.PU,@propparams_no_noise.propagatefcn,ps.propparams.stateCovariance,propparams_no_noise);
        vg.M_prior = vg.M_prior + ps.propparams.stateMean';
        [vg.M,vg.PU] = ukf_update1(vg.M_prior,vg.PP,z(:,tt)-ps.likeparams.totalGmmMean',ps.likeparams.h_func,ps.likeparams.totalGmmCovariance,ps.likeparams);
    end
    
    % Regularize the Kalman covariance matrix if necessary
    [~,regind] = chol(vg.PU);
    if regind
        vg.PU = cov_regularize(vg.PU);
    end;
    
    output.x_est(:,tt) = vg.M;
end

output.x = ps.x;
output.execution_time = toc;
if ps.setup.EKF
    alg_name = 'EKF';
else
    alg_name = 'UKF';
end
calculateErrors(output,ps,alg_name);
end