function [vg,output] = initializationFilter(ps)
%% initialize filter outputs
% inputs:
% ps: structure with filter and simulation parameters
%
% Output:
% vg: a struct that contains the filter output
% output: a struct that contains the filter estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = ps.setup.T; % number of time steps

output = struct('x_est',zeros(ps.setup.dimState_all,T));%,'M',zeros(ps.setup.dimState,T), 'PU',zeros(ps.setup.dimState,ps.setup.dimState,T));
% output.x_est_unweighted = output.x_est;

%  Kalman variables
 vg.M = zeros(ps.setup.dimState_all,1); % EKF1/EKF2/UKF estimate
 vg.PP = zeros(ps.setup.dimState_all,ps.setup.dimState_all); % EKF1/EKF2/UKF pred. covariance
 vg.PU = zeros(ps.setup.dimState_all,ps.setup.dimState_all); % EKF1/EKF2/UKF upd. covariance

[vg.xp,vg.M,vg.PU,vg.dist_xp] = ps.initparams.initfcn(ps,ps.setup.nParticle); % Initialize particles, KF estimate and cov
vg.xp_m = vg.M;
vg.logW = zeros(1,size(vg.xp,2));

end