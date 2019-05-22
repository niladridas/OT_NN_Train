function ps = NonStationaryGrowth_initialization(ps)
%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the simulation setup and filter parameter values for the
% "Septier16" example
% Input:
% ps: a structure containg the simulation setup and filter parameters.
% Output:
% ps: a structure containg the simulation setup and filter parameters.
%%%%%%%%%%%%%%%%%%%%%%%%

ps.setup.nTrack = 100;% the number of different states/measurements.
ps.setup.nAlg_per_track = 5;% the number of different execution for each algorithms on the same measurements.
ps.setup.nParticle = 50;% number of particles used in particle flow type algorithms.


%% Simulation area
ps.setup.T=50;
init_dist = 'delta';%delta means we know the true initial state and set it that way.
dim = 1;
ps.setup.dimState = dim;

%% Prior
x0 = zeros(ps.setup.dimState,1);

stateCovariance=3;

stateCovarianceSR=sqrt(stateCovariance);
stateCovarianceInv=inv(stateCovariance);

switch init_dist
    case 'delta'
        initCov = eps*eye(size(stateCovariance));
end

ps.initparams = struct(...
    'x0',x0,...
    'init_dist',init_dist,...
    'initCov',initCov,...
    'initfcn',@NonStationaryGrowth_init);

%% Dynamic model
ps.propparams = struct(...
    'stateCovariance',stateCovariance,...
    'stateCovarianceSR',stateCovarianceSR,...
    'stateCovarianceInv',stateCovarianceInv,...
    'propagatefcn',@NonStationaryGrowth_propagate,...
    'dpropagatefcn_dx',@NonStationaryGrowth_dpropagatefcn_dx);
%% Measurement model
%multimodal likelihood
sigma_z = 0.1;
gmmNumComp = 3;
gmmMixWeight = ones(1,gmmNumComp)/gmmNumComp;
gmmMeans = 5*[-1*ones(1,ps.setup.dimState);zeros(1,ps.setup.dimState);ones(1,ps.setup.dimState)];
gmmCovariance = repmat(sigma_z^2*eye(ps.setup.dimState,ps.setup.dimState),1,1,gmmNumComp);
%% heavy tail likelihood
% gmmNumComp = 2;
% gmmMixWeight = [0.9,0.1];
% gmmMeans = zeros(gmmNumComp,ps.setup.dimState);
% gmmCovariance = zeros(ps.setup.dimState,ps.setup.dimState,gmmNumComp);
% 
% gmmCovariance(:,:,1) = 0.01*eye(ps.setup.dimState,ps.setup.dimState);
% gmmCovariance(:,:,2) = 0.1*eye(ps.setup.dimState,ps.setup.dimState);
 %% single gaussian
% gmmNumComp = 1;
% gmmMixWeight = 1;
% gmmMeans = zeros(gmmNumComp,ps.setup.dimState);
% gmmCovariance = zeros(ps.setup.dimState,ps.setup.dimState);
% 
% gmmCovariance(:,:,1) = 0.01*eye(ps.setup.dimState,ps.setup.dimState);
%%
gmmNoiseModel = gmdistribution(gmmMeans,gmmCovariance,gmmMixWeight);
% parameters of GMM noise in measurement model
totalGmmMean = sum(gmmMeans.*repmat(gmmMixWeight',1,dim),1);
totalGmmCovariance = calculateTotalCovarianceGMM(gmmNoiseModel); % this particular step is required 
% exclusively for EKF comparison,EKF runs altogether differently in my
% setup. I do not need it for prior estimation % Required to Map supp(X_k) to supp(Y_k)
gmmNoiseModel_singlegauss = gmdistribution(zeros(1,dim),totalGmmCovariance); %for LEDH
ps.likeparams = struct(...
    'gmmNoiseModel', gmmNoiseModel,...
    'gmmNoiseModel_singlegauss',gmmNoiseModel_singlegauss,... %for ledh
    'totalGmmMean',totalGmmMean,...             %not required if result not compared to EKF/UKF
    'totalGmmCovariance',totalGmmCovariance,... %not required if result not compared to EKF/UKF
    'h_func',@NonStationaryGrowth_hfunc,...
    'dh_dx_func',@NonStationaryGrowth_dh_dxfunc,...
    'dH_dx_func',@NonStationaryGrowth_dH_dx,...
    'observation_noise','gmm');

switch ps.likeparams.observation_noise
    case 'gmm'
        ps.likeparams.llh = @GMM_llh;
end

%%
ps.setup.plotfcn = @NonStationaryGrowth_ParticlePlot;               % plotting function
