function ps = InteractingNonStationaryGrowth_initialization(ps)
%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the simulation setup and filter parameter values for the
% "Septier16" example
% Input:
% ps: a structure containg the simulation setup and filter parameters.
% Output:
% ps: a structure containg the simulation setup and filter parameters.
%%%%%%%%%%%%%%%%%%%%%%%%

ps.setup.nTrack = 200;% the number of different states/measurements.
ps.setup.nAlg_per_track = 1;% the number of different execution for each algorithms on the same measurements.
ps.setup.nParticle = 50;% number of particles used in particle flow type algorithms.


%% Simulation area
ps.setup.T=50;
init_dist = 'gaussian';%delta: known true initial state / gaussian : sampled from a gaussian
dim = 64;
ps.setup.dimState = dim;
%% Prior
x0 = zeros(ps.setup.dimState,1);

switch init_dist
    case 'delta'
        initCov = eps*eye(dim);
    case 'gaussian'
        initCov = 1*eye(dim);    
end

ps.initparams = struct(...
    'x0',x0,...
    'init_dist',init_dist,...
    'initCov',initCov,...
    'initfcn',@InteractingNonStationaryGrowth_init);
%% Dynamic model

sigma_x = 0.5;
gmmNumComp_dyn = 3;
gmmMixWeight_dyn = ones(1,gmmNumComp_dyn)/gmmNumComp_dyn;
gmmMeans_dyn = 1*[-1*ones(1,ps.setup.dimState);zeros(1,ps.setup.dimState);ones(1,ps.setup.dimState)];
gmmCovariance_dyn = repmat(sigma_x^2*eye(ps.setup.dimState,ps.setup.dimState),1,1,gmmNumComp_dyn);

gmmNoiseModel_dyn = gmdistribution(gmmMeans_dyn,gmmCovariance_dyn,gmmMixWeight_dyn);

% parameters of GMM noise in dynamic model

stateMean = sum(gmmMeans_dyn.*repmat(gmmMixWeight_dyn',1,size(gmmMeans_dyn,2)),1);
stateCovariance = calculateTotalCovarianceGMM(gmmNoiseModel_dyn); %

stateCovariance=(stateCovariance + stateCovariance')/2;

ps.propparams = struct(...
    'stateMean',stateMean,...
    'stateCovariance',stateCovariance,...
    'prop_type','gmm',...
    'gmmNoiseModel_dyn',gmmNoiseModel_dyn,...
    'propagatefcn',@InteractingNonStationaryGrowth_propagate,...
    'dpropagatefcn_dx',@InteractingNonStationaryGrowth_dpropagatefcn_dx,...
    'log_GMM_state_transition',@log_GMM_state_transition);

%% Measurement model
%multimodal likelihood
sigma_z = 0.1;
gmmNumComp = 3;
gmmMixWeight = ones(1,gmmNumComp)/gmmNumComp;
gmmMeans = 3*[-1*ones(1,ps.setup.dimState);zeros(1,ps.setup.dimState);ones(1,ps.setup.dimState)];
gmmCovariance = repmat(sigma_z^2*eye(ps.setup.dimState,ps.setup.dimState),1,1,gmmNumComp);
%% heavy tail likelihood
% gmmNumComp = 2;
% gmmMixWeight = [0.8,0.2];
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
totalGmmMean = sum(gmmMeans.*repmat(gmmMixWeight',1,size(gmmMeans,2)),1);
totalGmmCovariance = calculateTotalCovarianceGMM(gmmNoiseModel); % this particular step is required 
% exclusively for EKF comparison,EKF runs altogether differently in my
% setup. I do not need it for prior estimation 

gmmNoiseModel_singlegauss = gmdistribution(totalGmmMean,totalGmmCovariance);

ps.likeparams = struct(...
    'gmmNoiseModel', gmmNoiseModel,...
    'gmmNoiseModel_singlegauss',gmmNoiseModel_singlegauss,... %for ledh
    'totalGmmMean',totalGmmMean,...             %not required if result not compared to EKF/UKF
    'totalGmmCovariance',totalGmmCovariance,... %not required if result not compared to EKF/UKF
    'h_func',@InteractingNonStationaryGrowth_hfunc,...
    'dh_dx_func',@InteractingNonStationaryGrowth_dh_dxfunc,...
    'observation_noise','gmm');

switch ps.likeparams.observation_noise
    case 'gmm'
        ps.likeparams.llh = @GMM_llh;
end
%%
ps.setup.plotfcn = @NonStationaryGrowth_ParticlePlot;               % plotting function
end