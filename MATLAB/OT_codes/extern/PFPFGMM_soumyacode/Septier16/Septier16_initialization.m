function ps = Septier16_initialization(ps)
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
init_dist = 'gaussian';%delta means we know the true initial state and set it that way.
dim = 64;
ps.setup.dimState = dim;
% ps.setup.domain=[1 sqrt(ps.setup.dimState)];%[0 30];
%
% ps.setup.nSensor=ps.setup.dimState; % Dimension of the observation
% % Sensors placed on a regular grid
% [XSensors,YSensors]=meshgrid(linspace(ps.setup.domain(1),ps.setup.domain(2),sqrt(ps.setup.nSensor)));
% ps.setup.sensorPosition=[reshape(XSensors,1,[]);reshape(YSensors,1,[])];

%% Prior
x0 = zeros(ps.setup.dimState,1);

switch init_dist
    case 'delta'
        initCov = eps*eye(dim);
    case 'gaussian'
        initCov = 1*eye(dim);
end

alpha = 0.9;

ps.initparams = struct(...
    'x0',x0/alpha,...
    'init_dist',init_dist,...
    'initCov',initCov,...
    'initfcn',@Septier16_init);

%% Dynamic model

sigma_x = 1;
gmmNumComp_dyn = 3;
gmmMixWeight_dyn = ones(1,gmmNumComp_dyn)/gmmNumComp_dyn;
gmmMeans_dyn = 1*[-1*ones(1,ps.setup.dimState);zeros(1,ps.setup.dimState);ones(1,ps.setup.dimState)];
gmmCovariance_dyn = repmat(sigma_x^2*eye(ps.setup.dimState,ps.setup.dimState),1,1,gmmNumComp_dyn);

gmmNoiseModel_dyn = gmdistribution(gmmMeans_dyn,gmmCovariance_dyn,gmmMixWeight_dyn);

% parameters of GMM noise in dynamic model

stateMean = sum(gmmMeans_dyn.*repmat(gmmMixWeight_dyn',1,size(gmmMeans_dyn,2)),1);
stateCovariance = calculateTotalCovarianceGMM(gmmNoiseModel_dyn); %

stateCovariance=(stateCovariance + stateCovariance')/2;
transitionMatrix=alpha*eye(ps.setup.dimState);

ps.propparams = struct(...
    'transitionMatrix',transitionMatrix,...
    'stateMean',stateMean,...
    'stateCovariance',stateCovariance,...
    'prop_type','gmm',...
    'gmmNoiseModel_dyn',gmmNoiseModel_dyn,...
    'propagatefcn',@Septier16_propagate,...
    'dpropagatefcn_dx',@Septier16_dpropagatefcn_dx,...
    'log_GMM_state_transition',@log_GMM_state_transition);

%% Measurement model
%multimodal likehood
sigma_z = 0.1;
gmmNumComp = 3;
gmmMixWeight = ones(1,gmmNumComp)/gmmNumComp;
gmmMeans = 5*[-1*ones(1,ps.setup.dimState);zeros(1,ps.setup.dimState);ones(1,ps.setup.dimState)];
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
% gmmMeans = zeros(gmmNumComp,dim);
% gmmCovariance = zeros(dim,dim);
%
% gmmCovariance(:,:,1) = 0.01*eye(dim,dim);
%%
gmmNoiseModel = gmdistribution(gmmMeans,gmmCovariance,gmmMixWeight);
% parameters of GMM noise in measurement model
totalGmmMean = sum(gmmMeans.*repmat(gmmMixWeight',1,dim),1);
totalGmmCovariance = calculateTotalCovarianceGMM(gmmNoiseModel); % this particular step is required
% exclusively for EKF comparison,EKF runs altogether differently in my
% setup. I do not need it for prior estimation

%want to test how ledh works if single gaussian is assumed. so we overwrite
%the noise model by a gmm with 1 component , same mean and covariance as
%gmm. This is for comparison, not required in the original code.


gmmNoiseModel_singlegauss = gmdistribution(totalGmmMean,totalGmmCovariance); %for LEDH


observationTransition=eye(ps.setup.dimState); % Required to Map supp(X_k) to supp(Y_k)

ps.likeparams = struct('observationTransition', observationTransition,...
    'gmmNoiseModel', gmmNoiseModel,...
    'gmmNoiseModel_singlegauss',gmmNoiseModel_singlegauss,... %for ledh
    'totalGmmMean',totalGmmMean,...             %not required if result not compared to EKF/UKF
    'totalGmmCovariance',totalGmmCovariance,... %not required if result not compared to EKF/UKF
    'h_func',@Septier16_hfunc,...
    'dh_dx_func',@Septier16_dh_dxfunc,...
    'observation_noise','gmm');%'poisson','normal'

switch ps.likeparams.observation_noise
    case 'normal'
        ps.likeparams.llh = @Gaussian_llh;
    case 'poisson'
        ps.likeparams.llh = @Poisson_llh;
    case 'gmm'
        ps.likeparams.llh = @GMM_llh;
end

%%
ps.setup.plotfcn = @Septier16_ParticlePlot;               % plotting function
