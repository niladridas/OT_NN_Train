function ps = Acoustic_example_initialization(ps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes the parameter structure for the acoustic sensor example
%
% Input:
% ps: a structure containg the simulation setup and filter parameters.
%
% Output:
% ps: a structure containg the simulation setup and filter parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ps.setup.nTrack = 100;% the number of different states/measurements. It is set to 100 in [li2016].
ps.setup.nAlg_per_track = 5;% the number of different execution for each algorithms on the same measurements.
ps.setup.nParticle = 125;% number of particles used in particle flow type algorithms.

% Simulation and Measurement Parameters
ps.setup.T=40;
ps.setup.nTarget = 4;
ps.setup.dimState_per_target = 4;

% Error metrics
ps.setup.ospa_c = 40;
ps.setup.ospa_p = 1;

nTarget = ps.setup.nTarget;
load 'sensorsXY';
simAreaSize=40; %size of the area

sensorsPos = simAreaSize/40*sensorsXY; %physical positions of the sensors
nSensor = size(sensorsPos,1);

% Copy the sensor positions so that they are easier for per-target
% processing
sensorsPos = sensorsPos';
sensorsPos = kron(sensorsPos,ones(nTarget,1));

% Dimension of sensorsPos: (2*nTarget) x nSensor
% x-coordinate replicated nTarget timessigma_z = 0.1;

% y-coordinate replicated nTarget times

% See paper for measurement model description
Amp = 10; %amplitude of the sound source (target)
invPow = 1; % rate of decay
d0 = 0.1;

switch nTarget
    case 4
        % Initialization: 4 targets
        x0 = [12 6 0.001 0.001 32 32 -0.001 -0.005 20 13 -0.1 0.01 15 35 0.002 0.002]';
        survRegion = [0,0,40,40];
        trackBounds = [-10,-10,50,50];
end

if length(x0)~= 4* nTarget
    error('incorrect initial state dimension');
end

% Motion model parameters
% For each target x_{k,i} = Phi*x_{k-1,i} + Gamma*sqrt(gammavar)*normrnd(2,1);
Phi = [1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1];


Qii_real = 0.05*[1/3 0 0.5 0; 0 1/3 0 0.5; 0.5 0 1 0; 0 0.5 0 1];

Qii = [3 0 0.1 0; 0 3 0 0.1; 0.1 0 0.03 0; 0 0.1 0 0.03];

Q_ii_correction = 0.1*[1 0 0 0; 0 1 0 0;0 0 0.01 0;0 0 0 0.01];

Phi = kron(eye(nTarget),Phi);
Q = kron(eye(nTarget),Qii);
Q_real = kron(eye(nTarget),Qii_real);
Q_correction = kron(eye(nTarget),Q_ii_correction);

ps.initparams = struct(...
    'x0',x0,...
    'sigma0',[],...
    'nTarget',nTarget,...
    'survRegion',survRegion,...
    'simAreaSize',simAreaSize);

ps.initparams.sigma0 = repmat(10*[1;1;0.1;0.1],nTarget,1);
stateCovariance = Q;
stateCovarianceSR=real(sqrtm(stateCovariance));
stateCovarianceInv=real(inv(stateCovariance));

stateCovariance_real = Q_real;
stateCovarianceSR_real=real(sqrtm(stateCovariance_real));
stateCovarianceInv_real=real(inv(stateCovariance_real));
ps.propparams_real = struct(...
    'Phi',Phi,...
    'Q',Q_real,...
    'stateCovariance',stateCovariance_real,...
    'stateCovarianceSR',stateCovarianceSR_real,...
    'stateCovarianceInv',stateCovarianceInv_real,...
    'nTarget',nTarget);

ps.propparams = struct(...
    'Phi',Phi,...
    'Q',Q_real,...
    'stateCovariance',stateCovariance,...
    'stateCovarianceSR',stateCovarianceSR,...
    'stateCovarianceInv',stateCovarianceInv,...
    'Q_correction', Q_correction,...
    'Q_regularized', Q_correction,...
    'nTarget',nTarget,...
    'propagatefcn',@AcousticPropagate,...
    'dpropagatefcn_dx',@Acoustic_dpropagatefcn_dx,...
    'dimState_per_target',ps.setup.dimState_per_target);
%% Measurement model
%% heavy tail likelihood
% gmmNumComp = 2;
% gmmMixWeight = [0.8,0.2];
% gmmMeans = zeros(gmmNumComp,nSensor);
% gmmCovariance = zeros(nSensor,nSensor,gmmNumComp);
% 
% gmmCovariance(:,:,1) = 0.0025*eye(nSensor,nSensor);
% gmmCovariance(:,:,2) = 0.04*eye(nSensor,nSensor);
 %% single gaussian
gmmNumComp = 1;
gmmMixWeight = 1;
gmmMeans = zeros(gmmNumComp,nSensor);
gmmCovariance = zeros(nSensor,nSensor);

gmmCovariance(:,:,1) = 0.01*eye(nSensor,nSensor);
%%
gmmNoiseModel = gmdistribution(gmmMeans,gmmCovariance,gmmMixWeight);
% parameters of GMM noise in measurement model
totalGmmMean = sum(gmmMeans.*repmat(gmmMixWeight',1,nSensor),1);
totalGmmCovariance = calculateTotalCovarianceGMM(gmmNoiseModel); % this particular step is required
% exclusively for EKF comparison,EKF runs altogether differently in my
% setup. I do not need it for prior estimation

%want to test how ledh works if single gaussian is assumed. so we overwrite
%the noise model by a gmm with 1 component , same mean and covariance as
%gmm. This is for comparison, not required in the original code.


gmmNoiseModel_singlegauss = gmdistribution(totalGmmMean,totalGmmCovariance); %for LEDH


ps.likeparams = struct('sensorsPos',sensorsPos,...
    'Amp',Amp,...
    'd0',d0,...
    'invPow',invPow,...
    'gmmNoiseModel', gmmNoiseModel,...
    'gmmNoiseModel_singlegauss',gmmNoiseModel_singlegauss,... %for ledh
    'totalGmmMean',totalGmmMean,...             %not required if result not compared to EKF/UKF
    'totalGmmCovariance',totalGmmCovariance,... %not required if result not compared to EKF/UKF
    'observation_noise','gmm',...
    'simAreaSize',simAreaSize,...
    'nTarget',nTarget,...
    'nSensor',nSensor,...
    'dimMeasurement_per_target',nSensor,...
    'dimMeasurement_all',nSensor,...);
    'survRegion',survRegion,...
    'trackBounds',trackBounds);


ps.setup.plotfcn = @AcousticParticlePlot;               % plotting function

ps.likeparams.llh = @GMM_llh;
ps.likeparams.h_func = @Acoustic_hfunc;
ps.likeparams.dh_dx_func = @Acoustic_dh_dxfunc;
%ps.likeparams.dH_dx_func = @Acoustic_dH_dx;
ps.initparams.initfcn = @AcousticGaussInit;             % filter initialization function
