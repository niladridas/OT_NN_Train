function run(example_name,algs_executed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab codes associated to the paper
% 	[li2016] Y. Li and M. Coates, "Particle filtering with invertible particle flow",
%             arXiv: 1607.08799, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main function to run the filtering algorithms.
% We compare the proposed PF-PF (LEDH) and PF-PF (EDH) with
% other filtering algorithms in an acoustic tracking example
% or a large sensor field example, specified in initializePS.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function starts by initializating the simulation and
% algorithm parameters, then runs each filter with
% the required number of random trials.
% The outputs are saved in a mat file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% alg_executed: a cell that can include the following algorithm names.
%       PFPF_LEDH: the particle flow particle filter based on LEDH (PF-PF (LEDH))
%       PFPF_EDH: the particle flow particle filter based on EDH (PF-PF (EDH))
%       LEDH: the localized Daum and Huang filter (LEDH)
%       EDH: the exact Daum and Huang filter (EDH)
%       GPFIS: the Gaussian particle flow particle filter (GPFIS)
%       SmHMC: the Sequential Markov chain Monte Carlo based on the
%              Manifold Hamiltonian Monte Carlo kernel (SmHMC).
%              Note that it can be only applied to the Septier16 example,
%              as it requires the target distribution to be log-concave.
%       EKF: the extended Kalman filter (EKF)
%       BPF: the bootstrap particle filter (BPF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Yunpeng Li and Mark Coates
% Date: July 24, 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if the cell algs_executed is not specified, use the following default
% ones.
if nargin < 2
    %     algs_executed = {'EKF_GMM','UKF_GMM','PF_GMM_LEDH','PF_GMM_EDH','PFPF_GMM','GSPF',...
    %  'EKF','UKF','LEDH','EDH','PFPF_LEDH','PFPF_EDH','BPF_1e6'};
    algs_executed = {'EKF_GMM','PF_GMM_LEDH','PFPF_GMM','GSPF',...
      'UKF','LEDH','EDH','PFPF_LEDH','PFPF_EDH','BPF_1e6'};
end

clearvars -except algs_executed example_name
rng('default');

%% Initialize the simulation setup and filter parameter values
addpath('initialization/');
ps_initial = initializePS(algs_executed,example_name);

%% Generate states and measurements for all simulation trials.
ps_initial = generateTracksMeasurements(ps_initial);

%% Produce tracking results for each algorithm and each trial.
output = cell(ps_initial.setup.nTrial,1);

if ps_initial.setup.parallel_run
    ps_initial.doplot = false;
    gcp();
    parfor trial_ix=1:ps_initial.setup.nTrial
        ps_initial_new = ps_initial;
        ps_initial_new.trial_ix = trial_ix;
        output{trial_ix} = run_one_trial(ps_initial_new);
    end
else
    for trial_ix=1:ps_initial.setup.nTrial
        ps_initial.trial_ix = trial_ix;
        output{trial_ix} = run_one_trial(ps_initial);
    end
end

%% Save results in a mat file
saveResults(ps_initial, output);