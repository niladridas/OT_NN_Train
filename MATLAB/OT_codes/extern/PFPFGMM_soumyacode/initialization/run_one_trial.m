function tracking_output = run_one_trial(ps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the interface used to run different filtering algorithms.
%
% Inputs:
% ps: structure with filter and simulation parameters
%
% Output:
% tracking_output: a struct that contains outputs from different filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tracking_output = [];

trial_ix = ps.trial_ix;

ps.x = ps.x_all{ceil(trial_ix/ps.setup.nAlg_per_track)};
y = ps.y_all{ceil(trial_ix/ps.setup.nAlg_per_track)};


if ismember('EKF_GMM',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    tracking_output.EKF_GMM = EKF_UKF_GMM(ps_new,y);
end

if ismember('UKF_GMM',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.EKF = false;
    tracking_output.UKF_GMM = EKF_UKF_GMM(ps_new,y);
end

if ismember('PF_GMM_LEDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    tracking_output.PF_GMM_LEDH = DH_ExactFlow_Filter_GMM(ps_new,y);
end

if ismember('PF_GMM_EDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.LEDH = false;
    ps_new.setup.nParticle = 1e4;
    tracking_output.PF_GMM_EDH = DH_ExactFlow_Filter_GMM(ps_new,y);
end

if ismember('PFPF_GMM',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 200;
    tracking_output.PFPF_GMM = PFPF_GMM(ps_new,y);
end

if ismember('GSPF',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 1e4;
    tracking_output.GSPF = GSPF(ps_new,y);
end

if ismember('EKF',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    tracking_output.EKF = EKF_UKF_Filter(ps_new,y);
end

if ismember('UKF',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.EKF = false;
    tracking_output.UKF = EKF_UKF_Filter(ps_new,y);
end

if ismember('LEDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 500;
    tracking_output.LEDH = LEDH_EDH(ps_new,y);
end

if ismember('EDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.LEDH = false;
    ps_new.setup.nParticle = 500;
    tracking_output.EDH = LEDH_EDH(ps_new,y);
end

if ismember('PFPF_LEDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 500;
    tracking_output.PFPF_LEDH = PFPF_LEDH(ps_new,y);
end

if ismember('PFPF_EDH',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 1e5;
    tracking_output.PFPF_EDH = PFPF_EDH(ps_new,y);
end

if ismember('BPF_1e6',ps.setup.algs_executed)
    rng(ps.setup.random_seeds(trial_ix),'twister');
    ps_new = ps;
    ps_new.setup.nParticle = 1e6;
    tracking_output.BPF_1e6 = BootstrapParticleFilter(ps_new,y);
end

% if ismember('APF_1e6',ps.setup.algs_executed)
%     rng(ps.setup.random_seeds(trial_ix),'twister');
%     ps_new = ps;
%     ps_new.setup.nParticle = 1e6;
%     tracking_output.APF_1e6 = AuxiliaryParticleFilter(ps_new,y);
% end
%% Displaying estimation errors for all algorithms
% disp('------------------------------------------------------------')
% disp('------------------------------------------------------------')
% disp(['Displaying results for all algorithms, at Trial ', num2str(trial_ix), ':']);
% alg_field_names = fieldnames(tracking_output);
% for field_ix = 1:length(alg_field_names)
%     field_name = alg_field_names{field_ix};
%     calculateErrors(tracking_output.(field_name),ps,field_name);
% end
% disp('------------------------------------------------------------')
% disp('------------------------------------------------------------')
end