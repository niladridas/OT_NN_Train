function log_process_density = log_GMM_state_transition(xp_new,xp_prop_deterministic,propparams)

process_density = pdf(propparams.gmmNoiseModel_dyn,(xp_new-xp_prop_deterministic)');

log_process_density = log(process_density');
end