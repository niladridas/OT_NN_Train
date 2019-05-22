function log_process = log_process_density(vg,ps)

switch ps.setup.example_name
    case 'Septier16'
        log_process = loggausspdf(vg.xp,vg.xp_prop_deterministic,ps.propparams.stateCovariance);
%     otherwise
%         log_process = loggausspdf(vg.xp,vg.xp_prop_deterministic,ps.propparams.Q);
end
end