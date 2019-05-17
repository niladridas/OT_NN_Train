function xp_new = Septier16_propagate(xp,prop_params)

nParticle = size(xp,2);
%Linear stste dynamics
xp_new = prop_params.transitionMatrix*xp;
switch prop_params.prop_type
    case 'gaussian'
        xp_new = xp_new + mvnrnd(prop_params.stateMean,prop_params.stateCovariance,nParticle)';
    case 'gmm'
        xp_new = xp_new + random(prop_params.gmmNoiseModel_dyn,nParticle)';
end
end