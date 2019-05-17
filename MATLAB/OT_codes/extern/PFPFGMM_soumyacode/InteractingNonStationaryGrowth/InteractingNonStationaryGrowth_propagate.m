function xp_new = InteractingNonStationaryGrowth_propagate(xp,prop_params)
tt = prop_params.time_step;
[dimState,nParticle] = size(xp);
xp_new = zeros(size(xp));
a = 0.5;
b = 2.5;
c = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non Linear stste dynamics

xp_new(1,:) = a*xp(1,:)+ b*xp(2,:)./(1+xp(1,:).^2)+c*cos(1.2*(tt-1));

xp_new(2:dimState-1,:) = a*xp(2:dimState-1,:)+ b*xp(3:dimState,:)./(1+xp(1:dimState-2,:).^2)+c*cos(1.2*(tt-1));

xp_new(dimState,:) = a*xp(dimState,:)+ b*xp(dimState,:)./(1+xp(dimState-1,:).^2)+c*cos(1.2*(tt-1));

switch prop_params.prop_type
    case 'gaussian'
        xp_new = xp_new + mvnrnd(prop_params.stateMean,prop_params.stateCovariance,nParticle)';
    case 'gmm'
        xp_new = xp_new + random(prop_params.gmmNoiseModel_dyn,nParticle)';
end
end