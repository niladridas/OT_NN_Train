function xp_new = NonStationaryGrowth_propagate(xp,prop_params)
Sigma_sqrt = prop_params.stateCovarianceSR;
tt = prop_params.time_step;
dimState = size(xp,1);
nParticle = size(xp,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non Linear stste dynamics
xp_new = 0.5*xp+ 25*xp./(1+xp.^2)+8*cos(1.2*(tt-1))+Sigma_sqrt*randn(dimState,nParticle);
end