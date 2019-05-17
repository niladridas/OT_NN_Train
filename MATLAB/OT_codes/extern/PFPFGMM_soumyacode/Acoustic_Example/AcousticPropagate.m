function [xp] = AcousticPropagate(xp,prop_params)

[dim,nParticle] = size(xp);

Phi = prop_params.Phi;
Sigma_sqrt = prop_params.stateCovarianceSR;

xp = Phi*xp + Sigma_sqrt*randn(dim,nParticle);

end