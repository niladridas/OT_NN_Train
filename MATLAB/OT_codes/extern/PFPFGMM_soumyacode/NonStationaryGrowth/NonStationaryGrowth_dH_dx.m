function dH_dx = NonStationaryGrowth_dH_dx(xp, likeparams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hessian of measurement function for the poisson mean used in the
% "Septier16" example.
%
% measurement function for the poisson mean
%
% Inputs:
% xp: state value
%  
% likeparams: a structure containing the field 'ObservationTransition',
%             'MappingIntensityCoeff','MappingIntensityScale'.
%
% Output:
%
%  dhdx = MappingIntensityCoeff*exp(xp/MappingIntensityScale)/MappingIntensityScale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,nParticle] = size(xp);

if nParticle > 1
    error('This function is only defined for a single state')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non Linear Measurement model
dH_dx = 1/10;
end