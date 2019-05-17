function dhdx = NonStationaryGrowth_dh_dxfunc(xp , likeparams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of measurement function for the poisson mean used in the
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

 [dim,nParticle] = size(xp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Non Linear measurement model
dhdx = zeros(dim,dim,nParticle);
dhdx(1,1,:) = xp/10;
end