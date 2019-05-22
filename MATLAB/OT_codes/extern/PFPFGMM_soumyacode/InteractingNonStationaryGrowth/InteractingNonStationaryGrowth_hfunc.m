function z = InteractingNonStationaryGrowth_hfunc(xp , likeparams)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non Linear measurement model 
z =(xp.^2)/20;
end