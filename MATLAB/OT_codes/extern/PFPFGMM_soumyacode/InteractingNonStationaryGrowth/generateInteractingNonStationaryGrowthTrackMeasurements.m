function [x,y] = generateInteractingNonStationaryGrowthTrackMeasurements(ps)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates target trajectories for the Septier 16 example
%
% Input:
% ps: parameter structure, including initial states and propagation
% parameters
%
% Output: the state trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%

dimState = ps.setup.dimState;

x = zeros(dimState,ps.setup.T);
y = zeros(dimState,ps.setup.T);

for tt = 1:ps.setup.T
    ps.propparams.time_step = tt;
    if tt == 1
        x(:,tt)=ps.propparams.propagatefcn(ps.initparams.x0,ps.propparams);
    else
        x(:,tt)=ps.propparams.propagatefcn(x(:,tt-1),ps.propparams);
    end
    
    measurement_mean=ps.likeparams.h_func(x(:,tt),ps.likeparams);
    
    switch ps.likeparams.observation_noise
        case 'gmm'
            y(:,tt)=measurement_mean+ random(ps.likeparams.gmmNoiseModel,1)';
    end
end
end