function  [logW,xp,PU,eff] = resample_PFPF_GMM(logW,xp,PU,ps)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform resampling and add regularized noise if required.
%
% Inputs:
% vg: a struct that contains the filter output
% ps: structure with filter and simulation parameters
%
% Output:
% vg: a struct that contains the filter output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weights = exp(logW-max(logW));

if sum(weights)==0||any(isnan(weights))
    weights = ones(size(weights))/ps.setup.nParticle;
    eff = 1;
else
    weights = weights/sum(weights);
    eff = 1/sum(weights.^2);
end

if ps.setup.Resampling && eff < ps.setup.Neff_thresh_ratio * ps.setup.nParticle
    I = resample(ps.setup.nParticle,weights,'stratified');
    xp = xp(:,I);
    logW = zeros(size(I));
    PU = PU(:,:,I);
end
end