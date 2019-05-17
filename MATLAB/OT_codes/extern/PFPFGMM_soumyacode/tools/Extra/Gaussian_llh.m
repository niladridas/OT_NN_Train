function llh = Gaussian_llh(xp,z_current,likeparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

llh = loggausspdf(bsxfun(@minus,likeparams.h_func(xp,likeparams),z_current),zeros(size(z_current,1),1),likeparams.R);

end

