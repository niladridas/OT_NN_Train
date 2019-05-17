function llh = GMM_llh(xp,z_current,likeparams)

X = bsxfun(@minus,z_current,likeparams.h_func(xp,likeparams));
y = pdf(likeparams.gmmNoiseModel,X');

llh = log(y');
end