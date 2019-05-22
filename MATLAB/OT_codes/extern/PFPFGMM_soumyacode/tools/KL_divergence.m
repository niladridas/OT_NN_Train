function [kl_div] = KL_divergence(mu1,sigma1,mu2,sigma2)
kl_div = 0.5*(log(det(sigma2))-log(det(sigma1)))+...
    0.5*trace(sigma2\(sigma1-sigma2+(mu2-mu1)*(mu2-mu1)'));
end