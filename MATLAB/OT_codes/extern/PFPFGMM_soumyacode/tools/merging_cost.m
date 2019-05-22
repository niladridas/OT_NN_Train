function [cost] = merging_cost(alpha1,mu1,sigma1,...
    alpha2,mu2,sigma2)

alpha12 = alpha1 + alpha2;
mu12 = (alpha1*mu1 + alpha2*mu2)/alpha12;
sigma12 = (alpha1*sigma1 + alpha2*sigma2)/alpha12...
    + (alpha1*alpha2/(alpha12^2))*(mu1-mu2)*(mu1-mu2)';
cost = alpha12*log(alpha12) - alpha12*log(alpha1*exp(-1*KL_divergence(mu12,sigma12,mu1,sigma1))...
    +alpha2*exp(-1*KL_divergence(mu12,sigma12,mu2,sigma2)));
end

