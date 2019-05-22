function totalGmmCovariance = calculateTotalCovarianceGMM(gmmNoiseModel)
alpha = gmmNoiseModel.ComponentProportion;
mu = gmmNoiseModel.mu;
sigma = gmmNoiseModel.Sigma;
totalGmmCovariance  = zeros(size(mu,2));
overallMean = sum(mu.*repmat(alpha',1,size(mu,2)),1)';

for i =1:length(alpha)
  totalGmmCovariance = totalGmmCovariance + alpha(i)*(sigma(:,:,i) + mu(i,:)'*mu(i,:));
end  
totalGmmCovariance = totalGmmCovariance - overallMean*overallMean';
end

