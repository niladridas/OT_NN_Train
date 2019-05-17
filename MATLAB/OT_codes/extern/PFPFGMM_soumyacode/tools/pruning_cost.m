function [cost] = pruning_cost(alpha,mu,sigma,i)
criterion = zeros(1,length(alpha));
for j = 1:length(alpha)
    if j==i
        criterion(i) = Inf;
    else
         criterion(j) = -log(1-alpha(i)) - ...
        alpha(j)/(1-alpha(i))*log(1+alpha(i)/alpha(j)*exp(-1*KL_divergence(mu(:,j),sigma(:,:,j),mu(:,i),sigma(:,:,i)))); 
    end
end
cost = min(criterion);
end