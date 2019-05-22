function [dist]= omat(X,Y)
nTarget = size(X,2);
d = zeros(factorial(nTarget),1);
allperm = perms(1:nTarget);
for i = 1:length(d)
    d(i) =  sum(sqrt(sum((X-Y(:,allperm(i,:))).^2,1)));
end
dist = min(d)/nTarget;
end