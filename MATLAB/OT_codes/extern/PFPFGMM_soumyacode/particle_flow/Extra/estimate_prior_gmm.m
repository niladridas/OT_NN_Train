function [prior_dist_obj] = estimate_prior_gmm(xp_prop,logW,likeparams)
%%
% [priorNumComp,clusterid] = calculateNumComp(xp_prop');
%
% theta = zeros(1,priorNumComp);
% mu = zeros(priorNumComp,size(xp_prop,1));
% P = zeros(size(xp_prop,1),size(xp_prop,1),priorNumComp);
% for comp = 1:priorNumComp
%     theta(comp) = sum(clusterid==comp)/sum(clusterid~=0);
%     mu(comp,:) = mean(xp_prop(:,clusterid==comp),2)';
%     P(:,:,comp) = cov_regularize(cov(xp_prop(:,clusterid==comp)'));
% end
% only implementing 'retaining strong components', later will do merger of components
%numRetainedComp = min(likeparams.gmmNumComp,priorNumComp);
% [~,retainedCompID] = sort(theta,'descend');
% retainedCompID = retainedCompID(1:numRetainedComp);
%
% theta = theta(retainedCompID);
% theta = theta/sum(theta);
% mu = mu(retainedCompID,:);
% P = P(:,:,retainedCompID);
%
% Start = struct('mu',mu,'Sigma',P,'ComponentProportion',theta);

%prior_dist_obj = gmdistribution(mu,P,theta);
%%
numRetainedComp = likeparams.gmmNumComp;
if var(logW)>0
    weights = exp(logW);
    weights = weights/sum(weights);
    I = resample(length(weights),weights,'stratified');
    xp_prop = xp_prop(:,I);
end
    options = statset('MaxIter',5000);
    prior_dist_obj =  fitgmdist(xp_prop',numRetainedComp,'Start','plus','Replicates',5,'RegularizationValue',1e-9,'Options',options);
end
