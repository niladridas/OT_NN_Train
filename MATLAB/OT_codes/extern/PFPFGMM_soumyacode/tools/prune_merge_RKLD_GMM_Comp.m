function [CompMean,CompCov,proportion,idx] = prune_merge_RKLD_GMM_Comp(CompMean,CompCov,proportion,numComp)
proportion = proportion(proportion>0);
CompMean = CompMean(:,proportion>0);
CompCov = CompCov(:,:,proportion>0);

proportion1 = proportion;
CompMean1 = CompMean;
CompCov1 = CompCov;
while(length(proportion)>numComp)
    cost = zeros(length(proportion),length(proportion));
    for i = 1:length(proportion)
        for j = 1:i
            if i>j
                 cost(i,j) = merging_cost(proportion(i),CompMean(:,i),CompCov(:,:,i),...
                     proportion(j),CompMean(:,j),CompCov(:,:,j));
               % cost(i,j) = Inf;
                cost(j,i) = cost(i,j);
            else
                cost(i,i) = pruning_cost(proportion1,CompMean1,CompCov1,i);
            end
        end
    end
    [minrow,yid] = min(cost);
    [~, xid] = min(minrow);
    yid = yid(xid);
    if xid~=yid
        newCompProp = proportion(xid) + proportion(yid);
        
        newCompMean = (proportion(xid)*CompMean(:,xid) + proportion(yid)*CompMean(:,yid))/newCompProp;
        
        newCompCov = (proportion(xid)*CompCov(:,:,xid) + proportion(yid)*CompCov(:,:,yid))/newCompProp + ...
            proportion(xid)*proportion(yid)/(newCompProp^2)*(CompMean(:,yid)-CompMean(:,xid))*(CompMean(:,yid)-CompMean(:,xid))';
        newCompCov = cov_regularize(newCompCov);
        
        proportion([xid,yid]) =[];
        CompMean(:,[xid,yid]) = [];
        CompCov(:,:,[xid,yid]) = [];
        
        proportion =[proportion, newCompProp];
        CompMean = [CompMean,newCompMean];
        CompCov(:,:,end+1) = newCompCov;
    else
        proportion(xid) =[];
        proportion = proportion/sum(proportion);
        CompMean(:,xid) = [];
        CompCov(:,:,xid) = [];
    end
end
% idx = [];
% for i = 1:size(CompMean,2)
%     for j = 1:size(CompMean1,2)
%         if CompMean(:,i) == CompMean1(:,j)
%             idx = [idx,j];
%         end
%     end
% end
end