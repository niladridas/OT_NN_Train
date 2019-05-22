function [CompMean,CompCov,proportion] = merge_GMM_Comp(CompMean,CompCov,proportion,numComp)
proportion = proportion(proportion>0);
CompMean = CompMean(:,proportion>0);
CompCov = CompCov(:,:,proportion>0);

gmmModel = gmdistribution(CompMean',CompCov,proportion);
totalGmmCovariance = calculateTotalCovarianceGMM(gmmModel);
totalMean = sum(CompMean.*repmat(proportion,size(CompMean,1),1),2);
invCov = inv(totalGmmCovariance);

while(length(proportion)>numComp)
    dist = zeros(length(proportion));
    for i = 1: length(proportion)
        for j = 1:i
            if i>j
                dist(i,j) = proportion(i)* proportion(j)/(proportion(i)+proportion(j))*...
                    (CompMean(:,i)-CompMean(j))'*invCov*(CompMean(:,i)-CompMean(j));
                dist(j,i) = dist(i,j);
            else
                
                dist(i,i) = Inf;
            end
        end
    end
    [minrow,yid] = min(dist);
    [~, xid] = min(minrow);
    yid = yid(xid);
    newCompProp = proportion(xid) + proportion(yid);
    newCompMean = (proportion(xid)*CompMean(:,xid) + proportion(yid)*CompMean(:,yid))/newCompProp;
    newCompCov =  (proportion(xid)*CompCov(:,:,xid)+proportion(yid)*CompCov(:,:,yid))/newCompProp + ...
        proportion(xid)*proportion(yid)/(newCompProp^2)*(CompMean(:,xid)-CompMean(:,yid))*(CompMean(:,xid)-CompMean(:,yid))';
    newCompCov = cov_regularize(newCompCov);
    
    proportion([xid,yid]) =[];
    CompMean(:,[xid,yid]) = [];
    CompCov(:,:,[xid,yid]) = [];
    
    proportion =[proportion, newCompProp];
    CompMean = [CompMean,newCompMean];
    CompCov(:,:,end+1) = newCompCov;   
end
end