function [CompMean,CompCov,proportion,idx] = resample_GMM_Comp(CompMean,CompCov,proportion,numComp,proportion_th)

if numComp<length(proportion)
  if proportion(numComp)>proportion_th
        idx = 1:numComp;
        CompMean = CompMean(:,idx);
        CompCov = CompCov(:,:,idx);
        proportion = proportion(idx);    
  else
        idx = resample(numComp,proportion(1:numComp)/sum(proportion(1:numComp)),'stratified');
        CompMean = CompMean(:,idx);
        CompCov = CompCov(:,:,idx);
        proportion = ones(1,numComp);
   end
    
    proportion = proportion/sum(proportion);
else
    idx = 1:length(proportion);
end
end