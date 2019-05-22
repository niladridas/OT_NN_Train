function [CompMean,CompCov,proportion,idx] = prune_RKLD_GMM_Comp(CompMean,CompCov,proportion,numComp)
cost = zeros(1,length(proportion));

for i = 1:length(proportion)   
    cost(i) = pruning_cost(proportion,CompMean,CompCov,i);
end
[~,indices] = sort(cost,'descend');
idx = indices(1:numComp);
CompMean = CompMean(:,idx);
CompCov = CompCov(:,:,idx);
proportion = proportion(idx); 

%%
% CompMean1 = CompMean;
% 
% while length(proportion)>numComp
%     cost = zeros(1,length(proportion));
%     
%     for i = 1:length(proportion)
%         cost(i) = pruning_cost(proportion,CompMean,CompCov,i);
%     end
%     [~,index] = min(cost);
%     CompMean(:,index) = [];
%     CompCov(:,:,index) = [];
%     proportion(index) = [];
%     proportion = proportion/sum(proportion);
% end
% [~,indices] = sort( proportion,'descend');
% CompMean = CompMean(:,indices);
% CompCov = CompCov(:,:,indices);
% proportion = proportion(indices);
% 
% idx = [];
% for j = 1:size(CompMean,2)
%     for i = 1:size(CompMean1,2)
%         if CompMean(:,j)==CompMean1(:,i)
%             idx = [idx,i];
%         end
%     end
% end
%%
end