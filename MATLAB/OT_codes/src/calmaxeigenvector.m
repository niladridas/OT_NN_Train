% Author: *Niladri Das*
function [maxeigenvector]= calmaxeigenvector(PriorMAT)
% calmaxeigenvector  Finds the eigenvector corresponding to the largest
% eigen value
%   [maxeigenvector] = calmaxeigenvector(PriorMAT) 
%   PriorMAT is the positive semidefinite matrix, representing the state
%   error covariance matrix
[V,D] = eig(PriorMAT);
lam = diag(D);
[~,I] = max((lam));
maxeigenvector = V(:,I(1));
end