function cova = cov_regularize(cova)
% Regularize covariance matrices if the Cholesky factorization
% is not positive definite
% Avoids possible problems with the (pseudo-)inverse in later code
%
% Loop, adding identity matrix (multiplied by very small constant), until
% Cholesky factorization is positive definite
% Exception if too any iterations of the loop
%
dim = size(cova,1);
reg = eye(dim,dim)*1e-14;

if any(isinf(cova))
    cova =  reg;
elseif any(isnan(cova))
    cova =  reg;
elseif ~all(isreal(cova))
    cova =  reg;    
end
% Perform Cholesky decomposition
[~,indicator]=chol(cova);
count = 0; maxCount = 1e2;

% Check whether the factorization is positive definite
% If not add regularization matrix to covariance matrix
while indicator>0 && count < maxCount && any(eig(cova))<=eps
    cova=cova+reg;
    [~,indicator]=chol(cova);
    count=count+1;
end

% throw an exception if positive-definiteness cannot be achieved
if count == maxCount
    %     ME = MException('cov_regularize:TooManyIterations', ...
    %         'Could not regularize the covariance matrix');
    %     throw(ME);
    %   warning(['cov_regularize:TooManyIterations', ...
    %         'Could not regularize the covariance matrix']);
    cova = eye(dim,dim)*max(eig(cova)) + reg;
end
end