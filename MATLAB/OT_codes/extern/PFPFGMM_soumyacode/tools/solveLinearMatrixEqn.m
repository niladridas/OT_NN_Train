function [X] = solveLinearMatrixEqn(A,B,F)
p = size(A,3);

mu = 0;
for i = 1:p
    mu = mu + max(eig(A(:,:,i)*A(:,:,i)'))*max(eig(B(:,:,i)'*B(:,:,i)));
end

mu = 0.5/mu;

X = 1e-6*ones(size(A,1));

Xj = zeros([size(X),p]);

residual = calcResidual(A,B,F,X,p);

%while(1)
for ite = 1:100000
for j = 1:p
    Xj(:,:,j) = X +mu*A(:,:,j)'*(residual + A(:,:,j)*X*B(:,:,j))*B(:,:,j)';
end
X_prev = X;
X = mean(Xj,3);
residual = calcResidual(A,B,F,X,p);
if norm(X-X_prev)/norm(X)<1e-2 && norm(residual)/norm(F)<1e-2
    break;
end
disp(ite);
end

end

function residual = calcResidual(A,B,F,X,p)
residual = F;
for i = 1:p
    residual = residual -A(:,:,i)*X*B(:,:,i);
end
end

