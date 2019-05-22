function [A_GMM,b_GMM] = calculateAbGMM(zc,xp_linearization,H,ps,vg,lambda)
%particle flow equations for GMM prior and GMM likelihood

%% prior parameters calculation
x_dim = length(xp_linearization);

priorGMM = vg.dist_prior;

M =priorGMM.NumComponents;

mu = cell(M,1);
invP = cell(M,1);
%theta = priorGMM.ComponentProportion;
for i = 1:M
    mu{i} = priorGMM.mu(i,:)';
    invP{i} = inv(cov_regularize(priorGMM.Sigma(:,:,i)));
end
alpha = -1*posterior(priorGMM,xp_linearization');
a = zeros(x_dim,M);
common_term_a = zeros(x_dim,1);
for i = 1:M
     common_term_a = common_term_a + alpha(i)*invP{i}*(xp_linearization-mu{i});
end
for i = 1:M
    a(:,i) = -1*alpha(i)*(invP{i}*(xp_linearization-mu{i}) + common_term_a) ;
end

%% likelihood parameters calculation
z_dim = length(zc);

likelihoodGMM = ps.likeparams.gmmNoiseModel;
N = likelihoodGMM.NumComponents;
zeta = cell(N,1);
invR = cell(N,1);
%gamma = likelihoodGMM.ComponentProportion;
for j = 1:N
    zeta{j} = likelihoodGMM.mu(j,:)';
    invR{j} = inv(cov_regularize(likelihoodGMM.Sigma(:,:,j)));
end
beta = -1*posterior(likelihoodGMM,(zc -H*xp_linearization)');
b_dashed = zeros(z_dim,N);
common_term_b_dashed = zeros(z_dim,1);
for j = 1:N
     common_term_b_dashed = common_term_b_dashed + ...
         beta(j)*H'*invR{j}*(zc-H*xp_linearization-zeta{j});
end
for j = 1:N
    b_dashed(:,j) = beta(j)*(H'*invR{j}*(zc-H*xp_linearization-zeta{j}) ...
        +common_term_b_dashed);
end

%% calculation of Q, f, G and q 

Q1 = zeros(x_dim);
s = zeros(1,x_dim);
U = zeros(x_dim);
for i = 1:M
    Q1 = Q1 + ...
    ((alpha(i) - xp_linearization'*a(:,i))*eye(x_dim) + a(:,i)*(xp_linearization-mu{i})')*invP{i};
    s = s - (alpha(i) - xp_linearization'*a(:,i))*mu{i}'*invP{i};
    U = U + (alpha(i) - xp_linearization'*a(:,i))*invP{i} + invP{i}*mu{i}*a(:,i)';
end

Q2 = zeros(x_dim);
m = zeros(1,x_dim);
K_dashed = zeros(x_dim);
q = zeros(1,x_dim);
for j = 1:N
    Q2 = Q2 + ((beta(j)-xp_linearization'*b_dashed(:,j))*H' -...
       b_dashed(:,j)'*(zc-H*xp_linearization-zeta{j}))*invR{j}*H;
     m = m - (beta(j)-xp_linearization'*b_dashed(:,j))*(zc-zeta{j})'*invR{j}*H;
     K_dashed = K_dashed + H'*invR{j}*((beta(j)-xp_linearization'*b_dashed(:,j))*...
         H  - (zc-zeta{j})*b_dashed(:,j)');
     q = q - beta(j)*(zc-H*xp_linearization-zeta{j})'*invR{j}*H;
end
Q = Q1 + lambda*Q2;
f = s + lambda*m;
G = U + lambda*K_dashed;
%invG = inv(G);

%% calculation of E, A_GMM and b_GMM
% AA = zeros(x_dim,x_dim,M+N+1);
% BB = zeros(x_dim,x_dim,M+N+1);
% FF = zeros(x_dim,x_dim);
% 
% for l = 1:M+N+1
%  if l==1
%      AA(:,:,l) = Q;
%      BB(:,:,l) = eye(size(Q));
%  elseif l>=2 && l<=M+1
%      AA(:,:,l) = -a(:,l-1)*f;
%      BB(:,:,l) = invG*invP{l-1};
%      FF = FF + a(:,l-1)*q*invG*invP{l-1};
%  else
%      AA(:,:,l) = -lambda*b_dashed(:,l-M-1)*f;
%      BB(:,:,l) = invG*H'*invR{l-M-1}*H;
%      FF = FF + lambda*b_dashed(:,l-M-1)*q*invG*H'*invR{l-M-1}*H;
%  end    
% end
% 
% A_GMM = solveLinearMatrixEqn(AA,BB,FF);

E = kron(eye(x_dim),Q);
pvec = zeros(x_dim^2,1);
for i = 1:M
    E = E + kron((G\invP{i})',(a(:,i)*f));
    temp1 = a(:,i)*q*(G\invP{i});
    pvec = pvec + temp1(:);
end

for j = 1:N
    E = E + lambda*kron(((G\H')*invR{j}*H)',(b_dashed(:,j)*f));
    temp2 = lambda*b_dashed(:,j)*q*(G\H')*invR{j}*H;
    pvec = pvec + temp2(:);
end
preCondtioningMatrix = cov_regularize(diag(sum(abs(E))));
intermediateSoln  = (E/preCondtioningMatrix)\pvec;
Avec = preCondtioningMatrix\intermediateSoln;

%Avec = E\pvec;
A_GMM = reshape(Avec,x_dim,x_dim);
b_GMM = -1*(f*A_GMM+q)/G;
end