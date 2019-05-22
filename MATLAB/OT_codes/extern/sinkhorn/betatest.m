clc;clear;close;
M = 50;
a = 0;
b = 1;
alp = 2;
bet = 5;
B = gamma(alp)*gamma(bet)/gamma(alp+bet);
f = @(x) x.^(alp-1).*(1-x).^(bet-1)/B;
uni_samples = rand(M,1);
% uni_samples = zeros(M,1);
% for i = 1:M
%     uni_samples(i) = a+(b-a)*i/(M+1);
% end
weight_gen = arrayfun(f,uni_samples);
w = weight_gen/sum(weight_gen) ;
% Cost Matrix
D = abs(repmat(uni_samples,1,M)-repmat(uni_samples',M,1))*10;
% set lambda
lambda=100;
% the matrix to be scaled.
K=exp(-lambda*D);
% in practical situations it might be a good idea to do the following:
%K(K<1e-100)=1e-100;
% pre-compute matrix U, the Schur product of K and M.
U=K.*D;
% draw and normalize 1 point in the simplex with a few zeros (this is not a uniform sampling)
r=ones(M,1)/M;
% draw and normalize N points in the simplex with a few zeros (not uniform)
c=w;
[Dist,lowerEMD,l,m]=sinkhornTransport(r,c,K,U,lambda,'marginalDifference',...
    1,05e-7,5000,0); % running with VERBOSE
T=bsxfun(@times,m(:,1)',(bsxfun(@times,l(:,1),K))); % this is the optimal transport.
Pmat = T*M;
x_a = zeros(M,1);
for i= 1:M
    for j = 1:M
        x_a(i) = sum(Pmat(i,:)*uni_samples);
    end
end
h = histogram(x_a,10);%,'Normalization','probability');
hold on;
plot(0:0.05:1,f(0:0.05:1));