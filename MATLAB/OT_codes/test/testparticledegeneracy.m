% Author: Niladri Das
% About: Simulation showing how OT avoids the problem of particle
% degeneracy
clc;clear;close all;
prior_mean = 0;
prior_sig = 5;
lik_mean = 10;
lik_sig = 0.5;
%% PRIOR
x = -20:.1:20;
norm = normpdf(x,prior_mean,prior_sig);
figure;
% plot(x,norm);hold on;
area(x,norm);hold on;
%% LIKELIHOOD
x = -20:.1:20;
norm = normpdf(x,lik_mean,lik_sig);
% plot(x,norm);hold on;
area(x,norm);hold on;
%%
% % Sample from Gaussian centered around origin
rng('default');
rng(1);
N = 5;
Nsamples = normrnd(prior_mean,prior_sig,N,1);
scatter(Nsamples,zeros(N,1),'filled','MarkerFaceColor','black');
% % The likelihood function is also a Gaussian
%% OT
W_post = normpdf(Nsamples,10,1);
W_post = W_post./(sum(W_post));
W0 = ones(1,N)/N; % Prior is equally weighted
tic; OT_samplesX = eq_wsamples(Nsamples',W0,W_post'); toc;
%% DAUM-HUANG
P = prior_sig^2;%cov(Nsamples);% Prior sample covariance
mu_0 = 0;%mean(Nsamples);% Prior sample mean
R = 1;
y = lik_mean;
flowdyn = @(t,x) EDH(t,x,y,1,mu_0,P,R);
Daum_samplesX = zeros(1,N);
tic
for k = 1:N
    [~,xf] = ode23(flowdyn,[0 1],Nsamples(k,1));
    Daum_samplesX(:,k) = xf(end,:)';
end
toc
%% Posterior Samples
OT_mean = mean(OT_samplesX);
OT_cov = cov(OT_samplesX);

Daum_mean = mean(Daum_samplesX);
Daum_cov = cov(Daum_samplesX);
