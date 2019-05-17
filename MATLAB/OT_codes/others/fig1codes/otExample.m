clc; clear;
nsamp = 300;

mu1 = [-7;-5];
mu2 = [7;5];
mu3 = [4;1];
Sig1 = 2*eye(2,2);
Sig2 = 1*[2 .5;.5 3];
Sig3 = [0.5 0; 0 1];

% Priors
X = randn(2,nsamp);
Y1 = repmat(mu1,1,nsamp) + Sig1*X;
Y2 = repmat(mu2,1,nsamp) + Sig2*X;
Y = [Y1 Y2];

Y_enkf = Y;
% EnKF
EnKF_samples = enkf_samples(Y_enkf, mu3, Y_enkf, Sig3);


% Likelihood
parfor i=1:2*nsamp
    y = Y(:,i);
    W(1,i)  = exp(-0.5*(y-mu3)'*inv(Sig3)*(y-mu3)); % Weights from likelihood.
end
% COMMENTS: why *.3 at the end of exp(-0.5*(y-mu3)'*inv(Sig3)*(y-mu3)*.3)

%% Plot
pointSize = 10;
figure(1); clf;
subplot(1,3,1); scatter(Y(1,:),Y(2,:),pointSize,'filled'); 
title('Equally weighted prior');grid on;

subplot(1,3,2); scatter(Y(1,:),Y(2,:),pointSize,W/sum(W),'filled'); colorbar;colormap jet
title('Weighted posterior'); grid on;
ax = axis;

% Determine OT Map 
W0 = ones(1,size(Y,2))/size(Y,2);
tic
X = eq_wsamples(Y,W0,W/sum(W));
toc
subplot(1,3,3); scatter(X(1,:),X(2,:),pointSize,'filled','green'); hold on;
scatter(EnKF_samples(1,:),EnKF_samples(2,:),pointSize,'filled','blue');
legend('OT','EnKF');
title('Equally weighted posterior');grid on;
axis(ax);
% 
% print -depsc -r300 ot3.eps
% save otData X Y W W0 mu1 mu2 mu3 Sig1 Sig2 Sig3




