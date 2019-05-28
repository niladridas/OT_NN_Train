% Author: Vedang Deshpande
% Date: 13th May 2019
% Inputs: 
% xMean - mean of random vector x
% Px - covariance of x
% a,k,b - parameters for scaling
% Outputs: 
% xSigmaPts - sigma pts of x (each column is a sigma pt)
% Ref1: Eq.15 of https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf
% Ref2: Eq.5-Eq.9 https://www.cse.sc.edu/~terejanu/files/tutorialUKF.pdf
function [xSigmaPts,Wm,Wc] = getSigmaPts(xMean,Px,a,k,b)
    L = length(xMean); % dim of x
    nPt = 2*L+1; % number of sigma pts
    
    
    % Ref1: Sigma Pts 
    lambda = a^2*(L+k)-L;
%     TEMP = sqrtm((L+lambda)*Px)'; % Mind the transpose
    TEMP = sqrt(L+lambda)*chol(Px)';
    xSigmaPts = zeros(L,nPt);
    xSigmaPts(:,1) = xMean;
    for i = 1:L
        xSigmaPts(:,i+1) = xMean + TEMP(:,i);
        xSigmaPts(:,i+1+L) = xMean - TEMP(:,i);
    end
    % Ref1: Weights
    Wm = 0.5*ones(1,nPt)/(L+lambda); Wm(1,1) = lambda/(L+lambda);
    Wc = Wm; Wc(1,1) = lambda/(L+lambda)+(1-a^2+b);

    
    % Ref2: Sigma Pts
%     W0 = 2*rand-1;
%     xSigmaPts = zeros(L,nPt);
%     xSigmaPts(:,1) = xMean;
%     TEMP = sqrtm((L/(1-W0))*Px)';
% %     isreal(TEMP)
%     for i = 1:L
%         xSigmaPts(:,i+1) = xMean + TEMP(:,i);
%         xSigmaPts(:,i+1+L) = xMean - TEMP(:,i);
%     end
%     % Ref2: Weights
%     Wm = (1-W0)*ones(1,nPt)/(2*L); Wm(1,1) = W0;
%     Wc=Wm;
    
end