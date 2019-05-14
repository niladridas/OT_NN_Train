% Author: Vedang Deshpande
% Date: 13th May 2019
% Inputs: 
% xMean - mean of random vector x
% Px - covariance of x
% a,k,b - parameters for scaling
% Outputs: 
% xSigmaPts - sigma pts of x (each column is a sigma pt)
% Ref: Eq.15 of https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf
function [xSigmaPts,Wm,Wc] = getSigmaPts(xMean,Px,a,k,b)
    L = length(xMean); % dim of x
    nPt = 2*L+1; % number of sigma pts
    lambda = a^2*(L+k)-L;
    TEMP = sqrtm((L+lambda)*Px)'; % Mind the transpose
    
    % Sigma Pts
    xSigmaPts = zeros(L,nPt);
    xSigmaPts(:,1) = xMean;
    for i = 1:L
        xSigmaPts(:,i+1) = xMean + TEMP(:,i);
        xSigmaPts(:,i+1+L) = xMean - TEMP(:,i);
    end
    
    % Weights
    Wm = 0.5*ones(1,nPt)/(L+lambda); Wm(1,1) = lambda/(L+lambda);
    Wc = Wm; Wc(1,1) = lambda/(L+lambda)+(1-a^2+b);
end