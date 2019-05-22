% Author: Niladri Das
% Implementation of Exact Duam Huang filter
% Date:8 Aug 2018
function [slope] = EDH(lambda,xp,z,H,mu_0,P,R)
% computes gradient for individual particle location
% xp: particle locations
% z: measurement vector (zdim,1)
% lambda: particle flow time parameter
% 
% The computation is as follows
% A = -0.5*PP*H'*inv(lambda*H*PP*H'+R)*H;
% b = (I + 2*lambda*A)*((I+lambda*A)*PP*H'*inv(R)*(z-e) + A*mX)
%
% slope = Ax + b;

[dim,~] = size(xp); % N denoted the sample size

%%
A = -0.5*P*H'*((lambda*H*P*H'+R)\H);
b = (eye(dim)+2*lambda*A)...
        *((eye(dim)+lambda*A)*P*H'...
        *(R\z)+A*mu_0);
%%
slope = A*xp + b;
end