function [slope] = calculateSlope_LEDH(z,xp,mu_0,P,R,lambda)
% computes the update for each particle
% using individual gradients at each particle location
%
% Inputs:
% z: measurement
% lambda: the particle flow time step
% step_size: a scalar that shows the step size.
%
% The computation is as follows (see paper for details)
%
% e = h(x)-Hx;
%
% A = -0.5*PP*H'*inv(lambda*H*PP*H'+R)*H;
% b = (I + 2*lambda*A)*((I+lambda*A)*PP*H'*inv(R)*(z-e) + A*mX)
%
% slope = Ax + b;
%
% Outputs:
% slope: a struct contains the field real, which includes slopes for vg.xp
%       (the original particles) and the field auxiliary_individual which contains
%        slopes for vg.xp_auxiliary_individual (the slopes for particles used to calculate the slope,
%        used when clustering is not enanled or clustering is performed in each intermediate time step.)
%        and the field auxiliary_cluster: slopes for cluster centroid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dim,N] = size(xp);
zdim = size(z,1);

xp_linearization = xp;

H = repmat(eye(dim),1,1,N);
h = xp;
e = zeros(size(h));
for particle_ix = 1:size(h,2)
    e(:,particle_ix) = h(:,particle_ix)-H(:,:,particle_ix)*xp_linearization(:,particle_ix);
end

zc = bsxfun(@minus,z,e);
% Reshape to a [zdim,1,N] matrix to facilitate later computation
zc = reshape(zc,[zdim,1,size(zc,2)]);
%%
A = zeros(dim,dim,size(H,3));
b = zeros(dim,size(H,3));

slope = zeros(dim,N);

for particle_ix = 1:size(H,3)
    Hi = squeeze(H(:,:,particle_ix));
    
    A_i = -0.5*P*Hi'*((lambda*Hi*P*Hi'+R)\Hi);
    A(:,:,particle_ix) = A_i;
    b(:,particle_ix) = (eye(dim)+2*lambda*A_i)...
        *((eye(dim)+lambda*A_i)*P*Hi'...
        *(R\zc(:,1,particle_ix))+A_i*mu_0);
    slope(:,particle_ix) = A_i*xp(:,particle_ix) + b(:,particle_ix);
end
end