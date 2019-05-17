function [slope_deterministic,slope,logJacobianDet] = calculateSlope_PFPF_LEDH...
    (z,xp_deterministic,xp,mu_0,P,R,zeta,ps,lambda,step_size)
% computes the update for each particle
% using individual gradients at each particle location
%
% Inputs:
% z: measurement
% vg: a struct that contains the filter output
% ps: a struct with filter and simulation parameters
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
% log_jacobian_det: the log determinant of |I+step_size A^i|, used later in
% weight update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dim,N] = size(xp_deterministic);
zdim = size(z,1);


xp_linearization = xp_deterministic;

% Calculate the error term due to linearization
% and subtract it from the measurement vector
% This results in a [zdim,N] matrix, since the linearization and measurement evaluation is at each
% particle location.
%
H = ps.likeparams.dh_dx_func(xp_linearization,ps.likeparams);  % dh/dx at particle locations zdim x dim x N particles
h = ps.likeparams.h_func(xp_linearization,ps.likeparams);
e = zeros(size(h));
for particle_ix = 1:size(h,2)
    e(:,particle_ix) = h(:,particle_ix)-H(:,:,particle_ix)*xp_linearization(:,particle_ix);
end

zc = bsxfun(@minus,z,e+zeta);
% Reshape to a [zdim,1,N] matrix to facilitate later computation
zc = reshape(zc,[zdim,1,size(zc,2)]);


%%
A = zeros(dim,dim,size(H,3));
b = zeros(dim,size(H,3));

slope_deterministic = zeros(dim,N);
slope = zeros(dim,N);
logJacobianDet = zeros(1,N);

for particle_ix = 1:size(H,3)
    Hi = squeeze(H(:,:,particle_ix));
    
    A_i = -0.5*P(:,:,particle_ix)*Hi'*((lambda*Hi*P(:,:,particle_ix)*Hi'+R(:,:,particle_ix))\Hi);
    A(:,:,particle_ix) = A_i;
    b(:,particle_ix) = (eye(dim)+2*lambda*A_i)...
        *((eye(dim)+lambda*A_i)*P(:,:,particle_ix)*Hi'...
        *(R(:,:,particle_ix)\zc(:,1,particle_ix))+A_i*mu_0(:,particle_ix));
    
   
    slope_deterministic(:,particle_ix) = A_i*xp_deterministic(:,particle_ix) + b(:,particle_ix);
    slope(:,particle_ix) = A_i*xp(:,particle_ix) + b(:,particle_ix);
    logJacobianDet(particle_ix) = log(abs(det(eye(dim)+step_size*A_i)));
end
end