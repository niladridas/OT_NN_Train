function [slope] = calculateSlopeComponentwise_EDH(z,xp,mu_0,P,R,zeta,ps,lambda)
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

[dim,N] = size(xp);
zdim = size(z,1);


xp_linearization = mean(xp,2);

% Calculate the error term due to linearization
% and subtract it from the measurement vector

%
H = ps.likeparams.dh_dx_func(xp_linearization,ps.likeparams);  % dh/dx at particle locations zdim x dim x N particles
h = ps.likeparams.h_func(xp_linearization,ps.likeparams);
%e = zeros(size(h));
e = h-H*xp_linearization;

zc = z-zeta-e;

%%
A = -0.5*P*H'*((lambda*H*P*H'+R)\H);
b = (eye(dim)+2*lambda*A)...
        *((eye(dim)+lambda*A)*P*H'...
        *(R\zc)+A*mu_0);
%%
slope = A*xp + repmat(b,1,size(xp,2));

end