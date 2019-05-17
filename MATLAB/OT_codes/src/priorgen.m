function [Xfinal,x_truefinalcirc,z_circ,W_post,H,Sigmeas] = priorgen(nSamp)
d2r = pi/180;
%% LEO Parameters
Re = 6378.1363; %km radius of earth
h = 1000; % height of satellite km
R0 = (Re + h)/Re; % Normalized
th0 = 0*d2r;
V = 7.35; %km/s speed of satellite
Tp = 105*60; % Orbit Period in sec
Vb = Re/Tp;   % Velocity normalizing factor.
flagJ2 = 1;
rd = 0;  % normalized rate
thd = (V/Vb)/R0; % rad/s
% Gaussian uncertainty
mu = [R0; th0];
Sig = diag([0.0001*h/Re;0.2*d2r]);
% rng default  % For reproducibility
X = mvnrnd(mu,Sig,nSamp);
Rd = ones(nSamp,1)*rd;
Thd = ones(nSamp,1)*thd;
Xinit = [X(:,1) Rd X(:,2) Thd];
%% Final truth data after propagation 
x_true = [R0;0; th0;thd];
tf = 1;
[t1,x1] = ode45(@satelliteDynamics,[0 tf],x_true,[],flagJ2,Tp,Re);
x_truefinalcirc = x1(end,:);
% x1 = x_truefinalcirc(1,1)*cos(x_truefinalcirc(1,3));
% y1 = x_truefinalcirc(1,1)*sin(x_truefinalcirc(1,3));
% xv1 = x_truefinalcirc(1,2)*cos(x_truefinalcirc(1,3))-x_truefinalcirc(1,1)*sin(x_truefinalcirc(1,3))*x_truefinalcirc(1,4);
% yv1 = x_truefinalcirc(1,2)*sin(x_truefinalcirc(1,3))+x_truefinalcirc(1,1)*cos(x_truefinalcirc(1,3))*x_truefinalcirc(1,4);
% x_truefinalcar = [x1,y1,xv1,yv1];
%% Propagating the initial samples
Xfinal = zeros(nSamp,4);
for i=1:nSamp
    x0 = Xinit(i,:)';
    [t,x] = ode45(@satelliteDynamics,[0 tf],x0,[],flagJ2,Tp,Re);
    Xfinal(i,:) = x(end,:);
end
% % Propagated samples in XYVxVy
% xpos = Xfinal(:,1).*cos(Xfinal(:,3));
% ypos = Xfinal(:,1).*sin(Xfinal(:,3));
% xvel = Xfinal(:,2).*cos(Xfinal(:,3))-Xfinal(:,1).*sin(Xfinal(:,3)).*Xfinal(:,4);
% yvel = Xfinal(:,2).*sin(Xfinal(:,3))+Xfinal(:,1).*cos(Xfinal(:,3)).*Xfinal(:,4);
% X_posvel = [xpos,ypos,xvel,yvel];
%% Measurement value (only R and Theta measurements)
z_nonoise = [x_truefinalcirc(1,1);x_truefinalcirc(1,3)];
Sigmeas = diag([0.00001*h/Re;0.05*d2r]);
meas_noise = mvnrnd([0;0],Sigmeas,1);
z_circ = z_nonoise + meas_noise';
% z_cart = [z_circ(1,1)*cos(z_circ(2,1));z_circ(1,1)*sin(z_circ(2,1))];
%% Calculate weight of each samples based on likelihood function and the measurement
H = [1 0 0 0;
     0 0 1 0];
W_post = zeros(1,nSamp); 
for i=1:nSamp
    y = z_circ;
    W_post(1,i)  = exp(-0.5*(y-H*Xfinal(i,:)')'*(Sigmeas\(y-H*Xfinal(i,:)'))); % Weights from likelihood.
end
W_post = W_post./(sum(W_post));
%% Propagation Dynamics
function xdot = satelliteDynamics(~,x,J2Flag,Tp,Re)
    mu = 398600.4415; %(km^3/sec^2)
    mub = mu*Tp^2/Re^3; % Normalized
    J2 = 1.7555*10^10;
    J2b = J2*Tp^2/Re^5;
    % J2 effect
    if J2Flag
        r = x(1);
        th = x(3);
        Fr = J2b*(3/2)*(3*sin(th)^2-1)/r^4;
        Fth = -J2b*3*cos(th)*sin(th)/r^4;
    else
        Fr = 0;
        Fth = 0;
    end
    xdot(1,1) = x(2);
    xdot(2,1) = -mub/x(1)^2 + x(1)*x(4)^2 + Fr;
    xdot(3,1) = x(4);
    xdot(4,1) = -2*x(2)*x(4)/x(1) + Fth;
end
end