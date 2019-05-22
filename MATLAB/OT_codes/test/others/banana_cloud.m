clear; clc;
d2r = pi/180;

%% LEO Parameters
Re = 6378.1363; %km radius of earth
h = 1000; % height of satellite km
R0 = (Re + h)/Re; % Normalized
th0 = 0*d2r;
V = 7.35; %km/s speed of satellite
Tp = 105*60; % Orbit Period in sec
Vb = Re/Tp;   % Velocity normalizing factor.
flagJ2 = 0;

%%
rd = 0;  % normalized rate
thd = (V/Vb)/R0; % rad/s

%% Initial Condition Uncertainty
nSamp = 500;
% Gaussian uncertainty
mu = [R0; th0];
Sig = diag([0.00001*h/Re;0.05*d2r]);
rng default  % For reproducibility
X = mvnrnd(mu,Sig,nSamp);

fprintf('Range Rmin = %f km, thmin = %f deg \n',min((X(:,1)-1)*Re),min(X(:,2))/d2r);
fprintf('Range Rmax = %f km, thmax = %f deg \n',max((X(:,1)-1)*Re),max(X(:,2))/d2r);

%% MC Simulation
%T = [0 0.001 1 5 10 15]; % Orbit interval
T = [0 .000001 1];% 5 10 20 50];
dT = diff(T);
disp('Starting Monte-Carlo ...');
%X0 = X;
% figure(1); clf; hold on;

Rd = ones(nSamp,1)*rd;
Thd = ones(nSamp,1)*thd;
Xt = [X(:,1) Rd X(:,2) Thd];
clear r1 th1

%% Thruth data
x_true = [R0;0; th0;thd];
[t1,x1] = ode45(@satelliteDynamics,[0 1],x_true,[],flagJ2,Tp,Re);

x_actual = x1(end,:);

%%
for k=1:length(dT)
    fprintf('Simulating nOrbit = %f',T(k+1));
    tic;
    for i=1:nSamp
        x0 = Xt(i,:)';
        [t,x] = ode45(@satelliteDynamics,[0 dT(k)],x0,[],flagJ2,Tp,Re);
        Xt(i,:) = x(end,:);
    end
    toc;
%     subplot(2,3,k);
%     polar(Xt(:,3),Xt(:,1),'r.');
%     title(sprintf('Time = %.2f orbits (%.1f hr)',T(k+1),Tp*T(k+1)/3600));
end
disp('Done simulating ...');

% The last Xt is saved and plotted on XY plane and enlarged
X_pos = Xt(:,1).*cos(Xt(:,3));%rcostheta
Y_pos = Xt(:,1).*sin(Xt(:,3));%rcostheta
Y = [Xt(:,1)';Xt(:,3)'];
% figure(1);clf;
% scatter(X_pos,Y_pos,'filled');
r_mean = mean(Xt(:,1));
theta_mean = mean(Xt(:,3));

mu3 = [r_mean+0.1*sqrt(0.0001*h/Re);theta_mean+0.1*sqrt(0.0*d2r)];
Sig3 = diag([0.00001*h/Re;0.05*d2r]);


Y_enkf = Y;
% EnKF
EnKF_samples = enkf_samples(Y_enkf, mu3, Y_enkf, Sig3);

% Likelihood
parfor i=1:nSamp
    y = Y(:,i);
    W(1,i)  = exp(-0.5*(y-mu3)'*inv(Sig3)*(y-mu3).*3); % Weights from likelihood.
end
% COMMENTS: why *.3 at the end of exp(-0.5*(y-mu3)'*inv(Sig3)*(y-mu3)*.3)

%% Plot
pointSize = 10;
figure(1); clf;
subplot(2,2,1); scatter(Y(1,:).*cos(Y(2,:)),Y(1,:).*sin(Y(2,:)),pointSize,'filled'); 
xlabel('X');ylabel('Y');
title('Equally weighted prior');grid on;

subplot(2,2,2); scatter(Y(1,:).*cos(Y(2,:)),Y(1,:).*sin(Y(2,:)),pointSize,W/sum(W),'filled'); %colorbar;colormap jet
title('Weighted posterior'); grid on;
xlabel('X');ylabel('Y');
ax = axis;

% Determine OT Map 
W0 = ones(1,size(Y,2))/size(Y,2);
tic
OT_samplesX = eq_wsamples(Y,W0,W/sum(W));
toc
subplot(2,2,3); scatter(OT_samplesX(1,:).*cos(OT_samplesX(2,:)),OT_samplesX(1,:).*sin(OT_samplesX(2,:)),pointSize,'filled','green');hold on;
% Draw the mean 
ot_mean_x = mean(OT_samplesX(1,:).*cos(OT_samplesX(2,:)));
ot_mean_y = mean(OT_samplesX(1,:).*sin(OT_samplesX(2,:)));
plot(ot_mean_x, ot_mean_y, '.r', 'MarkerSize',30);hold on;
plot(x_actual(1)*cos(x_actual(3)), x_actual(1)*sin(x_actual(3)), '.g', 'MarkerSize',30);
title('Equally weighted posterior for OT');grid on;
xlabel('X');ylabel('Y');
axis(ax);
subplot(2,2,4); scatter(EnKF_samples(1,:).*cos(EnKF_samples(2,:)),EnKF_samples(1,:).*sin(EnKF_samples(2,:)),pointSize,'filled','blue');hold on
% Draw the mean 
enkf_mean_x = mean(EnKF_samples(1,:).*cos(EnKF_samples(2,:)));
enkf_mean_y = mean(EnKF_samples(1,:).*sin(EnKF_samples(2,:)));
plot(enkf_mean_x, enkf_mean_y, '.r', 'MarkerSize',30);hold on;
plot(x_actual(1)*cos(x_actual(3)), x_actual(1)*sin(x_actual(3)), '.r', 'MarkerSize',30);
title('Equally weighted posterior for EnKF');grid on;
xlabel('X');ylabel('Y');
axis(ax);


%% ALL FUNCTIONS USED WRITTEN BELOW
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