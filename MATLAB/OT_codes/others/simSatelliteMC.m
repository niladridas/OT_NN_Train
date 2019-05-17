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

rd = 0;  % normalized rate
thd = (V/Vb)/R0; % rad/s


%% Initial Condition Uncertainty
nSamp = 300;

% Gaussian uncertainty
% mu = [R0; th0];
% Sig = diag([0.01*h/Re;1*d2r]);
% rng default  % For reproducibility
% X = mvnrnd(mu,Sig,nSamp);

% Uniform Uncertainty
% -- Generate uniform in [0,1]
X = rand(nSamp,2);
eR = 5/Re; % Error is about 1.5 km
eTh = 10*d2r;  % +- 5 degree --> Angle errors are pretty bad.
% Scale uniform distribution in normalized r and theta
X(:,1) = R0;%+(-1 + 2*X(:,1))*eR;
X(:,2) = (-1 + 2*X(:,2))*eTh;

disp(sprintf('Range Rmin = %f km, thmin = %f deg',min((X(:,1)-1)*Re),min(X(:,2))/d2r));
disp(sprintf('Range Rmax = %f km, thmax = %f deg',max((X(:,1)-1)*Re),max(X(:,2))/d2r));

%% MC Simulation
%T = [0 0.001 1 5 10 15]; % Orbit interval
T = [0 .000001 1 5 10 20 50];
dT = diff(T);
disp('Starting Monte-Carlo ...');
%X0 = X;
figure(1); clf; hold on;

Rd = ones(nSamp,1)*rd;
Thd = ones(nSamp,1)*thd;
Xt = [X(:,1) Rd X(:,2) Thd];
clear r1 th1

for k=1:length(dT)
    disp(sprintf('Simulating nOrbit = %f',T(k+1)));
    tic;
    for i=1:nSamp
        x0 = Xt(i,:)';
        [t,x] = ode45(@satelliteDynamics,[0 dT(k)],x0,[],flagJ2,Tp,Re);
        Xt(i,:) = x(end,:);
        r1(i,1) = x(end,1);
        th1(i,1) = x(end,3);
    end
    toc;
    
    subplot(2,3,k);
    polar(Xt(:,3),Xt(:,1),'r.');
    title(sprintf('Time = %.2f orbits (%.1f hr)',T(k+1),Tp*T(k+1)/3600));
end
disp('Done simulating ...');

%% Sensor Model
% There are four states x = (r,rdot,th,thdot). We are interested in tracking
% (r,th) := (x(1),x(3)). Typically we can sense y = (r,th).
% Typically the errors in r are between plus minus 1.5/Re
% Errors in th are really bad, we can assume +- 10deg. 
% 
% Not sure how to convert them to Gaussian sensor noise. Maybe we
% can take the min max to be +- 3sigma and construct a noise model.
% 
% The objective is to keep estimate consistent with large sensing gaps and
% large sensor noise.
%