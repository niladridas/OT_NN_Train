% Parameters to generate the real x and the observations
simparams.burnperiod = 5;
% simparams.delt = 1;
simparams.delt = 1;
simparams.Nt = 5;
simparams.d2r = pi/180;
%% LEO Parameters
simparams.Re = 6378.1363; %km radius of earth
simparams.h = 1000; % height of satellite km
simparams.Tp = 105*60; % Orbit Period in sec
simparams.flagJ2 = 1;
simparams.R0 = (simparams.Re + simparams.h)/simparams.Re; % Normalized
simparams.th0 = 0*simparams.d2r;
simparams.V = 7.35; %km/s speed of satellite
simparams.Vb = simparams.Re/simparams.Tp;   % Velocity normalizing factor.
simparams.rd = 0;  % normalized rate
simparams.thd = (simparams.V/simparams.Vb)/simparams.R0; % rad/s
simparams.H = [1 0 0 0;
         0 0 1 0];
simparams.meas_sig = diag([0.00001*simparams.h/simparams.Re;0.05*simparams.d2r]);
% simparams.meas_sig = diag([0.001*simparams.h/simparams.Re;0.1*simparams.d2r]);
simparams.dyn = @(t,x) satelliteDynamics(t,x,simparams.flagJ2,simparams.Tp,simparams.Re);
simparams.likfun = @(y,x) likfun(y,x,simparams.meas_sig,simparams.H);

