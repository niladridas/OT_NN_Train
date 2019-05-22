filterparams.sample_sizes = 500;
filterparams.nRep = 1;
% Intial sample generation parameters
% Gaussian uncertainty
filterparams.mu = [simparams.R0; simparams.th0];
filterparams.Sig = diag([0.0001*simparams.h/simparams.Re;0.2*simparams.d2r]);
% MAP_Flag
% If MAP_Flag = 1 calculate MAP estimates in OT, Daum and EnKF
filterparams.MAP_Flag = 1;