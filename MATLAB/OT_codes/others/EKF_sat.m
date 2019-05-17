% Demonstration of 1st Order EKF on noisy orbital data
%
% A non-linear satelite model is provided in 'sateliteDynamics.m'
%
% Nildri Das, Department of Aerospace Engineering

clc;
clear all;
disp('Filtering orbital signal with EKF...');
addpath('C:\Users\niladridas\Documents\Kalman\ekfukf');
save_plots = 1;


d2r = pi/180;

%% Single Initial Condition
% LEO Parameters
Re = 6378.1363; %km radius of earth
h = 900; % height of satellite km
R0 = (Re + h)/Re; % Normalized
th0 = 5*d2r;
V = 7.4; %km/s speed of satellite
Tp = 102.8*60; % Orbit Period in sec
Vb = Re/Tp;   % Velocity normalizing factor.

Rd = 0;  % normalized rate
thd = (V/Vb)/R0; %rad/s

% Initial Conditions
x0 = [R0;Rd;th0;thd];



  
% Spectral power density of the white noise.
q1 = 0.02;
q2 = 0.001;
q3 = 0.004;
q4 = 0.2;

Qc = diag([q1 q2 q3 q4]);
  
% Generate the real signal.
% num_orbit : total number of orbits
% num_samp : number of samples per orbit
 
num_orbit = 5;

[t1,y1] = generate_orbit(num_orbit);

num_samples = 4;

% Generating equispaced orbital sample time instants
t_sample = linspace(0,num_orbit,num_orbit*num_samples);

samplePoints = {t1, 1:(size(y1,2)+1)};
F = griddedInterpolant(samplePoints,[t1,y1]);
queryPoints = {t_sample,1:(size(y1,2)+1)};
Vq = F(queryPoints);

% Generate real sample points with noise

X = Vq(:,2:end) + gauss_rnd([0 0 0 0]', Qc,size(t_sample,2))';
  
% % Generate the observations with Gaussian noise.
sd1 = 0.04;
sd2 = 0.4;
R = diag([sd1^2 sd2^2]);

Y = zeros(1,size(t_sample,2));
Y_real = X(:,[1,3]);     
Y = Y_real + gauss_rnd([0 0]',R,size(t_sample,2))';

figure(1)
plot(1:size(t_sample,2),Y,'.',1:size(t_sample,2),Y_real)
  
% Initial guesses for the state mean and covariance.
M = x0;
P = diag([3 3 3 3]);    
%   
% Reserve space for estimates.
MM = zeros(size(M,1),size(Y,2));
PP = zeros(size(M,1),size(M,1),size(Y,2));

J2Flag = 0 ;
% Estimate with EKF
for k=1:size(Y,1)
   % The LEO satelite model is non-linear 
   [A,H] = lin_param_sat(M,J2Flag);
   [M,P] = ekf_predict1(M,P,A,Qc);
   [M,P] = ekf_update1(M,P,Y(k,:)',H,R);
   MM(:,k)   = M;
   PP(:,:,k) = P;
end

figure(2)
% Project the estimates to measurement space
plot(t_sample,MM([1,3],:)','o',t_sample,Y,'*', t_sample, Vq(:,[2,4]))
xlabel('Orbital time');
legend('Estimated Radius','Estimated Theta','Measured Radius','Measured Theta','Actual Radius','Actual Theta')
title('Implementation of EKF on satellite model')
grid on
figure(3)
polar(y1(:,3),y1(:,1),'r--')