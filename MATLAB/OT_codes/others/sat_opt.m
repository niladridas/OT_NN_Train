clear; clc;
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
x0 = [R0;Rd;th0;thd];
nOrbits = 5;

J2flag = 0;
[t1,y1] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);
r = y1(:,1);
th = y1(:,3);
x0 = y1(end,:);


bx1 = y1(10,1);
bx2 = y1(10,2);
bx3 = y1(10,3);
bx4 = y1(10,4);

mu = 398600.4415; %(km^3/sec^2)
mub = mu*Tp^2/Re^3; % Normalized


A = [0 1 0 0;
     bx4^2+(2*mub/(bx1^3)) 0 0 2*bx4*bx1;
     0 0 0 1;
     (2*bx2*bx4/(bx1^2)) (-2*bx4/bx1) 0 -(2*bx2/bx1)];
 
C = eye(4);
sig_Q = 0.1*[bx1;bx2;bx3;bx4];
% Q = diag([0.0130;0.05;0.0685;0.3933]);
Q = diag(sig_Q.*sig_Q);
P0 = 2*Q;
t = 0.01; % propagation and then update
P1_pr = expm(A*t)*P0*(expm(A*t))'+Q;


cvx_solver mosek

cvx_begin sdp
    % Variable declaration
    variable X1(4,4) symmetric 
    variable lambda(4,1)
    variable X2(4,4) symmetric
    variable K(4,4)
    % Main Optimization problem
    minimize(norm(lambda,1))
    % Variable defining

    trace(X1+X2)<=0.1*trace(P1_pr);
    
    [X1 (eye(4)-K*C);
    (eye(4)-K*C)' inv(P1_pr)]>=0;

    [X2 K;
     K' diag(lambda)] >= 0
cvx_end

P_f = (eye(4)-K*C)*P1_pr;
trace(P_f)
 