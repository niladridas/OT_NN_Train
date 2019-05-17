function [t1,y1] = generate_orbit(nOrbits)

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

J2flag = 1;
[t1,y1] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);
end