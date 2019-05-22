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

J2flag = 1;
[t1,y1] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);
r = y1(:,1);
th = y1(:,3);
x0 = y1(end,:);
J2flag = 0;
[t2,y2] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);

% x0 = [R0;Rd;th0;thd];
% J2flag = 1;
% [t3,y3] = ode45(@satelliteDynamics,[0 2*nOrbits],x0,[],J2flag,Tp,Re);
% figure(1); clf; hold on;
% plot(t1,y1(:,1),'r',t2+t1(end),y2(:,1),'b',t3,y3(:,1),'k-.')


% figure(2); clf;
% J2flag = 0;
% x0 = [R0;Rd;th0;thd];
% [t1,y1] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);
% r = y1(:,1);
% th = y1(:,3);
% subplot(2,1,1); hold on; plot(t1,r,'b');
% subplot(2,1,2); hold on; plot(t1,th/d2r,'b');
% figure(3); clf; polar(th,r,'b');
% 
% J2flag = 1;
% [t,y] = ode45(@satelliteDynamics,[0 nOrbits],x0,[],J2flag,Tp,Re);
% r = y(:,1);
% th = y(:,3);
% 
% figure(2); 
% subplot(2,1,1); hold on; plot(t,r,'r');
% subplot(2,1,2); hold on; plot(t,th/d2r,'r');
% 
% figure(3); hold on; polar(th,r,'r--');
% break;

%% Initial Condition Uncertainty
% nSamp = 300;
% 
% % Gaussian uncertainty
% mu = [R0; th0];
% Sig = diag([0.01*h/Re;1*d2r]);
% rng default  % For reproducibility
% X = mvnrnd(mu,Sig,nSamp);
% 
% % Generate uniform in [0,1]
% X = rand(nSamp,2);
% eR = 0.10; % +- 5 percent of h
% eTh = 1;  % +- 1 degree
% % Scale uniform distribution in normalized r and theta
% X(:,1) = R0+(-1 + 2*X(:,1))*(eR*h/Re);
% X(:,2) = (-1 + 2*X(:,2))*eTh*d2r;
% 
% disp(sprintf('Range Rmin = %f km, thmin = %f deg',min((X(:,1)-1)*Re),min(X(:,2))/d2r));
% disp(sprintf('Range Rmax = %f km, thmax = %f deg',max((X(:,1)-1)*Re),max(X(:,2))/d2r));
% 
% %% MC Simulation
% T = [0 0.001 1 5 10 15 50]; % Orbit interval
% dT = diff(T);
% disp('Starting Monte-Carlo ...');
% clear th r
% X0 = X;
% figure(2); clf;
% for k=1:length(dT)
%     disp(sprintf('Simulating nOrbit = %f',T(k)));
%     for i=1:nSamp
%         x0 = [X0(i,1) Rd X0(i,2) thd]';
%         [t,y] = ode45(@satelliteDynamics,[0 dT(k)],x0,[],1,Tp,Re);
%         X0(i,1) = y(end,1); % r
%         X0(i,2) = y(end,2); % th
%     end
%     subplot(2,3,k);
%     polar(X0(:,2),X0(:,1),'r.');
%     title(sprintf('Time = %.2f orbits',T(k+1)));
%     error.r(k,:) = [min(X0(:,1))*Re max(X0(:,1))*Re];
%     error.th(k,:) = [min(X0(:,2)/d2r) max(X0(:,2)/d2r)];
% end
% disp('Done simulating ...');
% 
% % Plot
% figure(2); clf;
% for i=2:length(T)
%     subplot(2,3,i-1)
%     polar(th(:,i),r(:,i),'r.');
%     title(sprintf('Time = %.2f orbits',T(i)));
%     error.th(i,:) = [min(th(:,i)/d2r) max(th(:,i)/d2r)];
%     error.r(i,:) = [min(r(:,i))*Re max(r(:,i))*Re];
% end
% figure(3); clf;
% subplot(2,1,1);plot(T,error.th); legend('min','max'); xlabel('Orbit'); ylabel('th (deg)'); grid on;
% subplot(2,1,2);plot(T,error.r); legend('min','max'); xlabel('Orbit'); ylabel('r (km)'); grid on;
