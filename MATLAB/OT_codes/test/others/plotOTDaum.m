clc;clear;
load('OTDaumdata');
Xp = OTDaumdata.Xfinal;
OTx = OTDaumdata.OT_samplesX' ;
OTxbar = mean(OTx,1);
Dx = OTDaumdata.Daum_samplesX' ;
Dxbar = mean(Dx,1);
xtrue = OTDaumdata.real_state;
% First plot for prior samples
%% Plot
pointSize = 20;
% figure(1); clf;
% scatter(Xp(:,1).*cos(Xp(:,3)),Xp(:,1).*sin(Xp(:,3)),pointSize,'filled'); 
% xlabel('X');ylabel('Y');
% title('Equally weighted prior');grid on;



figure(2); clf;

ax1 = subplot(1,2,1);
scatter(OTx(:,1).*cos(OTx(:,3)),OTx(:,1).*sin(OTx(:,3)),pointSize,'filled'); hold on;
plot(xtrue(1,1)*cos(xtrue(1,3)),xtrue(1,1)*sin(xtrue(1,3)),'.r', 'MarkerSize',30);hold on;
plot(OTxbar(1,1)*cos(OTxbar(1,3)),OTxbar(1,1)*sin(OTxbar(1,3)),'.g', 'MarkerSize',30);
xlabel('X');ylabel('Y');
title('OT Posterior');grid on;
set(gca,'fontsize',10,'fontweight','bold')

ax2 = subplot(1,2,2);
scatter(Dx(:,1).*cos(Dx(:,3)),Dx(:,1).*sin(Dx(:,3)),pointSize,'filled'); hold on;
plot(xtrue(1,1)*cos(xtrue(1,3)),xtrue(1,1)*sin(xtrue(1,3)),'.r', 'MarkerSize',30);hold on;
plot(Dxbar(1,1)*cos(Dxbar(1,3)),Dxbar(1,1)*sin(Dxbar(1,3)),'.g', 'MarkerSize',30);
xlabel('X');ylabel('Y');
title('Daum-Huang Posterior');grid on;
set(gca,'fontsize',10,'fontweight','bold')


linkaxes([ax1,ax2],'xy')

%% Error in position
truex = xtrue(1,1)*cos(xtrue(1,3));
truey = xtrue(1,1)*sin(xtrue(1,3));
xOT = OTxbar(1,1)*cos(OTxbar(1,3));
yOT = OTxbar(1,1)*sin(OTxbar(1,3));
xDaum = Dxbar(1,1)*cos(Dxbar(1,3));
yDaum = Dxbar(1,1)*sin(Dxbar(1,3));

Re = 6378.1363;
errorOT = norm([truex-xOT;truey-yOT])*Re
errorDaum = norm([truex-xDaum;truey-yDaum])*Re

%% Covariance
OTcovcart = cov([OTx(:,1).*cos(OTx(:,3)),OTx(:,1).*sin(OTx(:,3))]);
covxxOT = OTcovcart(1,1)*Re^2;
covyyOT = OTcovcart(2,2)*Re^2;
%
Daumcovcart = cov([Dx(:,1).*cos(Dx(:,3)),Dx(:,1).*sin(Dx(:,3))]);
covxxDaum = Daumcovcart(1,1)*Re^2;
covyyDaum = Daumcovcart(2,2)*Re^2;

