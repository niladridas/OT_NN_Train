% plot error for low sensor noise and high sensor noise
clc;clear;
high = load('errorhigh');
low = load('errorlow');
load('sample_sizes');
figure;
pointSize = 25;

plot(sample_sizes,high.est_error(1,:),'r');hold on;
plot(sample_sizes,high.est_error(2,:),'b'); hold on;
plot(sample_sizes,low.est_error(1,:),'--r');hold on;
plot(sample_sizes,low.est_error(2,:),'--b'); hold on;

scatter(sample_sizes,high.est_error(1,:),pointSize,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
scatter(sample_sizes,high.est_error(2,:),pointSize,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
scatter(sample_sizes,low.est_error(1,:),pointSize,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
scatter(sample_sizes,low.est_error(2,:),pointSize,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');hold on;



xlabel('sample size');ylabel('Estimation error');
set(gca,'fontsize',10,'fontweight','bold');
title('Estimation Error vs Sample size');
legend('OT filter error (high R)','Daum filter error (high R)','OT filter error (low R)','Daum filter error (low R)');
export_fig(gcf,'FontMode','fixed','FontSize',25,'Color','Error.pdf','-nocrop');
hold off