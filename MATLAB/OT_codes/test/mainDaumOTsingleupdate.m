% Author: Niladri Das
% Date: Aug 9, 2018
% About: Daum-Huang filter update is compared with OT based update.
%% TO-DO
% 1. Real state gen can be done just one time outside the loops
clear; clc;close all;
sample_sizes = 100:100:400;
nRep = 3; % repitition of the experiment with 'nRep' different sets of samples each of size 'nSamp'
run_time = zeros(2,length(sample_sizes),nRep);
est_error = zeros(2,length(sample_sizes),nRep);
rng('default');
seeds = 1:1:nRep; % Generate nRep seeds
PDF_Flag = 0; % Flag to produce pdf plots
Show_Plot = 0;
Show_Plotruntime = 0;
Show_Ploterror = 1;
plotfolderhd = dir('plots');
plotfolder = plotfolderhd.folder;
%%
for k = 1:length(sample_sizes)
    for m = 1:nRep
        nSamp = sample_sizes(k);
        rng(seeds(m));
        [Xfinal,x_truefinalcirc,z_circ,W_post,H,Sigmeas] = priorgen(nSamp);
        %% OT Filter
        W0 = ones(1,nSamp)/nSamp; % Prior is equally weighted
        tic; OT_samplesX = eq_wsamples(Xfinal',W0,W_post); toc;
        run_time(1,k,m) = toc;
        %% Duam Huang Filter
        P = cov(Xfinal);% Prior sample covariance
        mu_0 = mean(Xfinal)';% Prior sample mean
        R = Sigmeas;
        flowdyn = @(t,x) EDH(t,x,z_circ,H,mu_0,P,R);
        Daum_samplesX = zeros(4,nSamp);
        tic
        for i = 1:nSamp
            [tf,xf] = ode45(flowdyn,[0 1],Xfinal(i,:)');
            Daum_samplesX(:,i) = xf(end,:)';
        end
        toc
        run_time(2,k,m) = toc;  
        %% Plotting
        Xp = Xfinal;
        OTx = OT_samplesX' ;
        OTxbar = mean(OTx,1);
        Dx = Daum_samplesX' ;
        Dxbar = mean(Dx,1);
        xtrue = x_truefinalcirc;
        %% Plot only first instance of the repitions
        if Show_Plot == 1
            if m == 1
                % Plot prior samples
                pointSize = 20;
                figure(k); clf;
                scatter(Xfinal(:,1).*cos(Xfinal(:,3)),Xfinal(:,1).*sin(Xfinal(:,3)),pointSize,'filled'); hold on;
                plot(Xfinal(1,1)*cos(Xfinal(1,3)),Xfinal(1,1)*sin(Xfinal(1,3)),'.g', 'MarkerSize',30);
                xlabel('X');ylabel('Y');
                title('Prior');grid on;
                set(gca,'fontsize',10,'fontweight','bold');
                leg1txt = strcat('prior samples (',num2str(nSamp),')');
                legend(leg1txt,'true location');
                filename = strcat(plotfolder,'/prior.pdf');
                if PDF_Flag ==1;export_fig(gcf,'FontMode','fixed','FontSize',25,'Color',filename,'-nocrop');end
                hold off
                % Plot posterior of OT filtering
                pointSize = 20;
                figure(k+1); clf;
                ax1 = subplot(1,2,1);
                scatter(OTx(:,1).*cos(OTx(:,3)),OTx(:,1).*sin(OTx(:,3)),pointSize,'filled'); hold on;
                plot(xtrue(1,1)*cos(xtrue(1,3)),xtrue(1,1)*sin(xtrue(1,3)),'.r', 'MarkerSize',30);hold on;
                plot(OTxbar(1,1)*cos(OTxbar(1,3)),OTxbar(1,1)*sin(OTxbar(1,3)),'.g', 'MarkerSize',30);
                xlabel('X');ylabel('Y');
                title('OT Posterior');grid on;
                set(gca,'fontsize',10,'fontweight','bold');
                leg1txt = strcat('posterior samples (',num2str(nSamp),')');
                legend(leg1txt,'true location','estimated location');
                % Plot posterior of Daum Huang filtering
                ax2 = subplot(1,2,2);
                scatter(Dx(:,1).*cos(Dx(:,3)),Dx(:,1).*sin(Dx(:,3)),pointSize,'filled'); hold on;
                plot(xtrue(1,1)*cos(xtrue(1,3)),xtrue(1,1)*sin(xtrue(1,3)),'.r', 'MarkerSize',30);hold on;
                plot(Dxbar(1,1)*cos(Dxbar(1,3)),Dxbar(1,1)*sin(Dxbar(1,3)),'.g', 'MarkerSize',30);
                xlabel('X');ylabel('Y');
                title('Daum-Huang Posterior');grid on;
                set(gca,'fontsize',10,'fontweight','bold');
                legend(leg1txt,'true location','estimated location');
                linkaxes([ax1,ax2],'xy')
                filename = strcat(plotfolder,'/postsamples', num2str(nSamp), '.pdf');
                if PDF_Flag ==1;export_fig(gcf,'FontMode','fixed','FontSize',25,'Color',filename,'-nocrop');end
                hold off
            end
        end
        %%
        % Error in position
        truex = xtrue(1,1)*cos(xtrue(1,3));
        truey = xtrue(1,1)*sin(xtrue(1,3));
        xOT = OTxbar(1,1)*cos(OTxbar(1,3));
        yOT = OTxbar(1,1)*sin(OTxbar(1,3));
        xDaum = Dxbar(1,1)*cos(Dxbar(1,3));
        yDaum = Dxbar(1,1)*sin(Dxbar(1,3));
        Re = 1;%6378.1363;
        errorOT = norm([truex-xOT;truey-yOT])*Re;
        errorDaum = norm([truex-xDaum;truey-yDaum])*Re;
        est_error(:,k,m) = [errorOT;errorDaum];
    end
end
%% Run-time plot
run_time = mean(run_time,3);
if Show_Plotruntime == 1
    figure(length(sample_sizes)+1);
    pointSize = 25;
    scatter(sample_sizes,run_time(1,:),pointSize,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
    scatter(sample_sizes,run_time(2,:),pointSize,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
    plot(sample_sizes,run_time(1,:),'r');hold on;
    plot(sample_sizes,run_time(2,:),'b'); hold on;
    xlabel('sample size');ylabel('Run time in seconds');
    set(gca,'fontsize',10,'fontweight','bold');
    title('Run time vs Sample size');
    legend('OT filter run time','Daum filter run time');
    filename = strcat(plotfolder,'/runtime.pdf');
    if PDF_Flag ==1;export_fig(gcf,'FontMode','fixed','FontSize',25,'Color',filename,'-nocrop');end
    hold off
end

%% Error plot
avg_error = mean(est_error,3);
cov_error = zeros(2,length(sample_sizes));
for i = 1:length(sample_sizes)
    cov_error(1,i) = cov(reshape(est_error(1,i,:),nRep,1));
    cov_error(2,i) = cov(reshape(est_error(2,i,:),nRep,1));
end
est_error = avg_error;
figure(length(sample_sizes)+2);
if Show_Ploterror == 1
    pointSize = 25;
    scatter(sample_sizes,est_error(1,:),pointSize,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');hold on;
    scatter(sample_sizes,est_error(2,:),pointSize,'filled','MarkerEdgeColor','b','MarkerFaceColor','b');hold on;
    plot(sample_sizes,est_error(1,:),'r');hold on;
    plot(sample_sizes,est_error(2,:),'b'); hold on;
    xlabel('sample size');ylabel('Estimation error');
    set(gca,'fontsize',10,'fontweight','bold');
    title('Average Estimation Error vs Sample size');
    legend('OT filter error','Daum filter error');
    filename = strcat(plotfolder,'/Error.pdf.pdf');
    export_fig(gcf,'FontMode','fixed','FontSize',25,'Color',filename,'-nocrop');
    hold off
end


