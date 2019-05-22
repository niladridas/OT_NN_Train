close all;clear;clc;
load('./data/fulldata');
%%
PDF_Flag = 1; % Flag to produce pdf plots
Show_Plotruntime = 0;
Show_Ploterror = 1;
plotfolderhd = dir('plots');
plotfolder = plotfolderhd.folder;
%% Data Preprocessing
% fulldata
%   1. samples1
%       1. samplesize
%       2. repitition1
%           1. initsamples
%           2. Otfilterdata
%               1. prior_samples
%               2. post_samples
%               3. run_time
%               4. MAP_estimate: if MAP_Flag is on
%           3. Daumfilterdata
%               1. prior_samples
%               2. post_samples
%               3. run_time
%               4. MAP_estimate: if MAP_Flag is on
%           4. EnKFfilterdata
%               1. prior_samples
%               2. post_samples
%               3. run_time
%               4. MAP_estimate: if MAP_Flag is on
%       3. repitition2
%       ....
%       M. repititionM
%   2. samples2
%   .....
%   N. realdata
%       1. x_true: true states at observation points
%       2. x_start: initial state
%       3. z_noise
%       4. fullorbit
%       5. delt
%       6. Nt
%   N+1: filterparams
%   N+2: simparams

%% First Plot
% Subplots showing evolution of the point cloud starting from the point
% cloud at T=0 which is same for all filters.
% After T=0 all the POSTERIOR POINT CLOUDS shown at different time points
% will have the true state, the mean estimate and the MAP point estimate
% shown in color.
% Data set is for multiple repitions and varying sample size.
% The plot will be shown for the input arguments [sample_size, rep_n]
% 
% Available sample sizes
xtruept = fulldata.realdata.x_true;
xstart = fulldata.realdata.x_start';
Nt = fulldata.realdata.Nt;
delt = fulldata.realdata.delt;
fullorbit = fulldata.realdata.fullorbit;
sample_sizes = fulldata.filterparams.sample_sizes;
nRep = fulldata.filterparams.nRep;
%%
% Display choise of sample sizes
disp("Choice of sample sizes:")
disp(sample_sizes);
sample_sizechoosen = 500; % (USER-DEFINED) -------------------------------------------
disp("Sample size chosen for plotting:")
disp(sample_sizechoosen);
% Display number of repititions available
disp("Number of repititions available:");
disp(nRep);
nRep_choosen = 1; % (USER-DEFINED)----------------------------------------------------
disp("Repititio number chosen for plotting:")
disp(nRep_choosen);
if(nRep_choosen>nRep);error("User defined repitition number exceeded");end
%%
% Extract data that correspond to sample size = sample_sizechoose
k = find(sample_sizes==sample_sizechoosen); % since sample_sizes is a vector
if(isempty(k));error("Choose sample size correctly");
end
% Data in fulldata.samplek needs to be extracted
dataplot = fulldata.(sprintf("samples%i",k));
pointcloud = dataplot.(sprintf("repitition%i",k));
% Verify sample size is correct
if(dataplot.samplesize ~= sample_sizechoosen);error("Sample sizes donot match");end
%
tpts = 4; % number of time point we show in the plot. tpoints <= Nt (USER-DEFINED)---
figure(1)
for i=1:tpts+1 % First one is T=0, initial samples 
ax1 = subplot(3,tpts+1,i);
%% OT
if(i==1)
    pointcloudplot(pointcloud.initsamples,xstart); % Just raw plot
    title('Initial Samples');
    pointcloudattributes('X','Y',[]);
else
   pointcloudplot(pointcloud.Otfilterdata.post_samples(:,:,i-1),xtruept(i-1,:));
   plotpoint(pointcloud.Otfilterdata.MAP_estimate(i-1,:),'*r')
   timepoint = delt*(i-1);
   texttitle = sprintf('Posterior at T = %0.2f',timepoint);
   title({'OT Filter',texttitle},'FontSize',10);
  pointcloudattributes('X','Y',[]);
end
ax2 = subplot(3,tpts+1,tpts+1+i);
%% EnKF
if(i==1)
    pointcloudplot(pointcloud.initsamples,xstart);% Just raw plot
    title('Initial Samples');
    pointcloudattributes('X','Y',[]);
else
    pointcloudplot(pointcloud.EnKFfilterdata.post_samples(:,:,i-1),xtruept(i-1,:)); 
    plotpoint(pointcloud.EnKFfilterdata.MAP_estimate(i-1,:),'*r')
    timepoint = delt*(i-1);
    texttitle = sprintf('Posterior at T = %0.2f',timepoint);
    title({'EnKF Filter',texttitle},'FontSize',10);
    pointcloudattributes('X','Y',[]);
end
ax3 = subplot(3,tpts+1,2*tpts+2+i);
%% Daum
if(i==1)
    pointcloudplot(pointcloud.initsamples,xstart); % Just raw plot
    title('Initial Samples');
    pointcloudattributes('X','Y',[]);
else
    pointcloudplot(pointcloud.Daumfilterdata.post_samples(:,:,i-1),xtruept(i-1,:)); 
    plotpoint(pointcloud.Daumfilterdata.MAP_estimate(i-1,:),'*r')
    timepoint = delt*(i-1);
    texttitle = sprintf('Posterior at T = %0.2f',timepoint);
    title({'Daum-Huang Filter',texttitle},'FontSize',10);
   if i == 3;pointcloudattributes('X','Y',{'samples','sample mean','true position','MAP Estimate'});
   else
      pointcloudattributes('X','Y',[]);
   end
end
linkaxes([ax1,ax2,ax3],'xy')
end

% filename = strcat(plotfolder,'/fullplot.pdf');
% fig = gcf;
% set(fig, 'PaperPositionMode', 'auto')
% if PDF_Flag ==1;export_fig(fig,'FontMode','fixed','FontSize',10,'Color',filename,'-nocrop');end

