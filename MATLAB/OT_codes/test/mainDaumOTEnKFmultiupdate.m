% Author: Niladri Das
% Date: Aug 9, 2018
% About: Daum-Huang filter is compared with OT filter for several measurements
clear; clc;close all;
init;
simulation_params;
filter_params;
[realdata] = gen_realdata(simparams);
sample_sizes = filterparams.sample_sizes;
nRep = filterparams.nRep; % repitition of the experiment with 'nRep' different sets of samples each of size 'nSamp'
run_time = zeros(2,length(sample_sizes),nRep);
est_error = zeros(2,length(sample_sizes),nRep);
rng('default');
seeds = 1:1:nRep; % Generate nRep seeds

datasavefolderhd = dir('data');
datafolder = datasavefolderhd.folder;
%%
for k = 1:length(sample_sizes)
    nSamp = sample_sizes(k);
    fulldata.(sprintf('samples%i', k)).samplesize = nSamp;
    for m = 1:nRep
        rng(seeds(m));
        [initsamples] = gen_initsamples(filterparams,simparams,nSamp);
        [Otfilterdata] = OTfilterandprop(initsamples,realdata.z_noise,simparams,filterparams);
%         [Otfilterdata] = OTsinkfilterandprop(initsamples,realdata.z_noise,simparams,filterparams);
        % ode23 is used in Daum. Explore more.
        % This is the linear formulation of Daum
        [Daumfilterdata] = Daumfilterandprop(initsamples,realdata.z_noise,simparams,filterparams);
        [EnKFfilterdata] = EnKFfilterandprop(initsamples,realdata.z_noise,simparams,filterparams);
        % save data for multiple repititions
        fulldata.(sprintf('samples%i', k)).(sprintf('repitition%i', m)).initsamples = initsamples;
        fulldata.(sprintf('samples%i', k)).(sprintf('repitition%i', m)).Otfilterdata = Otfilterdata;
        fulldata.(sprintf('samples%i', k)).(sprintf('repitition%i', m)).Daumfilterdata = Daumfilterdata;
        fulldata.(sprintf('samples%i', k)).(sprintf('repitition%i', m)).EnKFfilterdata = EnKFfilterdata;
    end
end
fulldata.realdata = realdata;
fulldata.filterparams = filterparams;
fulldata.simparams = simparams;
filename = strcat(datafolder,'/fulldata.mat');
save(filename,'fulldata');

