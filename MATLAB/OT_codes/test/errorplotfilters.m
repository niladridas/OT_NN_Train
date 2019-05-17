% Author: Niladri Das
% Date: Aug 16, 2018
% About: Error plots between OT, EnKF and Particle Flow.
close all;clear;clc;
load('./data/fulldata');
xtruept = fulldata.realdata.x_true;
xstart = fulldata.realdata.x_start';
Nt = fulldata.realdata.Nt;
delt = fulldata.realdata.delt;
fullorbit = fulldata.realdata.fullorbit;
sample_sizes = fulldata.filterparams.sample_sizes;
nRep = fulldata.filterparams.nRep;

% Accross each nRep data calculate the mean of the posrterior estimated 
% Calculate mean of this mean 
% The calculate the distance from the true position and plot them with
% respect to time for OT EnKF and Particle flow
data = fulldata.samples1;
%
OTpriormean = zeros(Nt,4);
EnKFpriormean = zeros(Nt,4);
Daumpriormean = zeros(Nt,4);
OTruntime = zeros(1,Nt);
EnKFruntime = zeros(1,Nt);
Daumruntime = zeros(1,Nt); 

for i = 1:nRep
    OTprior = data.(sprintf('repitition%i',i)).Otfilterdata.post_samples;
    OTpriormean = OTpriormean + reshape(mean(OTprior,1),4,Nt)'/nRep;    
    OTruntime = OTruntime + data.(sprintf('repitition%i',i)).Otfilterdata.run_time/nRep;
    EnKFprior = data.(sprintf('repitition%i',i)).EnKFfilterdata.post_samples;
    EnKFpriormean = EnKFpriormean + reshape(mean(EnKFprior,1),4,Nt)'/nRep;     
    EnKFruntime = EnKFruntime + data.(sprintf('repitition%i',i)).EnKFfilterdata.run_time/nRep;
    Daumprior = data.(sprintf('repitition%i',i)).Daumfilterdata.post_samples;
    Daumpriormean = Daumpriormean + reshape(mean(Daumprior,1),4,Nt)'/nRep;      
    Daumruntime = Daumruntime + data.(sprintf('repitition%i',i)).Daumfilterdata.run_time/nRep;
end

OTpriormeanxy = [OTpriormean(:,1).*cos(OTpriormean(:,3)),OTpriormean(:,1).*sin(OTpriormean(:,3))];
EnKFpriormeanxy = [EnKFpriormean(:,1).*cos(EnKFpriormean(:,3)),EnKFpriormean(:,1).*sin(EnKFpriormean(:,3))];
Daumpriormeanxy = [Daumpriormean(:,1).*cos(Daumpriormean(:,3)),Daumpriormean(:,1).*sin(Daumpriormean(:,3))];
truexy = [xtruept(:,1).*cos(xtruept(:,3)),xtruept(:,1).*sin(xtruept(:,3))];
disterrorOT = sqrt(sum((OTpriormeanxy - truexy).*(OTpriormeanxy - truexy),2));
disterrorEnKF = sqrt(sum((EnKFpriormeanxy - truexy).*(EnKFpriormeanxy - truexy),2));
disterrorDaum = sqrt(sum((Daumpriormeanxy - truexy).*(Daumpriormeanxy - truexy),2));

scatter(1:Nt,disterrorOT,'filled');hold on;
scatter(1:Nt,disterrorEnKF,'filled');hold on;
scatter(1:Nt,disterrorDaum,'filled');hold on;
legend('OT','EnKF','Daum');
plot(1:Nt,disterrorOT);hold on;
plot(1:Nt,disterrorEnKF);hold on;
plot(1:Nt,disterrorDaum);hold off;

