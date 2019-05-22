% Author: Niladri Das
% Test calculating k nearest neighbour of points in a set of N particles
clc;clear;close all;
rng('default'); rng(1);
d = 6; % dimension of the state
N = 20; % number of samples
Mvalue = 10;
kvalue = 3;
PriorMAT = 0.2*eye(6);
Nparticles = mvnrnd(zeros(6,1),0.2*eye(6),N);
% [maxeigenvector]= calmaxeigenvector(PriorMAT);
% [Nsortparticles] = Nsortprojection(Nparticles,maxeigenvector);
% for i = 1:N
%     [Mneighboursindex] = Mnearestneighbour(Nsortparticles,i,Mvalue);
%     [knearestindex] = knearestneighbour(Mneighboursindex,Nsortparticles,i,kvalue);
% end


[Nsortparticles] = kNN(Nparticles,Mvalue,kvalue,PriorMAT);