% Author: Niladri Das
% Date: 24th April 2019
clc;clear;close all;
format long 

load('data/X.mat');
X = X(1:400,1);
% Normalize the data
X = (X - min(X))./(max(X)-min(X));
del = 17; % delay parameter

for i=1:(size(X,1)-del)
    iP(i,:) = X(i:(i+del-1),1)';
    oP(i,1) = X(i+del,1);
end


ni = 17;
Ln = [10;10;1];
eta = 0.35;
maxitr = 2000;
NN = NNconstruct(ni,Ln); % Weights and Bias randomly initialized

for tau = 1: maxitr
    [Wnext,Bnext] = weightbiasup(NN,iP,oP,eta);
    NN.W = Wnext;
    NN.B = Bnext;
end


y1 = zeros(size(oP,1),size(oP,2));

for i = 1:size(iP,1)
    [A,Z] = forwardprop(NN,iP(i,:)');
    y1(i,1) = A{end};
end








