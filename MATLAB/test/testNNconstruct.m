% Author: Niladri Das
% Date: 24th April 2019
clc;clear;close all;
format long 

load('data/X.mat');
% Construct data set from X
% Input matrix (no. of samples x dimension of states)
del = 17;
iP = zeros(size(X,1), del);
oP = X;

for i = 1:size(X,1)
    for j = 1:del
        if (i-(del-j+1))<=0
            iP(i,j) = 0;
        else 
            iP(i,j) = X(i-(del-j+1),1);
        end
    end
end

ni = 17;
Ln = [10;10;1];
eta = 0.55;
maxitr = 1000;
NN = NNconstruct(ni,Ln); % Weights and Bias randomly initialized

% % Create data set
% x = linspace(0,4*pi,100)';
% y = 0.5*(sin(x)+1);
% plot(x,y);
for tau = 1: maxitr
    [Wnext,Bnext] = weightbiasup(NN,iP,oP,eta);
    NN.W = Wnext;
    NN.B = Bnexy;
end
y1 = zeros(size(y,1),size(y,2));
for i = 1:size(x,1)
    [A,Z] = forwardprop(NN,x(i,:)');
    y1(i,1) = A{end};
end
hold on; plot(x,y1);






