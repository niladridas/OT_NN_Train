% Author: Niladri Das
% Date: 24th April 2019
clc;clear;close all;
format long 

load('data/X.mat');
X_nn = X(1:600,1); % use this X to train NN
% Normalize the data
X_nn = (X_nn - min(X_nn))./(max(X_nn)-min(X_nn));
del = 17; % delay parameter

for i=1:(size(X_nn,1)-del)
    iP(i,:) = X_nn(i:(i+del-1),1)';
    oP(i,1) = X_nn(i+del,1);
end


ni = del;
Ln = [10;12;1];
eta = 0.25;
maxitr = 2000;
NN = NNconstruct(ni,Ln); % Weights and Bias randomly initialized
avgDelta = 1;
tau =1;
% while avgDelta > 1e-3
for tau = 1: maxitr
    clc
    fprintf('Itr = %d.\n', tau)
    [Wnext,Bnext,avgDelta] = weightbiasup(NN,iP,oP,eta);
    NN.W = Wnext;
    NN.B = Bnext;
%     tau =tau+1;
end
% end

y1 = zeros(size(oP,1),size(oP,2));

for i = 1:size(iP,1)
    [A,Z] = forwardprop(NN,iP(i,:)');
    y1(i,1) = A{end};
end

%% Plot
% load('data/oP.mat'); % Normalized real output data 
% load('data/y1.mat'); % NN output
% load('data/X.mat');
figure; clf; hold on;
plot(oP,'r','Linewidth',1)
plot(y1,'b','Linewidth',1)
legend('Real','NN')%, 'X')






