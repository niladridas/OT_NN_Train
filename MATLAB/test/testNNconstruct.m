% Author: Niladri Das
% Date: 24th April 2019
clc;clear;close all;
format long 

ni = 1;
Ln = [20;1];
eta = 0.55;
maxitr = 1000;
NN = NNconstruct(ni,Ln); % Weights and Bias randomly initialized

% Create data set
x = linspace(0,4*pi,100)';
y = 0.5*(sin(x)+1);
plot(x,y);
for tau = 1: maxitr
    Wnext = weightup(NN,x,y,eta);
    NN.W = Wnext;
end
y1 = zeros(size(y,1),size(y,2));
for i = 1:size(x,1)
    [A,Z] = forwardprop(NN,x(i,:)');
    y1(i,1) = A{end};
end
hold on; plot(x,y1);






