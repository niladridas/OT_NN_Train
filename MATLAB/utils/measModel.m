% Author: Vedang Deshpande
% Date: 11th May 2019
% Measurement model for neural network
% Inputs: 
% NN - neural network under training
% ip - inputs to the NN, column vector
% 
% Outputs: 
% y - prediction of output
function y = measModel(NN,ip)
        [A,~] = forwardprop(NN,ip);
        y = A{end};
end