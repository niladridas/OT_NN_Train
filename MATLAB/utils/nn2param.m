% Author: Vedang Deshpande
% Date: 11th May 2019
% Inputs: 
% NN -  a neural network
% Outputs: 
% param - Column of parameters (weights and biases)
function param = nn2param(NN)
    vecB=[];vecW=[];
    for i=1:length(NN.W)
        vecW = [vecW;vec(NN.W{i})];
    end
    for i=1:length(NN.B)
        vecB = [vecB;vec(NN.B{i})];
    end
    param = [vecW;vecB];
end