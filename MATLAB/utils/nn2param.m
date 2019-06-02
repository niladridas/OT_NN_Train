% Author: Vedang Deshpande
% Date: 11th May 2019
% Inputs: 
% NN -  a neural network
% Outputs: 
% param - Column of parameters (weights and biases)
function param = nn2param(NN)
    vecB=[];vecW=[];
    for i=1:length(NN.W)
%        vecW = [vecW;vec(NN.W{i})];
        % Matlab 2016
        tmp = NN.W{i};
        tmp = tmp';
        vecW = [vecW;tmp(:)];        
    end
    for i=1:length(NN.B)
        tmp1 = NN.B{i};
        tmp1 = tmp1';
        vecB = [vecB;tmp1(:)];
    end
    param = [vecW;vecB];
end