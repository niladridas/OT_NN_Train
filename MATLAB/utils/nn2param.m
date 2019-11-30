% Author: Vedang Deshpande
% Date: 11th May 2019
% Inputs: 
% NN -  a neural network
% Outputs: 
% param - Column of parameters (weights and biases)
function param = nn2param(NN)
    vecB=[];vecW=[];
    for i=1:length(NN.W)
        t1 = NN.W{i};
        vecW = [vecW;t1(:)];
    end
    for i=1:length(NN.B)
        t2 = NN.B{i};
        vecB = [vecB;t2(:)];
    end
    param = [vecW;vecB];
end