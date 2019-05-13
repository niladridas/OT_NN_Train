% Author: Vedang Deshpande
% Date: 11th May 2019
% Dynamics for the weights/parameters of neural network
% Inputs: 
% NN -  a neural network
% param - vector of parameters
% Outputs: 
% NN - neural network with updated weights and bias as per params
function NN = param2nn(NN,param)
    numParam = 0;
    for i=1:length(NN.W)
        rw = size(NN.W{i},1);
        col = size(NN.W{i},2);
        NN.W{i} =  reshape(param(numParam+1:numParam+rw*col),rw,col);
        numParam = numParam + rw*col;
    end
    for i=1:length(NN.B)
        rw = size(NN.B{i},1);
        col = size(NN.B{i},2);
        NN.B{i} =  reshape(param(numParam+1:numParam+rw*col),rw,col);
        numParam = numParam + rw*col;
    end
end