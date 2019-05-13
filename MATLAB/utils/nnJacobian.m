% Author: Vedang Deshpande
% Date: 13th May 2019
% Inputs: 
% NN -  a neural network
% param - vector of parameters (weights and biases)
% ip - inputs to the NN, column vector
% Outputs: 
% H - Jacobian of NN wrt parameters, evaluated at current parameters of NN
function H = nnJacobian(NN,ip)
    [A,Z] = forwardprop(NN,ip);
    for i = 1:length(Z)
        der{i} = sigmoidot(Z{i});
    end
    for i = length(NN.W):-1:1
        if i == length(NN.W)
            Delta{i} =  der{i};
        else
            Delta{i} = ((NN.W{i+1})'*Delta{i+1}).*der{i};
        end
    end
    dOPdW = CderW(NN.W,Delta,A,ip); % derivative of o/p wrt weights
    dOPdB = Delta; % derivative of o/p wrt biases (Bias correspond to input=1)
    
    % Derivatives are in matrix form, convert to a row vector
    tempNN = NN;
    tempNN.W = dOPdW;
    tempNN.B = dOPdB;
    H = nn2param(tempNN)';
end