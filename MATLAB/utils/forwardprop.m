function [A,Z] = forwardprop(NN,x)
    for i = 1:length(NN.W) % Propagating through each layer
        if i == 1
            a = x;
        else
            a = A{i-1};
        end
        Z{i} = NN.W{i}*a+NN.B{i};
        A{i} = sigmoid(Z{i});
    end
end