function [Delta] = backprop(E,Z,NN)
    for i = length(NN.W):-1:1
        if i == length(NN.W)
            Delta{i} = -E.*sigmoidot(Z{i});% delta_L
        else
            Delta{i} = ((NN.W{i+1})'*Delta{i+1}).*sigmoidot(Z{i});
        end
    end
end
