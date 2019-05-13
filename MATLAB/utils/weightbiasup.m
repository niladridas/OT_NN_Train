function [Wnext,Bnext,tempDelta] = weightbiasup(NN,x,y,eta,CdW)
    for i = 1:length(NN.W)
        CdW{i} = zeros(size(NN.W{i},1),size(NN.W{i},2));
    end
    tempDelta = zeros(size(NN.W{end},1));
    for i = 1:size(x,1) % Number of individual samples
        % Forward propagation
        [A,Z] = forwardprop(NN,x(i,:)');

        % Error signal in the final layer
        E = y(i,:)-A{end};

        % Backward propagation
        Delta = backprop(E,Z,NN);
        tempDelta = tempDelta + Delta{end};

        % Derivative of the cost with respect to weights
        tmpCdW = CderW(NN.W,Delta,A,x(i,:)') ;

        for t = 1:length(NN.W)
            CdW{t} = CdW{t}+tmpCdW{t}/(size(x,1));
        end  
    end
    tempDelta = tempDelta/size(x,1);
    % Weight update
    for i = 1:length(NN.W)
        Wnext{i} = NN.W{i} - eta*CdW{i};
    end
    % Bias update
    for i = 1:length(NN.B)
        Bnext{i} = NN.B{i} - eta*Delta{i};
    end
end