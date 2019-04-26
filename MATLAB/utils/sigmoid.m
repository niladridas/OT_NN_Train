function y = sigmoid(x)
    y = 1+exp(-x);
    y = 1./y;
end