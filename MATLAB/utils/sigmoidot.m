function y = sigmoidot(x)
    y = (1-sigmoid(x)).*sigmoid(x);
end