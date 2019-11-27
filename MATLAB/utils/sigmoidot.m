function y = sigmoidot(x)
%     y = (1-sigmoid(x)).*sigmoid(x);
    
%     y = (x>0);  ReLU

      y = 1-tanh(x).^2;
end