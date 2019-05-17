function [p] = likfun(y,x,meas_sig,H)
    p = exp(-0.5*(y-H*x)'*(meas_sig\(y-H*x)));
end