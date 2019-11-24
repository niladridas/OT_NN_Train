function J = nntool_Jacobian(net,ip)
    w0 = getParams(net);
    fun = @(w) nntool_eval(net,w,ip);
    del = 1e-4;
    delta_f = zeros(1,numel(w0));
    parfor i = 1:numel(w0)
        w_f = w0; w_f(i) = w0(i) + del;
        w_b = w0; w_b(i) = w0(i) - del;
        delta_f(i) = fun(w_f) - fun(w_b);
    end
    
    J = delta_f/(2*del);
%     J =  jacobianest(fun,w0);
    

end

