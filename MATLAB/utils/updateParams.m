function net = updateParams(net, params)
    w_tmp1 = net.IW{1,1};
    net.IW{1,1} = reshape(params(1:numel(w_tmp1)), size(w_tmp1));
    
    w_tmp2 = net.LW{2,1};
    net.LW{2,1} = reshape( params(numel(w_tmp1)+1:numel(w_tmp1)+numel(w_tmp2)), size(w_tmp2));
    
    w_tmp3 = net.LW{3,2};
    net.LW{3,2} = reshape( params(numel(w_tmp1)+numel(w_tmp2)+1:numel(w_tmp1)+numel(w_tmp2)+numel(w_tmp3)), size(w_tmp3));
    
    biases = params(numel(w_tmp1)+numel(w_tmp2)+numel(w_tmp3)+1:end);
    
    b_tmp1 = net.b{1,1};
    net.b{1,1} =  reshape(biases(1:numel(b_tmp1)), size(b_tmp1));
    
    b_tmp2 = net.b{2,1};
    net.b{2,1} =  reshape(biases(numel(b_tmp1)+1:numel(b_tmp1)+numel(b_tmp2)  ), size(b_tmp2));
    
    b_tmp3 = net.b{3,1};
    net.b{3,1} =  reshape(biases(numel(b_tmp1)+numel(b_tmp2)+1:end ), size(b_tmp3));
       
end