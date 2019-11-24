function params = getParams(net)
    w_tmp1 = net.IW{1,1};
    w_tmp2 = net.LW{2,1};
    w_tmp3 = net.LW{3,2};
    
    b_tmp1 = net.b{1,1};
    b_tmp2 = net.b{2,1};
    b_tmp3 = net.b{3,1};
    
    params = [w_tmp1(:); w_tmp2(:); w_tmp3(:); b_tmp1(:);  b_tmp2(:);  b_tmp3(:)]; 
    
    


end