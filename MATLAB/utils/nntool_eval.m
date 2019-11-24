function val =  nntool_eval(nn,w,ip1)
    nn = updateParams(nn,w);
    val = nn(ip1);
end