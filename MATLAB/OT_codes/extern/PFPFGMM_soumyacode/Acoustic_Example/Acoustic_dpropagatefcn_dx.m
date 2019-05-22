function dpropagatefcn_dx = Acoustic_dpropagatefcn_dx(xp,prop_params )
dpropagatefcn_dx  = repmat(prop_params.Phi,1,1,size(xp,2));
end