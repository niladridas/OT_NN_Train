function dpropagatefcn_dx = Septier16_dpropagatefcn_dx(xp,prop_params )
dpropagatefcn_dx  = repmat(prop_params.transitionMatrix,1,1,size(xp,2));
end

