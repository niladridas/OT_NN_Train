function dpropagatefcn_dx = NonStationaryGrowth_dpropagatefcn_dx(xp,prop_params )
dpropagatefcn_dx  = 0.5+25*(1-xp.^2)./((1+xp.^2).^2);
end