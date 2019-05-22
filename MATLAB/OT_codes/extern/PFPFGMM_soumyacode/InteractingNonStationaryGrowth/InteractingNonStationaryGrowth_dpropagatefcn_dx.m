function dpropagatefcn_dx = InteractingNonStationaryGrowth_dpropagatefcn_dx(xp,prop_params )
[dim,nParticle] = size(xp);
dpropagatefcn_dx = zeros(dim,dim,nParticle);
a = 0.5;
b = 2.5;
c = 8;

dpropagatefcn_dx(1,1,:) = a -2*b*xp(1,:).*xp(2,:)./((1+xp(1,:).^2).^2);
dpropagatefcn_dx(1,2,:) = b./(1+xp(1,:).^2);
for i = 2:dim-1
    dpropagatefcn_dx(i,i-1,:)  = -2*b*xp(i+1,:).*xp(i-1,:)./((1+xp(i-1,:).^2).^2);
    dpropagatefcn_dx(i,i,:) = repmat(a,1,nParticle);
    dpropagatefcn_dx(i,i+1,:) = repmat(b,1,nParticle)./(1+xp(i-1,:).^2);
end
dpropagatefcn_dx(dim,dim-1,:) = -2*b*xp(dim,:).*xp(dim-1,:)./((1+xp(dim-1,:).^2).^2);
dpropagatefcn_dx(dim,dim,:) = a + repmat(b,1,nParticle)./(1+xp(dim-1,:).^2);
end