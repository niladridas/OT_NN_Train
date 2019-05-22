function [knearestindex] = knearestneighbour(Mnearestindex,Nparticles,particlenumber,kvalue)
% knearestneighbour calculates the index of the k nearest neighbour in terms 
% of a distance out of the M nearest neighbour calculate in terms of
% projection distance.
% [knearestindex] = knearestneighbour(Mnearestindex,particle,kvalue)
if length(Mnearestindex) < kvalue
    error('k must be less or equal to M')
end
Mnearest = Nparticles.positions(Mnearestindex,:);
% Calculate distances
tmp = (Mnearest-Nparticles.positions(particlenumber,:));
tmp = sqrt(sum(tmp.*tmp,2));
[~,I] = sort(tmp);
I = I(1:kvalue);
knearestindex = Mnearestindex(I);
end