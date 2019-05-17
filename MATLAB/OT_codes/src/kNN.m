function [Nsortparticles] = kNN(Nparticles,Mvalue,kvalue,PriorMAT)
    [maxeigenvector]= calmaxeigenvector(PriorMAT);
    [Nsortparticles] = Nsortprojection(Nparticles,maxeigenvector);
    % 'Nparticles.projvalues'. The N particles are arranged in
    % 'Nparticles.positions' according to the descending value of 'Nparticles.projvalues'
    N = length(Nsortparticles.projvalues);
    kindexmat = zeros(kvalue,N);
    for i = 1:N
        [Mneighboursindex] = Mnearestneighbour(Nsortparticles,i,Mvalue);
        [knearestindex] = knearestneighbour(Mneighboursindex,Nsortparticles,i,kvalue);
        kindexmat(:,i) = knearestindex;
    end
    Nsortparticles.kindex = kindexmat;
end