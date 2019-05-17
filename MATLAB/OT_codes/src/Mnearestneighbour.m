function [Mneighboursindex] = Mnearestneighbour(Nparticles,particlenumber,Mvalue)
% Mnearestneighbour Calculates the M nearest neighbour to particle number
% 'particlenumber' out of N particles, with respect to the projection value
% 'Nparticles.projvalues'. The N particles are arranged in
% 'Nparticles.positions' according to the descending value of 'Nparticles.projvalues'
% [Mneighbours] = Mnearestneighbour(Nparticles,particlenumber,Mvalue)
% 

% First change origin of 'Nparticles.projvalues' with respect to the particlenumber.
Nparticles.projvalues = Nparticles.projvalues - Nparticles.projvalues(particlenumber);
% The 'particlenumber' entry in 'Nparticles.projvalues' is zero.
% Find out the closest M neighbours based on the new 'Nparticles.projvalues'
% First take abs() of 'Nparticles.projvalues'
Nparticles.projvalues = abs(Nparticles.projvalues);
count = 0;
N = length(Nparticles.projvalues);
% Removing the cases of 'particlenumber' equal to 1 or N.
if particlenumber == 1
    Mneighboursindex = 2:(Mvalue+1);
elseif particlenumber == N
    Mneighboursindex = (N-Mvalue-1):(N-1);
else
    [~,I] = sort(Nparticles.projvalues);
    Mneighboursindex = I(2:(Mvalue+1));
end 
end