function [Nsortparticles] = Nsortprojection(Nparticles,projvec)
% Nsortprojection  Sorts the input N particles with respect to their
% projection on the projvec vector in descending order
% [Nsortparticles] = Nsortprojection(Nparticles,projvec)
% Nparticles.positions : are the N number of input particles with dimension N x d
% Nparticles.projvalues : are the N number of projection values
% projvec: the vector on which the projections of the N particles will be
% calculated. It has a dimension of d x 1.

% First we calculate the projections
projvalues = Nparticles*projvec;
[M,I] = sort(projvalues,'descend');
Nsortparticles.positions = Nparticles(I,:);
Nsortparticles.projvalues = M(I);
end