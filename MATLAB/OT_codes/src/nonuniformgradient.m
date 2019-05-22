function [gradient] = nonuniformgradient(kpoints,priorfunc,evalpoint)
% nonuniformgradient calculated gradient at 'evalpoint' based on k
% nonuniform points 'kpoints'. The smooth function whose gradient is
% calculated is 'priorfunc'.
% [gradient] = nonuniformgradient(kpoints,priorfunc,evalpoint)
% INPUT: 
% 1. kpoints: neighboring k points with size k x d, where d is the size of
% each point
% 2. priorfunc: smooth function whose gradient is calculated
% 3. evalpoint: gradient is evaluated at this point. 'evalpoint' cannot be
% in the set 'kpoints'.
% Ref:
% Gradient Estimation from Irregularly Spaced Data Sets - Thomas H. Meyer,Marian Eriksson,and Robert C. Maggio

% First check that k>=d
[k,d1] = size(kpoints);
[~,d2] = size(evalpoint);
if k < d1 || d1~=d2
    error('Either k less than d or sizes of neighbor state and evaluation point donot match');
end
% Calculate all the connecting unit vectors from 'evalpoint' to each of
% 'kpoints'
tmp1 = kpoints - evalpoint;
tmp1norm = sqrt(sum(tmp1.*tmp1,2));
kvunitvec = tmp1./tmp1norm;
delS = priorfunc(kpoints)-priorfunc(evalpoint);
delSbyv = delS./tmp1norm;
% if rank(kvunitvec'*kvunitvec)<d1
%     error('Ill conditioned')
% end
gradient = pinv(kvunitvec)*delSbyv;
end