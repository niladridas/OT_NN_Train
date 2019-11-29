function pdfval = GaussBimodalPdf(A,B,W,x)
% GaussBimodalPdf  Gives PDF value from Gaussian Bimodal Distrbution.
%   pdfval = GAUSSBIMODAL(A) generates 1 sample from the mean location in A.
%
%   pdfval = GAUSSBIMODAL(A,n) generates n sample from the mean location in A.
%
%   pdfval = GAUSSBIMODAL(A,n,B) generates n sample from the mean location in A and SD in B.
%
%   pdfval = GAUSSBIMODAL(A,n,B,W) generates n sample from the mean location in A
%   and SD in B, with relative weighting in W
%
%   See also GMDISTRIBUTION, RANDOM.
switch nargin
    case 1
    % Only mean locations are given
        m = [A(1); A(2)];
        s = cat(3,1,1);
        w = ones(1,2) / 2;
        gmd = gmdistribution(m,s,w);
        pdfval = gmd.pdf(x);
    case 2
    % Only mean locations and sample number given
        m = [A(1); A(2)];
        s = cat(3,1,1);
        w = ones(1,2) / 2;
        gmd = gmdistribution(m,s,w);
        pdfval = gmd.pdf(x);
    case 3
    % Only mean locations, sample number, and sigma given
        m = [A(1); A(2)];
        s = cat(3,B(1),B(2));
        w = ones(1,2) / 2;
        gmd = gmdistribution(m,s,w);
        pdfval = gmd.pdf(x);
    case 4
    % Only mean locations, sample number, and sigma given
        m = [A(1); A(2)];
        s = cat(3,B(1),B(2));
        w = [W(1),W(2)];
        gmd = gmdistribution(m,s,w);
        pdfval = gmd.pdf(x);
    otherwise
    % Only mean locations are given
        m = [-1; 1];
        s = cat(3,1,1);
        w = ones(1,2) / 2;
        gmd = gmdistribution(m,s,w);
        pdfval = gmd.pdf(x);
end
end