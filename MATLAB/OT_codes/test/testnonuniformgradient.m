% Author: Niladri Das
% About: Testing correctness of the function 'nonuniformgradient.m'
clc;clear;close all;
% Problem Setup:
% The testing surface is a hyperboloid given by the equation:
% z = 115*((e-364560)/420)^2 - 174((n-3901140)/450)^2+2400
% e: assume it to be x co-ordinates
% n: assume it to be y co-ordinates
% e = [364140,364980]; n = [3900690,3901590]

% % We plot the surface
% x = 364140:364980;
% y = 3900690:3901590;
% [X,Y] = meshgrid(x,y);
% Z = 115*((X-364560)./420).^2 - 174*((Y-3901140)/450).^2+2400;
% surf(X,Y,Z);


% First we generate N number of random points on this space of (e,n)
rng('default');
rng(1);
N = 1000;
Npoints = [(364980-364140);(3901590-3900690)].*rand(2,N)+[364140;3900690];
Npoints = Npoints';
% Calculate gradient at each of these N points based on the k nearest
% neighbour
k = 50;
%%
gradientN = zeros(N,2);
gradientNreal = zeros(N,2);
for i = 1:N
    x = Npoints(i,:);
    normdelx = sqrt(sum((Npoints - x).*(Npoints - x),2));
    [~,I] = sort(normdelx);
    I = I(2:(k+1));
    xkneighboridx = I; 
    xkneighbor = Npoints(I,:);
    xkseval = zeros(k,1);
    % Calculate the gradient
    % Evaluate function values at these points
    for j = 1:k
        xkseval(j) = zcal(xkneighbor(j,:));
    end
    [gradient] = nonuniformgradient(xkneighbor,@zcal,x);
    gradientN(i,:) = gradient'; 
    gradientNreal(i,:) = gradrealcal(x);
end

%% Error in calculating gradient
errorangle = zeros(N,1);
dottmp = sum(gradientN.*gradientNreal,2);
absgradientN = sqrt(sum(gradientN.*gradientN,2));
absgradientNreal = sqrt(sum(gradientNreal.*gradientNreal,2));
costhetatmp = dottmp./(absgradientN.*absgradientNreal);
thetaerror = acosd(costhetatmp);
figure(1);plot(thetaerror);
magerror = (absgradientN - absgradientNreal);
figure(2);plot(magerror);

%%
function [val] = zcal(X)
    val = 115*((X(:,1)-364560)/420).^2 - 174*((X(:,2)-3901140)/450).^2+2400;
end

function [gradreal] = gradrealcal(X)
    gradrealx =  115*2*(X(:,1)-364560)/(420^2);
    gradrealy =  -174*2*(X(:,2)-3901140)/(450^2);
    gradreal = [gradrealx,gradrealy];
end

