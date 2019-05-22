function [A,H] = lin_param_sat(x,J2Flag)
    Tp = 102.8*60; % Orbit Period in sec
    Re = 6378.1363; %km radius of earth

    mu = 398600.4415; %(km^3/sec^2)
    mub = mu*Tp^2/Re^3; % Normalized

    J2 = 1.7555*10^10;
    J2b = J2*Tp^2/Re^5;

    % J2 effect
    if J2Flag
        
        F = [0,0,0,0;...
            -4*J2b*(3/2)*(3*sin(x(3))^2-1)/x(1)^5,0,J2b*(3/2)*(6*cos(x(3)))/x(1)^4,0;...
            0,0,0,0;...
            4*J2b*3*cos(x(3))*sin(x(3))/x(1)^5,0,-J2b*3*cos(2*x(3))/x(1)^4,0];
    else
        F = zeros(4,4);
    end
    
    A = [0,1,0,0;...
        2*mub/x(1)^3 + x(4)^2, 0, 0, 2*x(1)*x(4);...
        0,0,0,1;...
        2*x(2)*x(4)/x(1)^2,-2*x(4)/x(1),0,-2*x(2)/x(1)];
    A = A + F;
    H = [1,0,0,0;0,0,1,0];
end