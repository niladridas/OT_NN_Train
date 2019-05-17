function xdot = satelliteDynamics(~,x,J2Flag,Tp,Re)
    mu = 398600.4415; %(km^3/sec^2)
    mub = mu*Tp^2/Re^3; % Normalized
    J2 = 1.7555*10^10;
    J2b = J2*Tp^2/Re^5;
    % J2 effect
    if J2Flag
        r = x(1);
        th = x(3);
        Fr = J2b*(3/2)*(3*sin(th)^2-1)/r^4;
        Fth = -J2b*3*cos(th)*sin(th)/r^4;
    else
        Fr = 0;
        Fth = 0;
    end
    xdot(1,1) = x(2);
    xdot(2,1) = -mub/x(1)^2 + x(1)*x(4)^2 + Fr;
    xdot(3,1) = x(4);
    xdot(4,1) = -2*x(2)*x(4)/x(1) + Fth;
end