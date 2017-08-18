function [ u ] = rcs( dtheta, dthetadot, opts )
%RCS Summary of this function goes here
%   Detailed explanation goes here

c1 = opts.UpperLim;
c2 = opts.LowerLim;
a3 = -opts.Deadband/2;
a6 = opts.Deadband/2;
a2 = a3 - 0.1;
a7 = a6 + 0.1;

thrust = opts.Thrust;

u = [0 0 0];
for i = 1:3
    k1 = 1/(2*thrust(i));
    a1 = -k1 * c2^2 + a2;
    a8 = k1 * c2^2 + a7;
    theta = rad2deg(dtheta(i));
    thetadot = rad2deg(dthetadot(i));
    if thetadot >= c1
        u(i) = -thrust(i);
    elseif theta >= (-k1*thetadot^2 + a6) && thetadot > 0 && thetadot < c1
        u(i) = -thrust(i);
    elseif thetadot >= 0 && theta >= a6
        u(i) = -thrust(i);
    elseif theta >= (k1*thetadot^2 + a7) && thetadot > -c2 && thetadot < 0
        u(i) = -thrust(i);
    elseif theta >= a8 && thetadot > -c2
        u(i) = -thrust(i);
    elseif thetadot <= -c1
        u(i) = thrust(i);
    elseif thetadot <= c2 && theta <= a1
        u(i) = thrust(i);
    elseif theta <= (-k1*thetadot^2 + a2) && thetadot < c2 && thetadot > 0
        u(i) = thrust(i);
    elseif thetadot <= 0 && theta <= a3
        u(i) = thrust(i);
    elseif theta <= (k1*thetadot^2 + a3) && thetadot < 0 && thetadot > -c1
        u(i) = thrust(i);
    elseif thetadot <= -c1
        u(i) = thrust(i);
    end
end

end

