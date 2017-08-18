function [ b ] = magfield( r )
%MAGFIELD Calculate magnetic field vector.
%   MAGFIELD(R) calculates the magnetic field vector at the point R in the
%   inertial frame.

B0 = 3.12e-5;
rE = 6378e3;

a = norm(r);

el = acos(r(3) / a);
az = atan2(r(2), r(1));

R = [sin(el)*cos(az)    cos(el)*cos(az) -sin(az); ...
     sin(el)*sin(az)    cos(el)*sin(az) cos(az); ...
     cos(el)            -sin(el)        0];

br = -2 * B0 * (rE/a)^3 * cos(el);
bel = -B0 * (rE/a)^3 * sin(el);
      
b = reshape(R * [br bel 0]', size(r));

end

