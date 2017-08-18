function [ M ] = gravgrad( J, r, n )
%GRAVGRAD Gravity gradient in body frame, assuming circular orbit.
%   GRAVGRAD(J,R) computes the torque M acting on a satellite body
%   with inertia matrix J at a position R with attitude represented by the
%   quaternion Q. The position R must be in the body frame.

if nargin < 3
    mu = 3.986005e14;
    n = sqrt(mu / norm(r)^3);
end

rHat = r / norm(r);
M = 3 * n^2 * skew(rHat) * J * rHat;

end

