function [ s ] = mekf( s )
%MEKF Multiplicative Extended Kalman Filter.
%   Detailed explanation goes here

if isnan(s.x)
    s.x = s.H \ s.z;
    s.P = s.H \ s.R / s.H';
else
    % Predict step
    s.x = s.A*s.x + s.B*s.u;
    s.P = s.A * s.P * s.A' + s.Q;

    % Compute Kalman gain
    K = s.P * s.H' / (s.H * s.P * s.H' + s.R);

    % Update step
    s.x = s.x + K*(s.z - s.H*s.x);
    s.P = s.P - K*s.H*s.P;
end

end

