function [ wTilde, qTilde, R ] = sensors( w, q, r, bias, opts )
%SENSORS Simulate sensors.
%   [WTILDE,QTILDE,R] = SENSORS(W,Q,R,BIAS,OPTS) simulates a gyro sensor
%   measurement with bias BIAS of the angular velocity WTILDE and a quaternion 
%   attitude measurement QTILDE with measurement covariance matrix R given the
%   true (simulated) value of the angular velocity W, the true (simulated) 
%   attitude quaternion Q, a 3-by-N weighted set of reference measurement 
%   vectors R.  
%
%   OPTS must be a struct containing parameters that define the sensors.
%   The following parameters are required:
%       - 'MeasurementNoise' : a N-by-1 vector representing the standard
%                                 deviation of each corresponding
%                                 measurement
%       - 'GyroNoise'        : standard deviation of gyro noise

dt = 1/opts.ComputerFreq;

T = quat2rotm(q)';
b = zeros(size(r));

% Normalize measurements
r = normc(r);

% Calculate noisy measurements in body frame
sigma = opts.MeasurementNoise;
heta = bsxfun(@times, randn([3 length(sigma)]), sigma);
for ii = 1:size(r, 2)
    b(:, ii) = T * r(:, ii) + heta(:, ii);
end

% Normalize after adding noise
b = normc(b);

% Calculate attitude measurement
weights = 1./(opts.MeasurementNoise).^2;
[T, B] = svdatt(b, r, weights);

wTilde = w + bias + (opts.GyroNoise / sqrt(dt)) * randn([3 1]);
qTilde = rotm2quat(T');
R = inv(-B*T' + trace(B*T')*eye(3));

end

