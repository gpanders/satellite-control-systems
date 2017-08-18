function [ wHat, qHat, P ] = nav( wTilde, qTilde, R, qPrev, wPrev, betaPrev, PPrev, opts )
%NAV Estimate angular velocity and attitude.
%   [WHAT,QHAT,P] = NAV(WTILDE,QTILDE,R,QPREV,WPREV,BETAPREV,PPREV,OPTS) is a
%   Multiplicative Extended Kalman Filter (MEKF) that estimates the angular
%   velocity WHAT, attitude quaternion QHAT, and state covariance matrix P given
%   an angular velocity measurment WTILDE, measurement quaternion QTILDE and
%   measurement covariance R, the previous quaternion estimate PREV, the 
%   previous gyro bias estimate BETAPREV, and the previous error covariance 
%   matrix PPREV. OPTS contains options that define simulation and estimator
%   parameters.
%
%   See also SENSORS, PHIF, CONTROLLER.

% Integration time step
dt = 1/opts.ComputerFreq;

% Propagate quaternion reference from last time step to current time step
qHat = qpropagate(qPrev, wPrev, dt);

dq = qmult(qTilde, qconj(qHat));
da = 2*dq(1:3) / dq(4); % Use Gibbs vector parameterization

% Define sensitivity matrix
H = eye(3, 6);

% Define process noise
Q = blkdiag(opts.GyroNoise^2 * eye(3), opts.GyroBiasRateNoise^2 * eye(3));

% Predict P
G = blkdiag(-eye(3), eye(3));
F = phiF(wPrev, dt);
P = F*PPrev*F' + G*(Q*dt)*G';

% Calculate Kalman gain
S = H*P*H' + R;
K = P*H' / S;

dx = K*da;
da = dx(1:3);
dbeta = dx(4:6);
dq = 1/sqrt(4 + norm(da)^2) * [da; 2];
beta = betaPrev + dbeta;

wHat = wTilde - beta;
qHat = qnormalize(qmult(dq, qHat));
P = (eye(6) - K*H)*P;

end