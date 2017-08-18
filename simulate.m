%% Simulation

%% Setup
%
% Define constants and utility functions
%
J2 = 1.082e-3;  % J2 constant
rE = 6378e3;    % radius of earth
GM = 3.986004415e14;    % gravitational constant
a = rE + 400e3; % semi-major axis of satellite orbit
I = pi/6;   % orbit inclination
n = sqrt(GM / a^3) * (1 + (3/2) * (rE/a)^2 * J2 * (1 - 3*cos(I)^2));   % mean motion of sat
orbitPeriod = 2*pi / n; % satellite orbit period
J = diag([90 70 60]);  % satellite inertia matrix
fSim = 100;  % simulation sampling frequency
fCom = 10;  % flight computer sampling frequency
magmoment = skew([1 0 0]); % satellite magnetic moment
numTrials = 1;

dt = 1/fSim;
t0 = 0;                                 % initial time
tf = 0.5 * orbitPeriod;                   % final time
tVec = 0:dt:tf-dt;                      % time vector
Nt = length(tVec);

% Sun reference measurement (assumed constant)
rSun = [1 0 0]';

% Position as a function of time
x = @(t) a*cos(n*t);
y = @(t) a*sin(n*t)*cos(I);
z = @(t) a*sin(n*t)*sin(I);
r = [x(tVec') y(tVec') z(tVec')];

%
% Configure simulation components
%
opts = struct( ...
    ... % Simulation sampling frequency
    'SimulationFreq', fSim, ...
    ... % Flight computer sampling frequency
    'ComputerFreq', fCom, ... 
    ... % Standard deviation of measurement noise (sun, horizon, magnetometer)
    'MeasurementNoise', [deg2rad(0.05), deg2rad(0.015) deg2rad(0.5)], ...
    ... % Gyro measurement noise - angle random walk (ARW) (deg / sqrt(sec))
    'GyroNoise', deg2rad(0.45)/60, ...
    ... % Gyro bias rate (deg / sec / sec)
    'GyroBiasRateNoise', deg2rad(4)/3600/1000, ...
    ... % Initial attitude error covariance
    'AttitudeError', deg2rad(5), ...
    ... % Initial bias error covariance
    'GyroBiasError', deg2rad(0.02), ...
    ... % Actuator type
    'ActuatorType', 'rcs', ...
    ... % Actuator noise (Newton-meters)
    'ActuatorNoise', 0.05, ...
    ... % Configure RCS
    'Rcs', struct( ...
        'Thrust', 0.5 ./ diag(J), ...
        'UpperLim', 0.08, ...
        'LowerLim', 0.04, ...
        'Deadband', 1 ...
    ), ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-9 ...
);

%
% Initialize and preallocate variables
%

% Reference angular velocity
wRef = [0 0 -n]';

% Initalize true quaternion
qTrue = zeros(Nt, 4);
qTrue(1, :) = eul2quat(-I, pi, 0);

% Reference quaternion (pointing toward Earth)
[~, qRef] = ode45(@(~, q) qmult([wRef/2; 0], q), tVec, qTrue(1, :)', opts);

% Track time history of theta and theta dot for phase plane plot
dtheta = zeros(Nt, 3);
dthetadot = zeros(Nt, 3);

% Calculate matrices used in linearized dynamics
A = J \ (skew(J*wRef) - skew(wRef)*J);
phiA = @(t, t0) expm(A * (t - t0));
Ad = phiA(dt, 0);
Bd = integral(@(u) phiA(dt, u), 0, dt, 'ArrayValued', true) / J;

attEstError = zeros(Nt, 3, numTrials);
biasEstError = zeros(Nt, 3, numTrials);
onTimeRatio = zeros(numTrials, 1);

%
%% Main simulation loop
%
for jj = 1:numTrials
    % Initialize true angular velocity
    wTrue = zeros(Nt, 3);
    wTrue(1, :) = wRef;

    % Angular velocity estimate
    wEst = zeros(Nt, 3);
    wEst(1, :) = wRef;
    
    % Initialize gyro estimate (initial estimate is zero)
    betaEst = zeros(Nt, 3);
    betaEst(1, :) = [0 0 0];
    
    % Initialize quaternion estimate with some random error
    qEst = zeros(Nt, 4);
    dtheta0 = opts.AttitudeError * rand([3 1]);
    qEst(1, :) = qnormalize(qmult([dtheta0/2; 1], qTrue(1, :)'));

    % Initialize true gyro bias
    betaTrue = zeros(Nt, 3);
    betaTrue(1, :) = opts.GyroBiasError * randn([3 1]);
    
    % Track time history of control law
    u = zeros(Nt, 3);
    
    % Initialize covariance
    P = zeros(6, 6, Nt);
    P(:, :, 1) = blkdiag(opts.AttitudeError^2 * eye(3), opts.GyroBiasError^2 * eye(3));
    
    for ii = 1:Nt
        % True values
        t = tVec(ii);
        w = wTrue(ii, :)';
        q = qTrue(ii, :)';
        beta = betaTrue(ii, :)';
        rI = r(ii, :)';
        T = quat2rotm(q)';

        % Calculate magnetic field at time t
        rMag = magfield(rI);

        % Should the computer fire?
        if ~mod(t - t0 - dt, 1/opts.ComputerFreq)
            % Calculate reference horizon sensor measurement
            rEarth = -rI / norm(rI);

            % Put all measurements into a matrix and normalize
            ref = normc([rSun rEarth rMag]);

            % Get sensor measurements
            [wTilde, qTilde, R] = sensors(w, q, ref, beta, opts);

            % Estimate angular velocity and attitude
            [wHat, qHat, P(:, :, ii)] = nav( ...
                wTilde, ...
                qTilde, ...
                R, ...
                qEst(ii-1, :)', ...
                wEst(ii-1, :)', ...
                betaEst(ii-1, :)', ...
                P(:, :, ii-1), ...
                opts ...
            );

            % Calculate control
            [u(ii, :), dtheta(ii, :), dthetadot(ii, :)] = controller(wHat, qHat, wRef, qRef(ii, :)', opts);

            % Update estimate time history
            wEst(ii, :) = wHat';
            betaEst(ii, :) = (wTilde - wHat)';
            qEst(ii, :) = qHat';
        else
            % If computer doesn't fire, control is zero
            u(ii, :) = zeros(1, 3);

            if ii > 1
                % Computer doesn't update, so use same estimate as last epoch
                wEst(ii, :) = wEst(ii-1, :);
                betaEst(ii, :) = betaEst(ii-1, :);
                qEst(ii, :) = qEst(ii-1, :);
                P(:, :, ii) = P(:, :, ii-1);
            end
        end

        if ii < Nt
            % Convert control law into a torque
            tau = actuator(u(ii, :)', J, opts);

            % Calculate external torques in body frame
            M = sum([ ...
                gravgrad(J, -T*rI, n), ...           % gravity gradient
                magmoment * (T*rMag), ...           % magnetic torque
                tau ...                                     % control torque
                ], 2);

            % Propagate angular velocity using linearized dynamics
            dw = Ad*(w - wRef) + Bd*M;
            wTrue(ii+1, :) = (wRef + dw)';

            % Propagate attitude quaternion
            qTrue(ii+1, :) = qpropagate(q, w, dt);

            % Update bias random walk
            betaTrue(ii+1, :) = beta + (opts.GyroBiasRateNoise*dt) * randn([3 1]);
            
            % Calculate attitude estimate error
            qe = qmult(qTrue(ii, :), qconj(qEst(ii, :)));
            attEstError(ii, :, jj) = 2*qe(1:3) / qe(4);           
        end
    end
    
    % Thruster on-time percentage
    onTimeRatio(jj) = sum(any(u, 2)) / Nt * 100;
    
    % Bias estimate error for this MC run
    biasEstError(:, :, jj) = betaTrue - betaEst;
end

%% Plot results
close all;

%% Plot true bias and bias estimate
figure;
plot(tVec, rad2deg(betaTrue), tVec, rad2deg(betaEst), '--');
h = legend('$\beta_1$', '$\beta_2$', '$\beta_3$', '$\hat{\beta}_1$', '$\hat{\beta}_2$', '$\hat{\beta}_3$');
set(h, 'Interpreter', 'latex');
xlabel('seconds');
ylabel('deg / sec');
title('Gyro bias (true vs. estimated)');

% Plot bias estimate error
figure;
hold on;
plot(tVec, rad2deg(betaTrue - betaEst));
sigma = zeros(Nt, 3);
for ii = 1:Nt
    d = diag(P(4:6, 4:6, ii));
    sigma(ii, :) = sqrt(d);
end
plot(tVec, rad2deg(3*sigma), 'r--', tVec, rad2deg(-3*sigma), 'r--');
xlabel('seconds');
ylabel('deg / sec');
title('Bias estimate error');
legend('\beta_1', '\beta_2', '\beta_3');

%% Plot true angular velocity vs reference
figure;
plot(tVec, rad2deg(wTrue), [tVec(1) tVec(end)], rad2deg([wRef wRef]'), '--');
h = legend('$\omega_1$', '$\omega_2$', '$\omega_3$', '$\bar{\omega}_1$', '$\bar{\omega}_2$', '$\bar{\omega}_3$');
set(h, 'Interpreter', 'latex');
xlabel('seconds');
ylabel('deg / sec');
title('True angular velocity vs reference');

%% Plot true attitude vs reference
% Convert to angles
eulTrue = zeros(Nt, 3);
eulRef = zeros(Nt, 3);
attE = zeros(Nt, 3);
sigma = zeros(Nt, 3);
for ii = 1:Nt
    [eulTrue(ii, 1), eulTrue(ii, 2), eulTrue(ii, 3)] = quat2eul(qTrue(ii, :));
    [eulRef(ii, 1), eulRef(ii, 2), eulRef(ii, 3)] = quat2eul(qRef(ii, :));
    dq = qmult(qTrue(ii, :), qconj(qEst(ii, :)));
    attE(ii, :) = 2*dq(1:3) / dq(4);
    d = diag(P(1:3, 1:3, ii));
    sigma(ii, :) = sqrt(d);
end

% Plot true attitude vs reference
figure;
plot(tVec, rad2deg(unwrap(eulTrue)), tVec, rad2deg(unwrap(eulRef)), '--');
xlabel('seconds');
ylabel('degrees');
h = legend('$\theta_1$', '$\theta_2$', '$\theta_3$', '$\bar{\theta}_1$', '$\bar{\theta}_2$', '$\bar{\theta}_3$');
set(h, 'Interpreter', 'latex');
title('True attitude vs. reference attitude');

% Plot attitude error
figure;
for k = 1:3
    subplot(3, 1, k), ...
        hold on, ...
        plot(tVec, rad2deg(attE(:, k))), ...
        plot(tVec, rad2deg(3*sigma(:, k)), 'r--', tVec, rad2deg(-3*sigma(:, k)), 'r--'), ...
        legend(['\theta_' num2str(k)]), xlabel('seconds'), ylabel('degrees');
    if k == 1
        title('Attitude Error (true vs. estimate)');
    end
end

%% Plot phase plane
plotphaseplane(rad2deg(dtheta(:, 1)), rad2deg(dthetadot(:, 1)), opts.Rcs.UpperLim, opts.Rcs.LowerLim, opts.Rcs.Deadband, u(:, 1));
plotphaseplane(rad2deg(dtheta(:, 2)), rad2deg(dthetadot(:, 2)), opts.Rcs.UpperLim, opts.Rcs.LowerLim, opts.Rcs.Deadband, u(:, 2));
plotphaseplane(rad2deg(dtheta(:, 3)), rad2deg(dthetadot(:, 3)), opts.Rcs.UpperLim, opts.Rcs.LowerLim, opts.Rcs.Deadband, u(:, 3));

%% Plot thruster time history
figure;
subplot(3, 1, 1);
plot(tVec, u(:, 1));
subplot(3, 1, 2);
plot(tVec, u(:, 2));
subplot(3, 1, 3);
plot(tVec, u(:, 3));