function [ M ] = actuator( u, J, opts )
%ACTUATOR Actuate a control into a torque.
%   ACTUATOR(U,J,OPTS) creates a torque on a spacecraft with inertia matrix J
%   from the control U. OPTS is a struct containing actuator parameters.

if ~isstruct(opts)
    error('Invalid argument');
elseif ~isfield(opts, 'ActuatorType')
    error('Missing required field: ActuatorType');
elseif ~isfield(opts, 'ActuatorNoise')
    error('Missing required field: ActuatorNoise');
end

if any(u)
    switch (opts.ActuatorType)
        case 'rcs'
            % In a reaction control system, the control U represents the desired
            % angular acceleration. The control torque is therefore J*U, where J
            % is the inertia matrix of the spacecraft
            M = J*u + (u ~= 0) .* (opts.ActuatorNoise * randn([3 1]));
        otherwise
            M = zeros(3, 1);
    end
else
    M = zeros(3, 1);

end