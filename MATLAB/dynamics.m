function [ ydot ] = dynamics( ~, y, J, M )
%DYNAMICS Propagate dynamics.
%   y: state vector
%       y(1:3)      :   angular velocity
%   J: inertia matrix
%   M: external torque

ydot = zeros(size(y));

w = y(1:3);
q = quat(y(4:7));

ydot(1:3) = J \ (skew(-w)*J*w + M);
ydot(4:7) = quat2vec(qmult(quat(w/2), q));

end

