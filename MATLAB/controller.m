function [ u, dtheta, dthetadot ] = controller( w, q, wref, qref, opts )
%CONTROLLER Calculate control law.
%   U = CONTROLLER(W,Q,WREF,QREF,OPTS) calculates the control law U based on
%   the current angular velocity W and attitude Q and the reference angular
%   velocity WREF and reference attitude QREF. OPTS is a struct containing
%   controller parameters.
%
%   [U,DTHETA,DTHETADOT] = CONTROLLER(W,Q,WREF,QREF,OPTS) also returns the
%   calculated attitude error and its derivative.

dw = w - wref;
dq = qmult(q, qconj(qref));
dtheta = 2*dq(1:3);
dthetadot = cross(-wref, dtheta) + dw;
u = rcs(dtheta, dthetadot, opts.Rcs);

end