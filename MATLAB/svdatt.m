function [ A, B ] = svdatt( b, r, w )
%SVDATT Find optimal attitude matrix using SVD method.
%   A = SVDATT(B, R, W) finds the optimal 3-by-3 attitude matrix A from a 
%   N-by-3 set of body-frame measurements B and N-by-3 reference frame 
%   measurements R and a N-by-1 vector of weights W that minimizes the cost
%   function in Wahba's problem.
%
%   [A, Y] = SVDATT(B, R, W) also returns the computed measurement matrix Y.
%
%   See also SVD, DAVQ, QUEST, FOAM, ESOQ.

if size(b, 2) < 2
    error('Minimum of 2 measurements required');
end

B = bsxfun(@times, w, b) * r';
[U, ~, V] = svd(B);
Up = U * diag([1 1 det(U)]);
Vp = V * diag([1 1 det(V)]);

A = Up * Vp';

end

