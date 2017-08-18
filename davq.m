function [ A, B, t ] = davq( b, r, w )
%DAVQ Find optimal attitude matrix using Davenport's Q method.
%   A = DAVQ(B,R,W) finds the optimal 3-by-3 attitude matrix A from a 
%   N-by-3 set of body-frame measurements B and N-by-3 reference frame 
%   measurements R and a N-by-1 vector of weights W that minimizes the cost
%   function in Wahba's problem.
%
%   See also SVDATT, FOAM, QUEST, ESOQ.

if length(w) < 2
    error('Minimum of 2 measurements required');
end

tic;
B = bsxfun(@times, w, b) * r';
z = unskew(B' - B);

S = B + B';
K = [S - eye(3)*trace(B), z; z.', trace(B)];
[V, D] = eig(K);

q = V(:, max(D) == max(max(D)));

A = quat2dcm(q);
t = toc;

end

