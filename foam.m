function [ A, B, t ] = foam( b, r, w )
%FOAM Find optimal attitude matrix using FOAM method
%   A = FOAM(B,R,W) finds the optimal 3-by-3 attitude matrix A from a 
%   N-by-3 set of body-frame measurements B and N-by-3 reference frame 
%   measurements R and a N-by-1 vector of weights W that minimizes the cost
%   function in Wahba's problem.
%
%   See also DAVQ, SVDATT, QUEST, ESOQ.

if length(w) < 2
    error('Minimum of 2 measurements required');
end

tic;
B = bsxfun(@times, w, b) * r';

alpha = norm(B, 'fro')^2;
beta = det(B);
gamma = norm(adj(B), 'fro')^2;

k = @(x) (x^2 - alpha)/2;
zeta = @(x) k(x)*x - beta;
f = @(x) (x^2 - alpha)^2 - 8*x*beta - 4*gamma;
df = @(x) 8*zeta(x);

lambda = newtraph(f, df, sum(a), 1e-6);

A = ((k(lambda) + alpha)*B + lambda*adj(B') - (B*B')*B)/zeta(lambda);
t = toc;

end

