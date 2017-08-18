function [ A, B, t ] = quest( b, r, w )
%QUEST Find optimal attitude matrix using QUEST method
%   A = QUEST(B,R,W) finds the optimal 3-by-3 attitude matrix A from a 
%   N-by-3 set of body-frame measurements B and N-by-3 reference frame 
%   measurements R and a N-by-1 vector of weights W that minimizes the cost
%   function in Wahba's problem.
%
%   See also DAVQ, SVDATT, FOAM, ESOQ.

if length(w) < 2
    error('Minimum of 2 measurements required');
end

tic;
B = bsxfun(@times, w, b) * r';
S = B + B';
z = unskew(B' - B);

s = trace(S)/2;
k = trace(adj(S));
D = det(S);

a = s^2 - k;
b = s^2 + z'*z;
c = D + z' * S * z;
d = z' * S^2 * z;

f = @(x) x^4 - (a + b)*x^2 - c*x + (a*b + c*s - d);
df = @(x) 4*x^3 - 2*(a + b)*x - c;

lambda = newtraph(f, df, 1);

alpha = lambda^2 - s^2 + k;
beta = lambda - s;
gamma = (lambda + s)*alpha - D;

X = (alpha*eye(3) + beta*S + S^2) * z;
q = (gamma^2 + norm(X)^2)^(-1/2) * [X; gamma];

A = quat2dcm(q);
t = toc;

end

