function [ A, B, t ] = esoq( b, r, w )
%ESOQ Find optimal attitude matrix using ESOQ method.
%   A = ESOQ(B,R,W) finds the optimal 3-by-3 attitude matrix A from a 
%   N-by-3 set of body-frame measurements B and N-by-3 reference frame 
%   measurements R and a N-by-1 vector of weights W that minimizes the cost
%   function in Wahba's problem.
%
%   See also SVDATT, DAVQ, FOAM, QUEST.

if length(w) < 2
    error('Minimum of 2 measurements required');
end

tic;
B = bsxfun(@times, w, b) * r';
z = unskew(B' - B);
S = B + B';
K = [S - eye(3)*trace(B), z; z.', trace(B)];

b = -2 * trace(B)^2 + trace(adj(B + B')) - z'*z;
c = -trace(adj(K));
d = det(K);

p = (b/3)^2 + 4*d/3;
q = (b/3)^3 - 4*d*b/3 + c^2/2;

u1 = 2 * sqrt(p) * cos((1/3)*acos(q/p^(3/2))) + b/3;

if length(a) > 2
    g1 = sqrt(u1 - b);
    g2 = 2*sqrt(u1^2 - 4*d);

    lambda = (1/2) * (g1 + sqrt(-u1 - b - g2));
elseif length(a) == 2
    g3 = sqrt(2*sqrt(d) - b);
    g4 = sqrt(-2*sqrt(d) - b);
    lambda = (g3 + g4)/2;
end

H = (K - lambda*eye(4));

% Compute diagonal of cofactor matrix of H
q = zeros(4, 1);
for k = 1:4
    q(k) = cof(H, k, k);
end

% The index of the farthest-from-zero value corresponds to the maximum
% modulus vector cross product
[~, k] = max(abs(q));

% Calculate the rest of the cross product
for i = 1:4
    if i ~= k
        q(i) = cof(H, i, k);
    end
end

% Normalize and convert to A
A = quat2dcm(normc(q));
t = toc;

end

