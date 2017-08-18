function [ A ] = triad( b, r )
%TRIAD Find optimal attitude matrix using Davenport's Q method
%   A = TRIAD(B, R) estimates the optimal attitude matrix A given
%   a 2-by-3 matrix B and a 2-by-3 matrix R.

V = zeros(3);
W = zeros(3);

V(:, 1) = r(1, :);
V(:, 2) = cross(r(1, :), r(2, :)) / norm(cross(r(1, :), r(2, :)));
V(:, 3) = cross(V(:, 1), V(:, 2));

W(:, 1) = b(1, :);
W(:, 2) = cross(b(1, :), b(2, :)) / norm(cross(b(1, :), b(2, :)));
W(:, 3) = cross(W(:, 1), W(:, 2));

A = W * V';

end

