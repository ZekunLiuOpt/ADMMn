% Our algorithm for solving the NMC problem
%
% Input: 
%         M: the observation matrix
%         P: the location matrix 
%         r: the upper bound of rank
%         w: the F-norm paramter
%      beta: the penalty parameter
%
% Output: 
%   X, Y, Z: the reconstructed nonnegative low-rank matrix 
%       chg: relative change
%      iter: number of iterations
%      time: running time
%
% Written by Zekun Liu, 06/02/2024
%
% Reference:
% [1] Z. Liu. 
%     An Extended ADMM for 3-Block Nonconvex Nonseparable Problems With Applications. 
%     arXiv:2402.02193.
%
% Latest Revision: 17/10/2024


function [X, Y, Z, chg, iter, time] = NMC_ADMMn(M, P, r, w, beta)

m = size(M, 1);
n = size(M, 2);
X = zeros(m, n);
Y = zeros(m, n);
Z = zeros(m, n);
Lam = zeros(m, n);
E = ones(m, n);

invA = E ./ (2 * P + (w + beta) * E);
Q = 2 * P .* M;

eps = 1e-6;
MaxIter = 3000;

tic;
for k = 1 : MaxIter

    Y = 1 / (w + beta) * max(2 * beta * X + (w - beta) * Z - Lam, 0);

    Z = invA .* (Q + (w - beta) * Y + 2 * beta * X - Lam);

    [U, S, V] = svds((Y + Z + Lam / beta) / 2, r);
    X = U * S * V';

    Z = invA .* (Q + (w - beta) * Y + 2 * beta * X - Lam);

    Lam = Lam - beta * (2 * X - Y - Z);

    chg = norm(P .* (M - X), 'fro') / (norm(M, 'fro') + 1);

    iter = k;

    if chg < eps
        break
    end

end
toc;

time = toc;

end
