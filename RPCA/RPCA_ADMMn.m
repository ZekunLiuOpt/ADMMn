% Our algorithm for solving the RPCA problem
%
% Input: 
%         M: the observation matrix 
%     theta: the L1-norm parameter
%         w: the F-norm paramter
%      beta: the penalty parameter
%
% Output: 
%         S: the reconstructed sparse matrix
%         L: the low-rank matrix
%         T: the observation matrix
%       chg: relative change 
%      iter: number of iterations
%      time: running time
%
% Written by Zekun Liu, 17/12/2023
%
% Reference:
% [1] Z. Liu. 
%     An Extended ADMM for 3-Block Nonconvex Nonseparable Problems With Applications.
%     arXiv:2402.02193.
%
% Latest Revision: 17/10/2024


function [L, S, T, chg, iter, time] = RPCA_ADMMn(M, theta, w, beta)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

eps = 1e-7;
MaxIter = 3000;

tic;
for k = 1 : MaxIter

    Lold = L;
    Sold = S;
    Told = T;

    S = prox_L1(T + Lam / beta - L, theta / beta);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    L = prox_NNfrac12(T - S + Lam / beta, 1 / beta);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    Lam = Lam - beta * (L + S - T);

    chg = norm([L - Lold, S - Sold, T - Told], 'fro') / (norm([Lold, Sold, Told], 'fro') + 1);

    iter = k;

    if chg < eps
        break
    end

end
toc;

time = toc;

end
