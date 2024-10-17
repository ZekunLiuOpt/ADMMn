% Existed algorithm for solving the RPCA problem
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
% [1] X. Wang, H. Shao, P. Liu and T. Wu. 
%     An Inertial Proximal Partially Symmetric ADMM-Based Algorithm for 
%     Linearly Constrained Multi-Block Nonconvex Optimization Problems with Applications. 
%     Journal of Computational and Applied Mathematics, 420 (2023), 114821.
%
% Latest Revision: 17/10/2024


function [L, S, T, chg, iter, time] = RPCA_IPPS_ADMM(M, theta, w, beta)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

bL = L;
bS = S;

eps = 1e-7;
MaxIter = 3000;
r = 0.2;
s = 0.5;
alpha1 = 0.25;
alpha2 = 0.25;

tic;
for k = 1 : MaxIter

    Lold = L;
    Sold = S;
    Told = T;

    bL = L + alpha1 * (L - bL);
    bS = S + alpha2 * (S - bS);

    L = prox_NNfrac12((beta * (T - S) + Lam + bL) / (beta + 1), 1 / (beta + 1));

    S = prox_L1((beta * (T - L) + Lam + bS) / (beta + 1), theta / (beta + 1));

    Lam = Lam - r * beta * (L + S - T);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    Lam = Lam - s * beta * (L + S - T);

    chg = norm([L - Lold, S - Sold, T - Told], 'fro') / (norm([Lold, Sold, Told], 'fro') + 1);

    iter = k;

    if chg < eps
        break
    end

end
toc;

time = toc;

end
