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
% Written by Zekun Liu, 18/12/2023
%
% Reference:
% [1] J. Jian, P. Liu and X. Jiang. 
%     A Partially Symmetric Regularized Alternating Direction Method of Multipliers 
%     for Nonconvex Multi-block Optimization. 
%     Acta Mathematica Sinica, Chinese Series, 64 (2021), 1005–1026.
%
% Latest Revision: 17/10/2024


function [L, S, T, chg, iter, time] = RPCA_PRS_ADMM(M, theta, w, beta)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

eps = 1e-7;
MaxIter = 3000;

r = -0.1;
s = 1.05;

tic;
for k = 1 : MaxIter

    Lold = L;
    Sold = S;
    Told = T;

    L = prox_NNfrac12((beta * (T - S) + Lam + L) / (beta + 1), 1 / (beta + 0.07)); 

    S = prox_L1((beta * (T - L) + Lam + S) / (beta + 1), theta / (beta + 0.07));

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
