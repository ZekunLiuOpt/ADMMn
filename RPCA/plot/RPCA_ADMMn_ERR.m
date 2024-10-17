% Our algorithm for solving the RPCA problem
% used to record the change of the relative error during iterations
%
% Input: 
%         M: the observation matrix
%     theta: the L1-norm parameter
%         w: the F-norm paramter 
%      beta: the penalty parameter
%        S1: the ground truth sparse matrix
%        L1: the ground truth low-rank matrix
%        T1: = L1 + S1
%
% Output: 
%       err: the change of the relative error during iterations
%
% Written by Zekun Liu, 17/12/2023
%
% Reference:
% [1] Z. Liu. 
%     An Extended ADMM for 3-Block Nonconvex Nonseparable Problems With Applications. 
%     arXiv:2402.02193.
%
% Latest Revision: 17/10/2024


function err = RPCA_ADMMn_ERR(M, theta, w, beta, L1, S1, T1)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

MaxIter = 3000;
err = zeros(MaxIter, 1);

err(1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

for k = 1 : MaxIter

    S = prox_L1(T + Lam / beta - L, theta / beta);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    L = prox_NNfrac12(T + Lam / beta - S, 1 / beta);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    Lam = Lam - beta * (L + S - T);

    err(k + 1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

end

end
