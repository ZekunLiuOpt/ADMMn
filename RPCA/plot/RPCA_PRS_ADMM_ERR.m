% Existed algorithm for solving the RPCA problem
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
% Written by Zekun Liu, 18/12/2023
% 
% Reference:
% [1] J. Jian, P. Liu and X. Jiang. 
%     A Partially Symmetric Regularized Alternating Direction Method of Multipliers
%     for Nonconvex Multi-block Optimization. 
%     Acta Mathematica Sinica, Chinese Series, 64 (2021), 1005–1026.
% 
% Latest Revision: 17/10/2024


function err = RPCA_PRS_ADMM_ERR(M, theta, w, beta, L1, S1, T1)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

MaxIter = 3000;
err = zeros(MaxIter, 1);

r = -0.1;
s = 1.05;

err(1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

for k = 1 : MaxIter

    L = prox_NNfrac12((beta * (T - S) + Lam + L) / (beta + 1), 1 / (beta + 0.07)); 

    S = prox_L1((beta * (T - L) + Lam + S) / (beta + 1), theta / (beta + 0.07));

    Lam = Lam - r * beta * (L + S - T);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    Lam = Lam - s * beta * (L + S - T);

    err(k + 1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

end

end
