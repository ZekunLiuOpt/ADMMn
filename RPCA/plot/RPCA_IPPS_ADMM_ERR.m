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
% Written by Zekun Liu, 17/12/2023
%
% Reference:
% [1] X. Wang, H. Shao, P. Liu and T. Wu. 
%     An Inertial Proximal Partially Symmetric ADMM-Based Algorithm for 
%     Linearly Constrained Multi-Block Nonconvex Optimization Problems with Applications. 
%     Journal of Computational and Applied Mathematics, 420 (2023), 114821.
%
% Latest Revision: 17/10/2024


function err = RPCA_IPPS_ADMM_ERR(M, theta, w, beta, L1, S1, T1)

m = size(M, 1);
n = size(M, 2);
L = zeros(m, n);
S = zeros(m, n);
T = zeros(m, n);
Lam = zeros(m, n);

bL = L;
bS = S;

MaxIter = 3000;
err = zeros(MaxIter, 1);

r = 0.2;
s = 0.5;
alpha1 = 0.25;
alpha2 = 0.25;

err(1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

for k = 1 : MaxIter

    bL = L + alpha1 * (L - bL);
    bS = S + alpha2 * (S - bS);

    L = prox_NNfrac12((beta * (T - S) + Lam + bL) / (beta + 1), 1 / (beta + 1));

    S = prox_L1((beta * (T - L) + Lam + bS) / (beta + 1), theta / (beta + 1));

    Lam = Lam - r * beta * (L + S - T);

    T = 1 / (w + beta) * (w * M + beta * (L + S) - Lam);

    Lam = Lam - s * beta * (L + S - T);

    err(k + 1) = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L1, S1, T1], 'fro') + 1);

end

end
