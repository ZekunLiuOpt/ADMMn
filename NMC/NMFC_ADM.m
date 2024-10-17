% Existed algorithm for solving the NMFC problem
%
% Input: 
%          M: the observation matrix
%          P: the location matrix 
%          r: the upper bound of rank
%      alpha: the F-norm paramter
%       beta: the F-norm paramter
%      gamma: the penalty parameter
%
% Output: 
%          Z: the reconstructed nonnegative low-rank matrix 
% X, Y, U, V: fractorization matrices
%        chg: relative change
%       iter: number of iterations
%       time: running time
%
% Written by Zekun Liu 19/12/2023
%
% Reference:
% [1] Y. Xu, W. Yin, Z. Wen and Y. Zhang. 
%     An alternating direction algorithm for matrix completion with nonnegative factors. 
%     Frontiers of Mathematics in China, 7 (2012), 365–384.
% 
% Latest Revision: 17/10/2024


function [X, Y, Z, U, V, chg, iter, time] = NMFC_ADM(M, P, r, alpha, beta, gamma)

m = size(M, 1);
n = size(M, 2);
X = zeros(m, r);
Y = rand(r, n);
Z = M;
U = zeros(m, r);
V = zeros(r, n);
Lam1 = zeros(m, r);
Lam2 = zeros(r, n);
I = eye(r);

eps = 1e-6;
MaxIter = 3000;

tic;
for k = 1 : MaxIter

    X = (Z * Y' + alpha * U - Lam1) / (Y * Y' + alpha * I);

    Y = (X' * X + beta * I) \ (X' * Z + beta * V - Lam2);

    Z = X * Y + P .* (M - X * Y);

    U = max(X + Lam1 / alpha, 0);

    V = max(Y + Lam2 / beta, 0);

    Lam1 = Lam1 + gamma * alpha * (X - U);

    Lam2 = Lam2 + gamma * beta * (Y - V);

    chg = norm(P .* (M - X * Y), 'fro') / (norm(M, 'fro') + 1);

    iter = k;

    if chg < eps
        break
    end

end
toc;

time = toc;

end
