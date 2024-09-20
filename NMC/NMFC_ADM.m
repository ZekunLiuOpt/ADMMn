% Existed algorithm for solving the NMFC problem
% Input: the observation matrix M, the location matrix P, the upper bound of rank r, the F-norm paramters alpha and beta, the penalty parameter gamma
% Output: the reconstructed nonnegative low-rank matrix Z, fractorization matrices X, Y, U, V; relative change chg, number of iterations iter, running time
% Written by: Zekun Liu (19/12/2023)
% Reference:
% [1] Y. Xu, W. Yin, Z. Wen and Y. Zhang. An alternating direction algorithm for matrix completion with nonnegative factors. 
%     Frontiers of Mathematics in China, 7 (2012), 365–384.
% Latest Revision: 20/09/2024


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
