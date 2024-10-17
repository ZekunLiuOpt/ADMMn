% L1/2 shrinkage operator
%
% Input: 
%         X: a given matrix X
%    lambda: the shrinkage parameter 
%
% Output: 
%         Y: the L1/2-shrink matrix
%
% Written by Zekun Liu, 16/12/2023
%
% Reference:
% [1] Z. Xu, X. Chang, F. Xu and H. Zhang. 
%     L1/2 Regularization: A Thresholding Representation Theory and a Fast Solver. 
%     IEEE Transactions on Neural Networks and Learning Systems, 23 (2012), 1013–1027.
%
% Latest Revision: 17/10/2024


function Y = prox_Lfrac12(X, lambda)

t = 2 / 3;
m = size(X, 1);
n = size(X, 2);
Y = zeros(m, n);
T = find(X > 54^(1 / 3) / 4 * lambda^t);
Y(T) = t * X(T) .* (1 + cos(t * pi - t * acos((lambda / 8) ./ ((abs(X(T)) / 3).^(1 / t)))));

end
