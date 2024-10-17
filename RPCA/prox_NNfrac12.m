% Half shrinkage operator
%
% Input: 
%         X: a given matrix X
%         t: the shrinkage parameter 
%
% Output: 
%         Y: the half-shrink matrix
%
% Written by Zekun Liu, 16/12/2023
%
% Reference:
% [1] Z. Xu, X. Chang, F. Xu and H. Zhang. 
%     L1/2 Regularization: A Thresholding Representation Theory and a Fast Solver. 
%     IEEE Transactions on Neural Networks and Learning Systems, 23 (2012), 1013–1027.
%
% Latest Revision: 17/10/2024


function Y = prox_NNfrac12(X, t)

[U, S, V] = svd(X, 'econ');
s = diag(S);
s = prox_Lfrac12(s, t);
Y = U * diag(s) * V';

end

