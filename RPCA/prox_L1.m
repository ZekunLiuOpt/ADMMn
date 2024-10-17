% Solf shrinkage operator
%
% Input: 
%         X: a given matrix
%         t: the shrinkage parameter 
%
% Output: 
%         Y: the soft-shrink matrix 
%
% Written by Zekun Liu, 16/12/2023
%
% Reference:
% [1] I. Daubechies, M. Defrise and C. D. Mol. 
%     An iterative thresholding algorithm for linear inverse problems with a sparsity constraint. 
%     Communications on Pure and Applied Mathematics, 57 (2004), 1413–1457.
%
% Latest Revision: 17/10/2024


function Y = prox_L1(X, t)

Y = sign(X) .* max(abs(X) - t, 0);

end
