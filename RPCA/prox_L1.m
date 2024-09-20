% L1-regularization proximal operator
% Input: a given matrix X, with the parameter t
% Output: the L1-shrink matrix Y
% Written by: Zekun Liu (16/12/2023)
% Reference:
% [1] I. Daubechies, M. Defrise and C. D. Mol. An iterative thresholding algorithm for linear inverse problems with a sparsity constraint. 
%     Communications on Pure and Applied Mathematics, 57 (2004), 1413–1457.
% Latest Revision: 20/09/2024


function Y = prox_L1(X,t)

Y = sign(X) .* max(abs(X) - t, 0);

end
