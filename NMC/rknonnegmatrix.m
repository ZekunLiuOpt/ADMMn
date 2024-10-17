% Function that generates rank k matrix which belongs to R^{m*n}
%
% Input: 
%         k: the rank
%      m, n: the size of the nonnegative matrix
%
% Output: 
%         Y: the nonnegative rank k matrix 
%
% Written by Zekun Liu, 14/12/2023
%
% Latest Revision: 17/10/2024


function Y = rknonnegmatrix(m, n, k)

P = rand(m, k);
Q = rand(n, k);
Y = P * (Q');

end
