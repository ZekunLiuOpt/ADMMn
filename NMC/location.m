% Function that generates the matrix with p values 1 and others 0 which belongs to R^{m*n}
%
% Input: 
%         p: the number of nonzero entries
%      m, n: the size of location matrix
%
% Output: 
%         P: the location matrix 
%
% Written by Zekun Liu, 23/09/2023
%
% Latest Revision: 17/10/2024


function P = location(m, n, p)

P = zeros(m, n);
ind = randi(numel(P), 1, p);
P(ind) = 1;

end
