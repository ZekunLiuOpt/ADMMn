% Function that generates the matrix with p values 1 and others 0 which belongs to R^{m*n}
% Input: the number of nonzero entries p, the size of location matrix m and n
% Output: the location matrix P
% Written by: Zekun Liu (23/09/2023)
% Latest Revision: 20/09/2024


function P = location(m, n, p)

P = zeros(m, n);
ind = randi(numel(P), 1, p);
P(ind) = 1;

end
