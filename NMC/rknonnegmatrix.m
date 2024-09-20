% Function that generates rank k matrix which belongs to R^{m*n}
% Input: the rank k, the size of the nonnegative matrix m and n
% Output: the nonnegative rank k matrix Y
% Written by: Zekun Liu (14/12/2023)
% Latest Revision: 20/09/2024


function Y = rknonnegmatrix(m, n, k)

P = rand(m, k);
Q = rand(n, k);
Y = P * (Q');

end
