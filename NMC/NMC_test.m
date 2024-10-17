% Experiments that comparing the performance of different algorithms on NMC problem
%
% Written by Zekun Liu, 06/02/2024
%
% Latest Revision: 17/10/2024


clear
clc

m = 500;
n = 500;
RANK = [2; 10; 20];
SR = [0.7; 0.5; 0.3];
CNT = 20;

w = 1;
beta = 1;
gamma = 1.618;

Timesaver1 = zeros(3, 3, CNT);
Timesaver2 = zeros(3, 3, CNT);
Timesaver3 = zeros(3, 3, CNT);

Errsaver1 = zeros(3, 3, CNT);
Errsaver2 = zeros(3, 3, CNT);
Errsaver3 = zeros(3, 3, CNT);

Chgsaver1 = zeros(3, 3, CNT);
Chgsaver2 = zeros(3, 3, CNT);
Chgsaver3 = zeros(3, 3, CNT);

Itersaver1 = zeros(3, 3, CNT);
Itersaver2 = zeros(3, 3, CNT);
Itersaver3 = zeros(3, 3, CNT);

for i = 1:3
    for j = 1:3
        for k = 1:CNT
        
            r = RANK(i);
            d = ceil(m * n * SR(j));

            X = rknonnegmatrix(m, n, r);  
            P = location(m, n, d);   
            M = P .* X; 
        
            alpha = 1.91 * 10^(-4) * norm(M, 'fro') * max(m, n) / r;
            beta2 = n * alpha / m;
        
            [X1, Y1, Z1, chg1, iter1, time1] = NMC_ADMMn(M, P, r, w, beta);
        
            [X2, Y2, Z2, U, V, chg2, iter2, time2] = NMFC_ADM(M, P, r, alpha, beta2, gamma);

            [X3, Y3, Z3, chg3, iter3, time3] = NMC_pADMMz(M, P, r, w, beta);
        
            err1 = norm(X - X1, 'fro') / (norm(X, 'fro') + 1); 
            err2 = norm(X - X2 * Y2, 'fro') / (norm(X, 'fro') + 1);
            err3 = norm(X - X3, 'fro') / (norm(X, 'fro') + 1);

            Timesaver1(i, j, k) = time1;
            Timesaver2(i, j, k) = time2;
            Timesaver3(i, j, k) = time3;

            Errsaver1(i, j, k) = err1;
            Errsaver2(i, j, k) = err2;
            Errsaver3(i, j, k) = err3;
            
            Chgsaver1(i, j, k) = chg1;
            Chgsaver2(i, j, k) = chg2;
            Chgsaver3(i, j, k) = chg3;
            
            Itersaver1(i, j, k) = iter1;
            Itersaver2(i, j, k) = iter2;
            Itersaver3(i, j, k) = iter3;
        end
    end
end

AveTime1 = mean(Timesaver1, 3);
AveTime2 = mean(Timesaver2, 3);
AveTime3 = mean(Timesaver3, 3);

AveErr1 = mean(Errsaver1, 3);
AveErr2 = mean(Errsaver2, 3);
AveErr3 = mean(Errsaver3, 3);

AveChg1 = mean(Chgsaver1, 3);
AveChg2 = mean(Chgsaver2, 3);
AveChg3 = mean(Chgsaver3, 3);

AveIter1 = mean(Itersaver1, 3);
AveIter2 = mean(Itersaver2, 3);
AveIter3 = mean(Itersaver3, 3);

save test3
