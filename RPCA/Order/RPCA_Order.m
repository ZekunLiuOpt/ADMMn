% Experiments that comparing ADMMn and pADMMz with XY and YX-orders on RPCA
%
% Written by Zekun Liu, 18/04/2025
%
% Notes: RPCA_ADMMn_Order, RPCA_pADMMz_Order only change the update order of L and S,
%        hence we omit their codes for brevity
%
% Latest Revision: 18/07/2025


clear
clc

p = 100;
n = 100;
SPR = [0.05, 0.1];
RANK = [1, 5, 10, 20];
CNT = 20;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11; 

Timesaver1 = zeros(2, 4, CNT);
Timesaver2 = zeros(2, 4, CNT);
Timesaver3 = zeros(2, 4, CNT);
Timesaver4 = zeros(2, 4, CNT);

Errsaver1 = zeros(2, 4, CNT);
Errsaver2 = zeros(2, 4, CNT);
Errsaver3 = zeros(2, 4, CNT);
Errsaver4 = zeros(2, 4, CNT);

Chgsaver1 = zeros(2, 4, CNT);
Chgsaver2 = zeros(2, 4, CNT);
Chgsaver3 = zeros(2, 4, CNT);
Chgsaver4 = zeros(2, 4, CNT);

Itersaver1 = zeros(2, 4, CNT);
Itersaver2 = zeros(2, 4, CNT);
Itersaver3 = zeros(2, 4, CNT);
Itersaver4 = zeros(2, 4, CNT);

for i = 1:2
    spr = SPR(i);
    for j = 1:4
        rank = RANK(j);
        parfor k = 1:CNT
            L = randn(p, rank) * randn(rank, n); 
            q = randperm(p * n);
            S = zeros(p, n); 
            K = round(spr * p * n); 
            S(q(1:K)) = randn(K, 1);
            sigma = 0.01; % set sigma = 0 for noiseless case
            N = randn(p, n) * sigma;
            T = L + S; 
            M = T + N;
            
            [L1, S1, T1, chg1, iter1, time1] = RPCA_ADMMn(M, theta, w, beta1);
            
            [L2, S2, T2, chg2, iter2, time2] = RPCA_pADMMz(M, theta, w, beta1);
            
            [L3, S3, T3, chg3, iter3, time3] = RPCA_ADMMn_Order(M, theta, w, beta1);

            [L4, S4, T4, chg4, iter4, time4] = RPCA_pADMMz_Order(M, theta, w, beta1);
            
            err1 = norm([L - L1, S - S1, T - T1], 'fro') / (norm([L, S, T], 'fro') + 1);
            err2 = norm([L - L2, S - S2, T - T2], 'fro') / (norm([L, S, T], 'fro') + 1);
            err3 = norm([L - L3, S - S3, T - T3], 'fro') / (norm([L, S, T], 'fro') + 1);
            err4 = norm([L - L4, S - S4, T - T4], 'fro') / (norm([L, S, T], 'fro') + 1);

            Timesaver1(i, j, k) = time1;
            Timesaver2(i, j, k) = time2;
            Timesaver3(i, j, k) = time3;
            Timesaver4(i, j, k) = time4;
            
            Errsaver1(i, j, k) = err1;
            Errsaver2(i, j, k) = err2;
            Errsaver3(i, j, k) = err3;
            Errsaver4(i, j, k) = err4;
            
            Chgsaver1(i, j, k) = chg1;
            Chgsaver2(i, j, k) = chg2;
            Chgsaver3(i, j, k) = chg3;
            Chgsaver4(i, j, k) = chg4;
            
            Itersaver1(i, j, k) = iter1;
            Itersaver2(i, j, k) = iter2;
            Itersaver3(i, j, k) = iter3;
            Itersaver4(i, j, k) = iter4;
        end
    end
end

AveTime1 = mean(Timesaver1, 3);
AveTime2 = mean(Timesaver2, 3);
AveTime3 = mean(Timesaver3, 3);
AveTime4 = mean(Timesaver4, 3);

AveErr1 = mean(Errsaver1, 3);
AveErr2 = mean(Errsaver2, 3);
AveErr3 = mean(Errsaver3, 3);
AveErr4 = mean(Errsaver4, 3);

AveChg1 = mean(Chgsaver1, 3);
AveChg2 = mean(Chgsaver2, 3);
AveChg3 = mean(Chgsaver3, 3);
AveChg4 = mean(Chgsaver4, 3);

AveIter1 = mean(Itersaver1, 3);
AveIter2 = mean(Itersaver2, 3);
AveIter3 = mean(Itersaver3, 3);
AveIter4 = mean(Itersaver4, 3);

save testOrderNoisy
