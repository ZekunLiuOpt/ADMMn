% Experiments that comparing the performance of different algorithms on NMC problem
%
% Written by Zekun Liu, 06/02/2024
%
% Latest Revision: 18/07/2025


clear
clc

m = 500;
n = 500;
RANK = 0.01 * m : 0.01 * m : 0.1 * m;
SR = [0.9, 0.7, 0.5, 0.3, 0.1];
CNT = 20;

w = 1;
beta = 1;
gamma = 1.618;

Timesaver1 = zeros(length(RANK), length(SR), CNT);
Timesaver2 = zeros(length(RANK), length(SR), CNT);
Timesaver3 = zeros(length(RANK), length(SR), CNT);

Errsaver1 = zeros(length(RANK), length(SR), CNT);
Errsaver2 = zeros(length(RANK), length(SR), CNT);
Errsaver3 = zeros(length(RANK), length(SR), CNT);

Chgsaver1 = zeros(length(RANK), length(SR), CNT);
Chgsaver2 = zeros(length(RANK), length(SR), CNT);
Chgsaver3 = zeros(length(RANK), length(SR), CNT);

Itersaver1 = zeros(length(RANK), length(SR), CNT);
Itersaver2 = zeros(length(RANK), length(SR), CNT);
Itersaver3 = zeros(length(RANK), length(SR), CNT);

for i = 1:length(RANK)
    r = RANK(i);
    for j = 1:length(SR)
        d = ceil(m * n * SR(j));
        for k = 1:CNT
        
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

save NMC500Time

figure;
plot(RANK, log10(AveErr1(:, 1)), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr1(:, 2)), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr1(:, 3)), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr1(:, 4)), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr1(:, 5)), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([-7 0]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
fontsize(12, "points")
xlabel('rank');
ylabel('log_{10}(RelErr)');

figure;
plot(RANK, log10(AveErr2(:, 1)), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr2(:, 2)), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr2(:, 3)), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr2(:, 4)), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr2(:, 5)), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([-7 0]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
fontsize(12, "points")
xlabel('rank');
ylabel('log_{10}(RelErr)');

figure;
plot(RANK, log10(AveErr3(:, 1)), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr3(:, 2)), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr3(:, 3)), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr3(:, 4)), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, log10(AveErr3(:, 5)), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([-7 0]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
fontsize(12, "points")
xlabel('rank');
ylabel('log_{10}(RelErr)');


figure;
plot(RANK, AveTime1(:, 1), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, AveTime1(:, 2), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, AveTime1(:, 3), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, AveTime1(:, 4), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, AveTime1(:, 5), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([min([min(min(AveTime1)), min(min(AveTime2)), min(min(AveTime3))]), ...
      max([max(max(AveTime1)), max(max(AveTime2)), max(max(AveTime3))])]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
fontsize(12, "points")
xlabel('rank');
ylabel('CPU Time (Seconds)');

figure;
plot(RANK, AveTime2(:, 1), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, AveTime2(:, 2), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, AveTime2(:, 3), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, AveTime2(:, 4), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, AveTime2(:, 5), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([min([min(min(AveTime1)), min(min(AveTime2)), min(min(AveTime3))]), ...
      max([max(max(AveTime1)), max(max(AveTime2)), max(max(AveTime3))])]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
xlabel('rank');
ylabel('CPU Time (Seconds)');
fontsize(12, "points")

figure;
plot(RANK, AveTime3(:, 1), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(RANK, AveTime3(:, 2), '-*', 'color', '#f97306', 'Linewidth', 2);
hold on;
plot(RANK, AveTime3(:, 3), '-x', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(RANK, AveTime3(:, 4), '-o', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(RANK, AveTime3(:, 5), '-+', 'color', '[0.93 0.69 0.13]', 'Linewidth', 2);
hold off;
grid on;
xlim([min(RANK) max(RANK)])
ylim([min([min(min(AveTime1)), min(min(AveTime2)), min(min(AveTime3))]), ...
      max([max(max(AveTime1)), max(max(AveTime2)), max(max(AveTime3))])]);
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('sr 90%', 'sr 70%', 'sr 50%', 'sr 30%', 'sr 10%', 'Location', 'northwest');
fontsize(12, "points")
xlabel('rank');
ylabel('CPU Time (Seconds)');
