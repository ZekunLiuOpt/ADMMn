% Experiments that comparing the change of the relative error of different algorithms
% Written by: Zekun Liu (18/12/2023)
% Latest Revision: 20/09/2024


clear
clc

p = 100;
n = 100;

% noiseless case, (0.05，10)
spr = 0.05;
rank = 10;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11;
beta2 = 7.09; 

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0; 
N = randn(p, n) * sigma;
T = L + S; 
M = T + N;

err11 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err12 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err13 = RPCA_IPPS_ADMM_ERR(M, theta, w, beta1, L, S, T);

err14 = RPCA_PRS_ADMM_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err14), '-.', 'color', '#77AC30', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err12), '--', 'color', 'm', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err13), ':', 'color', '#EDB120', 'Linewidth', 1.75);
hold on;
plot(ind, log10(err11), '-', 'color', '#0072BD', 'Linewidth', 1.25);
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11);
h = legend('PSR-ADMM', 'pADMMz', 'IPPS-ADMM', 'ADMMn', 'Location', 'northeast');
set(h, 'FontSize', 8);
xlabel('Iterations');
ylabel('log_{10}(RelErr)');


% noiseless case, (0.05，20)
spr = 0.05;
rank = 20;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11; 
beta2 = 7.09; 

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0; 
N = randn(p, n) * sigma;
T = L + S; 
M = T + N;

err21 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err22 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err23 = RPCA_IPPS_ADMM_ERR(M, theta, w, beta1, L, S, T);

err24 = RPCA_PRS_ADMM_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err24), '-.', 'color', '#77AC30', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err22), '--', 'color', 'm', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err23), ':', 'color', '#EDB120', 'Linewidth', 1.75);
hold on;
plot(ind, log10(err21), '-', 'color', '#0072BD', 'Linewidth', 1.25);
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11);
h = legend('PSR-ADMM', 'pADMMz', 'IPPS-ADMM', 'ADMMn', 'Location', 'northeast');
set(h, 'FontSize', 8);
xlabel('Iterations');
ylabel('log_{10}(RelErr)');


% noisy case, (0.05，10)
spr = 0.05;
rank = 10;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11; 
beta2 = 7.09;

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0.01; 
N = randn(p, n) * sigma;
T = L + S; 
M = T + N;

err31 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err32 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err33 = RPCA_IPPS_ADMM_ERR(M, theta, w, beta1, L, S, T);

err34 = RPCA_PRS_ADMM_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err34), '-.', 'color', '#77AC30', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err32), '--', 'color', 'm', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err33), ':', 'color', '#EDB120', 'Linewidth', 1.75);
hold on;
plot(ind, log10(err31), '-', 'color', '#0072BD', 'Linewidth', 1.25);
hold off;
ylim([-3 0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11);
h = legend('PSR-ADMM', 'pADMMz', 'IPPS-ADMM', 'ADMMn', 'Location', 'northeast');
set(h, 'FontSize', 8);
xlabel('Iterations');
ylabel('log_{10}(RelErr)');


% noisy case, (0.05，20)
spr = 0.05;
rank = 20;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11;
beta2 = 7.09; 

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0.01; 
N = randn(p, n) * sigma;
T = L + S; 
M = T + N;

err41 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err42 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err43 = RPCA_IPPS_ADMM_ERR(M, theta, w, beta1, L, S, T);

err44 = RPCA_PRS_ADMM_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err44), '-.', 'color', '#77AC30', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err42), '--', 'color', 'm', 'Linewidth', 1.25);
hold on;
plot(ind, log10(err43), ':', 'color', '#EDB120', 'Linewidth', 1.75);
hold on;
plot(ind, log10(err41), '-', 'color', '#0072BD', 'Linewidth', 1.25);
hold off;
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11);
h = legend('PSR-ADMM', 'pADMMz', 'IPPS-ADMM', 'ADMMn', 'Location', 'northeast');
set(h, 'FontSize', 8);
xlabel('Iterations');
ylabel('log_{10}(RelErr)');

save test2
load test2
