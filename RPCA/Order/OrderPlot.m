% Experiments that recording the Relative Error of ADMMn and pADMMz with XY and YX-orders on RPCA
% 
% Written by Zekun Liu, 18/04/2025
%
% Notes: RPCA_ADMMn_Order_ERR, RPCA_pADMMz_Order_ERR are similar to RPCA_pADMMz_ERR, RPCA_pADMMz_ERR, 
%        with an additional "err" variance to record the change of relative error during iterations,
%        hence we omit their codes for simplicity
%
% Latest Revision: 18/07/2025


clear
clc

p = 100;
n = 100;

% noiseless
spr = 0.1;
rank = 10;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11; % PRS-ADMM uses it

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0; 
N = randn(p,n) * sigma;
T = L + S; 
M = T + N;

err11 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err12 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err13 = RPCA_ADMMn_Order_ERR(M, theta, w, beta1, L, S, T);

err14 = RPCA_pADMMz_Order_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err14), '-.', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(ind, log10(err12), '-.', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(ind, log10(err13), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(ind, log10(err11), '-', 'color', '#f97306', 'Linewidth', 2);
hold off;
grid on;
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('pADMMz-XY', 'pADMMz-YX', 'ADMMn-XY', 'ADMMn-YX', 'Location', 'east');
fontsize(12, "points")
xlabel('Iterations');
ylabel('log_{10}(RelErr)');


% noisy
spr = 0.05;
rank = 20;

theta = 0.1 / sqrt(p);
w = 1e3;
beta1 = 3.11; % PRS-ADMM uses it

L = randn(p, rank) * randn(rank, n); 
q = randperm(p * n);
S = zeros(p, n); 
K = round(spr * p * n); 
S(q(1:K)) = randn(K, 1);
sigma = 0.01; 
N = randn(p, n) * sigma;
T = L + S; 
M = T + N;

err21 = RPCA_ADMMn_ERR(M, theta, w, beta1, L, S, T);

err22 = RPCA_pADMMz_ERR(M, theta, w, beta1, L, S, T);

err23 = RPCA_ADMMn_Order_ERR(M, theta, w, beta1, L, S, T);

err24 = RPCA_pADMMz_Order_ERR(M, theta, w, beta1, L, S, T);

figure;
ind = 0:3000;
plot(ind, log10(err24), '-.', 'color', '#A066D3', 'Linewidth', 2);
hold on;
plot(ind, log10(err22), '-.', 'color', '#77AC30', 'Linewidth', 2);
hold on;
plot(ind, log10(err23), '-', 'color', '#0072BD', 'Linewidth', 2);
hold on;
plot(ind, log10(err21), '-', 'color', '#f97306', 'Linewidth', 2);
hold off;
grid on;
set(gca, 'FontName', 'Arial', 'Fontsize', 12);
legend('pADMMz-XY', 'pADMMz-YX', 'ADMMn-XY', 'ADMMn-YX', 'Location', 'east');
fontsize(12, "points")
xlabel('Iterations');
ylabel('log_{10}(RelErr)');

save Orderplot
