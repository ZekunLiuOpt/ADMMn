% Experiments that investigating the performance of tested algorithms under different parameter settings on RPCA
%
% Written by Zekun Liu, 17/04/2025
%
% Latest Revision: 18/07/2025


clear
clc

p = 100;
n = 100;
spr = 0.1;
rank = 10;
CNT = 20;

Theta = [0.01, 0.05, 0.1, 0.5, 1] / sqrt(p);
W = [1e2, 5e2, 1e3, 5e3, 1e4];
beta1 = 3.11; % PRS-ADMM uses it
beta2 = 7.09; % IPPS-ADMM uses it, and its paper uses it for PRS-ADMM

Timesaver1 = zeros(length(Theta), length(W), CNT);
Timesaver2 = zeros(length(Theta), length(W), CNT);
Timesaver3 = zeros(length(Theta), length(W), CNT);
Timesaver4 = zeros(length(Theta), length(W), CNT);

Errsaver1 = zeros(length(Theta), length(W), CNT);
Errsaver2 = zeros(length(Theta), length(W), CNT);
Errsaver3 = zeros(length(Theta), length(W), CNT);
Errsaver4 = zeros(length(Theta), length(W), CNT);

Chgsaver1 = zeros(length(Theta), length(W), CNT);
Chgsaver2 = zeros(length(Theta), length(W), CNT);
Chgsaver3 = zeros(length(Theta), length(W), CNT);
Chgsaver4 = zeros(length(Theta), length(W), CNT);

Itersaver1 = zeros(length(Theta), length(W), CNT);
Itersaver2 = zeros(length(Theta), length(W), CNT);
Itersaver3 = zeros(length(Theta), length(W), CNT);
Itersaver4 = zeros(length(Theta), length(W), CNT);

for i = 1:length(Theta)
    theta = Theta(i);
    for j = 1:length(W)
        w = W(j);
        parfor k = 1:CNT % use "parfor" to save time
            L = randn(p, rank) * randn(rank, n); 
            q = randperm(p * n);
            S = zeros(p, n); 
            K = round(spr * p * n); 
            S(q(1:K)) = randn(K, 1);
            sigma = 0.01; 
            N = randn(p, n) * sigma;
            T = L + S; 
            M = T + N;
            
            [L1, S1, T1, chg1, iter1, time1] = RPCA_ADMMn(M, theta, w, beta1);
            
            [L2, S2, T2, chg2, iter2, time2] = RPCA_ADMMz(M, theta, w, beta1);
            
            [L3, S3, T3, chg3, iter3, time3] = RPCA_IPPS_ADMM(M, theta, w, beta2);

            [L4, S4, T4, chg4, iter4, time4] = RPCA_PRS_ADMM(M, theta, w, beta1);
            
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

save testParam

% plot the heapmap of error and time

figure;
h = heatmap(W, Theta * sqrt(p), AveErr1, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveErr1(:)), min(AveErr2(:)), min(AveErr3(:)), min(AveErr4(:))]), ...
                 max([max(AveErr1(:)), max(AveErr2(:)), max(AveErr3(:)), max(AveErr4(:))])];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.4f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveErr2, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveErr1(:)), min(AveErr2(:)), min(AveErr3(:)), min(AveErr4(:))]), ...
                 max([max(AveErr1(:)), max(AveErr2(:)), max(AveErr3(:)), max(AveErr4(:))])];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.4f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveErr3, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveErr1(:)), min(AveErr2(:)), min(AveErr3(:)), min(AveErr4(:))]), ...
                 max([max(AveErr1(:)), max(AveErr2(:)), max(AveErr3(:)), max(AveErr4(:))])];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.4f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveErr4, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveErr1(:)), min(AveErr2(:)), min(AveErr3(:)), min(AveErr4(:))]), ...
                 max([max(AveErr1(:)), max(AveErr2(:)), max(AveErr3(:)), max(AveErr4(:))])];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.4f'; 


figure;
h = heatmap(W, Theta * sqrt(p), AveIter1, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveIter1(:)), min(AveIter2(:)), min(AveIter3(:)), min(AveIter4(:))]), ...
                 3e3];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.0f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveIter2, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveIter1(:)), min(AveIter2(:)), min(AveIter3(:)), min(AveIter4(:))]), ...
                 3e3];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.0f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveIter3, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveIter1(:)), min(AveIter2(:)), min(AveIter3(:)), min(AveIter4(:))]), ...
                 3e3];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.0f'; 

figure;
h = heatmap(W, Theta * sqrt(p), AveIter4, interpreter = 'latex');
h.XLabel = '$\omega$';
h.YLabel = '$\rho \cdot \sqrt{m}$';
h.Colormap = parula; 
h.ColorLimits = [min([min(AveIter1(:)), min(AveIter2(:)), min(AveIter3(:)), min(AveIter4(:))]), ...
                 3e3];
h.FontSize = 12;
h.FontName = 'Arial';
h.CellLabelFormat = '%.0f'; 
