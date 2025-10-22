%% Single and Dual movement representation

clear
clc

option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;  % 201;
option.fmax = 150;  % 100;
option.maxnff = 512;  % 256;
option.tps = [-1000 200];

fileuse = 8:53;
load('E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\readme.mat')
data_root = 'E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\';

[Spectra, Ylab, fb] = Step0_psd_dataset(Label, data_root, fileuse, option);

% Preprocess_funy = @(Ylab) Ylab;
% Preprocess_funx = @(Spectra) bsxfun(@minus, sqrt(Spectra), mean(sqrt(Spectra(:, :, Ylab==100)), 3));
Preprocess_funy = @(Ylab) Ylab(2:2:end);
Preprocess_funx = @(Spectra) sqrt(Spectra(:, :, 2:2:end)) - sqrt(Spectra(:, :, 1:2:end));


%%
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
Color = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980];

Title = {'LE', 'LH'};
indx = [101 102];

for conds = 1:length(indx)
    h_fig = figure('Position', [442,633,163,116]);

    imagesc(fb, 1:8, mean(X(:, :, ylab==indx(conds)), 3)', [-0.3, 0.3]);
    colormap(slanCM('coolwarm'));
    xlim([0 150])
    yticks([]);
    xticks([]);
    xlabel('Freq', 'FontSize', 10)
    ylabel('channel', 'FontSize', 10)

    title(Title{conds}, 'FontSize', 10, 'Color', Color(conds, :), 'FontWeight', 'bold');
end


%% 马氏距离矩阵
ylab = Preprocess_funy(Ylab);
X = permute(Preprocess_funx(Spectra), [3 1 2]);
X = X(:, :);

condu = [101 102 103 104 121 122 125 126];
Tbox = {'LE', 'LH', 'RE', 'RH', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};

X = X(ismember(ylab, condu), :);
% X = zscore(X, [], 1);
ylab = ylab(ismember(ylab, condu));

Xs = zeros(length(condu), size(X, 2));
for i = 1:1:length(condu)
    Xs(i, :) = mean(X(ylab==condu(i), :));
end

Xem = [Xs; X];

tic
D = pdist(Xem, 'mahalanobis');
D = squareform(D);
fig = figure('Position', [573,219,498,462]);
imagesc(D([7 2 5 1 6 3 8 4], [7 2 5 1 6 3 8 4]), [0 2.7]);
colormap(slanCM('gist_heat'))
box off
% xticklabels(Tbox([7 2 5 1 6 3 8 4]));
xticklabels([]);
yticklabels(Tbox([7 2 5 1 6 3 8 4]));
% 创建 colorbar
colorbar('Position',[0.920237617135192 0.109278499278499 0.025 0.2],...
    'Ticks',[0 2.7],...
    'Limits',[0 2.7], ...
    'Box','off');

% 创建 textbox
annotation(fig,'textbox',...
    [0.882707496653265 0.0491919191919188 0.1 0.0500000000000001],...
    'VerticalAlignment','baseline',...
    'String','Distance',...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% 进行层次聚类
Z = linkage(D(1:8, 1:8), 'ward'); % 使用完全连接聚类算法
figure('Position', [581,495.6,410,187.3]);
Z = dendrogram(Z,'ColorThreshold', 3.4, 'Reorder', optimalleaforder(Z,D(1:8, 1:8)));
set(Z, 'LineWidth', 1.2);
ylabel('Cluster Distance', 'FontSize',11);
xticks(1:10);
ylim([2.2 5]);
yticks([2.2 5]);
xticklabels(Tbox(str2num(xticklabels)'));
title('Hierarchical Clustering', 'FontSize',11);

