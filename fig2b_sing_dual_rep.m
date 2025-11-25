%% Single and Dual movement representation

clear
clc

load('data\PSDdata.mat');
Preprocess_funy = @(Ylab) Ylab(2:2:end);
Preprocess_funx = @(Spectra) sqrt(Spectra(:, :, 2:2:end)) - sqrt(Spectra(:, :, 1:2:end));


%% maha dist mat
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

xticklabels([]);
yticklabels(Tbox([7 2 5 1 6 3 8 4]));

colorbar('Position',[0.920237617135192 0.109278499278499 0.025 0.2],...
    'Ticks',[0 2.7],...
    'Limits',[0 2.7], ...
    'Box','off');

annotation(fig,'textbox',...
    [0.882707496653265 0.0491919191919188 0.1 0.0500000000000001],...
    'VerticalAlignment','baseline',...
    'String','Distance',...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Cluster
Z = linkage(D(1:8, 1:8), 'ward');  % ward
figure('Position', [581,495.6,410,187.3]);
Z = dendrogram(Z,'ColorThreshold', 3.4, 'Reorder', optimalleaforder(Z,D(1:8, 1:8)));
set(Z, 'LineWidth', 1.2);
ylabel('Cluster Distance', 'FontSize',11);
xticks(1:10);
ylim([2.2 5]);
yticks([2.2 5]);
xticklabels(Tbox(str2num(xticklabels)'));
title('Hierarchical Clustering', 'FontSize',11);


