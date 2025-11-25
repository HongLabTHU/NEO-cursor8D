%% Dual movement representation

clear; clc
load('data\PSDdata.mat');

freq_func = @(X, fb) X(fb>50, :, :);
Preprocess_funy = @(Ylab) Ylab;
Preprocess_funx = @(Spectra) bsxfun(@minus, sqrt(Spectra), mean(sqrt(Spectra(:, :, Ylab==100)), 3));


%% dPCA of single mode
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);

X = permute(X, [3 1 2]);
X = X(:, :);
[~, muX, stdX] = zscore(X);


condu = [101 102 103 104];
X = X(ismember(ylab, condu), :);
X = bsxfun(@rdivide, bsxfun(@minus, X, muX), stdX);
ylab = ylab(ismember(ylab, condu));
ylab_s = ylab;

[W_single, V_single, explVar_single] = NEO_single_DPCA(X', ylab, condu);


%% Project single mode on V_single
CMAP = slanCM('tab20c');
Cmp = [CMAP(1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(13*6+1, :)];
Cmp = [CMAP(13*6+1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(1, :)];
Y_s = X * V_single(:, 1:2);

figure('Position', [573,503.7,310,254]);
ax = axes('Position', [0.15 0.15 0.8 0.8]);

for i = 1:4
    scatter(ax, Y_s(ylab==condu(i), 1), Y_s(ylab==condu(i), 2), 10, Cmp(i, :), 'filled', ...
        'AlphaData', 0.9*ones(1, sum(ylab==condu(i))), 'MarkerFaceAlpha', 'flat');
    hold(ax, 'on');
end

legend('LE', 'LH', 'RE', 'RH', ...
    'Position',[0.722424051474638 0.167058674117348 0.220430107526882 0.242125984251968],...
    'FontName','Arial',...
    'FontSize', 8,...
    'EdgeColor','none',...
    'Color','none');

xticks([-18 0 30]); yticks([-10 0 10])
axis([-18 30 -10 10]);
xlabel(ax, 'Single PC1', 'FontName', 'Arial');
ylabel(ax, 'Single PC2', 'FontName', 'Arial');


%% Project dual mode on V_single
Cmps = repmat([225   162   164]/255, 4, 1);
Cmps = repmat([0.65 0.65 0.65], 4, 1);
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);
X = permute(X, [3 1 2]);
X = X(:, :);

condu = [121 122 125 126];
X = X(ismember(ylab, condu), :);
X = bsxfun(@rdivide, bsxfun(@minus, X, muX), stdX);
ylab = ylab(ismember(ylab, condu));


Tbox = {'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
fPos = {[0.06, 0.1, 0.225, 0.7], [0.29, 0.1, 0.225, 0.7], [0.52, 0.1, 0.225, 0.7], [0.75, 0.1, 0.225, 0.7]};

figure('Position', [447,409,900,240], 'Color', [1 1 1]);
for i = 1:4
    ax(i) = axes('Position', fPos{i});
    Y{i} = squeeze(X(ylab==condu(i), :) * V_single(:, 1:2));
    scatter(ax(i), Y{i}(:, 1), Y{i}(:, 2), 8, Cmps(i, :), 'filled', 'MarkerFaceAlpha', 0.85, 'LineWidth', 20);

    hold(ax(i), 'on');
    for j = [101 102 103 104]
        error_ellipse(Y_s(ylab_s==j, 1), Y_s(ylab_s==j, 2), ax(i), Cmp(j-100, :), 0.5);
    end
    
    axis(ax(i), [-12 30 -7 5]);
    xticks([-15 0 35]); yticks([-8 0 6]);
    axis off
    title(ax(i), Tbox{i}, 'FontName', 'Arial', 'FontSize', 11); %, 'Color', Cmps(i, :), 'FontWeight','bold');
end

