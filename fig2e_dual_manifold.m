%% dual movement regression
clear;clc
load('data\PSDdata.mat');

X = Spectra(logical((fb>50).*(fb<150)), :, :);
X = permute(X, [3 1 2]); 
X = sqrt(X(:, :));
baseline = mean(X(Ylab==100, :));
condu = [101 102 103 104 121 122 125 126];

ylab = Ylab;
X = X(ismember(ylab, condu), :);
ylab = ylab(ismember(ylab, condu));


%% linear combination
base = {[101 102], [101 103], [102 104], [103 104]};
move = {121, 122, 125, 126};

Param = cell(length(base)+2, length(move));
modelP = cell(length(base)+2, length(move));
Rsquared = cell(length(base)+2, length(move));


% baseline fitting
for j = 1:length(move)
    x_p = X(ismember(ylab, move{j}), :);
    [~, score,~,~,explained,~] = pca(x_p');
    % x_Self = [mean(x_p)' score(:, [1 2])];
    % x_Orth = [mean(x_p)' score(:, end-1:end)];
    x_Self = score(:, 1:2);
    x_Orth = score(:, end-1:end);
    
    for k = 1:size(x_p, 1)
        rg_xp = x_p(k, :)' - mean(x_p)';
        RSS_ = sum((x_p(k, :) - mean(x_p(k, :))).^2);
        RSS = sum((rg_xp - mean(rg_xp)).^2);
        mdl = fitlm(x_Self, rg_xp); %, 'Intercept', false);
        Param{1, j} = [Param{1, j}; [mdl.Coefficients{2,1} mdl.Coefficients{3,1}]];
        RS_ = 1 - (1 - mdl.Rsquared.Ordinary) * RSS / RSS_;
        Rsquared{1, j} = [Rsquared{1, j} RS_];
        modelP{1, j} = [modelP{1, j} mdl.ModelFitVsNullModel.Pvalue];

        mdl = fitlm(x_Orth, rg_xp); %, 'Intercept', false);
        RS_ = 1 - (1 - mdl.Rsquared.Ordinary) * RSS / RSS_;
        Rsquared{end, j} = [Rsquared{end, j} RS_];
        modelP{end, j} = [modelP{end, j} mdl.ModelFitVsNullModel.Pvalue];
    end
end


% different basis
for i = 1:length(base)
    x_a = mean(X(ismember(ylab, base{i}(1)), :))';
    x_b = mean(X(ismember(ylab, base{i}(2)), :))';

    for j = 1:length(move)
        x_p = X(ismember(ylab, move{j}), :);
        x_p = [mean(x_p); x_p];
        for k = 1:size(x_p, 1)
            mdl = fitlm([x_a, x_b], x_p(k, :)');%, 'Intercept', false);
            Param{i+1, j} = [Param{i+1, j}; [mdl.Coefficients{2,1} mdl.Coefficients{3,1}]];
            Rsquared{i+1, j} = [Rsquared{i+1, j} mdl.Rsquared.Ordinary]; % 获取拟合模型的 R-squared 值
            modelP{i+1, j} = [modelP{i+1, j} mdl.ModelFitVsNullModel.Pvalue];
        end
    end
    i
end

S = zeros(length(base)+2, length(move));
for i = 1:length(base)+2
    for j = 1:length(move)
    S(i, j) = mean(Rsquared{i, j});
    end
end


%% bar
Tbox = {'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
ind = [0.7 2 3 4.3] - 0.2;

Catch = [1 2 5 6; 
    1 3 4 6;
    1 3 4 6;
    1 2 5 6];
falpha = [0.1 0.5 0.1 0.1;
          0.1 0.5 0.1 0.1;
          0.1 0.1 0.5 0.1;
          0.1 0.1 0.5 0.1];
XLAB = {'Self', '\textbf{LE-LH}', 'RE-RH', 'Ortho';
        'Self', '\textbf{LE-RE}', 'LH-RH', 'Ortho';
        'Self', 'LE-RE', '\textbf{LH-RH}', 'Ortho';
        'Self', 'LE-LH', '\textbf{RE-RH}', 'Ortho'};

figure('Position', [177,502,1237,263]);
 for j = 1:size(Catch, 1)
    ax(j) = subplot(1, 4, j);
    for i = 1:size(Catch, 2)
        RRR = Rsquared{Catch(j, i), j};
        meanRR = mean(RRR);
        bar(ax(j), ind(i)+0.2, meanRR, 'FaceColor', [0.7 0.8 0.7], 'FaceAlpha', falpha(j, i), ...
            'EdgeColor', [0.6 0.8 0.6], 'LineWidth', 0.9);
        hold on
        V = RRR(modelP{Catch(j, i), j} <= 1);
        scatter(ax(j), 0.4*rand(1,length(V))+ind(i), V, 4, [0.6 0.6 0.6], 'filled', 'AlphaData', 0.9*ones(size(V)), 'MarkerFaceAlpha', 'flat');
        V = RRR(modelP{Catch(j, i), j} > 1);
        scatter(ax(j), 0.4*rand(1,length(V))+ind(i), V, 4, [0.6 0.6 0.6], 'filled', '>', 'AlphaData', 0.9*ones(size(V)), 'MarkerFaceAlpha', 'flat');
        ylim(ax(j), [0.6 1.03]);
        ylabel(ax(j), '$R^{2}$', "Interpreter", "latex");
        box(ax(j), 'off');
        xticks(ax(j), [0.7 2 3 4.3]);
        xticklabels(ax(j), XLAB(j ,:));
        set(gca, "TickLabelInterpreter", "latex", "FontSize", 10);
    end
    title(ax(j), Tbox{j}, 'FontName', 'Arial', 'FontSize', 11); %, "Interpreter", "latex");
 end


 %% test
TTmat = zeros(6, 4);
for j = 1:4
    for i = 1:6
        [~, TTmat(i, j), ~, ~] = ttest2(Rsquared{i, j}, Rsquared{j+1, j});
    end
end

plot(ax(1), [2 3], [0.96 0.96], 'k', 'LineWidth', 1);
text(ax(1), 2.5, 0.97, '***', 'HorizontalAlignment', 'center');
plot(ax(1), [2 4.3], [0.99 0.99], 'k', 'LineWidth', 1);
text(ax(1), 3.15, 1, '***', 'HorizontalAlignment', 'center');

plot(ax(2), [2 3], [0.96 0.96], 'k', 'LineWidth', 1);
text(ax(2), 2.5, 0.97, '***', 'HorizontalAlignment', 'center');
plot(ax(2), [2 4.3], [0.99 0.99], 'k', 'LineWidth', 1);
text(ax(2), 3.15, 1, '***', 'HorizontalAlignment', 'center');

plot(ax(3), [2 3], [0.96 0.96], 'k', 'LineWidth', 1);
text(ax(3), 2.5, 0.97, '***', 'HorizontalAlignment', 'center');
plot(ax(3), [3 4.3], [0.99 0.99], 'k', 'LineWidth', 1);
text(ax(3), 3.65, 1, '***', 'HorizontalAlignment', 'center');

plot(ax(4), [2 3], [0.96 0.96], 'k', 'LineWidth', 1);
text(ax(4), 2.5, 0.97, '***', 'HorizontalAlignment', 'center');
plot(ax(4), [3 4.3], [0.99 0.99], 'k', 'LineWidth', 1);
text(ax(4), 3.65, 1, '***', 'HorizontalAlignment', 'center');


%% scatter
Tbox = {'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'}; % 
A1 = {'LE', 'LE', 'LH', 'RE';
      'LH', 'RE', 'RH', 'RH'};

for i = 1:4
    figure('Position', [573,404,355,353])
    ax1 = axes("Position", [0.15 0.15 0.5 0.5]);
    ZZ = Rsquared{i+1, i} > 0.75;
    scatter(ax1, Param{i+1, i}(ZZ, 1), Param{i+1, i}(ZZ, 2), 4, Rsquared{i+1, i}(ZZ), 'filled');
    hold(ax1, "on")
    colormap([210 210 210; 190 190 190; 170 170 170; 150 150 150; 130 130 130]/255);
    scatter(ax1, Param{i+1, i}(1, 1), Param{i+1, i}(1, 2), 20, 'r', 'd', 'LineWidth', 1.2);
    grid on
    
    plot(ax1, [-3 3], [-3 3], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    text(2.2, 2.6, '$p<1\times10^{-6}$', 'HorizontalAlignment','center', "Interpreter", "latex");
    xlabel(ax1, ['$W_{' A1{1,i} '}$'], "Interpreter", "latex", 'FontSize', 11);
    ylabel(ax1, ['$W_{' A1{2,i} '}$'], "Interpreter", "latex", 'FontSize', 11);
    xline(ax1, 0);
    yline(ax1, 0);
    axis(ax1, [-1.5 3 -1.5 3]);
    title(Tbox{i}, 'FontName', 'Arial', 'FontSize', 10);

    ax2 = axes("Position", [0.29 0.29 0.5*sqrt(2) 0.5*sqrt(2)]);
    resp = (Param{i+1, i}(ZZ, 1) - Param{i+1, i}(ZZ, 2))/sqrt(2);
    histogram(ax2, resp, 20, 'FaceAlpha', 0.25, 'FaceColor', [0.7 0.7 0.7],...
        'EdgeColor','none', 'Normalization','probability');
    hold on;
    [kde, xi] = ksdensity(resp, 'Bandwidth', 0.1);
    plot(ax2, xi, kde/4, 'LineWidth', 0.8, 'Color', [0.6 0.6 0.6]);
    yline(0);
    axis(ax2, [-2.25 2.25 0 1.2]);
    axis(ax2, 'off')
    view(ax2, 45, 90);
end



