%% Center-out mapping vector compare
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Center');

idxf = [37 38 39;
     31 35 36;
     27 28 33];

ITR = zeros(3, 9);
Fitts_ITR = zeros(3, 9);

for j = 1:3

    kk = 1;
    for i = 1:3
        load(Flist{idxf(j, i)});
        Behave = Control.Behave;
        for u = 1:3
            for v = 1:8
                T(u, v) = Behave(u, v).timecost;
            end
            T_ = T(u, :); T_(T_>8) = 8;
            p = sum(T(u, :)<=8) / 8;
            if p < 1
                itr = log2(8) + p * log2(p) + (1 - p) * log2((1-p)/(8-1));
            else
                itr = log2(8);
            end
            ITR(j, kk) = itr / mean(T_) * 60;
            Fitts_ITR(j, kk) = mean(log2((0.4+0.11)/0.11)./T_);
            kk = kk + 1;
        end
    end
end



%% Create figure
figure('Position', [573,430.3,370.6,327.3]);
AX = axes('Position', [0.16 0.16 0.8 0.8]);
% Group names
formattedGroupNames = {'Sing. Move.\newline        4D',...
                       'Sing. & Dual.\newline        4D',...
                       'Sing. & Dual.\newline        8D'};

% Colors for each group
colors = [
    0.5 0.5 0.5; % Lighter blue for Single 4D
    0.65, 0.65, .65;  % Blue for Sin.-Dual. 4D
    0.8500 0.3250 0.098    % Lightest blue for Sin.-Dual. 8D
];

% Generate boxplot
boxplot(AX,ITR', 'Colors', colors, 'Symbol', '', 'Widths', 0.5, 'Labels', {}, 'Whisker', Inf);

% Hold the plot to add scatter points and means
hold on;

% Overlay individual data points and mean markers for each group
for i = 1:size(ITR, 1)
    scatter(repmat(i, 1, size(ITR, 2)), ITR(i, :), 20, 'filled', 'MarkerFaceColor', colors(i, :));
    % plot(i, mean(ITR(i, :)), 'd', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i, :));
end

% Customize plot
ylabel("Hitting ITR (bpm)");
ylim([10, 70]);

% Remove the box outline and add grid
box off;

% Add significance annotation (example between group 1 and group 3)
text(2, 67.5, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
line([1, 3], [67, 67], 'Color', 'k', 'LineWidth', 1.5);
text(2.5, 64, '***', 'FontSize', 14, 'HorizontalAlignment', 'center');
line([2, 3], [63.5, 63.5], 'Color', 'k', 'LineWidth', 1.5);

% 添加带换行符的X轴标签
ax = gca;
ax.XTickLabel = repmat({''}, 1, 3);
for i = 1:3
    text(i, ax.YLim(1), formattedGroupNames{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end


[p, h, stats] = ranksum(ITR(1, :), ITR(2, :))
[p, h, stats] = ranksum(ITR(1, :), ITR(3, :))
[p, h, stats] = ranksum(ITR(2, :), ITR(3, :))




%% 跨月解码
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Center');


idxf = [15:19 27 28 33  54:56 61:63 64:66 70:72 76:78];
FittsMat = [];

k = 1;
for i = 1:length(idxf)
    load(Flist{idxf(i)});
    Behave = Control.Behave';
    Behave = Behave(:);

    for h = 1:length(Behave)
        if Behave(h).timecost < 8
            FittsMat(k) = 60 * log2((0.4 + 0.22)/0.22) / Behave(h).timecost;
        else
            FittsMat(k) = NaN;
        end
        k = k + 1;
    end
end

idxf = {1:5*24, 5*24+1:8*24, 8*24+1:11*24, 11*24+1:14*24, 14*24+1:17*24, 17*24+1:19*24, 19*24+1:21*24};
for i = 1:length(idxf)
    P(i) = sum(~isnan(FittsMat(idxf{i}))) / length(idxf{i});
end

g = [zeros(1, 5*24) 1*ones(1, 3*24) 2*ones(1, 3*24) 3*ones(1, 3*24) 4*ones(1, 3*24)  5*ones(1, 3*24)  6*ones(1, 3*24)];
g = g(~isnan(FittsMat));
FittsMat = FittsMat(~isnan(FittsMat));


%%
figure('Position', [573,483.7,442,260])
qing = [0.65 0.65 0.65];  % %[0.1250 0.5781 0.8047];
hong = [0.8500 0.3250 0.0980];
colororder([qing; hong]);

H = [0 1 2 3 4 5 6];

% 左侧
yyaxis left
boxchart(g, FittsMat, 'BoxFaceColor', qing, 'JitterOutliers', 'off', 'BoxWidth', 0.35, ...
    'BoxFaceAlpha', 0.2, 'MarkerStyle', 'none');
hold on;

g1 = unique(g);
for i = 1:7 
    Z(i) = mean(FittsMat(g==g1(i)));
end
plot(H, Z, 'x--', 'LineWidth', 0.6, 'Color', [0.7 0.7 0.7]);
box off

ylabel('Fitts ITR (bpm)');
ylim([0 90]);
ax = gca;
ax.YAxis(1).LineWidth = 1;
ax.FontWeight = 'normal';

% 右侧
yyaxis right
plot(H, P, '.-', 'MarkerSize', 15, 'Color', hong, 'LineWidth', 1);
ylim([0 1]);
ylabel('Hit rate');
ax = gca;
ax.YAxis(2).LineWidth = 1;
ax.FontWeight = 'normal';

xlim([-0.6 6.7])
xticks(H);
xticklabels({'0M', '1M', '2M', '3M', '4M', '5M', '6M'});
xlabel('Decoder time span');



%%
function matFiles = findMatFilesByName(rootDir, searchStr)
    % 查找文件名包含 searchStr 的 .mat 文件
    files = dir(fullfile(rootDir, ['**/*' searchStr '*.mat']));  % 例如 *data*.mat
    matFiles = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);
    matFiles = matFiles(:);
end


