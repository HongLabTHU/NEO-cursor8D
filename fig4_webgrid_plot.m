%% Center-out Task performance
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Web');
useid = 1:22;


%% 计算每个session的命中时和命中率
TimeMat = NaN(24, length(useid));
PathMat = NaN(24, length(useid));

for i = 1:length(useid)
    try
        load(Flist{useid(i)});
        Control.Tarlist = Control.Tarlist;
        Behave = Control.Behave(1:24);
        for h = 1:length(Behave)
            if Behave(h).timecost < 14
                TimeMat(h, i) = Behave(h).timecost;
                PathMat(h, i) = check_path(Behave(h).Path);
            end
        end
    catch
    end
   
end
        
clc
for i = 1:length(useid)
    Flist{useid(i), 2} = mean(TimeMat(:, i), 'omitnan');
    Flist{useid(i), 3} = sum(~isnan(TimeMat(:, i)))/24;
end


%% webgird 跨月命中率
idxf = [1:5, 6:8, 12:22];
for i = 1:length(idxf)
    L(i) = Flist{idxf(i), 3};
end
g1 = repmat({'1mo'},5,1);
g2 = repmat({'2mo'},3,1);
g3 = repmat({'3mo'},2,1);
g4 = repmat({'4mo'},3,1);
g5 = repmat({'5mo'},3,1);
g6 = repmat({'6mo'},3,1);
g = [g1; g2; g3; g4; g5; g6];

figure('Position', [573,529,224,223])
H = [1 2 3 4 5 6];
boxplot(L, g, 'Colors', [0.5529 0.6275 0.7961], 'Symbol', '', 'Widths', 0.5, ...
    'Labels', {}, 'Whisker', Inf, 'Positions', H);
hold on;

G = {1:5, 6:8, 9:10, 11:13, 14:16, 17:19};
Z = [];
for i = 1:length(G)
    scatter(repmat(H(i), 1, length(G{i})), L(G{i}), 20, 'filled', 'MarkerFaceColor', [0.5529 0.6275 0.7961]);
    Z(i) = mean(L(G{i}));
end
plot(H, Z, 'x--', 'LineWidth', 0.6, 'Color', [0.7 0.7 0.7]);

ylabel('Hit rate');
xlabel('Decoder time span');
axis([-0 6.5 0.6 1]);
box off


%% webgrid 画轨迹图
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Web');
% load(Flist{4}); slice = 18:22;
% load(Flist{6}); slice = 20:24;
load(Flist{8}); slice = 6:10;
load(Flist{6}); slice = 1:5;
% load(Flist{17}); slice = 5:9;
L = Control.Behave(slice);

tarpos = [Control.Mesh1(Control.Tarlist(slice)); Control.Mesh2(Control.Tarlist(slice))]';

clc

cMap = slanCM('Set2');
cMap = cMap(1:32:150, :);
% cMap = repmat(cMap(3, :), 5, 1);

figure('Position', [573,331.6666666666666,436,426])
for i = 1:9
    plot([0.1*i-0.5 0.1*i-0.5], [-0.40 0.40], 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
    hold on
    plot([-0.40 0.40], [0.1*i-0.5 0.1*i-0.5], 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
end

plot([-0.42 0.42], [-0.42 -0.42], 'Color', [0 0 0], 'LineWidth', 1);
plot([-0.42 0.42], [0.42 0.42], 'Color', [0 0 0], 'LineWidth', 1);
plot([-0.42 -0.42], [-0.42 0.42], 'Color', [0 0 0], 'LineWidth', 1);
plot([0.42 0.42], [-0.42 0.42], 'Color', [0 0 0], 'LineWidth', 1);

for i = 1:length(slice)
    scatter(tarpos(i, 1), tarpos(i, 2), 780, 'filled', 'LineWidth', 1, ...
        'MarkerFaceColor', cMap(i, :), ...
        'MarkerEdgeColor', cMap(i, :), ...
        'MarkerFaceAlpha', 0.9);
    transP = linspace(0.1, 0.6, length(L(i).Path));
    for j = 1:length(L(i).Path)
        scatter(L(i).Path(j, 1), L(i).Path(j, 2), 100, 'filled', 'LineWidth', 1, ...
            'MarkerFaceColor', cMap(i, :), 'MarkerFaceAlpha', transP(j), ...
            'MarkerEdgeColor', cMap(i, :), 'MarkerEdgeAlpha', transP(j));
    end
end
for i = 1:length(slice)
    text(tarpos(i, 1), tarpos(i, 2), num2str(i), 'FontSize', 13, 'Color', [1 1 1], ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

axis equal
axis off


%% Fitts ITR 统计
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'web');

k = 1;
for i = 1:5
    load(Flist{i});
    Control.Tarlist = Control.Tarlist;
    Behave = Control.Behave(1:24);
    for h = 1:length(Behave)
        if Behave(h).timecost < 14
            tarpos = [Control.Mesh1(Control.Tarlist(h)) Control.Mesh2(Control.Tarlist(h))];
            D(k) = norm(tarpos - Control.Behave(h).Path(1, :));
            FittsMat(k) = 60 * log2((D(k) + 0.16)/0.16) / Behave(h).timecost;
        else
            FittsMat(k) = NaN;
            D(k) = NaN;
        end
        k = k + 1;
    end
end

figure('Position', [573,507,394,250.6666666666666])
axes('Position', [0.14 0.15 0.65 0.8])
scatter(D, log10(FittsMat), [], [0.5529 0.6275 0.7961], 'filled', 'MarkerFaceAlpha', 0.5);
ylim([0.5 3]);
yticks([0.5 1 2 3]);
yticklabels({'', '1e1', '1e2', '1e3'});
box off
xlabel('​​Start-Target dist.​', 'FontName', 'Arial');
ylabel('Fitts ITR (bpm)', 'FontName', 'Arial');

axes('Position', [0.83 0.15 0.13 0.8])
h = histogram(log10(FittsMat), 20, ...
    'Orientation', 'horizontal', ...
    'FaceColor', [0.5529 0.6275 0.7961], ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none', ...
    'Normalization', 'pdf'); % 归一化为概率密度

hold on;

[kde, xi] = ksdensity(log10(FittsMat), 'Bandwidth', 0.1); % 可调整带宽
plot(kde, xi, 'LineWidth', 1.5, 'Color', [0.5529 0.6275 0.7961]); % 红色KDE曲线
axis([0 2.2 0.5 3]);

% 标记均值（红色虚线）
yline(mean(log10(FittsMat), 'omitnan'), '--', 'Color', [0.8529 0.6275 0.701], 'LineWidth', 1.5);
axis off




%% 跨月解码
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'web');

idxf = [1:5, 6:8, 12:22];
FittsMat = [];

k = 1;
for i = 1:length(idxf)
    load(Flist{i});
    Control.Tarlist = Control.Tarlist;
    Behave = Control.Behave(1:24);
    for h = 1:length(Behave)
        if Behave(h).timecost < 14
            tarpos = [Control.Mesh1(Control.Tarlist(h)) Control.Mesh2(Control.Tarlist(h))];
            D = norm(tarpos - Control.Behave(h).Path(1, :));
            FittsMat(k) = 60 * log2((D+0.16)/0.16) / Behave(h).timecost;
        else
            FittsMat(k) = NaN;
        end
        k = k + 1;
    end
end


g1 = ones(1, 5*24,1);
g2 = 2*ones(1, 3*24,1);
g3 = 3*ones(1, 2*24,1);
g4 = 4*ones(1, 3*24,1);
g5 = 5*ones(1, 3*24,1);
g6 = 6*ones(1, 3*24,1);
g = [g1 g2 g3 g4 g5 g6];


g = g(~isnan(FittsMat));
FittsMat = FittsMat(~isnan(FittsMat));
for i = 1:6 
    Z(i) = median(FittsMat(g==g1(i)));
end


figure('Position', [573,483.7,420,260])
H = [1 2 3 4 5 6];
boxchart(g, FittsMat, 'BoxFaceColor', [0.5529 0.6275 0.7961], 'JitterOutliers', 'off', ...
    'BoxFaceAlpha', 0.2, 'MarkerStyle', 'none');
axis([0 7 0 100]);
hold on;
g1 = unique(g);
plot(H, Z, 'o--', 'LineWidth', 0.6, 'Color', [0.7 0.7 0.7]);
box off
xticks(H);
xticklabels({'1mo', '2mo', '3mo', '4mo', '5mo', '6mo'});
xlabel('Decoder time span');
ylabel('Fitts ITR (bpm)');



%%
figure('Position', [573,483.7,420,260])
qing = [0.5529 0.6275 0.7961];  % %[0.1250 0.5781 0.8047];
hong = [0.8500 0.3250 0.0980];
colororder([qing; hong]);

H = [1 2 3 4 5 6];

% 左侧
yyaxis left
boxchart(g, FittsMat, 'BoxFaceColor', qing, 'JitterOutliers', 'off', 'BoxWidth', 0.35, ...
    'BoxFaceAlpha', 0.2, 'MarkerStyle', 'none');
hold on;

g1 = unique(g);
for i = 1:6 
    Z(i) = median(FittsMat(g==g1(i)));
end
plot(H, Z, 'o--', 'LineWidth', 0.6, 'Color', [0.7 0.7 0.7]);
box off
ylabel('Fitts ITR (bpm)');
ylim([0 91]);
ax = gca;
ax.YAxis(1).LineWidth = 1;
ax.FontWeight = 'normal';

% 右侧
yyaxis right
plot(H, [0.9000    0.9444    0.8542    0.9306    0.9583    0.9028], '.-', 'MarkerSize', 15, 'Color', hong, 'LineWidth', 1);
ylim([0 1]);
ylabel('Hit rate');
ax = gca;
ax.YAxis(2).LineWidth = 1;
ax.FontWeight = 'normal';

xlim([0.2 6.8])
xticks(H);
xticklabels({'1M', '2M', '3M', '4M', '5M', '6M'});
xlabel('Decoder time span');


%%
function matFiles = findMatFilesByName(rootDir, searchStr)
    % 查找文件名包含 searchStr 的 .mat 文件
    files = dir(fullfile(rootDir, ['**/*' searchStr '*.mat']));  % 例如 *data*.mat
    matFiles = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);
    matFiles = matFiles(:);
end



%%
function r = check_path(X)
X = max(abs(X), [], 2);
r = sum(X==0.42)/length(X);
end
