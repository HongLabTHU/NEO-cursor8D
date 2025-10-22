%% 长期解码离线准确率
clear;clc

% 将tfr_for_mat.py添加到python目录
if count(py.sys.path, 'E:/pycharm/MyCode/NEO-BCI/neuracle-offline/dataloaders/') == 0
    insert(py.sys.path, int32(0), 'E:/pycharm/MyCode/NEO-BCI/neuracle-offline/dataloaders/');
end

% 将 python NEO_data 类导入matlab
import py.tfr_for_mat.NEO_data
NEO_loader = py.tfr_for_mat.NEO_data();
rootDir = 'E:\pycharm\NEOdata\TT01';


sM = Return_filelist(rootDir, 'single-');
dM = Return_filelist(rootDir, 'dual-');
Dirs = [sM(38:end); dM(36:end)];
% Dirs([5 6 9 14 25]) = [];
clear sM dM

for i = 1:length(Dirs)
    D(i) = days(datetime(Dirs{i}(51:58), "Format", "yyyyMMdd") - datetime('2024-11-11'));
end



%% 计算特征
option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;
option.fmax = 100;
option.maxnff = 256;
option.tps = [-2000, 100, 1100];

Spectra = [];
Ylab = [];

for i = 1:length(Dirs)
    try
        [Spect, ylab, fb] = PSD_dataSet(Dirs(i), NEO_loader, option);
        
        Dirs{i, 2} = (length(Ylab)+1):(length(Ylab)+length(ylab));
        Spectra = cat(3, Spectra, Spect);
        Ylab = [Ylab; ylab];
    catch
        disp('ERRRRRRRRRRRRRRRRRRRRROOOOORRRR');
        disp(Dirs{i})
    end
end


%% 解码验证
model = load('.\neo_vis\TT_LDA_1111.mat');
LBmat = [1 2 3 4 7 8 10 11 12 13 15:21;
    22:25 27:39];

for h = 1:17
    idx = [Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
    idx = idx(:);
    X = permute(Spectra(:, :, idx), [3 1 2]);
    X = log10(X(:, :));
    ylab = Ylab(idx);
    [macroF1(h), macroACC(h)] = longterm_vaild(model, X, ylab);
end


%%
% Days
idx = [];
for h = 3:6
    idx = [idx Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
end
X = permute(Spectra(:, :, idx), [3 1 2]);
X = log10(X(:, :));
ylab = Ylab(idx);
[~, macroACC1, C1] = longterm_vaild(model, X, ylab);

% Days
idx = [];
for h = 14:17
    idx = [idx Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
end
X = permute(Spectra(:, :, idx), [3 1 2]);
X = log10(X(:, :));
ylab = Ylab(idx);
[~, macroACC2, C2] = longterm_vaild(model, X, ylab);



%%
figure('Position', [573,428.3333333333333,386,329]);
idx = [1 5 9 4 7 2 6 3 8];
imagesc(C1(idx, idx), [0 0.8]);
colormap(slanCM('Greys'));
xticks(1:9); yticks(1:9);
Tbox = {'Rest', 'LE', 'LE' , 'RE', 'RE', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
xticklabels(Tbox(idx));
yticklabels(Tbox(idx));
colorbar('limits',  [0 0.8], 'Ticks', [0 0.8]);
annotation('textbox', [.89 .65 .1 .05], 'String', 'Confusion', 'EdgeColor', 'none', 'Rotation', 270, 'FontSize', 11);



%%
figure('Position', [573,520.3333333333333,380,237.3333333333333])
% 绘制散点图（自定义样式）
D1 = D(LBmat(1,:));
s = scatter(D1, macroF1, 20, 'filled', ...
    'MarkerFaceColor', [0.5 0.5 0.5], ...  % 青蓝色填充
    'MarkerEdgeColor', [0.5 0.5 0.5], ...  % 深蓝色边框
    'MarkerFaceAlpha', 0.4);               % 半透明效果
hold on;

% 计算线性回归
p = polyfit(D1, macroF1, 1);                % 1次多项式拟合
yfit = polyval(p, D1);                      % 拟合值

% 绘制回归线（红色实线）
plot(D1, yfit, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% 添加R²和方程标注
R2 = corr(D1', macroF1')^2;
text(0.05, 0.95, sprintf('R^2 = %.3f', R2), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', [1 1 1 0.7]);       % 半透明白底

box off
axis([-1 200 0 0.6]);
xticks([0 50 100 150 200]);
xlabel('Decoder time span (Days)', 'FontName', 'Arial');
ylabel('Macro F1-Score', 'FontName', 'Arial');



%%
[F1, macroACC1, ~] = longterm_vaild(model, X, ylab)


%%
function [macroF1, macroACC, C] = longterm_vaild(model, X, ylab)
test_X = bsxfun(@rdivide,  bsxfun(@minus, X, model.meanX), model.stdX);
[pred, ~] = predict(model.classer, test_X);
C = confusionmat(ylab, pred);
[macroF1, macroACC] = computeMacroF1AndACC(C);
C = C ./ repmat(sum(C, 2), 1, size(C, 1));
end

