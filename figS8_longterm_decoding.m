%% long-term decoding
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


%% 训练模型
classer = load('.\neo_vis\TT_LDA_1111.mat');
X = classer.classer.X;
ylab = classer.classer.Y;

X = bsxfun(@times, X, classer.stdX);
X = bsxfun(@plus, X, classer.meanX);
        
X = X(ismember(ylab, [100 101 102]), :);
ylab = ylab(ismember(ylab, [100 101 102]));

[X, meanX, stdX] = zscore(X, [], 1);
classer.model = fitcdiscr(X, ylab, 'DiscrimType', 'linear', 'Prior', 'uniform', 'Gamma', 0.5);
classer.meanX = meanX;
classer.stdX = stdX;
clear meanX  stdX X ylab


%% TFR 预测
clc
option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;
option.fmax = 100;
option.maxnff = 256;

DD = {'E:\pycharm\NEOdata\TT01\2024\11\26\Data\ECoG\TT01_20241126-single-MI',...
    'E:\pycharm\NEOdata\TT01\2025\02\21\Data\ECoG\TT01_20250221-single-MI',...
    'E:\pycharm\NEOdata\TT01\2025\03\14\Data\ECoG\TT01_20250314-single-MI',...
    'E:\pycharm\NEOdata\TT01\2025\05\23\Data\ECoG\TT01_20250523-single-MI'};

EE = {[31 111], [17 81], [7 109], [15 105]};


for ii = 1:4
    NEO_loader.neo_loadinto(DD{ii});
    events = int64(NEO_loader.events);
    data = double(NEO_loader.raw.get_data());
    data = NEO_reref(data/1e-6, 'average');

    for jj = 1:2
        Dots{ii, jj} = events(EE{ii}(jj):EE{ii}(jj)+5, 1) - events(EE{ii}(jj), 1) + 1000;
        Dots{ii, jj} = round(Dots{ii, jj}/200);
        NEO_loader.neo_get_tfr(events(EE{ii}(jj), 1), [int32(-1000), int32(17000)]);
        Power{ii, jj} = double(NEO_loader.power);

        Eve = events(EE{ii}(jj), 1)-1000:200:events(EE{ii}(jj), 1)+17000;
        Eve = [Eve' 0*Eve' 0*Eve'];
        [Espe, fb] = neo_calc_spectra(data, Eve, option);
        Espe = permute(Espe, [3 1 2]); Espe = log10(Espe(:, :));
        Espe = bsxfun(@minus, Espe, classer.meanX);
        Espe = bsxfun(@rdivide, Espe, classer.stdX);

        [~, Preds{ii, jj}] = predict(classer.model, Espe);
    end
    ii
end


%%

for i = 1:8
    subplot(8, 1, i);
    plot(Preds{i}(:, [2 3]), 'LineWidth', 1);
    hold on
    plot(Preds{i}(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    ylim([0 1.2]);
    box off
end


%%
Dx = 1;

Ax_x = [0.13 0.53;
        0.58 0.53;
        0.13 0.13;
        0.58 0.13];

AA = [];
for i = 1:4
    for j = 1:2
        AA = [AA squeeze(10*log10(Power{i, j}(3, :, :)))];
    end
end

meanA = mean(AA, 2); stdA = std(AA, [], 2);

marksC = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
figure('Position',[573,582.3333333333333,510,160])
for i = 1:4
    Ax(i) = axes('Position', [Ax_x(i, 1) Ax_x(i, 2) 0.4 0.35]);
    if ismember(i, [1 2])
        psD = squeeze(Power{Dx, i}(3, :, :));
        psD = 10*log10(psD);
        psD = (psD - meanA) ./ stdA;
        psD = conv2(psD, gausswin(200)'/sum(gausswin(200)), 'same');
        imagesc(Ax(i), flipud(psD), [-2 2]);
        colormap(Ax(i), slanCM('coolwarm'));
        xlim(Ax(i), [-50 18000]);
        yticks(Ax(i), [1 20 39]);
        yticklabels(Ax(i), {'200', '100', '0'});
        box(Ax(i), 'off')
        Ax(i).XAxis.Visible = 'off';
        ylabel(Ax(round(i/2)), 'freq (Hz)');
    end
    if ismember(i, [3 4])
        hold(Ax(i), 'on')
        plot(Ax(i), Preds{Dx, i-2}(:, 1), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        plot(Ax(i), Preds{Dx, i-2}(:, 2), 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
        plot(Ax(i), Preds{Dx, i-2}(:, 3), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1);
        scatter(Ax(i), Dots{Dx, i-2}(1:2:end), 1.3*ones(1, length(Dots{Dx, i-2}(1:2:end))), 9, marksC{i-2}, 'filled');
        scatter(Ax(i), Dots{Dx, i-2}(2:2:end), 1.3*ones(1, length(Dots{Dx, i-2}(2:2:end))), 9, [0.6 0.6 0.6], 'filled');
        axis(Ax(i), [0 91 -0.2 1.4]);
        box(Ax(i), 'off');
        yticks(Ax(i), [0 1]);
        yticklabels(Ax(i), {'0.0', '1.0'});
        xticks(Ax(i), 0:15:90);
        xticklabels(Ax(i), {'0', '3', '6', '9', '12', '15', '18'});
        ylabel(Ax(round(i/2)+1), 'Prob.');
    end
end



    
