%% 全数据分析 时频图
clear; clc;

% 添加 python 代码路径，这里要改成自己的tfr_for_mat.py文件的路径
if count(py.sys.path, 'E:/pycharm/MyCode/NEO-BCI/neuracle-offline/dataloaders/') == 0
    insert(py.sys.path, int32(0), 'E:/pycharm/MyCode/NEO-BCI/neuracle-offline/dataloaders/');
end


% 将 python NEO_data 类导入matlab
import py.tfr_for_mat.NEO_data
NEO_loader = py.tfr_for_mat.NEO_data();

% 定义根目录和要查找的字符串
rootDir = 'E:\pycharm\NEOdata\XB_data\';
searchString = 'sing';

% 使用dir递归查找所有文件夹
allDirs = dir(fullfile(rootDir, '**', '*'));

% 只保留文件夹，并检查是否包含目标字符串
matchingDirs = {};
for i = 1:length(allDirs)
    if allDirs(i).isdir
        if contains(allDirs(i).name, searchString, 'IgnoreCase', true)
            matchingDirs{end+1} = fullfile(allDirs(i).folder, allDirs(i).name);
        end
    end
end

for i = 1:length(matchingDirs)
    matchingDirs{i}(1:27) = [];
end
matchingDirs = matchingDirs(1:19)';

%%

TFR = [];
for ii = 1:length(matchingDirs)
    NEO_loader.neo_loadinto(matchingDirs{ii});    % 导入数据
    events = int64(NEO_loader.events);
    data = double(NEO_loader.raw.get_data());
    data = NEO_reref(data/1e-6, 'average');

    NEO_loader.neo_get_epotfr(-1, 4, '5');
    tfr = double(NEO_loader.tfr.get_data());
    times = double(NEO_loader.tfr.times);
    freqs = double(NEO_loader.tfr.freqs);

    % 删除被噪声严重污染的trial
    if ii == 3
        tfr(20, :, :, :) = [];
    end
    if ii == 13
        tfr(7, :, :, :) = [];
    end

    TFR = cat(1, TFR, tfr);
    ii
end
tfr = squeeze(mean(TFR, 1));


%% 噪声校验

for ii = [3 13]
    NEO_loader.neo_loadinto(matchingDirs{ii});    % 导入数据
    events = int64(NEO_loader.events);
    data = double(NEO_loader.raw.get_data());
    data = NEO_reref(data/1e-6, 'average');

    NEO_loader.neo_get_epotfr(-1, 4, '5');
    tfr = double(NEO_loader.tfr.get_data());
    times = double(NEO_loader.tfr.times);
    freqs = double(NEO_loader.tfr.freqs);
    figure
    for j = 1:20
        subplot(4, 5, j);
        imagesc(times, freqs, flipud(squeeze(tfr(j, 5, :, :))), [-1.5 1.5]);
        colormap(slanCM('coolwarm'));
    end

end



%%
C_lev = 1.5;
figure('Position', [513.7,147,662,724]);

APos = [0.06 ]

CHH = [1 3 5 7 2 4 6 8];
for i = 1:8
    subplot(4, 2, CHH(i));
    tfr_c = conv2(squeeze(tfr(i, :, :)), gausswin(100)'/sum(gausswin(100)), 'same');
    % TFR = squeeze(tfr(i, :, :));
    imagesc(times, freqs, flipud(tfr_c), [-C_lev C_lev]);
    colormap(slanCM('coolwarm'));
    xline(0, 'k--');
    xline(2.5, 'k--');
    axis([-1 4 0 200])  
    yticks([0 100 200]);
    yticklabels({'200', '100', '0'});
    title(['ch' num2str(i)]);
    xlabel('time (s)');
    ylabel('Freq (Hz)');
    box off
end
colorbar('Ticks', [-C_lev 0 C_lev], 'Position',[0.93 0.2 0.03 0.6], 'FontSize', 10);
sgtitle('Hand MI', 'FontSize', 12);
% annotation('textbox', [.96 .5 .1 .05], 'String', 'ERSP', 'EdgeColor', 'none', 'Rotation', 270);


