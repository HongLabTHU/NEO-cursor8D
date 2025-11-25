%% TFR
clear; clc;

if count(py.sys.path, 'your/python/path') == 0
    insert(py.sys.path, int32(0), 'your/python/path');
end

import py.tfr_for_mat.NEO_data
NEO_loader = py.tfr_for_mat.NEO_data();

rootDir = 'TT01_data';
matchingDirs = Return_filelist(rootDir, 'sing');


%% compute TFR 
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
    
    TFR = cat(1, TFR, tfr);
    ii
end
tfr = squeeze(mean(TFR, 1));


%% imshow TFR
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


%% 
function matchingDirs = Return_filelist(rootDir,searchString)
allDirs = dir(fullfile(rootDir, '**', '*'));
matchingDirs = {};
for i = 1:length(allDirs)
    if allDirs(i).isdir
        if contains(allDirs(i).name, searchString, 'IgnoreCase', true)
            matchingDirs{end+1} = fullfile(allDirs(i).folder, allDirs(i).name);
        end
    end
end

matchingDirs = matchingDirs';
end
