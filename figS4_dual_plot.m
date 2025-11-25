%% single dual movement
clear
clc

load('data\PSDdata.mat');
Preprocess_funy = @(Ylab) Ylab(2:2:end);
Preprocess_funx = @(Spectra) sqrt(Spectra(:, :, 2:2:end)) - sqrt(Spectra(:, :, 1:2:end));


%%
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);


Title = {'LE', 'LH', 'RE', 'RH', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
indx = [101 102 103 104 121 122 125 126];
subp = [1 8 9 4 2 3 5 6];

h_fig = figure('Position', [547,-23.6,572,461.3], 'Name', 'Fig3B');
for conds = 1:length(indx)
    subplot(3, 3, subp(conds));
    imagesc(fb, 1:8, mean(X(:, :, ylab==indx(conds)), 3)', [-0.3, 0.3]);

    colormap(slanCM('coolwarm'));
    xlim([0 150])
    xticks(0:50:150)
    yticks(1:8)
    box off

    title(Title{conds}, 'FontSize', 9.5);
end


colorbar('Ticks', [-0.3 -0.2 -0.1 0 0.1 0.2 0.3], 'Position',[0.93 0.2 0.03 0.6]);
annotation('textbox', [.93 .14 .05 .05], 'String', 'RMS Amp (a.u.)', 'FontSize', 9, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

annotation('line',[0.155 0.155], [0.288 0.136], 'LineWidth', 1);
annotation('line',[0.155 0.317], [0.136 0.136], 'LineWidth', 1);
annotation('textbox',...
    [0.16 0.0741445086705175 0.23 0.0549132947976878],...
    'String',{'frequence (Hz)'},...
    'FontName','Arial',...
    'EdgeColor','none');

annotation('textbox',...
    [0.142 0.142 0.12062937062937 0.0549132947976877],...
    'String',{'channel'},...
    'Rotation',90,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');




