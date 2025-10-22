%%
clear
clc

option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;  % 201;
option.fmax = 150;  % 100;
option.maxnff = 512;  % 256;
option.tps = [-1000 200];

fileuse = 8:53;
load('E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\readme.mat')
data_root = 'E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\';

[Spectra, Ylab, fb] = Step0_psd_dataset(Label, data_root, fileuse, option);
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
   
    % run('E:\MATLAB_softhub\brainstorm3\toolbox\misc\cmap_mandrill.m');
    colormap(slanCM('coolwarm'));
    xlim([0 150])
    xticks(0:50:150)
    yticks(1:8)
    box off

    title(Title{conds}, 'FontSize', 9.5);
end


colorbar('Ticks', [-0.3 -0.2 -0.1 0 0.1 0.2 0.3], 'Position',[0.93 0.2 0.03 0.6]);
annotation('textbox', [.93 .14 .05 .05], 'String', 'RMS Amp (a.u.)', 'FontSize', 9, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% 创建 line
annotation('line',[0.155 0.155], [0.288 0.136], 'LineWidth', 1);

% 创建 line
annotation('line',[0.155 0.317], [0.136 0.136], 'LineWidth', 1);

% 创建 textbox
annotation('textbox',...
    [0.16 0.0741445086705175 0.23 0.0549132947976878],...
    'String',{'frequence (Hz)'},...
    'FontName','Arial',...
    'EdgeColor','none');

% 创建 textbox
annotation('textbox',...
    [0.142 0.142 0.12062937062937 0.0549132947976877],...
    'String',{'channel'},...
    'Rotation',90,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');



%% feak freq distribution
ylab = Preprocess_funy(Ylab);
spectra = Preprocess_funx(Spectra);
condu = [101 102 103 104];
spectra = spectra(:, :, ismember(ylab, condu));
ylab = ylab(ismember(ylab, condu));


Color = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];

H = zeros(size(ylab, 1), 8);
for ch = 1:8
    for i = 1:length(ylab)
        [pks, locs] = findpeaks(squeeze(spectra(25:51, ch, i)));
        [pk, loc] = max(pks);
        H(i, ch) = fb(24 + locs(loc));
    end
    
end
H = {H(ylab==101, :), H(ylab==102, :), H(ylab==103, :), H(ylab==104, :)};


x1 = H{1}(:, 2);
x2 = H{2}(:, 2);

% Create a figure
figure1 = figure('Position', [573,539.6,225.3,217.86]);

ax1 = axes('Position', [0.13 0.14 0.39 0.73]);
h1 = histogram(ax1, x1, 50:2:80, 'Normalization', 'pdf', 'FaceColor', Color(1, :), 'EdgeColor', 'none'); % Blue color for x1 histogram
hold on;

% Fit a Gaussian distribution to the datad
[f1, xi1] = ksdensity(x1);
plot(ax1, xi1, f1, "Color", Color(1, :), 'LineWidth', 1.5); 
yticks([0 0.1]);
ylabel('Left elbow', 'FontSize', 10, 'Color', Color(1, :)); %, 'FontName', 'Arial');
xlabel('RMS Amp. peak freq (Hz)', 'FontSize', 10);  %, 'FontName', 'Arial');
axis(ax1, [50 85 0 0.1]);
view(-90, 90) %# Swap the axes
% set(gca, 'YDir', 'reverse');
box off
hold off;


% Create the histogram for x2
ax2 = axes('Position', [0.56 0.14 0.39 0.73]);
h2 = histogram(ax2, x2, 50:2:80, 'Normalization', 'pdf', 'FaceColor', Color(2, :), 'EdgeColor', 'none'); % Red color for x2 histogram
hold on;

% Fit a Gaussian distribution to the data
[f2, xi2] = ksdensity(x2);
plot(ax2, xi2, f2, "Color", Color(2, :), 'LineWidth', 1.5); % Blue color for Gaussian fit
hold off;
box off
axis(ax2, [50 85 0 0.1]);
yticks([0 0.1]);
xticklabels([]);
ylabel('Left hand', 'FontSize', 10, 'Color', Color(2, :)); %, 'FontName', 'Arial');
view(-90, 90) %# Swap the axes
set(gca, 'YDir', 'reverse');

% 创建 textbox
annotation(figure1,'textbox',...
    [0.4586,0.9244,0.1598,0.103],...
    'String', 'ch8',...
    'HorizontalAlignment','center',...
    'FontSize', 10,...
    'FontWeight', 'bold', ...
    'FontName', 'Arial', ...
    'FitBoxToText','off',...
    'EdgeColor','none');

% 创建 line
annotation(figure1,'line',[0.325 0.325], [0.854 0.89]);
annotation(figure1,'line',[0.755 0.755], [0.854 0.89]);
annotation(figure1,'line',[0.325 0.755], [0.89 0.89]);
annotation(figure1,'textbox',...
    [0.4586,0.896,0.1598,0.06],...
    'String', '***',...
    'HorizontalAlignment','center',...
    'FontSize', 10,...
    'FontWeight', 'bold', ...
    'FontName', 'Arial', ...
    'FitBoxToText','off',...
    'EdgeColor','none');



