%% single movement anayl
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
load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\readme.mat')
data_root = 'E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\';

[Spectra, Ylab, fb] = Step0_psd_dataset(Label, data_root, fileuse, option);

% Preprocess_funy = @(Ylab) Ylab;
% Preprocess_funx = @(Spectra) bsxfun(@minus, sqrt(Spectra), mean(sqrt(Spectra(:, :, Ylab==100)), 3));
Preprocess_funy = @(Ylab) Ylab(2:2:end);
Preprocess_funx = @(Spectra) sqrt(Spectra(:, :, 2:2:end)) - sqrt(Spectra(:, :, 1:2:end));


%% ERDS
ylab = Preprocess_funy(Ylab);
spectra = Preprocess_funx(Spectra);
condu = [101 102 103 104];
spectra = spectra(:, :, ismember(ylab, condu));
ylab = ylab(ismember(ylab, condu));


x = fb';

Bands = {x(logical((x>7).*(x<=16))),...
    x(logical((x>14).*(x<=32))), ...
    x(logical((x>30).*(x<=51))), ...
    x(logical((x>50).*(x<200)))};
     
BDColor = [0.3176    0.5882    0.8039;
        0.4902    0.7961    0.8;
        0.5765    0.7686    0.4902;
        0.9805 .7031 .6797];

CMAP = slanCM('tab20c');
Color = [CMAP(13*6+1, :); CMAP(13*2+1, :); CMAP(13*4+1, :); CMAP(1, :)];

bdname = {'$\alpha$', '$\beta$', 'low $\gamma$', 'high $\gamma$'};

Cod = {ylab==101, ylab==103, ylab==102, ylab==104};

titl = {'Left elbow', 'Right elbow', 'Left hand', 'Right hand'};

k = 1;

h_fig = figure('Position', [281.67,383.6,560,470.0], 'Name', 'Fig2A');

upb = 0.25;
for i = 1:4
    ax(i) = subplot(2, 2, i);
    X = squeeze(mean(spectra(:, 2, Cod{i}), 2))';
    for k = 1:length(fb)
        H(k) = ttest(X(:, k), 0, 'Alpha', 0.001);
    end
    H(H==1) = -0.9865*upb;
    H(H==0) = -1e8;
    % 
    yem = mean(X);  %mean(X1 ./ X2) - 1;
    yes = std(X) / sqrt(size(X, 1)) * 2.576;

    
    for jj = 1:4
        bands = Bands{jj};
        yu = upb*ones(1, length(bands));
        yl = (0.75*upb)*ones(1, length(bands));
        fill([bands fliplr(bands)], [yu fliplr(yl)], BDColor(jj, :), 'linestyle', 'none', 'FaceAlpha',0.4);
        text(mean(bands(bands<150))-3*jj, 0.91*upb, bdname{jj}, 'FontSize',10, 'Interpreter', 'latex');
        hold on
    end

    yline(0, 'LineWidth', 1.5);
    plot(x, yem, 'color', Color(i, :), 'LineWidth', 1.5)
    xticks([0 25 50 75 100 125 150])
    fill([x fliplr(x)], [yem+yes fliplr(yem-yes)], Color(i, :), 'linestyle', 'none', 'FaceAlpha',0.3);

    plot(x, H, 'color', [0 0 0], 'LineWidth', 2)
    box off
    axis([5 150 -upb upb]);
    yticks([-upb 0 upb]);
    % legend(titl{i});
    title(titl{i}, 'FontSize', 11, 'FontName', 'Arial');
    xlabel('Freq (Hz)');
    ylabel('RMS amplitude (a.u.)');
    % ylabel('Power (10lg $\mu V^{2}$)', 'interpreter', 'Latex');

    hold(ax(i), 'off');
end

% saveas(h_fig, h_fig.Name, 'svg');



%% 地形图
figure1 = figure('Colormap', slanCM('coolwarm'), 'Position', [573,311.6,120,446]);
% 创建 axes
axes1 = axes('Parent', figure1);
axis off

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[0.5 1.5]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[0.5 1.5]);
box(axes1,'on');
axis(axes1,'ij');
% 设置其余坐标区属性
set(axes1, 'CLim',[-0.075 0.075], 'Layer','top');
% 创建 colorbar
colorbar(axes1,'Position',[0.4 0.2 0.2 0.6],...
    'Ticks',[0 0.025 0.050 0.075],...,
    'TickLength', 0, ...
    'AxisLocation','in',...
    'FontSize', 10, ...
    'Box', 'off', ...
    'Limits',[0 0.075]);

% 创建 textbox
annotation(figure1,'textbox',...
    [0.238674033149171 0.0338164251207729 0.584530386740332 0.143301127214171],...
    'String',{'RMS','amplitude','(a.u.)'},...
    'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','off',...
    'EdgeColor','none');




%% 频域对比
% 定义不同的颜色 (RGB格式)
colors = [0 0.4470 0.7410;   % 蓝色
          0.8500 0.3250 0.0980;   % 橙色
          0.9290 0.6940 0.1250;   % 黄色
          0.4940 0.1840 0.5560];  % 紫色

for ch = 1:8
    subplot(2, 8, ch);
    bd1 = logical((fb>58).*(fb<62));
    boxplot(squeeze(mean(spectra(bd1, ch, :))), ylab, 'Colors', colors, 'Symbol', '');
    ylim([-0.2 0.4])
    % 找到图中的箱子对象并设置颜色
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(5-j, :), ...
              'FaceAlpha', 0.5, 'EdgeColor', colors(5-j, :), 'linewidth', 0.9); % 颜色设置和边界颜色一致
    end
    
    % 隐藏离群点标记
    set(findobj(gca, 'Tag', 'Outliers'), 'Visible', 'off');
    %%%%%%%%%%%%%%%
    subplot(2, 8, ch+8);
    bd2 = logical((fb>63).*(fb<67));
    boxplot(squeeze(mean(spectra(bd2, ch, :))), ylab, 'Colors', colors, 'Symbol', '');
    ylim([-0.2 0.4])
    % 找到图中的箱子对象并设置颜色
    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(5-j, :), ...
              'FaceAlpha', 0.5, 'EdgeColor', colors(5-j, :), 'linewidth', 0.9); % 颜色设置和边界颜色一致
    end
    
    % 隐藏离群点标记
    set(findobj(gca, 'Tag', 'Outliers'), 'Visible', 'off');
end



