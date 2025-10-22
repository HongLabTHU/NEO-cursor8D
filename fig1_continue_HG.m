%% 投影联合模式

clear;
clc

tmin = -1; tmax = 4;
option.fs = 1000;
option.fb = 55:5:95;
% option.fb = 30:4:50;
option.tmin = tmin;
option.tmax = tmax;

load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\readme.mat')
data_root = 'E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\';
fileuse = [8 9 10 19:2:27 30 31 33 35 36:2:40 43:3:52];

[Scalingp, Ylab, time] = Step2_continue_dataset(Label, data_root, fileuse, option);

Preprocess_funy = @(Ylab) Ylab;
% Preprocess_funx = @(Scalingp) sqrt(Scalingp);
Preprocess_funx = @(Scalingp) bsxfun(@minus, sqrt(Scalingp), mean(sqrt(Scalingp(:, 1:1000, :, :)), 2));


%%  信号平均
ylab = Preprocess_funy(Ylab);
X = mean(Preprocess_funx(Scalingp), 4);

% 定义边界和间距参数
leftMargin = 0.05;    % 左侧留白
rightMargin = 0.03;   % 右侧留白
topMargin = 0.09;     % 顶部留白
bottomMargin = 0.18;   % 底部留白
horizontalGap = 0.065; % 子图水平间距
verticalGap = 0.05;   % 子图垂直间距

% 计算每个子图的宽度和高度
axWidth = (1 - leftMargin - rightMargin - 3*horizontalGap)/4;
axHeight = 1 - topMargin - bottomMargin;

CMAP = slanCM('tab20c');
Color = [CMAP(13*6+1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(1, :)];

figure('Position', [331,460.3333333333333,966.6666666666665,205.9999999999999]);
for i = [3 4 7 8]
    j = floor(i/2);
    ax(j) = axes('Position', [leftMargin + (j-1)*(axWidth + horizontalGap), bottomMargin, axWidth, axHeight]);
    for c = [101 102 103 104]
        plot(ax(j), [c c], [c c], 'Color', Color(c-100, :), 'LineWidth', 1);
        hold(ax(j), 'on')
    end
    for c = [101 102 103 104]
        % baseline = mean(X(ylab==c, 1:1000, i, :), "all") * ones(1, length(time));
        x_m = squeeze(mean(X(ylab==c, :, i)));
        x_v = squeeze(std(X(ylab==c, :, i)) / sqrt(sum(ylab==c))) * 2.576;
        plot(ax(j), time, squeeze(mean(X(ylab==c, :, i))), 'Color', Color(c-100, :), 'LineWidth', 1);
        fill(ax(j), [time fliplr(time)], [x_m-x_v fliplr(x_m+x_v)], Color(c-100, :), 'linestyle', 'none', 'FaceAlpha',0.2);
    end
    axis(ax(j), [tmin tmax -3 12]);
    xline(ax(j), 0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    yline(ax(j), 0);
    xline(ax(j), 2.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    yticks(ax(j), [-3 0 12]);
    xlabel(ax(j), 'time (s)', 'FontSize', 10, 'FontName', 'Arial');
    ylabel(ax(j), 'HG RMS Amp (a.u.)', 'FontSize', 10, 'FontName', 'Arial');
    title(ax(j), ['CH' num2str(i)], 'FontName', 'Arial');
    box(ax(j), 'off');
end

% legend({'LE', 'LH', 'RE', 'RH'}, 'Orientation', 'horizontal');



%% 
clear;
clc;


load('E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\readme_raw.mat');
data_root = 'E:\pycharm\MyCode\NEO-BCI\Piplinecode\XB_data\';
fileuse = 2:48;

fs = 1000;
tmin = -1;
tmax= 4;

BB = zeros(8, floor((tmax-tmin)*fs), length(fileuse)*80);
Ylab = [];

k = 1;
for id = 1:length(fileuse)
    % 导入数据
    neo_data = load([data_root Label(fileuse(id)).name]);

    % 平均重参考
    data = NEO_reref(neo_data.data/1e-6, 'average');
    events = Label(fileuse(id)).events;

    for j = 2:3:size(events, 1)
        data_block = data(:,  (events(j, 1)+ fs*tmin - 199):(events(j, 1)+ fs*tmax - 200));
        BB(:, :, k) = data_block;
        Ylab = [Ylab; events(j, 3)];
        k = k + 1;
    end

    
    id
end


%% MRCP LE LH
tt = tmin+1/fs:1/fs:tmax;
Color = [0 0.4470 0.7410;       % LE
        0.8500 0.3250 0.0980;   % LH
        0.9290 0.6940 0.1250;   % RE
        0.4940 0.1840 0.5560];  % RH

figure('Position', [331,285, 901,189]);

CHs = [1 1 1 2 2 3 3 4];
for ch = [1 4 6 8]
    subplot(1, 4, CHs(ch)); 
    for c = [101 102 103 104]
        plot([c c], [c c], 'Color', Color(c-100, :), 'LineWidth', 1);
        hold on
    end
    for c = [101 102 103 104]
        % baseline = mean(X(ylab==c, 1:1000, i, :), "all") * ones(1, length(time));
        x_m = squeeze(mean(BB(ch, :, Ylab==c), 3));
        x_v = squeeze(std(BB(ch, :, Ylab==c), [], 3)) / sqrt(sum(Ylab==101));
        plot(tt, x_m, 'Color', Color(c-100, :), 'LineWidth', 1);
        hold on
        fill([tt fliplr(tt)], [x_m-x_v fliplr(x_m+x_v)], Color(c-100, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on
    end

    axis([-1 4 -7 10]);
    xline(0, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    yline(0);
    xline(2.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    xlabel('time (s)', 'FontSize', 10);
    ylabel('MRCP (\muV)', 'FontSize', 10);
    title(['ch' num2str(ch)], 'FontSize', 10);
    box off
end
% 
% legend1 = legend('LE', 'LH', 'RE', 'RH');
% set(legend1,...
%     'Position',[0.41335883733737 0.951612903225807 0.188252621928521 0.0346667533818946],...
%     'NumColumns',4,...
%     'EdgeColor','none');