%% continue_HG
clear;
clc

% tmin = -1; 
% tmax = 4;
% option.fs = 1000;
% option.fb = 55:5:95;
% option.tmin = tmin;
% option.tmax = tmax;

continueHG = load('data\continueHG.mat');
Scalingp = continueHG.Scalingp;
Ylab = continueHG.Ylab;
time = continueHG.time;

Preprocess_funy = @(Ylab) Ylab;
Preprocess_funx = @(Scalingp) bsxfun(@minus, sqrt(Scalingp), mean(sqrt(Scalingp(:, 1:1000, :, :)), 2));


%%  Trial-averaged
ylab = Preprocess_funy(Ylab);
X = mean(Preprocess_funx(Scalingp), 4);

leftMargin = 0.05;
rightMargin = 0.03;
topMargin = 0.09;
bottomMargin = 0.18;
horizontalGap = 0.065;
verticalGap = 0.05;
axWidth = (1 - leftMargin - rightMargin - 3*horizontalGap)/4;
axHeight = 1 - topMargin - bottomMargin;

CMAP = slanCM('tab20c');
Color = [CMAP(13*6+1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(1, :)];

figure('Position', [331,460.3,966.6,205.9]);
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

