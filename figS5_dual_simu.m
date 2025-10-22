%% 生成基础数据
clear
clc

load('data.mat');
% fb = fb(1:51);
% Spectra = Spectra(1:51, :, :);

option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;
option.fmax = 150;
option.maxnff = 512;
option.just_recat = 1;
option.tps = [-1000 200];
option.sfb = 50;


baseline = mean(sqrt(Spectra(fb>option.sfb, :, Ylab==100)), 3);
option.BL = 0; %baseline;
Fc = @(X, baseline) sqrt(X(fb>option.sfb, :)) - baseline;
option.Fc = Fc;

fileuse = 8:53;


[data_A, data_B, spt_A, spt_B] = Mixture_dataset(fileuse, [101 102], option);
[data_C, data_D, spt_C, spt_D] = Mixture_dataset(fileuse, [103 104], option);



%% 生成代理数据回归结果
[Param{1}, Rsquared{1}, modelP{1}] = check_surro_fit(data_A, data_B, spt_A, spt_B, fb, option);
[Param{2}, Rsquared{2}, modelP{2}] = check_surro_fit(data_A, data_C, spt_A, spt_C, fb, option);
[Param{3}, Rsquared{3}, modelP{3}] = check_surro_fit(data_B, data_D, spt_B, spt_D, fb, option);
[Param{4}, Rsquared{4}, modelP{4}] = check_surro_fit(data_C, data_D, spt_C, spt_D, fb, option);



%% 散点图
Tbox = {'LE plus LH', 'LE plus RE', 'LH plus RH', 'RE plus RH'}; % 
A1 = {'LE', 'LE', 'LH', 'RE';
      'LH', 'RE', 'RH', 'RH'};

for i = 1:4
    figure('Position', [573,404,355,353])
    ax1 = axes("Position", [0.15 0.15 0.5 0.5]);
    ccD = Rsquared{i};
    ccD = min(max(ccD, 0.8), 1);
    ccD(end-1) = 0.8;ccD(end) = 1;
    scatter(ax1, Param{i}(:, 1), Param{i}(:, 2), 4, ccD, 'filled');
    hold(ax1, "on");
    colormap([227 240 250; 195 224 244; 186 209 236; 164 191 223; 137 170 215]/255);
    scatter(ax1, Param{i}(1, 1), Param{i}(1, 2), 20, 'r', 'd', 'LineWidth', 1.2);
    grid on
    
    plot(ax1, [-3 3], [-3 3], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    text(2.2, 2.6, 'n.s.', 'HorizontalAlignment','center', "Interpreter", "latex");
    xlabel(ax1, ['$W_{' A1{1,i} '}$'], "Interpreter", "latex", 'FontSize', 11);
    ylabel(ax1, ['$W_{' A1{2,i} '}$'], "Interpreter", "latex", 'FontSize', 11);
    xline(ax1, 0);
    yline(ax1, 0);
    axis(ax1, [-1.5 3 -1.5 3]);
    title(Tbox{i}, 'FontName', 'Arial', 'FontSize', 10);

    % 统计图
    ax2 = axes("Position", [0.29 0.29 0.5*sqrt(2) 0.5*sqrt(2)]);
    resp = (Param{i}(:, 1) - Param{i}(:, 2))/sqrt(2);
    histogram(ax2, resp, 20, 'FaceAlpha', 0.25, 'FaceColor', [137 170 215]/255,...
        'EdgeColor','none', 'Normalization','probability');
    hold on;
    [kde, xi] = ksdensity(resp, 'Bandwidth', 0.1); % 可调整带宽
    plot(ax2, xi, kde/4, 'LineWidth', 0.8, 'Color', [164 191 223]/255);
    yline(0);
    axis(ax2, [-2.25 2.25 0 1.2]);
    axis(ax2, 'off')
    view(ax2, 45, 90);
end



%% 柱状图
ind = [1 4 7 10
    2 5 8 11];
figure('Position', [573,557.7,525,200]);
axes('Position', [0.05 0.13 0.92 0.85]);
for i = 1:4
    
    XX = Param{i}(:, 1);
    YY = Param{i}(:, 2);
    RR = Rsquared{i}(:);
    scatter(0.4*rand(1, length(XX))+ind(1, i), XX, 4, RR, ...
        'filled', 'AlphaData', 0.9*ones(1 ,length(XX)), 'MarkerFaceAlpha', 'flat');
    hold on
    scatter(0.4*rand(1, length(XX))+ind(2, i), YY, 4, RR, ...
        'filled', 'AlphaData', 0.9*ones(1, length(XX)), 'MarkerFaceAlpha', 'flat');
    plot([ind(1, i)-0.03 ind(1, i)+0.43], [mean(XX) mean(XX)], 'Color', [0.314	0.0145	0.0235], 'LineWidth', 2);
    plot([ind(2, i)-0.03 ind(2, i)+0.43], [mean(YY) mean(YY)], 'Color', [0.314	0.0145	0.0235], 'LineWidth', 2);
    plot([ind(1, i)+0.2 ind(2, i)+0.2], [3.8 3.8], 'k', 'LineWidth', 1);
    p = signrank(XX, YY)
    % [~, p] = ttest2(XX, YY)
    text((ind(1, i)+ind(2, i))/2+0.2, 4.15, 'n.s.', 'HorizontalAlignment', 'center');
end
colormap(slanCM('bilbao'));
axis([0.5 11.9 -2 4.3]);
yticks([-2 -1 0 1 2 3])
xticks(ind(:)+0.2);
xticklabels({'$W_{LE}$', '$W_{LH}$', '$W_{LE}$', '$W_{RE}$', '$W_{LH}$', '$W_{RH}$', '$W_{RE}$', '$W_{RH}$'})
set(gca, "TickLabelInterpreter", "latex");



%% Make surroguate dataset
function [data_A, data_B, spt_A, spt_B] = Mixture_dataset(fileuse, Yuse, option)

load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\readme.mat');
data_root = 'E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\';

bsize = length(option.fs*option.tmin + 1: option.fs*option.tmax);

h1 = 1;
h2 = 1;

for id = 1:length(fileuse)
    neo_data = load([data_root Label(fileuse(id)).name]);
    data = NEO_reref(neo_data.data/1e-6, 'average');
    [~, events, ~] = Step0_psd_dataset(Label, data_root, fileuse(id), option);

    for i = 1:size(events, 1)
        if events(i, 3) == Yuse(1)
            data_A{h1} = data(:,  (events(i, 1)+ 1):(events(i, 1)+ bsize));
            spt_A{h1} = neo_calc_spectra(data_A{h1}, [0 0 Yuse(1)], option);
            h1 = h1 + 1;
        end
        if events(i, 3) == Yuse(2)
            data_B{h2} = data(:,  (events(i, 1)+ 1):(events(i, 1)+ bsize));
            spt_B{h2} = neo_calc_spectra(data_B{h2}, [0 0 Yuse(2)], option);
            h2 = h2 + 1;
        end
    end
    id
end


end



%%
function [Param, Rsquared, modelP] = check_surro_fit(data_C, data_D, spt_C, spt_D, fb, option)
Fc = option.Fc;
num = min(length(spt_C), length(spt_D));
num = min(num, 320);
X = zeros(num, sum(fb>option.sfb)*8);
x_a = zeros(num, sum(fb>option.sfb), 8);
x_b = zeros(num, sum(fb>option.sfb), 8);

Ind = randperm(num);
for i = 1:num
    Q = 1*data_C{Ind(i)} + 1*data_D{i};
    Q = NEO_reref(Q, 'average');
    % Q = highpass(Q', 50, 1000)';
    [SS, fb] = neo_calc_spectra(Q, [0 0 0], option);
    SS = Fc(SS, sqrt(2)*option.BL);
    % SS = bsxfun(@minus, SS, mean(SS) - re_aline);
    X(i, :) = SS(:)';

    x_a(i, :, :) = Fc(spt_C{i}, option.BL);
    x_b(i, :, :) = Fc(spt_D{i}, option.BL);
end

X = [mean(X); X];
x_a = mean(x_a(:, :))';
x_b = mean(x_b(:, :))';


Param = []; Rsquared = []; modelP = [];
for k = 1:size(X, 1)
    mdl = fitlm([x_a, x_b], X(k, :)');
    Param = [Param; [mdl.Coefficients{2,1} mdl.Coefficients{3,1}]];
    Rsquared = [Rsquared mdl.Rsquared.Ordinary]; % 获取拟合模型的 R-squared 值
    modelP = [modelP mdl.ModelFitVsNullModel.Pvalue];
end

end




%% A and B
% num = min(length(spt_A), length(spt_B));
%
%
% Param = []; Rsquared = []; modelP = [];
% for u = 0.5:0.05:1.1
%     X = zeros(num, sum(fb>sfb)*8);
%     x_a = zeros(num, sum(fb>sfb), 8);
%     x_b = zeros(num, sum(fb>sfb), 8);
%     for i = 1:num
%         Q = u*data_A{i} + u*data_B{i};
%         Q = NEO_reref(Q, 'average');
%         % Q = highpass(Q', 50, 1000)';
%         [SS, fb] = neo_calc_spectra(Q, [0 0 0], option);
%         SS = sqrt(SS(fb>sfb, :)) - baseline;
%         % SS = bsxfun(@minus, SS, mean(SS) - re_aline);
%         X(i, :) = SS(:)';
%
%         x_a(i, :, :) = sqrt(spt_A{i}(fb>sfb, :));
%         x_b(i, :, :) = sqrt(spt_B{i}(fb>sfb, :));
%     end
%
%     X = [mean(X); X];
%     x_a = mean(x_a(:, :))' - baseline(:);
%     x_b = mean(x_b(:, :))' - baseline(:);
%
%     mdl = fitlm([x_a, x_b], X(1, :)', 'Intercept', false);
%     Param = [Param; [mdl.Coefficients{1,1} mdl.Coefficients{2,1}]];
%     Rsquared = [Rsquared mdl.Rsquared.Ordinary]; % 获取拟合模型的 R-squared 值
%     modelP = [modelP mdl.ModelFitVsNullModel.Pvalue];
%     u
% end
%
% %%
% scatter(Param(:, 1), Param(:, 2), 8, Rsquared, 'filled');
% hold on
% colormap(slanCM('Greys'));
% grid on
% cbar = colorbar(); cbar.Limits = [0.5 1];
% plot([-3 3], [-3 3], '--', 'Color', [0.5 0.5 0.5]);
% xlabel('$W_{LE}$', "Interpreter", "latex");
% ylabel('$W_{LH}$', "Interpreter", "latex");
% xline(0)
% yline(0)
% axis([-1.5 2.5 -2 3])
% axis equal