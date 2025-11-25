%% one channel check
clear;clc

load('data\PSDdata.mat')
freq_func = @(X, fb) X(logical((fb>50).*(fb<150)), :, :);
Preprocess_funy = @(Ylab) Ylab;
Preprocess_funx = @(Spectra) bsxfun(@minus, sqrt(Spectra), mean(sqrt(Spectra(:, :, Ylab==100)), 3));



%% linear comb
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);

X = permute(X, [3 1 2]);        % 注意这里

condu = [101 102 103 104 121 122 125 126];
X = X(ismember(ylab, condu), :, :);
ylab = ylab(ismember(ylab, condu));


base = {[101 102], [101 103], [102 104], [103 104]};
move = {121, 122, 125, 126};

Param = cell(length(move), 8);
modelP = cell(length(move), 8);
Rsquared = cell(length(move), 8);

m = 1;
% different basis
for i = 1:length(base)
    for ch = 1:8
        x_a = mean(X(ismember(ylab, base{i}(1)), :, ch))';
        x_b = mean(X(ismember(ylab, base{i}(2)), :, ch))';
        x_p = squeeze(X(ismember(ylab, move{i}), :, ch));
        x_p = [mean(x_p); x_p];
        subplot(4, 8, m);
        plot([x_a x_b x_p(1,:)']);

        for k = 1:size(x_p, 1)
            mdl = fitlm([x_a, x_b], x_p(k, :)'); %, 'Intercept', false);
            Param{i, ch} = [Param{i, ch}; [mdl.Coefficients{2,1} mdl.Coefficients{3,1}]];
            Rsquared{i, ch} = [Rsquared{i, ch} mdl.Rsquared.Ordinary]; % 获取拟合模型的 R-squared 值
            modelP{i, ch} = [modelP{i, ch} mdl.ModelFitVsNullModel.Pvalue];
        end
        m = m +1;
    end
    i
end



%%
k = 1;

for i = 1:length(base)
    for ch = 1:8
        subplot(4, 8, k);
        scatter(Param{i, ch}(:, 1), Param{i, ch}(:, 2), 7, Rsquared{i, ch}, 'filled');
        hold on
        colormap("sky");
        scatter(Param{i, ch}(1, 1), Param{i, ch}(1, 2), 20, 'k', 'd');
        plot([-3 3], [-3 3], '--', 'Color', [0.5 0.5 0.5]);
        xlabel('$W_{a}$', "Interpreter", "latex");
        ylabel('$W_{b}$', "Interpreter", "latex");
        xline(0)
        yline(0)
        axis equal
        axis([-1.5 2.5 -2 3]);
        k = k + 1;
    end
end


%%
TTmat = zeros(4, 8);
RR = zeros(4, 8);
UU = zeros(4, 8);
for i = 1:length(base)
    for ch = 1:8
        UU(i, ch) = mean(Param{i, ch}(:, 1)) > mean(Param{i, ch}(:, 2));
        UU(i, ch) = 2*UU(i, ch)-1;
        [~, TTmat(i, ch)] = ttest(Param{i, ch}(:, 1), Param{i, ch}(:, 2));
        RR(i, ch) = mean(Rsquared{i, ch});
    end
end

Tbox = {'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
FigP = [0.1 0.55 0.3 0.35;
        0.1 0.05 0.3 0.35;
        0.5 0.55 0.3 0.35;
        0.5 0.05 0.3 0.35];
Cmp = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560];


fig1 = figure('Position', [573,427,402.7,330.7]);
CC = colormap(slanCM('coolwarm'));
for i = 1:length(base)
    Ax(i) = axes('Position', FigP(i, :));

    Pos = [1.3 0;
        1.3 1;
        1.3 2;
        1.3 3;
        0 0;
        0 1;
        0 2;
        0 3] * 100;
    hold on

    tt_cc = -UU(i, :).*log10(TTmat(i, :));
    tt_cc(tt_cc > 3) = 3;
    tt_cc(tt_cc < -3) = -3;
    tt_cc = CC(floor((tt_cc+3)*255/6)+1, :);

    scatter(Ax(i), Pos(:, 1), Pos(:, 2),  70*RR(i, :), tt_cc, "filled");
    if i == 1
        text(Pos(1:4, 1)+30, Pos(1:4, 2), {'1', '2', '3', '4'});
        text(Pos(5:8, 1)-60, Pos(5:8, 2), {'5', '6', '7', '8'});
    end
    
    title(Tbox{i}, 'Color', Cmp(i, :), 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    xticks([]); yticks([]);
    axis([-50 180 -50 350]);
    axis equal
    box on
    
end


cbar = colorbar('Position',[0.84 0.2 0.03 0.6], ...
    'FontName', 'Arial', ...
    'FontSize', 9);
cbar.Ticks = linspace(cbar.Limits(1), cbar.Limits(2), 7);
cbar.TickLabels = {'.001', '.010', '.100', '1.00', '.100', '.010', '.001'};


annotation(fig1,'textbox',...
    [0.935 0.593698115109363 0.161410479264961 0.0766051809293418],...
    'String',{'p value'},...
    'Rotation',270,...
    'FontName','Arial',...
    'EdgeColor','none');

annotation(fig1,'textbox',...
    [0.83 0.82 0.223491432828408 0.0907166616268519],...
    'String',{'$W_{a}>W_{b}$'},...
    'FitBoxToText','off',...
    'Interpreter', 'latex',...
    'EdgeColor','none', ...
    'Color', [0.717435000000000	0.0511180000000000	0.158737000000000]);

annotation(fig1,'textbox',...
    [0.83 0.0675414776736203 0.223491432828408 0.090716661626852],...
    'String',{'$W_{a}<W_{b}$'},...
    'Interpreter', 'latex',...
    'FitBoxToText','off',...
    'EdgeColor','none', ...
    'Color', [0.238948000000000	0.312365000000000	0.765676000000000]);


