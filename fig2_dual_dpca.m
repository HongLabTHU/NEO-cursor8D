%%
clear
clc

option.fs = 1000;
option.tmin = 0;
option.tmax = 1;
option.fpoint = 201;  % 201;
option.fmax = 150;  % 100;
option.maxnff = 512;
option.tps = [-1000 200];

fileuse = 8:53;
load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\readme.mat')
data_root = 'E:\MATLAB_softhub\SPMdataset\NEO_TT01\XB_data\';

[Spectra, Ylab, fb] = Step0_psd_dataset(Label, data_root, fileuse, option);



%%
clear;clc
load('data.mat');

freq_func = @(X, fb) X(fb>50, :, :);
Preprocess_funy = @(Ylab) Ylab;
Preprocess_funx = @(Spectra) bsxfun(@minus, sqrt(Spectra), mean(sqrt(Spectra(:, :, Ylab==100)), 3));
% Preprocess_funy = @(Ylab) Ylab(2:2:end);
% Preprocess_funx = @(Spectra) sqrt(Spectra(:, :, 2:2:end)) - sqrt(Spectra(:, :, 1:2:end));


%% dPCA of single mode
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);

X = permute(X, [3 1 2]);
X = X(:, :);
[~, muX, stdX] = zscore(X);


condu = [101 102 103 104];
X = X(ismember(ylab, condu), :);
X = bsxfun(@rdivide, bsxfun(@minus, X, muX), stdX);
ylab = ylab(ismember(ylab, condu));
ylab_s = ylab;

% addpath 'E:\MATLAB_softhub\SPMdataset\NEO_BCI\dPCA'
[W_single, V_single, explVar_single] = NEO_single_DPCA(X', ylab, condu);

% [W_single, V_single, whM_single, explVar_single] = NEO_reduced_dPCA(X, ylab, 5);
% V_single = reshape(V_single, [1 size(V_single)]);
% V_single = cat(1, V_single,V_single);
X_mar = [mean(X(ylab==101, :)); mean(X(ylab==102, :)); mean(X(ylab==103, :)); mean(X(ylab==104, :))];
% [V_single,~,~] = eigs(X_mar'*X_mar, 24);
% V_single = - V_single;


%% Bootstrap resampling
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);

X = permute(X, [3 1 2]);
X = X(:, :);

condu = [101 102 103 104];
X_ = X(ismember(ylab, condu), :);
ylab = ylab(ismember(ylab, condu));

for N = 1:1000
    index = randsrc(length(ylab), 1, 1:length(ylab));
    index = sort(index);
    X = X_(index, :);
    ylab_s = ylab(index);
    [~, V_single(N, :, :), ~, ~] = NEO_reduced_dPCA(X, ylab_s, 5);
    if abs(max(squeeze(V_single(N, :, 2)))) > abs(min(squeeze(V_single(i, :, 2))))
        V_single(N, :, 2) = -V_single(N, :, 2);
    end
    N
end

%% View single PC
figure('Position', [573,378.3,349,365])
sgs = [1 3 5 7 2 4 6 8];
for i = 1:8
    subplot(4, 2, sgs(i));
    hold on

    PC1 = squeeze(V_single(:, 51*i-50:51*i, 1));
    plot(fb(fb>50), mean(PC1), 'Color', [155 109 96]/255, 'LineWidth', 1.2);
    fill([fb(fb>50)' fliplr(fb(fb>50)')], ...
        [prctile(PC1, 0.5, 1) fliplr(prctile(PC1, 99.5, 1))], ...
        [155 109 96]/255, 'linestyle', 'none', 'FaceAlpha', 0.3);

    PC2 = squeeze(V_single(:, 51*i-50:51*i, 2));
    plot(fb(fb>50), mean(PC2), 'Color', [203 180 172]/255,'LineWidth', 1.2);
    fill([fb(fb>50)' fliplr(fb(fb>50)')], ...
        [prctile(PC2, 0.5, 1) fliplr(prctile(PC2, 99.5, 1))], ...
        [203 180 172]/255, 'linestyle', 'none', 'FaceAlpha', 0.3);

    title(['ch' num2str(i)]);
    yline(0, '--');
    axis([50 150 -0.25 0.25])
    box off
    if i < 5
        yticks([-0.25 0 0.25]);
    else
        yticks([]);
    end
    if sgs(i) > 6
        xticks([50 100 150]);
        xlabel('freq (Hz)')
    else
        xticks([]);
    end
end


%% Project single mode on V_single
CMAP = slanCM('tab20c');
Cmp = [CMAP(1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(13*6+1, :)];
Cmp = [CMAP(13*6+1, :); CMAP(13*4+1, :); CMAP(13*2+1, :); CMAP(1, :)];
Y_s = X * V_single(:, 1:2);

% [Y_s, mean_Ys, std_Ys] = zscore(Y_s, [], 1);

figure('Position', [573,503.7,310,254]);
ax = axes('Position', [0.15 0.15 0.8 0.8]);

for i = 1:4
    scatter(ax, Y_s(ylab==condu(i), 1), Y_s(ylab==condu(i), 2), 10, Cmp(i, :), 'filled', ...
        'AlphaData', 0.9*ones(1, sum(ylab==condu(i))), 'MarkerFaceAlpha', 'flat');
    %     error_ellipse(Y_s(ylab==condu(i), 1), Y_s(ylab==condu(i), 2), ax, Cmp(i, :));
    hold(ax, 'on');
end


legend('LE', 'LH', 'RE', 'RH', ...
    'Position',[0.722424051474638 0.167058674117348 0.220430107526882 0.242125984251968],...
    'FontName','Arial',...
    'FontSize', 8,...
    'EdgeColor','none',...
    'Color','none');

% xticks([-1.2 0 1.5]); yticks([-0.8 0 0.8])
% axis([-1.2 1.5 -0.8 0.8]);
xticks([-18 0 30]); yticks([-10 0 10])
axis([-18 30 -10 10]);
xlabel(ax, 'Single PC1', 'FontName', 'Arial');
ylabel(ax, 'Single PC2', 'FontName', 'Arial');



%% Project dual mode on V_single
% Cmps = [CMAP(13*8+1, :); CMAP(13*1+1, :); CMAP(13*5+1, :); CMAP(13*10+1, :)];
Cmps = repmat([225   162   164]/255, 4, 1);
Cmps = repmat([0.65 0.65 0.65], 4, 1);
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);
X = permute(X, [3 1 2]);
X = X(:, :);

condu = [121 122 125 126];
X = X(ismember(ylab, condu), :);
X = bsxfun(@rdivide, bsxfun(@minus, X, muX), stdX);
ylab = ylab(ismember(ylab, condu));


Tbox = {'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
fPos = {[0.06, 0.1, 0.225, 0.7], [0.29, 0.1, 0.225, 0.7], [0.52, 0.1, 0.225, 0.7], [0.75, 0.1, 0.225, 0.7]};

figure('Position', [447,409,900,240], 'Color', [1 1 1]);
for i = 1:4
    ax(i) = axes('Position', fPos{i});
    Y{i} = squeeze(X(ylab==condu(i), :) * V_single(:, 1:2));
    % Y{i} = bsxfun(@rdivide, bsxfun(@minus, Y{i}, mean_Ys), std_Ys);
    scatter(ax(i), Y{i}(:, 1), Y{i}(:, 2), 8, Cmps(i, :), 'filled', 'MarkerFaceAlpha', 0.85, 'LineWidth', 20);

    hold(ax(i), 'on');
    for j = [101 102 103 104]
        error_ellipse(Y_s(ylab_s==j, 1), Y_s(ylab_s==j, 2), ax(i), Cmp(j-100, :), 0.5);
    end
    
    axis(ax(i), [-12 30 -7 5]);
    xticks([-15 0 35]); yticks([-8 0 6]);
    % axis(ax(i), [-1 1.2 -0.5 0.5]);
    % xticks([-1 0 1.2]);
    axis off
    % xlabel(ax(i), 'Single PC1', 'FontName', 'Arial', 'FontSize', 11);
    % ylabel(ax(i), 'Single PC2', 'FontName', 'Arial', 'FontSize', 11);
    title(ax(i), Tbox{i}, 'FontName', 'Arial', 'FontSize', 11); %, 'Color', Cmps(i, :), 'FontWeight','bold');
end



%% dPCA of single mode
ylab = Ylab;
X = permute(log10(Spectra), [3 1 2]);
X = X(:, :);

condu = [121 122 123 124 125 126];
X = X(ismember(Ylab, condu), :);
ylab = ylab(ismember(Ylab, condu));
X = bsxfun(@minus, X, mean(X));

[W_dual, V_dual, whM_dual, explVar_dual] = NEO_reduced_dPCA(X, ylab);

X = bsxfun(@minus, X, mean(X));
for i = 1:6
    X_ = X(ismember(ylab, condu(i)), :);
    for j  = 1:5
        Svar(i, j) = trace(X_ * V_single(:, 1:j) * V_single(:, 1:j)' * X_');
        Svar(i, j) = Svar(i, j) / trace(X_ * X_');
        Dvar(i, j) = trace(X_ * V_dual(:, 1:j) * V_dual(:, 1:j)' * X_');
        Dvar(i, j) = Dvar(i, j) / trace(X_ * X_');
    end

end

Tbox = {'LE-LH', 'LE-RE', 'LE-RH', 'LH-RE', 'LH-RH', 'RE-RH'};
for i = 1:6
    subplot(2, 3, i);
    plot([Svar(i, :); Dvar(i, :)]');
    ylim([0 0.4]);
    title(Tbox{i});
    hold on
end



%% dPCA of All mode 
ylab = Preprocess_funy(Ylab);
X = Preprocess_funx(Spectra);
X = freq_func(X, fb);

X = permute(X, [3 1 2]);
X = X(:, :);

condu = [101 102 103 104 121 122 125 126];
X = X(ismember(ylab, condu), :);
ylab = ylab(ismember(ylab, condu));
ylab_s = ylab;

addpath 'E:\MATLAB_softhub\SPMdataset\NEO_BCI\dPCA'
[W_dual, V_dual, whM_dual, explVar_dual] = NEO_reduced_dPCA(X, ylab, 10);


Cmp = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.2118 0.8392 1.00;
    0.6350 0.0780 0.1840;
    0.8353 0.2196 1.00];

Y_s = bsxfun(@minus, X, mean(X)) * V_dual(:, 1:2);

figure('Position', [573,503.7,310,254]);
ax = axes('Position', [0.15 0.15 0.8 0.8]);
hold(ax, 'on');
for i = 1:8
    scatter(ax, -mean(Y_s(ylab==condu(i), 1)), mean(Y_s(ylab==condu(i), 2)), 10, Cmp(i, :), 'filled', ...
        'AlphaData', 0.9, 'MarkerFaceAlpha', 'flat');
    error_ellipse(-Y_s(ylab==condu(i), 1), Y_s(ylab==condu(i), 2), ax, Cmp(i, :), 0.5);
end


