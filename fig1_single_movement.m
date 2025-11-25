%% single movement Neural response
clear
clc

load('data.mat');

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



%% Topomap
figure1 = figure('Colormap', slanCM('coolwarm'), 'Position', [573,311.6,120,446]);

axes1 = axes('Parent', figure1);
axis off

box(axes1,'on');
axis(axes1,'ij');

set(axes1, 'CLim',[-0.075 0.075], 'Layer','top');

% colorbar
colorbar(axes1,'Position',[0.4 0.2 0.2 0.6],...
    'Ticks',[0 0.025 0.050 0.075],...,
    'TickLength', 0, ...
    'AxisLocation','in',...
    'FontSize', 10, ...
    'Box', 'off', ...
    'Limits',[0 0.075]);

% textbox
annotation(figure1,'textbox',...
    [0.238674033149171 0.0338164251207729 0.584530386740332 0.143301127214171],...
    'String',{'RMS','amplitude','(a.u.)'},...
    'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','off',...
    'EdgeColor','none');





