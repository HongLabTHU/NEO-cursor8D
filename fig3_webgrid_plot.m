%% Center-out Task performance
clear; clc
rootDir = 'behave\webgrid';
Flist = findMatFilesByName(rootDir, 'Web');
useid = 1:5;


%% Fitts ITR

k = 1;
for i = 1:length(useid)
    load(Flist{i});
    Control.Tarlist = Control.Tarlist;
    Behave = Control.Behave(1:24);
    for h = 1:length(Behave)
        if Behave(h).timecost < 14
            tarpos = [Control.Mesh1(Control.Tarlist(h)) Control.Mesh2(Control.Tarlist(h))];
            D(k) = norm(tarpos - Control.Behave(h).Path(1, :));
            FittsMat(k) = 60 * log2((D(k) + 0.16)/0.16) / Behave(h).timecost;
        else
            FittsMat(k) = NaN;
            D(k) = NaN;
        end
        k = k + 1;
    end
end

figure('Position', [573,507,394,250.6])
axes('Position', [0.14 0.15 0.65 0.8])
scatter(D, log10(FittsMat), [], [0.5529 0.6275 0.7961], 'filled', 'MarkerFaceAlpha', 0.5);
ylim([0.5 3]);
yticks([0.5 1 2 3]);
yticklabels({'', '1e1', '1e2', '1e3'});
box off
xlabel('​​Start-Target dist.​', 'FontName', 'Arial');
ylabel('Fitts ITR (bpm)', 'FontName', 'Arial');

axes('Position', [0.83 0.15 0.13 0.8])
h = histogram(log10(FittsMat), 20, ...
    'Orientation', 'horizontal', ...
    'FaceColor', [0.5529 0.6275 0.7961], ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none', ...
    'Normalization', 'pdf'); 

hold on;

[kde, xi] = ksdensity(log10(FittsMat), 'Bandwidth', 0.1); 
plot(kde, xi, 'LineWidth', 1.5, 'Color', [0.5529 0.6275 0.7961]);
axis([0 2.2 0.5 3]);

% mean
yline(mean(log10(FittsMat), 'omitnan'), '--', 'Color', [0.8529 0.6275 0.701], 'LineWidth', 1.5);
axis off



%%
function matFiles = findMatFilesByName(rootDir, searchStr)
    files = dir(fullfile(rootDir, ['**/*' searchStr '*.mat'])); 
    matFiles = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);
    matFiles = matFiles(:);
end


