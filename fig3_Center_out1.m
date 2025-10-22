%%
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Center');

useid = 13:75;


TimeMat = zeros(24, length(useid));

for i = 1:length(useid)
    try
        load(Flist{useid(i)});
        Control.Tarlist = Control.Tarlist;
        Behave = Control.Behave(:);
        for h = 1:size(Behave, 1)
            TimeMat(h, i) = Behave(h).timecost < 8;
        end
    catch
    end
   
end
        
clc
for i = 1:length(useid)
    Flist{useid(i), 2} = sum(TimeMat(:, i))/24;
end


%% Center-out Task performance
clear;clc
rootDir = 'E:\pycharm\NEOdata\TT01';
Flist = findMatFilesByName(rootDir, 'Center');
useid = 15:19;  % 11月的可用数据
% useid = [37 38 39];     % Sing. Move. 4D
% useid = [31 35 36];     % Sing. & Dual. 4D

N = length(useid);

T(24*N) = 0;
Pick(24*N).Target = 0;
Pick(24*N).Targetpos = 0;
Pick(24*N).hit = 0;
Pick(24*N).time = 0;
Pick(24*N).V = [0 0];
Pick(24*N).X = [0 0];

for i = 1:N
    load(Flist{useid(i)});
    Control.Tarlist = Control.Tarlist;
    Behave = Control.Behave(:);
    for h = 1:size(Behave, 1)
        Pick(24*i-24+h).Target = Control.Tarlist(h);
        Pick(24*i-24+h).Targetpos = 0.4*[cos(double(Control.Tarlist(h)-1)*pi/4) sin(double(Control.Tarlist(h)-1)*pi/4)];
        Pick(24*i-24+h).hit = Behave(h).timecost < 8;
        Pick(24*i-24+h).time = Behave(h).timecost;
        Pick(24*i-24+h).X = Behave(h).Path;
        T(24*i-24+h) = Behave(h).timecost;
    end
end

clc
%%
Cmaps = [0.7305    0.8320    0.5898;
         0.6196    0.7922    0.8824;
         0.0902    0.7451    0.8118;
         0.1922    0.5098    0.7412;
         0.4660    0.6740    0.1880;
         0.9020    0.3333    0.0510;
         0.8398    0.1680    0.1641;
         0.9922    0.6824    0.4196];

kk = [6 3 2 1 4 7 8 9];

figure('Position', [773.7,255.6667,270.6667,264.6667]);
axes('Position', [0 0 1 1])

r = 0.6;
for i = 1:8
    hold on
    C = r*[cos((i-1)*pi/4) sin((i-1)*pi/4)];
    plot((r+0.016)*[cos((i-1)*pi/4+0.2) cos((i-1)*pi/4-0.2)], ...
        (r+0.016)*[sin((i-1)*pi/4+0.2) sin((i-1)*pi/4-0.2)], 'black', 'LineWidth', 1);
    target = 0.4*[cos((i-1)*pi/4) sin((i-1)*pi/4)];
    rectangle('Position',[target(1)+C(1)-0.11, target(2)+C(2)-0.11, 0.11*2, 0.11*2], 'LineWidth', 1.5, ...
        'Curvature',[1 1], 'EdgeColor', [0.00,0.45,0.74], 'FaceColor','none');%[0.30,0.75,0.93]);
    % scatter(target(1)+C(1), target(2)+C(2), 250, [0.2000 0.5569 0.7922], 'filled', 'MarkerFaceAlpha', 1);
    for h = 1:72
        if (Pick(h).Target == i)
            if Pick(h).hit
                plot(Pick(h).X(:, 1)+C(1), Pick(h).X(:, 2)+C(2), 'Color', [.7 .7 .7], 'LineWidth', 1);
            else
                % plot(Pick(h).X(:, 1)+C(1), Pick(h).X(:, 2)+C(2), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.3);
            end

        end
    end

    scatter(C(1), C(2), 12, 'black', 'filled');
    axis off;
end
axis([-0.52-r 0.52+r -0.52-r 0.52+r]);
xticks([]); yticks([]);
axis equal


%% 计算时间累计正确曲线

data = zeros(1, 8);
for t = 1:8
    data(t) = sum(logical((T < t).*(T>t-1))) / length(T);
end

% qing = [0.1250 0.5781 0.8047];
qing = [0.65 0.65 0.65];
hong = [0.8500 0.3250 0.0980];

colororder([qing; hong]);

set(gcf, 'Position', [573,481.7, 264,258.7]); %256,276]);
% 左侧坐标轴（柱状图）
yyaxis left
bar(1:8, data, 'FaceColor', qing, 'EdgeColor', 'none');
ylabel('Hit trial histogram', 'FontSize', 11);
ylim([0 0.4]);
ax = gca;
ax.YAxis(1).LineWidth = 1.5;
ax.FontWeight = 'normal';

% 右侧坐标轴（累计折线图）
yyaxis right
plot(1:8, cumsum(data), '-', 'LineWidth', 1.3, 'Color', hong);
hold on
plot(1:8, cumsum(data), '.', 'MarkerSize', 15, 'Color', hong);
ylabel('Hit rate', 'FontSize', 11);
ylim([0 1.0]);
ax = gca;
ax.YAxis(2).LineWidth = 1.5;
ax.FontWeight = 'normal';

box off
xlabel('time (s)', 'FontSize', 11);


%% 计算 ITR, trajectory distance ratio，trajectory error，trajectory variability
clc
fprintf('mean hit time %.2f+%.2f\n', mean(T(T<8)), std(T(T<8)))

T_ = T;
T_(T_ > 8) = NaN; %8;
p = cumsum(data);
p = p(end);
ITR = log2(8) + p * log2(p) + (1 - p) * log2((1-p)/(8-1));
ITR = ITR / mean(T_, 'omitnan') * 60;
fprintf('ITR %.2f bpm\n ', ITR);
ITR_per_ch = ITR / 8;
fprintf('ITR_per_ch %.2f bpm\n ', ITR_per_ch);

Fitts_ITR = log2((0.4+0.22)/0.22)./T_;
Fitts_ITR = mean(Fitts_ITR, 'omitnan') * 60;
fprintf('Fitts_ITR %.2f bpm\n ', Fitts_ITR);

for h = 1:72
    if Pick(h).time < 10
        i = Pick(h).Target;
        rot_mat = [cos((i-1)*pi/4) -sin((i-1)*pi/4); sin((i-1)*pi/4) cos((i-1)*pi/4)];
        X = Pick(h).X * rot_mat;
        Road(h) = sum(abs(diff(X(:,1) + 1j*X(:, 2))));
        TE(h) = mean(abs(X(:, 2)));
        TV(h) = std(abs(X(:, 2)));
    else
        Road(h) = nan;
        TE(h) = nan;
        TV(h) = nan;
    end
end


fprintf('TDR  %.2f +- %.2f\n', nanmean(Road) / 0.29, nanstd(Road) / 0.29)
fprintf('TE  %.2f +- %.2f\n', nanmean(TE) / 0.29, nanstd(TE) / 0.29)
fprintf('TV  %.2f +- %.2f\n', nanmean(TV) / 0.29, nanstd(TV) / 0.29)



%%
function matFiles = findMatFilesByName(rootDir, searchStr)
    % 查找文件名包含 searchStr 的 .mat 文件
    files = dir(fullfile(rootDir, ['**/*' searchStr '*.mat']));  % 例如 *data*.mat
    matFiles = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);
    matFiles = matFiles(:);
end