%% Arrow out 与 Center-out 离线仿真

clear;
clc

% 混淆矩阵，信道传输矩阵
C = load('C.mat');



% 种群向量
Pv = [0         0;
   -0.4845    0.8748
   -0.4574   -0.8892
    0.9643    0.2647
    0.8567   -0.5158
   -0.9877   -0.1564
   -0.3308    0.9437
   -0.1445   -0.9895
    0.9893   -0.1458]*2;
sPv = [0         0;
    -sqrt(2)  sqrt(2);
    -sqrt(2)    -sqrt(2);
    sqrt(2)     sqrt(2);
    sqrt(2)     -sqrt(2)];

Text = {'100', '101', '102', '103', '104', '121', '122', '125', '126'};

% rPV = Pv;
% for i = 2:9
%     hold on
%     text(rPV(i, 1), rPV(i, 2), Text{i}, 'Color', 'k');
%     axis([-3 3 -3 3])
% end


%% Arrow-out任务仿真
dphi = pi / 36;

% sing. & Dual. 8D
Comb = perms(3:9);
for j = 1:size(Comb, 1)
    rPV = Pv;
    rPV(3:9, :) = Pv(Comb(j, :), :);
    for tar = 1:72
        Tarpos = [cos((tar-1)*dphi) sin((tar-1)*dphi)];
        m_porb = gen_move_prob(Tarpos, rPV, 0);  %
        V1 = m_porb * C.C_sd8D * rPV;
        rErr(j, tar) = V1 * Tarpos' / 2;
    end
    if mod(j, 1000) == 0
        j
    end
end


% single 4D
for tar = 1:72
    Tarpos = [cos((tar-1)*dphi) sin((tar-1)*dphi)];
    m_porb = gen_move_prob(Tarpos, sPv, 0);  %
    V1 = m_porb * C.C_sg4D * sPv;
    aff(tar) = V1 * Tarpos' / 2;
end
  = [rErr; aff];

% sing. & Dual. 4D
for tar = 1:72
    Tarpos = [cos((tar-1)*dphi) sin((tar-1)*dphi)];
    m_porb = gen_move_prob(Tarpos, sPv, 0)  %
    V1 = m_porb * C.C_sd4D * sPv;
    aff(tar) = V1 * Tarpos' / 2;
end
rErr = [rErr; aff];

%%
figure('Position', [234,484,1292,259]);
subplot(1, 3, 1)
heatmap(Text, Text, C.C_sd8D, 'GridVisible', 'off', 'ColorLimits', [0 0.9]);

subplot(1, 3, 2)
figure('Position', [573,546.3,285.3,211.3])
histogram(mean(rErr, 2), 20, 'EdgeColor', 'none', 'FaceColor', [0.6055 0.8242 0.9336], 'Normalization', 'probability');
xline(mean(rErr(end-2, :)), '--', 'Color', [0, 0.4, 0.7], 'LineWidth', 1);
text(mean(rErr(end-2, :))+0.02, 0.14, 'Sing. & Dual. 8D', 'Color', [0, 0.4, 0.7], 'FontName', 'Arial', 'Rotation',270);

xline(mean(rErr(end, :)), ':', 'Color', [0.2, 0.65, 0.87], 'LineWidth', 1);
text(mean(rErr(end, :))+0.02, 0.14, 'Sing. & Dual. 4D', 'Color', [0.2, 0.65, 0.87], 'FontName', 'Arial', 'Rotation',270);

xline(mean(rErr(end-1, :)), '-.', 'Color', [0.4, 0.8, 1], 'LineWidth', 1);
text(mean(rErr(end-1, :))+0.02, 0.14, 'Sing. Move. 4D', 'Color', [0.4, 0.8, 1], 'FontName', 'Arial', 'Rotation',270);

xlabel('Direction affinity', 'FontName', 'Arial');
yticks([0 0.07 0.14]);
xticks([0.1 0.3 0.5 0.7])
axis([0.1 0.7 0 0.14]);
box off


subplot(1, 3, 3)
[~, ind] = max(mean(rErr, 2));
rPV(3:9, :) = Pv(Comb(ind, :), :);
for i = 2:9
    hold on
    text(rPV(i, 1), rPV(i, 2), Text{i}, 'Color', 'k');
    axis([-3 3 -3 3])
end


%% center-out 任务仿真
Comb = perms(3:9);
TTim = zeros(1, size(Comb, 1));
for j = 1:size(Comb, 1)
    rPV = Pv;
    rPV(3:9, :) = Pv(Comb(j, :), :);
    [tt, LoGer] = co_sim(rPV, C.C_sd8D, 9);
    % itr = log2((0.4+0.11)/0.11)./tt;
    TTim(j) = mean(tt(:));
    if mod(j, 200) == 0
        j
    end
end

% single 4D
[tt, LoGer] = co_sim(sPv, C.C_sg4D, 9);
TTim = [TTim mean(tt(:))];
% sing. & Dual. 4D
[tt, LoGer] = co_sim(sPv, C.C_sd4D, 9);
itr = log2((0.4+0.11)/0.11)./tt;
TTim = [TTim mean(tt(:))];

% Cmaps = [0.8353 0.2196 1.00;
%     0.9290 0.6940 0.1250;
%     0.2118 0.8392 1.00;
%     0 0.4470 0.7410;
%     0.4660 0.6740 0.1880;
%     0.8500 0.3250 0.0980;
%     0.6350 0.0780 0.1840;
%     0.4940 0.1840 0.5560];
% 
% for tar = 1:8
%     for j = 1:40
%         plot(LoGer(tar, j).X(:, 1), LoGer(tar, j).X(:, 2), 'Color', Cmaps(tar, :));
%         hold on
%     end
% end

%%
figure('Position', [573,546.3,285.3,211.3])
histogram(TTim, 20, 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'Normalization', 'probability');
xline(TTim(end-2), '--', 'Color', [0.8500 0.3250 0.098], 'LineWidth', 1);
text(TTim(end-2)+0.15, 0.14, 'Sing. & Dual. 8D', 'Color', [0.8500 0.3250 0.098], 'FontName', 'Arial', 'Rotation',270);

xline(TTim(end), ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 1);
text(TTim(end)+0.15, 0.14, 'Sing. & Dual. 4D', 'Color', [0.65 0.65 0.65], 'FontName', 'Arial', 'Rotation',270);

xline(TTim(end-1), '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(TTim(end-1)+0.15, 0.14, 'Sing. Move. 4D', 'Color', [0.5, 0.5, 0.5], 'FontName', 'Arial', 'Rotation',270);

% histogram(TTim, 20, 'EdgeColor', 'none', 'FaceColor', [0.6055 0.8242 0.9336], 'Normalization', 'probability');
% xline(TTim(end-2), '--', {'Sing. & Dual. 8D'}, 'Color', [0, 0.4, 0.7], 'LineWidth', 1);
% xline(TTim(end), ':', {'Sing. & Dual. 4D'}, 'Color', [0.2, 0.65, 0.87], 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
% xline(TTim(end-1), '-.', {'Sing. Move. 4D'}, 'Color', [0.4, 0.8, 1], 'LineWidth', 1);

xlabel('Hit time (s)', 'FontName', 'Arial');
yticks([0 0.07 0.14]);
axis([4 10 0 0.14]);
box off



%% center-out任务程序
function [TT,LoGer] = co_sim(Pv, M, times)

TT = zeros(8, times);

for tar = 1:8
    Tarpos = 0.4*[cos((tar-1)*pi/4) sin((tar-1)*pi/4)];
    for j = 1:times
        Loger.X = [0 0];
        Loger.T = 0;
        while (norm(Loger.X(end, :) - Tarpos) >= 0.11) && (Loger.T<10)
            P_tar = Tarpos - Loger.X(end, :);   % 指向目标的向量
            m_porb = gen_move_prob(P_tar, Pv, 0.1);  % 动作先验概率
            m_porb = gen_post_porb(m_porb, M);  % 动作后验概率
            V1 = m_porb * Pv;
            V2 = nfunc(P_tar, m_porb) .* V1;
            Loger.X = [Loger.X; Loger.X(end, :) + 0.025 * V2];
            Loger.T = Loger.T + 0.2;
        end
        LoGer(tar, j) = Loger;
        TT(tar, j) = Loger.T;
    end
end

end



%% 依据指向目标的向量生成信源概率
function prob = gen_move_prob(P_tar, PV, prest)

% 计算种群向量中的每一个到P_tar的夹角
theta = zeros(1, size(PV, 1));
prob = zeros(1, size(PV, 1));
for i = 2:size(PV, 1)
    in_pord = (PV(i, :)*P_tar') / (norm(PV(i, :)) * norm(P_tar));
    in_pord = min(in_pord, 1);
    in_pord = max(in_pord, -1);
    theta(i) = acos(in_pord);
end

theta(1) = pi;
[s_theta, ind] = sort(theta);

dphi = angle(PV(ind(1), 1)+1j*PV(ind(1), 2)) - angle(PV(ind(2), 1)+1j*PV(ind(2), 2));
dphi = abs(angle(exp(1j*dphi)));

prob(ind(1)) = s_theta(2) / dphi;
prob(ind(2)) = s_theta(1) / dphi;
prob(1) = prest;
% prob(2:end) = cos(theta(2:end)/2).^2;
prob = prob ./ sum(prob);
end


%% 依据属于某类的概率分布随机抽样一个分类后验概率
function prob = gen_post_porb(m_porb, M)
ind = randsrc(1, 1, [1:size(M, 1); m_porb]);
V = 0 * m_porb;
V(ind) = 1;
prob = V * M;
end


%% 依据到目标的误差向量生成模长
function u = nfunc(P_tar, m_porb)
x = 5*norm(P_tar) - 1;
x = x + normrnd(0, 0.4, 1);
u = 0.2*tanh(x) + 0.8;
u = (1-m_porb(1).^0.25) .* u;
end
