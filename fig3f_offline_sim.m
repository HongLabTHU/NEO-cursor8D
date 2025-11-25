%% Center-out offline simulation

clear;
clc

% confusion
C = load('data\Confusion.mat');


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


%% center-out task
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



%%
figure('Position', [573,546.3,285.3,211.3])
histogram(TTim, 20, 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'Normalization', 'probability');
xline(TTim(end-2), '--', 'Color', [0.8500 0.3250 0.098], 'LineWidth', 1);
text(TTim(end-2)+0.15, 0.14, 'Sing. & Dual. 8D', 'Color', [0.8500 0.3250 0.098], 'FontName', 'Arial', 'Rotation',270);

xline(TTim(end), ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 1);
text(TTim(end)+0.15, 0.14, 'Sing. & Dual. 4D', 'Color', [0.65 0.65 0.65], 'FontName', 'Arial', 'Rotation',270);

xline(TTim(end-1), '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(TTim(end-1)+0.15, 0.14, 'Sing. Move. 4D', 'Color', [0.5, 0.5, 0.5], 'FontName', 'Arial', 'Rotation',270);

xlabel('Hit time (s)', 'FontName', 'Arial');
yticks([0 0.07 0.14]);
axis([4 10 0 0.14]);
box off



%% center-out sim
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



function prob = gen_move_prob(P_tar, PV, prest)

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



function prob = gen_post_porb(m_porb, M)
ind = randsrc(1, 1, [1:size(M, 1); m_porb]);
V = 0 * m_porb;
V(ind) = 1;
prob = V * M;
end


function u = nfunc(P_tar, m_porb)
x = 5*norm(P_tar) - 1;
x = x + normrnd(0, 0.4, 1);
u = 0.2*tanh(x) + 0.8;
u = (1-m_porb(1).^0.25) .* u;
end


