%% 新数据混淆矩阵
clear;
clc
load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\code_0804\neo_vis\TT_LDA_1111.mat')
X = classer.X;
ylab = classer.Y;

ylab(ylab==122) = 101;
ylab(ylab==121) = 102;
ylab(ylab==125) = 102;
ylab(ylab==126) = 104;
condu = [100 101 102 103 104];
Tbox = {'Rest', '  LE\newlineLE&RE',  '  LH\newlineLE&LH\newlineLH&RH', 'R Elbow', '  RH\newlineRE&RH'};

% condu = [100 101 102 103 104];
% Tbox = {'Rest', 'L Elbow', 'L Hand' , 'R Elbow', 'R Elbow'};

% condu = [100 101 102 103 104 121 122 125 126];
% Tbox = {'Rest', 'L Elbow', 'L Hand' , 'R Elbow', 'R Hand', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
X = X(ismember(ylab, condu), :);
ylab = ylab(ismember(ylab, condu));
[macroF1, macroACC, AbM] = cross_vaild(X, ylab, 100);


%%
idx = [1 5 9 4 7 2 6 3 8];
idx = 1:5;
figure('Position', [573,428.3333333333333,386,329]); % 573,344.3,480,412.7
imagesc(AbM(idx, idx), [0 0.8]);
colormap(slanCM('Greys'));
xticks(1:length(condu)); yticks(1:length(condu));
xticklabels(Tbox(idx));
yticklabels(Tbox(idx));
colorbar('limits',  [0 0.8], 'Ticks', [0 0.8]);
title(['Cross-Vaildated Decoding ACC ' num2str(macroACC(1), '%.3f')]);
annotation('textbox', [.89 .6 .1 .05], 'String', 'Confusion', 'EdgeColor', 'none', 'Rotation', 270, 'FontSize', 11);



%% 迭代计算最优的编码
clear;
clc
load('E:\MATLAB_softhub\SPMdataset\NEO_TT01\code_0804\neo_vis\TT_LDA_1111.mat')
X = classer.X;
ylab = classer.Y;
condu = [101 102 103 104 121 122 125 126];
X = X(ismember(ylab, condu), :);
ylab = ylab(ismember(ylab, condu));
[~, ~, Conf_Mat] = cross_vaild(X, ylab, 100);

% Pv 初值
Pv = [-sqrt(2)  sqrt(2);
      -sqrt(2) -sqrt(2);
       sqrt(2)  sqrt(2);
       sqrt(2) -sqrt(2);
      -2        0;
       0        2;
       0       -2;
       2        0];

dx = pi/36;
P_tar = [cos(dx:dx:2*pi); sin(dx:dx:2*pi)]';
P_tar = [P_tar; [cos(pi*3/2:dx:2*pi); sin(pi*3/2:dx:2*pi)]'];


for epoch = 1
    Prob = P_tar * Pv' ./ (sqrt(sum(P_tar.^2, 2)) * sqrt(sum(Pv.^2, 2))');
    Prob = pi - abs(acos(Prob));
    Prob = Prob.^3;
    human_move = zeros(size(Prob));

    % 循环遍历每一行
    for i = 1:size(Prob, 1)
        [sorted_row, idx] = sort(Prob(i, :), 'descend');
        human_move(i, idx(1:3)) = sorted_row(1:3);
        human_move(i, :) = human_move(i, :) ./ sum(human_move(i, :));
    end
end

B = Conf_Mat' * human_move' * P_tar;
for i = 1:8
    Bv(i, :) = B(i, :) / norm(B(i, :));
end


Cmp = [0.9922    0.6824    0.4196; % LE
    0.9020    0.3333    0.0510;  % LH
    0.6196    0.7922    0.8824;  % RE
    0.1922    0.5098    0.7412;  % RH
    0.8398    0.1680    0.1641;  % LE&LH
    0.7305    0.8320    0.5898;
    0.4660    0.6740    0.1880;
    0.0902    0.7451    0.8118];



%%
Cmp(100, :) = [0.2 0.2 0.2];
Cmp([101:104 121 122 125 126], :) = [0.1922    0.5098    0.7412;
                                    0.9020    0.3333    0.0510;
                                    0.6196    0.7922    0.8824;
                                    0.9922    0.6824    0.4196;
                                    0.4660    0.6740    0.1880;
                                    0.0902    0.7451    0.8118;
                                    0.8398    0.1680    0.1641;
                                    0.7305    0.8320    0.5898];
X = classer.X(1:5000,:);
ylab = classer.Y(1:5000,:);
[~, posterior] = predict(classer, X);
V_turning = Bv(:, 1) + 1j*Bv(:, 2);
V = posterior(:, 2:end) * V_turning;
V1 = (1-posterior(:, 1).^0.25) .* V;
polarscatter(angle(V1), abs(V1), 8, Cmp(ylab, :), 'filled', 'MarkerFaceAlpha', 0.5);
thetaticks(0:45:315);
rticks(0:1/3:1); rticklabels([]);




%% Fig 7a
figure('Position', [403,395,238,242]);
Tbox = {'LE', 'LH' , 'RE', 'RH', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};

for i = 1:8
    Bv(i, :) = B(i, :) / norm(B(i, :));
    quiver(0, 0, Bv(i, 1), Bv(i, 2), 'Color', Cmp(i, :), 'LineWidth', 2, 'MaxHeadSize', 0.4);
    text(Bv(i, 1), Bv(i, 2), Tbox{i}, 'Color', Cmp(i, :), 'FontWeight', 'bold', 'HorizontalAlignment','center');
    hold on
end

axis equal 
axis([-1 1 -1 1]);
axis off 

% 
figure('Position', [403,395,238,242]);
Tbox = {'LE', 'LH' , 'RE', 'RH', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
for i = 1:4
    quiver(0, 0, Pv(i, 1)/2, Pv(i, 2)/2, 'Color', Cmp(i, :), 'LineWidth', 2);
    text(Pv(i, 1)/2, Pv(i, 2)/2, Tbox{i}, 'Color', Cmp(i, :), 'FontWeight', 'bold', 'HorizontalAlignment','center');
    text(Pv(i+4, 1)/2, Pv(i+4, 2)/2, Tbox{i+4}, 'Color', Cmp(i+4, :), 'FontWeight', 'bold', 'HorizontalAlignment','center');
    hold on
end

axis equal 
axis([-1 1 -1 1]);
axis off 


