%% long term decoding
model = load('data\TT_LDA.mat');
LBmat = [1 2 3 4 7 8 10 11 12 13 15:21; 22:25 27:39];

for h = 1:17
    idx = [Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
    idx = idx(:);
    X = permute(Spectra(:, :, idx), [3 1 2]);
    X = log10(X(:, :));
    ylab = Ylab(idx);
    [macroF1(h), macroACC(h)] = longterm_vaild(model, X, ylab);
end


%%
idx = [];
for h = 3:6
    idx = [idx Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
end
X = permute(Spectra(:, :, idx), [3 1 2]);
X = log10(X(:, :));
ylab = Ylab(idx);
[~, macroACC1, C1] = longterm_vaild(model, X, ylab);

idx = [];
for h = 14:17
    idx = [idx Dirs{LBmat(1, h), 2} Dirs{LBmat(2, h), 2}];
end
X = permute(Spectra(:, :, idx), [3 1 2]);
X = log10(X(:, :));
ylab = Ylab(idx);
[~, macroACC2, C2] = longterm_vaild(model, X, ylab);



%%
figure('Position', [573,428.3333333333333,386,329]);
idx = [1 5 9 4 7 2 6 3 8];
imagesc(C1(idx, idx), [0 0.8]);
colormap(slanCM('Greys'));
xticks(1:9); yticks(1:9);
Tbox = {'Rest', 'LE', 'LE' , 'RE', 'RE', 'LE&LH', 'LE&RE', 'LH&RH', 'RE&RH'};
xticklabels(Tbox(idx));
yticklabels(Tbox(idx));
colorbar('limits',  [0 0.8], 'Ticks', [0 0.8]);
annotation('textbox', [.89 .65 .1 .05], 'String', 'Confusion', 'EdgeColor', 'none', 'Rotation', 270, 'FontSize', 11);



%%
figure('Position', [573,520.3333333333333,380,237.3333333333333])
D1 = D(LBmat(1,:));
s = scatter(D1, macroF1, 20, 'filled', ...
    'MarkerFaceColor', [0.5 0.5 0.5], ...  
    'MarkerEdgeColor', [0.5 0.5 0.5], ... 
    'MarkerFaceAlpha', 0.4);               
hold on;


p = polyfit(D1, macroF1, 1);                
yfit = polyval(p, D1);                   
plot(D1, yfit, 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
R2 = corr(D1', macroF1')^2;
text(0.05, 0.95, sprintf('R^2 = %.3f', R2), ...
    'Units', 'normalized', 'FontSize', 10, ...
    'BackgroundColor', [1 1 1 0.7]);      
box off
axis([-1 200 0 0.6]);
xticks([0 50 100 150 200]);
xlabel('Decoder time span (Days)', 'FontName', 'Arial');
ylabel('Macro F1-Score', 'FontName', 'Arial');



%%
function [macroF1, macroACC, C] = longterm_vaild(model, X, ylab)
test_X = bsxfun(@rdivide,  bsxfun(@minus, X, model.meanX), model.stdX);
[pred, ~] = predict(model.classer, test_X);
C = confusionmat(ylab, pred);
[macroF1, macroACC] = computeMacroF1AndACC(C);
C = C ./ repmat(sum(C, 2), 1, size(C, 1));
end
