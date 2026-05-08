%% Step 1: analyze pointwise laws F(x,g) for fixed x
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage0_static_gap_data.mat'), 'gapList', 'xGrid', 'S');

nX = numel(xGrid);
r2LinearG = zeros(nX, 1);
r2LinearInvG = zeros(nX, 1);
r2QuadG = zeros(nX, 1);
r2ExpG = zeros(nX, 1);

for ix = 1:nX
    y = S(:, ix);
    r2LinearG(ix) = fit_r2(y, [ones(size(gapList)), gapList]);
    r2LinearInvG(ix) = fit_r2(y, [ones(size(gapList)), 1 ./ gapList]);
    r2QuadG(ix) = fit_r2(y, [ones(size(gapList)), gapList, gapList.^2]);
    y0 = max(y - min(y) + eps, eps);
    r2ExpG(ix) = fit_r2(log(y0), [ones(size(gapList)), gapList]);
end

[~, peakIdx] = max(max(S, [], 1));
[~, leftIdx] = min(abs(xGrid - (-2.0)));
[~, rightIdx] = min(abs(xGrid - (2.0)));
repIdx = unique([leftIdx, peakIdx, rightIdx]);

pointwiseStats = table(xGrid, r2LinearG, r2LinearInvG, r2QuadG, r2ExpG);
modelSummary = table( ...
    mean(r2LinearG), mean(r2LinearInvG), mean(r2QuadG), mean(r2ExpG), ...
    median(r2LinearG), median(r2LinearInvG), median(r2QuadG), median(r2ExpG), ...
    'VariableNames', {'mean_R2_linear_g', 'mean_R2_linear_inv_g', 'mean_R2_quad_g', 'mean_R2_logexp_g', ...
    'median_R2_linear_g', 'median_R2_linear_inv_g', 'median_R2_quad_g', 'median_R2_logexp_g'});

save(fullfile(outDir, 'stage1_pointwise_laws.mat'), ...
    'pointwiseStats', 'modelSummary', 'repIdx', '-v7.3');

figure('Name', 'Static Gap Step 1 - Pointwise Laws', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1.5, 1.5, 22, 12]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(xGrid, r2LinearG, 'DisplayName', 'linear in g');
plot(xGrid, r2LinearInvG, 'DisplayName', 'linear in 1/g');
plot(xGrid, r2QuadG, 'DisplayName', 'quadratic in g');
plot(xGrid, r2ExpG, 'DisplayName', 'log-linear exp(g)');
xlabel('Position x (mm)');
ylabel('R^2');
title('Pointwise fit quality across x');
legend('Location', 'best', 'Box', 'off');

nexttile;
bar(categorical({'linear g', 'linear 1/g', 'quad g', 'log-exp g'}), ...
    [mean(r2LinearG), mean(r2LinearInvG), mean(r2QuadG), mean(r2ExpG)]);
ylabel('Mean R^2 over x');
title('Average pointwise law quality');

nexttile; hold on;
colors = lines(numel(repIdx));
for k = 1:numel(repIdx)
    ix = repIdx(k);
    plot(gapList, S(:, ix), 'o-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('x = %.2f mm', xGrid(ix)));
end
xlabel('Gap g (mm)');
ylabel('F(x,g)');
title('Representative pointwise curves');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
for k = 1:numel(repIdx)
    ix = repIdx(k);
    y = S(:, ix);
    yHat = fit_predict(y, [ones(size(gapList)), 1 ./ gapList]);
    plot(gapList, y, 'o', 'Color', colors(k, :), 'HandleVisibility', 'off');
    plot(gapList, yHat, '-', 'Color', colors(k, :), ...
        'DisplayName', sprintf('1/g fit at x = %.2f', xGrid(ix)));
end
xlabel('Gap g (mm)');
ylabel('Fitted F(x,g)');
title('Representative linear-in-1/g fits');
legend('Location', 'best', 'Box', 'off');

fprintf('[Static Step 1] Mean R^2:\n');
disp(modelSummary);

function r2 = fit_r2(y, X)
    yHat = fit_predict(y, X);
    ssRes = sum((y - yHat).^2);
    ssTot = sum((y - mean(y)).^2);
    r2 = 1 - ssRes / max(ssTot, eps);
end

function yHat = fit_predict(y, X)
    beta = X \ y;
    yHat = X * beta;
end

