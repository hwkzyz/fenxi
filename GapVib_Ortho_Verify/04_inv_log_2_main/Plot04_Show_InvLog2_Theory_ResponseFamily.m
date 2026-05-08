%% Plot 04: from finite gap library to continuous inv_log_2 template family.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
gapList = ctx.gapList(:);
xCell = ctx.xCell;
yCell = ctx.yCell;
cfgAna = ctx.cfgAna;

templateBaseline = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);
templateInv = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);
templateInv.gapInterpMode = "field_basis";
templateInv.gapInterpBasisName = "inv_log_2";
templateInv.gapInterpBasisOrder = 3;

xGrid = templateInv.xGrid(:);
gQuery = linspace(min(gapList), max(gapList), 120);
gSample = [0.2, 0.4, 1.0, 1.4];
gSample = intersect(gSample, gapList.', 'stable');

[~, peakIdx] = max(mean(templateInv.S, 1));
xPick = [xGrid(round(0.25 * numel(xGrid))), xGrid(peakIdx), xGrid(round(0.75 * numel(xGrid)))];

fig = figure('Name', 'inv_log_2 response family', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
libraryColors = turbo(numel(gapList));
for ig = 1:numel(gapList)
    plot(ax1, xCell{ig}(:), yCell{ig}(:), '-', 'Color', libraryColors(ig, :), ...
        'LineWidth', 1.0, 'DisplayName', sprintf('g = %.1f mm', gapList(ig)));
end
title(ax1, '(a) Finite static gap library', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'eastoutside');

ax2 = nexttile;
hold(ax2, 'on');
sampleColors = lines(numel(gSample));
for ig = 1:numel(gSample)
    g = gSample(ig);
    yPred = eval_gap_template(templateInv, g, xGrid);
    plot(ax2, xGrid, yPred, '-', 'Color', sampleColors(ig, :), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('Fhat(g,x), g = %.1f mm', g));
end
title(ax2, '(b) Continuous inv\_log\_2 template family slices', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax2, 'Voltage (V)', 'FontSize', 9);
legend(ax2, 'Location', 'best');

ax3 = nexttile;
hold(ax3, 'on');
trainValues = zeros(numel(gapList), numel(xPick));
baselineFit = zeros(numel(gQuery), numel(xPick));
invFit = zeros(numel(gQuery), numel(xPick));
for ix = 1:numel(xPick)
    for ig = 1:numel(gapList)
        trainValues(ig, ix) = interp1(xCell{ig}(:), yCell{ig}(:), xPick(ix), 'pchip');
    end
    for ig = 1:numel(gQuery)
        baselineFit(ig, ix) = eval_gap_template(templateBaseline, gQuery(ig), xPick(ix));
        invFit(ig, ix) = eval_gap_template(templateInv, gQuery(ig), xPick(ix));
    end
end
markerSet = {'o', 's', '^'};
fitColors = [0.35 0.35 0.35; 0.00 0.45 0.74; 0.85 0.33 0.10];
for ix = 1:numel(xPick)
    plot(ax3, gapList, trainValues(:, ix), markerSet{ix}, ...
        'Color', fitColors(ix, :), 'MarkerFaceColor', fitColors(ix, :), ...
        'LineStyle', 'none', 'DisplayName', sprintf('library at x = %.3f mm', xPick(ix)));
    plot(ax3, gQuery, baselineFit(:, ix), '--', 'Color', fitColors(ix, :), 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax3, gQuery, invFit(:, ix), '-', 'Color', fitColors(ix, :), 'LineWidth', 1.4, ...
        'HandleVisibility', 'off');
end
title(ax3, '(c) Gap-direction response fit at representative x', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax3, 'Gap g (mm)', 'FontSize', 9);
ylabel(ax3, 'Voltage (V)', 'FontSize', 9);
legend(ax3, 'Location', 'best');
text(ax3, 0.05, 0.08, 'dashed: baseline [1,1/g], solid: inv\_log\_2 [1,1/g,log(g)]', ...
    'Units', 'normalized', 'FontSize', 8, 'FontName', 'Times New Roman');

ax4 = nexttile;
FMat = zeros(numel(gQuery), numel(xGrid));
for ig = 1:numel(gQuery)
    FMat(ig, :) = eval_gap_template(templateInv, gQuery(ig), xGrid);
end
imagesc(ax4, xGrid, gQuery, FMat);
set(ax4, 'YDir', 'normal');
hold(ax4, 'on');
for ig = 1:numel(gapList)
    yline(ax4, gapList(ig), 'w:', 'LineWidth', 0.7, 'HandleVisibility', 'off');
end
cb = colorbar(ax4);
cb.Label.String = 'Voltage (V)';
title(ax4, '(d) Continuous template family Fhat(g,x)', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax4, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax4, 'Gap g (mm)', 'FontSize', 9);

apply_inv_log_2_figure_style(ax1, 170, 120);
