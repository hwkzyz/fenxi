%% Step 1: build the leave-one-gap-out static template library
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'gapaware_varpro_holdout_results_modular');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
assert(any(abs(gapList - cfg.g_low) < 1e-12), 'g_low must exist in gapList.');
assert(any(abs(gapList - cfg.g_holdout) < 1e-12), 'g_holdout must exist in gapList.');

templateLib = build_gap_template_library(gapList, xCell, yCell, cfg.g_holdout, cfg.xGridN);
templateLibFull = build_gap_template_library(gapList, xCell, yCell, NaN, cfg.xGridN);

save(fullfile(outDir, 'stage1_template_library.mat'), ...
    'gapList', 'xCell', 'yCell', 'templateLib', 'templateLibFull', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Step 1 - Template Library', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 20, 10]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    for i = 1:numel(gapList)
        if abs(gapList(i) - cfg.g_holdout) < 1e-12
            plot(xCell{i}, yCell{i}, 'r--', 'LineWidth', 1.7, ...
                'DisplayName', sprintf('holdout %.1f mm', gapList(i)));
        elseif abs(gapList(i) - cfg.g_low) < 1e-12
            plot(xCell{i}, yCell{i}, 'b-', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('low ref %.1f mm', gapList(i)));
        else
            plot(xCell{i}, yCell{i}, '-', 'Color', [0.65 0.65 0.65], ...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    xlabel('Position x (mm)');
    ylabel('Capacitance');
    title('Static curves: training vs held-out gap');
    legend('Location', 'best', 'Box', 'off');

    nexttile;
    imagesc(templateLib.xGrid, templateLib.gapTrain, templateLib.S);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('Position x (mm)');
    ylabel('Training gap g (mm)');
    title('Training template surface F(x,g)');
end

fprintf('[Step 1] Training gaps: %s mm\n', mat2str(templateLib.gapTrain(:)', 3));
fprintf('[Step 1] Held-out gap %.3f mm is not used in the template library.\n', cfg.g_holdout);

function [gapList, xCell, yCell] = load_stacked_gap_curves(filePath)
data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
x = data(:, 1);
y = data(:, 2);
breakIdx = [find(diff(x) < 0); numel(x)];
startIdx = [1; breakIdx(1:end-1) + 1];
nCurve = numel(breakIdx);
if nCurve == 7
    gapList = (0.2:0.2:1.4)';
else
    gapList = (1:nCurve)';
end
xCell = cell(nCurve, 1);
yCell = cell(nCurve, 1);
for i = 1:nCurve
    idx = startIdx(i):breakIdx(i);
    [xUnique, ia] = unique(x(idx), 'stable');
    yUnique = y(idx);
    xCell{i} = xUnique(:);
    yCell{i} = yUnique(ia);
end
end

function templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN)
if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
trainIds = find(trainIdx);
G_train = gapList(trainIds);
xMin = -inf;
xMax = inf;
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(G_train), numel(xGrid));
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    S(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ii = 1:size(S, 1)
    FxMat(ii, :) = gradient(S(ii, :), dx);
end
templateLib = struct();
templateLib.gapTrain = G_train(:);
templateLib.xGrid = xGrid(:);
templateLib.S = S;
templateLib.FxMat = FxMat;
templateLib.F = griddedInterpolant({G_train(:), xGrid(:)}, S, 'linear', 'nearest');
templateLib.Fx = griddedInterpolant({G_train(:), xGrid(:)}, FxMat, 'linear', 'nearest');
templateLib.domain = [min(xGrid), max(xGrid)];
end

