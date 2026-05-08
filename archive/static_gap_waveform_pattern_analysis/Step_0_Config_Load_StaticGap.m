%% Step 0: load static gap waveforms and build a common representation
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

cfg = struct();
cfg.rootDir = rootDir;
cfg.scriptDir = scriptDir;
cfg.outDir = outDir;
cfg.xGridN = 2401;
cfg.levelList = (0.1:0.1:0.9)';

txtCandidates = dir(fullfile(rootDir, '**', '*2mm*.txt'));
txtNames = lower(string({txtCandidates.name}));
isTarget = ~startsWith(txtNames, "wave");
targetIdx = find(isTarget, 1, 'first');
if isempty(targetIdx)
    error('Could not find the stacked straight-blade 2 mm gap data file under: %s', rootDir);
end
cfg.dataFile = fullfile(txtCandidates(targetIdx).folder, txtCandidates(targetIdx).name);

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

xMin = max(cellfun(@min, xCell));
xMax = min(cellfun(@max, xCell));
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, cfg.xGridN)';

S = zeros(numel(gapList), numel(xGrid));
for i = 1:numel(gapList)
    S(i, :) = interp1(xCell{i}, yCell{i}, xGrid, 'pchip');
end

save(fullfile(outDir, 'stage0_static_gap_data.mat'), ...
    'cfg', 'gapList', 'xCell', 'yCell', 'xGrid', 'S', '-v7.3');

figure('Name', 'Static Gap Step 0 - Loaded Curves', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 9]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
colors = parula(numel(gapList));
for i = 1:numel(gapList)
    plot(xCell{i}, yCell{i}, 'Color', colors(i, :), 'DisplayName', sprintf('g = %.1f mm', gapList(i)));
end
xlabel('Position x (mm)');
ylabel('Capacitance');
title('Original static waveforms');
legend('Location', 'eastoutside', 'Box', 'off');

nexttile;
imagesc(xGrid, gapList, S);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Position x (mm)');
ylabel('Gap g (mm)');
title('Common-grid surface F(x,g)');

fprintf('[Static Step 0] Loaded %d gap curves from: %s\n', numel(gapList), cfg.dataFile);

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

