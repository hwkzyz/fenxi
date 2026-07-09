%% Step05G: Build per-blade static gap response surfaces
% Experimental branch. It does not overwrite Step05 outputs.
% Each static-library blade is fitted as an independent response family:
%
%   F_b(g,x) = B0_b(x) + B1_b(x)/g + B2_b(x)*log(g/g0)
%
% The downstream idea is to use the low-speed rotating template to weight
% these blade families, instead of averaging their waveforms directly.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step05g_per_blade_response_surface');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

step05File = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
if ~isfile(step05File)
    error('Run Step05 first. Missing file: %s', step05File);
end

S = load(step05File, 'responseSurface');
base = S.responseSurface;

minValidGapsForFit = 4;
effectiveWindowThreshold = 0.12;
bladeIds = base.bladeIds(:).';
xGrid = base.xGrid(:);
trueGaps = base.trueGapMm(:);
recordedGaps = base.recordedGapMm(:);
waveforms = base.waveforms;

fprintf('\n=== Step05G: per-blade response surfaces ===\n');
fprintf('Source Step05 file: %s\n', step05File);
fprintf('Blades: %s, gaps: %d, x points: %d\n', mat2str(bladeIds), numel(trueGaps), numel(xGrid));

blade = repmat(empty_blade_surface_local(), numel(bladeIds), 1);
fitRows = {};

for ib = 1:numel(bladeIds)
    Yall = squeeze(waveforms(:, ib, :));
    if isvector(Yall)
        Yall = reshape(Yall, numel(xGrid), []);
    end
    validGapMask = all(isfinite(Yall), 1).';
    gTrain = trueGaps(validGapMask);
    recordedGapTrain = recordedGaps(validGapMask);
    Yref = Yall(:, validGapMask);
    if numel(gTrain) < minValidGapsForFit
        warning('Blade %d skipped: only %d valid gaps.', bladeIds(ib), numel(gTrain));
        continue;
    end

    peakEnvelope = max(Yref, [], 2, 'omitnan');
    peakMax = max(peakEnvelope, [], 'omitnan');
    effectiveWindow = peakEnvelope >= effectiveWindowThreshold * peakMax;
    idxWin = find(effectiveWindow);
    if ~isempty(idxWin)
        effectiveWindow(idxWin(1):idxWin(end)) = true;
    end

    basis = [ones(numel(gTrain), 1), 1 ./ gTrain(:), log(gTrain(:) ./ base.g0Mm)];
    coeff = nan(numel(xGrid), 3);
    Yfit = nan(size(Yref));
    for ix = 1:numel(xGrid)
        if ~effectiveWindow(ix)
            continue;
        end
        y = Yref(ix, :).';
        if nnz(isfinite(y)) < minValidGapsForFit
            continue;
        end
        coeff(ix, :) = (basis \ y).';
        Yfit(ix, :) = (basis * coeff(ix, :).').';
    end

    rmseByGap = nan(numel(gTrain), 1);
    relRmseByGap = nan(numel(gTrain), 1);
    for ig = 1:numel(gTrain)
        keep = effectiveWindow & isfinite(Yfit(:, ig)) & isfinite(Yref(:, ig));
        residual = Yref(keep, ig) - Yfit(keep, ig);
        rmseByGap(ig) = sqrt(mean(residual.^2, 'omitnan'));
        relRmseByGap(ig) = rmseByGap(ig) / max(abs(Yref(keep, ig)), [], 'omitnan');
        fitRows{end+1, 1} = table(bladeIds(ib), recordedGapTrain(ig), gTrain(ig), ...
            rmseByGap(ig), relRmseByGap(ig), nnz(keep), ...
            'VariableNames', {'bladeId','recordedGapMm','trueReferenceGapMm', ...
            'rmseMv','relativeRmse','pointCount'}); %#ok<SAGROW>
    end

    dFdgGrid = nan(numel(xGrid), numel(gTrain));
    dFdxGrid = nan(numel(xGrid), numel(gTrain));
    for ig = 1:numel(gTrain)
        g = gTrain(ig);
        dFdgGrid(:, ig) = -coeff(:, 2) ./ (g .^ 2) + coeff(:, 3) ./ g;
        dFdxGrid(:, ig) = gradient(Yfit(:, ig), xGrid);
    end

    blade(ib) = struct('bladeId', bladeIds(ib), 'bladeIndex', ib, ...
        'xGrid', xGrid, 'g0Mm', base.g0Mm, 'gTrainMm', gTrain(:), ...
        'recordedGapTrainMm', recordedGapTrain(:), 'Yref', Yref, 'Yfit', Yfit, ...
        'coeff', coeff, 'effectiveWindow', effectiveWindow, ...
        'dFdgGrid', dFdgGrid, 'dFdxGrid', dFdxGrid, ...
        'rmseByGapMv', rmseByGap, 'relativeRmseByGap', relRmseByGap, ...
        'sourceCounts', squeeze(base.sourceCounts(ib, validGapMask)).');

    fprintf('  blade %d: mean fit RMSE %.2f mV, valid gaps %d\n', ...
        bladeIds(ib), mean(rmseByGap, 'omitnan'), numel(gTrain));
end

fitTable = vertcat(fitRows{:});

PerBladeResponseLibrary = struct();
PerBladeResponseLibrary.dataset = '20250527';
PerBladeResponseLibrary.method = 'per_blade_static_gap_response_surfaces';
PerBladeResponseLibrary.description = ['Each static-library blade is fitted independently. ' ...
    'No raw waveform averaging across blades is used.'];
PerBladeResponseLibrary.sourceStep05File = step05File;
PerBladeResponseLibrary.baseResponseSurface = base;
PerBladeResponseLibrary.bladeIds = bladeIds;
PerBladeResponseLibrary.xGrid = xGrid;
PerBladeResponseLibrary.trueGapMm = trueGaps;
PerBladeResponseLibrary.recordedGapMm = recordedGaps;
PerBladeResponseLibrary.g0Mm = base.g0Mm;
PerBladeResponseLibrary.blade = blade;
PerBladeResponseLibrary.fitTable = fitTable;

matFile = fullfile(outDir, 'Step05G_PerBlade_Response_Surface_20250527.mat');
csvFile = fullfile(outDir, 'Step05G_PerBlade_Response_Surface_Fit_20250527.csv');
figFile = fullfile(figDir, 'Step05G_PerBlade_Surface_Comparison_20250527.png');

save(matFile, 'PerBladeResponseLibrary', 'fitTable', '-v7.3');
writetable(fitTable, csvFile);
plot_per_blade_library_local(PerBladeResponseLibrary, figFile);

fprintf('\nStep05G complete.\nSaved:\n  %s\n  %s\n  %s\n', matFile, csvFile, figFile);

function blade = empty_blade_surface_local()
blade = struct('bladeId', NaN, 'bladeIndex', NaN, 'xGrid', [], 'g0Mm', NaN, ...
    'gTrainMm', [], 'recordedGapTrainMm', [], 'Yref', [], 'Yfit', [], ...
    'coeff', [], 'effectiveWindow', [], 'dFdgGrid', [], 'dFdxGrid', [], ...
    'rmseByGapMv', [], 'relativeRmseByGap', [], 'sourceCounts', []);
end

function plot_per_blade_library_local(Lib, figFile)
style = paper_style_local();
fig = figure('Name', 'Step05G per-blade response surface', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.5, 11.5]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(numel(Lib.blade));

nexttile; hold on; box on;
for ib = 1:numel(Lib.blade)
    B = Lib.blade(ib);
    if isempty(B.Yref)
        continue;
    end
    [~, ig] = min(abs(B.gTrainMm - median(B.gTrainMm, 'omitnan')));
    plot(B.xGrid, B.Yref(:, ig), '-', 'Color', colors(ib, :), ...
        'LineWidth', 1.0, 'DisplayName', sprintf('B%d', B.bladeId));
end
xlabel('x_{lib} (mm)', 'Interpreter', 'tex');
ylabel('mV');
title('Single-gap waveform family');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
for ib = 1:numel(Lib.blade)
    B = Lib.blade(ib);
    if isempty(B.rmseByGapMv)
        continue;
    end
    plot(B.gTrainMm, B.rmseByGapMv, '-o', 'Color', colors(ib, :), ...
        'LineWidth', 1.0, 'MarkerSize', 3.5, 'DisplayName', sprintf('B%d', B.bladeId));
end
xlabel('g (mm)');
ylabel('Fit RMSE (mV)');
title('Per-blade response-surface fit quality');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.tickFontSize = 8.2;
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on');
grid(ax, 'off');
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
end

function export_paper_figure_local(fig, pngFile)
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
try
    exportgraphics(fig, fullfile(folder, [name, '.pdf']), 'ContentType', 'vector');
catch
end
try
    print(fig, fullfile(folder, [name, '.emf']), '-dmeta', '-r300');
catch
end
end
