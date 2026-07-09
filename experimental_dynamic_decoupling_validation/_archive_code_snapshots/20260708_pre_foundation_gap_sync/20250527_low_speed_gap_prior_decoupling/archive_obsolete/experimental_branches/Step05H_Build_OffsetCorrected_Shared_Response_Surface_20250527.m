%% Step05H: Build offset-corrected shared static response surface
% Experimental branch. It does not overwrite Step05/Step05G outputs.
%
% Different static-library blades are assumed to share the same blade
% thickness and response physics. Blade-to-blade differences are first
% represented by an initial equivalent-clearance offset delta_g_b obtained
% from reference_blade_gap_analysis. The shared response surface is fitted
% using all blade waveforms after replacing each nominal gap g_j with
%
%   g_eff,bj = g_j + delta_g_b
%
% Sensor/blade tilt is not estimated in this step; it should be handled by
% Step06-style low-speed-template calibration.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(thisDir));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step05h_offset_corrected_shared_surface');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

step05File = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
offsetFile = fullfile(rootDir, 'reference_blade_gap_analysis', 'results', ...
    'analysis_04_blade_offset_summary_matched.csv');
if ~isfile(step05File)
    error('Run Step05 first. Missing file: %s', step05File);
end
if ~isfile(offsetFile)
    error('Missing blade offset table: %s', offsetFile);
end

S = load(step05File, 'responseSurface');
base = S.responseSurface;
offsetTable = readtable(offsetFile);

minValidSamplesForFit = 16;
effectiveWindowThreshold = 0.12;
g0Mm = base.g0Mm;
xGrid = base.xGrid(:);
bladeIds = base.bladeIds(:).';
trueGaps = base.trueGapMm(:);
recordedGaps = base.recordedGapMm(:);
waveforms = base.waveforms;

offsetByBlade = nan(size(bladeIds));
for ib = 1:numel(bladeIds)
    row = offsetTable.bladeId == bladeIds(ib);
    if ~any(row)
        error('No offset entry for bladeId %d in %s.', bladeIds(ib), offsetFile);
    end
    if ismember('meanOffsetZeroedMm', offsetTable.Properties.VariableNames)
        offsetByBlade(ib) = offsetTable.meanOffsetZeroedMm(find(row, 1));
    else
        offsetByBlade(ib) = offsetTable.meanOffsetMm(find(row, 1));
    end
end

fprintf('\n=== Step05H: offset-corrected shared response surface ===\n');
fprintf('Step05 source: %s\n', step05File);
fprintf('Offset source: %s\n', offsetFile);
for ib = 1:numel(bladeIds)
    fprintf('  blade %d: delta_g %.5f mm\n', bladeIds(ib), offsetByBlade(ib));
end

peakEnvelope = max(reshape(waveforms, numel(xGrid), []), [], 2, 'omitnan');
effectiveWindow = peakEnvelope >= effectiveWindowThreshold * max(peakEnvelope, [], 'omitnan');
idxWin = find(effectiveWindow);
if ~isempty(idxWin)
    effectiveWindow(idxWin(1):idxWin(end)) = true;
end

sampleGap = [];
sampleBladeId = [];
sampleRecordedGap = [];
Yall = [];
for ib = 1:numel(bladeIds)
    for ig = 1:numel(trueGaps)
        y = waveforms(:, ib, ig);
        if ~all(isfinite(y))
            continue;
        end
        sampleGap(end+1, 1) = trueGaps(ig) + offsetByBlade(ib); %#ok<SAGROW>
        sampleBladeId(end+1, 1) = bladeIds(ib); %#ok<SAGROW>
        sampleRecordedGap(end+1, 1) = recordedGaps(ig); %#ok<SAGROW>
        Yall(:, end+1) = y; %#ok<SAGROW>
    end
end

validSample = isfinite(sampleGap) & sampleGap > 0;
sampleGap = sampleGap(validSample);
sampleBladeId = sampleBladeId(validSample);
sampleRecordedGap = sampleRecordedGap(validSample);
Yall = Yall(:, validSample);

basis = [ones(numel(sampleGap), 1), 1 ./ sampleGap(:), log(sampleGap(:) ./ g0Mm)];
coeff = nan(numel(xGrid), 3);
Yfit = nan(size(Yall));
for ix = 1:numel(xGrid)
    if ~effectiveWindow(ix)
        continue;
    end
    y = Yall(ix, :).';
    valid = isfinite(y) & all(isfinite(basis), 2);
    if nnz(valid) < minValidSamplesForFit
        continue;
    end
    coeff(ix, :) = (basis(valid, :) \ y(valid)).';
    Yfit(ix, valid) = (basis(valid, :) * coeff(ix, :).').';
end

fitRows = cell(numel(sampleGap), 1);
for i = 1:numel(sampleGap)
    keep = effectiveWindow & isfinite(Yall(:, i)) & isfinite(Yfit(:, i));
    residual = Yall(keep, i) - Yfit(keep, i);
    rmseMv = sqrt(mean(residual.^2, 'omitnan'));
    relRmse = rmseMv / max(abs(Yall(keep, i)), [], 'omitnan');
    fitRows{i} = table(sampleBladeId(i), sampleRecordedGap(i), sampleGap(i), ...
        rmseMv, relRmse, nnz(keep), ...
        'VariableNames', {'bladeId','recordedGapMm','effectiveGapMm', ...
        'rmseMv','relativeRmse','pointCount'});
end
fitTable = vertcat(fitRows{:});

bladeSummary = groupsummary(fitTable, 'bladeId', {'mean','std'}, {'rmseMv','relativeRmse'});

dFdgGrid = nan(numel(xGrid), numel(sampleGap));
dFdxGrid = nan(numel(xGrid), numel(sampleGap));
for i = 1:numel(sampleGap)
    g = sampleGap(i);
    dFdgGrid(:, i) = -coeff(:, 2) ./ (g .^ 2) + coeff(:, 3) ./ g;
    dFdxGrid(:, i) = gradient(Yfit(:, i), xGrid);
end

OffsetCorrectedResponseSurface = base;
OffsetCorrectedResponseSurface.method = 'offset_corrected_shared_multi_blade_response_surface';
OffsetCorrectedResponseSurface.description = ['All static blade waveforms are fitted with one shared response surface ' ...
    'after applying blade-specific initial equivalent-clearance offsets from reference_blade_gap_analysis.'];
OffsetCorrectedResponseSurface.sourceStep05File = step05File;
OffsetCorrectedResponseSurface.offsetFile = offsetFile;
OffsetCorrectedResponseSurface.offsetTable = offsetTable;
OffsetCorrectedResponseSurface.offsetByBlade = table(bladeIds(:), offsetByBlade(:), ...
    'VariableNames', {'bladeId','deltaGapOffsetMm'});
OffsetCorrectedResponseSurface.gTrainMm = sampleGap(:);
OffsetCorrectedResponseSurface.recordedGapTrainMm = sampleRecordedGap(:);
OffsetCorrectedResponseSurface.sampleBladeId = sampleBladeId(:);
OffsetCorrectedResponseSurface.Yref = Yall;
OffsetCorrectedResponseSurface.Yfit = Yfit;
OffsetCorrectedResponseSurface.coeff = coeff;
OffsetCorrectedResponseSurface.effectiveWindow = effectiveWindow;
OffsetCorrectedResponseSurface.dFdgGrid = dFdgGrid;
OffsetCorrectedResponseSurface.dFdxGrid = dFdxGrid;
OffsetCorrectedResponseSurface.fitTable = fitTable;
OffsetCorrectedResponseSurface.bladeFitSummary = bladeSummary;
OffsetCorrectedResponseSurface.waveformAggregationLabel = 'offset-corrected all-blade shared surface';
OffsetCorrectedResponseSurface.targetBladeMode = 'offset_corrected_all_blades';

matFile = fullfile(outDir, 'Step05H_OffsetCorrected_Shared_Response_Surface_20250527.mat');
csvFile = fullfile(outDir, 'Step05H_OffsetCorrected_Shared_Response_Surface_Fit_20250527.csv');
summaryCsv = fullfile(outDir, 'Step05H_OffsetCorrected_Shared_Response_Surface_BladeSummary_20250527.csv');
figFile = fullfile(figDir, 'Step05H_OffsetCorrected_Shared_Surface_20250527.png');

save(matFile, 'OffsetCorrectedResponseSurface', 'fitTable', 'bladeSummary', '-v7.3');
writetable(fitTable, csvFile);
writetable(bladeSummary, summaryCsv);
plot_offset_corrected_surface_local(OffsetCorrectedResponseSurface, figFile);

fprintf('\nStep05H complete.\n');
disp(bladeSummary);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, csvFile, figFile);

function plot_offset_corrected_surface_local(S, figFile)
style = paper_style_local();
fig = figure('Name', 'Step05H offset-corrected shared surface', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.5, 11.5]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(numel(S.bladeIds));

targetGap = median(S.gTrainMm, 'omitnan');
nexttile; hold on; box on;
for ib = 1:numel(S.bladeIds)
    rows = S.sampleBladeId == S.bladeIds(ib);
    if ~any(rows)
        continue;
    end
    idx = find(rows);
    [~, k] = min(abs(S.gTrainMm(idx) - targetGap));
    i = idx(k);
    plot(S.xGrid, S.Yref(:, i), '.', 'Color', colors(ib, :), ...
        'MarkerSize', 3.0, 'DisplayName', sprintf('B%d obs', S.bladeIds(ib)));
    plot(S.xGrid, S.Yfit(:, i), '-', 'Color', colors(ib, :), ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
end
xlabel('x_{lib} (mm)', 'Interpreter', 'tex'); ylabel('mV');
title(sprintf('Observed/fitted waveforms near g = %.3f mm', targetGap));
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
summary = S.bladeFitSummary;
bar(summary.bladeId, summary.mean_rmseMv, 'FaceColor', [0.25 0.45 0.70]);
errorbar(summary.bladeId, summary.mean_rmseMv, summary.std_rmseMv, ...
    'k.', 'LineWidth', 0.8);
xlabel('Blade ID'); ylabel('Fit RMSE (mV)');
title('Shared-surface residual by blade');
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
