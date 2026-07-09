%% Step05I: Learn offset-and-tilt corrected shared static response surface
% Current shared-surface refinement after Step05.
%
% This step extends Step05H by keeping the waveform-inverted blade clearance
% offsets delta_g_b fixed and learning one blade-specific tilt coefficient:
%
%   g_eff,bj(x) = g_j + delta_g_b + mu_b * x
%
% All blades are then fitted by one shared response surface F(g,x).
% To keep the first feasibility test identifiable, x-shift/scale and
% per-blade voltage gains are not released here.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(thisDir));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step05i_offset_tilt_shared_surface');
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

cfg = struct();
cfg.effectiveWindowThreshold = 0.12;
cfg.minValidSamplesForFit = 16;
cfg.muBounds = [-0.20, 0.20];
cfg.muGrid = linspace(cfg.muBounds(1), cfg.muBounds(2), 161);
cfg.iterations = 4;
cfg.g0Mm = base.g0Mm;

xGrid = base.xGrid(:);
bladeIds = base.bladeIds(:).';
trueGaps = base.trueGapMm(:);
recordedGaps = base.recordedGapMm(:);
waveforms = base.waveforms;

offsetByBlade = load_offset_by_blade_local(offsetTable, bladeIds);
muByBlade = zeros(size(bladeIds));

peakEnvelope = max(reshape(waveforms, numel(xGrid), []), [], 2, 'omitnan');
effectiveWindow = peakEnvelope >= cfg.effectiveWindowThreshold * max(peakEnvelope, [], 'omitnan');
idxWin = find(effectiveWindow);
if ~isempty(idxWin)
    effectiveWindow(idxWin(1):idxWin(end)) = true;
end

fprintf('\n=== Step05I: offset-and-tilt shared response surface ===\n');
fprintf('Step05 source: %s\n', step05File);
fprintf('Offset source: %s\n', offsetFile);
fprintf('Initial offsets:\n');
for ib = 1:numel(bladeIds)
    fprintf('  blade %d: delta_g %.5f mm\n', bladeIds(ib), offsetByBlade(ib));
end

history = table();
coeff = [];
Yfit = [];
fitTable = table();
for iter = 1:cfg.iterations
    [coeff, Yfit, sample] = fit_shared_surface_local( ...
        waveforms, xGrid, trueGaps, recordedGaps, bladeIds, offsetByBlade, muByBlade, effectiveWindow, cfg);
    fitTable = build_fit_table_local(waveforms, Yfit, xGrid, trueGaps, recordedGaps, ...
        bladeIds, offsetByBlade, muByBlade, effectiveWindow);
    bladeSummary = groupsummary(fitTable, 'bladeId', {'mean','std'}, {'rmseMv','relativeRmse'});

    muOld = muByBlade;
    for ib = 1:numel(bladeIds)
        muByBlade(ib) = estimate_mu_for_blade_local( ...
            waveforms(:, ib, :), xGrid, trueGaps, offsetByBlade(ib), coeff, effectiveWindow, cfg);
    end
    meanRmse = mean(fitTable.rmseMv, 'omitnan');
    history = [history; table(iter, meanRmse, norm(muByBlade - muOld), ... %#ok<AGROW>
        'VariableNames', {'iteration','meanRmseMv','muStepNorm'})];
    fprintf('  iter %d: mean RMSE %.2f mV, mu %s\n', ...
        iter, meanRmse, mat2str(muByBlade, 4));
end

[coeff, Yfit, sample] = fit_shared_surface_local( ...
    waveforms, xGrid, trueGaps, recordedGaps, bladeIds, offsetByBlade, muByBlade, effectiveWindow, cfg);
fitTable = build_fit_table_local(waveforms, Yfit, xGrid, trueGaps, recordedGaps, ...
    bladeIds, offsetByBlade, muByBlade, effectiveWindow);
bladeSummary = groupsummary(fitTable, 'bladeId', {'mean','std'}, {'rmseMv','relativeRmse'});

tiltTable = table(bladeIds(:), offsetByBlade(:), muByBlade(:), atan(muByBlade(:)) * 180 / pi, ...
    'VariableNames', {'bladeId','deltaGapOffsetMm','muGapPerXMm','tiltAngleDeg'});

dFdgGrid = nan(numel(xGrid), numel(sample.gEff));
dFdxGrid = nan(numel(xGrid), numel(sample.gEff));
for i = 1:numel(sample.gEff)
    gEffAtX = sample.gEffByX(:, i);
    dFdgGrid(:, i) = -coeff(:, 2) ./ (gEffAtX .^ 2) + coeff(:, 3) ./ gEffAtX;
    dFdxGrid(:, i) = gradient(Yfit(:, i), xGrid);
end

OffsetTiltResponseSurface = base;
OffsetTiltResponseSurface.method = 'offset_tilt_corrected_shared_multi_blade_response_surface';
OffsetTiltResponseSurface.description = ['All static blade waveforms are fitted with one shared response surface. ' ...
    'Blade clearance offsets are taken from reference_blade_gap_analysis; blade tilt coefficients are learned by alternating fit.'];
OffsetTiltResponseSurface.sourceStep05File = step05File;
OffsetTiltResponseSurface.offsetFile = offsetFile;
OffsetTiltResponseSurface.offsetTable = offsetTable;
OffsetTiltResponseSurface.bladeTiltTable = tiltTable;
OffsetTiltResponseSurface.gTrainMm = sample.gEff(:);
OffsetTiltResponseSurface.recordedGapTrainMm = sample.recordedGap(:);
OffsetTiltResponseSurface.sampleBladeId = sample.bladeId(:);
OffsetTiltResponseSurface.gEffByX = sample.gEffByX;
OffsetTiltResponseSurface.Yref = sample.Yall;
OffsetTiltResponseSurface.Yfit = Yfit;
OffsetTiltResponseSurface.coeff = coeff;
OffsetTiltResponseSurface.effectiveWindow = effectiveWindow;
OffsetTiltResponseSurface.dFdgGrid = dFdgGrid;
OffsetTiltResponseSurface.dFdxGrid = dFdxGrid;
OffsetTiltResponseSurface.fitTable = fitTable;
OffsetTiltResponseSurface.bladeFitSummary = bladeSummary;
OffsetTiltResponseSurface.fitHistory = history;
OffsetTiltResponseSurface.cfg = cfg;
OffsetTiltResponseSurface.waveformAggregationLabel = 'offset-and-tilt corrected all-blade shared surface';
OffsetTiltResponseSurface.targetBladeMode = 'offset_tilt_corrected_all_blades';

matFile = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_20250527.mat');
fitCsv = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_Fit_20250527.csv');
tiltCsv = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_Tilt_20250527.csv');
summaryCsv = fullfile(outDir, 'Step05I_OffsetTilt_Shared_Response_Surface_BladeSummary_20250527.csv');
figFile = fullfile(figDir, 'Step05I_OffsetTilt_Shared_Surface_20250527.png');

save(matFile, 'OffsetTiltResponseSurface', 'fitTable', 'tiltTable', 'bladeSummary', 'history', '-v7.3');
writetable(fitTable, fitCsv);
writetable(tiltTable, tiltCsv);
writetable(bladeSummary, summaryCsv);
plot_offset_tilt_surface_local(OffsetTiltResponseSurface, figFile);

fprintf('\nStep05I complete.\n');
disp(tiltTable);
disp(bladeSummary);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, tiltCsv, figFile);

function offsetByBlade = load_offset_by_blade_local(offsetTable, bladeIds)
offsetByBlade = nan(size(bladeIds));
for ib = 1:numel(bladeIds)
    row = offsetTable.bladeId == bladeIds(ib);
    if ~any(row)
        error('No offset entry for bladeId %d.', bladeIds(ib));
    end
    if ismember('meanOffsetZeroedMm', offsetTable.Properties.VariableNames)
        offsetByBlade(ib) = offsetTable.meanOffsetZeroedMm(find(row, 1));
    else
        offsetByBlade(ib) = offsetTable.meanOffsetMm(find(row, 1));
    end
end
end

function [coeff, Yfit, sample] = fit_shared_surface_local(waveforms, xGrid, trueGaps, recordedGaps, bladeIds, offsetByBlade, muByBlade, effectiveWindow, cfg)
nX = numel(xGrid);
sampleGap = [];
sampleBladeId = [];
sampleRecordedGap = [];
Yall = [];
gEffByX = [];
for ib = 1:numel(bladeIds)
    for ig = 1:numel(trueGaps)
        y = waveforms(:, ib, ig);
        if ~all(isfinite(y))
            continue;
        end
        gEff = trueGaps(ig) + offsetByBlade(ib) + muByBlade(ib) .* xGrid;
        if any(gEff <= 0)
            continue;
        end
        sampleGap(end+1, 1) = mean(gEff(effectiveWindow), 'omitnan'); %#ok<AGROW>
        sampleBladeId(end+1, 1) = bladeIds(ib); %#ok<AGROW>
        sampleRecordedGap(end+1, 1) = recordedGaps(ig); %#ok<AGROW>
        Yall(:, end+1) = y; %#ok<AGROW>
        gEffByX(:, end+1) = gEff; %#ok<AGROW>
    end
end

coeff = nan(nX, 3);
Yfit = nan(size(Yall));
for ix = 1:nX
    if ~effectiveWindow(ix)
        continue;
    end
    g = gEffByX(ix, :).';
    basis = [ones(numel(g), 1), 1 ./ g, log(g ./ cfg.g0Mm)];
    y = Yall(ix, :).';
    valid = isfinite(y) & all(isfinite(basis), 2);
    if nnz(valid) < cfg.minValidSamplesForFit
        continue;
    end
    coeff(ix, :) = (basis(valid, :) \ y(valid)).';
    Yfit(ix, valid) = (basis(valid, :) * coeff(ix, :).').';
end

sample = struct('gEff', sampleGap, 'bladeId', sampleBladeId, ...
    'recordedGap', sampleRecordedGap, 'Yall', Yall, 'gEffByX', gEffByX);
end

function mu = estimate_mu_for_blade_local(YbladeRaw, xGrid, trueGaps, offsetMm, coeff, effectiveWindow, cfg)
Yblade = squeeze(YbladeRaw);
bestMu = 0;
bestObj = inf;
for muTry = cfg.muGrid
    obj = blade_mu_objective_local(Yblade, xGrid, trueGaps, offsetMm, coeff, effectiveWindow, cfg, muTry);
    if obj < bestObj
        bestObj = obj;
        bestMu = muTry;
    end
end
fun = @(muRaw) blade_mu_objective_local(Yblade, xGrid, trueGaps, offsetMm, coeff, effectiveWindow, cfg, ...
    min(max(muRaw, cfg.muBounds(1)), cfg.muBounds(2)));
muOpt = fminsearch(fun, bestMu, optimset('Display', 'off', 'MaxIter', 120, 'MaxFunEvals', 320, ...
    'TolX', 1e-5, 'TolFun', 1e-5));
mu = min(max(muOpt, cfg.muBounds(1)), cfg.muBounds(2));
end

function obj = blade_mu_objective_local(Yblade, xGrid, trueGaps, offsetMm, coeff, effectiveWindow, cfg, mu)
pred = nan(size(Yblade));
for ig = 1:numel(trueGaps)
    gEff = trueGaps(ig) + offsetMm + mu .* xGrid;
    if any(gEff(effectiveWindow) <= 0)
        obj = inf;
        return;
    end
    pred(:, ig) = eval_response_surface_by_x_local(xGrid, coeff, cfg.g0Mm, gEff);
end
valid = effectiveWindow & all(isfinite(Yblade), 2) & all(isfinite(pred), 2);
res = Yblade(valid, :) - pred(valid, :);
obj = sqrt(mean(res(:).^2, 'omitnan'));
end

function F = eval_response_surface_by_x_local(xGrid, coeff, g0Mm, gQuery)
F = nan(size(gQuery));
valid = isfinite(gQuery) & isfinite(xGrid) & all(isfinite(coeff), 2);
B0 = coeff(valid, 1);
B1 = coeff(valid, 2);
B2 = coeff(valid, 3);
g = gQuery(valid);
F(valid) = B0 + B1 ./ g + B2 .* log(g ./ g0Mm);
end

function fitTable = build_fit_table_local(waveforms, Yfit, xGrid, trueGaps, recordedGaps, bladeIds, offsetByBlade, muByBlade, effectiveWindow)
rows = {};
iSample = 0;
for ib = 1:numel(bladeIds)
    for ig = 1:numel(trueGaps)
        y = waveforms(:, ib, ig);
        if ~all(isfinite(y))
            continue;
        end
        iSample = iSample + 1;
        keep = effectiveWindow & isfinite(y) & isfinite(Yfit(:, iSample));
        residual = y(keep) - Yfit(keep, iSample);
        rmseMv = sqrt(mean(residual.^2, 'omitnan'));
        relRmse = rmseMv / max(abs(y(keep)), [], 'omitnan');
        rows{end+1, 1} = table(bladeIds(ib), recordedGaps(ig), trueGaps(ig), ...
            offsetByBlade(ib), muByBlade(ib), rmseMv, relRmse, nnz(keep), ...
            'VariableNames', {'bladeId','recordedGapMm','trueReferenceGapMm', ...
            'deltaGapOffsetMm','muGapPerXMm','rmseMv','relativeRmse','pointCount'}); %#ok<AGROW>
    end
end
fitTable = vertcat(rows{:});
end

function plot_offset_tilt_surface_local(S, figFile)
style = paper_style_local();
fig = figure('Name', 'Step05I offset-tilt shared surface', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.5, 13.0]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; box on;
plot(S.fitHistory.iteration, S.fitHistory.meanRmseMv, '-o', 'LineWidth', 1.0, 'MarkerSize', 4);
xlabel('Iteration'); ylabel('Mean RMSE (mV)');
title('Alternating-fit convergence');
format_axes_local(gca, style);

nexttile; hold on; box on;
bar(S.bladeTiltTable.bladeId, S.bladeTiltTable.tiltAngleDeg, 'FaceColor', [0.25 0.45 0.70]);
xlabel('Blade ID'); ylabel('Tilt angle (deg)');
title('Learned blade tilt');
format_axes_local(gca, style);

nexttile; hold on; box on;
summary = S.bladeFitSummary;
bar(summary.bladeId, summary.mean_rmseMv, 'FaceColor', [0.35 0.55 0.35]);
errorbar(summary.bladeId, summary.mean_rmseMv, summary.std_rmseMv, 'k.', 'LineWidth', 0.8);
xlabel('Blade ID'); ylabel('Fit RMSE (mV)');
title('Offset-tilt shared-surface residual by blade');
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
