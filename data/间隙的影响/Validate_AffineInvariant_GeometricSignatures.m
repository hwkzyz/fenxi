%% Validate affine-invariant geometric signatures using no-vibration waveforms
% Idea to validate:
%   For a short pulse under circumferential vibration, the observed waveform
%   in x0 can be approximated as a locally re-parameterized version of the
%   no-vibration waveform. Absolute widths are not invariant, but some
%   normalized geometric signatures may be much more stable.
%
% We compare three classes of features under synthetic vibration:
%   1) absolute single equal-level width      -> expected to be unstable
%   2) normalized multi-level width signature -> expected to be much stabler
%   3) normalized cumulative area signature   -> expected to be much stabler
%
% The script uses only no-vibration waveforms, imposes synthetic
% circumferential vibration, and then checks feature stability and gap
% separability.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'affine_invariant_signature_validation_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Figure defaults
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.15, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% Load no-vibration baseline waveforms
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
levels = (0.10:0.10:0.90)';
nLevel = numel(levels);
singleLevelRef = 0.50;
[~, singleLevelIdx] = min(abs(levels - singleLevelRef));

%% Synthetic vibration cases
V_tip = 3.0e5;  % mm/s
cases = struct( ...
    'name', {'Constant shift', 'Low-frequency single-mode', 'User multimode', 'Fast varying multimode'}, ...
    'tag', {'const_shift', 'low_single', 'user_multimode', 'fast_multimode'}, ...
    'A', {0.60, 0.80, [0.5, 0.4], [0.8, 0.5]}, ...
    'f', {0, 250, [500, 1300], [900, 2200]}, ...
    'phi', {0, pi/6, [pi/4, -pi/3], [pi/5, -pi/2]} ...
    );
nCase = numel(cases);

%% Reference signatures from no-vibration data
allWidthSigRef = zeros(nGap, nLevel);
allAreaSigRef = zeros(nGap, nLevel);

for ig = 1:nGap
    xRef = xCell{ig}(:);
    yRef = yCell{ig}(:);
    [baseRef, ampRef, ~] = normalize_waveform(yRef);
    yTargets = baseRef + ampRef * levels;
    [wRef, ~, ~] = level_widths_abs(xRef, yRef, yTargets);
    allWidthSigRef(ig, :) = normalize_signature(wRef);
    allAreaSigRef(ig, :) = cumulative_area_signature(levels, wRef);
end

%% Validation loop
rows(nGap * nCase, 1) = struct( ...
    'gap_mm', NaN, ...
    'case_name', "", ...
    'distortion_eta', NaN, ...
    'single_width_rel_err_percent', NaN, ...
    'width_signature_l2_err', NaN, ...
    'width_signature_cosine_similarity', NaN, ...
    'area_signature_l2_err', NaN, ...
    'area_signature_cosine_similarity', NaN, ...
    'gap_by_width_signature_mm', NaN, ...
    'gap_by_area_signature_mm', NaN, ...
    'width_signature_correct', NaN, ...
    'area_signature_correct', NaN);

repIdx = find(abs(gapList - 0.8) == min(abs(gapList - 0.8)), 1, 'first');
repCaseIdx = 3;  % user multimode
rep = struct();
rowId = 0;

for ig = 1:nGap
    xRef = xCell{ig}(:);
    yRef = yCell{ig}(:);
    [baseRef, ampRef, ~] = normalize_waveform(yRef);
    yTargets = baseRef + ampRef * levels;
    [wRef, ~, ~] = level_widths_abs(xRef, yRef, yTargets);
    sigWRef = normalize_signature(wRef);
    sigARef = cumulative_area_signature(levels, wRef);

    for ic = 1:nCase
        synth = synthesize_vibrating_waveform(xRef, yRef, cases(ic), V_tip);
        [wObs, ~, ~] = level_widths_abs(synth.x0, synth.yObs, yTargets);
        sigWObs = normalize_signature(wObs);
        sigAObs = cumulative_area_signature(levels, wObs);

        distWidth = zeros(nGap, 1);
        distArea = zeros(nGap, 1);
        for j = 1:nGap
            distWidth(j) = norm(sigWObs - allWidthSigRef(j, :));
            distArea(j) = norm(sigAObs - allAreaSigRef(j, :));
        end
        [~, bestWidthIdx] = min(distWidth);
        [~, bestAreaIdx] = min(distArea);

        rowId = rowId + 1;
        rows(rowId).gap_mm = gapList(ig);
        rows(rowId).case_name = cases(ic).name;
        rows(rowId).distortion_eta = synth.eta;
        rows(rowId).single_width_rel_err_percent = ...
            100 * abs(wObs(singleLevelIdx) - wRef(singleLevelIdx)) / max(abs(wRef(singleLevelIdx)), eps);
        rows(rowId).width_signature_l2_err = norm(sigWObs - sigWRef);
        rows(rowId).width_signature_cosine_similarity = ...
            dot(sigWObs, sigWRef) / max(norm(sigWObs) * norm(sigWRef), eps);
        rows(rowId).area_signature_l2_err = norm(sigAObs - sigARef);
        rows(rowId).area_signature_cosine_similarity = ...
            dot(sigAObs, sigARef) / max(norm(sigAObs) * norm(sigARef), eps);
        rows(rowId).gap_by_width_signature_mm = gapList(bestWidthIdx);
        rows(rowId).gap_by_area_signature_mm = gapList(bestAreaIdx);
        rows(rowId).width_signature_correct = double(bestWidthIdx == ig);
        rows(rowId).area_signature_correct = double(bestAreaIdx == ig);

        if ig == repIdx && ic == repCaseIdx
            rep.gap_mm = gapList(ig);
            rep.case_name = cases(ic).name;
            rep.xRef = xRef;
            rep.yRef = yRef;
            rep.x0 = synth.x0;
            rep.yObs = synth.yObs;
            rep.xi = synth.xi;
            rep.wRef = wRef;
            rep.wObs = wObs;
            rep.sigWRef = sigWRef;
            rep.sigWObs = sigWObs;
            rep.sigARef = sigARef;
            rep.sigAObs = sigAObs;
            rep.eta = synth.eta;
        end
    end
end

summaryTable = struct2table(rows);
writetable(summaryTable, fullfile(outDir, 'affine_signature_validation_summary.csv'));

caseSummary = groupsummary(summaryTable, "case_name", "mean", ...
    ["distortion_eta", "single_width_rel_err_percent", "width_signature_l2_err", ...
     "area_signature_l2_err", "width_signature_correct", "area_signature_correct"]);
writetable(caseSummary, fullfile(outDir, 'affine_signature_validation_grouped.csv'));

globalSummary = table( ...
    mean(summaryTable.single_width_rel_err_percent), ...
    mean(summaryTable.width_signature_l2_err), ...
    mean(summaryTable.area_signature_l2_err), ...
    mean(summaryTable.width_signature_correct), ...
    mean(summaryTable.area_signature_correct), ...
    'VariableNames', {'mean_single_width_rel_err_percent', ...
    'mean_width_signature_l2_err', 'mean_area_signature_l2_err', ...
    'width_signature_identification_accuracy', ...
    'area_signature_identification_accuracy'});
writetable(globalSummary, fullfile(outDir, 'affine_signature_validation_global.csv'));

%% Pairwise gap separability
nPairs = nchoosek(nGap, 2);
pairRows = repmat(struct( ...
    'gap_i_mm', 0, ...
    'gap_j_mm', 0, ...
    'ref_width_signature_distance', 0, ...
    'obs_width_signature_distance', 0, ...
    'ref_area_signature_distance', 0, ...
    'obs_area_signature_distance', 0), nPairs, 1);

obsWidthRep = zeros(nGap, nLevel);
obsAreaRep = zeros(nGap, nLevel);
for ig = 1:nGap
    xRef = xCell{ig}(:);
    yRef = yCell{ig}(:);
    synth = synthesize_vibrating_waveform(xRef, yRef, cases(repCaseIdx), V_tip);
    [baseRef, ampRef, ~] = normalize_waveform(yRef);
    yTargets = baseRef + ampRef * levels;
    [wObs, ~, ~] = level_widths_abs(synth.x0, synth.yObs, yTargets);
    obsWidthRep(ig, :) = normalize_signature(wObs);
    obsAreaRep(ig, :) = cumulative_area_signature(levels, wObs);
end

rowId = 0;
for i = 1:nGap
    for j = i+1:nGap
        rowId = rowId + 1;
        pairRows(rowId).gap_i_mm = gapList(i);
        pairRows(rowId).gap_j_mm = gapList(j);
        pairRows(rowId).ref_width_signature_distance = norm(allWidthSigRef(i, :) - allWidthSigRef(j, :));
        pairRows(rowId).obs_width_signature_distance = norm(obsWidthRep(i, :) - obsWidthRep(j, :));
        pairRows(rowId).ref_area_signature_distance = norm(allAreaSigRef(i, :) - allAreaSigRef(j, :));
        pairRows(rowId).obs_area_signature_distance = norm(obsAreaRep(i, :) - obsAreaRep(j, :));
    end
end
pairTable = struct2table(pairRows);
writetable(pairTable, fullfile(outDir, 'affine_signature_pairwise_distance.csv'));

%% Figure 1: representative waveform and signatures
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(rep.xRef, rep.yRef, 'k-', 'DisplayName', 'No vibration');
plot(rep.x0, rep.yObs, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', rep.case_name);
xlabel('Rigid coordinate x_0 (mm)');
ylabel('Capacitance V (pF)');
title('(a) Representative waveform');
legend('Location', 'best', 'FontSize', 7.0);

nexttile; hold on;
plot(levels, rep.wRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference width');
plot(levels, rep.wObs, '-s', 'MarkerSize', 4, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Observed width');
xlabel('Normalized level q');
ylabel('Width in x_0 (mm)');
title('(b) Absolute widths');
legend('Location', 'best', 'FontSize', 7.0);

nexttile; hold on;
plot(levels, rep.sigWRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference');
plot(levels, rep.sigWObs, '-s', 'MarkerSize', 4, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Observed');
xlabel('Normalized level q');
ylabel('Normalized width signature');
title(sprintf('(c) Width signature, \\eta = %.4f', rep.eta));
legend('Location', 'best', 'FontSize', 7.0);

nexttile; hold on;
plot(levels, rep.sigARef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference');
plot(levels, rep.sigAObs, '-d', 'MarkerSize', 4, 'Color', [0.20 0.60 0.40], 'DisplayName', 'Observed');
xlabel('Normalized level q');
ylabel('Normalized cumulative area signature');
title('(d) Area signature');
legend('Location', 'best', 'FontSize', 7.0);

title(tl1, sprintf('Affine-invariant signature validation at g = %.1f mm', rep.gap_mm), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_affine_signature_representative.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_affine_signature_representative.pdf'), 'ContentType', 'vector');

%% Figure 2: grouped case comparison
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
caseCats = categorical(caseSummary.case_name);

nexttile;
bar(caseCats, [caseSummary.mean_single_width_rel_err_percent, ...
    100 * caseSummary.mean_width_signature_l2_err, ...
    100 * caseSummary.mean_area_signature_l2_err], 'grouped');
ylabel('Error metric (%, scaled)');
title('(a) Stability comparison');
legend({'Single width error (%)', '100 x width-signature L2', '100 x area-signature L2'}, ...
    'Location', 'northwest', 'FontSize', 7.0);

nexttile;
bar(caseCats, [caseSummary.mean_width_signature_correct, caseSummary.mean_area_signature_correct], 'grouped');
ylabel('Gap identification accuracy');
ylim([0, 1.05]);
title('(b) Gap separability comparison');
legend({'Width signature', 'Area signature'}, 'Location', 'northwest', 'FontSize', 7.0);

exportgraphics(fig2, fullfile(outDir, 'fig2_affine_signature_grouped.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_affine_signature_grouped.pdf'), 'ContentType', 'vector');

%% Figure 3: pairwise distance preservation
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl3 = tiledlayout(fig3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(pairTable.ref_width_signature_distance, pairTable.obs_width_signature_distance, 'o', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75]);
lims = [0, 1.05 * max(max(pairTable.ref_width_signature_distance), max(pairTable.obs_width_signature_distance))];
plot(lims, lims, 'k--');
xlabel('Reference width-signature distance');
ylabel('Observed width-signature distance');
title('(a) Width-signature separability');

nexttile; hold on;
plot(pairTable.ref_area_signature_distance, pairTable.obs_area_signature_distance, 'o', ...
    'Color', [0.20 0.60 0.40], 'MarkerFaceColor', [0.20 0.60 0.40]);
lims = [0, 1.05 * max(max(pairTable.ref_area_signature_distance), max(pairTable.obs_area_signature_distance))];
plot(lims, lims, 'k--');
xlabel('Reference area-signature distance');
ylabel('Observed area-signature distance');
title('(b) Area-signature separability');

exportgraphics(fig3, fullfile(outDir, 'fig3_affine_signature_pairwise.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_affine_signature_pairwise.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Affine-invariant geometric signature validation ===\n');
fprintf('Data source: no-vibration waveforms in %s\n', dataFile);
fprintf('Tip speed V = %.3e mm/s\n', V_tip);
disp(globalSummary);
fprintf('\nGrouped by vibration case:\n');
disp(caseSummary(:, {'case_name', 'mean_distortion_eta', 'mean_single_width_rel_err_percent', ...
    'mean_width_signature_l2_err', 'mean_area_signature_l2_err', ...
    'mean_width_signature_correct', 'mean_area_signature_correct'}));
fprintf('\nResults saved to: %s\n', outDir);

%% Local functions
function synth = synthesize_vibrating_waveform(xRef, yRef, caseDef, V_tip)
    xRef = xRef(:);
    yRef = yRef(:);
    x0Center = mean(xRef);

    % Find a rigid-coordinate window whose mapped relative coordinate spans
    % the original no-vibration support as fully as possible.
    x0 = xRef;
    for iter = 1:4
        tTmp = (x0 - x0Center) / V_tip;
        [uTmp, ~] = multimode_vibration(tTmp, caseDef.A, caseDef.f, caseDef.phi);
        x0 = linspace(min(xRef) + max(uTmp), max(xRef) + min(uTmp), numel(xRef))';
    end

    t = (x0 - x0Center) / V_tip;
    [u, du_dt] = multimode_vibration(t, caseDef.A, caseDef.f, caseDef.phi);
    xi = x0 - u;
    if any(diff(xi) <= 0) || min(xi) < min(xRef) - 1e-8 || max(xi) > max(xRef) + 1e-8
        error('Window construction failed for case "%s".', caseDef.name);
    end

    synth = struct();
    synth.x0 = x0;
    synth.t = t;
    synth.u = u;
    synth.du_dt = du_dt;
    synth.xi = xi;
    synth.yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    synth.eta = max(abs(du_dt)) / V_tip;
end

function [u, du_dt] = multimode_vibration(t, A, f, phi)
    t = t(:);
    if isscalar(A)
        if f == 0
            u = A * ones(size(t));
            du_dt = zeros(size(t));
        else
            u = A * sin(2 * pi * f * t + phi);
            du_dt = 2 * pi * f * A * cos(2 * pi * f * t + phi);
        end
        return;
    end

    u = zeros(size(t));
    du_dt = zeros(size(t));
    for k = 1:numel(A)
        u = u + A(k) * sin(2 * pi * f(k) * t + phi(k));
        du_dt = du_dt + 2 * pi * f(k) * A(k) * cos(2 * pi * f(k) * t + phi(k));
    end
end

function [gapList, xCell, yCell] = load_stacked_curves(filePath)
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
        xCell{i} = x(idx);
        yCell{i} = y(idx);
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function [widths, xLeft, xRight] = level_widths_abs(x, y, yTargets)
    x = x(:);
    y = y(:);
    yTargets = yTargets(:);
    [~, iPeak] = max(y);

    xRise = x(1:iPeak);
    yRise = y(1:iPeak);
    xFall = x(iPeak:end);
    yFall = y(iPeak:end);

    [yRiseU, idxRiseU] = unique(yRise, 'stable');
    xRiseU = xRise(idxRiseU);
    [yFallU, idxFallU] = unique(yFall, 'stable');
    xFallU = xFall(idxFallU);

    nLevel = numel(yTargets);
    widths = nan(nLevel, 1);
    xLeft = nan(nLevel, 1);
    xRight = nan(nLevel, 1);
    for k = 1:nLevel
        yt = yTargets(k);
        if yt >= min(yRiseU) && yt <= max(yRiseU) && yt >= min(yFallU) && yt <= max(yFallU)
            xLeft(k) = interp1(yRiseU, xRiseU, yt, 'linear');
            xRight(k) = interp1(flipud(yFallU), flipud(xFallU), yt, 'linear');
            widths(k) = xRight(k) - xLeft(k);
        end
    end
end

function sig = normalize_signature(w)
    w = sanitize_widths(w);
    sig = (w / max(sum(w), eps)).';
end

function sig = cumulative_area_signature(~, w)
    w = sanitize_widths(w);
    acc = cumsum(w);
    sig = (acc / max(acc(end), eps)).';
end

function w = sanitize_widths(w)
    w = w(:);
    if all(isfinite(w))
        return;
    end

    idx = (1:numel(w)).';
    valid = isfinite(w);
    if ~any(valid)
        w(:) = 0;
        return;
    end

    w = interp1(idx(valid), w(valid), idx, 'linear', 'extrap');
    w = max(w, 0);
end
