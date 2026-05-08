%% Analyze stability of normalized multi-level width signatures
% We test whether multi-level equal-capacitance width signatures remain
% much more stable than single widths under multimode circumferential
% vibration, and whether they can still identify the underlying gap class.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'multilevel_signature_stability_results');
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
    'DefaultLineLineWidth', 1.2, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% Vibration settings from the user
u0 = 0.0;               % mm
A_vib = [0.5, 0.4];     % mm
f_vib = [500, 1300];    % Hz
phi0 = [pi/4, -pi/3];   % rad
V_tip = 3.0e5;          % mm/s

levels = (0.10:0.10:0.90)';
singleLevelRef = 0.50;
[~, singleLevelIdx] = min(abs(levels - singleLevelRef));

%% Load all no-vibration baseline waveforms
[gapList, xCell, yCell] = load_stacked_comsol_curves(dataFile);
nGap = numel(gapList);
nLevel = numel(levels);

uMax = abs(u0) + sum(abs(A_vib));
summaryRows(nGap, 1) = struct( ...
    'gap_mm', NaN, ...
    'single_width_abs_err_mm', NaN, ...
    'single_width_rel_err_percent', NaN, ...
    'signature_l2_err', NaN, ...
    'signature_cosine_similarity', NaN, ...
    'identified_gap_by_signature_mm', NaN, ...
    'signature_identification_correct', NaN);

repIdx = find(abs(gapList - 0.8) == min(abs(gapList - 0.8)), 1, 'first');
rep = struct();

allSignatureRef = zeros(nGap, nLevel);
allSignatureObs = zeros(nGap, nLevel);

for ig = 1:nGap
    xRef = xCell{ig}(:);
    yRef = yCell{ig}(:);
    [baseRef, ampRef, ~] = normalize_waveform(yRef);
    yTargets = baseRef + ampRef * levels;
    [wRef, ~, ~] = level_widths_abs(xRef, yRef, yTargets);
    sigRef = normalize_signature(wRef);

    x0Grid = linspace(min(xRef) + uMax, max(xRef) - uMax, numel(xRef))';
    tGrid = (x0Grid - mean(x0Grid)) / V_tip;
    uGrid = u0 + multimode_vibration(tGrid, A_vib, f_vib, phi0);
    xiGrid = x0Grid - uGrid;

    if any(diff(xiGrid) <= 0)
        error('Monotonicity lost for gap %.1f mm. Increase V_tip or reduce amplitudes.', gapList(ig));
    end

    yObs = interp1(xRef, yRef, xiGrid, 'pchip', 'extrap');
    [wObs, ~, ~] = level_widths_abs(x0Grid, yObs, yTargets);
    sigObs = normalize_signature(wObs);

    allSignatureRef(ig, :) = sigRef;
    allSignatureObs(ig, :) = sigObs;

    l2Err = norm(sigObs - sigRef);
    cosSim = dot(sigObs, sigRef) / max(norm(sigObs) * norm(sigRef), eps);

    % Identify the most likely gap using the normalized width signature.
    distToAll = zeros(nGap, 1);
    for j = 1:nGap
        distToAll(j) = norm(sigObs - normalize_signature(signature_from_gap(xCell{j}, yCell{j}, levels)));
    end
    [~, bestIdx] = min(distToAll);

    summaryRows(ig).gap_mm = gapList(ig);
    summaryRows(ig).single_width_abs_err_mm = abs(wObs(singleLevelIdx) - wRef(singleLevelIdx));
    summaryRows(ig).single_width_rel_err_percent = 100 * abs(wObs(singleLevelIdx) - wRef(singleLevelIdx)) / max(abs(wRef(singleLevelIdx)), eps);
    summaryRows(ig).signature_l2_err = l2Err;
    summaryRows(ig).signature_cosine_similarity = cosSim;
    summaryRows(ig).identified_gap_by_signature_mm = gapList(bestIdx);
    summaryRows(ig).signature_identification_correct = double(bestIdx == ig);

    if ig == repIdx
        rep.gap_mm = gapList(ig);
        rep.wRef = wRef;
        rep.wObs = wObs;
        rep.sigRef = sigRef;
        rep.sigObs = sigObs;
        rep.yRef = yRef;
        rep.yObs = yObs;
        rep.xRef = xRef;
        rep.x0Grid = x0Grid;
        rep.xiGrid = xiGrid;
    end
end

summaryTable = struct2table(summaryRows);
writetable(summaryTable, fullfile(outDir, 'signature_stability_summary.csv'));

globalSummary = table( ...
    mean(summaryTable.single_width_abs_err_mm), ...
    mean(summaryTable.single_width_rel_err_percent), ...
    mean(summaryTable.signature_l2_err), ...
    mean(summaryTable.signature_cosine_similarity), ...
    mean(summaryTable.signature_identification_correct), ...
    'VariableNames', {'mean_single_width_abs_err_mm', 'mean_single_width_rel_err_percent', ...
    'mean_signature_l2_err', 'mean_signature_cosine_similarity', ...
    'signature_identification_accuracy'});
writetable(globalSummary, fullfile(outDir, 'signature_stability_global_summary.csv'));

% Pairwise gap separability in signature space
nPairs = nchoosek(nGap, 2);
pairRows = repmat(struct( ...
    'gap_i_mm', 0, ...
    'gap_j_mm', 0, ...
    'ref_signature_distance', 0, ...
    'obs_signature_distance', 0), nPairs, 1);
rowId = 0;
for i = 1:nGap
    for j = i+1:nGap
        rowId = rowId + 1;
        pairRows(rowId).gap_i_mm = gapList(i);
        pairRows(rowId).gap_j_mm = gapList(j);
        pairRows(rowId).ref_signature_distance = norm(allSignatureRef(i, :) - allSignatureRef(j, :));
        pairRows(rowId).obs_signature_distance = norm(allSignatureObs(i, :) - allSignatureObs(j, :));
    end
end
pairTable = struct2table(pairRows);
writetable(pairTable, fullfile(outDir, 'signature_pairwise_distance.csv'));

%% Figure 1: representative waveform and widths
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 12]);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(rep.xRef, rep.yRef, 'k-', 'DisplayName', 'No vibration');
plot(rep.x0Grid, rep.yObs, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Observed in x_0');
xlabel('Rigid coordinate x_0 (mm)');
ylabel('Capacitance (pF)');
title('(a) Representative waveform');
legend('Location', 'northeast', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, rep.wRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference widths');
plot(levels, rep.wObs, '-s', 'MarkerSize', 4, 'Color', [0.10 0.45 0.75], 'DisplayName', 'Observed widths');
xlabel('Normalized level q');
ylabel('Width in x_0 (mm)');
title('(b) Widths at representative gap');
legend('Location', 'best', 'FontSize', 7.2);

nexttile; hold on;
plot(levels, rep.sigRef, 'k--o', 'MarkerSize', 4, 'DisplayName', 'Reference signature');
plot(levels, rep.sigObs, '-s', 'MarkerSize', 4, 'Color', [0.85 0.35 0.10], 'DisplayName', 'Observed signature');
xlabel('Normalized level q');
ylabel('Normalized width signature');
title('(c) Width signature at representative gap');
legend('Location', 'best', 'FontSize', 7.2);

nexttile;
yyaxis left;
bar(categorical(string(summaryTable.gap_mm)), summaryTable.single_width_rel_err_percent, 0.55, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Single-width error (%)');
yyaxis right;
plot(categorical(string(summaryTable.gap_mm)), summaryTable.signature_l2_err, '-o', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10]);
ylabel('Signature L2 error');
title('(d) Single width vs signature stability');
exportgraphics(fig1, fullfile(outDir, 'fig1_signature_stability.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_signature_stability.pdf'), 'ContentType', 'vector');

%% Figure 2: all gap signatures and pairwise distance
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(nGap);

nexttile; hold on;
for ig = 1:nGap
    plot(levels, allSignatureRef(ig, :), '--', 'Color', colors(ig, :), ...
        'DisplayName', sprintf('ref g=%.1f', gapList(ig)));
    plot(levels, allSignatureObs(ig, :), '-', 'Color', colors(ig, :), ...
        'HandleVisibility', 'off');
end
xlabel('Normalized level q');
ylabel('Normalized width signature');
title('(a) Reference and observed signatures');
legend('Location', 'eastoutside', 'FontSize', 7.0);

nexttile; hold on;
plot(pairTable.ref_signature_distance, pairTable.obs_signature_distance, 'o', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75]);
lims = [0, 1.05 * max(max(pairTable.ref_signature_distance), max(pairTable.obs_signature_distance))];
plot(lims, lims, 'k--');
xlabel('Reference pairwise signature distance');
ylabel('Observed pairwise signature distance');
title('(b) Gap separability in signature space');
exportgraphics(fig2, fullfile(outDir, 'fig2_signature_separability.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_signature_separability.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Multi-level width signature stability analysis ===\n');
fprintf('Vibration: u(t) = %.3f + sum A_m sin(2*pi*f_m*t + phi_m)\n', u0);
fprintf('A = [%s] mm\n', num2str(A_vib, '%.3f '));
fprintf('f = [%s] Hz\n', num2str(f_vib, '%.1f '));
fprintf('phi = [%s] rad\n', num2str(phi0, '%.3f '));
fprintf('Tip speed V = %.3e mm/s\n\n', V_tip);
disp(globalSummary);
fprintf('\nGap-by-gap summary:\n');
disp(summaryTable(:, {'gap_mm', 'single_width_rel_err_percent', 'signature_l2_err', ...
    'signature_cosine_similarity', 'identified_gap_by_signature_mm', ...
    'signature_identification_correct'}));
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function [gapList, xCell, yCell] = load_stacked_comsol_curves(filePath)
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

function u = multimode_vibration(t, A, f, phi)
    u = zeros(size(t));
    for k = 1:numel(A)
        u = u + A(k) * sin(2 * pi * f(k) * t + phi(k));
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function sig = signature_from_gap(x, y, levels)
    [base, amp, ~] = normalize_waveform(y);
    yTargets = base + amp * levels;
    [widths, ~, ~] = level_widths_abs(x, y, yTargets);
    sig = normalize_signature(widths);
end

function sig = normalize_signature(widths)
    widths = widths(:);
    sig = widths / max(sum(widths, 'omitnan'), eps);
end

function [widths, xLeft, xRight] = level_widths_abs(x, y, yTargets)
    x = x(:);
    y = y(:);
    yTargets = yTargets(:);
    [~, iPeak] = max(y);
    xLeft = nan(size(yTargets));
    xRight = nan(size(yTargets));
    for k = 1:numel(yTargets)
        q = yTargets(k);
        xLeft(k) = find_crossing_left(x(1:iPeak), y(1:iPeak), q);
        xRight(k) = find_crossing_right(x(iPeak:end), y(iPeak:end), q);
    end
    widths = xRight - xLeft;
end

function xq = find_crossing_left(x, y, q)
    idx = find(y >= q, 1, 'first');
    if isempty(idx)
        xq = NaN;
    elseif idx == 1
        xq = x(1);
    else
        xq = interp1(y(idx-1:idx), x(idx-1:idx), q, 'linear');
    end
end

function xq = find_crossing_right(x, y, q)
    idx = find(y < q, 1, 'first');
    if isempty(idx)
        xq = x(end);
    elseif idx == 1
        xq = x(1);
    else
        xq = interp1(y(idx-1:idx), x(idx-1:idx), q, 'linear');
    end
end
