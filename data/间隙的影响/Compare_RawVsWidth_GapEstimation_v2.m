%% Compare raw-waveform matching vs width-signature matching for gap estimation
% Clean rewrite to avoid encoding issues
clc; clear; close all;

scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'raw_vs_width_comparison_results');
if ~isfolder(outDir), mkdir(outDir); end

set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, 'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, 'DefaultLineLineWidth', 1.15, ...
    'DefaultAxesTickDir', 'in', 'DefaultAxesBox', 'on');

%% Load data
data = readmatrix(dataFile, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
xAll = data(:, 1);
yAll = data(:, 2);
breakIdx = [find(diff(xAll) < 0); numel(xAll)];
startIdx = [1; breakIdx(1:end-1) + 1];
nGap = numel(breakIdx);
gapList = (0.2:0.2:1.4)';
xCell = cell(nGap, 1);
yCell = cell(nGap, 1);
for i = 1:nGap
    idx = startIdx(i):breakIdx(i);
    xCell{i} = xAll(idx);
    yCell{i} = yAll(idx);
end

qLevels = (0.10:0.05:0.90)';
nQ = numel(qLevels);
V_tip = 3.0e5;
fprintf('Loaded %d static calibration waveforms (%.1f - %.1f mm)\n', nGap, min(gapList), max(gapList));

%% Build libraries
% Common grid for raw waveforms
xAllCat = cell2mat(xCell);
xcFine = linspace(min(xAllCat)-mean(xAllCat), max(xAllCat)-mean(xAllCat), 2000)';

% Library A: raw normalized waveforms
Y_lib_raw = zeros(nGap, numel(xcFine));
for ig = 1:nGap
    xg = xCell{ig}(:);
    yg = yCell{ig}(:);
    base = min(yg);
    yBc = yg - base;
    [~, iPeak] = max(yBc);
    xc = xg - xg(iPeak);
    yn = yBc / max(yBc);
    Y_lib_raw(ig, :) = interp1(xc, yn, xcFine, 'linear', 'extrap');
    Y_lib_raw(ig, :) = max(0, min(1, Y_lib_raw(ig, :)));
end

% Library B: normalized width signatures
W_lib = zeros(nGap, nQ);
for ig = 1:nGap
    xg = xCell{ig}(:);
    yg = yCell{ig}(:);
    base = min(yg);
    yBc = yg - base;
    [~, iPeak] = max(yBc);
    xc = xg - xg(iPeak);
    yn = yBc / max(yBc);
    w = extract_level_widths(xc, yn, qLevels);
    w(~isfinite(w)) = 0;
    W_lib(ig, :) = w;
end
Rho_lib = W_lib ./ max(sum(W_lib, 2), eps);

%% Test sweep
vibAmps = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0];
vibFreqs = [0, 200, 400, 600, 800, 1000];
nAmp = numel(vibAmps);
nFreq = numel(vibFreqs);

allResults = [];
fprintf('Testing %d amps x %d freqs x %d gaps = %d cases...\n', nAmp, nFreq, nGap, nAmp*nFreq*nGap);

rng(42);

for ig = 1:nGap
    gTrue = gapList(ig);
    xTrue = xCell{ig}(:);
    yTrue = yCell{ig}(:);

    for ia = 1:nAmp
        A = vibAmps(ia);
        for ifr = 1:nFreq
            f = vibFreqs(ifr);

            % Generate vibrating waveform
            x0 = xTrue;
            x0Center = mean(xTrue);
            for iter = 1:4
                t = (x0 - x0Center) / V_tip;
                u = A * sin(2*pi*f*t);
                if f == 0, u = A * ones(size(t)); end
                x0 = linspace(min(xTrue)+max(u), max(xTrue)+min(u), numel(xTrue))';
            end
            t = (x0 - x0Center) / V_tip;
            if f == 0
                u = A * ones(size(t));
                du_dt = zeros(size(t));
            else
                u = A * sin(2*pi*f*t);
                du_dt = 2*pi*f*A * cos(2*pi*f*t);
            end
            xi = x0 - u;

            % Validity check
            if any(diff(xi) <= 0) || min(xi) < min(xTrue)-1e-8 || max(xi) > max(xTrue)+1e-8
                continue;
            end

            yObs = interp1(xTrue, yTrue, xi, 'linear', 'extrap');
            eta = max(abs(du_dt)) / V_tip;

            % Normalize observed waveform
            baseObs = min(yObs);
            yObsBc = yObs - baseObs;
            [~, iPeakObs] = max(yObsBc);
            xcObs = x0 - x0(iPeakObs);
            ynObs = yObsBc / max(yObsBc);
            ynObs = max(0, min(1, ynObs));

            % Method 1: Raw waveform correlation
            yObsOnGrid = interp1(xcObs, ynObs, xcFine, 'linear', 'extrap');
            yObsOnGrid = max(0, min(1, yObsOnGrid));
            yObsOnGrid = yObsOnGrid(:)';  % row vector to match Y_lib_raw

            rawCorrs = zeros(nGap, 1);
            for j = 1:nGap
                valid = isfinite(yObsOnGrid) & isfinite(Y_lib_raw(j, :));
                nv = sum(valid);
                if nv > 50
                    a = yObsOnGrid(valid);
                    b = Y_lib_raw(j, valid);
                    if std(a) > 0 && std(b) > 0
                        rawCorrs(j) = pearson_corr(a(:), b(:));
                    end
                end
            end
            [~, rawIdx] = max(rawCorrs);
            gHat_raw = gapList(rawIdx);

            % Method 1b: Raw waveform RMSE
            rawRmse = zeros(nGap, 1);
            for j = 1:nGap
                valid = isfinite(yObsOnGrid) & isfinite(Y_lib_raw(j, :));
                if sum(valid) > 50
                    rawRmse(j) = sqrt(mean((yObsOnGrid(valid) - Y_lib_raw(j, valid)).^2));
                else
                    rawRmse(j) = inf;
                end
            end
            [~, rawRmseIdx] = min(rawRmse);
            gHat_rawRmse = gapList(rawRmseIdx);

            % Method 2: Width signature matching
            wObs = extract_level_widths(xcObs, ynObs, qLevels);
            wObs(~isfinite(wObs)) = 0;
            if sum(wObs) > eps
                rhoObs = wObs(:)' / sum(wObs);
            else
                rhoObs = zeros(1, nQ);
            end
            widthDists = zeros(nGap, 1);
            for j = 1:nGap
                widthDists(j) = norm(rhoObs - Rho_lib(j, :));
            end
            [~, widthIdx] = min(widthDists);
            gHat_width = gapList(widthIdx);

            % Store results
            r.gap_true_mm = gTrue;
            r.vib_A_mm = A;
            r.vib_f_Hz = f;
            r.eta = eta;
            r.gHat_raw_corr_mm = gHat_raw;
            r.gHat_raw_rmse_mm = gHat_rawRmse;
            r.gHat_width_mm = gHat_width;
            r.err_raw_corr_mm = gHat_raw - gTrue;
            r.err_raw_rmse_mm = gHat_rawRmse - gTrue;
            r.err_width_mm = gHat_width - gTrue;
            r.raw_corr_correct = double(gHat_raw == gTrue);
            r.raw_rmse_correct = double(gHat_rawRmse == gTrue);
            r.width_correct = double(gHat_width == gTrue);
            r.best_raw_corr = max(rawCorrs);
            r.best_width_dist = min(widthDists);

            [sortedDists, sortedIdx] = sort(widthDists);
            correctRank = find(sortedIdx == ig);
            r.width_correct_rank = correctRank;
            if correctRank == 1
                r.width_separability = sortedDists(2) / max(sortedDists(1), eps) - 1;
            else
                r.width_separability = 0;
            end

            allResults = [allResults; r]; %#ok<AGROW>
        end
    end
end

T = struct2table(allResults);
writetable(T, fullfile(outDir, 'raw_vs_width_all_cases.csv'));

%% Summary statistics
fprintf('\n========== OVERALL COMPARISON ==========\n');

acc_raw_corr = mean(T.raw_corr_correct) * 100;
acc_raw_rmse = mean(T.raw_rmse_correct) * 100;
acc_width = mean(T.width_correct) * 100;

fprintf('Overall gap estimation accuracy (exact match):\n');
fprintf('  Raw waveform (correlation):  %.1f%%\n', acc_raw_corr);
fprintf('  Raw waveform (RMSE):         %.1f%%\n', acc_raw_rmse);
fprintf('  Width signature:             %.1f%%\n', acc_width);

rmse_raw_corr = sqrt(mean(T.err_raw_corr_mm .^ 2));
rmse_raw_rmse = sqrt(mean(T.err_raw_rmse_mm .^ 2));
rmse_width = sqrt(mean(T.err_width_mm .^ 2));

fprintf('\nGap estimation RMSE (mm):\n');
fprintf('  Raw waveform (correlation):  %.4f mm\n', rmse_raw_corr);
fprintf('  Raw waveform (RMSE):         %.4f mm\n', rmse_raw_rmse);
fprintf('  Width signature:             %.4f mm\n', rmse_width);

% Breakdown by vibration amplitude
fprintf('\n--- By vibration amplitude ---\n');
fprintf('%-10s %10s %10s %10s\n', 'A (mm)', 'RawCorr%%', 'RawRmse%%', 'Width%%');
for ia = 1:nAmp
    mask = abs(T.vib_A_mm - vibAmps(ia)) < 1e-10;
    if sum(mask) == 0, continue; end
    fprintf('%-10.2f %10.1f %10.1f %10.1f\n', vibAmps(ia), ...
        mean(T.raw_corr_correct(mask))*100, ...
        mean(T.raw_rmse_correct(mask))*100, ...
        mean(T.width_correct(mask))*100);
end

fprintf('\n--- By vibration frequency ---\n');
fprintf('%-10s %10s %10s %10s\n', 'f (Hz)', 'RawCorr%%', 'RawRmse%%', 'Width%%');
for ifr = 1:nFreq
    mask = abs(T.vib_f_Hz - vibFreqs(ifr)) < 1e-10;
    if sum(mask) == 0, continue; end
    fprintf('%-10d %10.1f %10.1f %10.1f\n', vibFreqs(ifr), ...
        mean(T.raw_corr_correct(mask))*100, ...
        mean(T.raw_rmse_correct(mask))*100, ...
        mean(T.width_correct(mask))*100);
end

% Save summary
summaryTable = table(...
    ["Raw waveform (correlation)"; "Raw waveform (RMSE)"; "Width signature"], ...
    [acc_raw_corr; acc_raw_rmse; acc_width], ...
    [rmse_raw_corr; rmse_raw_rmse; rmse_width], ...
    'VariableNames', {'Method', 'Accuracy_percent', 'RMSE_mm'});
writetable(summaryTable, fullfile(outDir, 'raw_vs_width_summary.csv'));

%% FIGURES
colors = lines(14);

%% Figure 1: Diagnostic comparison
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 16]);
tl1 = tiledlayout(fig1, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Representative case
repGapIdx = 4;
repA = 0.5;
repF = 600;
gTrue_rep = gapList(repGapIdx);
xTrue_rep = xCell{repGapIdx}(:);
yTrue_rep = yCell{repGapIdx}(:);

x0 = xTrue_rep;
x0Center = mean(xTrue_rep);
for iter = 1:4
    t = (x0 - x0Center) / V_tip;
    u = repA * sin(2*pi*repF*t);
    x0 = linspace(min(xTrue_rep)+max(u), max(xTrue_rep)+min(u), numel(xTrue_rep))';
end
t = (x0 - x0Center) / V_tip;
u = repA * sin(2*pi*repF*t);
xi = x0 - u;
yObs_rep = interp1(xTrue_rep, yTrue_rep, xi, 'linear', 'extrap');
eta_rep = 2*pi*repF*repA / V_tip;

baseObs_rep = min(yObs_rep);
yObsBc_rep = yObs_rep - baseObs_rep;
[~, iPeakObs_rep] = max(yObsBc_rep);
xcObs_rep = x0 - x0(iPeakObs_rep);
ynObs_rep = yObsBc_rep / max(yObsBc_rep);
ynObs_rep = max(0, min(1, ynObs_rep));

yObsOnGrid_rep = interp1(xcObs_rep, ynObs_rep, xcFine, 'linear', 'extrap');
yObsOnGrid_rep = max(0, min(1, yObsOnGrid_rep));
yObsOnGrid_rep = yObsOnGrid_rep(:)';

rawCorrs_rep = zeros(nGap, 1);
for j = 1:nGap
    valid = isfinite(yObsOnGrid_rep) & isfinite(Y_lib_raw(j, :));
    if sum(valid) > 50
        a = yObsOnGrid_rep(valid); b = Y_lib_raw(j, valid);
        if std(a) > 0 && std(b) > 0
            rawCorrs_rep(j) = pearson_corr(a(:), b(:));
        end
    end
end
[~, rawBest_rep] = max(rawCorrs_rep);

wObs_rep = extract_level_widths(xcObs_rep, ynObs_rep, qLevels);
wObs_rep(~isfinite(wObs_rep)) = 0;
if sum(wObs_rep) > eps
    rhoObs_rep = wObs_rep(:)' / sum(wObs_rep);
else
    rhoObs_rep = zeros(1, nQ);
end
widthDists_rep = zeros(nGap, 1);
for j = 1:nGap
    widthDists_rep(j) = norm(rhoObs_rep - Rho_lib(j, :));
end
[~, widthBest_rep] = min(widthDists_rep);

% (a) Observed vs library waveforms
nexttile; hold on;
for j = 1:nGap
    plot(xcFine, Y_lib_raw(j, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end
plot(xcFine, Y_lib_raw(repGapIdx, :), 'k-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Static, g=%.1f (truth)', gTrue_rep));
plot(xcObs_rep, ynObs_rep, 'r-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Observed, A=%.1f, f=%d', repA, repF));
xlabel('x - x_{peak} (mm)');
ylabel('Normalized capacitance');
title(sprintf('(a) Observed vs library, \\eta=%.4f', eta_rep));
legend('Location', 'best', 'FontSize', 6.5);

% (b) Raw correlation scores
nexttile; hold on;
bar(1:nGap, rawCorrs_rep, 0.6, 'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
set(gca, 'XTick', 1:nGap, 'XTickLabel', string(gapList));
xline(repGapIdx, 'k-', 'LineWidth', 2.5);
xline(rawBest_rep, 'r--', 'LineWidth', 2.0);
xlabel('Library gap (mm)');
ylabel('Correlation with observed');
title(sprintf('(b) Raw matching -> g=%.1f, truth g=%.1f', gapList(rawBest_rep), gTrue_rep));

% (c) Width signature comparison
nexttile; hold on;
for j = 1:nGap
    plot(qLevels, Rho_lib(j, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end
plot(qLevels, Rho_lib(repGapIdx, :), 'k-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Library, g=%.1f', gTrue_rep));
plot(qLevels, rhoObs_rep, 'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 5, ...
    'DisplayName', 'Observed');
xlabel('Normalized level q');
ylabel('Normalized width \rho(q)');
title('(c) Width signatures: library vs observed');
legend('Location', 'best', 'FontSize', 7);

% (d) Width signature distances
nexttile; hold on;
bar(1:nGap, widthDists_rep, 0.6, 'FaceColor', [0.84 0.89 0.75], 'EdgeColor', 'k', 'LineWidth', 0.5);
set(gca, 'XTick', 1:nGap, 'XTickLabel', string(gapList));
xline(repGapIdx, 'k-', 'LineWidth', 2.5);
xline(widthBest_rep, 'r--', 'LineWidth', 2.0);
xlabel('Library gap (mm)');
ylabel('||\rho_{obs} - \rho_{lib}||');
title(sprintf('(d) Width matching -> g=%.1f, truth g=%.1f', gapList(widthBest_rep), gTrue_rep));

% (e) Raw correlation score vs eta
nexttile; hold on;
idx_correct = T.raw_corr_correct == 1;
idx_wrong = T.raw_corr_correct == 0;
scatter(T.eta(idx_correct), T.best_raw_corr(idx_correct), ...
    10, [0.10 0.45 0.75], 'o', 'DisplayName', 'Raw correct');
scatter(T.eta(idx_wrong), T.best_raw_corr(idx_wrong), ...
    18, [0.85 0.35 0.10], 'x', 'LineWidth', 1.5, 'DisplayName', 'Raw WRONG');
xlabel('Distortion \eta');
ylabel('Best raw correlation');
title('(e) Raw correlation score vs \eta');
legend('Location', 'best', 'FontSize', 7);

% (f) Width advantage heatmap
nexttile;
accWidth2d = zeros(nAmp, nFreq);
accRaw2d = zeros(nAmp, nFreq);
for ia = 1:nAmp
    for ifr = 1:nFreq
        mask = abs(T.vib_A_mm - vibAmps(ia)) < 1e-10 & abs(T.vib_f_Hz - vibFreqs(ifr)) < 1e-10;
        if sum(mask) > 0
            accWidth2d(ia, ifr) = mean(T.width_correct(mask)) * 100;
            accRaw2d(ia, ifr) = mean(T.raw_corr_correct(mask)) * 100;
        end
    end
end
diff2d = accWidth2d - accRaw2d;
imagesc(vibFreqs, vibAmps, diff2d);
colormap(gca, jet);
caxis([-100 100]);
colorbar;
xlabel('Frequency (Hz)');
ylabel('Amplitude (mm)');
title('(f) Width advantage: acc_{width} - acc_{raw} (%)');

% (g) Accuracy vs eta
nexttile; hold on;
etaBins = [0, 0.001, 0.003, 0.01, 0.03, 0.1];
etaBinCenters = (etaBins(1:end-1) + etaBins(2:end)) / 2;
accWidthEta = zeros(numel(etaBins)-1, 1);
accRawEta = zeros(numel(etaBins)-1, 1);
for ib = 1:numel(etaBins)-1
    mask = T.eta >= etaBins(ib) & T.eta < etaBins(ib+1);
    if sum(mask) > 0
        accWidthEta(ib) = mean(T.width_correct(mask)) * 100;
        accRawEta(ib) = mean(T.raw_corr_correct(mask)) * 100;
    end
end
plot(etaBinCenters, accWidthEta, 'o-', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.5, ...
    'DisplayName', 'Width signature');
plot(etaBinCenters, accRawEta, 's-', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.5, ...
    'DisplayName', 'Raw waveform (corr)');
xlabel('Distortion \eta = max|du/dt| / V');
ylabel('Gap estimation accuracy (%)');
title('(g) Accuracy vs distortion strength');
legend('Location', 'southwest', 'FontSize', 8);

% (h) RMSE vs eta
nexttile; hold on;
rmseWidthEta = zeros(numel(etaBins)-1, 1);
rmseRawEta = zeros(numel(etaBins)-1, 1);
for ib = 1:numel(etaBins)-1
    mask = T.eta >= etaBins(ib) & T.eta < etaBins(ib+1);
    if sum(mask) > 0
        rmseWidthEta(ib) = sqrt(mean(T.err_width_mm(mask).^2));
        rmseRawEta(ib) = sqrt(mean(T.err_raw_corr_mm(mask).^2));
    end
end
plot(etaBinCenters, rmseWidthEta, 'o-', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.5, ...
    'DisplayName', 'Width signature');
plot(etaBinCenters, rmseRawEta, 's-', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.5, ...
    'DisplayName', 'Raw waveform (corr)');
xlabel('Distortion \eta = max|du/dt| / V');
ylabel('Gap estimation RMSE (mm)');
title('(h) RMSE vs distortion strength');
legend('Location', 'northwest', 'FontSize', 8);

% (i) Separability margin histogram
nexttile;
validSep = T.width_separability(T.width_separability > 0);
histogram(validSep, 30, 'FaceColor', [0.10 0.45 0.75], 'EdgeColor', 'k', 'LineWidth', 0.5);
xline(0.1, 'r--', 'LineWidth', 1.0);
xlabel('Separability margin (d_2/d_1 - 1)');
ylabel('Count');
title(sprintf('(i) Width separability margin, mean=%.3f', mean(validSep)));

title(tl1, sprintf('Raw waveform vs width signature for gap estimation (g=%.1f mm, A=%.1f mm, f=%d Hz)', ...
    gTrue_rep, repA, repF), 'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_raw_vs_width_diagnostic.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_raw_vs_width_diagnostic.pdf'), 'ContentType', 'vector');

%% Figure 2: Comprehensive summary
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 12]);
tl2 = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Overall accuracy bar
nexttile;
accVals = [acc_raw_corr, acc_raw_rmse, acc_width];
bar(categorical({'Raw corr', 'Raw RMSE', 'Width sig'}), accVals, 0.6, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Accuracy (%)');
title(sprintf('(a) Overall accuracy (N=%d)', height(T)));
ylim([0 105]);
for i = 1:3
    text(i, accVals(i)+2, sprintf('%.1f%%', accVals(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end

% (b) Accuracy vs amplitude
nexttile; hold on;
accWidthAmp = zeros(nAmp, 1);
accRawAmp = zeros(nAmp, 1);
for ia = 1:nAmp
    mask = abs(T.vib_A_mm - vibAmps(ia)) < 1e-10;
    if sum(mask) > 0
        accWidthAmp(ia) = mean(T.width_correct(mask)) * 100;
        accRawAmp(ia) = mean(T.raw_corr_correct(mask)) * 100;
    end
end
plot(vibAmps, accWidthAmp, 'o-', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.5, ...
    'DisplayName', 'Width signature');
plot(vibAmps, accRawAmp, 's-', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.5, ...
    'DisplayName', 'Raw waveform (corr)');
xlabel('Vibration amplitude A (mm)');
ylabel('Accuracy (%)');
title('(b) Accuracy vs amplitude');
legend('Location', 'southwest', 'FontSize', 8);
ylim([0 105]);

% (c) Accuracy vs frequency
nexttile; hold on;
accWidthFrq = zeros(nFreq, 1);
accRawFrq = zeros(nFreq, 1);
for ifr = 1:nFreq
    mask = abs(T.vib_f_Hz - vibFreqs(ifr)) < 1e-10;
    if sum(mask) > 0
        accWidthFrq(ifr) = mean(T.width_correct(mask)) * 100;
        accRawFrq(ifr) = mean(T.raw_corr_correct(mask)) * 100;
    end
end
plot(vibFreqs, accWidthFrq, 'o-', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.5, ...
    'DisplayName', 'Width signature');
plot(vibFreqs, accRawFrq, 's-', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.5, ...
    'DisplayName', 'Raw waveform (corr)');
xlabel('Vibration frequency f (Hz)');
ylabel('Accuracy (%)');
title('(c) Accuracy vs frequency');
legend('Location', 'southwest', 'FontSize', 8);
ylim([0 105]);

% (d) RMSE vs amplitude
nexttile; hold on;
rmseWidthAmp = zeros(nAmp, 1);
rmseRawAmp = zeros(nAmp, 1);
for ia = 1:nAmp
    mask = abs(T.vib_A_mm - vibAmps(ia)) < 1e-10;
    if sum(mask) > 0
        rmseWidthAmp(ia) = sqrt(mean(T.err_width_mm(mask).^2));
        rmseRawAmp(ia) = sqrt(mean(T.err_raw_corr_mm(mask).^2));
    end
end
plot(vibAmps, rmseWidthAmp, 'o-', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.5, ...
    'DisplayName', 'Width signature');
plot(vibAmps, rmseRawAmp, 's-', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.5, ...
    'DisplayName', 'Raw waveform (corr)');
xlabel('Vibration amplitude A (mm)');
ylabel('Gap RMSE (mm)');
title('(d) RMSE vs amplitude');
legend('Location', 'northwest', 'FontSize', 8);
yline(0.2, 'k--', 'LineWidth', 0.8);

% (e) Raw confusion matrix at high eta
nexttile;
maskHi = T.eta > 0.003;
confRaw = zeros(nGap, nGap);
for i = 1:nGap
    for j = 1:nGap
        m = maskHi & abs(T.gap_true_mm - gapList(i)) < 1e-10 & abs(T.gHat_raw_corr_mm - gapList(j)) < 1e-10;
        confRaw(i, j) = sum(m);
    end
end
confRaw = confRaw ./ max(sum(confRaw, 2), 1);
imagesc(gapList, gapList, confRaw);
colormap(gca, flipud(gray));
caxis([0 1]); colorbar;
xlabel('Estimated gap (mm)');
ylabel('True gap (mm)');
title('(e) Raw conf matrix, \eta>0.003');

% (f) Width confusion matrix at high eta
nexttile;
confWidth = zeros(nGap, nGap);
for i = 1:nGap
    for j = 1:nGap
        m = maskHi & abs(T.gap_true_mm - gapList(i)) < 1e-10 & abs(T.gHat_width_mm - gapList(j)) < 1e-10;
        confWidth(i, j) = sum(m);
    end
end
confWidth = confWidth ./ max(sum(confWidth, 2), 1);
imagesc(gapList, gapList, confWidth);
colormap(gca, flipud(gray));
caxis([0 1]); colorbar;
xlabel('Estimated gap (mm)');
ylabel('True gap (mm)');
title('(f) Width conf matrix, \eta>0.003');

title(tl2, 'Raw waveform vs width signature: summary', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_raw_vs_width_summary.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_raw_vs_width_summary.pdf'), 'ContentType', 'vector');

%% Figure 3: Why raw fails
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 8]);
tl3 = tiledlayout(fig3, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

failMask = T.raw_corr_correct == 0 & T.width_correct == 1 & T.eta > 0.003;
if any(failMask)
    failIdx = find(failMask, 1, 'first');
    fCase = T(failIdx, :);

    igF = find(abs(gapList - fCase.gap_true_mm) < 1e-10, 1);
    xTrueF = xCell{igF}(:);
    yTrueF = yCell{igF}(:);
    x0CenterF = mean(xTrueF);
    x0F = xTrueF;
    for iter = 1:4
        tF = (x0F - x0CenterF) / V_tip;
        uF = fCase.vib_A_mm * sin(2*pi*fCase.vib_f_Hz*tF);
        x0F = linspace(min(xTrueF)+max(uF), max(xTrueF)+min(uF), numel(xTrueF))';
    end
    tF = (x0F - x0CenterF) / V_tip;
    uF = fCase.vib_A_mm * sin(2*pi*fCase.vib_f_Hz*tF);
    xiF = x0F - uF;
    yObsF = interp1(xTrueF, yTrueF, xiF, 'linear', 'extrap');

    baseOF = min(yObsF); yObsFbc = yObsF - baseOF;
    [~, iPeakOF] = max(yObsFbc);
    xcOF = x0F - x0F(iPeakOF); ynOF = yObsFbc / max(yObsFbc);
    ynOF = max(0, min(1, ynOF));

    % True static
    baseTF = min(yTrueF); yTFbc = yTrueF - baseTF;
    [~, iPeakTF] = max(yTFbc);
    xcTF = xTrueF - xTrueF(iPeakTF);
    ynTF = yTFbc / max(yTFbc);

    % Wrong match
    igW = find(abs(gapList - fCase.gHat_raw_corr_mm) < 1e-10, 1);
    xW = xCell{igW}(:); yW = yCell{igW}(:);
    baseW = min(yW); yWbc = yW - baseW;
    [~, iPeakW] = max(yWbc);
    xcW = xW - xW(iPeakW);
    ynW = yWbc / max(yWbc);

    % (a) Waveforms
    nexttile; hold on;
    plot(xcTF, ynTF, 'k-', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Static g=%.1f (truth)', fCase.gap_true_mm));
    plot(xcW, ynW, '--', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Static g=%.1f (raw match)', fCase.gHat_raw_corr_mm));
    plot(xcOF, ynOF, 'r-', 'LineWidth', 1.0, 'DisplayName', 'Observed (vibrating)');
    xlabel('x - x_{peak} (mm)');
    ylabel('Normalized capacitance');
    title(sprintf('(a) Waveforms, \\eta=%.4f', fCase.eta));
    legend('Location', 'best', 'FontSize', 6.5);

    % (b) Width signatures
    wOF = extract_level_widths(xcOF, ynOF, qLevels);
    wOF(~isfinite(wOF)) = 0;
    rhoOF = wOF(:)' / max(sum(wOF), eps);

    nexttile; hold on;
    plot(qLevels, Rho_lib(igF, :), 'k-', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Library g=%.1f (truth)', fCase.gap_true_mm));
    plot(qLevels, Rho_lib(igW, :), '--', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Library g=%.1f (raw match)', fCase.gHat_raw_corr_mm));
    plot(qLevels, rhoOF, 'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 5, ...
        'DisplayName', 'Observed');
    xlabel('Normalized level q');
    ylabel('Normalized width \rho(q)');
    title('(b) Width signatures');
    legend('Location', 'best', 'FontSize', 6.5);

    % (c) Distance comparison
    nexttile; hold on;
    yObsOnGridF = interp1(xcOF, ynOF, xcFine, 'linear', 'extrap');
    yObsOnGridF = max(0, min(1, yObsOnGridF));
    yObsOnGridF = yObsOnGridF(:)';
    rawDists = zeros(nGap, 1);
    widthDistsF = zeros(nGap, 1);
    for j = 1:nGap
        valid = isfinite(yObsOnGridF) & isfinite(Y_lib_raw(j, :));
        if sum(valid) > 50
            rawDists(j) = 1 - pearson_corr(yObsOnGridF(valid)', Y_lib_raw(j, valid)');
        else
            rawDists(j) = 1;
        end
        widthDistsF(j) = norm(rhoOF - Rho_lib(j, :));
    end
    rawDistsNorm = rawDists / max(rawDists);
    widthDistsNorm = widthDistsF / max(widthDistsF);

    plot(gapList, rawDistsNorm, 's-', 'Color', [0.85 0.35 0.10], ...
        'MarkerFaceColor', [0.85 0.35 0.10], 'LineWidth', 1.2, ...
        'DisplayName', 'Raw waveform distance');
    plot(gapList, widthDistsNorm, 'o-', 'Color', [0.10 0.45 0.75], ...
        'MarkerFaceColor', [0.10 0.45 0.75], 'LineWidth', 1.2, ...
        'DisplayName', 'Width signature distance');
    xline(fCase.gap_true_mm, 'k-', 'LineWidth', 2.0);
    xlabel('Library gap (mm)');
    ylabel('Normalized distance');
    title('(c) Distance to library entries');
    legend('Location', 'best', 'FontSize', 7);
end

title(tl3, 'Why raw matching fails while width signature succeeds', 'FontWeight', 'normal');
exportgraphics(fig3, fullfile(outDir, 'fig3_why_raw_fails.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_why_raw_fails.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n========================================\n');
fprintf('RAW WAVEFORM vs WIDTH SIGNATURE COMPARISON\n');
fprintf('========================================\n');
fprintf('Test cases: %d\n', height(T));
fprintf('Vibration range: A = %.2f-%.1f mm, f = %d-%d Hz\n', min(vibAmps), max(vibAmps), min(vibFreqs), max(vibFreqs));
fprintf('eta range: %.5f - %.4f\n\n', min(T.eta), max(T.eta));
fprintf('OVERALL ACCURACY:\n');
fprintf('  Raw waveform (correlation):  %5.1f%%  (RMSE = %.3f mm)\n', acc_raw_corr, rmse_raw_corr);
fprintf('  Raw waveform (RMSE):         %5.1f%%  (RMSE = %.3f mm)\n', acc_raw_rmse, rmse_raw_rmse);
fprintf('  Width signature:             %5.1f%%  (RMSE = %.3f mm)\n', acc_width, rmse_width);
fprintf('\nKEY INSIGHT:\n');
fprintf('  Vibration shifts the waveform laterally in x0.\n');
fprintf('  Raw point-by-point comparison sees EVERY point shifted.\n');
fprintf('  Width signature normalizes out the common shift, preserving\n');
fprintf('  only the shape information that encodes gap.\n');
fprintf('\nResults saved to: %s\n', outDir);

%% Local functions

function w = extract_level_widths(xCenter, yNorm, qLevels)
    xCenter = xCenter(:);
    yNorm = yNorm(:);
    qLevels = qLevels(:);
    [~, iPeak] = max(yNorm);
    nQ = numel(qLevels);
    w = nan(nQ, 1);
    for k = 1:nQ
        q = qLevels(k);
        yRise = yNorm(1:iPeak);
        xRise = xCenter(1:iPeak);
        yFall = yNorm(iPeak:end);
        xFall = xCenter(iPeak:end);
        if q >= min(yRise) && q <= max(yRise) && q >= min(yFall) && q <= max(yFall)
            [yRiseU, idxU] = unique(yRise, 'stable');
            [yFallU, idxFU] = unique(yFall, 'stable');
            xL = interp1(yRiseU, xRise(idxU), q, 'linear');
            xR = interp1(flipud(yFallU), flipud(xFall(idxFU)), q, 'linear');
            w(k) = xR - xL;
        end
    end
end

function c = pearson_corr(x, y)
    x = x(:); y = y(:);
    xm = x - mean(x);
    ym = y - mean(y);
    denom = sqrt(xm' * xm) * sqrt(ym' * ym);
    if denom > 0
        c = (xm' * ym) / denom;
    else
        c = 0;
    end
end
