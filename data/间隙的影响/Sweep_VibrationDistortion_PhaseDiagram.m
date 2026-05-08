%% Sweep vibration parameters and visualize waveform distortion
% This script helps explain when a vibrating waveform in t-V / x0-V looks
% like a simple translation, and when it starts to show real shape
% distortion. The true relative-coordinate waveform xi-V is also checked.

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'vibration_distortion_phase_diagram_results');
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

%% Load one baseline waveform
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
[~, refIdx] = min(abs(gapList - 0.8));
xRef = xCell{refIdx}(:);
yRef = yCell{refIdx}(:);
x0Center = mean(xRef);

%% Sweep settings
AList = 0:0.2:1.8;               % mm
fList = 100:100:3000;            % Hz
VtipFixed = 3.0e5;               % mm/s
phiFixed = pi / 4;               % rad

VtipList = [1.5e5, 2.0e5, 3.0e5, 5.0e5];   % mm/s
AFixed = 0.8;                    % mm

% User multimode case for reference
userCase.A = [0.5, 0.4];
userCase.f = [500, 1300];
userCase.phi = [pi/4, -pi/3];
userCase.Vtip = VtipFixed;

%% Sweep 1: amplitude-frequency map at fixed tip speed
nA = numel(AList);
nF = numel(fList);
rawRmseMap = nan(nA, nF);
shiftRmseMap = nan(nA, nF);
peakShiftMap = nan(nA, nF);
etaMap = nan(nA, nF);
bestShiftMap = nan(nA, nF);

rowCounter = 0;
rowsSweep1(nA * nF) = struct('A_mm', 0, 'f_Hz', 0, 'Vtip_mm_per_s', 0, ...
    'eta', 0, 'rmse_x0', 0, 'rmse_after_best_shift', 0, ...
    'best_shift_mm', 0, 'peak_shift_mm', 0);

for ia = 1:nA
    for jf = 1:nF
        metrics = analyze_single_mode_case(xRef, yRef, x0Center, ...
            AList(ia), fList(jf), phiFixed, VtipFixed);
        rawRmseMap(ia, jf) = metrics.rmseX0;
        shiftRmseMap(ia, jf) = metrics.rmseAfterBestShift;
        peakShiftMap(ia, jf) = metrics.peakShiftX0;
        etaMap(ia, jf) = metrics.eta;
        bestShiftMap(ia, jf) = metrics.bestShift;

        rowCounter = rowCounter + 1;
        rowsSweep1(rowCounter).A_mm = AList(ia);
        rowsSweep1(rowCounter).f_Hz = fList(jf);
        rowsSweep1(rowCounter).Vtip_mm_per_s = VtipFixed;
        rowsSweep1(rowCounter).eta = metrics.eta;
        rowsSweep1(rowCounter).rmse_x0 = metrics.rmseX0;
        rowsSweep1(rowCounter).rmse_after_best_shift = metrics.rmseAfterBestShift;
        rowsSweep1(rowCounter).best_shift_mm = metrics.bestShift;
        rowsSweep1(rowCounter).peak_shift_mm = metrics.peakShiftX0;
    end
end

sweep1Table = struct2table(rowsSweep1);
writetable(sweep1Table, fullfile(outDir, 'sweep1_A_f_fixedVtip.csv'));

%% Sweep 2: frequency-tip-speed map at fixed amplitude
nV = numel(VtipList);
rawRmseFV = nan(nV, nF);
shiftRmseFV = nan(nV, nF);
etaFV = nan(nV, nF);

rowCounter = 0;
rowsSweep2(nV * nF) = struct('A_mm', 0, 'f_Hz', 0, 'Vtip_mm_per_s', 0, ...
    'eta', 0, 'rmse_x0', 0, 'rmse_after_best_shift', 0);

for iv = 1:nV
    for jf = 1:nF
        metrics = analyze_single_mode_case(xRef, yRef, x0Center, ...
            AFixed, fList(jf), phiFixed, VtipList(iv));
        rawRmseFV(iv, jf) = metrics.rmseX0;
        shiftRmseFV(iv, jf) = metrics.rmseAfterBestShift;
        etaFV(iv, jf) = metrics.eta;

        rowCounter = rowCounter + 1;
        rowsSweep2(rowCounter).A_mm = AFixed;
        rowsSweep2(rowCounter).f_Hz = fList(jf);
        rowsSweep2(rowCounter).Vtip_mm_per_s = VtipList(iv);
        rowsSweep2(rowCounter).eta = metrics.eta;
        rowsSweep2(rowCounter).rmse_x0 = metrics.rmseX0;
        rowsSweep2(rowCounter).rmse_after_best_shift = metrics.rmseAfterBestShift;
    end
end

sweep2Table = struct2table(rowsSweep2);
writetable(sweep2Table, fullfile(outDir, 'sweep2_f_Vtip_fixedA.csv'));

%% User multimode reference case
userMetrics = analyze_multimode_case(xRef, yRef, x0Center, ...
    userCase.A, userCase.f, userCase.phi, userCase.Vtip);
userSummary = struct2table(userMetrics, 'AsArray', true);
writetable(userSummary, fullfile(outDir, 'user_multimode_reference.csv'));

%% Representative cases for intuition
targets = [0.01, 0.05, 0.10];
repIdx = zeros(size(targets));
for k = 1:numel(targets)
    [~, repIdx(k)] = min(abs(sweep1Table.eta - targets(k)));
end
repCases = sweep1Table(repIdx, :);

repData = cell(numel(targets), 1);
for k = 1:numel(targets)
    repData{k} = analyze_single_mode_case(xRef, yRef, x0Center, ...
        repCases.A_mm(k), repCases.f_Hz(k), phiFixed, repCases.Vtip_mm_per_s(k));
end

%% Global summary
globalSummary = table( ...
    min(sweep1Table.eta), max(sweep1Table.eta), ...
    mean(sweep1Table.rmse_x0), mean(sweep1Table.rmse_after_best_shift), ...
    userMetrics.eta, userMetrics.rmseX0, userMetrics.rmseAfterBestShift, ...
    'VariableNames', {'eta_min_sweep1', 'eta_max_sweep1', ...
    'mean_rmse_x0_sweep1', 'mean_rmse_after_best_shift_sweep1', ...
    'user_eta', 'user_rmse_x0', 'user_rmse_after_best_shift'});
writetable(globalSummary, fullfile(outDir, 'global_summary.csv'));

%% Figure 1: A-f heatmaps at fixed tip speed
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 14]);
tl1 = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(fList, AList, rawRmseMap);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Amplitude A (mm)');
title('(a) Raw RMSE in x_0-V');
cb = colorbar; cb.Label.String = 'RMSE (pF)';

nexttile;
imagesc(fList, AList, shiftRmseMap);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Amplitude A (mm)');
title('(b) Residual RMSE after best shift');
cb = colorbar; cb.Label.String = 'RMSE (pF)';

nexttile;
imagesc(fList, AList, etaMap);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Amplitude A (mm)');
title('(c) Distortion index \eta');
cb = colorbar; cb.Label.String = '\eta = max|du/dt| / V_{tip}';

nexttile;
imagesc(fList, AList, bestShiftMap);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Amplitude A (mm)');
title('(d) Best horizontal shift');
cb = colorbar; cb.Label.String = 'Shift (mm)';

title(tl1, sprintf('Amplitude-frequency sweep at V_{tip} = %.1e mm/s', VtipFixed), ...
    'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_A_f_heatmaps.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_A_f_heatmaps.pdf'), 'ContentType', 'vector');

%% Figure 2: frequency-tip-speed heatmaps at fixed amplitude
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tl2 = tiledlayout(fig2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(fList, VtipList / 1e5, rawRmseFV);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Tip speed V_{tip} (10^5 mm/s)');
title('(a) Raw RMSE in x_0-V');
cb = colorbar; cb.Label.String = 'RMSE (pF)';

nexttile;
imagesc(fList, VtipList / 1e5, shiftRmseFV);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Tip speed V_{tip} (10^5 mm/s)');
title('(b) Residual RMSE after best shift');
cb = colorbar; cb.Label.String = 'RMSE (pF)';

nexttile;
imagesc(fList, VtipList / 1e5, etaFV);
set(gca, 'YDir', 'normal');
xlabel('Frequency f (Hz)');
ylabel('Tip speed V_{tip} (10^5 mm/s)');
title('(c) Distortion index \eta');
cb = colorbar; cb.Label.String = '\eta';

title(tl2, sprintf('Frequency-tip-speed sweep at A = %.1f mm', AFixed), ...
    'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_f_Vtip_heatmaps.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_f_Vtip_heatmaps.pdf'), 'ContentType', 'vector');

%% Figure 3: collapse against distortion index
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 8.5]);
tl3 = tiledlayout(fig3, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
scatter(sweep1Table.eta, sweep1Table.rmse_x0, 18, [0.1 0.45 0.75], 'filled', ...
    'DisplayName', 'Raw RMSE');
scatter(sweep1Table.eta, sweep1Table.rmse_after_best_shift, 18, [0.85 0.35 0.10], 'filled', ...
    'DisplayName', 'After best shift');
scatter(userMetrics.eta, userMetrics.rmseX0, 54, 'kp', 'filled', 'DisplayName', 'User multimode raw');
scatter(userMetrics.eta, userMetrics.rmseAfterBestShift, 54, 'kd', 'filled', 'DisplayName', 'User multimode aligned');
xlabel('Distortion index \eta');
ylabel('RMSE (pF)');
title('(a) RMSE versus distortion index');
legend('Location', 'northwest', 'FontSize', 7.0);

nexttile; hold on;
scatter(sweep1Table.eta, abs(sweep1Table.best_shift_mm - sweep1Table.peak_shift_mm), 18, ...
    [0.20 0.60 0.40], 'filled');
xlabel('Distortion index \eta');
ylabel('|Best shift - peak shift| (mm)');
title('(b) Peak shift versus best-shift mismatch');

exportgraphics(fig3, fullfile(outDir, 'fig3_distortion_index_collapse.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_distortion_index_collapse.pdf'), 'ContentType', 'vector');

%% Figure 4: representative waveform examples
fig4 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 18]);
tl4 = tiledlayout(fig4, numel(repData), 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(repData)
    rep = repData{k};
    labelBase = char('a' + (k-1) * 3);

    nexttile; hold on;
    plot(rep.t * 1e6, rep.yRefOnGrid, 'k-', 'DisplayName', 'No vibration');
    plot(rep.t * 1e6, rep.yObs, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Vibrating');
    xlabel('Time t (\mus)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) t-V, \\eta = %.3f', labelBase, rep.eta));
    legend('Location', 'best', 'FontSize', 7.0);

    nexttile; hold on;
    plot(rep.x0, rep.yRefOnGrid, 'k-', 'DisplayName', 'No vibration');
    plot(rep.x0, rep.yObs, '-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Vibrating');
    plot(rep.x0, rep.yBestShiftAligned, '--', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Best-shift aligned');
    xlabel('Rigid coordinate x_0 (mm)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) x_0-V, A = %.1f mm, f = %d Hz', labelBase + 1, rep.A, rep.f));
    legend('Location', 'best', 'FontSize', 7.0);

    nexttile; hold on;
    plot(xRef, yRef, 'k-', 'DisplayName', 'No vibration');
    plot(rep.xiSorted, rep.yXi, '-', 'Color', [0.20 0.60 0.40], 'DisplayName', 'In \xi');
    xlabel('Relative coordinate \xi (mm)');
    ylabel('Capacitance V (pF)');
    title(sprintf('(%c) \\xi-V, residual = %.2e', labelBase + 2, rep.rmseXi));
    legend('Location', 'best', 'FontSize', 7.0);
end

title(tl4, 'Representative cases from weak to strong distortion', 'FontWeight', 'normal');
exportgraphics(fig4, fullfile(outDir, 'fig4_representative_cases.png'), 'Resolution', 300);
exportgraphics(fig4, fullfile(outDir, 'fig4_representative_cases.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Vibration distortion phase-diagram analysis ===\n');
fprintf('Reference gap: %.1f mm\n', gapList(refIdx));
fprintf('Sweep 1: A in [%.1f, %.1f] mm, f in [%d, %d] Hz, Vtip = %.2e mm/s\n', ...
    min(AList), max(AList), min(fList), max(fList), VtipFixed);
fprintf('Sweep 2: A = %.1f mm, f in [%d, %d] Hz, Vtip in [%s] mm/s\n\n', ...
    AFixed, min(fList), max(fList), num2str(VtipList, '%.1e '));
disp(globalSummary);
fprintf('\nRepresentative cases:\n');
disp(repCases);
fprintf('\nUser multimode reference:\n');
disp(userSummary);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
function metrics = analyze_single_mode_case(xRef, yRef, x0Center, A, f, phi, Vtip)
    xMargin = abs(A) + 1e-12;
    x0 = linspace(min(xRef) + xMargin, max(xRef) - xMargin, numel(xRef))';
    t = (x0 - x0Center) / Vtip;
    u = A * sin(2 * pi * f * t + phi);
    xi = x0 - u;

    if any(diff(xi) <= 0)
        error('Non-monotone xi(x0) encountered for A = %.3f, f = %.1f, Vtip = %.3e.', A, f, Vtip);
    end

    yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    yRefOnGrid = interp1(xRef, yRef, x0, 'pchip', NaN);
    [bestShift, rmseAfterBestShift, yBestShiftAligned] = best_shift_alignment(x0, yObs, xRef, yRef, max(abs(u)) + 0.2);

    [~, iPeakObs] = max(yObs);
    [~, iPeakRef] = max(yRefOnGrid);
    peakShiftX0 = x0(iPeakObs) - x0(iPeakRef);

    [xiSorted, idx] = sort(xi, 'ascend');
    yXi = yObs(idx);
    yRefOnXi = interp1(xRef, yRef, xiSorted, 'pchip', NaN);
    rmseXi = sqrt(mean((yXi - yRefOnXi).^2, 'omitnan'));

    du_dt = 2 * pi * f * A * cos(2 * pi * f * t + phi);
    eta = max(abs(du_dt)) / Vtip;
    rmseX0 = sqrt(mean((yObs - yRefOnGrid).^2, 'omitnan'));

    metrics = struct( ...
        'A', A, 'f', f, 'phi', phi, 'Vtip', Vtip, ...
        'x0', x0, 't', t, 'u', u, 'xi', xi, ...
        'yObs', yObs, 'yXi', yXi, 'xiSorted', xiSorted, ...
        'yRefOnGrid', yRefOnGrid, 'yBestShiftAligned', yBestShiftAligned, ...
        'bestShift', bestShift, 'peakShiftX0', peakShiftX0, ...
        'rmseX0', rmseX0, 'rmseAfterBestShift', rmseAfterBestShift, ...
        'rmseXi', rmseXi, 'eta', eta);
end

function metrics = analyze_multimode_case(xRef, yRef, x0Center, A, f, phi, Vtip)
    xMargin = sum(abs(A)) + 1e-12;
    x0 = linspace(min(xRef) + xMargin, max(xRef) - xMargin, numel(xRef))';
    t = (x0 - x0Center) / Vtip;
    u = zeros(size(t));
    du_dt = zeros(size(t));
    for k = 1:numel(A)
        u = u + A(k) * sin(2 * pi * f(k) * t + phi(k));
        du_dt = du_dt + 2 * pi * f(k) * A(k) * cos(2 * pi * f(k) * t + phi(k));
    end
    xi = x0 - u;

    yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    yRefOnGrid = interp1(xRef, yRef, x0, 'pchip', NaN);
    [bestShift, rmseAfterBestShift] = best_shift_alignment(x0, yObs, xRef, yRef, max(abs(u)) + 0.2);

    [~, iPeakObs] = max(yObs);
    [~, iPeakRef] = max(yRefOnGrid);
    peakShiftX0 = x0(iPeakObs) - x0(iPeakRef);

    [xiSorted, idx] = sort(xi, 'ascend');
    yXi = yObs(idx);
    yRefOnXi = interp1(xRef, yRef, xiSorted, 'pchip', NaN);
    rmseXi = sqrt(mean((yXi - yRefOnXi).^2, 'omitnan'));

    metrics = struct( ...
        'A_desc', mat2str(A), 'f_desc', mat2str(f), 'phi_desc', mat2str(phi), ...
        'Vtip', Vtip, 'eta', max(abs(du_dt)) / Vtip, ...
        'rmseX0', sqrt(mean((yObs - yRefOnGrid).^2, 'omitnan')), ...
        'rmseAfterBestShift', rmseAfterBestShift, ...
        'bestShift', bestShift, 'peakShiftX0', peakShiftX0, ...
        'rmseXi', rmseXi, 'maxAbsU', max(abs(u)), ...
        'peak_shift_t_us', (peakShiftX0 / Vtip) * 1e6);
end

function [bestShift, bestRmse, yAligned] = best_shift_alignment(x0, yObs, xRef, yRef, searchRadius)
    objective = @(delta) shift_rmse(delta, x0, yObs, xRef, yRef);
    if searchRadius <= 0
        searchRadius = 1;
    end
    bestShift = fminbnd(objective, -searchRadius, searchRadius);
    [bestRmse, yAligned] = shift_rmse(bestShift, x0, yObs, xRef, yRef);
end

function [rmse, yShifted] = shift_rmse(delta, x0, yObs, xRef, yRef)
    yShifted = interp1(xRef, yRef, x0 - delta, 'pchip', NaN);
    valid = isfinite(yShifted) & isfinite(yObs);
    rmse = sqrt(mean((yObs(valid) - yShifted(valid)).^2));
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
