%% Compare raw waveform pointwise matching vs width-signature matching under vibration
% Goal:
%   Show, on the same synthetic vibrating waveform, why direct pointwise
%   matching in x0-domain is sensitive to vibration-induced coordinate shift,
%   while normalized width signatures remain much more stable.

clc; clear; close all;

scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 9, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% Load static gap library
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
qLevels = (0.10:0.05:0.90)';
V_tip = 3.0e5;  % mm/s

curve = struct();
for ig = 1:nGap
    curve(ig).gap = gapList(ig);
    curve(ig).x = xCell{ig}(:);
    curve(ig).y = yCell{ig}(:);
    [curve(ig).base, curve(ig).amp, curve(ig).yNorm] = normalize_waveform(curve(ig).y);
    [~, iPeak] = max(curve(ig).yNorm);
    curve(ig).xCenter = curve(ig).x - curve(ig).x(iPeak);
    [w, ~, ~] = extract_level_widths(curve(ig).xCenter, curve(ig).yNorm, qLevels);
    curve(ig).rho = w(:)' / max(sum(w), eps);
end

%% Vibration cases
vibCases = struct( ...
    'name', {'No vibration', 'Weak', 'Moderate', 'Strong', 'Very strong'}, ...
    'A_mm', {0, 0.10, 0.30, 0.60, 1.00}, ...
    'f_Hz', {0, 300, 500, 800, 1200});
nVib = numel(vibCases);

%% 1) Global comparison across all gaps
rows = [];
for ig = 1:nGap
    for iv = 1:nVib
        synth = synthesize_vibrating_waveform(curve(ig).x, curve(ig).y, vibCases(iv), V_tip);
        if isempty(synth)
            continue;
        end

        % Method A: raw waveform pointwise matching on the observed x0-grid
        [~, ~, yObsNorm] = normalize_waveform(synth.yObs);
        rawDist = nan(nGap, 1);
        for j = 1:nGap
            yLibOnObs = interp1(curve(j).x, curve(j).y, synth.x0, 'pchip', NaN);
            [~, ~, yLibNorm] = normalize_waveform(yLibOnObs);
            valid = isfinite(yObsNorm) & isfinite(yLibNorm);
            rawDist(j) = sqrt(mean((yObsNorm(valid) - yLibNorm(valid)).^2));
        end
        [~, idxRaw] = min(rawDist);
        gHatRaw = gapList(idxRaw);

        % Method B: normalized width signature matching
        [~, ~, yObsNorm2] = normalize_waveform(synth.yObs);
        [~, iPeakObs] = max(yObsNorm2);
        xCenterObs = synth.x0 - synth.x0(iPeakObs);
        [wObs, ~, ~] = extract_level_widths(xCenterObs, yObsNorm2, qLevels);
        rhoObs = wObs(:)' / max(sum(wObs), eps);

        widthDist = nan(nGap, 1);
        for j = 1:nGap
            widthDist(j) = norm(rhoObs - curve(j).rho);
        end
        [~, idxWidth] = min(widthDist);
        gHatWidth = gapList(idxWidth);

        row.gap_true_mm = curve(ig).gap;
        row.vib_case = vibCases(iv).name;
        row.vib_A_mm = vibCases(iv).A_mm;
        row.vib_f_Hz = vibCases(iv).f_Hz;
        row.distortion_eta = synth.eta;
        row.gHat_raw_mm = gHatRaw;
        row.gHat_width_mm = gHatWidth;
        row.err_raw_mm = gHatRaw - curve(ig).gap;
        row.err_width_mm = gHatWidth - curve(ig).gap;
        rows = [rows; row]; %#ok<AGROW>
    end
end

resultTable = struct2table(rows);

summaryRows = [];
for iv = 1:nVib
    mask = strcmp(resultTable.vib_case, vibCases(iv).name);
    srow.vib_case = vibCases(iv).name;
    srow.mean_eta = mean(resultTable.distortion_eta(mask));
    srow.rmse_raw_mm = sqrt(mean(resultTable.err_raw_mm(mask).^2));
    srow.rmse_width_mm = sqrt(mean(resultTable.err_width_mm(mask).^2));
    srow.mean_abs_raw_mm = mean(abs(resultTable.err_raw_mm(mask)));
    srow.mean_abs_width_mm = mean(abs(resultTable.err_width_mm(mask)));
    summaryRows = [summaryRows; srow]; %#ok<AGROW>
end
summaryTable = struct2table(summaryRows);

fprintf('\nRaw waveform vs width-signature gap estimation:\n');
disp(summaryTable);

%% 2) Representative case for direct visualization
[~, repGapIdx] = min(abs(gapList - 0.8));
repVibIdx = 4;  % Strong
repSynth = synthesize_vibrating_waveform(curve(repGapIdx).x, curve(repGapIdx).y, vibCases(repVibIdx), V_tip);
[~, ~, yRepNorm] = normalize_waveform(repSynth.yObs);
[~, iPeakRep] = max(yRepNorm);
xCenterRep = repSynth.x0 - repSynth.x0(iPeakRep);
[wRep, ~, ~] = extract_level_widths(xCenterRep, yRepNorm, qLevels);
rhoRep = wRep(:)' / max(sum(wRep), eps);

rawDistRep = nan(nGap, 1);
widthDistRep = nan(nGap, 1);
for j = 1:nGap
    yLibOnObs = interp1(curve(j).x, curve(j).y, repSynth.x0, 'pchip', NaN);
    [~, ~, yLibNorm] = normalize_waveform(yLibOnObs);
    valid = isfinite(yRepNorm) & isfinite(yLibNorm);
    rawDistRep(j) = sqrt(mean((yRepNorm(valid) - yLibNorm(valid)).^2));
    widthDistRep(j) = norm(rhoRep - curve(j).rho);
end
[~, idxRawRep] = min(rawDistRep);
[~, idxWidthRep] = min(widthDistRep);

%% 3) Visualize the difference
figure('Name', 'Raw waveform vs width signature under vibration', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 24, 14]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(curve(repGapIdx).x, curve(repGapIdx).y, 'k-', 'LineWidth', 1.4, 'DisplayName', 'True static');
plot(repSynth.x0, repSynth.yObs, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.1, ...
    'DisplayName', sprintf('Observed (%s)', vibCases(repVibIdx).name));
xlabel('x_0 (mm)');
ylabel('Capacitance');
title('(a) Same gap, vibration changes pointwise waveform');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
bar(categorical(string(gapList)), rawDistRep, 0.65, 'FaceColor', [0.85 0.35 0.10]);
xlabel('Library gap (mm)');
ylabel('Pointwise RMSE');
title(sprintf('(b) Raw waveform matching -> gHat = %.1f mm', gapList(idxRawRep)));

nexttile; hold on;
bar(categorical(string(gapList)), widthDistRep, 0.65, 'FaceColor', [0.10 0.45 0.75]);
xlabel('Library gap (mm)');
ylabel('||rho_{obs} - rho_{lib}||');
title(sprintf('(c) Width signature matching -> gHat = %.1f mm', gapList(idxWidthRep)));

nexttile; hold on;
plot(qLevels, rhoRep, 'ko-', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed rho');
plot(qLevels, curve(repGapIdx).rho, '-', 'Color', [0.10 0.45 0.75], 'LineWidth', 1.3, ...
    'DisplayName', sprintf('True gap %.1f', curve(repGapIdx).gap));
plot(qLevels, curve(idxRawRep).rho, '--', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.3, ...
    'DisplayName', sprintf('Raw-picked gap %.1f', gapList(idxRawRep)));
xlabel('Normalized level q');
ylabel('Normalized width signature');
title('(d) Width signature is nearly preserved');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(summaryTable.mean_eta, summaryTable.rmse_raw_mm, 's-', 'LineWidth', 1.4, ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], ...
    'DisplayName', 'Raw waveform pointwise');
plot(summaryTable.mean_eta, summaryTable.rmse_width_mm, 'o-', 'LineWidth', 1.4, ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], ...
    'DisplayName', 'Width signature');
xlabel('Distortion \eta');
ylabel('Gap RMSE (mm)');
title('(e) RMSE under increasing vibration');
legend('Location', 'northwest', 'Box', 'off');

nexttile; hold on;
plot(summaryTable.mean_eta, summaryTable.mean_abs_raw_mm, 's-', 'LineWidth', 1.4, ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], ...
    'DisplayName', 'Raw waveform pointwise');
plot(summaryTable.mean_eta, summaryTable.mean_abs_width_mm, 'o-', 'LineWidth', 1.4, ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], ...
    'DisplayName', 'Width signature');
xlabel('Distortion \eta');
ylabel('Mean |gap error| (mm)');
title('(f) Mean absolute error under increasing vibration');
legend('Location', 'northwest', 'Box', 'off');

fprintf('\nRepresentative case:\n');
fprintf('  true gap = %.1f mm, vibration = %s, eta = %.4f\n', ...
    curve(repGapIdx).gap, vibCases(repVibIdx).name, repSynth.eta);
fprintf('  raw waveform pointwise -> gHat = %.1f mm\n', gapList(idxRawRep));
fprintf('  width signature        -> gHat = %.1f mm\n', gapList(idxWidthRep));

%% Local functions
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

function [widths, xLeft, xRight] = extract_level_widths(xCenter, yNorm, qLevels)
    xCenter = xCenter(:);
    yNorm = yNorm(:);
    qLevels = qLevels(:);
    [~, iPeak] = max(yNorm);
    widths = nan(numel(qLevels), 1);
    xLeft = nan(numel(qLevels), 1);
    xRight = nan(numel(qLevels), 1);
    for k = 1:numel(qLevels)
        q = qLevels(k);
        if q >= min(yNorm(1:iPeak)) && q <= max(yNorm(1:iPeak)) && ...
           q >= min(yNorm(iPeak:end)) && q <= max(yNorm(iPeak:end))
            xL = interp1(yNorm(1:iPeak), xCenter(1:iPeak), q, 'linear');
            xR = interp1(flipud(yNorm(iPeak:end)), flipud(xCenter(iPeak:end)), q, 'linear');
            xLeft(k) = xL;
            xRight(k) = xR;
            widths(k) = xR - xL;
        end
    end
end

function synth = synthesize_vibrating_waveform(xRef, yRef, vibCase, V_tip)
    xRef = xRef(:);
    yRef = yRef(:);
    x0Center = mean(xRef);
    x0 = xRef;
    for iter = 1:4
        tTmp = (x0 - x0Center) / V_tip;
        [uTmp, ~] = vib_displacement(tTmp, vibCase);
        x0 = linspace(min(xRef) + max(uTmp), max(xRef) + min(uTmp), numel(xRef))';
    end
    t = (x0 - x0Center) / V_tip;
    [u, du_dt] = vib_displacement(t, vibCase);
    xi = x0 - u;
    if any(diff(xi) <= 0) || min(xi) < min(xRef) - 1e-8 || max(xi) > max(xRef) + 1e-8
        synth = [];
        return;
    end
    synth.x0 = x0;
    synth.yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    synth.eta = max(abs(du_dt)) / V_tip;
end

function [u, du_dt] = vib_displacement(t, vibCase)
    t = t(:);
    if vibCase.A_mm == 0 || vibCase.f_Hz == 0
        u = zeros(size(t));
        du_dt = zeros(size(t));
        return;
    end
    omega = 2 * pi * vibCase.f_Hz;
    u = vibCase.A_mm * sin(omega * t);
    du_dt = vibCase.A_mm * omega * cos(omega * t);
end
