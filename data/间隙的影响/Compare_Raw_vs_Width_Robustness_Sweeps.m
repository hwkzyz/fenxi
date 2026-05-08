%% Compare raw waveform matching vs width-signature matching under harder conditions
% Harder tests:
%   1) phase sweep
%   2) additive noise sweep
%   3) dual-frequency vibration

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

%% Static library
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

%% 1) Phase sweep
phaseList = linspace(0, 2*pi, 25);
phaseRows = [];
for ig = 1:nGap
    for ip = 1:numel(phaseList)
        vib.kind = 'single';
        vib.A_mm = 0.60;
        vib.f_Hz = 800;
        vib.phi_rad = phaseList(ip);
        synth = synthesize_waveform(curve(ig).x, curve(ig).y, vib, V_tip, 0);
        if isempty(synth), continue; end
        [gHatRaw, gHatWidth] = estimate_gap_two_methods(synth, curve, gapList, qLevels);
        row.gap_true_mm = curve(ig).gap;
        row.phase_rad = phaseList(ip);
        row.eta = synth.eta;
        row.gHat_raw_mm = gHatRaw;
        row.gHat_width_mm = gHatWidth;
        row.err_raw_mm = gHatRaw - curve(ig).gap;
        row.err_width_mm = gHatWidth - curve(ig).gap;
        phaseRows = [phaseRows; row]; %#ok<AGROW>
    end
end
phaseTable = struct2table(phaseRows);

%% 2) Noise sweep
noiseList = [0, 0.002, 0.005, 0.01, 0.02];
noiseRows = [];
rng(42);
for ig = 1:nGap
    for in = 1:numel(noiseList)
        vib.kind = 'single';
        vib.A_mm = 0.60;
        vib.f_Hz = 800;
        vib.phi_rad = 0.9;
        synth = synthesize_waveform(curve(ig).x, curve(ig).y, vib, V_tip, noiseList(in));
        if isempty(synth), continue; end
        [gHatRaw, gHatWidth] = estimate_gap_two_methods(synth, curve, gapList, qLevels);
        row.gap_true_mm = curve(ig).gap;
        row.noise_std_ratio = noiseList(in);
        row.eta = synth.eta;
        row.gHat_raw_mm = gHatRaw;
        row.gHat_width_mm = gHatWidth;
        row.err_raw_mm = gHatRaw - curve(ig).gap;
        row.err_width_mm = gHatWidth - curve(ig).gap;
        noiseRows = [noiseRows; row]; %#ok<AGROW>
    end
end
noiseTable = struct2table(noiseRows);

%% 3) Dual-frequency cases
multiCases = {
    struct('name','single-strong','kind','single','A_mm',0.60,'f_Hz',800,'phi_rad',0.7), ...
    struct('name','dual-mild','kind','dual','A1_mm',0.35,'f1_Hz',500,'phi1_rad',0.4,'A2_mm',0.20,'f2_Hz',1200,'phi2_rad',1.1), ...
    struct('name','dual-strong','kind','dual','A1_mm',0.50,'f1_Hz',500,'phi1_rad',0.4,'A2_mm',0.35,'f2_Hz',1200,'phi2_rad',1.1)
    };
multiRows = [];
for ig = 1:nGap
    for ic = 1:numel(multiCases)
        vib = multiCases{ic};
        synth = synthesize_waveform(curve(ig).x, curve(ig).y, vib, V_tip, 0);
        if isempty(synth), continue; end
        [gHatRaw, gHatWidth] = estimate_gap_two_methods(synth, curve, gapList, qLevels);
        row.gap_true_mm = curve(ig).gap;
        row.case_name = string(vib.name);
        row.eta = synth.eta;
        row.gHat_raw_mm = gHatRaw;
        row.gHat_width_mm = gHatWidth;
        row.err_raw_mm = gHatRaw - curve(ig).gap;
        row.err_width_mm = gHatWidth - curve(ig).gap;
        multiRows = [multiRows; row]; %#ok<AGROW>
    end
end
multiTable = struct2table(multiRows);

%% Summaries
phaseSummary = summarize_by_group(phaseTable, "phase_rad");
noiseSummary = summarize_by_group(noiseTable, "noise_std_ratio");
multiSummary = summarize_by_group(multiTable, "case_name");

fprintf('\nPhase sweep summary:\n');
disp(phaseSummary);
fprintf('\nNoise sweep summary:\n');
disp(noiseSummary);
fprintf('\nMulti-frequency summary:\n');
disp(multiSummary);

%% Representative phase-sweep case
[~, repGapIdx] = min(abs(gapList - 1.0));
repPhase = pi/2;
vibRep.kind = 'single';
vibRep.A_mm = 0.60;
vibRep.f_Hz = 800;
vibRep.phi_rad = repPhase;
repSynth = synthesize_waveform(curve(repGapIdx).x, curve(repGapIdx).y, vibRep, V_tip, 0);
[gHatRawRep, gHatWidthRep, rawDistRep, widthDistRep, rhoRep] = estimate_gap_two_methods(repSynth, curve, gapList, qLevels);

%% Visualization
figure('Name', 'Raw vs width robustness sweeps', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 25, 15]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(phaseList, phase_rmse(phaseTable.err_raw_mm, phaseTable.phase_rad, phaseList), 's-', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], 'DisplayName', 'Raw waveform');
plot(phaseList, phase_rmse(phaseTable.err_width_mm, phaseTable.phase_rad, phaseList), 'o-', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Width signature');
xlabel('Phase \phi (rad)');
ylabel('Gap RMSE (mm)');
title('(a) Phase sweep');
legend('Location', 'northwest', 'Box', 'off');

nexttile; hold on;
plot(noiseSummary.group_value, noiseSummary.rmse_raw_mm, 's-', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], 'DisplayName', 'Raw waveform');
plot(noiseSummary.group_value, noiseSummary.rmse_width_mm, 'o-', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Width signature');
xlabel('Noise std / signal range');
ylabel('Gap RMSE (mm)');
title('(b) Noise sweep');
legend('Location', 'northwest', 'Box', 'off');

nexttile;
bar(categorical(string(multiSummary.group_value)), [multiSummary.rmse_raw_mm, multiSummary.rmse_width_mm], 'grouped');
ylabel('Gap RMSE (mm)');
title('(c) Multi-frequency cases');
legend({'Raw waveform', 'Width signature'}, 'Location', 'northwest', 'Box', 'off');

nexttile; hold on;
plot(curve(repGapIdx).x, curve(repGapIdx).y, 'k-', 'LineWidth', 1.3, 'DisplayName', 'True static');
plot(repSynth.x0, repSynth.yObs, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.1, ...
    'DisplayName', 'Observed');
xlabel('x_0 (mm)');
ylabel('Capacitance');
title('(d) Representative vibrating waveform');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
bar(categorical(string(gapList)), rawDistRep, 0.65, 'FaceColor', [0.85 0.35 0.10]);
xlabel('Library gap (mm)');
ylabel('Pointwise RMSE');
title(sprintf('(e) Raw waveform -> gHat = %.1f mm', gHatRawRep));

nexttile; hold on;
plot(qLevels, rhoRep, 'ko-', 'LineWidth', 1.1, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed rho');
plot(qLevels, curve(repGapIdx).rho, '-', 'Color', [0.10 0.45 0.75], 'LineWidth', 1.3, ...
    'DisplayName', sprintf('True gap %.1f', curve(repGapIdx).gap));
plot(qLevels, curve(find(gapList == gHatWidthRep, 1)).rho, '--', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.3, ...
    'DisplayName', sprintf('Width-picked gap %.1f', gHatWidthRep));
xlabel('Normalized level q');
ylabel('Normalized width signature');
title(sprintf('(f) Width signature -> gHat = %.1f mm', gHatWidthRep));
legend('Location', 'best', 'Box', 'off');

fprintf('\nRepresentative phase-sweep case:\n');
fprintf('  true gap = %.1f mm, strong vibration, phi = %.2f rad, eta = %.4f\n', curve(repGapIdx).gap, repPhase, repSynth.eta);
fprintf('  raw waveform -> gHat = %.1f mm\n', gHatRawRep);
fprintf('  width signature -> gHat = %.1f mm\n', gHatWidthRep);

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
            xLeft(k) = interp1(yNorm(1:iPeak), xCenter(1:iPeak), q, 'linear');
            xRight(k) = interp1(flipud(yNorm(iPeak:end)), flipud(xCenter(iPeak:end)), q, 'linear');
            widths(k) = xRight(k) - xLeft(k);
        end
    end
end

function synth = synthesize_waveform(xRef, yRef, vib, V_tip, noiseStdRatio)
    xRef = xRef(:);
    yRef = yRef(:);
    x0Center = mean(xRef);
    x0 = xRef;
    for iter = 1:4
        tTmp = (x0 - x0Center) / V_tip;
        [uTmp, ~] = displacement_model(tTmp, vib);
        x0 = linspace(min(xRef) + max(uTmp), max(xRef) + min(uTmp), numel(xRef))';
    end
    t = (x0 - x0Center) / V_tip;
    [u, du_dt] = displacement_model(t, vib);
    xi = x0 - u;
    if any(diff(xi) <= 0) || min(xi) < min(xRef) - 1e-8 || max(xi) > max(xRef) + 1e-8
        synth = [];
        return;
    end
    yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    sigRange = max(yObs) - min(yObs);
    if noiseStdRatio > 0
        yObs = yObs + noiseStdRatio * sigRange * randn(size(yObs));
    end
    synth.x0 = x0;
    synth.yObs = yObs;
    synth.eta = max(abs(du_dt)) / V_tip;
end

function [u, du_dt] = displacement_model(t, vib)
    t = t(:);
    switch vib.kind
        case 'single'
            omega = 2 * pi * vib.f_Hz;
            u = vib.A_mm * sin(omega * t + vib.phi_rad);
            du_dt = vib.A_mm * omega * cos(omega * t + vib.phi_rad);
        case 'dual'
            omega1 = 2 * pi * vib.f1_Hz;
            omega2 = 2 * pi * vib.f2_Hz;
            u = vib.A1_mm * sin(omega1 * t + vib.phi1_rad) + ...
                vib.A2_mm * sin(omega2 * t + vib.phi2_rad);
            du_dt = vib.A1_mm * omega1 * cos(omega1 * t + vib.phi1_rad) + ...
                    vib.A2_mm * omega2 * cos(omega2 * t + vib.phi2_rad);
        otherwise
            error('Unknown vibration kind.');
    end
end

function [gHatRaw, gHatWidth, rawDist, widthDist, rhoObs] = estimate_gap_two_methods(synth, curve, gapList, qLevels)
    nGap = numel(gapList);
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
end

function summaryTable = summarize_by_group(T, groupVar)
    groups = unique(T.(groupVar));
    rows = [];
    for i = 1:numel(groups)
        g = groups(i);
        mask = T.(groupVar) == g;
        row.group_value = g;
        row.rmse_raw_mm = sqrt(mean(T.err_raw_mm(mask).^2));
        row.rmse_width_mm = sqrt(mean(T.err_width_mm(mask).^2));
        row.mean_abs_raw_mm = mean(abs(T.err_raw_mm(mask)));
        row.mean_abs_width_mm = mean(abs(T.err_width_mm(mask)));
        rows = [rows; row]; %#ok<AGROW>
    end
    summaryTable = struct2table(rows);
end

function y = phase_rmse(err, phi, phiGrid)
    y = zeros(size(phiGrid));
    for i = 1:numel(phiGrid)
        mask = abs(phi - phiGrid(i)) < 1e-12;
        y(i) = sqrt(mean(err(mask).^2));
    end
end
