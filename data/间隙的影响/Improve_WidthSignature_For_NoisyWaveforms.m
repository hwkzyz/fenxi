%% Improve width-signature matching for noisy waveforms
% Strategy:
%   1) denoise each observed waveform with robust smoothing
%   2) enforce monotone left/right branches by isotonic regression
%   3) extract width signature + auxiliary geometric features
%   4) estimate gap by local weighted feature fusion

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

%% Config
qLevels = (0.10:0.05:0.90)';
noiseList = [0, 0.002, 0.005, 0.01, 0.02];
V_tip = 3.0e5;  % mm/s
nTrials = 20;

vib.kind = 'single';
vib.A_mm = 0.60;
vib.f_Hz = 800;
vib.phi_rad = 0.9;

opts.medWindow = 9;
opts.meanWindow = 15;
opts.edgeFrac = 0.08;

%% Build static library
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);

lib = struct();
for ig = 1:nGap
    lib(ig).gap = gapList(ig);
    lib(ig).x = xCell{ig}(:);
    lib(ig).y = yCell{ig}(:);
    sig = extract_signature_basic(lib(ig).x, lib(ig).y, qLevels);
    lib(ig).rho_basic = sig.rho;
    lib(ig).fwhm_basic = sig.fwhm;
    lib(ig).area_basic = sig.area_norm;

    sigR = extract_signature_robust(lib(ig).x, lib(ig).y, qLevels, opts);
    lib(ig).rho_robust = sigR.rho;
    lib(ig).fwhm_robust = sigR.fwhm;
    lib(ig).area_robust = sigR.area_norm;
end

rhoMat = vertcat(lib.rho_robust);
fwhmVec = vertcat(lib.fwhm_robust);
areaVec = vertcat(lib.area_robust);
scale.rho = max(std(rhoMat, 0, 1), 1e-6);
scale.fwhm = max(std(fwhmVec), 1e-6);
scale.area = max(std(areaVec), 1e-6);

%% Noise sweep
rng(42);
rows = [];
rep = struct();
repCaptured = false;

for in = 1:numel(noiseList)
    noiseStdRatio = noiseList(in);
    for trial = 1:nTrials
        for ig = 1:nGap
            synth = synthesize_waveform(lib(ig).x, lib(ig).y, vib, V_tip, noiseStdRatio);
            if isempty(synth), continue; end

            sigBasic = extract_signature_basic(synth.x0, synth.yObs, qLevels);
            [gHatBasicNN, distBasic] = match_basic_width(sigBasic, lib, gapList);

            sigRobust = extract_signature_robust(synth.x0, synth.yObs, qLevels, opts);
            [gHatRobustNN, distRobust] = match_robust_width(sigRobust, lib, gapList);
            gHatRobustFuse = match_robust_fused(sigRobust, lib, gapList, scale);

            row.gap_true_mm = lib(ig).gap;
            row.noise_std_ratio = noiseStdRatio;
            row.trial = trial;
            row.eta = synth.eta;
            row.gHat_basic_nn_mm = gHatBasicNN;
            row.gHat_robust_nn_mm = gHatRobustNN;
            row.gHat_robust_fused_mm = gHatRobustFuse;
            row.err_basic_nn_mm = gHatBasicNN - lib(ig).gap;
            row.err_robust_nn_mm = gHatRobustNN - lib(ig).gap;
            row.err_robust_fused_mm = gHatRobustFuse - lib(ig).gap;
            rows = [rows; row]; %#ok<AGROW>

            if ~repCaptured && abs(noiseStdRatio - 0.01) < 1e-12 && abs(lib(ig).gap - 1.0) < 1e-12 && trial == 1
                repCaptured = true;
                rep.x0 = synth.x0;
                rep.yObs = synth.yObs;
                rep.sigBasic = sigBasic;
                rep.sigRobust = sigRobust;
                rep.distBasic = distBasic;
                rep.distRobust = distRobust;
                rep.gHatBasicNN = gHatBasicNN;
                rep.gHatRobustNN = gHatRobustNN;
                rep.gHatRobustFuse = gHatRobustFuse;
                rep.trueGap = lib(ig).gap;
                rep.trueStaticX = lib(ig).x;
                rep.trueStaticY = lib(ig).y;
            end
        end
    end
end

resultTable = struct2table(rows);

summaryRows = [];
for in = 1:numel(noiseList)
    mask = abs(resultTable.noise_std_ratio - noiseList(in)) < 1e-12;
    srow.noise_std_ratio = noiseList(in);
    srow.rmse_basic_nn_mm = sqrt(mean(resultTable.err_basic_nn_mm(mask).^2));
    srow.rmse_robust_nn_mm = sqrt(mean(resultTable.err_robust_nn_mm(mask).^2));
    srow.rmse_robust_fused_mm = sqrt(mean(resultTable.err_robust_fused_mm(mask).^2));
    srow.mean_abs_basic_nn_mm = mean(abs(resultTable.err_basic_nn_mm(mask)));
    srow.mean_abs_robust_nn_mm = mean(abs(resultTable.err_robust_nn_mm(mask)));
    srow.mean_abs_robust_fused_mm = mean(abs(resultTable.err_robust_fused_mm(mask)));
    summaryRows = [summaryRows; srow]; %#ok<AGROW>
end
summaryTable = struct2table(summaryRows);

fprintf('\nNoise-robust width-signature comparison:\n');
disp(summaryTable);

%% Visualize
figure('Name', 'Improved width signature for noisy waveforms', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 25, 15]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
plot(summaryTable.noise_std_ratio, summaryTable.rmse_basic_nn_mm, 's-', ...
    'Color', [0.70 0.35 0.10], 'MarkerFaceColor', [0.70 0.35 0.10], 'DisplayName', 'Original width NN');
plot(summaryTable.noise_std_ratio, summaryTable.rmse_robust_nn_mm, 'o-', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Robust width NN');
plot(summaryTable.noise_std_ratio, summaryTable.rmse_robust_fused_mm, 'd-', ...
    'Color', [0.15 0.60 0.30], 'MarkerFaceColor', [0.15 0.60 0.30], 'DisplayName', 'Robust fused interp');
xlabel('Noise std / signal range');
ylabel('Gap RMSE (mm)');
title('(a) RMSE under noise');
legend('Location', 'northwest', 'Box', 'off');

nexttile; hold on;
plot(summaryTable.noise_std_ratio, summaryTable.mean_abs_basic_nn_mm, 's-', ...
    'Color', [0.70 0.35 0.10], 'MarkerFaceColor', [0.70 0.35 0.10], 'DisplayName', 'Original width NN');
plot(summaryTable.noise_std_ratio, summaryTable.mean_abs_robust_nn_mm, 'o-', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Robust width NN');
plot(summaryTable.noise_std_ratio, summaryTable.mean_abs_robust_fused_mm, 'd-', ...
    'Color', [0.15 0.60 0.30], 'MarkerFaceColor', [0.15 0.60 0.30], 'DisplayName', 'Robust fused interp');
xlabel('Noise std / signal range');
ylabel('Mean |gap error| (mm)');
title('(b) Mean absolute error under noise');
legend('Location', 'northwest', 'Box', 'off');

nexttile; hold on;
plot(rep.trueStaticX, rep.trueStaticY, 'k-', 'LineWidth', 1.2, 'DisplayName', 'True static');
plot(rep.x0, rep.yObs, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 0.9, 'DisplayName', 'Noisy observed');
plot(rep.sigRobust.x, rep.sigRobust.yFitRawScale, '-', 'Color', [0.10 0.45 0.75], 'LineWidth', 1.2, ...
    'DisplayName', 'Robust fitted waveform');
xlabel('x_0 (mm)');
ylabel('Capacitance');
title('(c) Noisy waveform and robust fit');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
plot(qLevels, rep.sigBasic.rho, 's--', 'Color', [0.70 0.35 0.10], 'DisplayName', 'Original width rho');
plot(qLevels, rep.sigRobust.rho, 'o-', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Robust width rho');
plot(qLevels, lib(find(gapList == rep.trueGap, 1)).rho_robust, 'k-', 'LineWidth', 1.2, 'DisplayName', 'True library rho');
xlabel('Normalized level q');
ylabel('Normalized width signature');
title('(d) Signature after denoise + monotone fit');
legend('Location', 'best', 'Box', 'off');

nexttile; hold on;
bar(categorical(string(gapList)), rep.distBasic, 0.65, 'FaceColor', [0.70 0.35 0.10]);
xlabel('Library gap (mm)');
ylabel('||rho_{obs} - rho_{lib}||');
title(sprintf('(e) Original width NN -> %.1f mm', rep.gHatBasicNN));

nexttile; hold on;
fusedScores = nan(nGap, 1);
for j = 1:nGap
    fusedScores(j) = fused_feature_distance(rep.sigRobust, lib(j), scale);
end
bar(categorical(string(gapList)), fusedScores, 0.65, 'FaceColor', [0.15 0.60 0.30]);
xlabel('Library gap (mm)');
ylabel('Fused feature distance');
title(sprintf('(f) Robust fused interp -> %.3f mm', rep.gHatRobustFuse));

fprintf('\nRepresentative case:\n');
fprintf('  true gap = %.1f mm, noise std ratio = 0.01, eta = %.4f\n', rep.trueGap, vib.A_mm * 2*pi*vib.f_Hz / V_tip);
fprintf('  original width NN  -> gHat = %.1f mm\n', rep.gHatBasicNN);
fprintf('  robust width NN    -> gHat = %.1f mm\n', rep.gHatRobustNN);
fprintf('  robust fused interp-> gHat = %.3f mm\n', rep.gHatRobustFuse);

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

function synth = synthesize_waveform(xRef, yRef, vib, V_tip, noiseStdRatio)
    xRef = xRef(:);
    yRef = yRef(:);
    x0Center = mean(xRef);
    x0 = xRef;
    for iter = 1:4
        tTmp = (x0 - x0Center) / V_tip;
        omega = 2 * pi * vib.f_Hz;
        uTmp = vib.A_mm * sin(omega * tTmp + vib.phi_rad);
        x0 = linspace(min(xRef) + max(uTmp), max(xRef) + min(uTmp), numel(xRef))';
    end
    t = (x0 - x0Center) / V_tip;
    omega = 2 * pi * vib.f_Hz;
    u = vib.A_mm * sin(omega * t + vib.phi_rad);
    du_dt = vib.A_mm * omega * cos(omega * t + vib.phi_rad);
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

function sig = extract_signature_basic(x, y, qLevels)
    [base, amp, yNorm] = normalize_waveform_minmax(y);
    [~, iPeak] = max(yNorm);
    xCenter = x(:) - x(iPeak);
    [w, xL, xR] = extract_level_widths_from_monotone(xCenter, yNorm, qLevels);
    sig.x = x(:);
    sig.base = base;
    sig.amp = amp;
    sig.yNorm = yNorm;
    sig.yFitRawScale = y(:);
    sig.rho = w(:)' / max(sum(w), eps);
    sig.fwhm = calc_width_at_level(xCenter, yNorm, 0.5);
    sig.area_norm = trapz(xCenter, max(yNorm, 0));
    sig.xL = xL;
    sig.xR = xR;
end

function sig = extract_signature_robust(x, y, qLevels, opts)
    x = x(:);
    y = y(:);
    yMed = smoothdata(y, 'movmedian', opts.medWindow);
    ySm = smoothdata(yMed, 'movmean', opts.meanWindow);

    n = numel(ySm);
    nEdge = max(5, round(opts.edgeFrac * n));
    base = median([ySm(1:nEdge); ySm(end-nEdge+1:end)]);
    yBc = ySm - base;
    yBc = max(yBc, 0);
    amp = max(yBc);
    yNorm0 = yBc / max(amp, eps);

    [~, iPeak] = max(yNorm0);
    leftFit = pava_increasing(yNorm0(1:iPeak));
    rightFit = pava_decreasing(yNorm0(iPeak:end));
    yFit = [leftFit; rightFit(2:end)];
    yFit = yFit / max(max(yFit), eps);

    xCenter = x - x(iPeak);
    [w, xL, xR] = extract_level_widths_from_monotone(xCenter, yFit, qLevels);

    sig.x = x;
    sig.base = base;
    sig.amp = amp;
    sig.yNorm = yFit;
    sig.yFitRawScale = base + amp * yFit;
    sig.rho = w(:)' / max(sum(w), eps);
    sig.fwhm = calc_width_at_level(xCenter, yFit, 0.5);
    sig.area_norm = trapz(xCenter, max(yFit, 0));
    sig.xL = xL;
    sig.xR = xR;
end

function [base, amp, yNorm] = normalize_waveform_minmax(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function [w, xL, xR] = extract_level_widths_from_monotone(xCenter, yMono, qLevels)
    xCenter = xCenter(:);
    yMono = yMono(:);
    qLevels = qLevels(:);
    [~, iPeak] = max(yMono);
    yLeft = yMono(1:iPeak);
    xLeftCurve = xCenter(1:iPeak);
    yRight = yMono(iPeak:end);
    xRightCurve = xCenter(iPeak:end);

    [yLUniq, iaL] = unique(yLeft, 'stable');
    xLUniq = xLeftCurve(iaL);
    [yRUniq, iaR] = unique(flipud(yRight), 'stable');
    xRUniq = flipud(xRightCurve);
    xRUniq = xRUniq(iaR);

    xL = nan(numel(qLevels), 1);
    xR = nan(numel(qLevels), 1);
    for k = 1:numel(qLevels)
        q = qLevels(k);
        xL(k) = interp1(yLUniq, xLUniq, q, 'linear', 'extrap');
        xR(k) = interp1(yRUniq, xRUniq, q, 'linear', 'extrap');
    end
    w = xR - xL;
end

function width = calc_width_at_level(xCenter, yMono, q)
    [w, ~, ~] = extract_level_widths_from_monotone(xCenter, yMono, q);
    width = w(1);
end

function yFit = pava_increasing(y)
    y = y(:);
    n = numel(y);
    lvl = y;
    wt = ones(n, 1);
    idxStart = (1:n).';
    idxEnd = (1:n).';
    m = n;
    i = 1;
    while i < m
        if lvl(i) <= lvl(i + 1)
            i = i + 1;
        else
            newWt = wt(i) + wt(i + 1);
            newLvl = (wt(i) * lvl(i) + wt(i + 1) * lvl(i + 1)) / newWt;
            lvl(i) = newLvl;
            wt(i) = newWt;
            idxEnd(i) = idxEnd(i + 1);
            lvl(i + 1:m - 1) = lvl(i + 2:m);
            wt(i + 1:m - 1) = wt(i + 2:m);
            idxStart(i + 1:m - 1) = idxStart(i + 2:m);
            idxEnd(i + 1:m - 1) = idxEnd(i + 2:m);
            m = m - 1;
            if i > 1
                i = i - 1;
            end
        end
    end
    yFit = zeros(n, 1);
    for j = 1:m
        yFit(idxStart(j):idxEnd(j)) = lvl(j);
    end
end

function yFit = pava_decreasing(y)
    yFit = -pava_increasing(-y(:));
end

function [gHat, distVec] = match_basic_width(sig, lib, gapList)
    nGap = numel(gapList);
    distVec = nan(nGap, 1);
    for j = 1:nGap
        distVec(j) = norm(sig.rho - lib(j).rho_basic);
    end
    [~, idx] = min(distVec);
    gHat = gapList(idx);
end

function [gHat, distVec] = match_robust_width(sig, lib, gapList)
    nGap = numel(gapList);
    distVec = nan(nGap, 1);
    for j = 1:nGap
        distVec(j) = norm(sig.rho - lib(j).rho_robust);
    end
    [~, idx] = min(distVec);
    gHat = gapList(idx);
end

function gHat = match_robust_fused(sig, lib, gapList, scale)
    nGap = numel(gapList);
    scores = nan(nGap, 1);
    for j = 1:nGap
        scores(j) = fused_feature_distance(sig, lib(j), scale);
    end
    [scoreSort, idxSort] = sort(scores, 'ascend');
    topK = min(3, nGap);
    weights = 1 ./ max(scoreSort(1:topK), 1e-9);
    weights = weights / sum(weights);
    gHat = sum(gapList(idxSort(1:topK)) .* weights);
end

function score = fused_feature_distance(sig, libItem, scale)
    drho = (sig.rho - libItem.rho_robust) ./ scale.rho;
    dfwhm = (sig.fwhm - libItem.fwhm_robust) / scale.fwhm;
    darea = (sig.area_norm - libItem.area_robust) / scale.area;
    score = sqrt(mean(drho.^2)) + 0.35 * abs(dfwhm) + 0.25 * abs(darea);
end
