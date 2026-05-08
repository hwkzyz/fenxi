%% Experimental clearance-vibration decoupling for BTT waveforms
% The RAW calibration files contain 16 clearance levels. Each row is a
% complete blade-passing waveform measured at one clearance. This script
% aligns the experimental waveforms, builds s(n,g), and evaluates whether a
% fixed-clearance template produces apparent timing shifts.

clc; clear; close all;

%% Paths
rootDir = fileparts(mfilename('fullpath'));
expDir = fullfile(rootDir, '标定原始数据-试验');
if ~isfolder(expDir)
    d = dir(rootDir);
    isHit = [d.isdir] & contains({d.name}, '试验');
    if any(isHit)
        expDir = fullfile(rootDir, d(find(isHit, 1, 'first')).name);
    else
        error('Experimental data folder was not found under: %s', rootDir);
    end
end

outDir = fullfile(rootDir, 'experimental_decoupling_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Figure defaults
fontName = 'Times New Roman';
set(groot, 'DefaultAxesFontName', fontName, ...
    'DefaultTextFontName', fontName, ...
    'DefaultLegendFontName', fontName, ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

colors = [0.85 0.35 0.10;
          0.85 0.65 0.10;
          0.20 0.60 0.60;
          0.10 0.25 0.50;
          0.36 0.36 0.36];

%% File discovery
rawFiles = dir(fullfile(expDir, '*RAW.csv'));
resultFiles = dir(fullfile(expDir, '*RESULT.csv'));
rawFiles = sort_by_name(rawFiles);
resultFiles = sort_by_name(resultFiles);

if isempty(rawFiles)
    error('No RAW files found in: %s', expDir);
end

allMetrics = table();
summaryRows = table();

for fileIdx = 1:numel(rawFiles)
    rawPath = fullfile(rawFiles(fileIdx).folder, rawFiles(fileIdx).name);
    sensorId = sprintf('sensor_%02d', fileIdx);
    if fileIdx <= numel(resultFiles)
        token = regexp(resultFiles(fileIdx).name, '_YWP_(240\d+)', 'tokens', 'once');
        if ~isempty(token)
            sensorId = token{1};
        end
    end

    [gapList, yCell] = read_experimental_raw(rawPath);
    nGap = numel(gapList);

    proc = repmat(struct('gap', [], 'x', [], 'y', [], 'baseline', [], 'peakIndex', []), nGap, 1);
    for k = 1:nGap
        [xRel, ySeg, baseline, iPeak] = preprocess_waveform(yCell{k});
        proc(k).gap = gapList(k);
        proc(k).x = xRel;
        proc(k).y = ySeg;
        proc(k).baseline = baseline;
        proc(k).peakIndex = iPeak;
    end

    xMin = max(arrayfun(@(w) min(w.x), proc));
    xMax = min(arrayfun(@(w) max(w.x), proc));
    xGrid = linspace(xMin, xMax, 701);
    S = zeros(nGap, numel(xGrid));
    for k = 1:nGap
        S(k, :) = interp1(proc(k).x, proc(k).y, xGrid, 'pchip');
    end

    % Remove any tiny negative baseline left after edge subtraction.
    S = S - min(S, [], 2);

    F = griddedInterpolant({gapList, xGrid}, S, 'linear', 'nearest');

    gRef = 1.10;
    [~, iRef] = min(abs(gapList - gRef));
    gRef = gapList(iRef);

    %% Normalized residuals and waveform features
    Sn = zeros(size(S));
    for k = 1:nGap
        Sn(k, :) = normalize_wave(S(k, :));
    end
    Sref = Sn(iRef, :);
    shapeRms = sqrt(mean((Sn - Sref).^2, 2));

    peakVal = max(S, [], 2);
    fwhmVal = zeros(nGap, 1);
    areaVal = zeros(nGap, 1);
    peakX = zeros(nGap, 1);
    for k = 1:nGap
        [~, imax] = max(S(k, :));
        peakX(k) = xGrid(imax);
        fwhmVal(k) = calc_fwhm(xGrid, S(k, :));
        areaVal(k) = trapz(xGrid, S(k, :));
    end

    %% Bias caused by a fixed reference-clearance template
    dxRange = [-35, 35];
    fitMask = xGrid >= min(xGrid) + max(abs(dxRange)) & ...
        xGrid <= max(xGrid) - max(abs(dxRange));
    dxBias = zeros(nGap, 1);
    biasRmse = zeros(nGap, 1);
    biasAtBound = false(nGap, 1);
    for k = 1:nGap
        yObs = S(k, :);
        cost = @(dx) scaled_template_rmse(yObs, response_eval(F, xGrid, gRef, dx), fitMask);
        [dxBias(k), biasRmse(k)] = fminbnd(cost, dxRange(1), dxRange(2));
        biasAtBound(k) = abs(dxBias(k) - dxRange(1)) < 0.5 || ...
            abs(dxBias(k) - dxRange(2)) < 0.5;
    end

    %% Sensitivity independence of timing shift and clearance
    dxStep = mean(diff(xGrid));
    dSdx = zeros(size(S));
    dSdg = zeros(size(S));
    rho = zeros(nGap, 1);
    condJ = zeros(nGap, 1);
    for k = 1:nGap
        dSdx(k, :) = gradient(S(k, :), dxStep);
        if k == 1
            dSdg(k, :) = (S(k+1, :) - S(k, :)) ./ (gapList(k+1) - gapList(k));
        elseif k == nGap
            dSdg(k, :) = (S(k, :) - S(k-1, :)) ./ (gapList(k) - gapList(k-1));
        else
            dSdg(k, :) = (S(k+1, :) - S(k-1, :)) ./ (gapList(k+1) - gapList(k-1));
        end

        active = S(k, :) > 0.08 * max(S(k, :));
        jx = dSdx(k, active)';
        jg = dSdg(k, active)';
        rho(k) = dot(jx, jg) / max(norm(jx) * norm(jg), eps);
        J = [jx, jg];
        condJ(k) = cond(J' * J);
    end

    %% Joint estimation example using the experimental response surface
    dxTrue = 8.0;                         % samples
    gTrue = gapList(round(0.55 * nGap));  % an internal clearance
    if abs(gTrue - gRef) < 1e-12 && nGap > 9
        gTrue = gapList(round(0.65 * nGap));
    end
    noiseRatio = 0.004;
    rng(10 + fileIdx);

    yClean = response_eval(F, xGrid, gTrue, dxTrue);
    yMeas = yClean + noiseRatio * range(yClean) * randn(size(yClean));

    costFixed = @(dx) scaled_template_rmse(yMeas, response_eval(F, xGrid, gRef, dx), fitMask);
    [dxFixed, rmseFixed] = fminbnd(costFixed, dxRange(1), dxRange(2));

    p0 = coarse_joint_initial_guess(yMeas, F, xGrid, gapList, dxRange, fitMask);
    costJoint = @(p) joint_cost(p, yMeas, F, xGrid, gapList, dxRange, fitMask);
    pHat = fminsearch(costJoint, p0, optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10));
    dxJoint = pHat(1);
    gJoint = pHat(2);
    rmseJoint = joint_cost(pHat, yMeas, F, xGrid, gapList, dxRange, fitMask);

    yFitFixed = scaled_template_fit(yMeas, response_eval(F, xGrid, gRef, dxFixed));
    yFitJoint = scaled_template_fit(yMeas, response_eval(F, xGrid, gJoint, dxJoint));

    %% Metrics
    sensorCol = repmat(string(sensorId), nGap, 1);
    metrics = table(sensorCol, gapList(:), peakVal, fwhmVal, areaVal, peakX, ...
        shapeRms, dxBias, biasRmse, rho, condJ, ...
        'VariableNames', {'sensor_id', 'gap_mm', 'peak_voltage', 'fwhm_samples', ...
        'area', 'peak_x_samples', 'normalized_shape_rms', 'apparent_dx_bias_samples', ...
        'fixed_template_rmse', 'sensitivity_correlation', 'normal_matrix_condition'});
    metrics.bias_at_search_bound = biasAtBound;
    allMetrics = [allMetrics; metrics]; %#ok<AGROW>

    oneSummary = table(string(sensorId), gRef, dxTrue, gTrue, dxFixed, dxJoint, ...
        gJoint, rmseFixed, rmseJoint, ...
        'VariableNames', {'sensor_id', 'reference_gap_mm', 'true_dx_samples', ...
        'true_gap_mm', 'fixed_dx_samples', 'joint_dx_samples', ...
        'joint_gap_mm', 'fixed_rmse', 'joint_rmse'});
    summaryRows = [summaryRows; oneSummary]; %#ok<AGROW>

    %% Figures
    make_sensor_figures(sensorId, outDir, xGrid, gapList, S, Sn, Sref, ...
        dxBias, rho, fwhmVal, gRef, yMeas, yFitFixed, yFitJoint, ...
        dxFixed, dxJoint, gJoint, colors);

    fprintf('\n=== Experimental sensor %s ===\n', sensorId);
    fprintf('RAW file: %s\n', rawFiles(fileIdx).name);
    fprintf('Reference clearance: %.3f mm\n', gRef);
    fprintf('Max normalized shape RMS: %.4g\n', max(shapeRms));
    fprintf('Max fixed-template apparent shift: %.4g samples\n', max(abs(dxBias)));
    fprintf('Joint test true:  dx = %.3f samples, g = %.3f mm\n', dxTrue, gTrue);
    fprintf('Fixed template:   dx = %.3f samples, RMSE = %.4g\n', dxFixed, rmseFixed);
    fprintf('Joint fit:        dx = %.3f samples, g = %.3f mm, RMSE = %.4g\n', ...
        dxJoint, gJoint, rmseJoint);
end

writetable(allMetrics, fullfile(outDir, 'experimental_decoupling_metrics.csv'));
writetable(summaryRows, fullfile(outDir, 'experimental_joint_fit_summary.csv'));

fprintf('\nAll experimental results saved to: %s\n', outDir);

%% Local functions
function files = sort_by_name(files)
    [~, idx] = sort({files.name});
    files = files(idx);
end

function [gapList, yCell] = read_experimental_raw(filePath)
    lines = readlines(filePath, 'EmptyLineRule', 'skip');
    labels = str2double(split(strtrim(lines(1)), ','));
    labels = labels(isfinite(labels));
    gapList = labels(:) / 1000; % um to mm
    nGap = numel(gapList);
    yCell = cell(nGap, 1);

    for k = 1:nGap
        vals = str2double(split(strtrim(lines(k+1)), ','));
        vals = vals(isfinite(vals));
        yCell{k} = vals(:);
    end
end

function [xRel, ySeg, baseline, iPeak] = preprocess_waveform(y)
    y = y(:);
    n = numel(y);
    edgeN = max(20, round(0.05 * n));
    edgeVals = [y(1:edgeN); y(end-edgeN+1:end)];
    baseline = median(edgeVals);
    y0 = y - baseline;

    if abs(min(y0)) > max(y0)
        y0 = -y0;
    end

    [peakVal, iPeak] = max(y0);
    threshold = max(0.025 * peakVal, median(abs(edgeVals - baseline)) * 6);
    active = find(y0 > threshold);
    if isempty(active)
        startIdx = 1;
        endIdx = n;
    else
        margin = max(50, round(0.12 * numel(active)));
        startIdx = max(1, active(1) - margin);
        endIdx = min(n, active(end) + margin);
    end

    idx = (startIdx:endIdx)';
    xRel = idx - iPeak;
    ySeg = y0(idx);
    ySeg = ySeg - min(ySeg);
end

function yn = normalize_wave(y)
    yn = (y - min(y)) ./ max(range(y), eps);
end

function width = calc_fwhm(x, y)
    y0 = min(y);
    halfLevel = y0 + 0.5 * range(y);
    above = y >= halfLevel;
    idx = find(diff(above) ~= 0);
    if numel(idx) < 2
        width = NaN;
        return;
    end
    xCross = zeros(numel(idx), 1);
    for ii = 1:numel(idx)
        i = idx(ii);
        xCross(ii) = interp1(y(i:i+1), x(i:i+1), halfLevel, 'linear', 'extrap');
    end
    width = xCross(end) - xCross(1);
end

function y = response_eval(F, xGrid, g, dx)
    y = F(g * ones(size(xGrid)), xGrid - dx);
end

function rmse = scaled_template_rmse(yObs, yTpl, fitMask)
    yFit = scaled_template_fit(yObs(fitMask), yTpl(fitMask));
    rmse = sqrt(mean((yObs(fitMask) - yFit).^2));
end

function yFit = scaled_template_fit(yObs, yTpl)
    A = [ones(numel(yTpl), 1), yTpl(:)];
    theta = A \ yObs(:);
    yFit = reshape(A * theta, size(yObs));
end

function p0 = coarse_joint_initial_guess(yObs, F, xGrid, gapList, dxRange, fitMask)
    dxGrid = linspace(dxRange(1), dxRange(2), 29);
    gGrid = linspace(min(gapList), max(gapList), 31);
    bestVal = Inf;
    p0 = [0, median(gapList)];

    for ig = 1:numel(gGrid)
        for id = 1:numel(dxGrid)
            yTpl = response_eval(F, xGrid, gGrid(ig), dxGrid(id));
            val = scaled_template_rmse(yObs, yTpl, fitMask);
            if val < bestVal
                bestVal = val;
                p0 = [dxGrid(id), gGrid(ig)];
            end
        end
    end
end

function val = joint_cost(p, yObs, F, xGrid, gapList, dxRange, fitMask)
    dx = p(1);
    g = p(2);
    if g < min(gapList) || g > max(gapList) || dx < dxRange(1) || dx > dxRange(2)
        val = 1e3 + 1e3 * (max(0, min(gapList)-g)^2 + ...
            max(0, g-max(gapList))^2 + max(0, dxRange(1)-dx)^2 + ...
            max(0, dx-dxRange(2))^2);
        return;
    end
    yTpl = response_eval(F, xGrid, g, dx);
    val = scaled_template_rmse(yObs, yTpl, fitMask);
end

function make_sensor_figures(sensorId, outDir, xGrid, gapList, S, Sn, Sref, ...
    dxBias, rho, fwhmVal, gRef, yMeas, yFitFixed, yFitJoint, ...
    dxFixed, dxJoint, gJoint, colors)

    nGap = numel(gapList);
    pick = unique(round(linspace(1, nGap, min(5, nGap))));

    fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
    tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    for ii = 1:numel(pick)
        k = pick(ii);
        plot(xGrid, S(k, :), 'Color', colors(mod(ii-1, size(colors,1))+1, :), ...
            'DisplayName', sprintf('g = %.2f mm', gapList(k)));
    end
    xlabel('Relative sample n');
    ylabel('Voltage (V)');
    legend('Location', 'northeast', 'Box', 'off');
    title('(a) Raw aligned waveforms');

    nexttile; hold on;
    for ii = 1:numel(pick)
        k = pick(ii);
        plot(xGrid, Sn(k, :), 'Color', colors(mod(ii-1, size(colors,1))+1, :));
    end
    xlabel('Relative sample n');
    ylabel('Normalized response');
    title('(b) Normalized shapes');

    nexttile; hold on;
    for ii = 1:numel(pick)
        k = pick(ii);
        plot(xGrid, Sn(k, :) - Sref, 'Color', colors(mod(ii-1, size(colors,1))+1, :));
    end
    yline(0, 'k-', 'LineWidth', 0.6);
    xlabel('Relative sample n');
    ylabel('Residual');
    title(sprintf('(c) Residual to %.2f mm template', gRef));

    nexttile; hold on;
    bar(gapList, dxBias, 0.85, 'FaceColor', [0.5529, 0.6941, 0.8863], ...
        'EdgeColor', 'k', 'LineWidth', 0.4);
    yline(0, 'k-', 'LineWidth', 0.6);
    xlabel('Clearance g (mm)');
    ylabel('Apparent \Deltan (samples)');
    title('(d) Fixed-template timing bias');

    exportgraphics(fig1, fullfile(outDir, sprintf('%s_fig1_shape_bias.png', sensorId)), 'Resolution', 300);
    exportgraphics(fig1, fullfile(outDir, sprintf('%s_fig1_shape_bias.pdf', sensorId)), 'ContentType', 'vector');

    fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 17, 10]);
    tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile; hold on;
    imagesc(xGrid, gapList, S);
    set(gca, 'YDir', 'normal');
    colormap(gca, parula);
    cb = colorbar;
    cb.Label.String = 'Voltage (V)';
    xlabel('Relative sample n');
    ylabel('Clearance g (mm)');
    title('(a) Experimental s(n,g)');

    nexttile; hold on;
    plot(gapList, fwhmVal, '-o', 'Color', [0.20, 0.60, 0.60], ...
        'MarkerFaceColor', [0.20, 0.60, 0.60]);
    xlabel('Clearance g (mm)');
    ylabel('FWHM (samples)');
    title('(b) Width variation');

    nexttile; hold on;
    plot(gapList, rho, '-s', 'Color', [0.85, 0.35, 0.10], ...
        'MarkerFaceColor', [0.85, 0.35, 0.10]);
    yline(0, 'k-', 'LineWidth', 0.6);
    ylim([-1, 1]);
    xlabel('Clearance g (mm)');
    ylabel('Correlation');
    title('(c) Sensitivity independence');

    nexttile; hold on;
    plot(xGrid, yMeas, 'ko', 'MarkerSize', 2.3, 'DisplayName', 'Measured');
    plot(xGrid, yFitFixed, '-', 'Color', [0.85, 0.35, 0.10], ...
        'DisplayName', sprintf('Fixed g: \\Deltan = %.2f', dxFixed));
    plot(xGrid, yFitJoint, '-', 'Color', [0.10, 0.25, 0.50], ...
        'DisplayName', sprintf('Joint: \\Deltan = %.2f, g = %.2f', dxJoint, gJoint));
    xlabel('Relative sample n');
    ylabel('Voltage (V)');
    legend('Location', 'northeast', 'Box', 'off');
    title('(d) Joint estimation example');

    exportgraphics(fig2, fullfile(outDir, sprintf('%s_fig2_surface_joint.png', sensorId)), 'Resolution', 300);
    exportgraphics(fig2, fullfile(outDir, sprintf('%s_fig2_surface_joint.pdf', sensorId)), 'ContentType', 'vector');
end
