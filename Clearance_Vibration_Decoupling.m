%% Clearance-vibration decoupling demo for field-effect BTT waveforms
% This script reads simulated waveforms at different tip clearances, builds
% an s(x,g) response surface, and demonstrates why amplitude normalization
% alone can create an apparent timing/displacement bias.

clc; clear; close all;

%% Paths
rootDir = fileparts(mfilename('fullpath'));
simDir = fullfile(rootDir, '间隙的影响');
if ~isfolder(simDir)
    simDir = rootDir;
end

outDir = fullfile(rootDir, 'decoupling_results');
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
          0.10 0.25 0.50];

%% Load densified waveforms
files = dir(fullfile(simDir, 'wave_*_加密.txt'));
if isempty(files)
    error('No densified waveform files were found in: %s', simDir);
end

wave = struct('gap', {}, 'x', {}, 'y', {}, 'name', {});
for i = 1:numel(files)
    token = regexp(files(i).name, '间隙([0-9.]+)mm', 'tokens', 'once');
    if isempty(token)
        continue;
    end

    fn = fullfile(files(i).folder, files(i).name);
    T = readtable(fn, 'FileType', 'text', 'VariableNamingRule', 'preserve');
    M = table2array(T);
    M = M(:, 1:2);
    M = M(all(isfinite(M), 2), :);
    [x, ia] = unique(M(:, 1), 'stable');
    y = M(ia, 2);

    wave(end+1).gap = str2double(token{1}); %#ok<SAGROW>
    wave(end).x = x(:);
    wave(end).y = y(:);
    wave(end).name = files(i).name;
end

[~, order] = sort([wave.gap]);
wave = wave(order);
gapList = [wave.gap]';
nGap = numel(gapList);

if nGap < 3
    error('At least three clearance waveforms are recommended for s(x,g) interpolation.');
end

%% Common grid and response surface S(g,x)
xMin = max(arrayfun(@(w) min(w.x), wave));
xMax = min(arrayfun(@(w) max(w.x), wave));
xGrid = linspace(xMin, xMax, 501);
S = zeros(nGap, numel(xGrid));

for k = 1:nGap
    S(k, :) = interp1(wave(k).x, wave(k).y, xGrid, 'pchip');
end

F = griddedInterpolant({gapList, xGrid}, S, 'linear', 'nearest');
gRef = 1.0;
if ~any(abs(gapList - gRef) < 1e-12)
    [~, iRef] = min(abs(gapList - gRef));
    gRef = gapList(iRef);
end

%% Normalized shape residuals
Sn = zeros(size(S));
for k = 1:nGap
    Sn(k, :) = normalize_wave(S(k, :));
end
Sref = interp1(gapList, Sn, gRef, 'linear');
shapeRms = sqrt(mean((Sn - Sref).^2, 2));

%% Feature trends
peakVal = max(S, [], 2);
baseVal = min(S, [], 2);
fwhmVal = zeros(nGap, 1);
areaVal = zeros(nGap, 1);
peakX = zeros(nGap, 1);

for k = 1:nGap
    [peakVal(k), imax] = max(S(k, :));
    peakX(k) = xGrid(imax);
    fwhmVal(k) = calc_fwhm(xGrid, S(k, :));
    areaVal(k) = trapz(xGrid, S(k, :) - min(S(k, :)));
end

%% Apparent displacement bias caused by using a fixed clearance template
dxBias = zeros(nGap, 1);
biasRmse = zeros(nGap, 1);
fitRange = [-0.8, 0.8];

for k = 1:nGap
    yObs = S(k, :);
    cost = @(dx) scaled_template_rmse(yObs, response_eval(F, xGrid, gRef, dx), xGrid, dx);
    [dxBias(k), biasRmse(k)] = fminbnd(cost, fitRange(1), fitRange(2));
end

%% Local sensitivity correlation between timing shift and clearance change
dxStep = mean(diff(xGrid));
dSdx = zeros(size(S));
dSdg = zeros(size(S));
rho = zeros(nGap, 1);

for k = 1:nGap
    dSdx(k, :) = gradient(S(k, :), dxStep);
    if k == 1
        dSdg(k, :) = (S(k+1, :) - S(k, :)) ./ (gapList(k+1) - gapList(k));
    elseif k == nGap
        dSdg(k, :) = (S(k, :) - S(k-1, :)) ./ (gapList(k) - gapList(k-1));
    else
        dSdg(k, :) = (S(k+1, :) - S(k-1, :)) ./ (gapList(k+1) - gapList(k-1));
    end

    active = S(k, :) > min(S(k, :)) + 0.10 * range(S(k, :));
    jx = dSdx(k, active);
    jg = dSdg(k, active);
    rho(k) = dot(jx, jg) / max(norm(jx) * norm(jg), eps);
end

%% Joint estimation example
dxTrue = 0.25;       % mm, equivalent circumferential vibration displacement
gTrue = 1.45;        % mm, unseen clearance
noiseRatio = 0.002;  % relative to waveform span
rng(2);

yClean = response_eval(F, xGrid, gTrue, dxTrue);
yMeas = yClean + noiseRatio * range(yClean) * randn(size(yClean));

costFixed = @(dx) scaled_template_rmse(yMeas, response_eval(F, xGrid, gRef, dx), xGrid, dx);
[dxFixed, rmseFixed] = fminbnd(costFixed, fitRange(1), fitRange(2));

p0 = [0, gRef];
costJoint = @(p) joint_cost(p, yMeas, F, xGrid, gapList);
pHat = fminsearch(costJoint, p0, optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10));
dxJoint = pHat(1);
gJoint = pHat(2);
rmseJoint = joint_cost(pHat, yMeas, F, xGrid, gapList);

yFitFixed = scaled_template_fit(yMeas, response_eval(F, xGrid, gRef, dxFixed), xGrid, dxFixed);
yFitJoint = scaled_template_fit(yMeas, response_eval(F, xGrid, gJoint, dxJoint), xGrid, dxJoint);

%% Export metrics
metrics = table(gapList, peakVal, baseVal, fwhmVal, areaVal, peakX, ...
    shapeRms, dxBias, biasRmse, rho, ...
    'VariableNames', {'gap_mm', 'peak', 'baseline', 'fwhm_mm', 'area', ...
    'peak_x_mm', 'normalized_shape_rms', 'apparent_dx_bias_mm', ...
    'fixed_template_rmse', 'sensitivity_correlation'});
writetable(metrics, fullfile(outDir, 'decoupling_metrics.csv'));

%% Figure 1: waveform shape and clearance-induced bias
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 10]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
for k = 1:nGap
    plot(xGrid, S(k, :), 'Color', colors(mod(k-1, size(colors,1))+1, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'northeast', 'Box', 'off');
title('(a) Raw waveforms');

nexttile; hold on;
for k = 1:nGap
    plot(xGrid, Sn(k, :), 'Color', colors(mod(k-1, size(colors,1))+1, :), ...
        'DisplayName', sprintf('g = %.1f mm', gapList(k)));
end
xlabel('Position x (mm)');
ylabel('Normalized response');
title('(b) Amplitude-normalized shapes');

nexttile; hold on;
for k = 1:nGap
    plot(xGrid, Sn(k, :) - Sref, 'Color', colors(mod(k-1, size(colors,1))+1, :));
end
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Position x (mm)');
ylabel('Residual');
title(sprintf('(c) Residual to %.1f mm template', gRef));

nexttile; hold on;
bar(gapList, dxBias, 0.55, 'FaceColor', [0.5529, 0.6941, 0.8863], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
yline(0, 'k-', 'LineWidth', 0.6);
xlabel('Actual clearance g (mm)');
ylabel('Apparent \Deltax (mm)');
title('(d) Fixed-template timing bias');

exportgraphics(fig1, fullfile(outDir, 'fig1_shape_and_bias.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_shape_and_bias.pdf'), 'ContentType', 'vector');

%% Figure 2: response surface and joint estimation
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 17, 10]);
tiledlayout(fig2, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
imagesc(xGrid, gapList, S);
set(gca, 'YDir', 'normal');
colormap(gca, parula);
cb = colorbar;
cb.Label.String = 'Capacitance (pF)';
xlabel('Position x (mm)');
ylabel('Clearance g (mm)');
title('(a) Response surface s(x,g)');

nexttile; hold on;
plot(gapList, fwhmVal, '-o', 'Color', [0.20, 0.60, 0.60], ...
    'MarkerFaceColor', [0.20, 0.60, 0.60]);
xlabel('Clearance g (mm)');
ylabel('FWHM (mm)');
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
plot(xGrid, yMeas, 'ko', 'MarkerSize', 2.5, 'DisplayName', 'Measured');
plot(xGrid, yFitFixed, '-', 'Color', [0.85, 0.35, 0.10], ...
    'DisplayName', sprintf('Fixed g: \\Deltax = %.3f', dxFixed));
plot(xGrid, yFitJoint, '-', 'Color', [0.10, 0.25, 0.50], ...
    'DisplayName', sprintf('Joint: \\Deltax = %.3f, g = %.3f', dxJoint, gJoint));
xlabel('Position x (mm)');
ylabel('Capacitance (pF)');
legend('Location', 'northeast', 'Box', 'off');
title('(d) Joint estimation example');

exportgraphics(fig2, fullfile(outDir, 'fig2_surface_and_joint_fit.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_surface_and_joint_fit.pdf'), 'ContentType', 'vector');

%% Console summary
fprintf('\n=== Clearance-vibration decoupling summary ===\n');
fprintf('Loaded %d simulated waveforms from: %s\n', nGap, simDir);
fprintf('Reference fixed-template clearance: %.3f mm\n', gRef);
fprintf('\n%-8s %-10s %-10s %-12s %-12s %-10s\n', ...
    'g(mm)', 'peak', 'FWHM', 'shapeRMS', 'bias_dx(mm)', 'rho');
for k = 1:nGap
    fprintf('%-8.3f %-10.4g %-10.4g %-12.4g %-12.4g %-10.4g\n', ...
        gapList(k), peakVal(k), fwhmVal(k), shapeRms(k), dxBias(k), rho(k));
end

fprintf('\nJoint-estimation test:\n');
fprintf('True:        dx = %.4f mm, g = %.4f mm\n', dxTrue, gTrue);
fprintf('Fixed g:     dx = %.4f mm, g = %.4f mm, RMSE = %.4g\n', ...
    dxFixed, gRef, rmseFixed);
fprintf('Joint fit:   dx = %.4f mm, g = %.4f mm, RMSE = %.4g\n', ...
    dxJoint, gJoint, rmseJoint);
fprintf('Results saved to: %s\n', outDir);

%% Local functions
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

function rmse = scaled_template_rmse(yObs, yTpl, xGrid, dx)
    valid = (xGrid - dx >= min(xGrid)) & (xGrid - dx <= max(xGrid));
    yFit = scaled_template_fit(yObs(valid), yTpl(valid), xGrid(valid), 0);
    rmse = sqrt(mean((yObs(valid) - yFit).^2));
end

function yFit = scaled_template_fit(yObs, yTpl, ~, ~)
    A = [ones(numel(yTpl), 1), yTpl(:)];
    theta = A \ yObs(:);
    yFit = reshape(A * theta, size(yObs));
end

function val = joint_cost(p, yObs, F, xGrid, gapList)
    dx = p(1);
    g = p(2);
    if g < min(gapList) || g > max(gapList) || abs(dx) > 1.5
        val = 1e3 + 1e3 * (max(0, min(gapList)-g)^2 + max(0, g-max(gapList))^2 + max(0, abs(dx)-1.5)^2);
        return;
    end
    yTpl = response_eval(F, xGrid, g, dx);
    val = scaled_template_rmse(yObs, yTpl, xGrid, dx);
end
