%% Step05J: Learn blade geometry modes from all static gap waveforms
% Experimental theory/diagnostic branch. It does not overwrite previous
% Step05/Step05H/Step05I outputs.
%
% Starting from the offset-corrected shared response surface F0(g,x), this
% step decomposes each blade's residual into first-order physical modes:
%
%   R_bj(x) = y_bj(x) - F0(g_j + delta_g_b, x)
%          ~= ddelta_b * F_g
%            + mu_b * x * F_g
%            - tau_b * F_x
%            + kappa_b * x * F_x
%            + c_b
%
% where F_g = dF/dg and F_x = dF/dx. The learned coefficients indicate
% whether remaining waveform differences are mainly clearance, tilt,
% x-registration, or width/scale effects.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step05j_blade_geometry_modes');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

surfaceFile = fullfile(outDir, 'Step05H_OffsetCorrected_Shared_Response_Surface_20251222.mat');
if ~isfile(surfaceFile)
    error('Run Step05H first. Missing file: %s', surfaceFile);
end

S = load(surfaceFile, 'OffsetCorrectedResponseSurface');
Surf = S.OffsetCorrectedResponseSurface;

cfg = struct();
cfg.robustSigma = 4.0;
cfg.ridgeLambda = 1e-3;
cfg.minPointsPerBlade = 300;
cfg.useRobustSecondPass = true;

xGrid = Surf.xGrid(:);
bladeIds = Surf.bladeIds(:).';
trueGaps = Surf.trueGapMm(:);
recordedGaps = Surf.recordedGapMm(:);
waveforms = Surf.waveforms;
effectiveWindow = Surf.effectiveWindow(:);
coeff = Surf.coeff;
g0Model = Surf.g0Mm;
offsetTable = Surf.offsetByBlade;

offsetByBlade = nan(size(bladeIds));
for ib = 1:numel(bladeIds)
    row = offsetTable.bladeId == bladeIds(ib);
    if ~any(row)
        error('No deltaGapOffsetMm for blade %d.', bladeIds(ib));
    end
    offsetByBlade(ib) = offsetTable.deltaGapOffsetMm(find(row, 1));
end

fprintf('\n=== Step05J: first-order blade geometry modes ===\n');
fprintf('Source surface: %s\n', surfaceFile);
fprintf('Model: R ~= ddelta*Fg + mu*xFg - tau*Fx + kappa*xFx + c\n');

geometryRows = cell(numel(bladeIds), 1);
fitRows = {};
bladeCache = repmat(struct('bladeId', NaN, 'x', [], 'resBefore', [], 'resAfter', [], ...
    'predBefore', [], 'predAfter', [], 'obs', []), numel(bladeIds), 1);

for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    deltaRef = offsetByBlade(ib);
    [A, yRes, obs, pred0, xAll, gapAll, gapIndexAll] = build_blade_regression_local( ...
        squeeze(waveforms(:, ib, :)), xGrid, trueGaps, coeff, g0Model, deltaRef, effectiveWindow);
    valid = all(isfinite(A), 2) & isfinite(yRes) & isfinite(obs) & isfinite(pred0);
    A = A(valid, :);
    yRes = yRes(valid);
    obs = obs(valid);
    pred0 = pred0(valid);
    xAll = xAll(valid);
    gapAll = gapAll(valid);
    gapIndexAll = gapIndexAll(valid);

    if nnz(valid) < cfg.minPointsPerBlade
        error('Blade %d has too few valid points for geometry-mode learning.', bladeId);
    end

    beta = scaled_ridge_solve_local(A, yRes, cfg.ridgeLambda);
    predCorr = pred0 + A * beta;
    resBefore = obs - pred0;
    resAfter = obs - predCorr;

    if cfg.useRobustSecondPass
        sigma = 1.4826 * mad(resAfter, 1);
        if ~isfinite(sigma) || sigma <= 0
            sigma = std(resAfter, 'omitnan');
        end
        keep = abs(resAfter) <= cfg.robustSigma * max(sigma, eps);
        beta = scaled_ridge_solve_local(A(keep, :), yRes(keep), cfg.ridgeLambda);
        predCorr = pred0 + A * beta;
        resBefore = obs - pred0;
        resAfter = obs - predCorr;
    end

    rmseBefore = sqrt(mean(resBefore.^2, 'omitnan'));
    rmseAfter = sqrt(mean(resAfter.^2, 'omitnan'));
    contribution = mode_contribution_local(A, beta);

    geometryRows{ib} = table(bladeId, deltaRef, beta(1), deltaRef + beta(1), ...
        beta(2), atan(beta(2)) * 180 / pi, beta(3), beta(4), 1 + beta(4), ...
        beta(5), rmseBefore, rmseAfter, 1 - rmseAfter / max(rmseBefore, eps), ...
        contribution(1), contribution(2), contribution(3), contribution(4), contribution(5), ...
        'VariableNames', {'bladeId','deltaPriorMm','deltaCorrectionMm','deltaTotalMm', ...
        'muGapPerXMm','tiltAngleDeg','tauMm','kappaScale','xScale', ...
        'offsetMv','rmseBeforeMv','rmseAfterMv','rmseReductionRatio', ...
        'contribDelta','contribTilt','contribTau','contribScale','contribOffset'});

    for ig = 1:numel(trueGaps)
        rows = gapIndexAll == ig;
        if ~any(rows)
            continue;
        end
        fitRows{end+1, 1} = table(bladeId, recordedGaps(ig), trueGaps(ig), ...
            deltaRef, beta(1), beta(2), beta(3), beta(4), ...
            sqrt(mean(resBefore(rows).^2, 'omitnan')), ...
            sqrt(mean(resAfter(rows).^2, 'omitnan')), nnz(rows), ...
            'VariableNames', {'bladeId','recordedGapMm','trueReferenceGapMm', ...
            'deltaPriorMm','deltaCorrectionMm','muGapPerXMm','tauMm','kappaScale', ...
            'rmseBeforeMv','rmseAfterMv','pointCount'}); %#ok<SAGROW>
    end

    bladeCache(ib).bladeId = bladeId;
    bladeCache(ib).x = xAll;
    bladeCache(ib).resBefore = resBefore;
    bladeCache(ib).resAfter = resAfter;
    bladeCache(ib).predBefore = pred0;
    bladeCache(ib).predAfter = predCorr;
    bladeCache(ib).obs = obs;

    fprintf('  B%d: RMSE %.2f -> %.2f mV | ddelta %.4f, mu %.4f, tau %.4f, kappa %.4f\n', ...
        bladeId, rmseBefore, rmseAfter, beta(1), beta(2), beta(3), beta(4));
end

geometryTable = vertcat(geometryRows{:});
fitTable = vertcat(fitRows{:});

BladeGeometryModes = struct();
BladeGeometryModes.dataset = '20251222';
BladeGeometryModes.method = 'first_order_blade_geometry_mode_learning';
BladeGeometryModes.description = ['Residuals around the offset-corrected shared response surface are decomposed ' ...
    'into clearance, tilt, x-registration, and width/scale modes using F_g and F_x.'];
BladeGeometryModes.sourceSurfaceFile = surfaceFile;
BladeGeometryModes.cfg = cfg;
BladeGeometryModes.geometryTable = geometryTable;
BladeGeometryModes.fitTable = fitTable;
BladeGeometryModes.bladeCache = bladeCache;
BladeGeometryModes.formula = 'R ~= ddelta*Fg + mu*xFg - tau*Fx + kappa*xFx + c';

matFile = fullfile(outDir, 'Step05J_BladeGeometry_Modes_20251222.mat');
geometryCsv = fullfile(outDir, 'Step05J_BladeGeometry_Modes_20251222.csv');
fitCsv = fullfile(outDir, 'Step05J_BladeGeometry_Modes_Fit_20251222.csv');
figFile = fullfile(figDir, 'Step05J_BladeGeometry_Modes_20251222.png');

save(matFile, 'BladeGeometryModes', 'geometryTable', 'fitTable', '-v7.3');
writetable(geometryTable, geometryCsv);
writetable(fitTable, fitCsv);
plot_geometry_modes_local(BladeGeometryModes, figFile);

fprintf('\nStep05J complete.\n');
disp(geometryTable(:, {'bladeId','deltaPriorMm','deltaCorrectionMm','muGapPerXMm','tauMm','kappaScale','rmseBeforeMv','rmseAfterMv','rmseReductionRatio'}));
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, geometryCsv, figFile);

function [A, yRes, obs, pred0, xAll, gapAll, gapIndexAll] = build_blade_regression_local(Yblade, xGrid, trueGaps, coeff, g0Model, deltaRef, effectiveWindow)
Yblade = squeeze(Yblade);
if isvector(Yblade)
    Yblade = reshape(Yblade, numel(xGrid), []);
end
dCoeff = zeros(size(coeff));
for j = 1:size(coeff, 2)
    dCoeff(:, j) = gradient(coeff(:, j), xGrid);
end

A = [];
yRes = [];
obs = [];
pred0 = [];
xAll = [];
gapAll = [];
gapIndexAll = [];
for ig = 1:numel(trueGaps)
    g = trueGaps(ig) + deltaRef;
    if g <= 0
        continue;
    end
    y = Yblade(:, ig);
    F0 = eval_surface_local(coeff, xGrid, g0Model, g);
    Fg = eval_gap_derivative_local(coeff, g);
    Fx = eval_x_derivative_local(dCoeff, g0Model, g);
    keep = effectiveWindow & isfinite(y) & isfinite(F0) & isfinite(Fg) & isfinite(Fx);
    x = xGrid(keep);
    basis = [Fg(keep), x .* Fg(keep), -Fx(keep), x .* Fx(keep), ones(nnz(keep), 1)];
    A = [A; basis]; %#ok<AGROW>
    obs = [obs; y(keep)]; %#ok<AGROW>
    pred0 = [pred0; F0(keep)]; %#ok<AGROW>
    yRes = [yRes; y(keep) - F0(keep)]; %#ok<AGROW>
    xAll = [xAll; x]; %#ok<AGROW>
    gapAll = [gapAll; repmat(g, nnz(keep), 1)]; %#ok<AGROW>
    gapIndexAll = [gapIndexAll; repmat(ig, nnz(keep), 1)]; %#ok<AGROW>
end
end

function beta = scaled_ridge_solve_local(A, y, lambda)
valid = all(isfinite(A), 2) & isfinite(y);
A = A(valid, :);
y = y(valid);
scale = sqrt(mean(A.^2, 1, 'omitnan'));
scale(~isfinite(scale) | scale <= 0) = 1;
As = A ./ scale;
betaScaled = (As.' * As + lambda * eye(size(As, 2))) \ (As.' * y);
beta = betaScaled(:) ./ scale(:);
end

function contribution = mode_contribution_local(A, beta)
terms = A .* beta(:).';
energy = sqrt(mean(terms.^2, 1, 'omitnan'));
if sum(energy) > 0
    contribution = energy ./ sum(energy);
else
    contribution = zeros(1, numel(beta));
end
end

function F = eval_surface_local(coeff, xGrid, g0Model, g)
F = coeff(:, 1) + coeff(:, 2) ./ g + coeff(:, 3) .* log(g ./ g0Model);
F(~all(isfinite(coeff), 2) | ~isfinite(xGrid)) = NaN;
end

function Fg = eval_gap_derivative_local(coeff, g)
Fg = -coeff(:, 2) ./ (g .^ 2) + coeff(:, 3) ./ g;
Fg(~all(isfinite(coeff), 2)) = NaN;
end

function Fx = eval_x_derivative_local(dCoeff, g0Model, g)
Fx = dCoeff(:, 1) + dCoeff(:, 2) ./ g + dCoeff(:, 3) .* log(g ./ g0Model);
Fx(~all(isfinite(dCoeff), 2)) = NaN;
end

function plot_geometry_modes_local(M, figFile)
style = paper_style_local();
T = M.geometryTable;
fig = figure('Name', 'Step05J blade geometry modes', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 18.0, 15.0]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; box on;
bar(T.bladeId, [T.deltaPriorMm, T.deltaCorrectionMm], 'stacked');
xlabel('Blade ID'); ylabel('\Delta g (mm)', 'Interpreter', 'tex');
title('Clearance offset');
legend({'Prior', 'Correction'}, 'Location', 'northwest', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
bar(T.bladeId, T.tiltAngleDeg, 'FaceColor', [0.25 0.45 0.70]);
xlabel('Blade ID'); ylabel('Tilt angle (deg)');
title('Tilt mode');
format_axes_local(gca, style);

nexttile; hold on; box on;
bar(T.bladeId, T.tauMm, 'FaceColor', [0.45 0.35 0.65]);
xlabel('Blade ID'); ylabel('\tau (mm)', 'Interpreter', 'tex');
title('x-registration mode');
format_axes_local(gca, style);

nexttile; hold on; box on;
bar(T.bladeId, T.xScale, 'FaceColor', [0.35 0.55 0.35]);
yline(1, 'k--', 'LineWidth', 0.8);
xlabel('Blade ID'); ylabel('1+\kappa', 'Interpreter', 'tex');
title('Width/scale mode');
format_axes_local(gca, style);

nexttile; hold on; box on;
bar(T.bladeId, [T.rmseBeforeMv, T.rmseAfterMv]);
xlabel('Blade ID'); ylabel('RMSE (mV)');
title('Residual reduction');
legend({'Before', 'After'}, 'Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
contrib = [T.contribDelta, T.contribTilt, T.contribTau, T.contribScale, T.contribOffset];
bar(T.bladeId, contrib, 'stacked');
xlabel('Blade ID'); ylabel('Relative mode energy');
title('Mode contribution');
legend({'\Delta g','\mu x','\tau','\kappa x','c'}, 'Location', 'eastoutside', 'Box', 'off', 'Interpreter', 'tex');
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

