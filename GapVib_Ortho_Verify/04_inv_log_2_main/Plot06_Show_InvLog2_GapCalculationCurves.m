%% Plot 06: dedicated gap calculation curves for the inv_log_2 method.
close all;
clc;

ctx = load_inv_log_2_project_context(fileparts(mfilename('fullpath')));
gapList = ctx.gapList;
xCell = ctx.xCell;
yCell = ctx.yCell;
cfgAna = ctx.cfgAna;

if ~exist('gTrue', 'var') || isempty(gTrue)
    if any(abs(gapList - 0.2) < 1e-12)
        gTrue = 0.2;
    else
        gTrue = ctx.cfg.g_holdout;
    end
end

cfgCase = cfgAna;
cfgCase.g_holdout = gTrue;
cfgCase.A_true = [0.25, 0.15];
cfgCase.f_true = [500, 1300];
cfgCase.phi_true = [pi/4, -pi/3];
cfgCase.noiseMode = 'noise_ratio';
cfgCase.snrDb = 10;

templateTruth = build_gap_template_library(gapList, xCell, yCell, NaN, cfgAna.xGridN);
cfgCase.noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, ...
    gapList, xCell, yCell, templateTruth.domain, cfgCase.snrDb);

idxTrue = find(abs(gapList - gTrue) < 1e-12, 1);
if isempty(idxTrue)
    error('Holdout gap %.6g mm not found in data file.', gTrue);
end
FTrue = griddedInterpolant(xCell{idxTrue}, yCell{idxTrue}, 'pchip', 'nearest');

rng(20261205 + round(1000 * gTrue), 'twister');
dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
    cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

modelInv = table("field_basis_inv_log_2", "field_basis", "inv_log_2", NaN, 3, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
trustOptions = make_gap_trust_options(cfgCase);
templateInv = make_response_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelInv, trustOptions);
highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
    templateInv.domain, cfgCase.fitActiveLevel);

staticInv = estimate_highspeed_static_gap_raw(highMap, templateInv, cfgCase);
gapStateInv = estimate_gap_init_vib_basis_projected(highMap, templateInv, cfgCase, staticInv);

coarseGapGrid = staticInv.coarse_gap_grid(:);
coarseDxGrid = staticInv.coarse_dx_grid(:);
fineGapGrid = staticInv.fine_gap_grid(:);
fineDxGrid = staticInv.fine_dx_grid(:);

sampleIdx = pick_evenly_spaced_indices(numel(highMap.x_v(:)), get_cfg_field_local(cfgCase, 'gapStaticTargetSamples', 2500));
xCoarse = highMap.x_v(sampleIdx);
yCoarse = highMap.V_a(sampleIdx);
JCoarse2D = eval_static_gap_cost_grid_local(xCoarse, yCoarse, templateInv, coarseGapGrid, coarseDxGrid);
JFine2D = eval_static_gap_cost_grid_local(highMap.x_v(:), highMap.V_a(:), templateInv, fineGapGrid, fineDxGrid);

JCoarseGap = min(JCoarse2D, [], 2);
JFineGap = min(JFine2D, [], 2);

gapGridDense = linspace(max(0.05, staticInv.gHat - 0.12), ...
    min(max(templateInv.gapTrain) + cfgCase.rawGapSearchMargin, staticInv.gHat + 0.12), 91);
[JStaticDense, JProjDense] = compute_gap_cost_curves_local(highMap, templateInv, cfgCase, staticInv, gapGridDense);

[~, idxStaticBest] = min(JStaticDense);
[~, idxProjBest] = min(JProjDense);
gStaticCurve = gapGridDense(idxStaticBest);
gProjCurve = gapGridDense(idxProjBest);

fig = figure('Name', 'inv_log_2 gap calculation curves', 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
trainGaps = templateInv.gapTrain(:);
trainColors = turbo(numel(trainGaps));
for ig = 1:numel(trainGaps)
    idxTrain = find(abs(gapList(:) - trainGaps(ig)) < 1e-12, 1);
    plot(ax1, xCell{idxTrain}(:), yCell{idxTrain}(:), '-', ...
        'Color', trainColors(ig, :), 'LineWidth', 0.9, ...
        'DisplayName', sprintf('library g = %.1f mm', trainGaps(ig)));
end
plot(ax1, xCell{idxTrue}(:), yCell{idxTrue}(:), 'k--', ...
    'LineWidth', 1.2, 'DisplayName', sprintf('held-out truth g = %.1f mm', gTrue));
xline(ax1, templateInv.domain(1), ':', 'Color', [0.15 0.15 0.15], ...
    'LineWidth', 0.9, 'DisplayName', 'trusted fitting window');
xline(ax1, templateInv.domain(2), ':', 'Color', [0.15 0.15 0.15], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
title(ax1, '(a) Static gap library', ...
    'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'eastoutside');

ax2 = nexttile;
hold(ax2, 'on');
xGrid = templateInv.xGrid(:);
[~, nearestIdx] = min(abs(trainGaps - gProjCurve));
nearestGap = trainGaps(nearestIdx);
yNearest = eval_gap_template(templateInv, nearestGap, xGrid);
yStatic = eval_gap_template(templateInv, staticInv.gHat, xGrid);
yProj = eval_gap_template(templateInv, gapStateInv.gHat, xGrid);
yTruthOnGrid = interp1(xCell{idxTrue}(:), yCell{idxTrue}(:), xGrid, 'pchip');
staticTemplateRmse = sqrt(mean((yStatic(:) - yTruthOnGrid(:)) .^ 2));
projectedTemplateRmse = sqrt(mean((yProj(:) - yTruthOnGrid(:)) .^ 2));
plot(ax2, xGrid, yNearest, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('nearest library g = %.1f mm', nearestGap));
plot(ax2, xGrid, yStatic, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, ...
    'DisplayName', sprintf('static-template solution g = %.4f mm', staticInv.gHat));
plot(ax2, xGrid, yProj, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, ...
    'DisplayName', sprintf('projected-cost solution g = %.4f mm', gapStateInv.gHat));
plot(ax2, xCell{idxTrue}(:), yCell{idxTrue}(:), 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'held-out truth');
title(ax2, '(b) Reconstructed non-library gap waveform', ...
    'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax2, 'Voltage (V)', 'FontSize', 9);
legend(ax2, 'Location', 'best');
text(ax2, 0.04, 0.08, sprintf('static RMSE = %.3g V\nprojected RMSE = %.3g V', ...
    staticTemplateRmse, projectedTemplateRmse), ...
    'Units', 'normalized', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.75 0.75 0.75]);

ax3 = nexttile;
hold(ax3, 'on');
xAligned = highMap.x_v(:) - staticInv.dx0;
plot(ax3, xAligned, highMap.V_a(:), 'k.', 'MarkerSize', 5, ...
    'DisplayName', 'high-speed samples');
[xSorted, sortIdx] = sort(xAligned, 'ascend');
VStatic0 = eval_gap_template(templateInv, staticInv.gHat, xAligned);
VProj0 = eval_gap_template(templateInv, gapStateInv.gHat, xAligned);
VNearest0 = eval_gap_template(templateInv, nearestGap, xAligned);
plot(ax3, xSorted, VNearest0(sortIdx), '-', 'Color', [0.65 0.65 0.65], ...
    'LineWidth', 1.0, 'DisplayName', 'nearest library waveform');
plot(ax3, xSorted, VStatic0(sortIdx), '-', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.3, 'DisplayName', 'static-template solved waveform');
plot(ax3, xSorted, VProj0(sortIdx), '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.4, 'DisplayName', 'projected-cost solved waveform');
title(ax3, '(c) High-speed data and candidate waveforms', ...
    'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax3, 'Aligned spatial coordinate x - dx (mm)', 'FontSize', 9);
ylabel(ax3, 'Voltage (V)', 'FontSize', 9);
legend(ax3, 'Location', 'best');

ax4 = nexttile;
plot(ax4, gapGridDense, JStaticDense, '--', 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.1, 'DisplayName', 'static cost J_static(g)');
hold(ax4, 'on');
plot(ax4, gapGridDense, JProjDense, '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.5, 'DisplayName', 'projected cost J_proj(g)');
xline(ax4, gTrue, 'k--', 'LineWidth', 0.9, ...
    'DisplayName', sprintf('true g = %.4f mm', gTrue));
xline(ax4, gStaticCurve, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1, ...
    'DisplayName', sprintf('static-template minimum = %.4f mm', gStaticCurve));
xline(ax4, gProjCurve, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.1, ...
    'DisplayName', sprintf('projected-cost minimum = %.4f mm', gProjCurve));
title(ax4, '(d) Gap calculation curve', ...
    'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax4, 'Candidate gap g (mm)', 'FontSize', 9);
ylabel(ax4, 'Residual cost', 'FontSize', 9);
legend(ax4, 'Location', 'best');

apply_inv_log_2_figure_style(ax1, 190, 125);

function J2D = eval_static_gap_cost_grid_local(x, y, templateLib, gapGrid, dxGrid)
J2D = zeros(numel(gapGrid), numel(dxGrid));
for ig = 1:numel(gapGrid)
    for id = 1:numel(dxGrid)
        yPred = eval_gap_template(templateLib, gapGrid(ig), x - dxGrid(id));
        valid = isfinite(yPred) & isfinite(y);
        res = y(valid) - yPred(valid);
        J2D(ig, id) = mean(res .^ 2);
    end
end
end

function [JStatic, JProj] = compute_gap_cost_curves_local(highMap, templateLib, cfgCase, staticState, gapGrid)
x = highMap.x_v(:);
t = highMap.t_v(:);
V = highMap.V_a(:);
dx0 = staticState.dx0;
freqState = estimate_gap_init_vib_basis_projected(highMap, templateLib, cfgCase, staticState);
freqList = freqState.freq_list(:);

JStatic = zeros(numel(gapGrid), 1);
JProj = zeros(numel(gapGrid), 1);
for ig = 1:numel(gapGrid)
    g = gapGrid(ig);
    F = eval_gap_template(templateLib, g, x - dx0);
    valid = isfinite(F) & isfinite(V);
    r = V(valid) - F(valid);
    JStatic(ig) = mean(r .^ 2);

    Fx = eval_gap_derivative(templateLib, g, x(valid) - dx0);
    B = build_vib_basis_local(Fx, t(valid), freqList);
    coef = B \ r;
    rProj = r - B * coef;
    JProj(ig) = mean(rProj .^ 2);
end
end

function B = build_vib_basis_local(Fx, t, freqList)
cols = {Fx};
for ifr = 1:numel(freqList)
    w = 2 * pi * freqList(ifr);
    cols{end + 1} = Fx .* sin(w * t); %#ok<AGROW>
    cols{end + 1} = Fx .* cos(w * t); %#ok<AGROW>
end
B = [cols{:}];
end

function noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, gapList, xCell, yCell, domain, snrDb)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
dataClean = simulate_rotating_waveform_from_template(Fx, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    0, domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, 'noise_ratio');
signal = dataClean.V_clean(:) - min(dataClean.V_clean(:));
signalRms = sqrt(mean(signal .^ 2));
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, ...
    linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatio = signalRms ./ (templateRange .* 10^(snrDb / 20));
end

function val = get_cfg_field_local(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
