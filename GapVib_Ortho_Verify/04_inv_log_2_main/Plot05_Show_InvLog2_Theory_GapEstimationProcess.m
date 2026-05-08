%% Plot 05: how the gap library is used in high-speed gap estimation.
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

rng(20261175 + round(1000 * gTrue), 'twister');
dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
    cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

modelBaseline = table("baseline_1_over_g", "inv_g_linear", "", NaN, 1, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});
modelInv = table("field_basis_inv_log_2", "field_basis", "inv_log_2", NaN, 3, ...
    'VariableNames', {'model_label','response_model','basis_name','parameter','basis_order'});

trustOptions = make_gap_trust_options(cfgCase);
templateBaseline = make_response_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelBaseline, trustOptions);
templateInv = make_response_template_library(gapList, xCell, yCell, ...
    gTrue, cfgAna.xGridN, modelInv, trustOptions);
highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
    templateInv.domain, cfgCase.fitActiveLevel);

staticBaseline = estimate_highspeed_static_gap_raw(highMap, templateBaseline, cfgCase);
staticInv = estimate_highspeed_static_gap_raw(highMap, templateInv, cfgCase);
gapStateBaseline = estimate_gap_init_vib_basis_projected(highMap, templateBaseline, cfgCase, staticBaseline);
gapStateInv = estimate_gap_init_vib_basis_projected(highMap, templateInv, cfgCase, staticInv);

gapGrid = linspace(max(0.05, min(templateInv.gapTrain) - cfgCase.rawGapSearchMargin), ...
    max(templateInv.gapTrain) + cfgCase.rawGapSearchMargin, 121);
[JStaticBase, JProjBase] = compute_gap_cost_curves(highMap, templateBaseline, cfgCase, staticBaseline, gapGrid);
[JStaticInv, JProjInv] = compute_gap_cost_curves(highMap, templateInv, cfgCase, staticInv, gapGrid);

gDisplay = unique([gTrue, staticInv.gHat, gapStateInv.gHat]);
gDisplay = gDisplay(isfinite(gDisplay));
displayColors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

fig = figure('Name', 'inv_log_2 gap estimation process', 'Color', 'w');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
plot(ax1, highMap.x_v(:), highMap.V_a(:), 'k.', 'MarkerSize', 6, 'DisplayName', 'high-speed waveform');
hold(ax1, 'on');
for ig = 1:numel(gDisplay)
    g = gDisplay(ig);
    yPred = eval_gap_template(templateInv, g, highMap.x_v(:) - staticInv.dx0);
    [xSorted, sortIdx] = sort(highMap.x_v(:), 'ascend');
    plot(ax1, xSorted, yPred(sortIdx), '-', 'Color', displayColors(ig, :), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('V_0(t; g = %.4g mm)', g));
end
title(ax1, '(a) Candidate no-vibration waveforms from the library', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax1, 'Spatial coordinate x (mm)', 'FontSize', 9);
ylabel(ax1, 'Voltage (V)', 'FontSize', 9);
legend(ax1, 'Location', 'best');
xline(ax1, templateInv.domain(1), ':', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
xline(ax1, templateInv.domain(2), ':', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');

ax2 = nexttile;
plot(ax2, gapGrid, JStaticBase, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.1, 'DisplayName', 'baseline static cost');
hold(ax2, 'on');
plot(ax2, gapGrid, JStaticInv, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'inv_log_2 static cost');
xline(ax2, gTrue, 'k--', 'LineWidth', 0.9, 'DisplayName', 'true g');
xline(ax2, staticInv.gHat, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1, 'DisplayName', 'inv_log_2 static g0');
title(ax2, '(b) Static gap-matching cost J_{static}(g)', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax2, 'Candidate gap g (mm)', 'FontSize', 9);
ylabel(ax2, 'Mean squared residual', 'FontSize', 9);
legend(ax2, 'Location', 'best');

ax3 = nexttile;
plot(ax3, gapGrid, JProjBase, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.1, 'DisplayName', 'baseline projected cost');
hold(ax3, 'on');
plot(ax3, gapGrid, JProjInv, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'inv_log_2 projected cost');
xline(ax3, gTrue, 'k--', 'LineWidth', 0.9, 'DisplayName', 'true g');
xline(ax3, gapStateInv.gHat, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.1, 'DisplayName', 'inv_log_2 corrected g');
title(ax3, '(c) Projected gap cost J_{proj}(g)', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax3, 'Candidate gap g (mm)', 'FontSize', 9);
ylabel(ax3, 'Residual after vibration-basis projection', 'FontSize', 9);
legend(ax3, 'Location', 'best');

ax4 = nexttile;
freqList = gapStateInv.freq_list(:);
freqList = freqList(1:min(2, numel(freqList)));
Fx = eval_gap_derivative(templateInv, gapStateInv.gHat, highMap.x_v(:) - gapStateInv.dx0);
plot(ax4, highMap.t_v(:), Fx ./ max(abs(Fx)), '-', 'Color', [0.25 0.25 0.25], ...
    'LineWidth', 1.0, 'DisplayName', 'normalized F_x''(g,x(t))');
hold(ax4, 'on');
for ifr = 1:numel(freqList)
    w = 2 * pi * freqList(ifr);
    basisCol = Fx .* sin(w * highMap.t_v(:));
    basisCol = basisCol ./ max(abs(basisCol));
    plot(ax4, highMap.t_v(:), basisCol, '-', 'Color', displayColors(ifr, :), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('F_x'' sin(2\\pi %.0f t)', freqList(ifr)));
end
title(ax4, '(d) Vibration-shaped basis built from F_x''(g,x(t))', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax4, 'Time t (s)', 'FontSize', 9);
ylabel(ax4, 'Normalized amplitude', 'FontSize', 9);
legend(ax4, 'Location', 'best');

ax5 = nexttile;
rStatic = highMap.V_a(:) - eval_gap_template(templateInv, staticInv.gHat, highMap.x_v(:) - staticInv.dx0);
rProj = compute_projected_residual(highMap, templateInv, cfgCase, gapStateInv.gHat, staticInv.dx0, gapStateInv.freq_list(:));
plot(ax5, highMap.t_v(:), rStatic, '.', 'Color', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', 'residual before projection');
hold(ax5, 'on');
plot(ax5, highMap.t_v(:), rProj, '.', 'Color', [0.00 0.45 0.74], 'MarkerSize', 5, ...
    'DisplayName', 'residual after projection');
yline(ax5, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
title(ax5, '(e) Residual before and after vibration-basis projection', 'FontWeight', 'normal', 'FontSize', 9);
xlabel(ax5, 'Time t (s)', 'FontSize', 9);
ylabel(ax5, 'Residual (V)', 'FontSize', 9);
legend(ax5, 'Location', 'best');

ax6 = nexttile;
bar(ax6, [ ...
    mean(rStatic .^ 2), ...
    mean(rProj .^ 2), ...
    mean((highMap.V_a(:) - eval_gap_template(templateInv, gapStateInv.gHat, highMap.x_v(:) - staticInv.dx0)) .^ 2)]);
xticklabels(ax6, {'static residual', 'projected residual', 'corrected-gap residual'});
ylabel(ax6, 'Mean squared residual', 'FontSize', 9);
title(ax6, '(f) Residual-energy reduction along the process', 'FontWeight', 'normal', 'FontSize', 9);
xtickangle(ax6, 15);

apply_inv_log_2_figure_style(ax1, 170, 160);

function [JStatic, JProj] = compute_gap_cost_curves(highMap, templateLib, cfgCase, staticState, gapGrid)
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
    B = build_vib_basis(Fx, t(valid), freqList);
    if isempty(B)
        rProj = r;
    else
        coef = B \ r;
        rProj = r - B * coef;
    end
    JProj(ig) = mean(rProj .^ 2);
end
end

function rProj = compute_projected_residual(highMap, templateLib, cfgCase, g, dx0, freqList)
x = highMap.x_v(:);
t = highMap.t_v(:);
V = highMap.V_a(:);
F = eval_gap_template(templateLib, g, x - dx0);
Fx = eval_gap_derivative(templateLib, g, x - dx0);
valid = isfinite(F) & isfinite(Fx) & isfinite(V);
r = V(valid) - F(valid);
B = build_vib_basis(Fx(valid), t(valid), freqList);
if isempty(B)
    rFit = r;
else
    coef = B \ r;
    rFit = r - B * coef;
end
rProj = nan(size(V));
rProj(valid) = rFit;
end

function B = build_vib_basis(Fx, t, freqList)
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
