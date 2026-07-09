%% Single waveform comparison: direct bundle vs no-direct bundle
% A compact figure for visually comparing how the direct bundle point set
% affects the fitted Step07J gap_tilt waveform.
%
% Optional environment variables:
%   STEP07J_BUNDLE_SINGLE_WINDOW   window index, default direct best window
%   STEP07J_BUNDLE_SINGLE_SENSOR   sensor id, default first analysis sensor

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07j_bundle_waveform_compare');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

directFile = fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', ...
    C0.dataset, C0.caseTag));
nodirectFile = fullfile(outDir, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt_nodirect.mat', ...
    C0.dataset, C0.caseTag));
Sd = load(directFile, 'Result');
Sn = load(nodirectFile, 'Result');
Direct = Sd.Result;
NoDirect = Sn.Result;

windowId = parse_positive_integer_env_local('STEP07J_BUNDLE_SINGLE_WINDOW', 0);
if windowId <= 0
    windowId = Direct.BestWindowIndex;
end
sensorId = parse_positive_integer_env_local('STEP07J_BUNDLE_SINGLE_SENSOR', 0);
if sensorId <= 0
    sensorId = Direct.WindowResult(windowId).bundle.sensorIds(1);
end

wrD = Direct.WindowResult(windowId);
wrN = NoDirect.WindowResult(windowId);
PD = prepare_panel_local(wrD.bundle, wrD.modelFits.gap_tilt, sensorId);
PN = prepare_panel_local(wrN.bundle, wrN.modelFits.gap_tilt, sensorId);

style = paper_style_local();
fig = figure('Name', sprintf('Step07J bundle waveform single W%d CH%d', windowId, sensorId), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.0, 10.5]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

yWave = padded_limits_local([PD.vObs; PD.vPred; PN.vObs; PN.vPred], 0.06);
yResid = padded_limits_local([PD.residual; PN.residual], 0.12);
xLimits = padded_limits_local([PD.x; PN.x], 0.04);

draw_wave_panel_local(nexttile(1), PD, style.blue, 'Direct bundle', xLimits, yWave);
draw_wave_panel_local(nexttile(2), PN, style.orange, 'No direct bundle', xLimits, yWave);
draw_residual_panel_local(nexttile(3), PD, style.blue, xLimits, yResid);
draw_residual_panel_local(nexttile(4), PN, style.orange, xLimits, yResid);

sgtitle(sprintf('%s %s W%02d CH%d: direct bundle vs no-direct bundle', ...
    C0.dataset, C0.caseTag, windowId, sensorId), ...
    'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 11);

figBase = sprintf('Step07J_BundleWaveformSingle_%s_%s_W%02d_CH%d', ...
    C0.dataset, C0.caseTag, windowId, sensorId);
pngFile = fullfile(figDir, [figBase, '.png']);
pdfFile = fullfile(figDir, [figBase, '.pdf']);
csvFile = fullfile(figDir, [figBase, '_Curves.csv']);
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
CurveTable = struct2table([curve_rows_local(PD, "direct_bundle"); ...
    curve_rows_local(PN, "no_direct_bundle")]);
writetable(CurveTable, csvFile);

fprintf('\nSaved single waveform comparison:\n  %s\n  %s\n  %s\n', pngFile, pdfFile, csvFile);
fprintf('Direct bundle: %d points, sensor RMSE %.2f mV\n', PD.pointCount, PD.sensorRmseMv);
fprintf('No-direct bundle: %d points, sensor RMSE %.2f mV\n', PN.pointCount, PN.sensorRmseMv);

function P = prepare_panel_local(bundle, fit, sensorId)
idxSensor = find(bundle.sensorIds == sensorId, 1, 'first');
if isempty(idxSensor)
    error('CH%d not found in bundle.', sensorId);
end
idx = find(bundle.sensorIndex == idxSensor);
x = bundle.X(idx) - fit.dxMm - fit.uMm(idx);
vObs = bundle.V(idx);
vPred = fit.VPred(idx);
w = bundle.W(idx);
valid = isfinite(x) & isfinite(vObs) & isfinite(vPred) & isfinite(w);
x = x(valid);
vObs = vObs(valid);
vPred = vPred(valid);
w = w(valid);
[x, order] = sort(x(:));
vObs = vObs(order);
vPred = vPred(order);
w = w(order);
residual = vObs - vPred;
P = struct();
P.x = x;
P.vObs = vObs;
P.vPred = vPred;
P.residual = residual;
P.pointCount = numel(x);
P.sensorRmseMv = sqrt(mean(residual.^2, 'omitnan'));
P.sensorWeightedRmseMv = sqrt(sum(max(w, 0.05) .* residual.^2, 'omitnan') / max(numel(residual), 1));
P.fitEO = fit.EO;
P.fitFreqHz = fit.freqHz;
P.fitAmpMm = fit.amplitudeMm;
end

function draw_wave_panel_local(ax, P, color, titlePrefix, xLimits, yLimits)
axes(ax); hold on; grid on; box on;
scatter(P.x, P.vObs, 6, [0.42, 0.42, 0.42], 'filled', ...
    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25, ...
    'DisplayName', sprintf('measured, n=%d', P.pointCount));
plot(P.x, P.vPred, '-', 'Color', color, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('gap\\_tilt prediction, EO%d', P.fitEO));
xlim(xLimits); ylim(yLimits);
title(sprintf('%s | A %.3f mm | sensor RMSE %.2f mV', ...
    titlePrefix, P.fitAmpMm, P.sensorRmseMv), 'Interpreter', 'none');
ylabel('Voltage (mV)');
legend('Location', 'south', 'Box', 'off', 'Interpreter', 'none');
format_axes_local(ax);
end

function draw_residual_panel_local(ax, P, color, xLimits, yLimits)
axes(ax); hold on; grid on; box on;
plot(xLimits, [0, 0], '-', 'Color', [0.15, 0.15, 0.15], 'LineWidth', 0.8);
scatter(P.x, P.residual, 6, color, 'filled', ...
    'MarkerFaceAlpha', 0.28, 'MarkerEdgeAlpha', 0.28);
xlim(xLimits); ylim(yLimits);
xlabel('Equivalent x after vibration compensation (mm)', 'Interpreter', 'none');
ylabel('Residual (mV)');
format_axes_local(ax);
end

function limits = padded_limits_local(v, padFrac)
v = v(:);
v = v(isfinite(v));
if isempty(v)
    limits = [-1, 1];
    return;
end
lo = min(v); hi = max(v);
if hi <= lo
    span = max(abs(lo), 1);
    limits = [lo - 0.1 * span, hi + 0.1 * span];
else
    span = hi - lo;
    limits = [lo - padFrac * span, hi + padFrac * span];
end
end

function rows = curve_rows_local(P, caseName)
n = numel(P.x);
rows = repmat(struct('case_name', string(caseName), 'point_index', NaN, ...
    'x_mm', NaN, 'measured_mV', NaN, 'predicted_mV', NaN, ...
    'residual_mV', NaN), n, 1);
for i = 1:n
    rows(i).point_index = i;
    rows(i).x_mm = P.x(i);
    rows(i).measured_mV = P.vObs(i);
    rows(i).predicted_mV = P.vPred(i);
    rows(i).residual_mV = P.residual(i);
end
end

function style = paper_style_local()
style = struct();
style.blue = [0.000, 0.447, 0.698];
style.orange = [0.835, 0.369, 0.000];
end

function format_axes_local(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 8.5, ...
    'LineWidth', 0.8, 'TickDir', 'out');
end

function value = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value) || value < 0
    error('%s must be a nonnegative number.', name);
end
value = floor(value);
end
