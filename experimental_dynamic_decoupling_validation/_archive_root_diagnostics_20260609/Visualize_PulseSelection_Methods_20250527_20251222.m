%% Visualize_PulseSelection_Methods_20250527_20251222
% Read-only diagnostic figure for Step03 pulse/domain point selection.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'diagnostics_original_config_visuals', 'pulse_selection_methods');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

cases = build_cases_local(rootDir);
for ic = 1:numel(cases)
    C = cases(ic);
    Template = load_var_local(C.templateFile, 'Template');
    DynamicMap = load_var_local(C.dynamicFile, 'DynamicMap');
    Wmap = DynamicMap.Window(C.windowId);
    for sid = C.sensorIds
        fig = plot_one_sensor_selection_local(C, Template, Wmap, sid);
        stem = sprintf('%s_CH%d_W%02d_pulse_selection', C.tag, sid, C.windowId);
        exportgraphics(fig, fullfile(outDir, [stem, '.png']), 'Resolution', 300);
        exportgraphics(fig, fullfile(outDir, [stem, '.pdf']), 'ContentType', 'vector');
        fprintf('Saved %s\n', fullfile(outDir, [stem, '.png']));
    end
end

function cases = build_cases_local(rootDir)
cases = repmat(struct(), 1, 2);
cases(1).tag = '20250527_B1_S136';
cases(1).label = '20250527: single + hard';
cases(1).sensorIds = [1 3 6];
cases(1).windowId = 1;
cases(1).pulseMode = 'single';
cases(1).domainMode = 'hard';
cases(1).softMarginMm = 0;
cases(1).templateFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S136_GradientXRange030_20250527.mat');
cases(1).dynamicFile = fullfile(rootDir, '20250527_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S136_SlidingWindows_20250527.mat');

cases(2).tag = '20251222_B1_S123';
cases(2).label = '20251222: all + soft';
cases(2).sensorIds = [1 2 3];
cases(2).windowId = 1;
cases(2).pulseMode = 'all';
cases(2).domainMode = 'soft';
cases(2).softMarginMm = 0.20;
cases(2).templateFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'templates', 'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
cases(2).dynamicFile = fullfile(rootDir, '20251222_low_speed_rotating_calibration', ...
    'output', 'dynamic_maps', 'DynamicMap_B1_S123_SlidingWindows_20251222.mat');
end

function val = load_var_local(file, varName)
if exist(file, 'file') ~= 2
    error('Missing file: %s', file);
end
S = load(file, varName);
val = S.(varName);
end

function fig = plot_one_sensor_selection_local(C, Template, Wmap, sid)
Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
if isempty(Tpl) || isempty(D)
    error('Missing CH%d in %s.', sid, C.tag);
end

method = default_method_local(C);
v = D.V(:);
x = D.x_rel(:);
t = D.t(:);
theta = D.theta(:);
threshold = get_threshold_local(Tpl, method);
f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x, 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x, 'pchip', NaN);
xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];

[maskAll, nAll] = isolate_all_pulses_local(v, threshold);
[maskSingle, nSingle] = isolate_main_pulse_local(v, threshold);
if strcmpi(C.pulseMode, 'all')
    maskPulse = maskAll;
else
    maskPulse = maskSingle;
end
maskEffective = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method);
maskDomain = x >= Tpl.x_domain(1) + method.domain_margin_mm & ...
             x <= Tpl.x_domain(2) - method.domain_margin_mm;
maskQuerySafe = x >= xDomain(1) + method.query_guard_mm & ...
                x <= xDomain(2) - method.query_guard_mm;
maskFinal = maskPulse & maskEffective & maskDomain & isfinite(f0) & isfinite(fx0) & isfinite(v) & isfinite(theta);
if strcmpi(C.domainMode, 'hard')
    safeMask = maskFinal & maskQuerySafe;
    if nnz(safeMask) >= 8
        maskFinal = safeMask;
    end
end

wEdge = build_edge_weight_local(t(maskFinal), v(maskFinal), method.weight_floor);
wDomain = build_domain_soft_weight_local(x(maskFinal), xDomain, method.domain_soft_margin_mm, method.weight_floor);
wQuery = build_query_guard_soft_weight_local(x(maskFinal), xDomain, method.query_guard_mm, method.domain_soft_margin_mm, method.weight_floor);
wGradient = build_template_gradient_weight_local(Tpl, x(maskFinal), method.domain_selection_mode, method.weight_floor);
wTotal = max(method.weight_floor, wEdge .* wDomain .* wQuery .* wGradient);
if max(wTotal) > 0
    wTotal = max(method.weight_floor, wTotal ./ max(wTotal));
end

fig = figure('Name', sprintf('%s CH%d selection', C.tag, sid), 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 20, 12]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t - min(t), v, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7); hold on;
plot(t(maskAll) - min(t), v(maskAll), '.', 'Color', [0.15 0.45 0.9], 'MarkerSize', 5);
plot(t(maskSingle) - min(t), v(maskSingle), '.', 'Color', [0.9 0.15 0.1], 'MarkerSize', 5);
yline(threshold, 'k--', 'LineWidth', 0.8);
grid on; box on;
title(sprintf('Pulse mask: all=%d seg, single=%d group', nAll, nSingle), 'FontWeight', 'normal');
xlabel('Time in window (s)');
ylabel('Voltage (V)');
legend({'raw', 'all pulse', 'single pulse', 'threshold'}, 'Location', 'best');

nexttile;
plot(t - min(t), v, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7); hold on;
plot(t(maskEffective) - min(t), v(maskEffective), '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 4);
plot(t(maskFinal) - min(t), v(maskFinal), 'r.', 'MarkerSize', 5);
grid on; box on;
title(sprintf('Final selected: %d/%d points', nnz(maskFinal), numel(maskFinal)), 'FontWeight', 'normal');
xlabel('Time in window (s)');
ylabel('Voltage (V)');
legend({'raw', 'dynamic effective', 'final selected'}, 'Location', 'best');

nexttile;
plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.0); hold on;
scatter(x(maskPulse), v(maskPulse), 8, [0.2 0.5 0.9], 'filled', 'MarkerFaceAlpha', 0.20);
scatter(x(maskFinal), v(maskFinal), 12, [0.9 0.1 0.1], 'filled', 'MarkerFaceAlpha', 0.40);
xline(xDomain(1) + method.query_guard_mm, 'k:', 'LineWidth', 0.8);
xline(xDomain(2) - method.query_guard_mm, 'k:', 'LineWidth', 0.8);
grid on; box on;
title('x-V selected crop', 'FontWeight', 'normal');
xlabel('x_{rel} (mm)');
ylabel('Voltage (V)');
legend({'template', 'pulse mask', 'final selected', 'query-safe boundary'}, 'Location', 'best');

nexttile;
scatter(x(maskFinal), v(maskFinal), 14, wTotal, 'filled'); hold on;
plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 0.8);
colorbar;
grid on; box on;
title(sprintf('Final weights: %s + %s', C.pulseMode, C.domainMode), 'FontWeight', 'normal');
xlabel('x_{rel} (mm)');
ylabel('Voltage (V)');

sgtitle(sprintf('%s | CH%d | window %d', C.label, sid, C.windowId), 'Interpreter', 'none');
end

function method = default_method_local(C)
method = struct();
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.amplitude_limit_mm = 0.50;
method.dx_c_limit_mm = 0.20;
method.query_guard_mm = method.amplitude_limit_mm + method.dx_c_limit_mm + method.coverage_safety_margin_mm;
method.default_sensor_threshold = 0.5;
method.dynamic_effective_mode = 'gradient';
method.dynamic_template_gradient_min_ratio = 0.08;
method.dynamic_time_gradient_min_ratio = 0.15;
method.dynamic_peak_quantile = 85;
method.pulse_selection_mode = C.pulseMode;
method.domain_selection_mode = C.domainMode;
method.domain_soft_margin_mm = C.softMarginMm;
end

function threshold = get_threshold_local(Tpl, method)
if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = method.default_sensor_threshold;
end
end

function [mask, segment_count] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segment_count = numel(starts);
end

function [mask, segment_count] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segment_count = 1;
    return;
end
segment_gap = idx(starts(2:end)) - idx(ends(1:end-1));
large_gap_threshold = max(10, round(0.02 * numel(v)));
pulse_breaks = find(segment_gap > large_gap_threshold);
group_starts = [1; pulse_breaks + 1];
group_ends = [pulse_breaks; numel(starts)];
for ig = 1:numel(group_starts)
    seg_ids = group_starts(ig):group_ends(ig);
    best_seg = seg_ids(1);
    best_peak = -inf;
    for iseg = seg_ids
        seg = idx(starts(iseg):ends(iseg));
        peak_val = max(v(seg));
        if peak_val > best_peak
            best_peak = peak_val;
            best_seg = iseg;
        end
    end
    mask(idx(starts(best_seg)):idx(ends(best_seg))) = true;
end
segment_count = numel(group_starts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end
gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gTplMax = max(gTpl(finite), [], 'omitnan');
if ~isfinite(gTplMax) || gTplMax <= 0
    maskTpl = finite;
else
    maskTpl = gTpl >= method.dynamic_template_gradient_min_ratio * gTplMax;
end
gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if ~isfinite(gTimeMax) || gTimeMax <= 0
    maskTime = false(size(v));
else
    maskTime = gTime >= method.dynamic_time_gradient_min_ratio * gTimeMax;
end
peakQ = min(max(method.dynamic_peak_quantile, 0), 100);
peakLevel = prctile(v(finite), peakQ);
maskPeak = v >= max(threshold, peakLevel);
mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function wEdge = build_edge_weight_local(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_local(x, xDomain, marginMm, floorW)
if marginMm <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMm, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMm, marginMm, floorW)
if queryGuardMm <= 0 || marginMm <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMm - marginMm, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMm, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_local(Tpl, x, domainMode, floorW)
if ~strcmpi(domainMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end
