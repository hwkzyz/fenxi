%% Visualize_Step03_DataSelection_Comparison
% Compare Step03 data-selection rules for the two low-speed rotating cases.
% Raw dynamic data are plotted in light gray; selected points are overlaid.

clear; clc; close all;

rootDir = fileparts(mfilename('fullpath'));
Options.saveFigures = false;  % Set true only when exported PNG/PDF files are needed.
Options.saveSummary = false;  % Set true only when a CSV summary file is needed.
Options.exportPDF = false;    % PNG is much faster for dense selected-point diagnostics.

Cases = [
    struct( ...
        'caseDate', '20250527', ...
        'folder', '20250527_low_speed_rotating_calibration', ...
        'sensorTag', 'S136', ...
        'analysisSensors', [1 3 6], ...
        'windowID', 9)
    struct( ...
        'caseDate', '20251222', ...
        'folder', '20251222_low_speed_rotating_calibration', ...
        'sensorTag', 'S123', ...
        'analysisSensors', [1 2 3], ...
        'windowID', 9)
    ];

Strategies = build_selection_strategies_local();

for ic = 1:numel(Cases)
    C = Cases(ic);
    routeDir = fullfile(rootDir, C.folder);
    templateFile = fullfile(routeDir, 'output', 'templates', sprintf( ...
        'Template_LowSpeedRotating_B1_%s_GradientXRange030_OPRCenterStd_%s.mat', ...
        C.sensorTag, C.caseDate));
    dynamicFile = fullfile(routeDir, 'output', 'dynamic_maps', sprintf( ...
        'DynamicMap_B1_%s_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_%s.mat', ...
        C.sensorTag, C.caseDate));
    resultFile = fullfile(routeDir, 'output', 'identification', sprintf( ...
        'Result_Step03_Main_VPTop3SynchronousWaveform_B1_%s_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_%s.mat', ...
        C.sensorTag, C.caseDate));

    if exist(templateFile, 'file') ~= 2
        error('Template file not found: %s', templateFile);
    end
    if exist(dynamicFile, 'file') ~= 2
        error('Dynamic map file not found: %s', dynamicFile);
    end

    loadedTemplate = load(templateFile, 'Template');
    loadedDynamic = load(dynamicFile, 'DynamicMap');
    Template = loadedTemplate.Template;
    DynamicMap = loadedDynamic.DynamicMap;
    PhaseRef = [];
    if exist(resultFile, 'file') == 2
        loadedResult = load(resultFile, 'Result');
        if C.windowID <= numel(loadedResult.Result.WindowResult)
            PhaseRef = loadedResult.Result.WindowResult(C.windowID).Result;
        end
    end

    if C.windowID > numel(DynamicMap.Window)
        error('Requested window %d exceeds available windows (%d) for %s.', ...
            C.windowID, numel(DynamicMap.Window), C.caseDate);
    end

    figureDir = fullfile(routeDir, 'output', 'figures', 'data_selection_compare');
    summaryDir = fullfile(routeDir, 'output', 'identification');
    if Options.saveFigures && exist(figureDir, 'dir') ~= 7; mkdir(figureDir); end
    if Options.saveSummary && exist(summaryDir, 'dir') ~= 7; mkdir(summaryDir); end

    Wmap = DynamicMap.Window(C.windowID);
    Summary = table();
    Selection = cell(numel(C.analysisSensors), numel(Strategies));
    for is = 1:numel(C.analysisSensors)
        sid = C.analysisSensors(is);
        Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
        D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);
        for im = 1:numel(Strategies)
            method = Strategies(im).method;
            S = select_points_for_strategy_local(D, Tpl, method, PhaseRef, is);
            Selection{is, im} = S;
            Summary = [Summary; make_summary_row_local(C, Wmap, sid, Strategies(im), S)]; %#ok<AGROW>
        end
    end

    summaryFile = fullfile(summaryDir, sprintf( ...
        'Step03_DataSelectionCompare_W%02d_%s_%s.csv', C.windowID, C.sensorTag, C.caseDate));
    fprintf('\n%s %s W%02d selected-point summary:\n', C.caseDate, C.sensorTag, C.windowID);
    disp(Summary(:, {'SensorID','Label','SelectedCount','SelectedFraction','QueryGuardMM'}));
    if Options.saveSummary
        writetable(Summary, summaryFile);
    end

    plot_time_selection_figure_local(C, Wmap, Strategies, Selection, figureDir, Options);
    plot_x_selection_figure_local(C, Wmap, Template, Strategies, Selection, figureDir, Options);
    plot_combined_selection_figure_local(C, Wmap, Template, Strategies, Selection, figureDir, Options);

    if Options.saveFigures || Options.saveSummary
        fprintf('Finished data-selection comparison for %s W%02d:\n', C.caseDate, C.windowID);
        if Options.saveFigures; fprintf('  Figures: %s\n', figureDir); end
        if Options.saveSummary; fprintf('  Summary: %s\n', summaryFile); end
    else
        fprintf('Displayed data-selection comparison for %s W%02d. No files were saved.\n', ...
            C.caseDate, C.windowID);
    end
end

function Strategies = build_selection_strategies_local()
base = default_method_local();
Strategies = repmat(struct('name', '', 'label', '', 'color', [], 'method', base), 6, 1);

Strategies(1).name = 'Main_GradientSingleHardAdaptive';
Strategies(1).label = 'Main';
Strategies(1).color = [0.000, 0.447, 0.741];
Strategies(1).method = base;

Strategies(2).name = 'AllPulse_GradientAllHardAdaptive';
Strategies(2).label = 'All pulses';
Strategies(2).color = [0.850, 0.325, 0.098];
Strategies(2).method = base;
Strategies(2).method.pulse_selection_mode = 'all';

Strategies(3).name = 'SoftDomain_GradientSingleSoftAdaptive';
Strategies(3).label = 'Soft domain';
Strategies(3).color = [0.466, 0.674, 0.188];
Strategies(3).method = base;
Strategies(3).method.domain_selection_mode = 'soft';
Strategies(3).method.domain_soft_margin_mm = 0.20;

Strategies(4).name = 'LegacyEffective_SingleHardAdaptive';
Strategies(4).label = 'No gradient';
Strategies(4).color = [0.494, 0.184, 0.556];
Strategies(4).method = base;
Strategies(4).method.dynamic_effective_mode = 'legacy';

Strategies(5).name = 'FixedGuard_GradientSingleHardFixed';
Strategies(5).label = 'Fixed guard';
Strategies(5).color = [0.929, 0.694, 0.125];
Strategies(5).method = base;
Strategies(5).method.query_guard_mode = 'fixed';
Strategies(5).method.query_guard_mm = 0.20;

Strategies(6).name = 'PhaseSafe_GradientSingleHardAdaptive';
Strategies(6).label = 'Phase-safe';
Strategies(6).color = [0.301, 0.745, 0.933];
Strategies(6).method = base;
Strategies(6).method.phase_safe_expansion = true;
Strategies(6).method.phase_safe_margin_mm = 0.03;
end

function method = default_method_local()
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.amplitude_limit_mm = 0.50;
method.dx_c_limit_mm = 0.35;
method.query_guard_mm = method.amplitude_limit_mm + method.dx_c_limit_mm + method.coverage_safety_margin_mm;
method.query_guard_mode = 'adaptive';
method.query_guard_quantile = 95;
method.query_guard_safety_mm = 0.05;
method.query_guard_min_mm = 0.12;
method.query_guard_max_mm = method.query_guard_mm;
method.domain_selection_mode = 'hard';
method.domain_soft_margin_mm = 0;
method.overshoot_penalty_weight = 100;
method.default_sensor_threshold = 0.5;
method.pulse_selection_mode = 'single';
method.dynamic_effective_mode = 'gradient';
method.dynamic_template_gradient_min_ratio = 0.08;
method.dynamic_time_gradient_min_ratio = 0.15;
method.dynamic_peak_quantile = 85;
method.sensor_eta_limit_mm = 0;
method.sensor_eta_reg_weight_v_per_mm = 0;
method.phase_safe_expansion = false;
method.phase_safe_margin_mm = 0.03;
end

function S = select_points_for_strategy_local(D, Tpl, method, phase_ref, sensor_local_index)
if nargin < 4
    phase_ref = [];
end
if nargin < 5
    sensor_local_index = 1;
end
v_raw = D.V(:);
x_raw = D.x_rel(:);
t_raw = D.t(:);
theta_raw = D.theta(:);
f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_raw, 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x_raw, 'pchip', NaN);

if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = method.default_sensor_threshold;
end

x_domain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
if isfield(Tpl, 'x_domain') && ~isempty(Tpl.x_domain)
    strict_domain = Tpl.x_domain;
else
    strict_domain = x_domain;
end

if strcmpi(method.pulse_selection_mode, 'all')
    [mask_pulse, pulse_segment_count] = isolate_all_pulses_local(v_raw, threshold);
else
    [mask_pulse, pulse_segment_count] = isolate_main_pulse_local(v_raw, threshold);
end
mask_effective = build_dynamic_effective_mask_local(x_raw, v_raw, t_raw, Tpl, threshold, method);
mask_domain = x_raw >= strict_domain(1) + method.domain_margin_mm & ...
              x_raw <= strict_domain(2) - method.domain_margin_mm;
mask_finite = isfinite(f0) & isfinite(fx0) & isfinite(v_raw) & isfinite(theta_raw);
mask_base = mask_pulse & mask_effective & mask_domain & mask_finite;
query_guard_mm = resolve_query_guard_mm_local(Tpl, x_raw, v_raw, mask_base, method);
mask_query_safe = x_raw >= x_domain(1) + query_guard_mm & ...
                  x_raw <= x_domain(2) - query_guard_mm;
mask_phase_safe = false(size(mask_query_safe));
if isfield(method, 'phase_safe_expansion') && method.phase_safe_expansion && ~isempty(phase_ref)
    [mask_phase_safe, x_query] = build_phase_safe_query_mask_local( ...
        x_raw, theta_raw, x_domain, sensor_local_index, phase_ref, method);
    if nnz(mask_base & mask_phase_safe) >= 8
        mask_query_safe = mask_phase_safe;
    end
else
    x_query = NaN(size(x_raw));
end

mask = mask_base;
fallback_level = 0;
if strcmpi(method.domain_selection_mode, 'hard')
    safe_mask = mask & mask_query_safe;
    if nnz(safe_mask) >= 8
        mask = safe_mask;
    end
end
if nnz(mask) < 8
    fallback_level = 1;
    mask = mask_effective & mask_domain & mask_finite;
    if strcmpi(method.domain_selection_mode, 'hard')
        safe_mask = mask & mask_query_safe;
        if nnz(safe_mask) >= 8
            mask = safe_mask;
        end
    end
end

if nnz(mask) > 0
    w_edge = build_edge_weight_local(t_raw(mask), v_raw(mask), method.weight_floor);
    if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
        w_template = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x_raw(mask), 'linear', method.weight_floor);
    else
        w_template = ones(nnz(mask), 1);
    end
    w_domain = build_domain_soft_weight_local(x_raw(mask), x_domain, method.domain_soft_margin_mm, method.weight_floor);
    w_query = build_query_guard_soft_weight_local(x_raw(mask), x_domain, query_guard_mm, method.domain_soft_margin_mm, method.weight_floor);
    w_gradient = build_template_gradient_weight_local(Tpl, x_raw(mask), method.domain_selection_mode, method.weight_floor);
    w_total = max(method.weight_floor, w_edge(:) .* w_template(:) .* w_domain(:) .* w_query(:) .* w_gradient(:));
    if max(w_total) > 0
        w_total = max(method.weight_floor, w_total ./ max(w_total));
    end
else
    w_total = [];
end

S = struct();
S.t = t_raw;
S.x = x_raw;
S.v = v_raw;
S.mask = mask;
S.mask_pulse = mask_pulse;
S.mask_effective = mask_effective;
S.mask_domain = mask_domain;
S.mask_query_safe = mask_query_safe;
S.mask_phase_safe = mask_phase_safe;
S.x_query = x_query;
S.mask_base = mask_base;
S.query_guard_mm = query_guard_mm;
S.threshold = threshold;
S.x_domain = x_domain;
S.selected_count = nnz(mask);
S.raw_count = numel(v_raw);
S.base_count = nnz(mask_base);
S.pulse_count = nnz(mask_pulse);
S.effective_count = nnz(mask_effective);
S.domain_count = nnz(mask_domain);
S.query_safe_selected_count = nnz(mask_query_safe(mask));
S.phase_safe_count = nnz(mask_base & mask_phase_safe);
S.pulse_segment_count = pulse_segment_count;
S.fallback_level = fallback_level;
S.weights = w_total;
end

function [mask_safe, x_query] = build_phase_safe_query_mask_local(x, theta, x_domain, sensor_local_index, phase_ref, method)
A = phase_ref.A_id;
eo = phase_ref.EO_id;
phi = phase_ref.phi_id_wrapped;
dx_c = phase_ref.dx_c_id;
eta = 0;
if isfield(phase_ref, 'sensor_eta_id') && numel(phase_ref.sensor_eta_id) >= sensor_local_index
    eta = phase_ref.sensor_eta_id(sensor_local_index);
end
u_est = A .* sin(eo .* theta(:) + phi);
x_query = x(:) - dx_c - eta - u_est;
margin = method.phase_safe_margin_mm;
mask_safe = x_query >= x_domain(1) + margin & x_query <= x_domain(2) - margin;
mask_safe = reshape(mask_safe, size(x));
x_query = reshape(x_query, size(x));
end

function row = make_summary_row_local(C, Wmap, sid, Strategy, S)
row = table();
row.CaseDate = string(C.caseDate);
row.WindowID = Wmap.window_id;
row.SensorID = sid;
row.Strategy = string(Strategy.name);
row.Label = string(Strategy.label);
row.RawCount = S.raw_count;
row.SelectedCount = S.selected_count;
row.SelectedFraction = S.selected_count / max(S.raw_count, 1);
row.BaseCount = S.base_count;
row.PulseCount = S.pulse_count;
row.EffectiveCount = S.effective_count;
row.DomainCount = S.domain_count;
row.QuerySafeSelectedCount = S.query_safe_selected_count;
row.PhaseSafeCount = S.phase_safe_count;
row.QueryGuardMM = S.query_guard_mm;
row.ThresholdV = S.threshold;
row.PulseSegmentCount = S.pulse_segment_count;
row.FallbackLevel = S.fallback_level;
if isempty(S.weights)
    row.WeightMin = NaN;
    row.WeightMedian = NaN;
    row.WeightMax = NaN;
else
    row.WeightMin = min(S.weights, [], 'omitnan');
    row.WeightMedian = median(S.weights, 'omitnan');
    row.WeightMax = max(S.weights, [], 'omitnan');
end
end

function plot_time_selection_figure_local(C, Wmap, Strategies, Selection, figureDir, Options)
fig = figure('Name', sprintf('%s W%02d time selection', C.caseDate, C.windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 16]);
tiledlayout(fig, numel(C.analysisSensors), numel(Strategies), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(C.analysisSensors)
    sid = C.analysisSensors(is);
    for im = 1:numel(Strategies)
        nexttile;
        S = Selection{is, im};
        t_ms = 1e3 * (S.t - min(S.t));
        plot(t_ms, S.v, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.7); hold on;
        scatter(t_ms(S.mask), S.v(S.mask), 9, Strategies(im).color, 'filled', ...
            'MarkerFaceAlpha', 0.82, 'MarkerEdgeAlpha', 0.82);
        yline(S.threshold, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
        style_axes_local();
        if is == 1
            title(sprintf('%s\n%d/%d pts', Strategies(im).label, S.selected_count, S.raw_count), ...
                'FontWeight', 'normal');
        else
            title(sprintf('%d/%d pts', S.selected_count, S.raw_count), 'FontWeight', 'normal');
        end
        if im == 1
            ylabel(sprintf('CH%d V (V)', sid));
        end
        if is == numel(C.analysisSensors)
            xlabel('Time (ms)');
        end
    end
end

sgtitle(sprintf('%s, %s, window %02d, lap [%d %d]: time-domain selected points', ...
    C.caseDate, C.sensorTag, C.windowID, Wmap.lap_range(1), Wmap.lap_range(end)), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

baseName = sprintf('Step03_DataSelection_Time_W%02d_%s_%s', C.windowID, C.sensorTag, C.caseDate);
if Options.saveFigures
    exportgraphics(fig, fullfile(figureDir, [baseName, '.png']), 'Resolution', 300);
    if Options.exportPDF
        exportgraphics(fig, fullfile(figureDir, [baseName, '.pdf']), 'ContentType', 'vector');
    end
    close(fig);
end
end

function plot_x_selection_figure_local(C, Wmap, Template, Strategies, Selection, figureDir, Options)
fig = figure('Name', sprintf('%s W%02d x selection', C.caseDate, C.windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 16]);
tiledlayout(fig, numel(C.analysisSensors), numel(Strategies), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(C.analysisSensors)
    sid = C.analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    for im = 1:numel(Strategies)
        nexttile;
        S = Selection{is, im};
        scatter(S.x, S.v, 5, [0.80 0.80 0.80], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45); hold on;
        plot(Tpl.x_grid(:), Tpl.v_grid(:), '-', 'Color', [0.20 0.20 0.20], 'LineWidth', 0.8);
        scatter(S.x(S.mask), S.v(S.mask), 11, Strategies(im).color, 'filled', ...
            'MarkerFaceAlpha', 0.82, 'MarkerEdgeAlpha', 0.82);
        xline(S.x_domain(1), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
        xline(S.x_domain(2), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);
        xline(S.x_domain(1) + S.query_guard_mm, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.7);
        xline(S.x_domain(2) - S.query_guard_mm, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.7);
        if any(S.mask)
            xline(min(S.x(S.mask), [], 'omitnan'), '-', 'Color', Strategies(im).color, 'LineWidth', 0.9);
            xline(max(S.x(S.mask), [], 'omitnan'), '-', 'Color', Strategies(im).color, 'LineWidth', 0.9);
        end
        style_axes_local();
        if is == 1
            title(make_panel_title_local(Strategies(im), S, true), 'FontWeight', 'normal');
        else
            title(make_panel_title_local(Strategies(im), S, false), 'FontWeight', 'normal');
        end
        if im == 1
            ylabel(sprintf('CH%d V (V)', sid));
        end
        if is == numel(C.analysisSensors)
            xlabel('x relative to center (mm)');
        end
    end
end

sgtitle(sprintf('%s, %s, window %02d, lap [%d %d]: x-domain selected points', ...
    C.caseDate, C.sensorTag, C.windowID, Wmap.lap_range(1), Wmap.lap_range(end)), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

baseName = sprintf('Step03_DataSelection_XDomain_W%02d_%s_%s', C.windowID, C.sensorTag, C.caseDate);
if Options.saveFigures
    exportgraphics(fig, fullfile(figureDir, [baseName, '.png']), 'Resolution', 300);
    if Options.exportPDF
        exportgraphics(fig, fullfile(figureDir, [baseName, '.pdf']), 'ContentType', 'vector');
    end
    close(fig);
end
end

function plot_combined_selection_figure_local(C, Wmap, Template, Strategies, Selection, figureDir, Options)
fig = figure('Name', sprintf('%s W%02d combined selection', C.caseDate, C.windowID), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 31, 17]);
tiledlayout(fig, numel(C.analysisSensors), numel(Strategies), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(C.analysisSensors)
    sid = C.analysisSensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    for im = 1:numel(Strategies)
        nexttile;
        S = Selection{is, im};
        plot_method_x_panel_local(S, Tpl, Strategies(im));
        if is == 1
            title(make_panel_title_local(Strategies(im), S, true), 'FontWeight', 'normal');
        else
            title(make_panel_title_local(Strategies(im), S, false), 'FontWeight', 'normal');
        end
        if im == 1
            ylabel(sprintf('CH%d V (V)', sid));
        end
        if is == numel(C.analysisSensors)
            xlabel('x relative to center (mm)');
        end
    end
end

sgtitle(sprintf('%s, %s, window %02d, lap [%d %d]: combined selected-point comparison', ...
    C.caseDate, C.sensorTag, C.windowID, Wmap.lap_range(1), Wmap.lap_range(end)), ...
    'FontName', 'Times New Roman', 'FontSize', 9, 'FontWeight', 'normal');

baseName = sprintf('Step03_DataSelection_Combined_W%02d_%s_%s', C.windowID, C.sensorTag, C.caseDate);
if Options.saveFigures
    exportgraphics(fig, fullfile(figureDir, [baseName, '.png']), 'Resolution', 300);
    if Options.exportPDF
        exportgraphics(fig, fullfile(figureDir, [baseName, '.pdf']), 'ContentType', 'vector');
    end
    close(fig);
end
end

function plot_method_x_panel_local(S, Tpl, Strategy)
is_phase_safe = isfield(Strategy.method, 'phase_safe_expansion') && Strategy.method.phase_safe_expansion;
scatter(S.x, S.v, 4, [0.82 0.82 0.82], 'filled', ...
    'MarkerFaceAlpha', 0.28, 'MarkerEdgeAlpha', 0.28); hold on;
plot(Tpl.x_grid(:), Tpl.v_grid(:), '-', 'Color', [0.20 0.20 0.20], 'LineWidth', 0.75);

y_min = min(S.v, [], 'omitnan');
y_max = max(S.v, [], 'omitnan');
y_span = max(y_max - y_min, eps);
y_patch_min = y_min - 0.05 * y_span;
y_patch_max = y_max + 0.05 * y_span;
if any(S.mask) && ~is_phase_safe
    x_sel_min = min(S.x(S.mask), [], 'omitnan');
    x_sel_max = max(S.x(S.mask), [], 'omitnan');
    patch([x_sel_min x_sel_max x_sel_max x_sel_min], ...
        [y_patch_min y_patch_min y_patch_max y_patch_max], Strategy.color, ...
        'FaceAlpha', 0.075, 'EdgeColor', Strategy.color, 'EdgeAlpha', 0.70, ...
        'LineWidth', 1.0);
end

xline(S.x_domain(1), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.65);
xline(S.x_domain(2), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.65);
xline(S.x_domain(1) + S.query_guard_mm, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.75);
xline(S.x_domain(2) - S.query_guard_mm, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.75);
if any(S.mask) && ~is_phase_safe
    xline(min(S.x(S.mask), [], 'omitnan'), '-', 'Color', Strategy.color, 'LineWidth', 1.0);
    xline(max(S.x(S.mask), [], 'omitnan'), '-', 'Color', Strategy.color, 'LineWidth', 1.0);
elseif any(S.mask) && is_phase_safe
    margin = Strategy.method.phase_safe_margin_mm;
    xline(S.x_domain(1) + margin, '-', 'Color', Strategy.color, 'LineWidth', 0.9);
    xline(S.x_domain(2) - margin, '-', 'Color', Strategy.color, 'LineWidth', 0.9);
end
scatter(S.x(S.mask), S.v(S.mask), 7, Strategy.color, 'filled', ...
    'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45);

ylim([y_min - 0.08 * y_span, y_max + 0.08 * y_span]);
xlim([min(S.x_domain(1) - 1.5, min(S.x, [], 'omitnan')), ...
      max(S.x_domain(2) + 1.5, max(S.x, [], 'omitnan'))]);
style_axes_local();
end

function txt = make_panel_title_local(Strategy, S, include_label)
if isfield(Strategy.method, 'phase_safe_expansion') && Strategy.method.phase_safe_expansion
    detail = sprintf('%d/%d pts, xq safe=%d', ...
        S.selected_count, S.raw_count, S.phase_safe_count);
else
    detail = sprintf('%d/%d pts, q=%.3f', ...
        S.selected_count, S.raw_count, S.query_guard_mm);
end
if include_label
    txt = sprintf('%s\n%s', Strategy.label, detail);
else
    txt = detail;
end
end

function style_axes_local()
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8, 'TickDir', 'in', 'LineWidth', 0.7);
end

function x_stat = invert_template_voltage_local(Tpl, v, x_ref)
x_grid = Tpl.x_grid(:);
v_grid = Tpl.v_grid(:);
x_stat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diff_v = v_grid - vv;
    crossing_x = [];
    exact_idx = find(abs(diff_v) <= 1e-10);
    if ~isempty(exact_idx)
        crossing_x = x_grid(exact_idx);
    end
    for k = 1:numel(diff_v)-1
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k+1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k+1) > 0
            continue;
        end
        denom = v_grid(k+1) - v_grid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - v_grid(k)) / denom;
        crossing_x(end+1, 1) = x_grid(k) + alpha * (x_grid(k+1) - x_grid(k)); %#ok<AGROW>
    end
    if isempty(crossing_x)
        [~, idx] = min(abs(diff_v));
        x_stat(i) = x_grid(idx);
    else
        [~, idx] = min(abs(crossing_x - x_ref(i)));
        x_stat(i) = crossing_x(idx);
    end
end
end

function query_guard_mm = resolve_query_guard_mm_local(Tpl, x, v, base_mask, method)
query_guard_mm = method.query_guard_mm;
if ~strcmpi(method.query_guard_mode, 'adaptive') || nnz(base_mask) < 8
    return;
end
x_sel = x(base_mask);
v_sel = v(base_mask);
x_stat = invert_template_voltage_local(Tpl, v_sel, x_sel);
u_app = abs(x_sel(:) - x_stat(:));
u_app = u_app(isfinite(u_app));
if isempty(u_app)
    return;
end
q = min(max(method.query_guard_quantile, 0), 100);
adaptive_guard = prctile(u_app, q) + method.query_guard_safety_mm;
query_guard_mm = min(max(adaptive_guard, method.query_guard_min_mm), method.query_guard_max_mm);
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

g_tpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    g_tpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
g_tpl_max = max(g_tpl(finite), [], 'omitnan');
if ~isfinite(g_tpl_max) || g_tpl_max <= 0
    mask_tpl = finite;
else
    mask_tpl = g_tpl >= method.dynamic_template_gradient_min_ratio * g_tpl_max;
end

g_time = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    g_time(finite) = abs(gradient(v(finite), t(finite)));
end
g_time_max = max(g_time(finite), [], 'omitnan');
if ~isfinite(g_time_max) || g_time_max <= 0
    mask_time = false(size(v));
else
    mask_time = g_time >= method.dynamic_time_gradient_min_ratio * g_time_max;
end

peak_q = min(max(method.dynamic_peak_quantile, 0), 100);
peak_level = prctile(v(finite), peak_q);
mask_peak = v >= max(threshold, peak_level);

mask = finite & (mask_tpl | mask_time | mask_peak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function w_edge = build_edge_weight_local(t, v, floor_w)
if numel(v) < 3 || range(t) <= 0
    w_edge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    w_edge = dv ./ max(dv);
else
    w_edge = ones(size(dv));
end
w_edge = max(floor_w, w_edge);
end

function w_domain = build_domain_soft_weight_local(x, x_domain, margin_mm, floor_w)
if margin_mm <= 0
    w_domain = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
ratio = min(max(dist_to_edge ./ margin_mm, 0), 1);
w_domain = floor_w + (1 - floor_w) .* ratio;
w_domain(~isfinite(w_domain)) = floor_w;
end

function w_query = build_query_guard_soft_weight_local(x, x_domain, query_guard_mm, margin_mm, floor_w)
if query_guard_mm <= 0 || margin_mm <= 0
    w_query = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
soft_start = max(query_guard_mm - margin_mm, 0);
ratio = min(max((dist_to_edge - soft_start) ./ max(margin_mm, eps), 0), 1);
w_query = floor_w + (1 - floor_w) .* ratio;
w_query(~isfinite(w_query)) = floor_w;
end

function w_gradient = build_template_gradient_weight_local(Tpl, x, domain_selection_mode, floor_w)
if ~strcmpi(domain_selection_mode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    w_gradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
g_max = max(g, [], 'omitnan');
if ~isfinite(g_max) || g_max <= 0
    w_gradient = ones(size(x));
    return;
end
w_gradient = floor_w + (1 - floor_w) .* g ./ g_max;
w_gradient(~isfinite(w_gradient)) = floor_w;
end
