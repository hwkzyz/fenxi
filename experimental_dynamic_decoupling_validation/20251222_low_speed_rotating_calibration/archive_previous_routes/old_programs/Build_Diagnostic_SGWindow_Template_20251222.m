%% Build_Diagnostic_SGWindow_Template_20251222
% Build a non-parametric low-speed template from the same static windows used
% by the 20251222 super-Gaussian Step4 calibration.
%
% This is a diagnostic artifact: it tests whether Step03D fails because the
% current Step01 template uses a different/poorer low-speed pulse selection.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
sg_dir = fullfile('E:\0小论文+程序\0博士期间小论文+程序\7超高斯模型-权重-瞬态', ...
    '程序\直叶片验证\实验验证\超高斯\20251222适配-3号传感间隙变化-改');

target_blade = 1;
analysis_sensors = [1, 2, 3];
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
template_grid_n = 1201;
spline_smoothing = 0.995;

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
if exist(template_dir, 'dir') ~= 7; mkdir(template_dir); end

template_file = fullfile(template_dir, sprintf( ...
    'Template_LowSpeedRotating_B%d_%s_SGWindowDiag_20251222.mat', target_blade, sensor_tag));

Template = struct();
Template.Route = 'diagnostic_sg_step4_window_nonparametric_template';
Template.TargetBlade = target_blade;
Template.SensorIDs = analysis_sensors;
Template.SensorTag = sensor_tag;
Template.LowSpeedCase = '1000rpm无振动';
Template.LowSpeedDataDir = '';
Template.SensorConfigFile = '';
Template.CreatedBy = mfilename;
Template.Settings = struct( ...
    'template_grid_n', template_grid_n, ...
    'spline_smoothing', spline_smoothing, ...
    'source', 'SGCalib_20251222_B1_S*.mat Diagnostics.Sensor.x_raw/v_raw');
Template.Sensor = repmat(struct( ...
    'sensor_id', NaN, ...
    'baseline', NaN, ...
    'threshold', NaN, ...
    'xc', NaN, ...
    'x_grid', [], ...
    'v_grid', [], ...
    'dv_dx', [], ...
    'bin_weight', [], ...
    'x_domain', [], ...
    'selection_mode', '', ...
    'point_count', NaN, ...
    'wide_point_count', NaN, ...
    'lap_count', NaN, ...
    'selected_pulse_indices', [], ...
    'fit_method', '', ...
    'spline_smoothing', NaN), numel(analysis_sensors), 1);

fprintf('\n=== Diagnostic SG-window non-parametric template ===\n');
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    sg_file = fullfile(sg_dir, sprintf('SGCalib_20251222_B%d_S%d.mat', target_blade, sid));
    if exist(sg_file, 'file') ~= 2
        error('Missing SG calibration file: %s', sg_file);
    end
    loaded = load(sg_file, 'Calib', 'Diagnostics');
    C = loaded.Calib.Sensor(1);
    D = loaded.Diagnostics.Sensor(1);

    x_rel = D.x_raw(:) - C.xc;
    v_raw = D.v_raw(:);
    finite_mask = isfinite(x_rel) & isfinite(v_raw);
    x_rel = x_rel(finite_mask);
    v_raw = v_raw(finite_mask);

    x_left = min(x_rel);
    x_right = max(x_rel);
    x_grid = linspace(x_left, x_right, template_grid_n).';
    edges = linspace(x_left, x_right, template_grid_n + 1).';

    bin_id = discretize(x_rel, edges);
    valid_bin = ~isnan(bin_id);
    v_med = accumarray(bin_id(valid_bin), v_raw(valid_bin), [template_grid_n, 1], @median, NaN);
    bin_count = accumarray(bin_id(valid_bin), 1, [template_grid_n, 1], @sum, 0);
    valid_grid = isfinite(v_med) & bin_count >= 3;
    if nnz(valid_grid) < 30
        error('Too few valid SG-window bins for CH%d.', sid);
    end

    x_fit = x_grid(valid_grid);
    v_fit = v_med(valid_grid);
    if exist('csaps', 'file') == 2
        v_grid = csaps(x_fit, v_fit, spline_smoothing, x_grid);
        fit_method = 'csaps smoothing spline over SG Step4 x_raw/v_raw';
    else
        v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
        v_grid = smoothdata(v_fill, 'sgolay', 41);
        fit_method = 'median bins + sgolay fallback over SG Step4 x_raw/v_raw';
    end
    v_grid = v_grid(:);

    Template.Sensor(is).sensor_id = sid;
    Template.Sensor(is).baseline = C.baseline;
    Template.Sensor(is).threshold = D.threshold_strict;
    Template.Sensor(is).xc = C.xc;
    Template.Sensor(is).x_grid = x_grid;
    Template.Sensor(is).v_grid = v_grid;
    Template.Sensor(is).dv_dx = gradient(v_grid, x_grid);
    Template.Sensor(is).bin_weight = max(bin_count(:), 1) ./ max(bin_count);
    Template.Sensor(is).x_domain = [x_left, x_right];
    Template.Sensor(is).selection_mode = 'SG Step4 fixed stable-window strict pulse points';
    Template.Sensor(is).point_count = numel(x_rel);
    Template.Sensor(is).wide_point_count = numel(D.x_wide);
    Template.Sensor(is).lap_count = numel(unique(D.peak_row_indices));
    Template.Sensor(is).selected_pulse_indices = unique(D.peak_row_indices(:));
    Template.Sensor(is).fit_method = fit_method;
    Template.Sensor(is).spline_smoothing = spline_smoothing;

    fprintf('CH%d: baseline %.6f V, xc %.6f mm, rows %d..%d, points %d, domain [%.3f, %.3f] mm, vmax %.4f V\n', ...
        sid, C.baseline, C.xc, Template.Sensor(is).selected_pulse_indices(1), ...
        Template.Sensor(is).selected_pulse_indices(end), numel(x_rel), x_left, x_right, max(v_grid));
end

save(template_file, 'Template', '-v7.3');
fprintf('Saved diagnostic template: %s\n', template_file);
