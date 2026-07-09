function Summary = Diagnose_OPRFlow_SG_vs_Foundation_20251222()
%DIAGNOSE_OPRFLOW_SG_VS_FOUNDATION_20251222
% Compare the OPR/angle/timing reference used by the current foundation
% Step04/Step05 route against the super-Gaussian route that identifies EO14.
%
% The script is read-only. It writes diagnostic CSV tables under
% output/diagnostics/opr_flow_sg_vs_foundation_20251222.

cfg = BTTProjectConfig_20251222();
route_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'opr_flow_sg_vs_foundation_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

sg_root = resolve_sg_root_local(route_dir);
addpath(sg_root);

target_blade = 1;
sensors = [1 2 3];

F = load(fullfile(cfg.step01_output_dir, 'Sensor_Config_20251222.mat'), ...
    'Sensor_Config');
FConfig = F.Sensor_Config;

source_file = fullfile(cfg.step04_output_dir, ...
    'LowSpeed_Template_SourceData_20251222.mat');
if ~isfile(source_file)
    error('Missing foundation Step04 source cache: %s', source_file);
end
loaded_source = load(source_file, 'LowSpeedTemplateSourceData', 'point_cloud');
FCase = loaded_source.LowSpeedTemplateSourceData;
FPointCloud = loaded_source.point_cloud;

template_file = fullfile(cfg.step04_output_dir, ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');
if isfile(template_file)
    loaded_template = load(template_file, 'Template');
    FTemplate = loaded_template.Template;
else
    FTemplate = [];
end

sg_config_file = fullfile(sg_root, 'output', ...
    '1000rpm无振动_reference', 'Sensor_Config_20251222.mat');
if ~isfile(sg_config_file)
    error('Missing SG Sensor_Config: %s', sg_config_file);
end
S = load(sg_config_file, 'Sensor_Config');
SGConfig = S.Sensor_Config;

base_cfg = Get_20251222_BTT_Config();
SGCase = Extract_BTT_Features_20251222( ...
    base_cfg.low_speed_case, base_cfg, struct('sensor_ids', sensors));

AngleTable = compare_angle_and_opr_offsets_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, sensors, target_blade);
EventOffsetTable = build_opr_event_offset_table_local(cfg, FCase, SGCase);
PointCloudTable = compare_point_cloud_centers_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, FPointCloud, FTemplate, ...
    sensors, target_blade);
ArrivalTable = compare_arrival_x_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, sensors, target_blade);
DynamicOPRTable = compare_dynamic_opr_products_local(cfg, sg_root);
DynamicEventOffsetTable = build_dynamic_opr_event_offset_table_local(cfg, sg_root);

writetable(AngleTable, fullfile(out_dir, ...
    'LowSpeed_OPR_Angle_Offset_B1_S123_20251222.csv'));
writetable(EventOffsetTable, fullfile(out_dir, ...
    'LowSpeed_OPR_EventCenterMinusRising_20251222.csv'));
writetable(PointCloudTable, fullfile(out_dir, ...
    'LowSpeed_PointCloud_Remap_By_OPRFlow_B1_S123_20251222.csv'));
writetable(ArrivalTable, fullfile(out_dir, ...
    'LowSpeed_ArrivalX_By_OPRFlow_B1_S123_20251222.csv'));
writetable(DynamicOPRTable, fullfile(out_dir, ...
    'Dynamic_OPR_Product_Offset_20251222.csv'));
if ~isempty(DynamicEventOffsetTable)
    writetable(DynamicEventOffsetTable, fullfile(out_dir, ...
        'Dynamic_OPR_EventCenterMinusRising_20251222.csv'));
end

Summary = struct();
Summary.OutputDir = out_dir;
Summary.SGRoot = sg_root;
Summary.AngleTable = AngleTable;
Summary.EventOffsetTable = EventOffsetTable;
Summary.PointCloudTable = PointCloudTable;
Summary.ArrivalTable = ArrivalTable;
Summary.DynamicOPRTable = DynamicOPRTable;
Summary.DynamicEventOffsetTable = DynamicEventOffsetTable;

fprintf('\n=== OPR-flow diagnostic: foundation vs SG success route ===\n');
fprintf('SG root:\n  %s\n', sg_root);
fprintf('Output dir:\n  %s\n\n', out_dir);
disp(AngleTable);
disp(EventOffsetTable);
disp(PointCloudTable);
disp(DynamicOPRTable);
if ~isempty(DynamicEventOffsetTable)
    disp(DynamicEventOffsetTable);
end
print_diagnostic_findings_local(AngleTable, PointCloudTable, DynamicOPRTable);
end


function sg_root = resolve_sg_root_local(route_dir)
env_root = strtrim(getenv('SG_SUCCESS_ROOT'));
if ~isempty(env_root) && isfolder(env_root)
    sg_root = env_root;
    return;
end

validation_root = fileparts(route_dir);
repo_root = fileparts(validation_root);
common_root = fileparts(repo_root);
candidate = fullfile(common_root, ...
    '7超高斯模型-权重-瞬态', '程序', '直叶片验证', '实验验证', ...
    '超高斯', '20251222适配-3号传感间隙变化-改');
if isfolder(candidate)
    sg_root = candidate;
    return;
end

error(['Cannot locate SG success route. Set SG_SUCCESS_ROOT to the folder ' ...
    'that contains SGCalib_20251222_B1_S1.mat.']);
end


function T = compare_angle_and_opr_offsets_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, sensors, blade_id)
f_opr = FCase.opr_times(:);
sg_opr = SGCase.opr_times(:);
n = min(numel(f_opr), numel(sg_opr));
epr = cfg.blades_num;
spd_deg_s = 360 ./ max(sg_opr((epr + 1):n) - sg_opr(1:(n - epr)), eps);
dt = f_opr(1:(n - epr)) - sg_opr(1:(n - epr));
opr_shift_deg = dt .* spd_deg_s;
event_offset = nan(epr, 1);
event_offset_iqr = nan(epr, 1);
for event_id = 1:epr
    idx = event_id:epr:numel(opr_shift_deg);
    vals = opr_shift_deg(idx);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    event_offset(event_id) = median(vals, 'omitnan');
    event_offset_iqr(event_id) = iqr(vals);
end

rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'dominant_previous_opr_event', NaN, ...
    'foundation_angle_deg', NaN, ...
    'sg_angle_deg', NaN, ...
    'foundation_minus_sg_angle_deg', NaN, ...
    'event_center_minus_rising_deg_median', NaN, ...
    'event_center_minus_rising_deg_iqr', NaN, ...
    'event_center_minus_rising_mm_median', NaN, ...
    'event_specific_compensation_residual_deg', NaN, ...
    'event_specific_compensation_residual_mm', NaN), numel(sensors), 1);

for i = 1:numel(sensors)
    sid = sensors(i);
    f_angle = FConfig.Standard_Relative_Angles_OPRCenter(sid, blade_id);
    sg_angle = SGConfig.Standard_Relative_Angles(sid, blade_id);
    delta_angle = f_angle - sg_angle;
    prev_event = dominant_previous_opr_event_local( ...
        FCase, FConfig, sid, blade_id, cfg.blades_num);
    event_shift = event_offset(prev_event);
    event_shift_iqr = event_offset_iqr(prev_event);
    residual_deg = delta_angle + event_shift;

    rows(i).sensor_id = sid;
    rows(i).blade_id = blade_id;
    rows(i).dominant_previous_opr_event = prev_event;
    rows(i).foundation_angle_deg = f_angle;
    rows(i).sg_angle_deg = sg_angle;
    rows(i).foundation_minus_sg_angle_deg = delta_angle;
    rows(i).event_center_minus_rising_deg_median = event_shift;
    rows(i).event_center_minus_rising_deg_iqr = event_shift_iqr;
    rows(i).event_center_minus_rising_mm_median = ...
        event_shift * pi / 180 * cfg.r_tip_mm;
    rows(i).event_specific_compensation_residual_deg = residual_deg;
    rows(i).event_specific_compensation_residual_mm = ...
        residual_deg * pi / 180 * cfg.r_tip_mm;
end
T = struct2table(rows);
end


function T = build_opr_event_offset_table_local(cfg, FCase, SGCase)
f_opr = FCase.opr_times(:);
sg_opr = SGCase.opr_times(:);
n = min(numel(f_opr), numel(sg_opr));
epr = cfg.blades_num;
spd_deg_s = 360 ./ max(sg_opr((epr + 1):n) - sg_opr(1:(n - epr)), eps);
dt = f_opr(1:(n - epr)) - sg_opr(1:(n - epr));
opr_shift_deg = dt .* spd_deg_s;

rows = repmat(struct( ...
    'opr_event_in_rev', NaN, ...
    'center_minus_rising_time_us_median', NaN, ...
    'center_minus_rising_time_us_iqr', NaN, ...
    'center_minus_rising_deg_median', NaN, ...
    'center_minus_rising_deg_iqr', NaN, ...
    'center_minus_rising_mm_median', NaN, ...
    'sample_count', NaN), epr, 1);
for event_id = 1:epr
    idx = event_id:epr:numel(opr_shift_deg);
    dt_vals = dt(idx);
    deg_vals = opr_shift_deg(idx);
    finite = isfinite(dt_vals) & isfinite(deg_vals);
    dt_vals = dt_vals(finite);
    deg_vals = deg_vals(finite);
    rows(event_id).opr_event_in_rev = event_id;
    rows(event_id).sample_count = numel(deg_vals);
    if isempty(deg_vals)
        continue;
    end
    rows(event_id).center_minus_rising_time_us_median = median(dt_vals) * 1e6;
    rows(event_id).center_minus_rising_time_us_iqr = iqr(dt_vals) * 1e6;
    rows(event_id).center_minus_rising_deg_median = median(deg_vals);
    rows(event_id).center_minus_rising_deg_iqr = iqr(deg_vals);
    rows(event_id).center_minus_rising_mm_median = ...
        rows(event_id).center_minus_rising_deg_median * pi / 180 * cfg.r_tip_mm;
end
T = struct2table(rows);
end


function T = compare_point_cloud_centers_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, point_cloud, Template, ...
    sensors, blade_id)
combo = struct( ...
    'label', {'foundation_center_angle_arrival', ...
              'sg_rising_angle_arrival', ...
              'foundation_center_sg_angle', ...
              'sg_rising_foundation_angle', ...
              'foundation_center_angle_sg_arrival'}, ...
    'opr_source', {'foundation_opr_center', 'sg_opr_rising', ...
                   'foundation_opr_center', 'sg_opr_rising', ...
                   'foundation_opr_center'}, ...
    'angle_source', {'foundation_angle', 'sg_angle', ...
                     'sg_angle', 'foundation_angle', ...
                     'foundation_angle'}, ...
    'arrival_source', {'foundation_probe_half_area', 'sg_probe_fit_arrival', ...
                       'foundation_probe_half_area', ...
                       'foundation_probe_half_area', ...
                       'sg_probe_fit_arrival'});

rows = repmat(make_point_cloud_row_local(), numel(sensors) * numel(combo), 1);
row_idx = 0;
for i = 1:numel(sensors)
    sid = sensors(i);
    pc = get_point_cloud_entry_local(point_cloud, sid, blade_id);
    if isempty(pc.x_mm)
        continue;
    end
    threshold = get_sensor_threshold_local(cfg, sid);
    template_xc = NaN;
    if ~isempty(Template)
        tpl = get_template_entry_local(Template, sid, blade_id);
        if ~isempty(tpl) && isfield(tpl, 'xc')
            template_xc = tpl.xc;
        end
    end

    for j = 1:numel(combo)
        row_idx = row_idx + 1;
        c = combo(j);
        [opr_times, angle_deg, arrivals] = resolve_combo_inputs_local( ...
            c, FCase, SGCase, FConfig, SGConfig, sid, blade_id);
        x_abs = remap_point_cloud_times_local( ...
            pc, opr_times, arrivals, angle_deg, cfg);
        center_mm = weighted_high_voltage_center_local(x_abs, pc.v, threshold);

        rows(row_idx).sensor_id = sid;
        rows(row_idx).blade_id = blade_id;
        rows(row_idx).combo = string(c.label);
        rows(row_idx).opr_source = string(c.opr_source);
        rows(row_idx).angle_source = string(c.angle_source);
        rows(row_idx).arrival_source = string(c.arrival_source);
        rows(row_idx).center_mm = center_mm;
        rows(row_idx).template_xc_mm = template_xc;
        rows(row_idx).delta_from_template_xc_mm = center_mm - template_xc;
        rows(row_idx).point_count = nnz(isfinite(x_abs) & isfinite(pc.v));
        rows(row_idx).pulse_count = numel(unique(pc.pulse_index(isfinite(x_abs))));
    end
end
rows = rows(1:row_idx);
T = struct2table(rows);

base = T(strcmp(T.combo, 'foundation_center_angle_arrival'), :);
if ~isempty(base)
    T.delta_from_foundation_flow_mm = nan(height(T), 1);
    for i = 1:height(base)
        mask = T.sensor_id == base.sensor_id(i);
        T.delta_from_foundation_flow_mm(mask) = ...
            T.center_mm(mask) - base.center_mm(i);
    end
end
end


function row = make_point_cloud_row_local()
row = struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'combo', "", ...
    'opr_source', "", ...
    'angle_source', "", ...
    'arrival_source', "", ...
    'center_mm', NaN, ...
    'template_xc_mm', NaN, ...
    'delta_from_template_xc_mm', NaN, ...
    'point_count', NaN, ...
    'pulse_count', NaN);
end


function T = compare_arrival_x_local( ...
    cfg, FCase, SGCase, FConfig, SGConfig, sensors, blade_id)
combo = struct( ...
    'label', {'foundation_center_angle_arrival', ...
              'sg_rising_angle_arrival', ...
              'foundation_center_sg_angle', ...
              'sg_rising_foundation_angle'}, ...
    'opr_source', {'foundation_opr_center', 'sg_opr_rising', ...
                   'foundation_opr_center', 'sg_opr_rising'}, ...
    'angle_source', {'foundation_angle', 'sg_angle', ...
                     'sg_angle', 'foundation_angle'}, ...
    'arrival_source', {'foundation_probe_half_area', 'sg_probe_fit_arrival', ...
                       'foundation_probe_half_area', ...
                       'foundation_probe_half_area'});
rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'blade_id', NaN, ...
    'combo', "", ...
    'median_arrival_x_mm', NaN, ...
    'mean_arrival_x_mm', NaN, ...
    'std_arrival_x_mm', NaN, ...
    'pulse_count', NaN), numel(sensors) * numel(combo), 1);
row_idx = 0;
for i = 1:numel(sensors)
    sid = sensors(i);
    for j = 1:numel(combo)
        row_idx = row_idx + 1;
        c = combo(j);
        [opr_times, angle_deg, arrivals] = resolve_combo_inputs_local( ...
            c, FCase, SGCase, FConfig, SGConfig, sid, blade_id);
        if strcmp(c.arrival_source, 'sg_probe_fit_arrival')
            start_idx = SGConfig.Target_Indices(sid);
        else
            start_idx = FConfig.Target_Indices(sid);
        end
        x = map_blade_arrivals_to_x_local( ...
            arrivals, opr_times, angle_deg, start_idx, cfg, blade_id);
        rows(row_idx).sensor_id = sid;
        rows(row_idx).blade_id = blade_id;
        rows(row_idx).combo = string(c.label);
        rows(row_idx).median_arrival_x_mm = median(x, 'omitnan');
        rows(row_idx).mean_arrival_x_mm = mean(x, 'omitnan');
        rows(row_idx).std_arrival_x_mm = std(x, 'omitnan');
        rows(row_idx).pulse_count = nnz(isfinite(x));
    end
end
T = struct2table(rows);
end


function T = compare_dynamic_opr_products_local(cfg, sg_root)
case_name = cfg.dynamic_cases{1};
f_file = fullfile(cfg.step02_output_dir, case_name, 'jiluOPR.mat');
sg_file = fullfile(sg_root, 'output', case_name, 'jiluOPR.mat');
if ~isfile(f_file) || ~isfile(sg_file)
    T = table();
    return;
end
F = load(f_file, 'jiluOPR');
S = load(sg_file, 'jiluOPR');
f = F.jiluOPR;
sg = S.jiluOPR;
n = min(size(f, 1), size(sg, 1));
f = f(1:n, :);
sg = sg(1:n, :);
if size(f, 2) < 3 || size(sg, 2) < 2
    T = table();
    return;
end

dt_start = f(:, 2) - sg(:, 1);
dt_end = f(:, 3) - sg(:, 2);
dt_center_vs_rising = f(:, 1) - sg(:, 1);
dt_center_vs_mid = f(:, 1) - 0.5 .* (sg(:, 1) + sg(:, 2));
width_f = f(:, 3) - f(:, 2);
width_sg = sg(:, 2) - sg(:, 1);

epr = cfg.blades_num;
valid_n = n - epr;
speed_deg_s = nan(n, 1);
if valid_n > 0
    speed_deg_s(1:valid_n) = 360 ./ max(sg((epr + 1):n, 1) - sg(1:valid_n, 1), eps);
end
center_shift_deg = dt_center_vs_rising .* speed_deg_s;
center_shift_mm = center_shift_deg * pi / 180 * cfg.r_tip_mm;

row = struct();
row.case_name = string(case_name);
row.row_count_foundation = size(F.jiluOPR, 1);
row.row_count_sg = size(S.jiluOPR, 1);
row.compare_count = n;
row.foundation_start_minus_sg_start_median_us = median(dt_start, 'omitnan') * 1e6;
row.foundation_end_minus_sg_end_median_us = median(dt_end, 'omitnan') * 1e6;
row.foundation_center_minus_sg_rising_median_us = median(dt_center_vs_rising, 'omitnan') * 1e6;
row.foundation_center_minus_sg_mid_median_us = median(dt_center_vs_mid, 'omitnan') * 1e6;
row.foundation_opr_width_median_us = median(width_f, 'omitnan') * 1e6;
row.sg_opr_width_median_us = median(width_sg, 'omitnan') * 1e6;
row.center_vs_rising_shift_deg_median = median(center_shift_deg, 'omitnan');
row.center_vs_rising_shift_mm_median = median(center_shift_mm, 'omitnan');
T = struct2table(row);
end


function T = build_dynamic_opr_event_offset_table_local(cfg, sg_root)
case_name = cfg.dynamic_cases{1};
f_file = fullfile(cfg.step02_output_dir, case_name, 'jiluOPR.mat');
sg_file = fullfile(sg_root, 'output', case_name, 'jiluOPR.mat');
if ~isfile(f_file) || ~isfile(sg_file)
    T = table();
    return;
end
F = load(f_file, 'jiluOPR');
S = load(sg_file, 'jiluOPR');
f = F.jiluOPR;
sg = S.jiluOPR;
n = min(size(f, 1), size(sg, 1));
if size(f, 2) < 3 || size(sg, 2) < 2 || n <= cfg.blades_num
    T = table();
    return;
end
f = f(1:n, :);
sg = sg(1:n, :);

epr = cfg.blades_num;
rev_period = sg((epr + 1):n, 1) - sg(1:(n - epr), 1);
speed_deg_s = 360 ./ max(rev_period, eps);
dt = f(1:(n - epr), 1) - sg(1:(n - epr), 1);
shift_deg = dt .* speed_deg_s;

rows = repmat(struct( ...
    'case_name', "", ...
    'opr_event_in_rev', NaN, ...
    'center_minus_rising_time_us_median', NaN, ...
    'center_minus_rising_time_us_iqr', NaN, ...
    'center_minus_rising_deg_median', NaN, ...
    'center_minus_rising_deg_iqr', NaN, ...
    'center_minus_rising_mm_median', NaN, ...
    'sample_count', NaN), epr, 1);

for event_id = 1:epr
    idx = event_id:epr:numel(shift_deg);
    dt_vals = dt(idx);
    deg_vals = shift_deg(idx);
    finite = isfinite(dt_vals) & isfinite(deg_vals);
    dt_vals = dt_vals(finite);
    deg_vals = deg_vals(finite);

    rows(event_id).case_name = string(case_name);
    rows(event_id).opr_event_in_rev = event_id;
    rows(event_id).sample_count = numel(deg_vals);
    if isempty(deg_vals)
        continue;
    end
    rows(event_id).center_minus_rising_time_us_median = median(dt_vals) * 1e6;
    rows(event_id).center_minus_rising_time_us_iqr = iqr(dt_vals) * 1e6;
    rows(event_id).center_minus_rising_deg_median = median(deg_vals);
    rows(event_id).center_minus_rising_deg_iqr = iqr(deg_vals);
    rows(event_id).center_minus_rising_mm_median = ...
        rows(event_id).center_minus_rising_deg_median * pi / 180 * cfg.r_tip_mm;
end
T = struct2table(rows);
end


function [opr_times, angle_deg, arrivals] = resolve_combo_inputs_local( ...
    combo, FCase, SGCase, FConfig, SGConfig, sid, blade_id)
switch combo.opr_source
    case 'foundation_opr_center'
        opr_times = FCase.opr_times(:);
    case 'sg_opr_rising'
        opr_times = SGCase.opr_times(:);
    otherwise
        error('Unsupported OPR source: %s', combo.opr_source);
end

switch combo.angle_source
    case 'foundation_angle'
        angle_deg = FConfig.Standard_Relative_Angles_OPRCenter(sid, blade_id);
    case 'sg_angle'
        angle_deg = SGConfig.Standard_Relative_Angles(sid, blade_id);
    otherwise
        error('Unsupported angle source: %s', combo.angle_source);
end

switch combo.arrival_source
    case 'foundation_probe_half_area'
        arrivals = FCase.channels(sid).arrival_times(:);
    case 'sg_probe_fit_arrival'
        arrivals = SGCase.channels(sid).arrival_times(:);
    otherwise
        error('Unsupported arrival source: %s', combo.arrival_source);
end
end


function event_id = dominant_previous_opr_event_local(FCase, FConfig, sid, blade_id, epr)
event_id = NaN;
if sid > numel(FCase.channels) || isempty(FCase.channels(sid).arrival_times)
    return;
end
start_idx = FConfig.Target_Indices(sid);
pulses = start_idx + (blade_id - 1):epr:numel(FCase.channels(sid).arrival_times);
if isempty(pulses)
    return;
end
prev = nan(size(pulses));
opr_times = FCase.opr_times(:);
for i = 1:numel(pulses)
    t = FCase.channels(sid).arrival_times(pulses(i));
    idx = find(opr_times < t, 1, 'last');
    if ~isempty(idx)
        prev(i) = idx;
    end
end
prev = prev(isfinite(prev));
if isempty(prev)
    return;
end
event_id = mode(mod(prev - 1, epr) + 1);
end


function x_abs = remap_point_cloud_times_local(pc, opr_times, arrivals, angle_deg, cfg)
t = pc.t_s(:);
pulse_id = pc.pulse_index(:);
x_abs = nan(size(t));
if isempty(t)
    return;
end
F = build_speed_interpolant_local(opr_times, cfg.blades_num);
ids = unique(pulse_id(isfinite(pulse_id))).';
for pulse = ids
    if pulse < 1 || pulse > numel(arrivals) || ~isfinite(arrivals(pulse))
        continue;
    end
    idx_ref = find(opr_times < arrivals(pulse), 1, 'last');
    if isempty(idx_ref)
        continue;
    end
    mask = pulse_id == pulse;
    x_abs(mask) = map_segment_times_to_x_local( ...
        opr_times(idx_ref), t(mask), F, angle_deg, cfg.r_tip_mm);
end
end


function x = map_blade_arrivals_to_x_local( ...
    arrivals, opr_times, angle_deg, start_idx, cfg, blade_id)
pulses = start_idx + (blade_id - 1):cfg.blades_num:numel(arrivals);
x = nan(numel(pulses), 1);
F = build_speed_interpolant_local(opr_times, cfg.blades_num);
for i = 1:numel(pulses)
    pulse = pulses(i);
    t = arrivals(pulse);
    idx_ref = find(opr_times < t, 1, 'last');
    if isempty(idx_ref)
        continue;
    end
    theta_deg = integrate_angle_local(opr_times(idx_ref), t, F);
    theta_diff = wrap_to_180_local(theta_deg - angle_deg);
    x(i) = theta_diff * pi / 180 * cfg.r_tip_mm;
end
end


function x = map_segment_times_to_x_local(t_ref, t, F, angle_deg, r_tip_mm)
t = t(:);
x = nan(size(t));
valid = isfinite(t);
if ~any(valid)
    return;
end
idx = find(valid);
[t_sort, order] = sort(t(valid));
theta = nan(size(t_sort));
theta0 = integrate_angle_local(t_ref, t_sort(1), F);
theta(:) = theta0 + cumtrapz(t_sort, F(t_sort));
theta = theta(order_inverse_local(order));
x(idx) = wrap_to_180_local(theta - angle_deg) .* pi / 180 .* r_tip_mm;
end


function inv = order_inverse_local(order)
inv = zeros(size(order));
inv(order) = 1:numel(order);
end


function F = build_speed_interpolant_local(opr_times, epr)
opr_times = opr_times(:);
if numel(opr_times) <= epr
    error('Not enough OPR events to build speed interpolant.');
end
speed_time = opr_times(1:(end - epr));
speed_deg_s = 360 ./ max(opr_times((epr + 1):end) - opr_times(1:(end - epr)), eps);
valid = isfinite(speed_time) & isfinite(speed_deg_s) & speed_deg_s > 0;
F = griddedInterpolant(speed_time(valid), speed_deg_s(valid), 'linear', 'nearest');
end


function theta_deg = integrate_angle_local(t0, t1, F)
t_grid = linspace(t0, t1, 10);
theta_deg = trapz(t_grid, F(t_grid));
end


function center = weighted_high_voltage_center_local(x, v, threshold)
x = x(:);
v = v(:);
finite = isfinite(x) & isfinite(v);
if nnz(finite) < 10
    center = NaN;
    return;
end
vf = v(finite);
xf = x(finite);
baseline = prctile(vf, 5);
weight = max(vf - baseline, 0);
high = vf >= threshold & weight > 0;
if nnz(high) < 8
    high = weight > 0;
end
if nnz(high) < 8 || sum(weight(high), 'omitnan') <= 0
    center = NaN;
else
    center = sum(xf(high) .* weight(high), 'omitnan') ./ ...
        sum(weight(high), 'omitnan');
end
end


function pc = get_point_cloud_entry_local(point_cloud, sid, blade_id)
idx = find([point_cloud.sensor_id] == sid & [point_cloud.blade_id] == blade_id, 1);
if isempty(idx)
    error('Missing point_cloud entry CH%d B%d.', sid, blade_id);
end
pc = point_cloud(idx);
end


function tpl = get_template_entry_local(Template, sid, blade_id)
tpl = [];
if ~isfield(Template, 'SensorBlade')
    return;
end
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id
        tpl = Template.SensorBlade(i);
        return;
    end
end
end


function threshold = get_sensor_threshold_local(cfg, sid)
threshold = cfg.sensor_threshold_default;
if isfield(cfg, 'sensor_thresholds') && isa(cfg.sensor_thresholds, 'containers.Map')
    if isKey(cfg.sensor_thresholds, sid)
        threshold = cfg.sensor_thresholds(sid);
    end
end
end


function y = wrap_to_180_local(x)
y = mod(x + 180, 360) - 180;
end


function print_diagnostic_findings_local(AngleTable, PointCloudTable, DynamicOPRTable)
fprintf('\nKey diagnostic checks:\n');
if ~isempty(AngleTable)
    r = AngleTable(1, :);
    fprintf(['  1) Foundation-SG standard-angle delta is a common %.6f deg, ' ...
        'but OPR center-rising offset is event-specific.\n'], ...
        r.foundation_minus_sg_angle_deg);
    for i = 1:height(AngleTable)
        fprintf(['     CH%d B%d uses previous OPR event %d: offset %.6f deg, ' ...
            'residual %.4f mm.\n'], ...
            AngleTable.sensor_id(i), AngleTable.blade_id(i), ...
            AngleTable.dominant_previous_opr_event(i), ...
            AngleTable.event_center_minus_rising_deg_median(i), ...
            AngleTable.event_specific_compensation_residual_mm(i));
    end
end

if ~isempty(PointCloudTable)
    base = PointCloudTable(strcmp(PointCloudTable.combo, ...
        'foundation_center_angle_arrival'), :);
    sg = PointCloudTable(strcmp(PointCloudTable.combo, ...
        'sg_rising_angle_arrival'), :);
    if ~isempty(base) && ~isempty(sg)
        fprintf('  2) Point-cloud high-voltage centers, foundation flow vs SG flow:\n');
        for i = 1:height(base)
            j = find(sg.sensor_id == base.sensor_id(i), 1);
            if isempty(j)
                continue;
            end
            fprintf('     CH%d: foundation %.4f mm, SG-flow %.4f mm, shift %.4f mm.\n', ...
                base.sensor_id(i), base.center_mm(i), sg.center_mm(j), ...
                sg.center_mm(j) - base.center_mm(i));
        end
    end
end

if ~isempty(DynamicOPRTable)
    fprintf(['  3) Dynamic Step02 OPR products share the same start/end as SG, ' ...
        'but foundation Step05 uses center. Median center-rising offset = %.3f us ' ...
        '(%.3f mm at local speed).\n'], ...
        DynamicOPRTable.foundation_center_minus_sg_rising_median_us, ...
        DynamicOPRTable.center_vs_rising_shift_mm_median);
end
end
