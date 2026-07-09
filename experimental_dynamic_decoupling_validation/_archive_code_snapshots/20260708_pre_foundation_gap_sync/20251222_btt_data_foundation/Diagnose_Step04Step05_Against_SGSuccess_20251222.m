function Summary = Diagnose_Step04Step05_Against_SGSuccess_20251222()
%DIAGNOSE_STEP04STEP05_AGAINST_SGSUCCESS_20251222
% Compare the current foundation Step04/Step05 route with the older
% super-Gaussian 20251222 route that identifies EO14 robustly.

cfg = BTTProjectConfig_20251222();
route_dir = fileparts(mfilename('fullpath'));
parent_root = fileparts(fileparts(route_dir));

sg_root = strtrim(getenv('SG_SUCCESS_ROOT'));
if isempty(sg_root)
    sg_root = fullfile(parent_root, ...
        '7超高斯模型-权重-瞬态', '程序', '直叶片验证', '实验验证', ...
        '超高斯', '20251222适配-3号传感间隙变化-改');
end
if ~isfolder(sg_root)
    error(['Cannot find SG success route. Set SG_SUCCESS_ROOT to the folder ' ...
        'that contains SGCalib_20251222_B1_S1.mat. Tried: %s'], sg_root);
end

target_blade = 1;
sensors = [1 2 3];
sensor_tag = ['S', sprintf('%d', sensors)];
case_name = cfg.dynamic_cases{1};

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'foundation_vs_supergaussian_success_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

Step04Compare = compare_step04_calibration_local(cfg, sg_root, ...
    target_blade, sensors, sensor_tag);
MethodSummary = compare_step05_results_local(cfg, sg_root, ...
    case_name, target_blade, sensor_tag);
EtaSummary = summarize_eta_preview_local(cfg, target_blade, sensor_tag);

writetable(Step04Compare, fullfile(out_dir, ...
    'Step04_Foundation_vs_SG_Calibration_B1_S123_20251222.csv'));
writetable(MethodSummary, fullfile(out_dir, ...
    'Step05_Foundation_vs_SG_ResultSummary_B1_S123_20251222.csv'));
if ~isempty(EtaSummary)
    writetable(EtaSummary, fullfile(out_dir, ...
        'Step05_StaticEtaPreview_B1_S123_20251222.csv'));
end

Summary = struct();
Summary.SGRoot = sg_root;
Summary.OutputDir = out_dir;
Summary.Step04Compare = Step04Compare;
Summary.MethodSummary = MethodSummary;
Summary.EtaSummary = EtaSummary;

fprintf('\n=== 20251222 foundation vs SG-success diagnostic ===\n');
fprintf('SG success root:\n  %s\n', sg_root);
fprintf('Diagnostic CSV output:\n  %s\n\n', out_dir);
disp(Step04Compare);
disp(MethodSummary);
if ~isempty(EtaSummary)
    disp(EtaSummary);
end

report_findings_local(Step04Compare, MethodSummary, EtaSummary);
end


function T = compare_step04_calibration_local(cfg, sg_root, blade_id, sensors, sensor_tag)
template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', ...
    sensor_tag, cfg.dataset));
if ~isfile(template_file)
    error('Missing foundation Step04 template: %s', template_file);
end
loaded = load(template_file, 'Template');
Template = loaded.Template;

sg_plan = table();
sg_plan_file = fullfile(sg_root, 'output', 'step4_window_scan', ...
    'Step4_StableWindowPlan_20251222.csv');
if isfile(sg_plan_file)
    sg_plan = readtable(sg_plan_file);
end

rows = repmat(make_empty_step04_row_local(), numel(sensors), 1);
for i = 1:numel(sensors)
    sid = sensors(i);
    tpl = get_template_entry_local(Template, sid, blade_id);
    sg = load_sg_calibration_local(sg_root, blade_id, sid);
    plan_row = get_sg_plan_row_local(sg_plan, blade_id, sid);

    rows(i).sensor_id = sid;
    rows(i).foundation_angle_deg = get_standard_angle_local(Template, sid, blade_id);
    rows(i).sg_angle_deg = sg.alpha_ref;
    rows(i).angle_delta_deg = rows(i).foundation_angle_deg - rows(i).sg_angle_deg;
    rows(i).angle_delta_equiv_mm = rows(i).angle_delta_deg * pi / 180 * cfg.r_tip_mm;

    rows(i).foundation_xc_mm = tpl.xc;
    rows(i).foundation_xc_detected_mm = getfield_default_local(tpl, 'xc_detected', NaN);
    rows(i).sg_xc_mm = sg.xc;
    rows(i).xc_delta_foundation_minus_sg_mm = tpl.xc - sg.xc;

    rows(i).foundation_stable_lap_start = getfield_default_local(tpl, ...
        'stable_window_start_lap', NaN);
    rows(i).foundation_stable_lap_end = getfield_default_local(tpl, ...
        'stable_window_end_lap', NaN);
    rows(i).foundation_pulse_count = getfield_default_local(tpl, ...
        'pulse_count', NaN);
    rows(i).foundation_x_domain_left_mm = tpl.x_domain(1);
    rows(i).foundation_x_domain_right_mm = tpl.x_domain(2);
    rows(i).foundation_query_left_mm = tpl.x_query_safe_domain(1);
    rows(i).foundation_query_right_mm = tpl.x_query_safe_domain(2);
    rows(i).foundation_quality = string(getfield_default_local(tpl, ...
        'quality_status', ''));

    rows(i).sg_lap_count = sg.lap_count;
    rows(i).sg_trust_left_mm = sg.dx_trust_left;
    rows(i).sg_trust_right_mm = sg.dx_trust_right;
    rows(i).sg_rmse_v = sg.rmse;
    rows(i).sg_window_start_time_s = getfield_default_local(plan_row, ...
        'StartTime', NaN);
    rows(i).sg_window_end_time_s = getfield_default_local(plan_row, ...
        'EndTime', NaN);
    rows(i).sg_start_offset_lap = getfield_default_local(plan_row, ...
        'StartOffsetLap', NaN);
end
T = struct2table(rows);
end


function T = compare_step05_results_local(cfg, sg_root, case_name, blade_id, sensor_tag)
rows = {};

base_dir = fullfile(cfg.step05_direct_template_output_dir, case_name);
rows{end + 1, 1} = summarize_trend_file_local( ...
    'foundation_noeta', ...
    fullfile(base_dir, sprintf( ...
    'Trend_Step05_FoundationMainPulseAdaptiveNoEta_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset)));
rows{end + 1, 1} = summarize_trend_file_local( ...
    'foundation_fixed_dynamic_residual_eta', ...
    fullfile(base_dir, sprintf( ...
    'Trend_Step05_FoundationMainPulseAdaptiveFixedJointEta_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset)));
rows{end + 1, 1} = summarize_trend_file_local( ...
    'foundation_old_vp_improved_free_eta', ...
    fullfile(base_dir, sprintf( ...
    'Trend_Step05_FoundationFullRawVPScreenTop3RMSE_B%d_%s_%s.csv', ...
    blade_id, sensor_tag, cfg.dataset)));

sg_trend_50 = fullfile(sg_root, sprintf( ...
    'Result_single_sync_20251222_B%d_%s_%s_Start50p0s_trend.csv', ...
    blade_id, sensor_tag, case_name));
rows{end + 1, 1} = summarize_trend_file_local( ...
    'sg_success_start50p0', sg_trend_50);

sg_trend_51 = fullfile(sg_root, sprintf( ...
    'Result_single_sync_20251222_B%d_%s_%s_Start51p0s_trend.csv', ...
    blade_id, sensor_tag, case_name));
rows{end + 1, 1} = summarize_trend_file_local( ...
    'sg_success_start51p0', sg_trend_51);

T = struct2table(vertcat(rows{:}));
end


function T = summarize_eta_preview_local(cfg, blade_id, sensor_tag)
eta_file = fullfile(cfg.output_root, 'step05_joint_static_eta_preview', ...
    sprintf('JointStaticEtaPreview_B%d_%s.mat', blade_id, sensor_tag));
if ~isfile(eta_file)
    T = table();
    return;
end
loaded = load(eta_file, 'Summary');
S = loaded.Summary;
C = S.Consensus;

eta = C.eta_median_mm(:).';
iqr = nan(size(eta));
if isfield(C, 'eta_plan_iqr_mm')
    iqr = C.eta_plan_iqr_mm(:).';
end
diag_status = "";
if isfield(C, 'diagnostic_status')
    diag_status = string(C.diagnostic_status);
elseif isfield(C, 'status')
    diag_status = string(C.status);
end
policy = "";
if isfield(C, 'eta_source_policy')
    policy = string(C.eta_source_policy);
end

rows = repmat(struct( ...
    'sensor_index', NaN, ...
    'eta_median_mm', NaN, ...
    'eta_plan_iqr_mm', NaN, ...
    'diagnostic_status', "", ...
    'eta_source_policy', "", ...
    'eta_file', ""), numel(eta), 1);
for i = 1:numel(eta)
    rows(i).sensor_index = i;
    rows(i).eta_median_mm = eta(i);
    rows(i).eta_plan_iqr_mm = iqr(i);
    rows(i).diagnostic_status = diag_status;
    rows(i).eta_source_policy = policy;
    rows(i).eta_file = string(eta_file);
end
T = struct2table(rows);
end


function report_findings_local(Step04Compare, MethodSummary, EtaSummary)
fprintf('\nKey diagnostic flags:\n');

angle_delta = Step04Compare.angle_delta_deg;
if all(isfinite(angle_delta)) && std(angle_delta, 'omitnan') < 1e-5
    fprintf(['  1) Foundation and SG Sensor_Config differ by a common ' ...
        'angle offset of %.6f deg (%.4f mm at the tip). This mostly ' ...
        'cancels after centering, but confirms the routes are not using ' ...
        'identical Step1 products.\n'], ...
        mean(angle_delta, 'omitnan'), ...
        mean(Step04Compare.angle_delta_equiv_mm, 'omitnan'));
else
    fprintf('  1) Foundation and SG sensor angle differences are not a pure common offset.\n');
end

max_xc_delta = max(abs(Step04Compare.xc_delta_foundation_minus_sg_mm), [], 'omitnan');
fprintf('  2) Max |foundation xc - SG xc| = %.4f mm.\n', max_xc_delta);
bad_xc = Step04Compare(abs(Step04Compare.xc_delta_foundation_minus_sg_mm) > 0.10, :);
if ~isempty(bad_xc)
    fprintf('     Sensors with >0.10 mm center mismatch: %s\n', ...
        strjoin("CH" + string(bad_xc.sensor_id(:).'), ', '));
end

if ~isempty(EtaSummary)
    eta_iqr = EtaSummary.eta_plan_iqr_mm;
    eta_iqr = eta_iqr(isfinite(eta_iqr));
    bad_eta = any(abs(eta_iqr) > 0.05) || ...
        any(strcmpi(EtaSummary.diagnostic_status, "unstable_plan_spread"));
    if bad_eta
        fprintf(['  3) Fixed dynamic-residual eta is unstable. Formal ' ...
            'fixed-eta Step05 should not report a final EO for this run.\n']);
    end
end

sg_ok = MethodSummary(strcmpi(MethodSummary.method, "sg_success_start50p0"), :);
if ~isempty(sg_ok)
    fprintf('  4) SG success route Start50p0 EO sequence: %s.\n', ...
        char(sg_ok.eo_sequence(1)));
end
foundation = MethodSummary(startsWith(MethodSummary.method, "foundation"), :);
if ~isempty(foundation)
    fprintf('     Foundation EO sequences:\n');
    for i = 1:height(foundation)
        fprintf('       %s: %s\n', char(string(foundation.method(i))), ...
            char(string(foundation.eo_sequence(i))));
    end
end
end


function row = summarize_trend_file_local(method, file)
row = struct( ...
    'method', string(method), ...
    'file_exists', isfile(file), ...
    'window_count', NaN, ...
    'dominant_eo', NaN, ...
    'unique_eo_count', NaN, ...
    'eo_sequence', "", ...
    'median_rmse_v', NaN, ...
    'mean_rmse_v', NaN, ...
    'median_A_mm', NaN, ...
    'time_start_s', NaN, ...
    'time_end_s', NaN, ...
    'file', string(file));
if ~isfile(file)
    return;
end

T = readtable(file);
if isempty(T) || ~ismember('EO_id', T.Properties.VariableNames)
    return;
end

eo = T.EO_id;
ok = isfinite(eo);
if ismember('status', T.Properties.VariableNames)
    ok = ok & strcmpi(string(T.status), "ok");
end
row.window_count = nnz(ok);
if any(ok)
    row.dominant_eo = mode(eo(ok));
    row.unique_eo_count = numel(unique(eo(ok)));
    row.eo_sequence = strjoin(compose('%d', eo(ok).'), '|');
end
if ismember('weighted_voltage_rmse', T.Properties.VariableNames)
    row.median_rmse_v = median(T.weighted_voltage_rmse(ok), 'omitnan');
    row.mean_rmse_v = mean(T.weighted_voltage_rmse(ok), 'omitnan');
end
if ismember('A_id', T.Properties.VariableNames)
    row.median_A_mm = median(T.A_id(ok), 'omitnan');
end
if ismember('window_center_time_s', T.Properties.VariableNames)
    t = T.window_center_time_s(ok);
elseif ismember('window_center_time', T.Properties.VariableNames)
    t = T.window_center_time(ok);
else
    t = [];
end
if ~isempty(t)
    row.time_start_s = min(t, [], 'omitnan');
    row.time_end_s = max(t, [], 'omitnan');
end
end


function row = make_empty_step04_row_local()
row = struct( ...
    'sensor_id', NaN, ...
    'foundation_angle_deg', NaN, ...
    'sg_angle_deg', NaN, ...
    'angle_delta_deg', NaN, ...
    'angle_delta_equiv_mm', NaN, ...
    'foundation_xc_mm', NaN, ...
    'foundation_xc_detected_mm', NaN, ...
    'sg_xc_mm', NaN, ...
    'xc_delta_foundation_minus_sg_mm', NaN, ...
    'foundation_stable_lap_start', NaN, ...
    'foundation_stable_lap_end', NaN, ...
    'foundation_pulse_count', NaN, ...
    'foundation_x_domain_left_mm', NaN, ...
    'foundation_x_domain_right_mm', NaN, ...
    'foundation_query_left_mm', NaN, ...
    'foundation_query_right_mm', NaN, ...
    'foundation_quality', "", ...
    'sg_lap_count', NaN, ...
    'sg_trust_left_mm', NaN, ...
    'sg_trust_right_mm', NaN, ...
    'sg_rmse_v', NaN, ...
    'sg_window_start_time_s', NaN, ...
    'sg_window_end_time_s', NaN, ...
    'sg_start_offset_lap', NaN);
end


function sg = load_sg_calibration_local(sg_root, blade_id, sid)
file = fullfile(sg_root, sprintf('SGCalib_20251222_B%d_S%d.mat', blade_id, sid));
if ~isfile(file)
    error('Missing SG calibration file: %s', file);
end
loaded = load(file, 'Calib', 'summary_table');
sensor = loaded.Calib.Sensor;
sg = struct();
sg.alpha_ref = sensor.alpha_ref;
sg.xc = sensor.xc;
sg.dx_trust_left = sensor.dx_trust_left;
sg.dx_trust_right = sensor.dx_trust_right;
sg.rmse = NaN;
sg.lap_count = NaN;
if isfield(loaded, 'summary_table') && ~isempty(loaded.summary_table)
    sg.rmse = loaded.summary_table.RMSE(1);
    sg.lap_count = loaded.summary_table.LapCount(1);
end
end


function tpl = get_template_entry_local(Template, sid, blade_id)
tpl = [];
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id
        tpl = Template.SensorBlade(i);
        return;
    end
end
error('Template has no SensorBlade entry for CH%d B%d.', sid, blade_id);
end


function angle = get_standard_angle_local(Template, sid, blade_id)
angles = Template.Standard_Relative_Angles;
angle = angles(sid, blade_id);
end


function row = get_sg_plan_row_local(plan, blade_id, sid)
row = struct();
if isempty(plan)
    return;
end
idx = find(plan.BladeID == blade_id & plan.SensorID == sid, 1, 'first');
if isempty(idx)
    return;
end
names = plan.Properties.VariableNames;
for i = 1:numel(names)
    row.(names{i}) = plan.(names{i})(idx);
end
end


function value = getfield_default_local(s, field, default_value)
if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
    value = s.(field);
else
    value = default_value;
end
end
