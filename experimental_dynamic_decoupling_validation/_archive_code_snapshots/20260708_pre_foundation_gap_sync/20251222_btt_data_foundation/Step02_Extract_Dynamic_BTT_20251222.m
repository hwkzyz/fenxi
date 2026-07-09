clc; clear; close all;

%STEP02_EXTRACT_DYNAMIC_BTT_20251222 Extract dynamic BTT timings from raw data.
%
% This step is the dynamic data layer of the 20251222 foundation route.
% It reads raw channel MAT files directly from cfg.dataset_root/case_name,
% computes the configured OPR reference times, recomputes speed, extracts probe arrivals,
% and assigns blade IDs from the Step01 low-speed fingerprints.
%
% It does not read old dynamic outputs from newflow/legacy folders.

%% 1. Run settings
cfg = BTTDataConfig_20251222();
target_cases = cfg_cellstr_local(cfg.dynamic_cases);
show_tables = true;
show_figures = true;
force_rebuild = true;

%% 2. Algorithm settings kept local to Step02
Step02_Setting = struct();
Step02_Setting.source_mode = 'raw_dynamic_extraction';
Step02_Setting.allow_old_dynamic_outputs_as_source = false;
Step02_Setting.raw_file_pattern = '4-<channel>-<file_id>.mat';
Step02_Setting.initial_gap_points = 1e4;
Step02_Setting.opr_center_method = 'multi_threshold_width_center';
Step02_Setting.opr_center_level_ratios = [0.30 0.40 0.50 0.60 0.70];
Step02_Setting.opr_timing_reference = get_field_or_default_local( ...
    cfg, 'opr_timing_reference', 'rising_edge');
Step02_Setting.opr_reference_time_column = get_field_or_default_local( ...
    cfg, 'opr_reference_time_column', 'opr_start_time_s');
Step02_Setting.opr_events_per_revolution_candidates = unique([1 cfg.blades_num], 'stable');
Step02_Setting.expected_rpm_range = [500 4500];
Step02_Setting.fingerprint_min_corr = 0.85;
Step02_Setting.fingerprint_min_corr_gap = 0.02;
Step02_Setting.fill_missing_blade_ids = false;
Step02_Setting.max_plot_points = 6000;
Step02_Setting.max_plot_revolutions = 250;
Step02_Setting.save_figures = cfg.save_figures;

%% 3. Load Step01 reference from the new foundation folder
sensor_config_file = fullfile(cfg.step01_output_dir, 'Sensor_Config_20251222.mat');
if ~isfile(sensor_config_file)
    error('Missing Step01 output. Run Step01_Build_LowSpeed_Reference_20251222 first:\n  %s', sensor_config_file);
end
loaded_step01 = load(sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_step01.Sensor_Config;
validate_step01_reference_local(Sensor_Config, cfg, sensor_config_file);

%% 4. Main loop
for iCase = 1:numel(target_cases)
    case_name = char(target_cases{iCase});
    case_dir = fullfile(cfg.dataset_root, case_name);
    out_dir = fullfile(cfg.step02_output_dir, case_name);
    fig_dir = fullfile(cfg.step02_figure_dir, case_name);
    ensure_dir_local(out_dir);
    if Step02_Setting.save_figures
        ensure_dir_local(fig_dir);
    end

    feature_out = fullfile(out_dir, 'DynamicBTTFeature_20251222.mat');
    if ~force_rebuild && isfile(feature_out)
        fprintf('Step02 output already exists:\n  %s\n', feature_out);
        continue;
    end

    fprintf('\n=== Step02 raw dynamic extraction [%s] ===\n', case_name);
    fprintf('Raw case dir:\n  %s\n', case_dir);
    if ~isfolder(case_dir)
        error('Raw dynamic case directory not found:\n  %s', case_dir);
    end

    SourceInventoryTable = inspect_raw_case_sources_local(case_dir, cfg, case_name);
    file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
    if isempty(file_ids)
        error('No OPR raw files found in:\n  %s', case_dir);
    end

    OPRRaw = extract_opr_from_raw_local(case_dir, file_ids, cfg, Step02_Setting);
    [OPRSemanticsTable, opr_events_per_revolution] = infer_opr_events_per_revolution_local( ...
        OPRRaw.reference_time_s, cfg, Step02_Setting);
    [OPRTable, Omega, OmegaTable] = build_opr_and_omega_tables_local( ...
        OPRRaw, opr_events_per_revolution, cfg, case_name);
    Step01Step02Consistency = check_step01_step02_consistency_local( ...
        Sensor_Config, cfg, case_name, SourceInventoryTable, opr_events_per_revolution);

    all_sensor_tables = cell(numel(cfg.sensor_ids), 1);
    all_sensor_quality = cell(numel(cfg.sensor_ids), 1);
    all_rev_quality = cell(numel(cfg.sensor_ids), 1);
    ProbeDiagnostics = struct();

    for iSensor = 1:numel(cfg.sensor_ids)
        sid = cfg.sensor_ids(iSensor);
        fprintf('>>> [%s] CH%d extracting raw probe pulses...\n', case_name, sid);
        ProbeRaw = extract_probe_from_raw_local(case_dir, file_ids, sid, cfg, Step02_Setting);
        [PulseTable, SensorQuality, RevQuality, ProbeDiag] = assign_blades_from_step01_local( ...
            ProbeRaw, OPRTable, Sensor_Config, cfg, Step02_Setting, case_name, sid);
        all_sensor_tables{iSensor} = PulseTable;
        all_sensor_quality{iSensor} = SensorQuality;
        all_rev_quality{iSensor} = RevQuality;
        ProbeDiagnostics.(sprintf('CH%d', sid)) = ProbeDiag;
    end

    DynamicPulseTable = vertcat(all_sensor_tables{:});
    SensorNumberingQuality = vertcat(all_sensor_quality{:});
    RevolutionNumberingQuality = vertcat(all_rev_quality{:});
    Step02Summary = build_step02_summary_local(cfg, case_name, DynamicPulseTable, ...
        SensorNumberingQuality, RevolutionNumberingQuality, OPRTable, OmegaTable);
    Step02Meta = build_step02_metadata_local(cfg, Step02_Setting, case_name, case_dir, ...
        sensor_config_file, Sensor_Config, OPRRaw, opr_events_per_revolution, ProbeDiagnostics);

    DynamicBTTFeature = struct();
    DynamicBTTFeature.Dataset = cfg.dataset;
    DynamicBTTFeature.CaseName = case_name;
    DynamicBTTFeature.SourceMode = Step02_Setting.source_mode;
    DynamicBTTFeature.RawCaseDir = case_dir;
    DynamicBTTFeature.OPRTable = OPRTable;
    DynamicBTTFeature.DynamicPulseTable = DynamicPulseTable;
    DynamicBTTFeature.SensorNumberingQuality = SensorNumberingQuality;
    DynamicBTTFeature.RevolutionNumberingQuality = RevolutionNumberingQuality;
    DynamicBTTFeature.Step02Summary = Step02Summary;
    DynamicBTTFeature.OmegaTable = OmegaTable;
    DynamicBTTFeature.OPRSemanticsTable = OPRSemanticsTable;
    DynamicBTTFeature.SourceInventoryTable = SourceInventoryTable;
    DynamicBTTFeature.Step01Step02Consistency = Step01Step02Consistency;
    DynamicBTTFeature.Metadata = Step02Meta;

    save(feature_out, 'DynamicBTTFeature', 'OPRTable', 'DynamicPulseTable', ...
        'SensorNumberingQuality', 'RevolutionNumberingQuality', 'Step02Summary', ...
        'OmegaTable', 'OPRSemanticsTable', 'SourceInventoryTable', ...
        'Step01Step02Consistency', 'Step02Meta', '-v7.3');
    save(fullfile(out_dir, 'Step02_Metadata_20251222.mat'), 'Step02Meta');
    writetable(DynamicPulseTable, fullfile(out_dir, 'DynamicPulseTable_20251222.csv'));
    writetable(SensorNumberingQuality, fullfile(out_dir, 'SensorNumberingQuality_20251222.csv'));
    writetable(RevolutionNumberingQuality, fullfile(out_dir, 'RevolutionNumberingQuality_20251222.csv'));
    writetable(Step02Summary, fullfile(out_dir, 'Step02_Summary_20251222.csv'));
    writetable(OPRTable, fullfile(out_dir, 'OPRTable_20251222.csv'));
    writetable(OmegaTable, fullfile(out_dir, 'OmegaTable_20251222.csv'));
    writetable(OPRSemanticsTable, fullfile(out_dir, 'OPRSemanticsCheck_20251222.csv'));
    writetable(SourceInventoryTable, fullfile(out_dir, 'Step02_SourceInventory_20251222.csv'));
    writetable(Step01Step02Consistency, fullfile(out_dir, 'Step01_Step02_Consistency_20251222.csv'));

    save_compatibility_outputs_local(out_dir, OPRRaw, Omega, DynamicPulseTable, cfg, Step02Meta);

    if show_figures || Step02_Setting.save_figures
        plot_step02_diagnostics_local(cfg, Step02_Setting, case_name, OPRTable, OmegaTable, ...
            DynamicPulseTable, SensorNumberingQuality, RevolutionNumberingQuality, ...
            fig_dir, show_figures, Step02_Setting.save_figures);
    end

    if show_tables
        disp(Step02Summary);
        disp(SensorNumberingQuality);
    end
end

function cases = cfg_cellstr_local(cases)
if ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end
end

function ensure_dir_local(folder)
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end

function validate_step01_reference_local(Sensor_Config, cfg, sensor_config_file)
if ~isfield(Sensor_Config, 'Fingerprints')
    error('Step01 Sensor_Config has no Fingerprints field:\n  %s', sensor_config_file);
end
for sid = cfg.sensor_ids(:).'
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Step01 Sensor_Config.Fingerprints is missing sensor %d:\n  %s', sid, sensor_config_file);
    end
    fp = Sensor_Config.Fingerprints(sid);
    if numel(fp) ~= cfg.blades_num
        error('Step01 fingerprint length for sensor %d is %d, expected cfg.blades_num=%d.', ...
            sid, numel(fp), cfg.blades_num);
    end
end
end

function T = inspect_raw_case_sources_local(case_dir, cfg, case_name)
rows = [];
channels = [cfg.sensor_ids(:).', cfg.opr_id];
for channel_id = channels
    d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
    file_ids = nan(numel(d), 1);
    bytes_total = 0;
    for i = 1:numel(d)
        tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
        if ~isempty(tok)
            file_ids(i) = str2double(tok{1});
        end
        bytes_total = bytes_total + d(i).bytes;
    end
    file_ids = sort(file_ids(isfinite(file_ids)));
    first_file_id = NaN;
    last_file_id = NaN;
    if ~isempty(file_ids)
        first_file_id = file_ids(1);
        last_file_id = file_ids(end);
    end
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'source_mode', "raw_dynamic_extraction", ...
        'channel_id', channel_id, ...
        'file_count', numel(file_ids), ...
        'first_file_id', first_file_id, ...
        'last_file_id', last_file_id, ...
        'bytes_total', bytes_total, ...
        'source_dir', string(case_dir), ...
        'file_pattern', string(sprintf('4-%d-*.mat', channel_id)))];
end
T = struct2table(rows);
end

function file_ids = list_case_file_ids_local(case_dir, channel_id)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    if isempty(tok)
        continue;
    end
    file_ids(i) = str2double(tok{1});
end
file_ids = sort(unique(file_ids(isfinite(file_ids))));
end

function OPRRaw = extract_opr_from_raw_local(case_dir, file_ids, cfg, Step02_Setting)
tail = [];
gap_points = Step02_Setting.initial_gap_points;
rows = [];

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_case_channel_local(case_dir, cfg.opr_id, file_id);
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gap_points] = segment_signal_local(raw, cfg.opr_threshold, gap_points, true);
    rows = [rows; build_opr_rows_from_segments_local(raw, segments, file_id, cfg, Step02_Setting)]; %#ok<AGROW>
end

if ~isempty(tail)
    [segments, ~, ~] = segment_signal_local(tail, cfg.opr_threshold, gap_points, false);
    rows = [rows; build_opr_rows_from_segments_local(tail, segments, file_ids(end), cfg, Step02_Setting)]; %#ok<AGROW>
end

if isempty(rows)
    error('No OPR pulses extracted from raw dynamic data.');
end

[~, order] = sort([rows.center_time_s]);
rows = rows(order);
OPRRaw = struct();
OPRRaw.center_time_s = [rows.center_time_s].';
OPRRaw.start_time_s = [rows.start_time_s].';
OPRRaw.end_time_s = [rows.end_time_s].';
OPRRaw.width_s = [rows.width_s].';
OPRRaw.center_offset_s = [rows.center_offset_s].';
OPRRaw.peak_value = [rows.peak_value].';
OPRRaw.file_id = [rows.file_id].';
OPRRaw.center_quality = [rows.center_quality].';
OPRRaw.center_method = Step02_Setting.opr_center_method;
OPRRaw.reference_time_s = select_opr_reference_times_local(OPRRaw, Step02_Setting);
OPRRaw.reference_method = resolve_opr_reference_method_local(Step02_Setting);
fprintf('>>> OPR raw pulses extracted: %d\n', numel(OPRRaw.center_time_s));
fprintf('>>> OPR timing reference: %s\n', OPRRaw.reference_method);
end

function rows = build_opr_rows_from_segments_local(raw, segments, file_id, cfg, Step02_Setting)
rows = [];
for iSeg = 1:numel(segments.start_idx)
    a = segments.start_idx(iSeg);
    b = segments.end_idx(iSeg);
    [center_t, peak_v, q] = compute_opr_center_local(raw, a, b, cfg, Step02_Setting);
    start_t = raw(a, 1) / cfg.pinlv;
    end_t = raw(b, 1) / cfg.pinlv;
    rows = [rows; struct( ... %#ok<AGROW>
        'center_time_s', center_t, ...
        'start_time_s', start_t, ...
        'end_time_s', end_t, ...
        'width_s', end_t - start_t, ...
        'center_offset_s', center_t - 0.5 * (start_t + end_t), ...
        'peak_value', peak_v, ...
        'file_id', file_id, ...
        'center_quality', q)];
end
end

function [center_t, peak_v, quality] = compute_opr_center_local(raw, a, b, cfg, Step02_Setting)
expand_pts = max(3, floor((b - a + 1) * 0.4));
a2 = max(1, a - expand_pts);
b2 = min(size(raw, 1), b + expand_pts);
t = raw(a2:b2, 1) / cfg.pinlv;
v = raw(a2:b2, 2);
try
    v_smooth = smoothdata(v, 'movmean', min(21, max(3, 2 * floor(numel(v) / 8) + 1)));
catch
    v_smooth = v;
end

peak_v = max(v_smooth);
baseline = min(v_smooth);
amp = peak_v - baseline;
centers = nan(numel(Step02_Setting.opr_center_level_ratios), 1);
if amp > 0 && numel(t) >= 3
    for i = 1:numel(Step02_Setting.opr_center_level_ratios)
        level = baseline + Step02_Setting.opr_center_level_ratios(i) * amp;
        centers(i) = pulse_width_center_at_level_local(t, v_smooth, level);
    end
end
centers = centers(isfinite(centers));
if isempty(centers)
    center_t = 0.5 * (raw(a, 1) + raw(b, 1)) / cfg.pinlv;
    quality = 0;
else
    center_t = median(centers);
    quality = numel(centers) / numel(Step02_Setting.opr_center_level_ratios);
end
end

function reference_time_s = select_opr_reference_times_local(OPRRaw, Step02_Setting)
ref = lower(strtrim(string(get_field_or_default_local( ...
    Step02_Setting, 'opr_timing_reference', 'rising_edge'))));
switch ref
    case {"rising_edge", "start_edge", "start", "threshold_rising_edge"}
        reference_time_s = OPRRaw.start_time_s(:);
    case {"center", "opr_center", "pulse_center", "multi_threshold_width_center"}
        reference_time_s = OPRRaw.center_time_s(:);
    case {"falling_edge", "end_edge", "end"}
        reference_time_s = OPRRaw.end_time_s(:);
    otherwise
        error('Unsupported OPR timing reference: %s.', ref);
end
end

function method = resolve_opr_reference_method_local(Step02_Setting)
ref = lower(strtrim(string(get_field_or_default_local( ...
    Step02_Setting, 'opr_timing_reference', 'rising_edge'))));
switch ref
    case {"rising_edge", "start_edge", "start", "threshold_rising_edge"}
        method = 'threshold_rising_edge';
    case {"center", "opr_center", "pulse_center", "multi_threshold_width_center"}
        method = char(string(Step02_Setting.opr_center_method));
    case {"falling_edge", "end_edge", "end"}
        method = 'threshold_falling_edge';
    otherwise
        method = char(ref);
end
end

function center_t = pulse_width_center_at_level_local(t, v, level)
above = find(v >= level);
if isempty(above)
    center_t = NaN;
    return;
end
i1 = above(1);
i2 = above(end);
t_rise = interp_crossing_local(t, v, i1 - 1, i1, level);
t_fall = interp_crossing_local(t, v, i2, i2 + 1, level);
if ~isfinite(t_rise) || ~isfinite(t_fall) || t_fall < t_rise
    center_t = NaN;
else
    center_t = 0.5 * (t_rise + t_fall);
end
end

function tc = interp_crossing_local(t, v, i_left, i_right, level)
if i_left < 1 || i_right > numel(t)
    tc = NaN;
    return;
end
v1 = v(i_left);
v2 = v(i_right);
if abs(v2 - v1) < eps
    tc = 0.5 * (t(i_left) + t(i_right));
else
    alpha = (level - v1) / (v2 - v1);
    alpha = min(max(alpha, 0), 1);
    tc = t(i_left) + alpha * (t(i_right) - t(i_left));
end
end

function [SemanticsTable, selected_events_per_rev] = infer_opr_events_per_revolution_local(opr_times, cfg, Step02_Setting)
rows = [];
for events_per_rev = Step02_Setting.opr_events_per_revolution_candidates(:).'
    if events_per_rev < 1 || numel(opr_times) <= events_per_rev
        continue;
    end
    period_s = opr_times((events_per_rev + 1):end) - opr_times(1:(end - events_per_rev));
    rpm = 60 ./ period_s;
    rpm = rpm(isfinite(rpm) & rpm > 0);
    if isempty(rpm)
        continue;
    end
    rpm_median = median(rpm);
    in_range = rpm_median >= Step02_Setting.expected_rpm_range(1) && ...
        rpm_median <= Step02_Setting.expected_rpm_range(2);
    rows = [rows; struct( ... %#ok<AGROW>
        'candidate_events_per_revolution', events_per_rev, ...
        'rpm_median', rpm_median, ...
        'rpm_min', min(rpm), ...
        'rpm_max', max(rpm), ...
        'rpm_std', std(rpm), ...
        'is_within_expected_rpm_range', in_range, ...
        'expected_rpm_min', Step02_Setting.expected_rpm_range(1), ...
        'expected_rpm_max', Step02_Setting.expected_rpm_range(2))];
end
if isempty(rows)
    error('Unable to infer OPR events per revolution from %d OPR events.', numel(opr_times));
end
SemanticsTable = struct2table(rows);
idx = find(SemanticsTable.is_within_expected_rpm_range, 1, 'first');
if isempty(idx)
    disp(SemanticsTable);
    error('No OPR events-per-revolution candidate gives RPM in expected range [%g, %g].', ...
        Step02_Setting.expected_rpm_range(1), Step02_Setting.expected_rpm_range(2));
end
selected_events_per_rev = SemanticsTable.candidate_events_per_revolution(idx);
SemanticsTable.selected = false(height(SemanticsTable), 1);
SemanticsTable.selected(idx) = true;
fprintf('>>> OPR events per revolution inferred as %d.\n', selected_events_per_rev);
end

function [OPRTable, Omega, OmegaTable] = build_opr_and_omega_tables_local(OPRRaw, events_per_rev, cfg, case_name)
opr_ref_time_s = OPRRaw.reference_time_s(:);
event_id = (1:numel(opr_ref_time_s)).';
rev_id = floor((event_id - 1) / events_per_rev) + 1;
event_in_rev = mod(event_id - 1, events_per_rev) + 1;
rev_period_s = nan(size(event_id));
rpm = nan(size(event_id));
valid = event_id + events_per_rev <= numel(event_id);
rev_period_s(valid) = opr_ref_time_s(event_id(valid) + events_per_rev) - opr_ref_time_s(valid);
rpm(valid) = 60 ./ rev_period_s(valid);

OPRTable = table( ...
    repmat(string(case_name), numel(event_id), 1), ...
    event_id, rev_id, event_in_rev, ...
    opr_ref_time_s, repmat(string(OPRRaw.reference_method), numel(event_id), 1), ...
    OPRRaw.center_time_s, OPRRaw.start_time_s, OPRRaw.end_time_s, OPRRaw.width_s, ...
    OPRRaw.center_offset_s, OPRRaw.peak_value, OPRRaw.center_quality, OPRRaw.file_id, ...
    rev_period_s, rpm, repmat(events_per_rev, numel(event_id), 1), ...
    'VariableNames', {'case_name','event_id','rev_id','event_in_rev', ...
    'opr_reference_time_s','opr_timing_reference', ...
    'opr_center_time_s','opr_start_time_s','opr_end_time_s','opr_width_s', ...
    'center_offset_s','opr_peak_value','center_quality','source_file_id', ...
    'rev_period_s','rpm','opr_events_per_revolution'});

[omega_time_s, omega_rad_s, omega_rpm, omega_period_s] = compute_speed_from_opr_local( ...
    opr_ref_time_s, events_per_rev);
Omega = struct();
Omega.omega_time_s = omega_time_s;
Omega.omega_rad_s = omega_rad_s;
Omega.omega_rpm = omega_rpm;
Omega.omega_zong = omega_rad_s;
Omega.metadata = struct( ...
    'source_mode', 'raw_dynamic_extraction', ...
    'opr_center_method', OPRRaw.center_method, ...
    'opr_timing_reference', OPRRaw.reference_method, ...
    'opr_events_per_revolution', events_per_rev, ...
    'rpm_recomputed_from_raw_opr_reference', true);

OmegaTable = table( ...
    repmat(string(case_name), numel(omega_time_s), 1), ...
    (1:numel(omega_time_s)).', omega_time_s(:), omega_period_s(:), ...
    omega_rad_s(:), omega_rpm(:), repmat(events_per_rev, numel(omega_time_s), 1), ...
    'VariableNames', {'case_name','omega_index','omega_time_s','rev_period_s', ...
    'omega_rad_s','omega_rpm','opr_events_per_revolution'});
end

function [omega_time_s, omega_rad_s, omega_rpm, rev_period_s] = compute_speed_from_opr_local(opr_times, events_per_rev)
if numel(opr_times) <= events_per_rev
    omega_time_s = [];
    omega_rad_s = [];
    omega_rpm = [];
    rev_period_s = [];
    return;
end
rev_period_s = opr_times((events_per_rev + 1):end) - opr_times(1:(end - events_per_rev));
omega_time_s = opr_times(1:(end - events_per_rev));
omega_rad_s = 2 * pi ./ rev_period_s;
omega_rpm = 60 ./ rev_period_s;
end

function T = check_step01_step02_consistency_local(Sensor_Config, cfg, case_name, SourceInventoryTable, events_per_rev)
rows = [];
for sid = cfg.sensor_ids(:).'
    fp = Sensor_Config.Fingerprints(sid);
    inv_idx = find(SourceInventoryTable.channel_id == sid, 1, 'first');
    file_count = NaN;
    if ~isempty(inv_idx)
        file_count = SourceInventoryTable.file_count(inv_idx);
    end
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'sensor_id', sid, ...
        'step01_fingerprint_length', numel(fp), ...
        'cfg_blades_num', cfg.blades_num, ...
        'fingerprint_length_matches_blades', numel(fp) == cfg.blades_num, ...
        'raw_file_count_for_sensor', file_count, ...
        'opr_events_per_revolution', events_per_rev, ...
        'sensor_in_cfg', ismember(sid, cfg.sensor_ids))];
end
T = struct2table(rows);
end

function ProbeRaw = extract_probe_from_raw_local(case_dir, file_ids, sid, cfg, Step02_Setting)
tail = [];
gap_points = Step02_Setting.initial_gap_points;
rows = [];
first_file_raw = [];
threshold = cfg.sensor_thresholds(sid);

for iFile = 1:numel(file_ids)
    file_id = file_ids(iFile);
    raw = load_raw_case_channel_local(case_dir, sid, file_id);
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if isempty(first_file_raw)
        first_file_raw = raw;
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gap_points] = segment_signal_local(raw, threshold, gap_points, true);
    rows = [rows; build_probe_rows_from_segments_local(raw, segments, file_id, sid, cfg, threshold)]; %#ok<AGROW>
end

if ~isempty(tail)
    [segments, ~, ~] = segment_signal_local(tail, threshold, gap_points, false);
    rows = [rows; build_probe_rows_from_segments_local(tail, segments, file_ids(end), sid, cfg, threshold)]; %#ok<AGROW>
end

if isempty(rows)
    warning('No probe pulses extracted for CH%d.', sid);
    ProbeRaw = empty_probe_raw_local(sid, first_file_raw);
    return;
end

[~, order] = sort([rows.arrival_time_s]);
rows = rows(order);
ProbeRaw = struct();
ProbeRaw.sensor_id = sid;
ProbeRaw.start_time_s = [rows.start_time_s].';
ProbeRaw.end_time_s = [rows.end_time_s].';
ProbeRaw.arrival_time_s = [rows.arrival_time_s].';
ProbeRaw.pulse_peak_value = [rows.pulse_peak_value].';
ProbeRaw.pulse_feature_value = [rows.pulse_feature_value].';
ProbeRaw.raw_peak_value = [rows.raw_peak_value].';
ProbeRaw.source_file_id = [rows.source_file_id].';
ProbeRaw.source_row = (1:numel(rows)).';
ProbeRaw.raw_first_file = first_file_raw;
fprintf('>>> CH%d raw probe pulses extracted: %d\n', sid, numel(ProbeRaw.arrival_time_s));
end

function ProbeRaw = empty_probe_raw_local(sid, first_file_raw)
ProbeRaw = struct();
ProbeRaw.sensor_id = sid;
ProbeRaw.start_time_s = [];
ProbeRaw.end_time_s = [];
ProbeRaw.arrival_time_s = [];
ProbeRaw.pulse_peak_value = [];
ProbeRaw.pulse_feature_value = [];
ProbeRaw.raw_peak_value = [];
ProbeRaw.source_file_id = [];
ProbeRaw.source_row = [];
ProbeRaw.raw_first_file = first_file_raw;
end

function rows = build_probe_rows_from_segments_local(raw, segments, file_id, sid, cfg, threshold)
rows = [];
for iSeg = 1:numel(segments.start_idx)
    a = segments.start_idx(iSeg);
    b = segments.end_idx(iSeg);
    [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = ...
        compute_half_area_arrival_local(raw, a, b, cfg.pinlv, threshold);
    rows = [rows; struct( ... %#ok<AGROW>
        'sensor_id', sid, ...
        'start_time_s', t_start, ...
        'end_time_s', t_end, ...
        'arrival_time_s', t_arrival, ...
        'pulse_peak_value', feat_raw, ...
        'pulse_feature_value', feat, ...
        'raw_peak_value', feat_fitted, ...
        'source_file_id', file_id)];
end
end

function raw = load_raw_case_channel_local(case_dir, channel_id, file_id)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    raw = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = [];
for i = 1:numel(fn)
    value = loaded.(fn{i});
    if isnumeric(value) && ismatrix(value) && size(value, 2) >= 2
        raw = value(:, 1:2);
        break;
    end
end
if isempty(raw)
    return;
end
raw(raw(:, 1) == 0, :) = [];
raw = raw(all(isfinite(raw), 2), :);
end

function [segments, tail, next_gap] = segment_signal_local(raw, threshold, gap_points, keep_last_as_tail)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];
next_gap = gap_points;
if isempty(raw)
    return;
end

sig = raw(:, 2);
try
    sig_smooth = smooth(sig, 16);
catch
    sig_smooth = smoothdata(sig, 'movmean', 16);
end

chase = find(sig_smooth > threshold);
if isempty(chase)
    return;
end

sample_order = raw(chase, 1);
seg_starts = 1;
seg_ends = [];
for ii = 1:(numel(sample_order) - 1)
    c = sample_order(ii + 1) - sample_order(ii);
    if c > next_gap
        seg_ends(end + 1, 1) = ii; %#ok<AGROW>
        seg_starts(end + 1, 1) = ii + 1; %#ok<AGROW>
        next_gap = 0.6 * c;
    end
end
seg_ends(end + 1, 1) = numel(sample_order);

if keep_last_as_tail
    tail_point = chase(seg_starts(end)) - floor(next_gap / 2);
    if tail_point > 0 && tail_point < size(raw, 1)
        tail = raw(tail_point:end, :);
    end
end

if keep_last_as_tail
    if numel(seg_starts) < 2
        return;
    end
    comp_starts = seg_starts(1:end-1);
    comp_ends = seg_ends(1:end-1);
else
    comp_starts = seg_starts;
    comp_ends = seg_ends;
end

segments.start_idx = chase(comp_starts);
segments.end_idx = chase(comp_ends);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end

function [t_start, t_end, t_arrival, feat, feat_raw, feat_fitted] = compute_half_area_arrival_local(raw, a, b, pinlv, threshold)
expand_pts = floor((b - a + 1) * 0.3);
a = max(1, a - expand_pts);
b = min(size(raw, 1), b + expand_pts);
t = raw(a:b, 1) / pinlv;
v = raw(a:b, 2);
t_start = t(1);
t_end = t(end);
feat_raw = max(v);
feat_fitted = fit_peak_feature_local(t, v, threshold);
if isfinite(feat_fitted)
    feat = feat_fitted;
else
    feat = feat_raw;
end
if numel(t) < 2
    t_arrival = t_start;
    return;
end
dt = median(diff(t));
area = cumsum(max(v, 0) * dt);
target = 0.5 * area(end);
idx = find(area >= target, 1, 'first');
if isempty(idx)
    t_arrival = t_start;
elseif idx == 1
    t_arrival = t(1);
else
    t_arrival = 0.5 * (t(idx - 1) + t(idx));
end
end

function peak_feature = fit_peak_feature_local(t_seg, v_seg, threshold)
peak_feature = NaN;
if isempty(t_seg) || isempty(v_seg)
    return;
end
smooth_seg = smooth_signal_for_peak_local(v_seg);
valid_mask = smooth_seg > threshold;
if sum(valid_mask) < 4
    peak_feature = max(smooth_seg);
    return;
end
t_fit = t_seg(valid_mask);
v_fit = smooth_seg(valid_mask);
mu = mean(t_fit);
order = min(3, numel(unique(t_fit)) - 1);
if order < 1
    peak_feature = max(v_fit);
    return;
end
try
    warn_state_1 = warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    warn_state_2 = warning('off', 'MATLAB:polyfit:PolyNotUnique');
    warn_state_3 = warning('off', 'MATLAB:singularMatrix');
    warn_state_4 = warning('off', 'MATLAB:nearlySingularMatrix');
    cleanup_obj = onCleanup(@() restore_polyfit_warnings_local( ...
        warn_state_1, warn_state_2, warn_state_3, warn_state_4)); %#ok<NASGU>
    p = polyfit(t_fit - mu, v_fit, order);
    tt = linspace(min(t_fit) - mu, max(t_fit) - mu, 100);
    peak_feature = max(polyval(p, tt));
catch
    peak_feature = max(smooth_seg);
end
end

function smooth_v = smooth_signal_for_peak_local(v)
try
    if numel(v) >= 21
        smooth_v = sgolayfilt(v, 3, 21);
    elseif numel(v) >= 5
        span = max(3, 2 * floor(numel(v) / 4) + 1);
        smooth_v = smoothdata(v, 'movmean', span);
    else
        smooth_v = v;
    end
catch
    smooth_v = v;
end
end

function restore_polyfit_warnings_local(w1, w2, w3, w4)
warning(w1);
warning(w2);
warning(w3);
warning(w4);
end

function [PulseTable, SensorQuality, RevQuality, ProbeDiag] = assign_blades_from_step01_local( ...
    ProbeRaw, OPRTable, Sensor_Config, cfg, Step02_Setting, case_name, sid)
n = numel(ProbeRaw.arrival_time_s);
if n == 0
    PulseTable = empty_pulse_table_local();
    SensorQuality = empty_sensor_quality_row_local(case_name, sid);
    RevQuality = empty_rev_quality_table_local();
    ProbeDiag = struct();
    return;
end

[prev_opr_event_id, rev_id, event_in_rev, relative_angle_deg] = assign_revolution_by_opr_local( ...
    ProbeRaw.arrival_time_s, OPRTable);
ref_fp = Sensor_Config.Fingerprints(sid);
ref_fp = ref_fp(:).';
[selected_shift, sensor_match, rev_match_map] = choose_sensor_shift_local( ...
    ProbeRaw, rev_id, ref_fp, cfg, Step02_Setting);

assigned_blade_id = nan(n, 1);
assignment_quality = nan(n, 1);
assignment_validity_flag = false(n, 1);
quality_flags = repmat({''}, n, 1);
assignment_source = repmat({'unassigned'}, n, 1);
quality_source = repmat({'same_revolution_fingerprint_corr'}, n, 1);

valid_rev_ids = unique(rev_id(isfinite(rev_id)));
for iRev = 1:numel(valid_rev_ids)
    rid = valid_rev_ids(iRev);
    idx = find(rev_id == rid);
    [~, order] = sort(ProbeRaw.arrival_time_s(idx));
    idx = idx(order);
    map_idx = find(rev_match_map.rev_id == rid, 1, 'first');
    if isempty(map_idx) || numel(idx) ~= cfg.blades_num || ~isfinite(selected_shift)
        quality_flags(idx) = {'incomplete_or_unmatched_revolution'};
        continue;
    end
    blade_seq = mod((0:(cfg.blades_num - 1)) - selected_shift, cfg.blades_num) + 1;
    assigned_blade_id(idx) = blade_seq(:);
    assignment_quality(idx) = rev_match_map.best_corr(map_idx);
    is_valid = rev_match_map.best_corr(map_idx) >= Step02_Setting.fingerprint_min_corr && ...
        rev_match_map.corr_gap(map_idx) >= Step02_Setting.fingerprint_min_corr_gap && ...
        rev_match_map.best_shift(map_idx) == selected_shift;
    assignment_validity_flag(idx) = is_valid;
    if is_valid
        quality_flags(idx) = {'ok'};
        assignment_source(idx) = {'selected_shift_from_step01_fingerprint'};
    else
        quality_flags(idx) = {'low_corr_or_shift_mismatch'};
        assignment_source(idx) = {'assigned_but_not_validated'};
    end
end

is_valid = isfinite(assigned_blade_id) & assignment_validity_flag;
PulseTable = table( ...
    repmat(string(case_name), n, 1), repmat(sid, n, 1), ProbeRaw.source_row(:), ...
    rev_id(:), event_in_rev(:), prev_opr_event_id(:), ...
    ProbeRaw.start_time_s(:), ProbeRaw.end_time_s(:), ProbeRaw.arrival_time_s(:), ...
    relative_angle_deg(:), assigned_blade_id(:), ProbeRaw.pulse_peak_value(:), ...
    ProbeRaw.pulse_feature_value(:), ProbeRaw.raw_peak_value(:), assignment_quality(:), ...
    assignment_validity_flag(:), is_valid(:), string(quality_flags(:)), ...
    string(assignment_source(:)), string(quality_source(:)), ProbeRaw.source_file_id(:), ...
    repmat(string('raw_dynamic_extraction'), n, 1), ...
    'VariableNames', {'case_name','sensor_id','source_row','rev_id','event_in_rev', ...
    'prev_opr_event_id','pulse_start_time_s','pulse_end_time_s','arrival_time_s', ...
    'relative_angle_deg','assigned_blade_id','pulse_peak_value','pulse_feature_value', ...
    'raw_peak_value','assignment_quality','assignment_validity_flag','is_valid', ...
    'quality_flags','assignment_source','assignment_quality_source','source_file_id','source_mode'});

SensorQuality = table( ...
    string(case_name), sid, n, numel(valid_rev_ids), sensor_match.complete_revolution_count, ...
    selected_shift, sensor_match.selected_shift_source, sensor_match.best_corr_max, ...
    sensor_match.best_corr_median, sensor_match.second_best_corr_median, ...
    sensor_match.corr_gap_median, sensor_match.valid_revolution_count, ...
    sum(isfinite(assigned_blade_id)), sum(is_valid), sensor_match.quality_status, ...
    'VariableNames', {'case_name','sensor_id','pulse_count','revolution_count', ...
    'complete_revolution_count','selected_shift','selected_shift_source','best_corr_max', ...
    'best_corr_median','second_best_corr_median','corr_gap_median','valid_revolution_count', ...
    'assigned_pulse_count','valid_assigned_pulse_count','quality_status'});

RevQuality = build_rev_quality_table_local(case_name, sid, rev_id, assigned_blade_id, ...
    assignment_validity_flag, rev_match_map, cfg);

ProbeDiag = struct();
ProbeDiag.reference_fingerprint = ref_fp;
ProbeDiag.selected_shift = selected_shift;
ProbeDiag.sensor_match = sensor_match;
ProbeDiag.rev_match_map = rev_match_map;
end

function [prev_opr_event_id, rev_id, event_in_rev, relative_angle_deg] = assign_revolution_by_opr_local(arrival_times, OPRTable)
n = numel(arrival_times);
prev_opr_event_id = nan(n, 1);
rev_id = nan(n, 1);
event_in_rev = nan(n, 1);
relative_angle_deg = nan(n, 1);
opr_times = resolve_opr_reference_times_from_table_local(OPRTable);
events_per_rev = OPRTable.opr_events_per_revolution(1);
for i = 1:n
    t = arrival_times(i);
    idx = find(opr_times <= t, 1, 'last');
    if isempty(idx)
        continue;
    end
    prev_opr_event_id(i) = OPRTable.event_id(idx);
    rev_id(i) = OPRTable.rev_id(idx);
    event_in_rev(i) = OPRTable.event_in_rev(idx);
    idx_next_rev = idx + events_per_rev;
    if idx_next_rev <= numel(opr_times)
        rev_period = opr_times(idx_next_rev) - opr_times(idx);
        if rev_period > 0
            relative_angle_deg(i) = 360 * (t - opr_times(idx)) / rev_period;
        end
    end
end
end

function opr_times = resolve_opr_reference_times_from_table_local(OPRTable)
if ismember('opr_reference_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_reference_time_s;
elseif ismember('opr_start_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_start_time_s;
else
    opr_times = OPRTable.opr_center_time_s;
end
end

function [selected_shift, sensor_match, rev_match_map] = choose_sensor_shift_local(ProbeRaw, rev_id, ref_fp, cfg, Step02_Setting)
valid_rev_ids = unique(rev_id(isfinite(rev_id)));
rows = [];
for iRev = 1:numel(valid_rev_ids)
    rid = valid_rev_ids(iRev);
    idx = find(rev_id == rid);
    [~, order] = sort(ProbeRaw.arrival_time_s(idx));
    idx = idx(order);
    if numel(idx) ~= cfg.blades_num
        rows = [rows; struct('rev_id', rid, 'pulse_count', numel(idx), 'best_shift', NaN, ... %#ok<AGROW>
            'best_corr', NaN, 'second_best_corr', NaN, 'corr_gap', NaN, 'is_complete', false)];
        continue;
    end
    feats = ProbeRaw.pulse_feature_value(idx).';
    corr_by_shift = nan(cfg.blades_num, 1);
    for shift = 0:(cfg.blades_num - 1)
        corr_by_shift(shift + 1) = safe_corr_local(ref_fp(:), circshift(feats(:), -shift));
    end
    [best_corr, best_pos] = max(corr_by_shift);
    second_best = max(corr_by_shift((1:numel(corr_by_shift)) ~= best_pos));
    rows = [rows; struct('rev_id', rid, 'pulse_count', numel(idx), 'best_shift', best_pos - 1, ... %#ok<AGROW>
        'best_corr', best_corr, 'second_best_corr', second_best, ...
        'corr_gap', best_corr - second_best, 'is_complete', true)];
end

if isempty(rows)
    rev_match_map = table();
    selected_shift = NaN;
    sensor_match = sensor_match_default_local(0, 0, selected_shift, "no_revolutions");
    return;
end
rev_match_map = struct2table(rows);
complete = rev_match_map.is_complete & isfinite(rev_match_map.best_corr);
quality_mask = complete & rev_match_map.best_corr >= Step02_Setting.fingerprint_min_corr & ...
    rev_match_map.corr_gap >= Step02_Setting.fingerprint_min_corr_gap;

if any(quality_mask)
    selected_shift = robust_mode_shift_local(rev_match_map.best_shift(quality_mask), ...
        rev_match_map.best_corr(quality_mask));
    selected_source = "quality_revolution_vote";
    quality_status = "good";
elseif any(complete)
    [~, idx] = max(rev_match_map.best_corr);
    selected_shift = rev_match_map.best_shift(idx);
    selected_source = "best_available_revolution_below_threshold";
    quality_status = "weak";
else
    selected_shift = NaN;
    selected_source = "not_available";
    quality_status = "failed_no_complete_revolution";
end

sensor_match = sensor_match_default_local(sum(complete), sum(quality_mask), selected_shift, quality_status);
sensor_match.selected_shift_source = selected_source;
if any(complete)
    sensor_match.best_corr_max = max(rev_match_map.best_corr(complete));
    sensor_match.best_corr_median = median(rev_match_map.best_corr(complete));
    sensor_match.second_best_corr_median = median(rev_match_map.second_best_corr(complete));
    sensor_match.corr_gap_median = median(rev_match_map.corr_gap(complete));
end
end

function selected_shift = robust_mode_shift_local(shifts, corr_values)
unique_shifts = unique(shifts(isfinite(shifts)));
score = nan(size(unique_shifts));
for i = 1:numel(unique_shifts)
    mask = shifts == unique_shifts(i);
    score(i) = sum(mask) + 0.01 * median(corr_values(mask));
end
[~, idx] = max(score);
selected_shift = unique_shifts(idx);
end

function s = sensor_match_default_local(complete_count, valid_count, selected_shift, quality_status)
s = struct();
s.complete_revolution_count = complete_count;
s.valid_revolution_count = valid_count;
s.selected_shift = selected_shift;
s.selected_shift_source = "not_available";
s.best_corr_max = NaN;
s.best_corr_median = NaN;
s.second_best_corr_median = NaN;
s.corr_gap_median = NaN;
s.quality_status = quality_status;
end

function RevQuality = build_rev_quality_table_local(case_name, sid, rev_id, assigned_blade_id, assignment_validity_flag, rev_match_map, cfg)
valid_rev_ids = unique(rev_id(isfinite(rev_id)));
rows = [];
for i = 1:numel(valid_rev_ids)
    rid = valid_rev_ids(i);
    idx = find(rev_id == rid);
    map_idx = [];
    if ~isempty(rev_match_map)
        map_idx = find(rev_match_map.rev_id == rid, 1, 'first');
    end
    best_shift = NaN; best_corr = NaN; second_best = NaN; corr_gap = NaN;
    if ~isempty(map_idx)
        best_shift = rev_match_map.best_shift(map_idx);
        best_corr = rev_match_map.best_corr(map_idx);
        second_best = rev_match_map.second_best_corr(map_idx);
        corr_gap = rev_match_map.corr_gap(map_idx);
    end
    assigned_count = sum(isfinite(assigned_blade_id(idx)));
    valid_count = sum(assignment_validity_flag(idx));
    duplicate_count = assigned_count - numel(unique(assigned_blade_id(idx(isfinite(assigned_blade_id(idx))))));
    rows = [rows; struct( ... %#ok<AGROW>
        'case_name', string(case_name), ...
        'sensor_id', sid, ...
        'rev_id', rid, ...
        'pulse_count', numel(idx), ...
        'assigned_count', assigned_count, ...
        'valid_assigned_count', valid_count, ...
        'missing_count', cfg.blades_num - assigned_count, ...
        'duplicate_blade_count', duplicate_count, ...
        'best_shift', best_shift, ...
        'best_corr', best_corr, ...
        'second_best_corr', second_best, ...
        'corr_gap', corr_gap, ...
        'is_complete_revolution', numel(idx) == cfg.blades_num)];
end
if isempty(rows)
    RevQuality = empty_rev_quality_table_local();
else
    RevQuality = struct2table(rows);
end
end

function T = empty_pulse_table_local()
T = table(string.empty(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), false(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), zeros(0,1), string.empty(0,1), ...
    'VariableNames', {'case_name','sensor_id','source_row','rev_id','event_in_rev', ...
    'prev_opr_event_id','pulse_start_time_s','pulse_end_time_s','arrival_time_s', ...
    'relative_angle_deg','assigned_blade_id','pulse_peak_value','pulse_feature_value', ...
    'raw_peak_value','assignment_quality','assignment_validity_flag','is_valid', ...
    'quality_flags','assignment_source','assignment_quality_source','source_file_id','source_mode'});
end

function T = empty_sensor_quality_row_local(case_name, sid)
T = table(string(case_name), sid, 0, 0, 0, NaN, string("not_available"), NaN, NaN, NaN, NaN, 0, 0, 0, string("no_pulses"), ...
    'VariableNames', {'case_name','sensor_id','pulse_count','revolution_count', ...
    'complete_revolution_count','selected_shift','selected_shift_source','best_corr_max', ...
    'best_corr_median','second_best_corr_median','corr_gap_median','valid_revolution_count', ...
    'assigned_pulse_count','valid_assigned_pulse_count','quality_status'});
end

function T = empty_rev_quality_table_local()
T = table(string.empty(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
    'VariableNames', {'case_name','sensor_id','rev_id','pulse_count','assigned_count', ...
    'valid_assigned_count','missing_count','duplicate_blade_count','best_shift', ...
    'best_corr','second_best_corr','corr_gap','is_complete_revolution'});
end

function c = safe_corr_local(a, b)
if numel(a) ~= numel(b) || numel(a) < 2
    c = -inf;
    return;
end
if all(abs(a - a(1)) < eps) || all(abs(b - b(1)) < eps)
    c = -inf;
    return;
end
r = corrcoef(a, b);
if numel(r) < 4 || ~isfinite(r(1, 2))
    c = -inf;
else
    c = r(1, 2);
end
end

function Summary = build_step02_summary_local(cfg, case_name, DynamicPulseTable, SensorNumberingQuality, RevolutionNumberingQuality, OPRTable, OmegaTable)
valid_rpm = OmegaTable.omega_rpm(isfinite(OmegaTable.omega_rpm));
if isempty(valid_rpm)
    rpm_min = NaN; rpm_median = NaN; rpm_max = NaN;
else
    rpm_min = min(valid_rpm);
    rpm_median = median(valid_rpm);
    rpm_max = max(valid_rpm);
end
Summary = table( ...
    string(case_name), string('raw_dynamic_extraction'), height(OPRTable), ...
    OPRTable.opr_events_per_revolution(1), height(DynamicPulseTable), ...
    sum(DynamicPulseTable.is_valid), rpm_min, rpm_median, rpm_max, ...
    sum(SensorNumberingQuality.quality_status == "good"), ...
    height(RevolutionNumberingQuality), cfg.blades_num, ...
    'VariableNames', {'case_name','source_mode','opr_event_count', ...
    'opr_events_per_revolution','pulse_count','valid_pulse_count', ...
    'rpm_min','rpm_median','rpm_max','good_sensor_count', ...
    'revolution_quality_rows','blades_num'});
end

function Meta = build_step02_metadata_local(cfg, Step02_Setting, case_name, case_dir, sensor_config_file, Sensor_Config, OPRRaw, events_per_rev, ProbeDiagnostics)
Meta = struct();
Meta.Dataset = cfg.dataset;
Meta.CaseName = case_name;
Meta.CreatedAt = string(datetime('now'));
Meta.SourceMode = Step02_Setting.source_mode;
Meta.RawCaseDir = case_dir;
Meta.SensorConfigFile = sensor_config_file;
Meta.SensorIds = cfg.sensor_ids;
Meta.OPRId = cfg.opr_id;
Meta.BladesNum = cfg.blades_num;
Meta.SampleRateHz = cfg.sample_rate_hz;
Meta.OPRCenterMethod = Step02_Setting.opr_center_method;
Meta.OPRTimingReference = OPRRaw.reference_method;
Meta.OPRReferenceTimeColumn = 'opr_reference_time_s';
Meta.OPREventsPerRevolution = events_per_rev;
Meta.AllowOldDynamicOutputsAsSource = Step02_Setting.allow_old_dynamic_outputs_as_source;
Meta.SensorConfigReference = get_field_or_default_local(Sensor_Config, ...
    'Standard_Relative_Angles_Reference', 'unknown');
Meta.OPRRawPulseCount = numel(OPRRaw.center_time_s);
Meta.OPRCenterQualityMedian = median(OPRRaw.center_quality);
Meta.ProbeDiagnostics = ProbeDiagnostics;
Meta.Note = ['Independent Step02 raw extraction. OPR reference time is the configured ', ...
    'timing reference, while center/start/end are retained for audit. Old dynamic jiluOPR/omega/jilublade ', ...
    'files are not read as sources; compatibility files saved here are newly generated.'];
end

function save_compatibility_outputs_local(out_dir, OPRRaw, Omega, DynamicPulseTable, cfg, Step02Meta)
jiluOPR = [OPRRaw.reference_time_s(:), OPRRaw.start_time_s(:), ...
    OPRRaw.end_time_s(:), OPRRaw.center_time_s(:)]; %#ok<NASGU>
metadata = struct();
metadata.method = Step02Meta.OPRTimingReference;
metadata.center_method = Step02Meta.OPRCenterMethod;
metadata.source_mode = Step02Meta.SourceMode;
metadata.opr_events_per_revolution = Step02Meta.OPREventsPerRevolution;
metadata.columns = {'reference_time_s','start_time_s','end_time_s','center_time_s'};
metadata.raw_extraction = true;
save(fullfile(out_dir, 'jiluOPR.mat'), 'jiluOPR', 'metadata');

omega_time_s = Omega.omega_time_s; %#ok<NASGU>
omega_rad_s = Omega.omega_rad_s; %#ok<NASGU>
omega_rpm = Omega.omega_rpm; %#ok<NASGU>
omega_zong = Omega.omega_zong; %#ok<NASGU>
metadata = Omega.metadata; %#ok<NASGU>
save(fullfile(out_dir, 'omega.mat'), 'omega_zong', 'omega_time_s', 'omega_rpm', 'omega_rad_s', 'metadata');

for sid = cfg.sensor_ids(:).'
    rows = DynamicPulseTable(DynamicPulseTable.sensor_id == sid, :);
    jilublade = [rows.pulse_start_time_s, rows.pulse_end_time_s, rows.arrival_time_s, rows.assigned_blade_id]; %#ok<NASGU>
    jilublade_standard = [rows.arrival_time_s, rows.pulse_peak_value, rows.pulse_feature_value, ...
        rows.assigned_blade_id, rows.rev_id, rows.assignment_quality]; %#ok<NASGU>
    save(fullfile(out_dir, sprintf('jilublade_probe%d.mat', sid)), 'jilublade');
    save(fullfile(out_dir, sprintf('jilublade_probe%d_standard.mat', sid)), 'jilublade_standard');
end
end

function plot_step02_diagnostics_local(cfg, Step02_Setting, case_name, OPRTable, OmegaTable, DynamicPulseTable, SensorNumberingQuality, RevolutionNumberingQuality, fig_dir, show_figures, save_figures)
visibility = 'off';
if show_figures
    visibility = 'on';
end

fig = figure('Name', 'Step02 RPM Trend', 'Color', 'w', 'Visible', visibility, 'Position', [80 80 1200 420]);
ax = axes(fig);
plot(ax, OmegaTable.omega_time_s, OmegaTable.omega_rpm, 'k-', 'LineWidth', 1.1);
grid(ax, 'on');
xlabel(ax, 'Time (s)');
ylabel(ax, 'RPM');
title(ax, sprintf('Raw OPR-derived RPM - %s', case_name), 'Interpreter', 'none');
save_figure_local(fig, fig_dir, 'Step02_RPM_Trend_20251222', save_figures);

fig = figure('Name', 'Step02 OPR Center Check', 'Color', 'w', 'Visible', visibility, 'Position', [100 100 1200 520]);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
plot(ax, OPRTable.opr_center_time_s, OPRTable.center_offset_s * 1e6, '.', 'MarkerSize', 4);
grid(ax, 'on');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Center offset (us)');
title(ax, 'OPR center minus midpoint(start,end)');
ax = nexttile;
plot(ax, OPRTable.opr_center_time_s, OPRTable.opr_width_s * 1e6, '.', 'MarkerSize', 4);
grid(ax, 'on');
xlabel(ax, 'Time (s)');
ylabel(ax, 'OPR width (us)');
title(ax, 'OPR pulse width from raw waveform');
save_figure_local(fig, fig_dir, 'Step02_OPR_Center_Check_20251222', save_figures);

fig = figure('Name', 'Step02 Blade ID Sequence', 'Color', 'w', 'Visible', visibility, 'Position', [120 120 1350 780]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
colors = lines(cfg.blades_num);
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = DynamicPulseTable(DynamicPulseTable.sensor_id == sid, :);
    rows_show = rows(1:min(height(rows), Step02_Setting.max_plot_points), :);
    ax = nexttile;
    hold(ax, 'on');
    for bid = 1:cfg.blades_num
        mask = rows_show.assigned_blade_id == bid;
        scatter(ax, rows_show.arrival_time_s(mask), rows_show.assigned_blade_id(mask), ...
            8, colors(bid, :), 'filled');
    end
    grid(ax, 'on');
    ylim(ax, [0.5 cfg.blades_num + 0.5]);
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, sprintf('CH%d blade', sid));
    title(ax, sprintf('CH%d first %d raw-extracted labeled pulses', sid, height(rows_show)));
end
save_figure_local(fig, fig_dir, 'Step02_BladeID_Sequence_20251222', save_figures);

fig = figure('Name', 'Step02 Revolution Coverage', 'Color', 'w', 'Visible', visibility, 'Position', [140 140 1250 720]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    q = RevolutionNumberingQuality(RevolutionNumberingQuality.sensor_id == sid, :);
    q_show = q(1:min(height(q), Step02_Setting.max_plot_revolutions), :);
    ax = nexttile;
    bar(ax, q_show.rev_id, q_show.assigned_count, 0.8);
    yline(ax, cfg.blades_num, 'r--', 'LineWidth', 1.1);
    grid(ax, 'on');
    xlabel(ax, 'Revolution ID');
    ylabel(ax, 'Assigned pulses');
    title(ax, sprintf('CH%d raw extraction revolution coverage', sid));
end
save_figure_local(fig, fig_dir, 'Step02_Revolution_Coverage_20251222', save_figures);

fig = figure('Name', 'Step02 Fingerprint Match Summary', 'Color', 'w', 'Visible', visibility, 'Position', [160 160 1000 450]);
ax = axes(fig);
bar(ax, SensorNumberingQuality.sensor_id, SensorNumberingQuality.best_corr_median, 0.6);
grid(ax, 'on');
set(ax, 'XTick', cfg.sensor_ids);
xlabel(ax, 'Sensor ID');
ylabel(ax, 'Median complete-rev correlation');
ylim(ax, [0 1.05]);
title(ax, 'Step01 fingerprint match quality on raw dynamic pulses');
save_figure_local(fig, fig_dir, 'Step02_Fingerprint_Match_20251222', save_figures);

fig = figure('Name', 'Step02 Pulse Feature Trend', 'Color', 'w', 'Visible', visibility, 'Position', [180 180 1350 780]);
tiledlayout(fig, numel(cfg.sensor_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    rows = DynamicPulseTable(DynamicPulseTable.sensor_id == sid, :);
    rows_show = rows(1:min(height(rows), Step02_Setting.max_plot_points), :);
    ax = nexttile;
    plot(ax, rows_show.arrival_time_s, rows_show.pulse_feature_value, '.', 'MarkerSize', 4);
    grid(ax, 'on');
    xlabel(ax, 'Arrival time (s)');
    ylabel(ax, 'Pulse feature');
    title(ax, sprintf('CH%d raw pulse feature trend', sid));
end
save_figure_local(fig, fig_dir, 'Step02_Pulse_Feature_Trend_20251222', save_figures);
end

function save_figure_local(fig, fig_dir, base_name, save_figures)
if ~save_figures
    return;
end
ensure_dir_local(fig_dir);
png_file = fullfile(fig_dir, [base_name, '.png']);
pdf_file = fullfile(fig_dir, [base_name, '.pdf']);
try
    exportgraphics(fig, png_file, 'Resolution', 300);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    saveas(fig, png_file);
    saveas(fig, pdf_file);
end
end

function value = get_field_or_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
    value = S.(field_name);
else
    value = default_value;
end
end
