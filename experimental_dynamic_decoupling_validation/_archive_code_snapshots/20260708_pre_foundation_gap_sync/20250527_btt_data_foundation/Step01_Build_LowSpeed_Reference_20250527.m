clc; clear; close all;

%STEP01_BUILD_LOWSPEED_REFERENCE_20250527 Build low-speed BTT reference.
%
% Step01 output is a method-neutral reference for later dynamic extraction:
%   - per-sensor blade fingerprints
%   - target start indices
%   - standard relative angle from previous OPR to each blade passage
%
% 2026-07 revision:
%   - The internal Blade 1 convention is explicit. By default Blade 1 is
%     defined as the largest reference-sensor pulse in the first few low-speed
%     revolutions; this is valid for relative BTT numbering but not a physical
%     blade label unless a separate offset is supplied.
%   - OPR/probe timing metadata from the extractor is preserved.
%   - Fingerprint matching reports best/second-best correlation and gap.
%   - Standard relative angles are saved with mean/median/std/IQR/count.

%% Step01 run settings
show_plots = true;

cfg = BTTDataConfig_20250527();
cfg = fill_missing_step01_defaults_local(cfg);

%% Step01 reference-building settings
low_speed_blade1_rule = cfg.low_speed_blade1_rule;
reference_search_revs = cfg.reference_search_revs;
max_alignment_search_starts = cfg.max_alignment_search_starts;
fingerprint_min_corr = cfg.fingerprint_min_corr;
fingerprint_min_corr_gap = cfg.fingerprint_min_corr_gap;
fingerprint_quality_policy = cfg.fingerprint_quality_policy;
standard_angle_value = cfg.standard_angle_value;

cfg.low_speed_blade1_rule = low_speed_blade1_rule;
cfg.reference_search_revs = reference_search_revs;
cfg.max_alignment_search_starts = max_alignment_search_starts;
cfg.fingerprint_min_corr = fingerprint_min_corr;
cfg.fingerprint_min_corr_gap = fingerprint_min_corr_gap;
cfg.fingerprint_quality_policy = fingerprint_quality_policy;
cfg.standard_angle_value = standard_angle_value;

case_data = extract_low_speed_btt_features_20250527(cfg, cfg.sensor_ids);

if ~exist(cfg.step01_output_dir, 'dir')
    mkdir(cfg.step01_output_dir);
end
if cfg.save_figures && ~exist(cfg.step01_figure_dir, 'dir')
    mkdir(cfg.step01_figure_dir);
end

anchor_sid = cfg.reference_sensor_id;
if ~ismember(anchor_sid, cfg.sensor_ids)
    error('Reference sensor CH%d is not included in cfg.sensor_ids.', anchor_sid);
end
anchor_feat = case_data.channels(anchor_sid).peak_features(:);
if numel(anchor_feat) < cfg.blades_num
    error('Reference sensor CH%d does not have enough pulses.', anchor_sid);
end

[anchor_start_idx, anchor_peak_idx, anchor_peak_value, blade_id_definition] = ...
    choose_reference_blade1_start_local(anchor_feat, cfg);
ref_fp = anchor_feat(anchor_start_idx:(anchor_start_idx + cfg.blades_num - 1));

fingerprints = containers.Map('KeyType', 'double', 'ValueType', 'any');
target_indices = containers.Map('KeyType', 'double', 'ValueType', 'double');
summary_rows = repmat(struct( ...
    'sensor_id', NaN, ...
    'start_index', NaN, ...
    'best_corr', NaN, ...
    'second_best_corr', NaN, ...
    'corr_gap', NaN, ...
    'quality_status', '', ...
    'pulse_count', NaN), numel(cfg.sensor_ids), 1);

for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    if numel(feats) < cfg.blades_num
        error('CH%d does not have enough low-speed pulses for fingerprint alignment.', sid);
    end

    [best_idx, best_corr, second_best_corr, second_best_idx, corr_gap] = ...
        match_reference_fingerprint_local(feats, ref_fp, cfg.blades_num, cfg.max_alignment_search_starts); %#ok<ASGLU>

    [quality_status, quality_msg] = evaluate_fingerprint_quality_local( ...
        sid, best_corr, corr_gap, cfg);
    if ~isempty(quality_msg)
        switch lower(strtrim(cfg.fingerprint_quality_policy))
            case 'error'
                error('%s', quality_msg);
            case 'warn'
                warning('%s', quality_msg);
        end
    end

    fingerprints(sid) = feats(best_idx:(best_idx + cfg.blades_num - 1)).';
    target_indices(sid) = best_idx;

    summary_rows(i).sensor_id = sid;
    summary_rows(i).start_index = best_idx;
    summary_rows(i).best_corr = best_corr;
    summary_rows(i).second_best_corr = second_best_corr;
    summary_rows(i).corr_gap = corr_gap;
    summary_rows(i).quality_status = quality_status;
    summary_rows(i).pulse_count = numel(feats);
end

[std_angles, angle_stats, angle_quality_table] = build_standard_relative_angle_stats_local( ...
    cfg, case_data, target_indices);

Sensor_Config = struct();
Sensor_Config.Dataset = cfg.dataset;
Sensor_Config.ReferenceCase = cfg.low_speed_case;
Sensor_Config.Sensor_IDs = cfg.sensor_ids;
Sensor_Config.Candidate_Sensor_IDs = cfg.candidate_sensor_ids;
Sensor_Config.Capacitance_IDs = cfg.capacitance_ids;
Sensor_Config.Eddy_Current_IDs = cfg.eddy_current_ids;
Sensor_Config.OPR_ID = cfg.opr_id;
Sensor_Config.Blades_Num = cfg.blades_num;
Sensor_Config.Sample_Rate_Hz = cfg.sample_rate_hz;
Sensor_Config.Pinlv = cfg.sample_rate_hz;
Sensor_Config.R_Tip_mm = cfg.r_tip_mm;
Sensor_Config.Gap_Points = cfg.gap_points;
Sensor_Config.OPR_Threshold = cfg.opr_threshold;
Sensor_Config.Sensor_Thresholds = cfg.sensor_thresholds;
Sensor_Config.Reference_Sensor_ID = cfg.reference_sensor_id;
Sensor_Config.Fingerprints = fingerprints;
Sensor_Config.Reference_Fingerprint = ref_fp(:).';
Sensor_Config.Target_Indices = target_indices;

% Reference-frame and timing metadata.
Sensor_Config.OPR_Timing_Method = case_data.opr_timing_method;
Sensor_Config.Probe_Arrival_Method = case_data.probe_arrival_method;
Sensor_Config.Time_Index_Mode = case_data.time_index_mode;
Sensor_Config.Time_Index_Mode_Detected = case_data.time_index_mode_detected;
Sensor_Config.Standard_Relative_Angles_Reference = case_data.opr_timing_method;
Sensor_Config.Standard_Relative_Angles_SelectedStatistic = cfg.standard_angle_value;

% Blade-numbering metadata.
Sensor_Config.Blade_ID_Definition = blade_id_definition;
Sensor_Config.Blade_ID_Rule = cfg.low_speed_blade1_rule;
Sensor_Config.Blade_ID_Offset = cfg.blade_id_offset;
Sensor_Config.Reference_Blade1_StartIndex = anchor_start_idx;
Sensor_Config.Reference_Blade1_PeakIndexInSearch = anchor_peak_idx;
Sensor_Config.Reference_Blade1_PeakFeature = anchor_peak_value;

% Compatibility plus full angle statistics.
Sensor_Config.Standard_Relative_Angles = std_angles;
Sensor_Config.Standard_Relative_Angles_Mean = angle_stats.mean;
Sensor_Config.Standard_Relative_Angles_Median = angle_stats.median;
Sensor_Config.Standard_Relative_Angles_Std = angle_stats.std;
Sensor_Config.Standard_Relative_Angles_IQR = angle_stats.iqr;
Sensor_Config.Standard_Relative_Angles_Count = angle_stats.count;
Sensor_Config.Standard_Relative_Angles_Min = angle_stats.min;
Sensor_Config.Standard_Relative_Angles_Max = angle_stats.max;
Sensor_Config.CreatedBy = mfilename;
Sensor_Config.CreatedOn = datestr(now, 31);

summary_table = struct2table(summary_rows);
channel_stats = build_channel_stats_table_local(cfg, case_data, summary_table);
disp(summary_table);

save(fullfile(cfg.step01_output_dir, 'Sensor_Config_20250527.mat'), 'Sensor_Config');
save(fullfile(cfg.step01_output_dir, 'LowSpeed_Features_20250527.mat'), 'case_data', '-v7.3');
writetable(summary_table, fullfile(cfg.step01_output_dir, 'Sensor_Config_Summary_20250527.csv'));
writetable(channel_stats, fullfile(cfg.step01_output_dir, 'LowSpeed_Channel_Stats_20250527.csv'));
writetable(angle_quality_table, fullfile(cfg.step01_output_dir, 'Standard_Relative_Angles_Quality_20250527.csv'));

fprintf('Saved Step01 low-speed reference to:\n  %s\n', ...
    fullfile(cfg.step01_output_dir, 'Sensor_Config_20250527.mat'));
fprintf('OPR timing method: %s\n', Sensor_Config.OPR_Timing_Method);
fprintf('Probe arrival method: %s\n', Sensor_Config.Probe_Arrival_Method);
fprintf('Blade-1 definition: %s\n', Sensor_Config.Blade_ID_Definition);

if show_plots || cfg.save_figures
    plot_step01_reference_local(case_data, Sensor_Config, cfg, summary_table, channel_stats, angle_quality_table, show_plots);
end

function cfg = fill_missing_step01_defaults_local(cfg)
if ~isfield(cfg, 'low_speed_blade1_rule') || isempty(cfg.low_speed_blade1_rule)
    cfg.low_speed_blade1_rule = 'reference_sensor_max_peak_in_initial_revs';
end
if ~isfield(cfg, 'blade_id_offset') || isempty(cfg.blade_id_offset)
    cfg.blade_id_offset = 0;
end
if ~isfield(cfg, 'fingerprint_min_corr') || isempty(cfg.fingerprint_min_corr)
    cfg.fingerprint_min_corr = 0.85;
end
if ~isfield(cfg, 'fingerprint_min_corr_gap') || isempty(cfg.fingerprint_min_corr_gap)
    cfg.fingerprint_min_corr_gap = 0.05;
end
if ~isfield(cfg, 'fingerprint_quality_policy') || isempty(cfg.fingerprint_quality_policy)
    cfg.fingerprint_quality_policy = 'warn';
end
if ~isfield(cfg, 'standard_angle_value') || isempty(cfg.standard_angle_value)
    cfg.standard_angle_value = 'mean';
end
end

function [anchor_start_idx, max_idx, max_val, definition] = choose_reference_blade1_start_local(anchor_feat, cfg)
limit = min(numel(anchor_feat), cfg.blades_num * cfg.reference_search_revs);
search_feat = anchor_feat(1:limit);
[max_val, max_idx] = max(search_feat);

switch lower(strtrim(cfg.low_speed_blade1_rule))
    case 'reference_sensor_max_peak_in_initial_revs'
        % Valid relative numbering convention: the largest reference-sensor
        % pulse in the initial low-speed search region is defined as Blade 1.
        anchor_start_idx = max_idx;
    otherwise
        error('Unsupported cfg.low_speed_blade1_rule: %s', cfg.low_speed_blade1_rule);
end

% Leave enough points for a full six-blade fingerprint. The previous version
% shifted early maxima by one revolution; here the max peak itself is allowed
% to be Blade 1, then clipped only if too close to the file boundary.
anchor_start_idx = min(max(anchor_start_idx, 1), numel(anchor_feat) - cfg.blades_num + 1);

if cfg.blade_id_offset ~= 0
    definition = sprintf(['Internal Blade 1 starts at reference sensor CH%d pulse %d, ' ...
        'selected as the maximum pulse feature in the first %d revolution(s); ' ...
        'reported physical blade IDs use offset %+d.'], ...
        cfg.reference_sensor_id, anchor_start_idx, cfg.reference_search_revs, cfg.blade_id_offset);
else
    definition = sprintf(['Internal Blade 1 starts at reference sensor CH%d pulse %d, ' ...
        'selected as the maximum pulse feature in the first %d revolution(s). ' ...
        'This is a relative blade-numbering convention, not an independently verified physical blade label.'], ...
        cfg.reference_sensor_id, anchor_start_idx, cfg.reference_search_revs);
end
end

function channel_stats = build_channel_stats_table_local(cfg, case_data, summary_table)
all_ids = [cfg.sensor_ids(:).' cfg.opr_id];
rows = repmat(struct( ...
    'ChannelID', NaN, ...
    'Role', '', ...
    'PulseCount', NaN, ...
    'FirstTimeS', NaN, ...
    'LastTimeS', NaN, ...
    'MeanFeature', NaN, ...
    'StdFeature', NaN, ...
    'FingerprintStartIndex', NaN, ...
    'FingerprintCorr', NaN, ...
    'FingerprintSecondCorr', NaN, ...
    'FingerprintCorrGap', NaN, ...
    'FingerprintQuality', ''), numel(all_ids), 1);

for i = 1:numel(all_ids)
    ch = all_ids(i);
    rows(i).ChannelID = ch;
    if ch == cfg.opr_id
        rows(i).Role = 'OPR';
        times = case_data.opr_times(:);
        features = [];
    else
        rows(i).Role = 'sensor';
        times = case_data.channels(ch).arrival_times(:);
        features = case_data.channels(ch).peak_features(:);
        j = find(summary_table.sensor_id == ch, 1);
        if ~isempty(j)
            rows(i).FingerprintStartIndex = summary_table.start_index(j);
            rows(i).FingerprintCorr = summary_table.best_corr(j);
            rows(i).FingerprintSecondCorr = summary_table.second_best_corr(j);
            rows(i).FingerprintCorrGap = summary_table.corr_gap(j);
            rows(i).FingerprintQuality = char(summary_table.quality_status(j));
        end
    end

    rows(i).PulseCount = numel(times);
    if ~isempty(times)
        rows(i).FirstTimeS = min(times);
        rows(i).LastTimeS = max(times);
    end
    if ~isempty(features)
        rows(i).MeanFeature = mean(features, 'omitnan');
        rows(i).StdFeature = std(features, 'omitnan');
    end
end

channel_stats = struct2table(rows);
end

function [std_angles, angle_stats, angle_quality_table] = build_standard_relative_angle_stats_local(cfg, case_data, target_indices)
max_id = max([cfg.candidate_sensor_ids, cfg.opr_id]);
angle_stats = struct();
angle_stats.mean = nan(max_id, cfg.blades_num);
angle_stats.median = nan(max_id, cfg.blades_num);
angle_stats.std = nan(max_id, cfg.blades_num);
angle_stats.iqr = nan(max_id, cfg.blades_num);
angle_stats.count = zeros(max_id, cfg.blades_num);
angle_stats.min = nan(max_id, cfg.blades_num);
angle_stats.max = nan(max_id, cfg.blades_num);

if numel(case_data.opr_times) <= cfg.blades_num
    std_angles = nan(max_id, cfg.blades_num);
    angle_quality_table = table();
    return;
end

spd_t = case_data.opr_times(1:(end - cfg.blades_num));
spd_v = 360 ./ max(case_data.opr_times((cfg.blades_num + 1):end) - ...
    case_data.opr_times(1:(end - cfg.blades_num)), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');

rows = repmat(struct( ...
    'SensorID', NaN, ...
    'BladeID', NaN, ...
    'AngleMeanDeg', NaN, ...
    'AngleMedianDeg', NaN, ...
    'AngleStdDeg', NaN, ...
    'AngleIQRDeg', NaN, ...
    'AngleMinDeg', NaN, ...
    'AngleMaxDeg', NaN, ...
    'Count', NaN), numel(cfg.sensor_ids) * cfg.blades_num, 1);
row_id = 0;

for sid = cfg.sensor_ids
    meas_times = case_data.channels(sid).arrival_times(:);
    start_idx = target_indices(sid);

    for blade_id = 1:cfg.blades_num
        rel_angles = [];
        for lap = 0:(cfg.max_laps_process - 1)
            pulse_idx = start_idx + lap * cfg.blades_num + (blade_id - 1);
            if pulse_idx > numel(meas_times)
                break;
            end

            t_meas = meas_times(pulse_idx);
            idx_prev_opr = find(case_data.opr_times < t_meas, 1, 'last');
            if isempty(idx_prev_opr)
                continue;
            end

            t_ref = case_data.opr_times(idx_prev_opr);
            t_grid = linspace(t_ref, t_meas, 10);
            rel_angles(end + 1) = trapz(t_grid, F_omega_deg(t_grid)); %#ok<AGROW>
        end

        row_id = row_id + 1;
        rows(row_id).SensorID = sid;
        rows(row_id).BladeID = blade_id;
        if ~isempty(rel_angles)
            rel_angles = rel_angles(:);
            angle_stats.mean(sid, blade_id) = mean(rel_angles, 'omitnan');
            angle_stats.median(sid, blade_id) = median(rel_angles, 'omitnan');
            angle_stats.std(sid, blade_id) = std(rel_angles, 'omitnan');
            angle_stats.iqr(sid, blade_id) = iqr_local(rel_angles);
            angle_stats.count(sid, blade_id) = sum(isfinite(rel_angles));
            angle_stats.min(sid, blade_id) = min(rel_angles, [], 'omitnan');
            angle_stats.max(sid, blade_id) = max(rel_angles, [], 'omitnan');

            rows(row_id).AngleMeanDeg = angle_stats.mean(sid, blade_id);
            rows(row_id).AngleMedianDeg = angle_stats.median(sid, blade_id);
            rows(row_id).AngleStdDeg = angle_stats.std(sid, blade_id);
            rows(row_id).AngleIQRDeg = angle_stats.iqr(sid, blade_id);
            rows(row_id).AngleMinDeg = angle_stats.min(sid, blade_id);
            rows(row_id).AngleMaxDeg = angle_stats.max(sid, blade_id);
            rows(row_id).Count = angle_stats.count(sid, blade_id);
        end
    end
end

switch lower(strtrim(cfg.standard_angle_value))
    case 'mean'
        std_angles = angle_stats.mean;
    case 'median'
        std_angles = angle_stats.median;
    otherwise
        error('cfg.standard_angle_value must be mean or median. Current: %s', cfg.standard_angle_value);
end
angle_quality_table = struct2table(rows(1:row_id));
end

function val = iqr_local(x)
x = sort(x(isfinite(x)));
if isempty(x)
    val = NaN;
    return;
end
try
    val = iqr(x);
catch
    val = prctile(x, 75) - prctile(x, 25);
end
end

function [best_idx, best_corr, second_best_corr, second_best_idx, corr_gap] = match_reference_fingerprint_local(feats, ref_fp, blades_num, max_starts)
max_start = min(max_starts, numel(feats) - blades_num + 1);
score = repmat(struct('idx', NaN, 'corr', -inf), max_start, 1);
for idx = 1:max_start
    seg = feats(idx:(idx + blades_num - 1));
    score(idx).idx = idx;
    score(idx).corr = safe_corr_local(ref_fp(:), seg(:));
end
corr_values = [score.corr];
[sorted_corr, order] = sort(corr_values, 'descend');
best_idx = score(order(1)).idx;
best_corr = sorted_corr(1);
if numel(sorted_corr) >= 2
    second_best_corr = sorted_corr(2);
    second_best_idx = score(order(2)).idx;
else
    second_best_corr = NaN;
    second_best_idx = NaN;
end
if ~isfinite(best_corr)
    best_corr = NaN;
end
if ~isfinite(second_best_corr)
    second_best_corr = NaN;
end
corr_gap = best_corr - second_best_corr;
if ~isfinite(corr_gap)
    corr_gap = NaN;
end
end

function [quality_status, msg] = evaluate_fingerprint_quality_local(sid, best_corr, corr_gap, cfg)
quality_status = 'pass';
msg = '';
if ~isfinite(best_corr) || best_corr < cfg.fingerprint_min_corr
    quality_status = 'low_corr';
    msg = sprintf('CH%d fingerprint correlation %.3g is below threshold %.3g.', ...
        sid, best_corr, cfg.fingerprint_min_corr);
elseif ~isfinite(corr_gap) || corr_gap < cfg.fingerprint_min_corr_gap
    quality_status = 'ambiguous';
    msg = sprintf('CH%d fingerprint corr gap %.3g is below threshold %.3g.', ...
        sid, corr_gap, cfg.fingerprint_min_corr_gap);
end
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

function plot_step01_reference_local(case_data, Sensor_Config, cfg, summary_table, channel_stats, angle_quality_table, show_plots)
sensor_ids = Sensor_Config.Sensor_IDs;
blades_num = Sensor_Config.Blades_Num;
n_show = min(4 * blades_num, min_pulse_count_local(case_data, sensor_ids));

fig0 = figure('Name', 'Step01 low-speed raw waveform overview', ...
    'Color', 'w', 'Position', [60, 60, 1450, 850], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig0, numel(sensor_ids) + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
plot_raw_preview_channel_local(case_data, cfg.opr_id, sprintf('CH%d OPR', cfg.opr_id), cfg.opr_threshold);
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    plot_raw_preview_channel_local(case_data, sid, sprintf('CH%d sensor', sid), cfg.sensor_thresholds(sid));
end
save_step01_figure_local(fig0, cfg, 'Step01_LowSpeed_RawWaveform_Overview');

fig1 = figure('Name', 'Step01 low-speed fingerprint alignment', ...
    'Color', 'w', 'Position', [80, 80, 1400, 820], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig1, numel(sensor_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
anchor_sid = Sensor_Config.Reference_Sensor_ID;
ref_fp = Sensor_Config.Fingerprints(anchor_sid);
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    feats = case_data.channels(sid).peak_features(:);
    start_idx = Sensor_Config.Target_Indices(sid);
    pulse_idx = 1:n_show;

    nexttile;
    hold on; grid on; box on;
    plot(pulse_idx, feats(pulse_idx), '-', 'Color', [0.70 0.70 0.70], 'LineWidth', 1.0);
    sel_idx = start_idx:(start_idx + blades_num - 1);
    plot(sel_idx, feats(sel_idx), 'o-', 'Color', [0 0.45 0.74], ...
        'LineWidth', 1.5, 'MarkerFaceColor', [0 0.45 0.74]);
    if sid ~= anchor_sid
        plot(sel_idx, ref_fp(:), 's--', 'Color', [0.85 0.33 0.10], ...
            'LineWidth', 1.2, 'MarkerFaceColor', [0.85 0.33 0.10]);
    end
    xlabel('Pulse index');
    ylabel('Pulse feature');
    title(sprintf('CH%d fingerprint alignment', sid));
end
save_step01_figure_local(fig1, cfg, 'Step01_Fingerprint_Alignment');

fig2 = figure('Name', 'Step01 standard relative angles', ...
    'Color', 'w', 'Position', [120, 120, 1200, 520], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
angle_mat = Sensor_Config.Standard_Relative_Angles(sensor_ids, :);
imagesc(angle_mat);
axis tight; grid on; colorbar;
set(gca, 'XTick', 1:blades_num, 'YTick', 1:numel(sensor_ids), ...
    'YTickLabel', compose('CH%d', sensor_ids));
xlabel('Blade ID');
ylabel('Sensor');
title(sprintf('Standard relative angles (%s, deg)', Sensor_Config.Standard_Relative_Angles_SelectedStatistic));

nexttile;
hold on; grid on; box on;
colors = lines(numel(sensor_ids));
for i = 1:numel(sensor_ids)
    plot(1:blades_num, angle_mat(i, :), 'o-', 'LineWidth', 1.5, ...
        'Color', colors(i, :), 'MarkerFaceColor', colors(i, :), ...
        'DisplayName', sprintf('CH%d', sensor_ids(i)));
end
xlabel('Blade ID');
ylabel('Angle relative to previous OPR (deg)');
title('Per-sensor blade angle curves');
legend('Location', 'best');
save_step01_figure_local(fig2, cfg, 'Step01_Standard_Relative_Angles');

fig3 = figure('Name', 'Step01 low-speed RPM and summary', ...
    'Color', 'w', 'Position', [160, 160, 1200, 620], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
tiledlayout(fig3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(case_data.omega_time_s, case_data.omega_rpm, 'k-', 'LineWidth', 1.2);
grid on; box on;
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('Low-speed RPM: %s', cfg.low_speed_case), 'Interpreter', 'none');
nexttile;
uitable('Data', table2cell(summary_table), ...
    'ColumnName', summary_table.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]);
save_step01_figure_local(fig3, cfg, 'Step01_RPM_And_Fingerprint_Summary');

fig4 = figure('Name', 'Step01 low-speed channel statistics', ...
    'Color', 'w', 'Position', [180, 180, 1100, 300], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
uitable('Data', table2cell(channel_stats), ...
    'ColumnName', channel_stats.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]);
save_step01_figure_local(fig4, cfg, 'Step01_Channel_Statistics');

fig5 = figure('Name', 'Step01 standard angle quality table', ...
    'Color', 'w', 'Position', [200, 200, 1300, 360], 'NumberTitle', 'off', ...
    'Visible', visibility_state_local(show_plots));
uitable('Data', table2cell(angle_quality_table), ...
    'ColumnName', angle_quality_table.Properties.VariableNames, ...
    'Units', 'normalized', 'Position', [0 0 1 1]);
save_step01_figure_local(fig5, cfg, 'Step01_Standard_Angle_Quality_Table');

if ~show_plots
    close([fig0 fig1 fig2 fig3 fig4 fig5]);
end
end

function plot_raw_preview_channel_local(case_data, channel_id, label_text, threshold)
nexttile;
hold on; grid on; box on;
if channel_id <= numel(case_data.raw_preview) && ~isempty(case_data.raw_preview(channel_id).t)
    t = case_data.raw_preview(channel_id).t(:);
    v = case_data.raw_preview(channel_id).v(:);
    max_points = 200000;
    if numel(t) > max_points
        idx = round(linspace(1, numel(t), max_points));
        t = t(idx);
        v = v(idx);
    end
    plot(t - t(1), v, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.6);
end
yline(threshold, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
xlabel('Time from file start (s)');
ylabel('Voltage (V)');
title(label_text);
end

function save_step01_figure_local(fig, cfg, tag)
if ~cfg.save_figures
    return;
end
if ~exist(cfg.step01_figure_dir, 'dir')
    mkdir(cfg.step01_figure_dir);
end
png_file = fullfile(cfg.step01_figure_dir, [tag, '.png']);
pdf_file = fullfile(cfg.step01_figure_dir, [tag, '.pdf']);
try
    exportgraphics(fig, png_file, 'Resolution', 300);
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
catch
    saveas(fig, png_file);
    saveas(fig, pdf_file);
end
end

function state = visibility_state_local(show_plots)
if show_plots
    state = 'on';
else
    state = 'off';
end
end

function n = min_pulse_count_local(case_data, sensor_ids)
counts = zeros(numel(sensor_ids), 1);
for i = 1:numel(sensor_ids)
    counts(i) = numel(case_data.channels(sensor_ids(i)).peak_features);
end
n = min(counts);
end
