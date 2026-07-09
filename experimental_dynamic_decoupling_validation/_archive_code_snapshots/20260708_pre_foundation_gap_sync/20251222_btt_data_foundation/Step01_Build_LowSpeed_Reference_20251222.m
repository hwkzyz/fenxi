clc; clear; close all;

%STEP01_BUILD_LOWSPEED_REFERENCE_20251222 Bootstrap low-speed reference.
%
% First-round role:
%   Adapt an existing 20251222 low-speed Sensor_Config artifact into the
%   btt_data_foundation interface without re-running raw data extraction.
%
% This script does not build LowSpeed_Features_20251222.mat yet. It writes a
% transparent metadata file stating which existing artifact was used.

%% Step01 run settings
show_tables = true;
show_figures = true;
force_rebuild = true;
validate_data_paths = true;

cfg = BTTDataConfig_20251222();

out_dir = cfg.step01_output_dir;
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
fig_dir = cfg.step01_figure_dir;
if cfg.save_figures && exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

sensor_config_out = fullfile(out_dir, 'Sensor_Config_20251222.mat');
summary_out = fullfile(out_dir, 'Sensor_Config_Summary_20251222.csv');
angle_quality_out = fullfile(out_dir, 'Standard_Relative_Angles_Quality_20251222.csv');
path_audit_out = fullfile(out_dir, 'Step01_Path_Audit_20251222.csv');
metadata_out = fullfile(out_dir, 'LowSpeed_Reference_Metadata_20251222.mat');

if ~force_rebuild && isfile(sensor_config_out)
    fprintf('Step01 output already exists:\n  %s\n', sensor_config_out);
    return;
end

[source_file, source_label] = resolve_step01_source_local(cfg);
loaded = load(source_file);
if ~isfield(loaded, 'Sensor_Config')
    error('Source file does not contain Sensor_Config:\n  %s', source_file);
end

Sensor_Config = standardize_sensor_config_local(loaded.Sensor_Config, cfg, source_file, source_label);
validate_sensor_config_local(Sensor_Config, cfg, source_file);

SensorConfigSummary = build_sensor_config_summary_local(Sensor_Config, cfg, source_file);
OPRCenterAngleTable = build_opr_center_angle_table_local(Sensor_Config, cfg);

AngleQualityTable = build_angle_quality_table_local(Sensor_Config, cfg, source_file, source_label);
PathAuditTable = build_path_audit_table_local(cfg, source_file, source_label);

Step01Meta = struct();
Step01Meta.Dataset = cfg.dataset;
Step01Meta.CreatedBy = mfilename;
Step01Meta.CreatedOn = datestr(now, 31);
Step01Meta.Mode = 'bootstrap_from_existing_sensor_config';
Step01Meta.SourceFile = source_file;
Step01Meta.SourceLabel = source_label;
Step01Meta.OutputFile = sensor_config_out;
Step01Meta.PathAuditFile = path_audit_out;
Step01Meta.FigureDir = fig_dir;
Step01Meta.Note = ['First-round adapter only. It standardizes existing ', ...
    'low-speed reference metadata; raw low-speed feature extraction will ', ...
    'be added in a later round.'];

LowSpeedReference = build_low_speed_reference_local(Sensor_Config, cfg, source_file, source_label);

save(sensor_config_out, 'Sensor_Config', 'SensorConfigSummary', ...
    'OPRCenterAngleTable', 'AngleQualityTable', 'PathAuditTable', 'Step01Meta');
save(metadata_out, 'LowSpeedReference', 'Step01Meta');
writetable(SensorConfigSummary, summary_out);
writetable(AngleQualityTable, angle_quality_out);
writetable(PathAuditTable, path_audit_out);

if show_figures || cfg.save_figures
    plot_step01_diagnostics_local(Sensor_Config, cfg, fig_dir, show_figures, cfg.save_figures);
end

fprintf('\n=== Step01: low-speed reference bootstrap ===\n');
fprintf('Source [%s]:\n  %s\n', source_label, source_file);
fprintf('Saved Sensor_Config:\n  %s\n', sensor_config_out);
fprintf('Saved summary:\n  %s\n', summary_out);
fprintf('Saved angle quality table:\n  %s\n', angle_quality_out);
fprintf('Saved path audit:\n  %s\n', path_audit_out);
fprintf('Saved metadata:\n  %s\n', metadata_out);
if cfg.save_figures
    fprintf('Saved figures under:\n  %s\n', fig_dir);
end

if show_tables
    disp(SensorConfigSummary);
    disp(AngleQualityTable);
    if validate_data_paths
        disp(PathAuditTable);
    end
end

function [source_file, source_label] = resolve_step01_source_local(cfg)
candidates = { ...
    cfg.source.newflow_sensor_config_file, 'newflow_sensor_config'; ...
    cfg.source.legacy_sensor_config_file, 'legacy_reference_sensor_config'};

for i = 1:size(candidates, 1)
    if isfile(candidates{i, 1})
        source_file = candidates{i, 1};
        source_label = candidates{i, 2};
        return;
    end
end

msg = sprintf('No usable Step01 bootstrap source found. Checked:\n');
for i = 1:size(candidates, 1)
    msg = sprintf('%s  - %s\n', msg, candidates{i, 1}); %#ok<AGROW>
end
error('%s', msg);
end

function Sensor_Config = standardize_sensor_config_local(Sensor_Config, cfg, source_file, source_label)
Sensor_Config.Dataset = get_field_or_default_local(Sensor_Config, 'Dataset', cfg.dataset);
Sensor_Config.ReferenceCase = get_field_or_default_local(Sensor_Config, 'ReferenceCase', cfg.low_speed_case);
Sensor_Config.Sensor_IDs = get_field_or_default_local(Sensor_Config, 'Sensor_IDs', cfg.sensor_ids);
Sensor_Config.Candidate_Sensor_IDs = get_field_or_default_local( ...
    Sensor_Config, 'Candidate_Sensor_IDs', cfg.candidate_sensor_ids);
Sensor_Config.Capacitance_IDs = get_field_or_default_local(Sensor_Config, 'Capacitance_IDs', cfg.capacitance_ids);
Sensor_Config.Eddy_Current_IDs = get_field_or_default_local(Sensor_Config, 'Eddy_Current_IDs', cfg.eddy_current_ids);
Sensor_Config.Unused_IDs = get_field_or_default_local(Sensor_Config, 'Unused_IDs', cfg.unused_ids);
Sensor_Config.OPR_ID = get_field_or_default_local(Sensor_Config, 'OPR_ID', cfg.opr_id);
Sensor_Config.Blades_Num = get_field_or_default_local(Sensor_Config, 'Blades_Num', cfg.blades_num);

if isfield(Sensor_Config, 'Sample_Rate_Hz') && ~isempty(Sensor_Config.Sample_Rate_Hz)
    sample_rate_hz = Sensor_Config.Sample_Rate_Hz;
elseif isfield(Sensor_Config, 'Pinlv') && ~isempty(Sensor_Config.Pinlv)
    sample_rate_hz = Sensor_Config.Pinlv;
else
    sample_rate_hz = cfg.sample_rate_hz;
end
Sensor_Config.Sample_Rate_Hz = sample_rate_hz;
Sensor_Config.Pinlv = sample_rate_hz;

Sensor_Config.R_Tip_mm = get_field_or_default_local(Sensor_Config, 'R_Tip_mm', cfg.r_tip_mm);
Sensor_Config.Gap_Points = get_field_or_default_local(Sensor_Config, 'Gap_Points', cfg.gap_points);
Sensor_Config.OPR_Threshold = get_field_or_default_local(Sensor_Config, 'OPR_Threshold', cfg.opr_threshold);
Sensor_Config.Sensor_Thresholds = get_field_or_default_local(Sensor_Config, 'Sensor_Thresholds', cfg.sensor_thresholds);
Sensor_Config.Reference_Sensor_ID = get_field_or_default_local( ...
    Sensor_Config, 'Reference_Sensor_ID', cfg.reference_sensor_id);

if ~isfield(Sensor_Config, 'Standard_Relative_Angles') && ...
        isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    Sensor_Config.Standard_Relative_Angles = Sensor_Config.Standard_Relative_Angles_OPRCenter;
end
if ~isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter') && ...
        isfield(Sensor_Config, 'Standard_Relative_Angles')
    Sensor_Config.Standard_Relative_Angles_OPRCenter = Sensor_Config.Standard_Relative_Angles;
end
if isfield(cfg, 'standard_angle_reference_field') && ...
        isfield(Sensor_Config, cfg.standard_angle_reference_field)
    Sensor_Config.Standard_Relative_Angles = ...
        Sensor_Config.(cfg.standard_angle_reference_field);
end

Sensor_Config.Standard_Relative_Angles_Reference = get_field_or_default_local( ...
    cfg, 'standard_angle_reference_name', ...
    get_field_or_default_local(Sensor_Config, ...
    'Standard_Relative_Angles_Reference', 'opr_pulse_center'));
Sensor_Config.Standard_Relative_Angles_SelectedStatistic = get_field_or_default_local( ...
    Sensor_Config, 'Standard_Relative_Angles_SelectedStatistic', 'bootstrap_existing');
Sensor_Config.OPR_Timing_Method = get_field_or_default_local( ...
    cfg, 'opr_timing_method', ...
    get_field_or_default_local(Sensor_Config, ...
    'OPR_Timing_Method', 'existing_opr_center_reference'));
Sensor_Config.Probe_Arrival_Method = get_field_or_default_local( ...
    Sensor_Config, 'Probe_Arrival_Method', 'existing_sensor_config');
Sensor_Config.Time_Index_Mode = get_field_or_default_local( ...
    Sensor_Config, 'Time_Index_Mode', 'unknown_bootstrap');

Sensor_Config.Bootstrap_Source_File = source_file;
Sensor_Config.Bootstrap_Source_Label = source_label;
end

function validate_sensor_config_local(Sensor_Config, cfg, source_file)
required_fields = {'Fingerprints', 'Target_Indices', 'Standard_Relative_Angles_OPRCenter'};
for i = 1:numel(required_fields)
    name = required_fields{i};
    if ~isfield(Sensor_Config, name)
        error('Sensor_Config missing required field %s:\n  %s', name, source_file);
    end
end
if ~isa(Sensor_Config.Fingerprints, 'containers.Map')
    error('Sensor_Config.Fingerprints must be containers.Map:\n  %s', source_file);
end
if ~isa(Sensor_Config.Target_Indices, 'containers.Map')
    error('Sensor_Config.Target_Indices must be containers.Map:\n  %s', source_file);
end
if size(Sensor_Config.Standard_Relative_Angles_OPRCenter, 2) < cfg.blades_num
    error('Standard_Relative_Angles_OPRCenter has fewer than %d blade columns.', cfg.blades_num);
end
for sid = cfg.sensor_ids(:).'
    if ~isKey(Sensor_Config.Fingerprints, sid)
        error('Sensor_Config.Fingerprints does not contain CH%d.', sid);
    end
    fp = Sensor_Config.Fingerprints(sid);
    if numel(fp) < cfg.blades_num
        error('CH%d fingerprint has %d values; expected at least %d.', ...
            sid, numel(fp), cfg.blades_num);
    end
end
end

function T = build_sensor_config_summary_local(Sensor_Config, cfg, source_file)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'HasTargetIndex', false, ...
    'TargetIndexB1', NaN, ...
    'AngleReference', "", ...
    'FingerprintB1', NaN, ...
    'FingerprintB2', NaN, ...
    'FingerprintB3', NaN, ...
    'FingerprintB4', NaN, ...
    'FingerprintB5', NaN, ...
    'FingerprintB6', NaN, ...
    'AngleB1Deg', NaN, ...
    'AngleB2Deg', NaN, ...
    'AngleB3Deg', NaN, ...
    'AngleB4Deg', NaN, ...
    'AngleB5Deg', NaN, ...
    'AngleB6Deg', NaN, ...
    'SourceFile', ""), numel(cfg.sensor_ids), 1);

for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    fp = force_row_vector_local(Sensor_Config.Fingerprints(sid));
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles(sid, 1:cfg.blades_num));
    rows(i).SensorID = sid;
    rows(i).AngleReference = string(Sensor_Config.Standard_Relative_Angles_Reference);
    rows(i).SourceFile = string(source_file);
    if isKey(Sensor_Config.Target_Indices, sid)
        rows(i).HasTargetIndex = true;
        rows(i).TargetIndexB1 = Sensor_Config.Target_Indices(sid);
    end
    for blade_id = 1:cfg.blades_num
        rows(i).(sprintf('FingerprintB%d', blade_id)) = fp(blade_id);
        rows(i).(sprintf('AngleB%dDeg', blade_id)) = angles(blade_id);
    end
end
T = struct2table(rows);
end

function T = build_opr_center_angle_table_local(Sensor_Config, cfg)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'BladeID', NaN, ...
    'OPRCenterAngleDeg', NaN, ...
    'FingerprintPeak', NaN), numel(cfg.sensor_ids) * cfg.blades_num, 1);

row_idx = 0;
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    fp = force_row_vector_local(Sensor_Config.Fingerprints(sid));
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles(sid, 1:cfg.blades_num));
    for blade_id = 1:cfg.blades_num
        row_idx = row_idx + 1;
        rows(row_idx).SensorID = sid;
        rows(row_idx).BladeID = blade_id;
        rows(row_idx).OPRCenterAngleDeg = angles(blade_id);
        rows(row_idx).FingerprintPeak = fp(blade_id);
    end
end
T = struct2table(rows);
end

function T = build_angle_quality_table_local(Sensor_Config, cfg, source_file, source_label)
rows = repmat(struct( ...
    'SensorID', NaN, ...
    'BladeID', NaN, ...
    'StandardAngleDeg', NaN, ...
    'AngleSource', "", ...
    'QualityStatus', "", ...
    'SourceLabel', "", ...
    'SourceFile', ""), numel(cfg.sensor_ids) * cfg.blades_num, 1);

row_idx = 0;
for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles(sid, 1:cfg.blades_num));
    for blade_id = 1:cfg.blades_num
        row_idx = row_idx + 1;
        rows(row_idx).SensorID = sid;
        rows(row_idx).BladeID = blade_id;
        rows(row_idx).StandardAngleDeg = angles(blade_id);
        rows(row_idx).AngleSource = string(Sensor_Config.Standard_Relative_Angles_Reference);
        rows(row_idx).QualityStatus = "bootstrap_from_existing_sensor_config";
        rows(row_idx).SourceLabel = string(source_label);
        rows(row_idx).SourceFile = string(source_file);
    end
end
T = struct2table(rows);
end

function T = build_path_audit_table_local(cfg, source_file, source_label)
dynamic_case = '';
if ~isempty(cfg.dynamic_cases)
    dynamic_case = char(cfg.dynamic_cases{1});
end

items = { ...
    'dataset_root', cfg.dataset_root, 'dir', 'raw sensor data root'; ...
    'low_speed_case_dir', fullfile(cfg.dataset_root, cfg.low_speed_case), 'dir', 'low-speed calibration case'; ...
    'dynamic_case_dir', fullfile(cfg.dataset_root, dynamic_case), 'dir', 'default dynamic case'; ...
    'strain_root', cfg.strain_root, 'dir', 'optional strain data root'; ...
    'newflow_sensor_config', cfg.source.newflow_sensor_config_file, 'file', 'preferred Step01 bootstrap source'; ...
    'legacy_sensor_config', cfg.source.legacy_sensor_config_file, 'file', 'fallback Step01 bootstrap source'; ...
    'selected_source', source_file, 'file', ['selected bootstrap source: ', source_label]};

rows = repmat(struct( ...
    'Name', "", ...
    'Path', "", ...
    'ExpectedType', "", ...
    'Exists', false, ...
    'Role', ""), size(items, 1), 1);

for i = 1:size(items, 1)
    expected_type = items{i, 3};
    path_text = items{i, 2};
    if strcmp(expected_type, 'dir')
        exists_tf = exist(path_text, 'dir') == 7;
    else
        exists_tf = exist(path_text, 'file') == 2;
    end
    rows(i).Name = string(items{i, 1});
    rows(i).Path = string(path_text);
    rows(i).ExpectedType = string(expected_type);
    rows(i).Exists = exists_tf;
    rows(i).Role = string(items{i, 4});
end

T = struct2table(rows);
end

function plot_step01_diagnostics_local(Sensor_Config, cfg, fig_dir, show_figures, save_figures)
[fingerprint_mat, angle_mat, target_idx] = collect_sensor_blade_matrices_local(Sensor_Config, cfg);
sensor_labels = compose('CH%d', cfg.sensor_ids);
blade_labels = compose('B%d', 1:cfg.blades_num);

visibility = 'off';
if show_figures
    visibility = 'on';
end

fig = figure('Name', 'Step01 Fingerprint Heatmap', 'Color', 'w', 'Visible', visibility);
tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
imagesc(ax, fingerprint_mat);
axis(ax, 'tight');
set(ax, 'XTick', 1:cfg.blades_num, 'XTickLabel', blade_labels, ...
    'YTick', 1:numel(cfg.sensor_ids), 'YTickLabel', sensor_labels);
xlabel(ax, 'Blade ID');
ylabel(ax, 'Sensor ID');
title(ax, 'Low-Speed Fingerprint Peaks');
cb = colorbar(ax);
cb.Label.String = 'Peak voltage feature';
annotate_matrix_local(ax, fingerprint_mat, '%.3g');
save_figure_local(fig, fig_dir, 'Step01_LowSpeed_Fingerprint_Heatmap_20251222', save_figures);
if ~show_figures
    close(fig);
end

fig = figure('Name', 'Step01 Standard Angle Heatmap', 'Color', 'w', 'Visible', visibility);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
imagesc(ax, angle_mat);
axis(ax, 'tight');
set(ax, 'XTick', 1:cfg.blades_num, 'XTickLabel', blade_labels, ...
    'YTick', 1:numel(cfg.sensor_ids), 'YTickLabel', sensor_labels);
xlabel(ax, 'Blade ID');
ylabel(ax, 'Sensor ID');
title(ax, 'OPR-Center Standard Angles');
cb = colorbar(ax);
cb.Label.String = 'Angle (deg)';
annotate_matrix_local(ax, angle_mat, '%.2f');

ax = nexttile;
angle_spread = max(angle_mat, [], 1, 'omitnan') - min(angle_mat, [], 1, 'omitnan');
bar(ax, 1:cfg.blades_num, angle_spread, 0.65);
grid(ax, 'on');
xlim(ax, [0.5 cfg.blades_num + 0.5]);
set(ax, 'XTick', 1:cfg.blades_num, 'XTickLabel', blade_labels);
xlabel(ax, 'Blade ID');
ylabel(ax, 'Across-sensor range (deg)');
title(ax, 'Standard-Angle Spread');
save_figure_local(fig, fig_dir, 'Step01_Standard_Angle_Quality_20251222', save_figures);
if ~show_figures
    close(fig);
end

fig = figure('Name', 'Step01 Numbering Overview', 'Color', 'w', 'Visible', visibility);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;
plot(ax, 1:cfg.blades_num, normalize_rows_local(fingerprint_mat).', '-o', 'LineWidth', 1.2);
grid(ax, 'on');
xlim(ax, [0.8 cfg.blades_num + 0.2]);
set(ax, 'XTick', 1:cfg.blades_num, 'XTickLabel', blade_labels);
xlabel(ax, 'Blade ID');
ylabel(ax, 'Normalized fingerprint');
title(ax, 'Sensor Fingerprint Profiles');
legend(ax, sensor_labels, 'Location', 'best');

ax = nexttile;
bar(ax, cfg.sensor_ids, target_idx, 0.65);
grid(ax, 'on');
set(ax, 'XTick', cfg.sensor_ids, 'XTickLabel', sensor_labels);
xlabel(ax, 'Sensor ID');
ylabel(ax, 'Target index for B1');
title(ax, 'Blade-1 Target Index by Sensor');
save_figure_local(fig, fig_dir, 'Step01_Numbering_Overview_20251222', save_figures);
if ~show_figures
    close(fig);
end
end

function [fingerprint_mat, angle_mat, target_idx] = collect_sensor_blade_matrices_local(Sensor_Config, cfg)
fingerprint_mat = nan(numel(cfg.sensor_ids), cfg.blades_num);
angle_mat = nan(numel(cfg.sensor_ids), cfg.blades_num);
target_idx = nan(numel(cfg.sensor_ids), 1);

for i = 1:numel(cfg.sensor_ids)
    sid = cfg.sensor_ids(i);
    fp = force_row_vector_local(Sensor_Config.Fingerprints(sid));
    angles = force_row_vector_local(Sensor_Config.Standard_Relative_Angles(sid, 1:cfg.blades_num));
    fingerprint_mat(i, :) = fp(1:cfg.blades_num);
    angle_mat(i, :) = angles(1:cfg.blades_num);
    if isKey(Sensor_Config.Target_Indices, sid)
        target_idx(i) = Sensor_Config.Target_Indices(sid);
    end
end
end

function annotate_matrix_local(ax, M, fmt)
for r = 1:size(M, 1)
    for c = 1:size(M, 2)
        if isfinite(M(r, c))
            text(ax, c, r, sprintf(fmt, M(r, c)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'w', ...
                'FontSize', 8, ...
                'FontWeight', 'bold');
        end
    end
end
end

function X = normalize_rows_local(X)
for i = 1:size(X, 1)
    row = X(i, :);
    lo = min(row, [], 'omitnan');
    hi = max(row, [], 'omitnan');
    if isfinite(lo) && isfinite(hi) && hi > lo
        X(i, :) = (row - lo) ./ (hi - lo);
    end
end
end

function save_figure_local(fig, fig_dir, base_name, save_figures)
if ~save_figures
    return;
end
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end
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

function LowSpeedReference = build_low_speed_reference_local(Sensor_Config, cfg, source_file, source_label)
LowSpeedReference = struct();
LowSpeedReference.Dataset = cfg.dataset;
LowSpeedReference.Mode = 'bootstrap_from_existing_sensor_config';
LowSpeedReference.SourceFile = source_file;
LowSpeedReference.SourceLabel = source_label;
LowSpeedReference.ReferenceCase = Sensor_Config.ReferenceCase;
LowSpeedReference.AnalysisSensors = cfg.sensor_ids;
LowSpeedReference.BladeCount = cfg.blades_num;
LowSpeedReference.Fingerprints = Sensor_Config.Fingerprints;
LowSpeedReference.TargetIndices = Sensor_Config.Target_Indices;
LowSpeedReference.StandardAnglesOPRCenter = Sensor_Config.Standard_Relative_Angles_OPRCenter;
LowSpeedReference.StandardAngleReference = Sensor_Config.Standard_Relative_Angles_Reference;
LowSpeedReference.OPRTimingMethod = Sensor_Config.OPR_Timing_Method;
LowSpeedReference.ProbeArrivalMethod = Sensor_Config.Probe_Arrival_Method;
LowSpeedReference.Note = ['This is a first-round metadata layer. ', ...
    'LowSpeed_Features_20251222.mat will be produced when raw low-speed ', ...
    'feature extraction is migrated.'];
end

function value = get_field_or_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
    value = S.(field_name);
else
    value = default_value;
end
end

function x = force_row_vector_local(x)
x = x(:).';
end
