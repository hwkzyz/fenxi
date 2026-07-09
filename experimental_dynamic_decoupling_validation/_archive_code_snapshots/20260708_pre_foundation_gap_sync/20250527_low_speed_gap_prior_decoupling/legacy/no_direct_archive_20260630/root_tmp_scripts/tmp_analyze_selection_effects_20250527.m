%% Audit data-selection effects for 20250527 B1 S136
% Read-only audit across:
%   1) rotating_calibration direct-template, no gap prior
%   2) gap_prior_decoupling Step07J gap_tilt

clear; clc;

rootDir = pwd;
rotDir = fullfile(rootDir, 'experimental_dynamic_decoupling_validation', ...
    '20250527_low_speed_rotating_calibration');
gapDir = fullfile(rootDir, 'experimental_dynamic_decoupling_validation', ...
    '20250527_low_speed_gap_prior_decoupling');

rotOut = fullfile(rotDir, 'output', 'identification');
gapOut = fullfile(gapDir, 'outputs');
auditDir = fullfile(gapOut, 'selection_effect_audit_20250527');
if exist(auditDir, 'dir') ~= 7
    mkdir(auditDir);
end

files = struct();
files.rot_core = fullfile(rotOut, ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_20250527.mat');
files.rot_phase = fullfile(rotOut, ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20250527.mat');
files.gap_direct = fullfile(gapOut, ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt.mat');
files.gap_nodirect = fullfile(gapOut, ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt_nodirect.mat');

assert_file_local(files.rot_core);
assert_file_local(files.rot_phase);
assert_file_local(files.gap_direct);
assert_file_local(files.gap_nodirect);

TrotCore = extract_rotating_result_local(files.rot_core, ...
    "Step03 direct-template core", "rotating_calibration", "no_gap", "no_phase_safe");
TrotPhase = extract_rotating_result_local(files.rot_phase, ...
    "Step03 direct-template PrevWinPhaseSafe", "rotating_calibration", "no_gap", "phase_safe");
TgapDirect = extract_gap_result_local(files.gap_direct, ...
    "Step07J gap_tilt direct-bundle", "gap_prior_decoupling", "gap_tilt", "direct_bundle");
TgapNo = extract_gap_result_local(files.gap_nodirect, ...
    "Step07J gap_tilt no-direct", "gap_prior_decoupling", "gap_tilt", "self_selected");

Tall = [TrotCore; TrotPhase; TgapDirect; TgapNo];
summary = summarize_case_table_local(Tall);

rotCompare = compare_pair_local(TrotCore, TrotPhase, ...
    "Step03 core -> PrevWinPhaseSafe");
gapCompare = compare_pair_local(TgapNo, TgapDirect, ...
    "Step07J no-direct -> direct-bundle");

writetable(Tall, fullfile(auditDir, 'SelectionEffect_WindowMetrics_20250527_B1_S136.csv'));
writetable(summary, fullfile(auditDir, 'SelectionEffect_Summary_20250527_B1_S136.csv'));
writetable(rotCompare, fullfile(auditDir, 'SelectionEffect_Rotating_PairCompare_20250527_B1_S136.csv'));
writetable(gapCompare, fullfile(auditDir, 'SelectionEffect_GapPrior_PairCompare_20250527_B1_S136.csv'));

fprintf('\n=== Selection effect summary: 20250527 B1 S136 ===\n');
disp(summary);
fprintf('\n=== Rotating calibration pair: core -> PrevWinPhaseSafe ===\n');
disp(pair_summary_local(rotCompare));
fprintf('\n=== Gap prior pair: no-direct -> direct-bundle ===\n');
disp(pair_summary_local(gapCompare));
fprintf('\nSaved audit tables under:\n  %s\n', auditDir);

function assert_file_local(f)
if exist(f, 'file') ~= 2
    error('Missing required file: %s', f);
end
end

function T = extract_rotating_result_local(fileName, label, folderName, modelName, selectionName)
S = load(fileName, 'Result');
R = S.Result;
n = numel(R.WindowResult);
rows = repmat(empty_row_local(), n, 1);
for iw = 1:n
    wr = R.WindowResult(iw);
    rr = wr.Result;
    b = wr.bundle;
    rows(iw).folder = string(folderName);
    rows(iw).case_label = string(label);
    rows(iw).model = string(modelName);
    rows(iw).selection = string(selectionName);
    rows(iw).window_id = iw;
    rows(iw).lap_start = wr.lap_range(1);
    rows(iw).lap_end = wr.lap_range(end);
    rows(iw).point_count = get_num_local(rr, 'point_count', get_num_local(b, 'point_count', NaN));
    rows(iw).core_point_count = get_trend_value_local(R.Trend, iw, 'core_point_count');
    rows(iw).added_point_count = get_trend_value_local(R.Trend, iw, 'added_point_count');
    rows(iw).selection_pass = string(get_str_local(b, 'selection_pass', 'unknown'));
    rows(iw).EO = rr.EO_id;
    rows(iw).frequency_hz = rr.fn_id;
    rows(iw).amplitude_mm = rr.A_id;
    rows(iw).rmse_mV = 1000 * rr.weighted_voltage_rmse;
    rows(iw).plain_rmse_mV = 1000 * rr.plain_voltage_rmse;
    rows(iw).sensor_count = numel(b.sensor_ids);
    rows(iw).mean_points_per_sensor = rows(iw).point_count / max(rows(iw).sensor_count, 1);
    rows(iw).query_guard_mm = get_num_local(b, 'query_guard_mm', NaN);
end
T = struct2table(rows);
T.source_file = repmat(string(fileName), height(T), 1);
end

function T = extract_gap_result_local(fileName, label, folderName, modelName, selectionName)
S = load(fileName, 'Result');
R = S.Result;
n = numel(R.WindowResult);
rows = repmat(empty_row_local(), n, 1);
for iw = 1:n
    wr = R.WindowResult(iw);
    fit = wr.modelFits.gap_tilt;
    b = wr.bundle;
    rows(iw).folder = string(folderName);
    rows(iw).case_label = string(label);
    rows(iw).model = string(modelName);
    rows(iw).selection = string(selectionName);
    rows(iw).window_id = iw;
    rows(iw).lap_start = wr.lapRange(1);
    rows(iw).lap_end = wr.lapRange(end);
    rows(iw).point_count = get_num_local(b, 'pointCount', NaN);
    rows(iw).core_point_count = get_selection_num_local(wr, 'phaseSafeCorePointCount');
    rows(iw).added_point_count = get_selection_num_local(wr, 'phaseSafeSelectedPointCount') - ...
        get_selection_num_local(wr, 'phaseSafeCorePointCount');
    rows(iw).selection_pass = string(get_str_local(b, 'selectionPass', selectionName));
    rows(iw).EO = fit.EO;
    rows(iw).frequency_hz = fit.freqHz;
    rows(iw).amplitude_mm = fit.amplitudeMm;
    rows(iw).rmse_mV = fit.weightedRmseMv;
    rows(iw).plain_rmse_mV = fit.plainRmseMv;
    rows(iw).sensor_count = numel(b.sensorIds);
    rows(iw).mean_points_per_sensor = rows(iw).point_count / max(rows(iw).sensor_count, 1);
    rows(iw).query_guard_mm = get_num_local(b, 'queryGuardMm', NaN);
end
T = struct2table(rows);
T.source_file = repmat(string(fileName), height(T), 1);
end

function row = empty_row_local()
row = struct('folder', "", 'case_label', "", 'model', "", 'selection', "", ...
    'window_id', NaN, 'lap_start', NaN, 'lap_end', NaN, ...
    'point_count', NaN, 'core_point_count', NaN, 'added_point_count', NaN, ...
    'selection_pass', "", 'EO', NaN, 'frequency_hz', NaN, 'amplitude_mm', NaN, ...
    'rmse_mV', NaN, 'plain_rmse_mV', NaN, 'sensor_count', NaN, ...
    'mean_points_per_sensor', NaN, 'query_guard_mm', NaN);
end

function summary = summarize_case_table_local(T)
cases = unique(T.case_label, 'stable');
rows = repmat(struct('case_label', "", 'folder', "", 'model', "", 'selection', "", ...
    'n_windows', NaN, 'dominant_EO', NaN, 'EO_consistency', NaN, ...
    'mean_points', NaN, 'median_points', NaN, 'mean_added_points', NaN, ...
    'mean_frequency_hz', NaN, 'std_frequency_hz', NaN, ...
    'mean_amplitude_mm', NaN, 'std_amplitude_mm', NaN, ...
    'mean_rmse_mV', NaN, 'median_rmse_mV', NaN, 'std_rmse_mV', NaN), numel(cases), 1);
for i = 1:numel(cases)
    idx = T.case_label == cases(i);
    Ti = T(idx, :);
    rows(i).case_label = cases(i);
    rows(i).folder = Ti.folder(1);
    rows(i).model = Ti.model(1);
    rows(i).selection = Ti.selection(1);
    rows(i).n_windows = height(Ti);
    rows(i).dominant_EO = mode(Ti.EO);
    rows(i).EO_consistency = mean(Ti.EO == rows(i).dominant_EO);
    rows(i).mean_points = mean(Ti.point_count, 'omitnan');
    rows(i).median_points = median(Ti.point_count, 'omitnan');
    rows(i).mean_added_points = mean(Ti.added_point_count, 'omitnan');
    rows(i).mean_frequency_hz = mean(Ti.frequency_hz, 'omitnan');
    rows(i).std_frequency_hz = std(Ti.frequency_hz, 'omitnan');
    rows(i).mean_amplitude_mm = mean(Ti.amplitude_mm, 'omitnan');
    rows(i).std_amplitude_mm = std(Ti.amplitude_mm, 'omitnan');
    rows(i).mean_rmse_mV = mean(Ti.rmse_mV, 'omitnan');
    rows(i).median_rmse_mV = median(Ti.rmse_mV, 'omitnan');
    rows(i).std_rmse_mV = std(Ti.rmse_mV, 'omitnan');
end
summary = struct2table(rows);
end

function C = compare_pair_local(A, B, label)
C = table();
C.compare_label = repmat(string(label), height(A), 1);
C.window_id = A.window_id;
C.points_A = A.point_count;
C.points_B = B.point_count;
C.point_ratio_B_over_A = B.point_count ./ A.point_count;
C.added_points_B_minus_A = B.point_count - A.point_count;
C.EO_A = A.EO;
C.EO_B = B.EO;
C.frequency_A_hz = A.frequency_hz;
C.frequency_B_hz = B.frequency_hz;
C.frequency_delta_hz = B.frequency_hz - A.frequency_hz;
C.amplitude_A_mm = A.amplitude_mm;
C.amplitude_B_mm = B.amplitude_mm;
C.amplitude_delta_mm = B.amplitude_mm - A.amplitude_mm;
C.rmse_A_mV = A.rmse_mV;
C.rmse_B_mV = B.rmse_mV;
C.rmse_delta_B_minus_A_mV = B.rmse_mV - A.rmse_mV;
C.rmse_ratio_B_over_A = B.rmse_mV ./ A.rmse_mV;
end

function S = pair_summary_local(C)
S = table();
S.compare_label = C.compare_label(1);
S.mean_points_A = mean(C.points_A, 'omitnan');
S.mean_points_B = mean(C.points_B, 'omitnan');
S.mean_point_ratio_B_over_A = mean(C.point_ratio_B_over_A, 'omitnan');
S.mean_rmse_A_mV = mean(C.rmse_A_mV, 'omitnan');
S.mean_rmse_B_mV = mean(C.rmse_B_mV, 'omitnan');
S.mean_rmse_delta_B_minus_A_mV = mean(C.rmse_delta_B_minus_A_mV, 'omitnan');
S.mean_rmse_ratio_B_over_A = mean(C.rmse_ratio_B_over_A, 'omitnan');
S.mean_amplitude_delta_mm = mean(C.amplitude_delta_mm, 'omitnan');
S.EO_changed_windows = sum(C.EO_A ~= C.EO_B);
end

function value = get_trend_value_local(T, row, name)
if istable(T) && ismember(name, T.Properties.VariableNames) && height(T) >= row
    value = T.(name)(row);
else
    value = NaN;
end
end

function value = get_selection_num_local(wr, fieldName)
value = NaN;
if isfield(wr, 'EOSelectionInfo') && isfield(wr.EOSelectionInfo, fieldName)
    value = wr.EOSelectionInfo.(fieldName);
end
end

function value = get_num_local(s, fieldName, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, fieldName)
    raw = s.(fieldName);
    if isnumeric(raw) && ~isempty(raw)
        value = raw(1);
    end
end
end

function value = get_str_local(s, fieldName, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, fieldName)
    raw = s.(fieldName);
    if isstring(raw) || ischar(raw)
        value = char(raw);
    end
end
end
