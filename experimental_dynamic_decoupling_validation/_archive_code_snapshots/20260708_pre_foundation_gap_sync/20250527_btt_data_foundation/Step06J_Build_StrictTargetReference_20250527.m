%% Step06J_Build_StrictTargetReference_20250527.m
% Build a strict target-window reference from Step03 BTT and strain basis.
%
% This script is deliberately narrower than Step06I:
%   1) only the target uniform resonance window is used for the main fit;
%   2) the main model has no free sensor bias, avoiding synchronous EO /
%      sensor-offset confounding;
%   3) free and ridge sensor-bias models are diagnostic only;
%   4) uncertainty is estimated by revolution-block bootstrap.
%
% Main reference:
%   d_i = B + Re{ C_14 h_14(t_i) } + e_i
%   u_ref(t) = Re{ C_14 h_14(t) }

clc; clear; close all;

cfg = BTTDataConfig_20250527();
cfg = apply_step06j_defaults_local(cfg);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06J strict target-window reference: %s ===\n', caseName);

    files = resolve_step06j_files_local(cfg, caseName);
    assert_step06j_files_local(files);

    S06I = load(files.step06i_mat, 'FitPointsEO', 'ReferenceTimeSeries', ...
        'RegionTable', 'files');
    FitPoints = S06I.FitPointsEO;
    TrefAll = S06I.ReferenceTimeSeries;
    ObsS = load(files.step03_bundle_mat, 'bundle');
    Obs = ObsS.bundle.Observation_Table;
    oprTimes = load_opr_times_step06j_local(files.jiluopr_mat);
    EvidenceS = load(files.step06a_evidence_mat, 'Evidence');
    RegionTableAll = get_evidence_region_table_step06j_local(EvidenceS.Evidence);
    LowBias = estimate_low_vibration_sensor_bias_local(Obs, RegionTableAll, cfg);

    [TargetPoints, TrefTarget, targetEO] = select_target_window_local( ...
        FitPoints, TrefAll, Obs, oprTimes, cfg);

    MainFit = fit_target_reference_model_local(TargetPoints, TrefTarget, ...
        cfg, "M0_no_sensor_bias", [], []);
    FixedBiasFit = fit_fixed_sensor_bias_reference_model_local( ...
        TargetPoints, TrefTarget, cfg, LowBias);
    RidgeFit = fit_target_reference_model_local(TargetPoints, TrefTarget, ...
        cfg, "M2_ridge_sensor_bias", MainFit.C_complex_mm_per_microstrain, []);
    FreeFit = fit_target_reference_model_local(TargetPoints, TrefTarget, ...
        cfg, "M3_free_sensor_bias", MainFit.C_complex_mm_per_microstrain, []);
    FixedFit = fit_fixed_k_reference_model_local(TargetPoints, TrefTarget, ...
        cfg, MainFit.C_complex_mm_per_microstrain);

    ModelTable = build_model_table_local( ...
        MainFit, FixedBiasFit, RidgeFit, FreeFit, FixedFit, targetEO, cfg);
    LeaveOneSensor = build_leave_one_sensor_table_local( ...
        TargetPoints, TrefTarget, targetEO, cfg);
    BootstrapTable = bootstrap_target_reference_local( ...
        TargetPoints, TrefTarget, targetEO, cfg);
    StrictTimeSeries = build_strict_time_series_local( ...
        TrefTarget, MainFit, BootstrapTable);
    StrictPoints = attach_strict_point_fit_local(TargetPoints, MainFit);
    StrictSummary = build_strict_summary_local(caseName, targetEO, ...
        TargetPoints, StrictTimeSeries, ModelTable, LeaveOneSensor, ...
        BootstrapTable, cfg);

    outDir = fullfile(cfg.step06j_output_dir, caseName);
    figDir = fullfile(cfg.step06j_figure_dir, caseName);
    ensure_dir_step06j_local(outDir);
    ensure_dir_step06j_local(figDir);

    summaryCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_Summary_20250527.csv');
    modelCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_ModelComparison_20250527.csv');
    pointsCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_FitPoints_20250527.csv');
    tsCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_TimeSeries_20250527.csv');
    losCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_LeaveOneSensor_20250527.csv');
    bootCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_BlockBootstrap_20250527.csv');
    lowBiasCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_LowVibrationSensorBias_20250527.csv');
    lowBiasWindowCsv = fullfile(outDir, ...
        'Step06J_StrictTargetReference_LowVibrationWindows_20250527.csv');
    matFile = fullfile(outDir, ...
        'Step06J_StrictTargetReference_20250527.mat');

    writetable(StrictSummary, summaryCsv);
    writetable(ModelTable, modelCsv);
    writetable(StrictPoints, pointsCsv);
    writetable(StrictTimeSeries, tsCsv);
    writetable(LeaveOneSensor, losCsv);
    writetable(BootstrapTable, bootCsv);
    writetable(LowBias.SensorBiasTable, lowBiasCsv);
    writetable(LowBias.WindowTable, lowBiasWindowCsv);
    save(matFile, 'StrictSummary', 'ModelTable', 'StrictPoints', ...
        'StrictTimeSeries', 'LeaveOneSensor', 'BootstrapTable', ...
        'LowBias', 'MainFit', 'FixedBiasFit', 'RidgeFit', 'FreeFit', ...
        'FixedFit', 'targetEO', ...
        'files', '-v7.3');

    figFile = plot_strict_reference_local(StrictTimeSeries, StrictPoints, ...
        StrictSummary, ModelTable, LeaveOneSensor, BootstrapTable, ...
        figDir, cfg);

    fprintf('\nStep06J strict target reference summary:\n');
    disp(StrictSummary);
    fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        summaryCsv, modelCsv, pointsCsv, tsCsv, losCsv, bootCsv, ...
        lowBiasCsv, lowBiasWindowCsv, matFile, figFile);
end


function cfg = apply_step06j_defaults_local(cfg)
if ~isfield(cfg, 'step06i_target_blade') || isempty(cfg.step06i_target_blade)
    cfg.step06i_target_blade = 1;
end
if ~isfield(cfg, 'step06i_sensor_ids') || isempty(cfg.step06i_sensor_ids)
    cfg.step06i_sensor_ids = cfg.sensor_ids(:).';
end
if ~isfield(cfg, 'step06i_target_time_range_s') || ...
        isempty(cfg.step06i_target_time_range_s)
    cfg.step06i_target_time_range_s = [1.0 5.0];
end
if ~isfield(cfg, 'step06i_fixed_strain_to_mm') || ...
        isempty(cfg.step06i_fixed_strain_to_mm)
    cfg.step06i_fixed_strain_to_mm = 1 / 1.266 / 1000;
end
if ~isfield(cfg, 'step06i_huber_k') || isempty(cfg.step06i_huber_k)
    cfg.step06i_huber_k = 1.345;
end
if ~isfield(cfg, 'step06i_irls_iterations') || ...
        isempty(cfg.step06i_irls_iterations)
    cfg.step06i_irls_iterations = 12;
end
if ~isfield(cfg, 'step06i_random_seed') || isempty(cfg.step06i_random_seed)
    cfg.step06i_random_seed = 20250527;
end
if ~isfield(cfg, 'step06j_bootstrap_count') || ...
        isempty(cfg.step06j_bootstrap_count)
    cfg.step06j_bootstrap_count = 200;
end
if ~isfield(cfg, 'step06j_sensor_bias_prior_mm') || ...
        isempty(cfg.step06j_sensor_bias_prior_mm)
    cfg.step06j_sensor_bias_prior_mm = 0.02;
end
if ~isfield(cfg, 'step06j_weight_sigma_floor_mm') || ...
        isempty(cfg.step06j_weight_sigma_floor_mm)
    cfg.step06j_weight_sigma_floor_mm = 0.008;
end
if ~isfield(cfg, 'step06j_output_dir') || isempty(cfg.step06j_output_dir)
    cfg.step06j_output_dir = fullfile(cfg.output_root, ...
        'step06j_strict_target_reference');
end
if ~isfield(cfg, 'step06j_figure_dir') || isempty(cfg.step06j_figure_dir)
    cfg.step06j_figure_dir = fullfile(cfg.figure_root, ...
        'step06j_strict_target_reference');
end
if ~isfield(cfg, 'step06j_fig_visible') || isempty(cfg.step06j_fig_visible)
    if usejava('desktop')
        cfg.step06j_fig_visible = 'on';
    else
        cfg.step06j_fig_visible = 'off';
    end
end
end


function files = resolve_step06j_files_local(cfg, caseName)
files = struct();
files.step06i_mat = fullfile(cfg.output_root, ...
    'step06i_eo_dependent_joint_reference', caseName, ...
    'Step06I_EODependentReference_20250527.mat');
files.step03_bundle_mat = fullfile(cfg.step03_output_dir, caseName, ...
    'BTT_Observation_Bundle_20250527.mat');
files.jiluopr_mat = fullfile(cfg.step02_output_dir, caseName, 'jiluOPR.mat');
files.step06a_evidence_mat = fullfile(cfg.output_root, ...
    'step06a_strain_rpm_resonance_evidence', caseName, ...
    'Step06A_StrainRPM_ResonanceEvidence_20250527.mat');
end


function assert_step06j_files_local(files)
names = fieldnames(files);
for i = 1:numel(names)
    if ~isfile(files.(names{i}))
        error('Missing Step06J input %s:\n  %s', names{i}, files.(names{i}));
    end
end
end


function R = get_evidence_region_table_step06j_local(Evidence)
R = table();
if isfield(Evidence, 'regionTable') && istable(Evidence.regionTable)
    R = Evidence.regionTable;
elseif isfield(Evidence, 'RegionTable') && istable(Evidence.RegionTable)
    R = Evidence.RegionTable;
end
end


function LowBias = estimate_low_vibration_sensor_bias_local(Obs, R, cfg)
twTarget = cfg.step06i_target_time_range_s(:).';
T = Obs(Obs.is_valid & Obs.blade_id == cfg.step06i_target_blade & ...
    ismember(Obs.sensor_id, cfg.step06i_sensor_ids) & ...
    isfinite(Obs.arrival_time_s) & isfinite(Obs.displacement_mm), :);
if isempty(T)
    error('No valid Step03 observations available for low-vibration bias estimation.');
end

tMin = min(T.arrival_time_s);
tMax = max(T.arrival_time_s);
winLen = 1.0;
winStep = 0.5;
starts = (ceil(tMin / winStep) * winStep):winStep:(tMax - winLen);
rows = repmat(make_low_window_row_local(cfg.step06i_sensor_ids), numel(starts), 1);
rowCount = 0;
for i = 1:numel(starts)
    t0 = starts(i);
    t1 = t0 + winLen;
    if intervals_overlap_local([t0 t1], twTarget)
        continue;
    end
    if overlaps_evidence_region_local([t0 t1], R)
        continue;
    end
    m = T.arrival_time_s >= t0 & T.arrival_time_s < t1;
    Tw = T(m, :);
    if height(Tw) < 30
        continue;
    end
    sids = cfg.step06i_sensor_ids(:).';
    medBySensor = nan(1, numel(sids));
    madBySensor = nan(1, numel(sids));
    cntBySensor = zeros(1, numel(sids));
    for k = 1:numel(sids)
        ms = Tw.sensor_id == sids(k);
        cntBySensor(k) = nnz(ms);
        if cntBySensor(k) > 0
            y = Tw.displacement_mm(ms);
            medBySensor(k) = median(y, 'omitnan');
            madBySensor(k) = robust_scale_step06j_local(y);
        end
    end
    if any(cntBySensor < 6) || any(~isfinite(medBySensor))
        continue;
    end
    delta = medBySensor - median(medBySensor, 'omitnan');
    spread = median(madBySensor, 'omitnan');
    rowCount = rowCount + 1;
    rows(rowCount).window_id = rowCount;
    rows(rowCount).time_start_s = t0;
    rows(rowCount).time_end_s = t1;
    rows(rowCount).point_count = height(Tw);
    rows(rowCount).btt_mad_mm = spread;
    rows(rowCount).rpm_std = std(Tw.rpm, 'omitnan');
    for k = 1:numel(sids)
        rows(rowCount).(sprintf('n_S%d', sids(k))) = cntBySensor(k);
        rows(rowCount).(sprintf('bias_S%d_mm', sids(k))) = delta(k);
    end
end

if rowCount == 0
    warning('No low-vibration windows found. Fixed-bias model will use zero bias.');
    WindowTable = struct2table(rows([]));
    SensorBiasTable = make_zero_sensor_bias_table_local(cfg.step06i_sensor_ids);
else
    WindowTable = struct2table(rows(1:rowCount));
    WindowTable = sortrows(WindowTable, 'btt_mad_mm');
    maxUse = min(height(WindowTable), max(3, ceil(0.30 * height(WindowTable))));
    Used = WindowTable(1:maxUse, :);
    sids = cfg.step06i_sensor_ids(:).';
    brow = repmat(struct('sensor_id', NaN, 'bias_mm', NaN, ...
        'bias_std_across_windows_mm', NaN, 'window_count', NaN, ...
        'status', ""), numel(sids), 1);
    for k = 1:numel(sids)
        col = sprintf('bias_S%d_mm', sids(k));
        vals = Used.(col);
        brow(k).sensor_id = sids(k);
        brow(k).bias_mm = median(vals, 'omitnan');
        brow(k).bias_std_across_windows_mm = std(vals, 'omitnan');
        brow(k).window_count = height(Used);
        if brow(k).bias_std_across_windows_mm <= 0.02
            brow(k).status = "stable";
        else
            brow(k).status = "diagnostic_unstable";
        end
    end
    SensorBiasTable = struct2table(brow);
end

LowBias = struct();
LowBias.WindowTable = WindowTable;
LowBias.SensorBiasTable = SensorBiasTable;
end


function row = make_low_window_row_local(sensorIds)
row = struct('window_id', NaN, 'time_start_s', NaN, 'time_end_s', NaN, ...
    'point_count', NaN, 'btt_mad_mm', NaN, 'rpm_std', NaN);
for sid = sensorIds(:).'
    row.(sprintf('n_S%d', sid)) = NaN;
    row.(sprintf('bias_S%d_mm', sid)) = NaN;
end
end


function T = make_zero_sensor_bias_table_local(sensorIds)
rows = repmat(struct('sensor_id', NaN, 'bias_mm', 0, ...
    'bias_std_across_windows_mm', NaN, 'window_count', 0, ...
    'status', "no_low_vibration_window"), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    rows(i).sensor_id = sensorIds(i);
end
T = struct2table(rows);
end


function tf = intervals_overlap_local(a, b)
tf = a(1) < b(2) && b(1) < a(2);
end


function tf = overlaps_evidence_region_local(win, R)
tf = false;
if isempty(R) || height(R) == 0
    return;
end
startT = get_table_numeric_step06j_local(R, {'bttStartSec', 'time_start_s', 'regionStart'});
endT = get_table_numeric_step06j_local(R, {'bttEndSec', 'time_end_s', 'regionEnd'});
if all(~isfinite(startT)) || all(~isfinite(endT))
    return;
end
tf = any(startT < win(2) & endT > win(1));
end


function x = get_table_numeric_step06j_local(T, names)
x = nan(height(T), 1);
for i = 1:numel(names)
    if ismember(names{i}, T.Properties.VariableNames)
        v = T.(names{i});
        if isnumeric(v)
            x = double(v(:));
        else
            x = str2double(string(v(:)));
        end
        return;
    end
end
end


function [Target, TrefTarget, targetEO] = select_target_window_local( ...
    FitPoints, TrefAll, Obs, oprTimes, cfg)
tw = cfg.step06i_target_time_range_s(:).';
maskTarget = FitPoints.arrival_time_s >= tw(1) & ...
    FitPoints.arrival_time_s <= tw(2) & ...
    FitPoints.blade_id == cfg.step06i_target_blade & ...
    ismember(FitPoints.sensor_id, cfg.step06i_sensor_ids) & ...
    isfinite(FitPoints.displacement_mm) & ...
    isfinite(FitPoints.strain_basis_real_microstrain) & ...
    isfinite(FitPoints.strain_basis_imag_microstrain);
if ~any(maskTarget)
    error('No Step06I fit points found in target window [%.3f, %.3f] s.', tw(1), tw(2));
end
targetEO = mode(FitPoints.dominant_order(maskTarget));
Target = FitPoints(maskTarget & FitPoints.dominant_order == targetEO, :);
Target = sortrows(Target, 'arrival_time_s');
Target.row_id = (1:height(Target)).';

TrefMask = TrefAll.time_s >= tw(1) & TrefAll.time_s <= tw(2) & ...
    TrefAll.dominant_order == targetEO & ...
    isfinite(TrefAll.strain_basis_real_microstrain) & ...
    isfinite(TrefAll.strain_basis_imag_microstrain);
TrefTarget = sortrows(TrefAll(TrefMask, :), 'time_s');
if height(TrefTarget) < 10
    error('Too few target time-series samples for EO%d in [%.3f, %.3f] s.', ...
        targetEO, tw(1), tw(2));
end

Target = attach_quality_from_step03_local(Target, Obs, cfg);
if ~ismember('revolution_index', Target.Properties.VariableNames) || ...
        all(~isfinite(Target.revolution_index))
    Target.revolution_index = compute_revolution_index_local( ...
        oprTimes, Target.arrival_time_s, cfg.blades_num);
end
end


function Target = attach_quality_from_step03_local(Target, Obs, cfg)
qVars = {'revolution_index', 'std_angle_std_deg', 'std_angle_iqr_deg', ...
    'std_angle_count', 'step02_best_corr', 'step02_corr_gap', ...
    'step02_shift_consistency', 'step02_label_quality', 'quality_flag'};
qVars = qVars(ismember(qVars, Obs.Properties.VariableNames));
ObsQ = Obs(Obs.blade_id == cfg.step06i_target_blade & ...
    ismember(Obs.sensor_id, cfg.step06i_sensor_ids), ...
    [{'sensor_id', 'blade_id', 'arrival_time_s'}, qVars]);
Target.join_key = int64(round(Target.arrival_time_s * 1e9));
ObsQ.join_key = int64(round(ObsQ.arrival_time_s * 1e9));
ObsQ.arrival_time_s = [];
Target = innerjoin(Target, ObsQ, ...
    'Keys', {'join_key', 'sensor_id', 'blade_id'});
Target = sortrows(Target, 'row_id');
Target.join_key = [];
end


function Fit = fit_target_reference_model_local(T, Tref, cfg, modelName, mainC, idxRows)
if nargin < 6 || isempty(idxRows)
    idxRows = (1:height(T)).';
end
Ts = T(idxRows, :);
h = Ts.strain_basis_real_microstrain + 1i * Ts.strain_basis_imag_microstrain;
y = Ts.displacement_mm(:);
X = [ones(height(Ts), 1), real(h), -imag(h)];
sensorIds = unique(T.sensor_id(:)).';
sensorBiasCols = [];
if modelName == "M2_ridge_sensor_bias" || modelName == "M3_free_sensor_bias"
    refSensor = min(sensorIds);
    for sid = sensorIds
        if sid == refSensor
            continue;
        end
        X(:, end+1) = double(Ts.sensor_id == sid); %#ok<AGROW>
        sensorBiasCols(end+1) = size(X, 2); %#ok<AGROW>
    end
end

baseWeights = build_quality_weights_local(Ts, cfg);
penalty = zeros(1, size(X, 2));
if modelName == "M2_ridge_sensor_bias" && ~isempty(sensorBiasCols)
    penalty(sensorBiasCols) = 1 / cfg.step06j_sensor_bias_prior_mm;
end
[beta, weights] = weighted_irls_ridge_local(X, y, baseWeights, penalty, cfg);

C = beta(2) + 1i * beta(3);
href = Tref.strain_basis_real_microstrain + 1i * Tref.strain_basis_imag_microstrain;
uRef = real(C .* href);
yFit = X * beta;
resid = y - yFit;

Fit = make_fit_struct_local(modelName, target_eo_from_tref_local(Tref), beta, C, ...
    Tref.time_s, uRef, yFit, resid, weights, baseWeights, Ts, cfg);
Fit.sensor_bias_mm = sensor_bias_from_beta_local(beta, sensorIds, sensorBiasCols);
if nargin >= 5 && ~isempty(mainC)
    Fit.phase_delta_to_main_rad = angle(C / mainC);
else
    Fit.phase_delta_to_main_rad = 0;
end
end


function Fit = fit_fixed_sensor_bias_reference_model_local(T, Tref, cfg, LowBias)
bias = zeros(height(T), 1);
if isfield(LowBias, 'SensorBiasTable') && ~isempty(LowBias.SensorBiasTable)
    BT = LowBias.SensorBiasTable;
    for i = 1:height(BT)
        m = T.sensor_id == BT.sensor_id(i);
        bias(m) = BT.bias_mm(i);
    end
end
Tc = T;
Tc.displacement_mm = Tc.displacement_mm - bias;
Fit = fit_target_reference_model_local(Tc, Tref, cfg, ...
    "M1_fixed_low_vibration_sensor_bias", [], []);
Fit.fixed_sensor_bias_mm_at_points = bias;
Fit.low_vibration_sensor_bias_table = LowBias.SensorBiasTable;
end


function Fit = fit_fixed_k_reference_model_local(T, Tref, cfg, mainC)
h = T.strain_basis_real_microstrain + 1i * T.strain_basis_imag_microstrain;
y = T.displacement_mm(:);
C = cfg.step06i_fixed_strain_to_mm * exp(1i * angle(mainC));
u = real(C .* h);
w = build_quality_weights_local(T, cfg);
B = sum(w .* (y - u)) / sum(w);
yFit = B + u;
resid = y - yFit;
href = Tref.strain_basis_real_microstrain + 1i * Tref.strain_basis_imag_microstrain;
uRef = real(C .* href);
beta = [B; real(C); imag(C)];
Fit = make_fit_struct_local("M4_fixed_K_main_phase", target_eo_from_tref_local(Tref), ...
    beta, C, Tref.time_s, uRef, yFit, resid, w, w, T, cfg);
Fit.sensor_bias_mm = zeros(numel(unique(T.sensor_id)), 1);
Fit.phase_delta_to_main_rad = angle(C / mainC);
end


function Fit = make_fit_struct_local(modelName, targetEO, beta, C, tRef, uRef, ...
    yFit, resid, weights, baseWeights, T, cfg)
uRms = std(uRef, 'omitnan');
Fit = struct();
Fit.model_name = string(modelName);
Fit.target_eo = targetEO;
Fit.beta = beta;
Fit.C_complex_mm_per_microstrain = C;
Fit.K_mm_per_microstrain = abs(C);
Fit.K_um_per_microstrain = 1000 * abs(C);
Fit.K_ratio_to_fixed = abs(C) / cfg.step06i_fixed_strain_to_mm;
Fit.phase_lag_rad = angle(C);
Fit.phase_lag_deg = angle(C) * 180 / pi;
Fit.time_s = tRef;
Fit.u_ref_mm = uRef;
Fit.u_rms_mm = uRms;
Fit.A_peak_equiv_mm = sqrt(2) * uRms;
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - min(uRef, [], 'omitnan'));
Fit.y_fit_mm = yFit;
Fit.residual_mm = resid;
Fit.weights = weights;
Fit.base_weights = baseWeights;
Fit.rmse_mm = sqrt(mean(resid.^2, 'omitnan'));
Fit.weighted_rmse_mm = sqrt(sum(weights .* resid.^2) / sum(weights));
Fit.mean_abs_residual_mm = mean(abs(resid), 'omitnan');
Fit.point_count = height(T);
Fit.sensor_count = numel(unique(T.sensor_id));
Fit.phase_coverage_fraction = phase_coverage_step06j_local(angle( ...
    T.strain_basis_real_microstrain + 1i * T.strain_basis_imag_microstrain)) / (2*pi);
end


function sensorBias = sensor_bias_from_beta_local(beta, sensorIds, sensorBiasCols)
sensorBias = zeros(numel(sensorIds), 1);
if isempty(sensorBiasCols)
    return;
end
refSensor = min(sensorIds);
k = 0;
for i = 1:numel(sensorIds)
    if sensorIds(i) == refSensor
        sensorBias(i) = 0;
    else
        k = k + 1;
        sensorBias(i) = beta(sensorBiasCols(k));
    end
end
end


function w = build_quality_weights_local(T, cfg)
n = height(T);
w = ones(n, 1);
thetaSigmaDeg = nan(n, 1);
if ismember('std_angle_std_deg', T.Properties.VariableNames)
    thetaStd = T.std_angle_std_deg(:);
    thetaSigmaDeg = thetaStd;
    if ismember('std_angle_count', T.Properties.VariableNames)
        cnt = max(T.std_angle_count(:), 1);
        thetaSigmaDeg = thetaStd ./ sqrt(cnt);
    end
end
if ismember('std_angle_iqr_deg', T.Properties.VariableNames)
    iqrSigma = T.std_angle_iqr_deg(:) ./ 1.349;
    thetaSigmaDeg = min_finite_local(thetaSigmaDeg, iqrSigma);
end
thetaSigmaDeg(~isfinite(thetaSigmaDeg)) = median(thetaSigmaDeg, 'omitnan');
if all(~isfinite(thetaSigmaDeg))
    thetaSigmaDeg = zeros(n, 1);
end
sigmaAngleMm = abs(thetaSigmaDeg) * pi / 180 * cfg.r_tip_mm;
sigmaMm = sqrt(cfg.step06j_weight_sigma_floor_mm.^2 + sigmaAngleMm.^2);
w = w .* (1 ./ max(sigmaMm.^2, eps));
w = w ./ median(w, 'omitnan');

if ismember('step02_best_corr', T.Properties.VariableNames)
    corrW = clamp_local((T.step02_best_corr(:) - 0.85) / 0.10, 0.20, 1.0);
    w = w .* corrW;
end
if ismember('step02_corr_gap', T.Properties.VariableNames)
    gapW = clamp_local(T.step02_corr_gap(:) / 0.20, 0.20, 1.0);
    w = w .* gapW;
end
if ismember('step02_shift_consistency', T.Properties.VariableNames)
    shiftW = clamp_local(T.step02_shift_consistency(:), 0.20, 1.0);
    w = w .* shiftW;
end
w(~isfinite(w) | w <= 0) = 1;
w = clamp_local(w ./ median(w, 'omitnan'), 0.10, 5.0);
end


function [beta, weights] = weighted_irls_ridge_local(X, y, baseWeights, penalty, cfg)
valid = all(isfinite(X), 2) & isfinite(y) & isfinite(baseWeights) & baseWeights > 0;
Xv = X(valid, :);
yv = y(valid);
baseW = baseWeights(valid);
robW = ones(size(yv));
p = size(X, 2);
beta = zeros(p, 1);
penaltyRows = zeros(0, p);
penCols = find(penalty > 0);
for i = 1:numel(penCols)
    row = zeros(1, p);
    row(penCols(i)) = penalty(penCols(i));
    penaltyRows(end+1, :) = row; %#ok<AGROW>
end
for iter = 1:cfg.step06i_irls_iterations
    w = baseW .* robW;
    sw = sqrt(w);
    Xa = Xv .* sw;
    ya = yv .* sw;
    if ~isempty(penaltyRows)
        Xa = [Xa; penaltyRows]; %#ok<AGROW>
        ya = [ya; zeros(size(penaltyRows, 1), 1)]; %#ok<AGROW>
    end
    beta = Xa \ ya;
    r = yv - Xv * beta;
    s = robust_scale_step06j_local(r);
    if ~isfinite(s) || s <= eps
        break;
    end
    cutoff = cfg.step06i_huber_k * s;
    robW = min(1, cutoff ./ max(abs(r), eps));
end
weights = zeros(size(y));
weights(valid) = baseW .* robW;
end


function ModelTable = build_model_table_local(varargin)
cfg = varargin{end};
targetEO = varargin{end-1};
fits = varargin(1:end-2);
rows = repmat(make_model_row_local(), numel(fits), 1);
for i = 1:numel(fits)
    F = fits{i};
    rows(i).model_name = F.model_name;
    rows(i).target_eo = targetEO;
    rows(i).point_count = F.point_count;
    rows(i).K_mm_per_microstrain = F.K_mm_per_microstrain;
    rows(i).K_um_per_microstrain = F.K_um_per_microstrain;
    rows(i).K_ratio_to_fixed = F.K_ratio_to_fixed;
    rows(i).phase_lag_rad = F.phase_lag_rad;
    rows(i).phase_lag_deg = F.phase_lag_deg;
    rows(i).u_rms_mm = F.u_rms_mm;
    rows(i).A_peak_equiv_mm = F.A_peak_equiv_mm;
    rows(i).A_half_range_mm = F.A_half_range_mm;
    rows(i).rmse_mm = F.rmse_mm;
    rows(i).weighted_rmse_mm = F.weighted_rmse_mm;
    rows(i).mean_abs_residual_mm = F.mean_abs_residual_mm;
    rows(i).phase_coverage_fraction = F.phase_coverage_fraction;
    rows(i).fixed_K_um_per_microstrain = 1000 * cfg.step06i_fixed_strain_to_mm;
    rows(i).phase_delta_to_main_rad = F.phase_delta_to_main_rad;
end
ModelTable = struct2table(rows);
end


function row = make_model_row_local()
row = struct('model_name', "", 'target_eo', NaN, 'point_count', NaN, ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'K_ratio_to_fixed', NaN, 'phase_lag_rad', NaN, ...
    'phase_lag_deg', NaN, 'u_rms_mm', NaN, ...
    'A_peak_equiv_mm', NaN, 'A_half_range_mm', NaN, ...
    'rmse_mm', NaN, 'weighted_rmse_mm', NaN, ...
    'mean_abs_residual_mm', NaN, 'phase_coverage_fraction', NaN, ...
    'fixed_K_um_per_microstrain', NaN, 'phase_delta_to_main_rad', NaN);
end


function LOS = build_leave_one_sensor_table_local(T, Tref, targetEO, cfg)
sids = unique(T.sensor_id(:)).';
rows = repmat(make_los_row_local(), numel(sids), 1);
for i = 1:numel(sids)
    sid = sids(i);
    keep = find(T.sensor_id ~= sid);
    rows(i).held_out_sensor_id = sid;
    rows(i).target_eo = targetEO;
    rows(i).train_point_count = numel(keep);
    if numel(keep) < 30 || numel(unique(T.sensor_id(keep))) < 2
        rows(i).status = "too_few_points";
        continue;
    end
    try
        F = fit_target_reference_model_local(T, Tref, cfg, ...
            "M0_no_sensor_bias", [], keep);
        rows(i).K_um_per_microstrain = F.K_um_per_microstrain;
        rows(i).phase_lag_rad = F.phase_lag_rad;
        rows(i).u_rms_mm = F.u_rms_mm;
        rows(i).A_peak_equiv_mm = F.A_peak_equiv_mm;
        rows(i).rmse_train_mm = F.rmse_mm;
        rows(i).status = "ok";
    catch ME
        rows(i).status = "failed_" + string(ME.identifier);
    end
end
LOS = struct2table(rows);
end


function row = make_los_row_local()
row = struct('held_out_sensor_id', NaN, 'target_eo', NaN, ...
    'train_point_count', NaN, 'K_um_per_microstrain', NaN, ...
    'phase_lag_rad', NaN, 'u_rms_mm', NaN, 'A_peak_equiv_mm', NaN, ...
    'rmse_train_mm', NaN, 'status', "");
end


function Boot = bootstrap_target_reference_local(T, Tref, targetEO, cfg)
rng(cfg.step06i_random_seed);
blocks = unique(T.revolution_index(isfinite(T.revolution_index))).';
if isempty(blocks)
    blocks = 1:height(T);
    T.revolution_index = blocks(:);
end
nBlocks = numel(blocks);
rows = repmat(make_boot_row_local(), cfg.step06j_bootstrap_count, 1);
for b = 1:cfg.step06j_bootstrap_count
    picked = blocks(randi(nBlocks, nBlocks, 1));
    idx = [];
    for k = 1:numel(picked)
        idx = [idx; find(T.revolution_index == picked(k))]; %#ok<AGROW>
    end
    rows(b).bootstrap_id = b;
    rows(b).target_eo = targetEO;
    rows(b).block_count = nBlocks;
    rows(b).point_count = numel(idx);
    try
        F = fit_target_reference_model_local(T, Tref, cfg, ...
            "M0_no_sensor_bias", [], idx);
        rows(b).K_mm_per_microstrain = F.K_mm_per_microstrain;
        rows(b).K_um_per_microstrain = F.K_um_per_microstrain;
        rows(b).phase_lag_rad = F.phase_lag_rad;
        rows(b).u_rms_mm = F.u_rms_mm;
        rows(b).A_peak_equiv_mm = F.A_peak_equiv_mm;
        rows(b).rmse_mm = F.rmse_mm;
        rows(b).C_real = real(F.C_complex_mm_per_microstrain);
        rows(b).C_imag = imag(F.C_complex_mm_per_microstrain);
        rows(b).status = "ok";
    catch ME
        rows(b).status = "failed_" + string(ME.identifier);
    end
end
Boot = struct2table(rows);
end


function row = make_boot_row_local()
row = struct('bootstrap_id', NaN, 'target_eo', NaN, 'block_count', NaN, ...
    'point_count', NaN, 'K_mm_per_microstrain', NaN, ...
    'K_um_per_microstrain', NaN, 'phase_lag_rad', NaN, ...
    'u_rms_mm', NaN, 'A_peak_equiv_mm', NaN, 'rmse_mm', NaN, ...
    'C_real', NaN, 'C_imag', NaN, 'status', "");
end


function TrefOut = build_strict_time_series_local(Tref, MainFit, Boot)
TrefOut = table();
TrefOut.time_s = Tref.time_s;
TrefOut.region_id = Tref.region_id;
TrefOut.dominant_order = Tref.dominant_order;
TrefOut.strain_basis_real_microstrain = Tref.strain_basis_real_microstrain;
TrefOut.strain_basis_imag_microstrain = Tref.strain_basis_imag_microstrain;
TrefOut.u_ref_strict_mm = MainFit.u_ref_mm;
TrefOut.u_ref_strict_mean_removed_mm = MainFit.u_ref_mm - ...
    mean(MainFit.u_ref_mm, 'omitnan');
if ismember('u_step05_mm', Tref.Properties.VariableNames)
    TrefOut.u_step05_mm = Tref.u_step05_mm;
end
ok = strcmp(string(Boot.status), "ok") & isfinite(Boot.C_real) & isfinite(Boot.C_imag);
if any(ok)
    h = Tref.strain_basis_real_microstrain + 1i * Tref.strain_basis_imag_microstrain;
    Cboot = Boot.C_real(ok) + 1i * Boot.C_imag(ok);
    U = real(h(:) * transpose(Cboot(:)));
    TrefOut.u_ref_boot_low_mm = prctile(U, 2.5, 2);
    TrefOut.u_ref_boot_high_mm = prctile(U, 97.5, 2);
else
    TrefOut.u_ref_boot_low_mm = nan(height(TrefOut), 1);
    TrefOut.u_ref_boot_high_mm = nan(height(TrefOut), 1);
end
end


function P = attach_strict_point_fit_local(T, MainFit)
P = T;
h = P.strain_basis_real_microstrain + 1i * P.strain_basis_imag_microstrain;
u = real(MainFit.C_complex_mm_per_microstrain .* h);
P.u_ref_strict_at_btt_mm = u;
P.strict_total_fit_mm = MainFit.beta(1) + u;
P.strict_residual_mm = P.displacement_mm - P.strict_total_fit_mm;
P.strict_dynamic_observed_mm = P.displacement_mm - MainFit.beta(1);
P.strict_fit_weight = MainFit.weights;
end


function Summary = build_strict_summary_local(caseName, targetEO, T, TS, ModelTable, LOS, Boot, cfg)
main = ModelTable(ModelTable.model_name == "M0_no_sensor_bias", :);
fixedBias = ModelTable(ModelTable.model_name == "M1_fixed_low_vibration_sensor_bias", :);
free = ModelTable(ModelTable.model_name == "M3_free_sensor_bias", :);
fixedK = ModelTable(ModelTable.model_name == "M4_fixed_K_main_phase", :);
okBoot = strcmp(string(Boot.status), "ok");
Kci = [NaN NaN];
Aci = [NaN NaN];
if any(okBoot)
    Kci = prctile(Boot.K_um_per_microstrain(okBoot), [2.5 97.5]);
    Aci = prctile(Boot.A_peak_equiv_mm(okBoot), [2.5 97.5]);
end
Summary = table();
Summary.case_name = string(caseName);
Summary.reference_definition = "strict_target_step03_btt_strain_shape_no_sensor_bias";
Summary.reference_uses_step05 = false;
Summary.step05_usage = "comparison_only";
Summary.target_time_start_s = cfg.step06i_target_time_range_s(1);
Summary.target_time_end_s = cfg.step06i_target_time_range_s(2);
Summary.target_blade = cfg.step06i_target_blade;
Summary.target_eo = targetEO;
Summary.sensor_ids = string(mat2str(unique(T.sensor_id(:)).'));
Summary.point_count = height(T);
Summary.revolution_block_count = numel(unique(T.revolution_index(isfinite(T.revolution_index))));
Summary.K_um_per_microstrain = main.K_um_per_microstrain;
Summary.K_um_per_microstrain_ci_low = Kci(1);
Summary.K_um_per_microstrain_ci_high = Kci(2);
Summary.K_ratio_to_fixed = main.K_ratio_to_fixed;
Summary.fixed_K_um_per_microstrain = 1000 * cfg.step06i_fixed_strain_to_mm;
Summary.phase_lag_rad = main.phase_lag_rad;
Summary.u_rms_mm = main.u_rms_mm;
Summary.A_peak_equiv_mm = main.A_peak_equiv_mm;
Summary.A_peak_equiv_ci_low_mm = Aci(1);
Summary.A_peak_equiv_ci_high_mm = Aci(2);
Summary.A_half_range_mm = main.A_half_range_mm;
Summary.rmse_mm = main.rmse_mm;
Summary.weighted_rmse_mm = main.weighted_rmse_mm;
Summary.fixed_low_vibration_bias_A_peak_equiv_mm = fixedBias.A_peak_equiv_mm;
Summary.fixed_low_vibration_bias_K_um_per_microstrain = fixedBias.K_um_per_microstrain;
Summary.free_sensor_bias_A_peak_equiv_mm = free.A_peak_equiv_mm;
Summary.free_sensor_bias_K_um_per_microstrain = free.K_um_per_microstrain;
Summary.fixed_K_A_peak_equiv_mm = fixedK.A_peak_equiv_mm;
Summary.leave_one_sensor_A_peak_min_mm = min(LOS.A_peak_equiv_mm, [], 'omitnan');
Summary.leave_one_sensor_A_peak_max_mm = max(LOS.A_peak_equiv_mm, [], 'omitnan');
if ismember('u_step05_mm', TS.Properties.VariableNames)
    mask = isfinite(TS.u_step05_mm) & isfinite(TS.u_ref_strict_mm);
    Summary.step05_time_rmse_mm = sqrt(mean((TS.u_step05_mm(mask) - ...
        TS.u_ref_strict_mm(mask)).^2, 'omitnan'));
    Summary.step05_time_corr = corr(TS.u_step05_mm(mask), ...
        TS.u_ref_strict_mm(mask), 'Rows', 'complete');
    Summary.step05_A_peak_equiv_mm = sqrt(2) * std(TS.u_step05_mm(mask), 'omitnan');
    Summary.step05_minus_strict_A_percent = 100 * ...
        (Summary.step05_A_peak_equiv_mm - Summary.A_peak_equiv_mm) / ...
        Summary.A_peak_equiv_mm;
else
    Summary.step05_time_rmse_mm = NaN;
    Summary.step05_time_corr = NaN;
    Summary.step05_A_peak_equiv_mm = NaN;
    Summary.step05_minus_strict_A_percent = NaN;
end
end


function figFile = plot_strict_reference_local(TS, P, Summary, ModelTable, LOS, Boot, figDir, cfg)
fig = figure('Name', 'Step06J strict target reference', ...
    'Color', 'w', 'Position', [80 60 1500 860], ...
    'Visible', cfg.step06j_fig_visible);
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(numel(unique(P.sensor_id)));
sids = unique(P.sensor_id(:)).';

ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
for i = 1:numel(sids)
    m = P.sensor_id == sids(i);
    scatter(ax1, P.arrival_time_s(m), P.displacement_mm(m), 10, ...
        colors(i, :), 'filled', 'DisplayName', sprintf('S%d raw BTT', sids(i)));
end
plot(ax1, TS.time_s, Summary.K_um_per_microstrain * 0 + ...
    mean(P.displacement_mm, 'omitnan') + TS.u_ref_strict_mm, ...
    'k-', 'LineWidth', 1.1, 'DisplayName', 'B + u_{ref}');
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'BTT displacement (mm)');
title(ax1, '1 Raw target-window BTT');
legend(ax1, 'Location', 'best');

ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
for i = 1:numel(sids)
    m = P.sensor_id == sids(i);
    scatter(ax2, P.arrival_time_s(m), P.strict_dynamic_observed_mm(m), 10, ...
        colors(i, :), 'filled', 'DisplayName', sprintf('S%d B removed', sids(i)));
end
plot(ax2, TS.time_s, TS.u_ref_strict_mm, 'k-', 'LineWidth', 1.2, ...
    'DisplayName', 'strict u_{ref}');
plot(ax2, TS.time_s, TS.u_ref_boot_low_mm, 'k--', 'LineWidth', 0.8, ...
    'DisplayName', 'bootstrap CI');
plot(ax2, TS.time_s, TS.u_ref_boot_high_mm, 'k--', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Dynamic displacement (mm)');
title(ax2, '2 Main reference without free sensor bias');
legend(ax2, 'Location', 'best');

ax3 = nexttile;
bar(ax3, categorical(ModelTable.model_name), ModelTable.A_peak_equiv_mm);
grid(ax3, 'on'); box(ax3, 'on');
ylabel(ax3, 'A_{peak,eq} (mm)');
title(ax3, '3 Model sensitivity');

ax4 = nexttile;
ok = strcmp(string(Boot.status), "ok");
histogram(ax4, Boot.A_peak_equiv_mm(ok), 18, 'FaceColor', [0.2 0.5 0.8]);
grid(ax4, 'on'); box(ax4, 'on');
xline(ax4, Summary.A_peak_equiv_mm, 'k-', 'LineWidth', 1.2);
xlabel(ax4, 'A_{peak,eq} (mm)'); ylabel(ax4, 'Count');
title(ax4, '4 Revolution-block bootstrap');

ax5 = nexttile;
yyaxis(ax5, 'left');
bar(ax5, categorical(string(LOS.held_out_sensor_id)), LOS.A_peak_equiv_mm);
ylabel(ax5, 'A_{peak,eq} (mm)');
yyaxis(ax5, 'right');
plot(ax5, categorical(string(LOS.held_out_sensor_id)), ...
    LOS.K_um_per_microstrain, 'ko-', 'LineWidth', 1.1);
ylabel(ax5, 'K (\mum/microstrain)');
grid(ax5, 'on'); box(ax5, 'on');
xlabel(ax5, 'Held-out sensor');
title(ax5, '5 Leave-one-sensor stability');

ax6 = nexttile;
hold(ax6, 'on'); grid(ax6, 'on'); box(ax6, 'on');
plot(ax6, TS.time_s, TS.u_ref_strict_mm, 'k-', 'LineWidth', 1.2, ...
    'DisplayName', 'strict reference');
if ismember('u_step05_mm', TS.Properties.VariableNames)
    plot(ax6, TS.time_s, TS.u_step05_mm, 'r--', 'LineWidth', 1.0, ...
        'DisplayName', 'Step05 comparison');
end
xlabel(ax6, 'Time (s)'); ylabel(ax6, 'Tip vibration (mm)');
title(ax6, sprintf('6 Strict reference vs Step05 | K=%.3f, A=%.3f', ...
    Summary.K_um_per_microstrain, Summary.A_peak_equiv_mm));
legend(ax6, 'Location', 'best');

sgtitle(fig, sprintf('Step06J strict target reference | EO%d | %.2f-%.2f s', ...
    Summary.target_eo, Summary.target_time_start_s, Summary.target_time_end_s));
figFile = fullfile(figDir, 'Step06J_StrictTargetReference_20250527.png');
save_figure_step06j_local(fig, figFile);
end


function targetEO = target_eo_from_tref_local(Tref)
targetEO = mode(Tref.dominant_order);
end


function x = min_finite_local(a, b)
x = a;
m = isfinite(b) & (~isfinite(a) | b < a);
x(m) = b(m);
end


function y = clamp_local(x, lo, hi)
y = min(max(x, lo), hi);
end


function s = robust_scale_step06j_local(r)
r = r(isfinite(r));
if isempty(r)
    s = NaN;
else
    s = 1.4826 * median(abs(r - median(r, 'omitnan')), 'omitnan');
    if ~isfinite(s) || s <= eps
        s = std(r, 'omitnan');
    end
end
end


function rev = compute_revolution_index_local(oprTimes, t, bladesNum)
prev = nan(size(t));
for i = 1:numel(t)
    idx = find(oprTimes < t(i), 1, 'last');
    if ~isempty(idx)
        prev(i) = idx;
    end
end
rev = floor((prev - 1) / bladesNum) + 1;
end


function oprTimes = load_opr_times_step06j_local(file)
S = load(file, 'jiluOPR');
if ~isfield(S, 'jiluOPR') || isempty(S.jiluOPR)
    error('Invalid jiluOPR file:\n  %s', file);
end
oprTimes = S.jiluOPR(:, 1);
oprTimes = oprTimes(isfinite(oprTimes));
end


function span = phase_coverage_step06j_local(phi)
phi = phi(isfinite(phi));
if numel(phi) < 2
    span = 0;
    return;
end
p = sort(mod(phi(:), 2*pi));
gaps = diff([p; p(1) + 2*pi]);
span = 2*pi - max(gaps);
end


function ensure_dir_step06j_local(d)
if ~isfolder(d)
    mkdir(d);
end
end


function save_figure_step06j_local(fig, pngFile)
[folder, base, ~] = fileparts(pngFile);
if ~isfolder(folder)
    mkdir(folder);
end
exportgraphics(fig, pngFile, 'Resolution', 220);
try
    exportgraphics(fig, fullfile(folder, [base, '.pdf']), 'ContentType', 'vector');
catch
end
try
    savefig(fig, fullfile(folder, [base, '.fig']));
catch
end
end
