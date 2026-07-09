%% Step06K_Build_SameFrequencyRobustReference_20250527.m
% Build a same-frequency, EO-aware reference for the target window.
%
% Purpose:
%   The 1-5 s EO14 target window has sparse BTT samples and poor phase
%   coverage.  Fitting only that window can make the gain absorb sensor
%   layering or sampling geometry.  This step keeps the target window for
%   the final comparison, but estimates the target EO gain from all
%   same-natural-frequency EO14 windows selected by Step06A/Step06I.
%
% Main target reference:
%   d_i = B_r + Re{ C_EO14 h_EO14(t_i) } + e_i
%   u_ref(t) = Re{ C_EO14 h_EO14(t) },  t in the target window
%
% Diagnostics:
%   1) per-region no-sensor-bias fits;
%   2) per-EO no-sensor-bias fits;
%   3) all-EO one-complex-gain fit, diagnostic only;
%   4) target-window-only fit, diagnostic only;
%   5) asynchronous amplitude-ratio checks that do not use strain/BTT
%      point-to-point phase;
%   6) sensor residuals after the main EO14 fit.

clc; clear; close all;

cfg = BTTDataConfig_20250527();
cfg = apply_step06k_defaults_local(cfg);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06K same-frequency robust reference: %s ===\n', caseName);

    files = resolve_step06k_files_local(cfg, caseName);
    assert_step06k_files_local(files);

    S06I = load(files.step06i_mat, 'FitPointsEO', 'ReferenceTimeSeries', ...
        'RegionTable', 'EOGainTable', 'Summary', 'files');
    ObsS = load(files.step03_bundle_mat, 'bundle');
    Obs = ObsS.bundle.Observation_Table;
    oprTimes = load_opr_times_step06k_local(files.jiluopr_mat);

    FitPointsAll = S06I.FitPointsEO;
    TrefAll = S06I.ReferenceTimeSeries;
    RegionTable = S06I.RegionTable;
    FitPointsAll = attach_quality_from_step03_step06k_local( ...
        FitPointsAll, Obs, cfg);

    [TargetPoints, TrefTarget, targetEO] = select_target_points_step06k_local( ...
        FitPointsAll, TrefAll, cfg);
    SameFreqPoints = select_same_frequency_points_step06k_local( ...
        FitPointsAll, RegionTable, cfg);
    SameFreqTref = select_same_frequency_tref_step06k_local( ...
        TrefAll, RegionTable);
    TargetEOPoints = SameFreqPoints(SameFreqPoints.dominant_order == targetEO, :);
    TargetEOTref = SameFreqTref(SameFreqTref.dominant_order == targetEO, :);

    RegionFits = build_region_fit_table_step06k_local( ...
        SameFreqPoints, SameFreqTref, TrefTarget, targetEO, cfg);
    EOFits = build_eo_fit_table_step06k_local( ...
        SameFreqPoints, SameFreqTref, cfg);
    AsyncRegionFits = build_async_region_fit_table_step06k_local( ...
        SameFreqPoints, SameFreqTref, oprTimes, targetEO, cfg);

    SameModeFit = fit_complex_reference_step06k_local( ...
        SameFreqPoints, SameFreqTref, cfg, "diagnostic_same_mode_all_regions");
    TargetEOFit = fit_complex_reference_step06k_local( ...
        TargetEOPoints, TargetEOTref, cfg, "diagnostic_target_EO_all_windows");
    TargetOnlyFit = fit_complex_reference_step06k_local( ...
        TargetPoints, TrefTarget, cfg, "diagnostic_target_window_only");
    [FEFit, FEReference, FEWidthSensitivity] = build_fe_reference_step06k_local( ...
        cfg, TrefTarget, TargetEOFit.C_complex_mm_per_microstrain);
    [PrimaryFit, PrimaryGrid] = build_fe_btt_grid_reference_step06k_local( ...
        TargetEOPoints, TrefTarget, TargetEOFit, FEReference, RegionTable, cfg);
    FEBandSensitivity = build_fe_band_sensitivity_step06k_local( ...
        TargetEOPoints, TrefTarget, TargetEOFit, FEReference, RegionTable, cfg);
    BttOnlyFit = fit_btt_harmonic_step06k_local( ...
        TargetPoints, TrefTarget, oprTimes, targetEO, cfg, ...
        "diagnostic_BTT_only_target_window");
    AsyncAmpFit = build_async_amp_scaled_fit_step06k_local( ...
        TrefTarget, oprTimes, BttOnlyFit, AsyncRegionFits, targetEO, cfg);

    BootstrapTable = bootstrap_primary_grid_reference_step06k_local( ...
        TargetEOPoints, TargetEOTref, TrefTarget, TargetEOFit, ...
        FEReference, RegionTable, targetEO, cfg);
    TargetTimeSeries = build_target_time_series_step06k_local( ...
        TrefTarget, PrimaryFit, SameModeFit, TargetEOFit, TargetOnlyFit, ...
        FEFit, BttOnlyFit, AsyncAmpFit, BootstrapTable, cfg);
    TargetFitPoints = attach_target_point_predictions_step06k_local( ...
        TargetPoints, PrimaryFit, SameModeFit, TargetOnlyFit, TargetEOFit, ...
        FEFit, BttOnlyFit, cfg);
    SensorResiduals = build_sensor_residuals_step06k_local( ...
        SameFreqPoints, SameModeFit, cfg);
    ModelComparison = build_model_comparison_step06k_local( ...
        PrimaryFit, SameModeFit, TargetEOFit, TargetOnlyFit, FEFit, ...
        BttOnlyFit, AsyncAmpFit, TargetTimeSeries, targetEO, cfg);
    Summary = build_summary_step06k_local(caseName, targetEO, ...
        RegionTable, RegionFits, EOFits, AsyncRegionFits, ModelComparison, ...
        BootstrapTable, TargetTimeSeries, TargetFitPoints, ...
        SensorResiduals, FEReference, FEWidthSensitivity, S06I, cfg);

    outDir = fullfile(cfg.step06k_output_dir, caseName);
    figDir = fullfile(cfg.step06k_figure_dir, caseName);
    ensure_dir_step06k_local(outDir);
    ensure_dir_step06k_local(figDir);

    summaryCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_Summary_20250527.csv');
    regionCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_RegionFits_20250527.csv');
    eoCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_EOFits_20250527.csv');
    modelCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_ModelComparison_20250527.csv');
    asyncRegionCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_AsyncRegionFits_20250527.csv');
    targetPointsCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_TargetFitPoints_20250527.csv');
    targetTsCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_TargetTimeSeries_20250527.csv');
    bootCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_PrimaryGridBootstrap_20250527.csv');
    sensorCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_SensorResiduals_20250527.csv');
    feCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_FEReference_20250527.csv');
    feWidthCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_FEWidthSensitivity_20250527.csv');
    gridCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_PrimaryGridSearch_20250527.csv');
    feBandSensitivityCsv = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_FEBandSensitivity_20250527.csv');
    matFile = fullfile(outDir, ...
        'Step06K_SameFrequencyRobustReference_20250527.mat');

    writetable(Summary, summaryCsv);
    writetable(RegionFits, regionCsv);
    writetable(EOFits, eoCsv);
    writetable(AsyncRegionFits, asyncRegionCsv);
    writetable(ModelComparison, modelCsv);
    writetable(TargetFitPoints, targetPointsCsv);
    writetable(TargetTimeSeries, targetTsCsv);
    writetable(BootstrapTable, bootCsv);
    writetable(SensorResiduals, sensorCsv);
    writetable(FEReference, feCsv);
    writetable(FEWidthSensitivity, feWidthCsv);
    writetable(PrimaryGrid, gridCsv);
    writetable(FEBandSensitivity, feBandSensitivityCsv);
    save(matFile, 'Summary', 'RegionFits', 'EOFits', 'AsyncRegionFits', ...
        'ModelComparison', 'TargetFitPoints', 'TargetTimeSeries', ...
        'BootstrapTable', 'SensorResiduals', 'FEReference', ...
        'FEWidthSensitivity', 'PrimaryGrid', 'FEBandSensitivity', ...
        'PrimaryFit', 'SameModeFit', ...
        'TargetEOFit', 'TargetOnlyFit', 'FEFit', 'BttOnlyFit', ...
        'AsyncAmpFit', 'RegionTable', ...
        'files', '-v7.3');

    [figFile, feFigFile] = plot_step06k_reference_local(TargetTimeSeries, ...
        TargetFitPoints, SameFreqPoints, RegionTable, RegionFits, ...
        EOFits, AsyncRegionFits, ModelComparison, SensorResiduals, ...
        BootstrapTable, FEWidthSensitivity, Summary, figDir, cfg);

    fprintf('\nStep06K same-frequency robust reference summary:\n');
    disp(Summary);
    fprintf(['Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        '  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n'], ...
        summaryCsv, regionCsv, eoCsv, asyncRegionCsv, modelCsv, targetPointsCsv, ...
        targetTsCsv, bootCsv, sensorCsv, feCsv, feWidthCsv, gridCsv, ...
        feBandSensitivityCsv, matFile, figFile, feFigFile);
end


function cfg = apply_step06k_defaults_local(cfg)
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
if ~isfield(cfg, 'step06k_bootstrap_count') || ...
        isempty(cfg.step06k_bootstrap_count)
    cfg.step06k_bootstrap_count = 0;
end
if ~isfield(cfg, 'step06k_weight_sigma_floor_mm') || ...
        isempty(cfg.step06k_weight_sigma_floor_mm)
    cfg.step06k_weight_sigma_floor_mm = 0.008;
end
if ~isfield(cfg, 'step06k_output_dir') || isempty(cfg.step06k_output_dir)
    cfg.step06k_output_dir = fullfile(cfg.output_root, ...
        'step06k_same_frequency_robust_reference');
end
if ~isfield(cfg, 'step06k_figure_dir') || isempty(cfg.step06k_figure_dir)
    cfg.step06k_figure_dir = fullfile(cfg.figure_root, ...
        'step06k_same_frequency_robust_reference');
end
if ~isfield(cfg, 'step06k_fig_visible') || isempty(cfg.step06k_fig_visible)
    if usejava('desktop')
        cfg.step06k_fig_visible = 'on';
    else
        cfg.step06k_fig_visible = 'off';
    end
end
if ~isfield(cfg, 'step06k_zoom_time_range_s') || ...
        isempty(cfg.step06k_zoom_time_range_s)
    cfg.step06k_zoom_time_range_s = [1.50 1.56];
end
if ~isfield(cfg, 'step06k_fe_transmissibility_per_m') || ...
        isempty(cfg.step06k_fe_transmissibility_per_m)
    cfg.step06k_fe_transmissibility_per_m = 0.9451;
end
if ~isfield(cfg, 'step06k_fe_strain_to_mm') || ...
        isempty(cfg.step06k_fe_strain_to_mm)
    cfg.step06k_fe_strain_to_mm = 1e-3 / cfg.step06k_fe_transmissibility_per_m;
end
if ~isfield(cfg, 'step06k_fe_source') || isempty(cfg.step06k_fe_source)
    cfg.step06k_fe_source = "ANSYS_580Hz_surface_neighborhood_interpolation";
end
if ~isfield(cfg, 'step06k_fe_point_label') || isempty(cfg.step06k_fe_point_label)
    cfg.step06k_fe_point_label = "strain_gauge_surface_neighborhood_Z2p82_Z22p18";
end
if ~isfield(cfg, 'step06k_fe_z_positions_mm') || ...
        isempty(cfg.step06k_fe_z_positions_mm)
    cfg.step06k_fe_z_positions_mm = [2.82 22.18];
end
if ~isfield(cfg, 'step06k_fe_nearest_transmissibility_per_m') || ...
        isempty(cfg.step06k_fe_nearest_transmissibility_per_m)
    cfg.step06k_fe_nearest_transmissibility_per_m = [0.9497 0.9497];
end
if ~isfield(cfg, 'step06k_fe_interpolated_transmissibility_per_m') || ...
        isempty(cfg.step06k_fe_interpolated_transmissibility_per_m)
    cfg.step06k_fe_interpolated_transmissibility_per_m = [0.9451 0.9451];
end
if ~isfield(cfg, 'step06k_fe_width_sensitivity_file') || ...
        isempty(cfg.step06k_fe_width_sensitivity_file)
    cfg.step06k_fe_width_sensitivity_file = fullfile( ...
        fileparts(fileparts(cfg.route_dir)), ...
        'blade_580_sweep', 'run_1p3127', 'z_width_sensitivity_580.json');
end
if ~isfield(cfg, 'step06k_fe_nominal_z_mm') || ...
        isempty(cfg.step06k_fe_nominal_z_mm)
    cfg.step06k_fe_nominal_z_mm = 2.82;
end
if ~isfield(cfg, 'step06k_fe_z_uncertainty_mm') || ...
        isempty(cfg.step06k_fe_z_uncertainty_mm)
    cfg.step06k_fe_z_uncertainty_mm = 2.5;
end
if ~isfield(cfg, 'step06k_fe_use_width_uncertainty_as_reference') || ...
        isempty(cfg.step06k_fe_use_width_uncertainty_as_reference)
    cfg.step06k_fe_use_width_uncertainty_as_reference = true;
end
if ~isfield(cfg, 'step06k_grid_frequency_half_width_hz') || ...
        isempty(cfg.step06k_grid_frequency_half_width_hz)
    cfg.step06k_grid_frequency_half_width_hz = 1.0;
end
if ~isfield(cfg, 'step06k_grid_frequency_count') || ...
        isempty(cfg.step06k_grid_frequency_count)
    cfg.step06k_grid_frequency_count = 17;
end
if ~isfield(cfg, 'step06k_grid_K_count') || isempty(cfg.step06k_grid_K_count)
    cfg.step06k_grid_K_count = 21;
end
if ~isfield(cfg, 'step06k_grid_phase_count') || ...
        isempty(cfg.step06k_grid_phase_count)
    cfg.step06k_grid_phase_count = 73;
end
if ~isfield(cfg, 'step06k_sensor_layering_penalty') || ...
        isempty(cfg.step06k_sensor_layering_penalty)
    cfg.step06k_sensor_layering_penalty = 0.0;
end
if ~isfield(cfg, 'step06k_fe_band_sensitivity_factors') || ...
        isempty(cfg.step06k_fe_band_sensitivity_factors)
    cfg.step06k_fe_band_sensitivity_factors = [1.0 1.5 2.0 3.0];
end
if ~isfield(cfg, 'step06k_sensitivity_frequency_count') || ...
        isempty(cfg.step06k_sensitivity_frequency_count)
    cfg.step06k_sensitivity_frequency_count = 13;
end
if ~isfield(cfg, 'step06k_sensitivity_K_count') || ...
        isempty(cfg.step06k_sensitivity_K_count)
    cfg.step06k_sensitivity_K_count = 21;
end
if ~isfield(cfg, 'step06k_sensitivity_phase_count') || ...
        isempty(cfg.step06k_sensitivity_phase_count)
    cfg.step06k_sensitivity_phase_count = 37;
end
end


function files = resolve_step06k_files_local(cfg, caseName)
files = struct();
files.step06i_mat = fullfile(cfg.output_root, ...
    'step06i_eo_dependent_joint_reference', caseName, ...
    'Step06I_EODependentReference_20250527.mat');
files.step03_bundle_mat = fullfile(cfg.step03_output_dir, caseName, ...
    'BTT_Observation_Bundle_20250527.mat');
files.jiluopr_mat = fullfile(cfg.step02_output_dir, caseName, 'jiluOPR.mat');
end


function assert_step06k_files_local(files)
names = fieldnames(files);
for i = 1:numel(names)
    if ~isfile(files.(names{i}))
        error('Missing Step06K input %s:\n  %s', names{i}, files.(names{i}));
    end
end
end


function T = attach_quality_from_step03_step06k_local(T, Obs, cfg)
T = T(T.blade_id == cfg.step06i_target_blade & ...
    ismember(T.sensor_id, cfg.step06i_sensor_ids), :);
T.row_id_step06k = (1:height(T)).';

qVars = {'revolution_index', 'std_angle_std_deg', 'std_angle_iqr_deg', ...
    'std_angle_count', 'step02_best_corr', 'step02_corr_gap', ...
    'step02_shift_consistency', 'step02_label_quality', 'quality_flag'};
qVars = qVars(ismember(qVars, Obs.Properties.VariableNames));
ObsQ = Obs(Obs.blade_id == cfg.step06i_target_blade & ...
    ismember(Obs.sensor_id, cfg.step06i_sensor_ids), ...
    [{'sensor_id', 'blade_id', 'arrival_time_s'}, qVars]);

T.join_key_step06k = int64(round(T.arrival_time_s * 1e9));
ObsQ.join_key_step06k = int64(round(ObsQ.arrival_time_s * 1e9));
ObsQ.arrival_time_s = [];
Tj = innerjoin(T, ObsQ, ...
    'Keys', {'join_key_step06k', 'sensor_id', 'blade_id'});
if height(Tj) < height(T)
    warning('Step06K quality join kept %d/%d fit points.', ...
        height(Tj), height(T));
end
Tj = sortrows(Tj, 'row_id_step06k');
Tj.join_key_step06k = [];
T = Tj;
end


function [Target, TrefTarget, targetEO] = select_target_points_step06k_local( ...
    FitPoints, TrefAll, cfg)
tw = cfg.step06i_target_time_range_s(:).';
mask = FitPoints.arrival_time_s >= tw(1) & ...
    FitPoints.arrival_time_s <= tw(2) & ...
    FitPoints.blade_id == cfg.step06i_target_blade & ...
    ismember(FitPoints.sensor_id, cfg.step06i_sensor_ids) & ...
    isfinite(FitPoints.displacement_mm) & ...
    isfinite(FitPoints.strain_basis_real_microstrain) & ...
    isfinite(FitPoints.strain_basis_imag_microstrain);
if ~any(mask)
    error('No target-window BTT fit points in [%.3f, %.3f] s.', ...
        tw(1), tw(2));
end
targetEO = mode(FitPoints.dominant_order(mask));
Target = FitPoints(mask & FitPoints.dominant_order == targetEO, :);
Target = sortrows(Target, 'arrival_time_s');

tmask = TrefAll.time_s >= tw(1) & TrefAll.time_s <= tw(2) & ...
    TrefAll.dominant_order == targetEO & ...
    isfinite(TrefAll.strain_basis_real_microstrain) & ...
    isfinite(TrefAll.strain_basis_imag_microstrain);
TrefTarget = sortrows(TrefAll(tmask, :), 'time_s');
if height(TrefTarget) < 10
    error('Too few target time samples for Step06K.');
end
end


function T = select_same_frequency_points_step06k_local(FitPoints, RegionTable, cfg)
regionIds = RegionTable.region_id(:);
T = FitPoints(ismember(FitPoints.region_id, regionIds) & ...
    FitPoints.blade_id == cfg.step06i_target_blade & ...
    ismember(FitPoints.sensor_id, cfg.step06i_sensor_ids) & ...
    isfinite(FitPoints.displacement_mm) & ...
    isfinite(FitPoints.strain_basis_real_microstrain) & ...
    isfinite(FitPoints.strain_basis_imag_microstrain), :);
T = sortrows(T, {'arrival_time_s', 'sensor_id'});
if height(T) < 40
    error('Too few same-frequency BTT points for Step06K.');
end
end


function Tref = select_same_frequency_tref_step06k_local(TrefAll, RegionTable)
regionIds = RegionTable.region_id(:);
Tref = TrefAll(ismember(TrefAll.region_id, regionIds) & ...
    isfinite(TrefAll.strain_basis_real_microstrain) & ...
    isfinite(TrefAll.strain_basis_imag_microstrain), :);
Tref = sortrows(Tref, 'time_s');
end


function RegionFits = build_region_fit_table_step06k_local( ...
    T, Tref, TrefTarget, targetEO, cfg)
regionIds = unique(T.region_id(:)).';
rows = repmat(make_region_fit_row_step06k_local(), numel(regionIds), 1);
for i = 1:numel(regionIds)
    rid = regionIds(i);
    Tp = T(T.region_id == rid, :);
    Tr = Tref(Tref.region_id == rid, :);
    rows(i).region_id = rid;
    rows(i).dominant_order = mode(Tp.dominant_order);
    rows(i).time_start_s = min(Tp.arrival_time_s, [], 'omitnan');
    rows(i).time_end_s = max(Tp.arrival_time_s, [], 'omitnan');
    rows(i).point_count = height(Tp);
    rows(i).sensor_ids = string(mat2str(unique(Tp.sensor_id(:)).'));
    try
        F = fit_complex_reference_step06k_local(Tp, Tr, cfg, ...
            "diagnostic_per_region");
        rows(i) = fill_region_fit_row_step06k_local(rows(i), F, Tp, Tr, ...
            TrefTarget, targetEO, cfg);
        rows(i).status = "ok";
    catch ME
        rows(i).status = "failed_" + string(ME.identifier);
    end
end
RegionFits = struct2table(rows);
RegionFits = sortrows(RegionFits, 'time_start_s');
end


function row = fill_region_fit_row_step06k_local(row, F, T, Tref, ...
    TrefTarget, targetEO, cfg)
row.K_mm_per_microstrain = F.K_mm_per_microstrain;
row.K_um_per_microstrain = F.K_um_per_microstrain;
row.K_ratio_to_fixed = F.K_ratio_to_fixed;
row.phase_lag_rad = F.phase_lag_rad;
row.phase_lag_deg = F.phase_lag_deg;
row.u_rms_mm = F.u_rms_mm;
row.A_peak_equiv_mm = F.A_peak_equiv_mm;
row.A_half_range_mm = F.A_half_range_mm;
row.rmse_mm = F.rmse_mm;
row.weighted_rmse_mm = F.weighted_rmse_mm;
row.mean_abs_residual_mm = F.mean_abs_residual_mm;
row.phase_coverage_fraction = F.phase_coverage_fraction;
row.rank_X = F.rank_X;
row.cond_X = F.cond_X;
row.strain_basis_peak_microstrain = max(Tref.strain_basis_abs_microstrain, [], ...
    'omitnan');
row.dynamic_BTT_std_mm = std(T.displacement_mm - ...
    median(T.displacement_mm, 'omitnan'), 'omitnan');
if row.dominant_order == targetEO
    hTarget = TrefTarget.strain_basis_real_microstrain + ...
        1i*TrefTarget.strain_basis_imag_microstrain;
    uTarget = real(F.C_complex_mm_per_microstrain .* hTarget);
    row.target_projection_A_peak_equiv_mm = sqrt(2) * ...
        std(uTarget, 'omitnan');
    row.target_projection_A_half_range_mm = 0.5 * ...
        (max(uTarget, [], 'omitnan') - min(uTarget, [], 'omitnan'));
else
    row.target_projection_A_peak_equiv_mm = NaN;
    row.target_projection_A_half_range_mm = NaN;
end
end


function row = make_region_fit_row_step06k_local()
row = struct('region_id', NaN, 'dominant_order', NaN, ...
    'time_start_s', NaN, 'time_end_s', NaN, ...
    'point_count', NaN, 'sensor_ids', "", ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'K_ratio_to_fixed', NaN, 'phase_lag_rad', NaN, ...
    'phase_lag_deg', NaN, 'u_rms_mm', NaN, ...
    'A_peak_equiv_mm', NaN, 'A_half_range_mm', NaN, ...
    'rmse_mm', NaN, 'weighted_rmse_mm', NaN, ...
    'mean_abs_residual_mm', NaN, 'phase_coverage_fraction', NaN, ...
    'rank_X', NaN, 'cond_X', NaN, ...
    'strain_basis_peak_microstrain', NaN, ...
    'dynamic_BTT_std_mm', NaN, ...
    'target_projection_A_peak_equiv_mm', NaN, ...
    'target_projection_A_half_range_mm', NaN, 'status', "");
end


function EOFits = build_eo_fit_table_step06k_local(T, Tref, cfg)
eoIds = unique(T.dominant_order(:)).';
rows = repmat(make_eo_fit_row_step06k_local(), numel(eoIds), 1);
for i = 1:numel(eoIds)
    eo = eoIds(i);
    Tp = T(T.dominant_order == eo, :);
    Tr = Tref(Tref.dominant_order == eo, :);
    rows(i).dominant_order = eo;
    rows(i).region_ids = string(mat2str(unique(Tp.region_id(:)).'));
    rows(i).region_count = numel(unique(Tp.region_id));
    rows(i).point_count = height(Tp);
    rows(i).sensor_ids = string(mat2str(unique(Tp.sensor_id(:)).'));
    try
        F = fit_complex_reference_step06k_local(Tp, Tr, cfg, ...
            "same_frequency_EO_fit");
        rows(i) = fill_eo_fit_row_step06k_local(rows(i), F, cfg);
        rows(i).status = "ok";
    catch ME
        rows(i).status = "failed_" + string(ME.identifier);
    end
end
EOFits = struct2table(rows);
EOFits = sortrows(EOFits, 'dominant_order');
end


function row = fill_eo_fit_row_step06k_local(row, F, cfg)
row.K_mm_per_microstrain = F.K_mm_per_microstrain;
row.K_um_per_microstrain = F.K_um_per_microstrain;
row.K_ratio_to_fixed = F.K_ratio_to_fixed;
row.phase_lag_rad = F.phase_lag_rad;
row.phase_lag_deg = F.phase_lag_deg;
row.u_rms_mm = F.u_rms_mm;
row.A_peak_equiv_mm = F.A_peak_equiv_mm;
row.A_half_range_mm = F.A_half_range_mm;
row.rmse_mm = F.rmse_mm;
row.weighted_rmse_mm = F.weighted_rmse_mm;
row.mean_abs_residual_mm = F.mean_abs_residual_mm;
row.phase_coverage_fraction = F.phase_coverage_fraction;
row.rank_X = F.rank_X;
row.cond_X = F.cond_X;
end


function row = make_eo_fit_row_step06k_local()
row = struct('dominant_order', NaN, 'region_ids', "", ...
    'region_count', NaN, 'point_count', NaN, 'sensor_ids', "", ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'K_ratio_to_fixed', NaN, 'phase_lag_rad', NaN, ...
    'phase_lag_deg', NaN, 'u_rms_mm', NaN, ...
    'A_peak_equiv_mm', NaN, 'A_half_range_mm', NaN, ...
    'rmse_mm', NaN, 'weighted_rmse_mm', NaN, ...
    'mean_abs_residual_mm', NaN, 'phase_coverage_fraction', NaN, ...
    'rank_X', NaN, 'cond_X', NaN, 'status', "");
end


function Async = build_async_region_fit_table_step06k_local( ...
    T, Tref, oprTimes, targetEO, cfg)
regionIds = unique(T.region_id(:)).';
rows = repmat(make_async_region_row_step06k_local(), numel(regionIds), 1);
for i = 1:numel(regionIds)
    rid = regionIds(i);
    Tp = T(T.region_id == rid, :);
    Tr = Tref(Tref.region_id == rid, :);
    eo = mode(Tp.dominant_order);
    rows(i).region_id = rid;
    rows(i).dominant_order = eo;
    rows(i).time_start_s = min(Tp.arrival_time_s, [], 'omitnan');
    rows(i).time_end_s = max(Tp.arrival_time_s, [], 'omitnan');
    rows(i).point_count = height(Tp);
    rows(i).sensor_ids = string(mat2str(unique(Tp.sensor_id(:)).'));
    try
        Bfit = fit_btt_harmonic_step06k_local(Tp, Tr, oprTimes, eo, cfg, ...
            "async_region_BTT_harmonic");
        [Astrain, strainRms] = strain_region_amplitude_step06k_local(Tr);
        rows(i).btt_A_peak_equiv_mm = Bfit.A_peak_equiv_mm;
        rows(i).btt_u_rms_mm = Bfit.u_rms_mm;
        rows(i).strain_A_peak_equiv_microstrain = Astrain;
        rows(i).strain_rms_microstrain = strainRms;
        rows(i).K_amp_um_per_microstrain = 1000 * ...
            Bfit.A_peak_equiv_mm / max(Astrain, eps);
        rows(i).phase_coverage_fraction = Bfit.phase_coverage_fraction;
        rows(i).rmse_mm = Bfit.rmse_mm;
        rows(i).weighted_rmse_mm = Bfit.weighted_rmse_mm;
        if eo == targetEO
            rows(i).target_eo_support = true;
        else
            rows(i).target_eo_support = false;
        end
        rows(i).status = "ok";
    catch ME
        rows(i).status = "failed_" + string(ME.identifier);
    end
end
Async = struct2table(rows);
Async = sortrows(Async, 'time_start_s');
end


function row = make_async_region_row_step06k_local()
row = struct('region_id', NaN, 'dominant_order', NaN, ...
    'time_start_s', NaN, 'time_end_s', NaN, ...
    'point_count', NaN, 'sensor_ids', "", ...
    'btt_A_peak_equiv_mm', NaN, 'btt_u_rms_mm', NaN, ...
    'strain_A_peak_equiv_microstrain', NaN, ...
    'strain_rms_microstrain', NaN, ...
    'K_amp_um_per_microstrain', NaN, ...
    'phase_coverage_fraction', NaN, 'rmse_mm', NaN, ...
    'weighted_rmse_mm', NaN, ...
    'target_eo_support', false, 'status', "");
end


function [ApeakEq, rmsVal] = strain_region_amplitude_step06k_local(Tref)
h = Tref.strain_basis_real_microstrain + 1i*Tref.strain_basis_imag_microstrain;
u = real(h);
rmsVal = std(u, 'omitnan');
ApeakEq = sqrt(2) * rmsVal;
if ~isfinite(ApeakEq) || ApeakEq <= eps
    ApeakEq = sqrt(2) * std(abs(h), 'omitnan');
end
end


function Fit = fit_complex_reference_step06k_local(T, Tref, cfg, modelName, idxRows)
if nargin < 5 || isempty(idxRows)
    idxRows = (1:height(T)).';
end
Ts = T(idxRows, :);
h = Ts.strain_basis_real_microstrain + 1i * Ts.strain_basis_imag_microstrain;
y = Ts.displacement_mm(:);
[X, meta] = build_design_step06k_local(Ts);
baseWeights = build_quality_weights_step06k_local(Ts, cfg);
[beta, weights] = weighted_irls_step06k_local(X, y, baseWeights, cfg);

C = beta(1) + 1i * beta(2);
yFit = X * beta;
resid = y - yFit;
offsetAtPoints = compute_region_offsets_step06k_local(Ts, beta, meta);

href = Tref.strain_basis_real_microstrain + 1i * ...
    Tref.strain_basis_imag_microstrain;
uRef = real(C .* href);

Fit = make_fit_struct_step06k_local(modelName, Ts, Tref, beta, C, ...
    yFit, resid, offsetAtPoints, weights, baseWeights, meta, uRef, cfg);
end


function [X, meta] = build_design_step06k_local(T)
h = T.strain_basis_real_microstrain + 1i * T.strain_basis_imag_microstrain;
X = [real(h), -imag(h)];
meta = struct();
meta.regionIds = unique(T.region_id(:)).';
meta.regionCols = zeros(size(meta.regionIds));
for i = 1:numel(meta.regionIds)
    rid = meta.regionIds(i);
    X(:, end+1) = double(T.region_id == rid); %#ok<AGROW>
    meta.regionCols(i) = size(X, 2);
end
end


function off = compute_region_offsets_step06k_local(T, beta, meta)
off = zeros(height(T), 1);
for i = 1:numel(meta.regionIds)
    m = T.region_id == meta.regionIds(i);
    off(m) = beta(meta.regionCols(i));
end
end


function Fit = make_fit_struct_step06k_local(modelName, T, Tref, beta, C, ...
    yFit, resid, offsetAtPoints, weights, baseWeights, meta, uRef, cfg)
uRms = std(uRef, 'omitnan');
Fit = struct();
Fit.model_name = string(modelName);
Fit.beta = beta;
Fit.C_complex_mm_per_microstrain = C;
Fit.K_mm_per_microstrain = abs(C);
Fit.K_um_per_microstrain = 1000 * abs(C);
Fit.K_ratio_to_fixed = abs(C) / cfg.step06i_fixed_strain_to_mm;
Fit.phase_lag_rad = angle(C);
Fit.phase_lag_deg = angle(C) * 180 / pi;
Fit.y_fit_mm = yFit;
Fit.residual_mm = resid;
Fit.offset_at_points_mm = offsetAtPoints;
Fit.u_at_points_mm = yFit - offsetAtPoints;
Fit.displacement_dynamic_mm = T.displacement_mm - offsetAtPoints;
Fit.weights = weights;
Fit.base_weights = baseWeights;
Fit.rmse_mm = sqrt(mean(resid.^2, 'omitnan'));
Fit.weighted_rmse_mm = sqrt(sum(weights .* resid.^2) / sum(weights));
Fit.mean_abs_residual_mm = mean(abs(resid), 'omitnan');
Fit.u_ref_mm = uRef;
Fit.u_rms_mm = uRms;
Fit.A_peak_equiv_mm = sqrt(2) * uRms;
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - ...
    min(uRef, [], 'omitnan'));
Fit.point_count = height(T);
Fit.region_ids = meta.regionIds;
Fit.region_offsets_mm = beta(meta.regionCols);
Fit.region_count = numel(meta.regionIds);
Fit.sensor_count = numel(unique(T.sensor_id));
Fit.phase_coverage_fraction = phase_coverage_step06k_local(angle( ...
    T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain)) ...
    / (2*pi);
Fit.rank_X = rank(build_design_step06k_local(T));
Fit.cond_X = cond(build_design_step06k_local(T));
if istable(Tref) && ismember('time_s', Tref.Properties.VariableNames)
    Fit.time_s = Tref.time_s;
else
    Fit.time_s = T.arrival_time_s;
end
if istable(Tref) && ismember('region_id', Tref.Properties.VariableNames)
    Fit.region_id_tref = Tref.region_id;
else
    Fit.region_id_tref = T.region_id;
end
end


function Fit = fit_fixed_k_step06k_local(T, Tref, cfg, phaseSourceC)
h = T.strain_basis_real_microstrain + 1i * T.strain_basis_imag_microstrain;
y = T.displacement_mm(:);
C = cfg.step06i_fixed_strain_to_mm * exp(1i * angle(phaseSourceC));
u = real(C .* h);
regionIds = unique(T.region_id(:)).';
Xoff = zeros(height(T), numel(regionIds));
for i = 1:numel(regionIds)
    Xoff(:, i) = double(T.region_id == regionIds(i));
end
w = build_quality_weights_step06k_local(T, cfg);
sw = sqrt(w);
B = (Xoff .* sw) \ ((y - u) .* sw);
off = Xoff * B;
yFit = u + off;
resid = y - yFit;
href = Tref.strain_basis_real_microstrain + 1i * ...
    Tref.strain_basis_imag_microstrain;
uRef = real(C .* href);
beta = [real(C); imag(C); B(:)];
meta = struct('regionIds', regionIds, ...
    'regionCols', 3:(2+numel(regionIds)));
Fit = make_fit_struct_step06k_local("diagnostic_fixed_K_main_phase", ...
    T, Tref, beta, C, yFit, resid, off, w, w, meta, uRef, cfg);
end


function Fit = build_fe_btt_constrained_reference_step06k_local( ...
    T, Tref, BttFit, FEReference, cfg)
phaseC = exp(1i * angle(BttFit.C_complex_mm_per_microstrain));
kExp = BttFit.K_mm_per_microstrain;
kLo = FEReference.band_K_um_per_microstrain_min(1) / 1000;
kHi = FEReference.band_K_um_per_microstrain_max(1) / 1000;
if ~isfinite(kLo) || ~isfinite(kHi) || kLo > kHi
    kLo = FEReference.K_mm_per_microstrain(1);
    kHi = FEReference.K_mm_per_microstrain(1);
end
kConstrained = min(max(kExp, kLo), kHi);
C = kConstrained * phaseC;

h = T.strain_basis_real_microstrain + 1i * T.strain_basis_imag_microstrain;
y = T.displacement_mm(:);
u = real(C .* h);
regionIds = unique(T.region_id(:)).';
Xoff = zeros(height(T), numel(regionIds));
for i = 1:numel(regionIds)
    Xoff(:, i) = double(T.region_id == regionIds(i));
end
w = build_quality_weights_step06k_local(T, cfg);
sw = sqrt(w);
B = (Xoff .* sw) \ ((y - u) .* sw);
off = Xoff * B;
yFit = u + off;
resid = y - yFit;

href = Tref.strain_basis_real_microstrain + 1i * ...
    Tref.strain_basis_imag_microstrain;
uRef = real(C .* href);
beta = [real(C); imag(C); B(:)];
meta = struct('regionIds', regionIds, ...
    'regionCols', 3:(2+numel(regionIds)));
Fit = make_fit_struct_step06k_local( ...
    "primary_FE_band_constrained_EO14_BTT", T, Tref, beta, C, yFit, ...
    resid, off, w, w, meta, uRef, cfg);
Fit.unconstrained_BTT_K_um_per_microstrain = 1000 * kExp;
Fit.FE_band_K_um_per_microstrain_min = 1000 * kLo;
Fit.FE_band_K_um_per_microstrain_max = 1000 * kHi;
Fit.constraint_status = "BTT_K_inside_FE_band";
if kExp < kLo
    Fit.constraint_status = "BTT_K_below_FE_band_clipped_to_lower_bound";
elseif kExp > kHi
    Fit.constraint_status = "BTT_K_above_FE_band_clipped_to_upper_bound";
end
Fit.source = "EO14 BTT phase and gain constrained by FE gauge-position band";
end


function [Fit, Grid] = build_fe_btt_grid_reference_step06k_local( ...
    T, Tref, BttFit, FEReference, RegionTable, cfg)
kLo = FEReference.band_K_um_per_microstrain_min(1) / 1000;
kHi = FEReference.band_K_um_per_microstrain_max(1) / 1000;
if ~isfinite(kLo) || ~isfinite(kHi) || kLo > kHi
    kLo = FEReference.K_mm_per_microstrain(1);
    kHi = FEReference.K_mm_per_microstrain(1);
end

freq0 = estimate_target_frequency_step06k_local(T, Tref, RegionTable);
freqGrid = linspace(freq0 - cfg.step06k_grid_frequency_half_width_hz, ...
    freq0 + cfg.step06k_grid_frequency_half_width_hz, ...
    cfg.step06k_grid_frequency_count);
freqGrid = freqGrid(isfinite(freqGrid) & freqGrid > 0);
if isempty(freqGrid)
    freqGrid = freq0;
end
kGrid = linspace(kLo, kHi, cfg.step06k_grid_K_count);
phase0 = angle(BttFit.C_complex_mm_per_microstrain);
phaseGrid = phase0 + linspace(-pi, pi, cfg.step06k_grid_phase_count);

t0 = mean(T.arrival_time_s, 'omitnan');
rows = repmat(make_primary_grid_row_step06k_local(), ...
    numel(freqGrid) * numel(kGrid) * numel(phaseGrid), 1);
rowIdx = 0;
for iF = 1:numel(freqGrid)
    fHz = freqGrid(iF);
    phaseAtBtt = 2*pi*fHz*(T.arrival_time_s(:) - t0);
    hAtBtt = abs(T.strain_basis_real_microstrain + ...
        1i*T.strain_basis_imag_microstrain) .* exp(1i*phaseAtBtt);
    phaseAtRef = 2*pi*fHz*(Tref.time_s(:) - t0);
    hAtRef = Tref.strain_basis_abs_microstrain .* exp(1i*phaseAtRef);
    for iK = 1:numel(kGrid)
        kVal = kGrid(iK);
        for iP = 1:numel(phaseGrid)
            rowIdx = rowIdx + 1;
            phaseVal = phaseGrid(iP);
            C = kVal * exp(1i*phaseVal);
            [score, rmse, wrmse, layer, phaseCov, Aref, resid] = ...
                score_constrained_candidate_step06k_local(T, hAtBtt, ...
                hAtRef, C, cfg);
            rows(rowIdx).frequency_hz = fHz;
            rows(rowIdx).K_mm_per_microstrain = kVal;
            rows(rowIdx).K_um_per_microstrain = 1000*kVal;
            rows(rowIdx).phase_lag_rad = wrap_to_pi_step06k_local(phaseVal);
            rows(rowIdx).objective = score;
            rows(rowIdx).rmse_mm = rmse;
            rows(rowIdx).weighted_rmse_mm = wrmse;
            rows(rowIdx).sensor_layering_score_mm = layer;
            rows(rowIdx).phase_coverage_fraction = phaseCov;
            rows(rowIdx).target_A_peak_equiv_mm = Aref;
            rows(rowIdx).residual_median_abs_mm = median(abs(resid), 'omitnan');
        end
    end
end
Grid = struct2table(rows);
Grid = sortrows(Grid, 'objective');

best = Grid(1, :);
phaseBestAtBtt = 2*pi*best.frequency_hz(1)*(T.arrival_time_s(:) - t0);
hBestAtBtt = abs(T.strain_basis_real_microstrain + ...
    1i*T.strain_basis_imag_microstrain) .* exp(1i*phaseBestAtBtt);
phaseBestAtRef = 2*pi*best.frequency_hz(1)*(Tref.time_s(:) - t0);
hBestAtRef = Tref.strain_basis_abs_microstrain .* exp(1i*phaseBestAtRef);
Cbest = best.K_mm_per_microstrain(1) * exp(1i*best.phase_lag_rad(1));
Fit = build_fixed_complex_fit_from_basis_step06k_local(T, Tref, ...
    hBestAtBtt, hBestAtRef, Cbest, cfg, ...
    "primary_FE_prior_frequency_grid_BTT_reference");
Fit.frequency_hz = best.frequency_hz(1);
Fit.grid_t0_s = t0;
Fit.grid_basis_type = "abs_strain_envelope_with_frequency_phase";
Fit.frequency_prior_center_hz = freq0;
Fit.frequency_prior_half_width_hz = cfg.step06k_grid_frequency_half_width_hz;
Fit.objective = best.objective(1);
Fit.sensor_layering_score_mm = best.sensor_layering_score_mm(1);
Fit.unconstrained_BTT_K_um_per_microstrain = BttFit.K_um_per_microstrain;
Fit.FE_band_K_um_per_microstrain_min = 1000*kLo;
Fit.FE_band_K_um_per_microstrain_max = 1000*kHi;
Fit.constraint_status = "grid_search_inside_FE_band";
Fit.source = "FE amplitude band + strain frequency grid + BTT residual search";
end


function Sens = build_fe_band_sensitivity_step06k_local( ...
    T, Tref, BttFit, FEReference, RegionTable, cfg)
factors = cfg.step06k_fe_band_sensitivity_factors(:).';
kLo0 = FEReference.band_K_um_per_microstrain_min(1);
kHi0 = FEReference.band_K_um_per_microstrain_max(1);
kCenter = 0.5 * (kLo0 + kHi0);
kHalf = 0.5 * (kHi0 - kLo0);
rows = repmat(make_fe_band_sensitivity_row_step06k_local(), numel(factors), 1);
sensCfg = cfg;
sensCfg.step06k_grid_frequency_count = cfg.step06k_sensitivity_frequency_count;
sensCfg.step06k_grid_K_count = cfg.step06k_sensitivity_K_count;
sensCfg.step06k_grid_phase_count = cfg.step06k_sensitivity_phase_count;
for i = 1:numel(factors)
    fac = factors(i);
    FEtmp = FEReference;
    kLo = max(eps, kCenter - fac * kHalf);
    kHi = kCenter + fac * kHalf;
    FEtmp.band_K_um_per_microstrain_min(1) = kLo;
    FEtmp.band_K_um_per_microstrain_max(1) = kHi;
    [Fit, ~] = build_fe_btt_grid_reference_step06k_local( ...
        T, Tref, BttFit, FEtmp, RegionTable, sensCfg);
    tol = 0.5 * (kHi - kLo) / max(sensCfg.step06k_grid_K_count - 1, 1);
    rows(i).band_expansion_factor = fac;
    rows(i).K_low_um_per_microstrain = kLo;
    rows(i).K_high_um_per_microstrain = kHi;
    rows(i).best_K_um_per_microstrain = Fit.K_um_per_microstrain;
    rows(i).best_A_peak_equiv_mm = Fit.A_peak_equiv_mm;
    rows(i).best_frequency_hz = Fit.frequency_hz;
    rows(i).objective = Fit.objective;
    rows(i).weighted_rmse_mm = Fit.weighted_rmse_mm;
    rows(i).sensor_layering_score_mm = Fit.sensor_layering_score_mm;
    rows(i).is_at_lower_bound = abs(Fit.K_um_per_microstrain - kLo) <= tol;
    rows(i).is_at_upper_bound = abs(Fit.K_um_per_microstrain - kHi) <= tol;
end
Sens = struct2table(rows);
end


function row = make_fe_band_sensitivity_row_step06k_local()
row = struct('band_expansion_factor', NaN, ...
    'K_low_um_per_microstrain', NaN, ...
    'K_high_um_per_microstrain', NaN, ...
    'best_K_um_per_microstrain', NaN, ...
    'best_A_peak_equiv_mm', NaN, ...
    'best_frequency_hz', NaN, ...
    'objective', NaN, ...
    'weighted_rmse_mm', NaN, ...
    'sensor_layering_score_mm', NaN, ...
    'is_at_lower_bound', false, ...
    'is_at_upper_bound', false);
end


function row = make_primary_grid_row_step06k_local()
row = struct('frequency_hz', NaN, ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'phase_lag_rad', NaN, 'objective', NaN, ...
    'rmse_mm', NaN, 'weighted_rmse_mm', NaN, ...
    'sensor_layering_score_mm', NaN, ...
    'phase_coverage_fraction', NaN, ...
    'target_A_peak_equiv_mm', NaN, ...
    'residual_median_abs_mm', NaN);
end


function fHz = estimate_target_frequency_step06k_local(T, Tref, RegionTable)
fHz = NaN;
rid = mode(T.region_id);
if istable(RegionTable) && ismember('peak_freq_hz', RegionTable.Properties.VariableNames)
    m = RegionTable.region_id == rid;
    if any(m)
        fHz = RegionTable.peak_freq_hz(find(m, 1));
    end
end
if ~isfinite(fHz)
    tt = Tref.time_s(:);
    h = Tref.strain_basis_real_microstrain(:);
    valid = isfinite(tt) & isfinite(h);
    tt = tt(valid);
    h = h(valid) - mean(h(valid), 'omitnan');
    if numel(tt) > 16
        dt = median(diff(tt), 'omitnan');
        fs = 1 / dt;
        n = numel(h);
        Y = abs(fft(h));
        f = (0:n-1).' * fs / n;
        keep = f > 100 & f < fs/2;
        [~, ix] = max(Y(keep));
        ff = f(keep);
        fHz = ff(ix);
    end
end
if ~isfinite(fHz)
    fHz = 580;
end
end


function [score, rmse, wrmse, layer, phaseCov, Aref, resid] = ...
    score_constrained_candidate_step06k_local(T, hAtBtt, hAtRef, C, cfg)
y = T.displacement_mm(:);
u = real(C .* hAtBtt(:));
regionIds = unique(T.region_id(:)).';
Xoff = zeros(height(T), numel(regionIds));
for i = 1:numel(regionIds)
    Xoff(:, i) = double(T.region_id == regionIds(i));
end
w = build_quality_weights_step06k_local(T, cfg);
sw = sqrt(w);
B = (Xoff .* sw) \ ((y - u) .* sw);
resid = y - (u + Xoff * B);
rmse = sqrt(mean(resid.^2, 'omitnan'));
wrmse = sqrt(sum(w .* resid.^2) / sum(w));
layer = sensor_layering_score_step06k_local(T, resid);
phaseCov = phase_coverage_step06k_local(angle(C .* hAtBtt)) / (2*pi);
uRef = real(C .* hAtRef);
Aref = sqrt(2) * std(uRef, 'omitnan');
score = wrmse + cfg.step06k_sensor_layering_penalty * layer;
end


function Fit = build_fixed_complex_fit_from_basis_step06k_local( ...
    T, Tref, hAtBtt, hAtRef, C, cfg, modelName)
y = T.displacement_mm(:);
u = real(C .* hAtBtt(:));
regionIds = unique(T.region_id(:)).';
Xoff = zeros(height(T), numel(regionIds));
for i = 1:numel(regionIds)
    Xoff(:, i) = double(T.region_id == regionIds(i));
end
w = build_quality_weights_step06k_local(T, cfg);
sw = sqrt(w);
B = (Xoff .* sw) \ ((y - u) .* sw);
off = Xoff * B;
yFit = u + off;
resid = y - yFit;
uRef = real(C .* hAtRef(:));
beta = [real(C); imag(C); B(:)];
meta = struct('regionIds', regionIds, ...
    'regionCols', 3:(2+numel(regionIds)));
Fit = make_fit_struct_step06k_local(modelName, T, Tref, beta, C, ...
    yFit, resid, off, w, w, meta, uRef, cfg);
end


function score = sensor_layering_score_step06k_local(T, resid)
sids = unique(T.sensor_id(:)).';
med = nan(numel(sids), 1);
for i = 1:numel(sids)
    med(i) = median(resid(T.sensor_id == sids(i)), 'omitnan');
end
score = max(abs(med - median(med, 'omitnan')), [], 'omitnan');
if ~isfinite(score)
    score = 0;
end
end


function [Fit, FEReference, FEWidthSensitivity] = build_fe_reference_step06k_local( ...
    cfg, Tref, phaseSourceC)
h = Tref.strain_basis_real_microstrain + 1i * Tref.strain_basis_imag_microstrain;
FEWidthSensitivity = build_fe_width_sensitivity_step06k_local(cfg, Tref);
Kmm = cfg.step06k_fe_strain_to_mm;
if cfg.step06k_fe_use_width_uncertainty_as_reference && ...
        ~isempty(FEWidthSensitivity) && ...
        ismember('in_nominal_uncertainty_band', FEWidthSensitivity.Properties.VariableNames)
    mBand = FEWidthSensitivity.in_nominal_uncertainty_band & ...
        isfinite(FEWidthSensitivity.K_mm_per_microstrain);
    if any(mBand)
        Kmm = median(FEWidthSensitivity.K_mm_per_microstrain(mBand), ...
            'omitnan');
    end
end
C = Kmm * exp(1i * angle(phaseSourceC));
uRef = real(C .* h);

Fit = struct();
Fit.model_name = "diagnostic_FE_mode_shape_K_main_phase";
Fit.C_complex_mm_per_microstrain = C;
Fit.K_mm_per_microstrain = Kmm;
Fit.K_um_per_microstrain = 1000 * Kmm;
Fit.K_ratio_to_fixed = Kmm / cfg.step06i_fixed_strain_to_mm;
Fit.phase_lag_rad = angle(C);
Fit.phase_lag_deg = angle(C) * 180 / pi;
Fit.rmse_mm = NaN;
Fit.weighted_rmse_mm = NaN;
Fit.mean_abs_residual_mm = NaN;
Fit.phase_coverage_fraction = NaN;
Fit.point_count = 0;
Fit.region_count = 0;
Fit.region_ids = [];
Fit.region_offsets_mm = [];
Fit.u_ref_mm = uRef;
Fit.u_rms_mm = std(uRef, 'omitnan');
Fit.A_peak_equiv_mm = sqrt(2) * Fit.u_rms_mm;
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - ...
    min(uRef, [], 'omitnan'));
Fit.source = string(cfg.step06k_fe_source);
Fit.transmissibility_per_m = cfg.step06k_fe_transmissibility_per_m;

FEReference = table();
FEReference.source = string(cfg.step06k_fe_source);
FEReference.point_label = string(cfg.step06k_fe_point_label);
FEReference.z_positions_mm = string(mat2str(cfg.step06k_fe_z_positions_mm));
FEReference.nominal_z_mm = cfg.step06k_fe_nominal_z_mm;
FEReference.z_uncertainty_mm = cfg.step06k_fe_z_uncertainty_mm;
FEReference.nearest_node_transmissibility_per_m = ...
    string(mat2str(cfg.step06k_fe_nearest_transmissibility_per_m));
FEReference.surface_interpolated_transmissibility_per_m = ...
    string(mat2str(cfg.step06k_fe_interpolated_transmissibility_per_m));
FEReference.transmissibility_strain_per_m = cfg.step06k_fe_transmissibility_per_m;
FEReference.K_mm_per_microstrain = Kmm;
FEReference.K_um_per_microstrain = 1000 * Kmm;
if ~isempty(FEWidthSensitivity) && ...
        ismember('in_nominal_uncertainty_band', FEWidthSensitivity.Properties.VariableNames)
    mBand = FEWidthSensitivity.in_nominal_uncertainty_band;
    FEReference.width_sweep_file = string(cfg.step06k_fe_width_sensitivity_file);
    FEReference.width_sweep_point_count = height(FEWidthSensitivity);
    FEReference.band_point_count = nnz(mBand);
    FEReference.band_K_um_per_microstrain_min = ...
        min(FEWidthSensitivity.K_um_per_microstrain(mBand), [], 'omitnan');
    FEReference.band_K_um_per_microstrain_median = ...
        median(FEWidthSensitivity.K_um_per_microstrain(mBand), 'omitnan');
    FEReference.band_K_um_per_microstrain_max = ...
        max(FEWidthSensitivity.K_um_per_microstrain(mBand), [], 'omitnan');
    FEReference.band_A_peak_equiv_mm_min = ...
        min(FEWidthSensitivity.target_A_peak_equiv_mm(mBand), [], 'omitnan');
    FEReference.band_A_peak_equiv_mm_median = ...
        median(FEWidthSensitivity.target_A_peak_equiv_mm(mBand), 'omitnan');
    FEReference.band_A_peak_equiv_mm_max = ...
        max(FEWidthSensitivity.target_A_peak_equiv_mm(mBand), [], 'omitnan');
else
    FEReference.width_sweep_file = "";
    FEReference.width_sweep_point_count = 0;
    FEReference.band_point_count = 0;
    FEReference.band_K_um_per_microstrain_min = NaN;
    FEReference.band_K_um_per_microstrain_median = NaN;
    FEReference.band_K_um_per_microstrain_max = NaN;
    FEReference.band_A_peak_equiv_mm_min = NaN;
    FEReference.band_A_peak_equiv_mm_median = NaN;
    FEReference.band_A_peak_equiv_mm_max = NaN;
end
FEReference.phase_lag_source = "target_EO_BTT_phase";
FEReference.phase_lag_rad = Fit.phase_lag_rad;
FEReference.target_A_peak_equiv_mm = Fit.A_peak_equiv_mm;
FEReference.target_A_half_range_mm = Fit.A_half_range_mm;
FEReference.note = "FE gives independent mode-shape magnitude. Width-position uncertainty is reported as sensitivity; phase is copied from the EO14 BTT-strain fit for time-series comparison.";
end


function FEWidthSensitivity = build_fe_width_sensitivity_step06k_local(cfg, Tref)
h = Tref.strain_basis_real_microstrain + 1i * Tref.strain_basis_imag_microstrain;
strainReal = real(h);
strainRms = std(strainReal, 'omitnan');
strainA = sqrt(2) * strainRms;
FEWidthSensitivity = table();
if isfield(cfg, 'step06k_fe_width_sensitivity_file') && ...
        isfile(cfg.step06k_fe_width_sensitivity_file)
    txt = fileread(cfg.step06k_fe_width_sensitivity_file);
    raw = jsondecode(txt);
    if isstruct(raw)
        z = reshape([raw.z_mm], [], 1);
        kidw = reshape([raw.k_idw_1_per_m], [], 1);
        knear = reshape([raw.k_nearest_1_per_m], [], 1);
    else
        z = [];
        kidw = [];
        knear = [];
    end
else
    z = cfg.step06k_fe_z_positions_mm(:);
    kidw = cfg.step06k_fe_interpolated_transmissibility_per_m(:);
    knear = cfg.step06k_fe_nearest_transmissibility_per_m(:);
end
if isempty(z)
    return;
end
[z, ord] = sort(z(:));
kidw = kidw(ord);
knear = knear(ord);
Kmm = 1e-3 ./ kidw;
Kum = 1000 .* Kmm;
nominalZ = cfg.step06k_fe_nominal_z_mm;
zUnc = cfg.step06k_fe_z_uncertainty_mm;
inBand = abs(z - nominalZ) <= zUnc;
targetA = Kmm .* strainA;
targetRms = Kmm .* strainRms;
FEWidthSensitivity = table();
FEWidthSensitivity.z_mm = z;
FEWidthSensitivity.k_idw_1_per_m = kidw;
FEWidthSensitivity.k_nearest_1_per_m = knear;
FEWidthSensitivity.K_mm_per_microstrain = Kmm;
FEWidthSensitivity.K_um_per_microstrain = Kum;
FEWidthSensitivity.target_A_peak_equiv_mm = targetA;
FEWidthSensitivity.target_u_rms_mm = targetRms;
FEWidthSensitivity.nominal_z_mm = repmat(nominalZ, numel(z), 1);
FEWidthSensitivity.z_uncertainty_mm = repmat(zUnc, numel(z), 1);
FEWidthSensitivity.in_nominal_uncertainty_band = inBand;
FEWidthSensitivity.source_file = repmat(string( ...
    cfg.step06k_fe_width_sensitivity_file), numel(z), 1);
end


function Fit = fit_btt_harmonic_step06k_local(T, Tref, oprTimes, eo, cfg, modelName)
theta = map_time_to_rotor_phase_step06k_local(oprTimes, T.arrival_time_s, ...
    cfg.opr_events_per_revolution);
phi = eo .* theta;
y = T.displacement_mm(:);
regionIds = unique(T.region_id(:)).';
X = [cos(phi), sin(phi)];
for i = 1:numel(regionIds)
    X(:, end+1) = double(T.region_id == regionIds(i)); %#ok<AGROW>
end
w = build_quality_weights_step06k_local(T, cfg);
[beta, weights] = weighted_irls_step06k_local(X, y, w, cfg);
a = beta(1);
b = beta(2);
off = zeros(height(T), 1);
for i = 1:numel(regionIds)
    m = T.region_id == regionIds(i);
    off(m) = beta(2+i);
end
yFit = X * beta;
resid = y - yFit;

thetaRef = map_time_to_rotor_phase_step06k_local(oprTimes, Tref.time_s, ...
    cfg.opr_events_per_revolution);
phiRef = eo .* thetaRef;
uRef = a .* cos(phiRef) + b .* sin(phiRef);

Fit = struct();
Fit.model_name = string(modelName);
Fit.eo = eo;
Fit.beta = beta;
Fit.a_cos_mm = a;
Fit.b_sin_mm = b;
Fit.K_mm_per_microstrain = NaN;
Fit.K_um_per_microstrain = NaN;
Fit.K_ratio_to_fixed = NaN;
Fit.phase_lag_rad = atan2(b, a);
Fit.phase_lag_deg = Fit.phase_lag_rad * 180 / pi;
Fit.A_peak_equiv_mm = hypot(a, b);
Fit.u_rms_mm = std(uRef, 'omitnan');
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - ...
    min(uRef, [], 'omitnan'));
Fit.time_s = Tref.time_s;
Fit.u_ref_mm = uRef;
Fit.y_fit_mm = yFit;
Fit.residual_mm = resid;
Fit.offset_at_points_mm = off;
Fit.u_at_points_mm = yFit - off;
Fit.displacement_dynamic_mm = y - off;
Fit.weights = weights;
Fit.base_weights = w;
Fit.rmse_mm = sqrt(mean(resid.^2, 'omitnan'));
Fit.weighted_rmse_mm = sqrt(sum(weights .* resid.^2) / sum(weights));
Fit.mean_abs_residual_mm = mean(abs(resid), 'omitnan');
Fit.point_count = height(T);
Fit.region_ids = regionIds;
Fit.region_offsets_mm = beta(3:end);
Fit.region_count = numel(regionIds);
Fit.phase_coverage_fraction = phase_coverage_step06k_local(phi) / (2*pi);
Fit.rank_X = rank(X);
Fit.cond_X = cond(X);
end


function Fit = build_async_amp_scaled_fit_step06k_local( ...
    TrefTarget, oprTimes, BttOnlyFit, AsyncRegionFits, targetEO, cfg)
ok = strcmp(string(AsyncRegionFits.status), "ok") & ...
    isfinite(AsyncRegionFits.K_amp_um_per_microstrain);
if any(ok)
    Kmm = median(AsyncRegionFits.K_amp_um_per_microstrain(ok), ...
        'omitnan') / 1000;
else
    [AstrainTarget, ~] = strain_region_amplitude_step06k_local(TrefTarget);
    Kmm = BttOnlyFit.A_peak_equiv_mm / max(AstrainTarget, eps);
end
[AstrainTarget, ~] = strain_region_amplitude_step06k_local(TrefTarget);
targetAmp = Kmm * AstrainTarget;
scale = targetAmp / max(BttOnlyFit.A_peak_equiv_mm, eps);
uRef = scale .* BttOnlyFit.u_ref_mm;

Fit = struct();
Fit.model_name = "diagnostic_async_amp_ratio_BTT_phase";
Fit.K_mm_per_microstrain = Kmm;
Fit.K_um_per_microstrain = 1000 * Kmm;
Fit.K_ratio_to_fixed = Kmm / cfg.step06i_fixed_strain_to_mm;
Fit.phase_lag_rad = NaN;
Fit.phase_lag_deg = NaN;
Fit.fit_rmse_mm = BttOnlyFit.rmse_mm;
Fit.rmse_mm = BttOnlyFit.rmse_mm;
Fit.weighted_rmse_mm = BttOnlyFit.weighted_rmse_mm;
Fit.phase_coverage_fraction = BttOnlyFit.phase_coverage_fraction;
Fit.point_count = BttOnlyFit.point_count;
Fit.region_count = BttOnlyFit.region_count;
Fit.region_ids = BttOnlyFit.region_ids;
Fit.u_ref_mm = uRef;
Fit.u_rms_mm = std(uRef, 'omitnan');
Fit.A_peak_equiv_mm = sqrt(2) * Fit.u_rms_mm;
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - ...
    min(uRef, [], 'omitnan'));
Fit.target_strain_A_peak_equiv_microstrain = AstrainTarget;
Fit.source = "BTT harmonic phase + asynchronous amplitude K";
Fit.eo = targetEO;
end


function w = build_quality_weights_step06k_local(T, cfg)
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
    thetaSigmaDeg = min_finite_step06k_local(thetaSigmaDeg, iqrSigma);
end
thetaSigmaDeg(~isfinite(thetaSigmaDeg)) = median(thetaSigmaDeg, 'omitnan');
if all(~isfinite(thetaSigmaDeg))
    thetaSigmaDeg = zeros(n, 1);
end
sigmaAngleMm = abs(thetaSigmaDeg) * pi / 180 * cfg.r_tip_mm;
sigmaMm = sqrt(cfg.step06k_weight_sigma_floor_mm.^2 + sigmaAngleMm.^2);
w = w .* (1 ./ max(sigmaMm.^2, eps));
w = w ./ median(w, 'omitnan');

if ismember('step02_best_corr', T.Properties.VariableNames)
    corrW = clamp_step06k_local((T.step02_best_corr(:) - 0.85) / 0.10, ...
        0.20, 1.0);
    w = w .* corrW;
end
if ismember('step02_corr_gap', T.Properties.VariableNames)
    gapW = clamp_step06k_local(T.step02_corr_gap(:) / 0.20, 0.20, 1.0);
    w = w .* gapW;
end
if ismember('step02_shift_consistency', T.Properties.VariableNames)
    shiftW = clamp_step06k_local(T.step02_shift_consistency(:), 0.20, 1.0);
    w = w .* shiftW;
end
w(~isfinite(w) | w <= 0) = 1;
w = clamp_step06k_local(w ./ median(w, 'omitnan'), 0.10, 5.0);
end


function [beta, weights] = weighted_irls_step06k_local(X, y, baseWeights, cfg)
valid = all(isfinite(X), 2) & isfinite(y) & ...
    isfinite(baseWeights) & baseWeights > 0;
Xv = X(valid, :);
yv = y(valid);
baseW = baseWeights(valid);
robW = ones(size(yv));
beta = zeros(size(X, 2), 1);
for iter = 1:cfg.step06i_irls_iterations
    w = baseW .* robW;
    sw = sqrt(w);
    betaV = (Xv .* sw) \ (yv .* sw);
    r = yv - Xv * betaV;
    s = robust_scale_step06k_local(r);
    if ~isfinite(s) || s <= eps
        beta = betaV;
        break;
    end
    cutoff = cfg.step06i_huber_k * s;
    robW = min(1, cutoff ./ max(abs(r), eps));
    beta = betaV;
end
weights = zeros(size(y));
weights(valid) = baseW .* robW;
end


function Boot = bootstrap_primary_grid_reference_step06k_local(T, TrefEO, ...
    TrefTarget, BttFit, FEReference, RegionTable, targetEO, cfg)
rng(cfg.step06i_random_seed);
if ismember('revolution_index', T.Properties.VariableNames)
    blocks = T.region_id * 1e7 + T.revolution_index;
    blocks = unique(blocks(isfinite(blocks))).';
else
    blocks = 1:height(T);
end
if isempty(blocks)
    blocks = 1:height(T);
end
nBlocks = numel(blocks);
rows = repmat(make_boot_row_step06k_local(), cfg.step06k_bootstrap_count, 1);
for b = 1:cfg.step06k_bootstrap_count
    picked = blocks(randi(nBlocks, nBlocks, 1));
    idx = [];
    if ismember('revolution_index', T.Properties.VariableNames)
        blockId = T.region_id * 1e7 + T.revolution_index;
        for k = 1:numel(picked)
            idx = [idx; find(blockId == picked(k))]; %#ok<AGROW>
        end
    else
        idx = picked(:);
    end
    rows(b).bootstrap_id = b;
    rows(b).target_eo = targetEO;
    rows(b).block_count = nBlocks;
    rows(b).point_count = numel(idx);
    try
        F = build_fe_btt_grid_reference_step06k_local( ...
            T(idx, :), TrefTarget, BttFit, FEReference, RegionTable, cfg);
        href = TrefTarget.strain_basis_real_microstrain + ...
            1i*TrefTarget.strain_basis_imag_microstrain;
        if isfield(F, 'u_ref_mm') && numel(F.u_ref_mm) == height(TrefTarget)
            uTarget = F.u_ref_mm;
        else
            uTarget = real(F.C_complex_mm_per_microstrain .* href);
        end
        rows(b).K_mm_per_microstrain = F.K_mm_per_microstrain;
        rows(b).K_um_per_microstrain = F.K_um_per_microstrain;
        rows(b).unconstrained_K_um_per_microstrain = ...
            F.unconstrained_BTT_K_um_per_microstrain;
        rows(b).constraint_status = F.constraint_status;
        rows(b).frequency_hz = F.frequency_hz;
        rows(b).objective = F.objective;
        rows(b).sensor_layering_score_mm = F.sensor_layering_score_mm;
        rows(b).phase_lag_rad = F.phase_lag_rad;
        rows(b).phase_lag_deg = F.phase_lag_deg;
        rows(b).target_u_rms_mm = std(uTarget, 'omitnan');
        rows(b).target_A_peak_equiv_mm = sqrt(2) * rows(b).target_u_rms_mm;
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


function row = make_boot_row_step06k_local()
row = struct('bootstrap_id', NaN, 'target_eo', NaN, ...
    'block_count', NaN, 'point_count', NaN, ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'unconstrained_K_um_per_microstrain', NaN, ...
    'constraint_status', "", ...
    'frequency_hz', NaN, 'objective', NaN, ...
    'sensor_layering_score_mm', NaN, ...
    'phase_lag_rad', NaN, 'phase_lag_deg', NaN, ...
    'target_u_rms_mm', NaN, 'target_A_peak_equiv_mm', NaN, ...
    'rmse_mm', NaN, 'C_real', NaN, 'C_imag', NaN, 'status', "");
end


function TS = build_target_time_series_step06k_local(Tref, PrimaryFit, ...
    SameModeFit, TargetEOFit, TargetOnlyFit, FEFit, BttOnlyFit, ...
    AsyncAmpFit, Boot, cfg)
h = Tref.strain_basis_real_microstrain + 1i*Tref.strain_basis_imag_microstrain;
TS = table();
TS.time_s = Tref.time_s;
TS.region_id = Tref.region_id;
TS.dominant_order = Tref.dominant_order;
TS.strain_basis_real_microstrain = Tref.strain_basis_real_microstrain;
TS.strain_basis_imag_microstrain = Tref.strain_basis_imag_microstrain;
TS.strain_basis_abs_microstrain = abs(h);
if isfield(PrimaryFit, 'u_ref_mm') && numel(PrimaryFit.u_ref_mm) == height(TS)
    TS.u_ref_primary_FE_BTT_constrained_mm = PrimaryFit.u_ref_mm(:);
else
    TS.u_ref_primary_FE_BTT_constrained_mm = ...
        real(PrimaryFit.C_complex_mm_per_microstrain .* h);
end
TS.u_ref_same_mode_all_regions_mm = ...
    real(SameModeFit.C_complex_mm_per_microstrain .* h);
TS.u_ref_main_same_mode_mm = TS.u_ref_primary_FE_BTT_constrained_mm;
TS.u_ref_target_EO_all_windows_mm = real(TargetEOFit.C_complex_mm_per_microstrain .* h);
TS.u_ref_target_only_mm = real(TargetOnlyFit.C_complex_mm_per_microstrain .* h);
TS.u_ref_FE_mode_shape_mm = real(FEFit.C_complex_mm_per_microstrain .* h);
TS.u_ref_btt_only_mm = BttOnlyFit.u_ref_mm;
TS.u_ref_async_amp_ratio_mm = AsyncAmpFit.u_ref_mm;
TS.u_ref_main_EOall_mm = TS.u_ref_main_same_mode_mm;
TS.u_ref_global_allEO_mm = TS.u_ref_main_same_mode_mm;
TS.u_ref_main_EOall_mean_removed_mm = TS.u_ref_main_same_mode_mm - ...
    mean(TS.u_ref_main_same_mode_mm, 'omitnan');
if ismember('u_step05_mm', Tref.Properties.VariableNames)
    TS.u_step05_mm = Tref.u_step05_mm;
end
if ismember('u_ref_eo_mm', Tref.Properties.VariableNames)
    TS.u_step06i_with_sensor_bias_mm = Tref.u_ref_eo_mm;
end
ok = strcmp(string(Boot.status), "ok") & isfinite(Boot.C_real) & ...
    isfinite(Boot.C_imag);
if any(ok)
    Cboot = Boot.C_real(ok) + 1i*Boot.C_imag(ok);
    U = real(h(:) * transpose(Cboot(:)));
    TS.u_ref_boot_low_mm = prctile(U, 2.5, 2);
    TS.u_ref_boot_high_mm = prctile(U, 97.5, 2);
else
    TS.u_ref_boot_low_mm = nan(height(TS), 1);
    TS.u_ref_boot_high_mm = nan(height(TS), 1);
end
end


function P = attach_target_point_predictions_step06k_local(T, PrimaryFit, ...
    SameModeFit, TargetOnlyFit, TargetEOFit, FEFit, BttOnlyFit, cfg)
P = T;
h = P.strain_basis_real_microstrain + 1i*P.strain_basis_imag_microstrain;
if isfield(PrimaryFit, 'grid_t0_s') && isfield(PrimaryFit, 'frequency_hz')
    hPrimary = abs(h) .* exp(1i * 2*pi*PrimaryFit.frequency_hz .* ...
        (P.arrival_time_s - PrimaryFit.grid_t0_s));
else
    hPrimary = h;
end
P.u_ref_primary_FE_BTT_constrained_at_btt_mm = ...
    real(PrimaryFit.C_complex_mm_per_microstrain .* hPrimary);
P.u_ref_same_mode_all_regions_at_btt_mm = ...
    real(SameModeFit.C_complex_mm_per_microstrain .* h);
P.u_ref_main_same_mode_at_btt_mm = P.u_ref_primary_FE_BTT_constrained_at_btt_mm;
P.u_ref_main_EOall_at_btt_mm = P.u_ref_main_same_mode_at_btt_mm;
P.u_ref_target_EO_all_windows_at_btt_mm = ...
    real(TargetEOFit.C_complex_mm_per_microstrain .* h);
P.u_ref_target_only_at_btt_mm = real(TargetOnlyFit.C_complex_mm_per_microstrain .* h);
P.u_ref_global_allEO_at_btt_mm = P.u_ref_main_same_mode_at_btt_mm;
P.u_ref_FE_mode_shape_at_btt_mm = real(FEFit.C_complex_mm_per_microstrain .* h);
P.u_ref_btt_only_at_btt_mm = BttOnlyFit.u_at_points_mm;
P.target_offset_main_mm = point_region_offset_step06k_local(P, PrimaryFit);
P.target_offset_same_mode_mm = point_region_offset_step06k_local(P, SameModeFit);
P.target_offset_target_only_mm = point_region_offset_step06k_local(P, TargetOnlyFit);
P.target_offset_btt_only_mm = BttOnlyFit.offset_at_points_mm;
P.dynamic_BTT_main_offset_removed_mm = P.displacement_mm - P.target_offset_main_mm;
P.dynamic_BTT_btt_only_offset_removed_mm = P.displacement_mm - ...
    P.target_offset_btt_only_mm;
P.total_fit_main_EOall_mm = P.target_offset_main_mm + ...
    P.u_ref_main_EOall_at_btt_mm;
P.residual_main_EOall_mm = P.displacement_mm - P.total_fit_main_EOall_mm;
P.total_fit_target_only_mm = P.target_offset_target_only_mm + ...
    P.u_ref_target_only_at_btt_mm;
P.residual_target_only_mm = P.displacement_mm - P.total_fit_target_only_mm;
P.total_fit_btt_only_mm = P.target_offset_btt_only_mm + ...
    P.u_ref_btt_only_at_btt_mm;
P.residual_btt_only_mm = P.displacement_mm - P.total_fit_btt_only_mm;
end


function off = point_region_offset_step06k_local(P, Fit)
off = zeros(height(P), 1);
for i = 1:numel(Fit.region_ids)
    m = P.region_id == Fit.region_ids(i);
    off(m) = Fit.region_offsets_mm(i);
end
end


function S = build_sensor_residuals_step06k_local(T, Fit, cfg)
h = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
u = real(Fit.C_complex_mm_per_microstrain .* h);
off = point_region_offset_step06k_local(T, Fit);
resid = T.displacement_mm - off - u;
sids = unique(T.sensor_id(:)).';
rows = repmat(make_sensor_row_step06k_local(), numel(sids), 1);
for i = 1:numel(sids)
    sid = sids(i);
    m = T.sensor_id == sid;
    rows(i).sensor_id = sid;
    rows(i).point_count = nnz(m);
    rows(i).mean_residual_mm = mean(resid(m), 'omitnan');
    rows(i).median_residual_mm = median(resid(m), 'omitnan');
    rows(i).std_residual_mm = std(resid(m), 'omitnan');
    rows(i).mean_abs_residual_mm = mean(abs(resid(m)), 'omitnan');
    rows(i).max_abs_median_relative_to_A = ...
        abs(rows(i).median_residual_mm) / max(Fit.A_peak_equiv_mm, eps);
    rows(i).status = "diagnostic_only_no_bias_used";
end
S = struct2table(rows);
end


function row = make_sensor_row_step06k_local()
row = struct('sensor_id', NaN, 'point_count', NaN, ...
    'mean_residual_mm', NaN, 'median_residual_mm', NaN, ...
    'std_residual_mm', NaN, 'mean_abs_residual_mm', NaN, ...
    'max_abs_median_relative_to_A', NaN, ...
    'status', "");
end


function Model = build_model_comparison_step06k_local(PrimaryFit, SameModeFit, ...
    TargetEOFit, TargetOnlyFit, FEFit, BttOnlyFit, AsyncAmpFit, TS, targetEO, cfg)
fits = {PrimaryFit, TargetEOFit, SameModeFit, TargetOnlyFit, FEFit, ...
    BttOnlyFit, AsyncAmpFit};
rows = repmat(make_model_row_step06k_local(), numel(fits), 1);
for i = 1:numel(fits)
    F = fits{i};
    rows(i).model_name = F.model_name;
    rows(i).target_eo = targetEO;
    rows(i).fit_point_count = F.point_count;
    rows(i).fit_region_count = F.region_count;
    rows(i).K_mm_per_microstrain = F.K_mm_per_microstrain;
    rows(i).K_um_per_microstrain = F.K_um_per_microstrain;
    rows(i).K_ratio_to_fixed = F.K_ratio_to_fixed;
    rows(i).phase_lag_rad = F.phase_lag_rad;
    rows(i).phase_lag_deg = F.phase_lag_deg;
    rows(i).fit_rmse_mm = F.rmse_mm;
    rows(i).fit_weighted_rmse_mm = F.weighted_rmse_mm;
    rows(i).fit_phase_coverage_fraction = F.phase_coverage_fraction;
    if isfield(F, 'unconstrained_BTT_K_um_per_microstrain')
        rows(i).unconstrained_BTT_K_um_per_microstrain = ...
            F.unconstrained_BTT_K_um_per_microstrain;
    end
    if isfield(F, 'FE_band_K_um_per_microstrain_min')
        rows(i).FE_band_K_um_per_microstrain_min = ...
            F.FE_band_K_um_per_microstrain_min;
        rows(i).FE_band_K_um_per_microstrain_max = ...
            F.FE_band_K_um_per_microstrain_max;
    end
    if isfield(F, 'constraint_status')
        rows(i).constraint_status = F.constraint_status;
    end
    if isfield(F, 'frequency_hz')
        rows(i).frequency_hz = F.frequency_hz;
        rows(i).frequency_prior_center_hz = F.frequency_prior_center_hz;
        rows(i).frequency_prior_half_width_hz = F.frequency_prior_half_width_hz;
        rows(i).objective = F.objective;
        rows(i).sensor_layering_score_mm = F.sensor_layering_score_mm;
    end
    u = target_series_for_model_step06k_local(TS, F.model_name);
    rows(i).target_u_rms_mm = std(u, 'omitnan');
    rows(i).target_A_peak_equiv_mm = sqrt(2) * rows(i).target_u_rms_mm;
    rows(i).target_A_half_range_mm = 0.5 * ...
        (max(u, [], 'omitnan') - min(u, [], 'omitnan'));
end
Model = struct2table(rows);
end


function u = target_series_for_model_step06k_local(TS, modelName)
    switch string(modelName)
    case {"primary_FE_prior_frequency_grid_BTT_reference", ...
            "primary_FE_band_constrained_EO14_BTT"}
        u = TS.u_ref_primary_FE_BTT_constrained_mm;
    case {"diagnostic_same_mode_all_regions", "main_same_mode_all_regions"}
        u = TS.u_ref_same_mode_all_regions_mm;
    case "diagnostic_target_EO_all_windows"
        u = TS.u_ref_target_EO_all_windows_mm;
    case "diagnostic_target_window_only"
        u = TS.u_ref_target_only_mm;
    case "diagnostic_FE_mode_shape_K_main_phase"
        u = TS.u_ref_FE_mode_shape_mm;
    case "diagnostic_BTT_only_target_window"
        u = TS.u_ref_btt_only_mm;
    case "diagnostic_async_amp_ratio_BTT_phase"
        u = TS.u_ref_async_amp_ratio_mm;
    otherwise
        u = nan(height(TS), 1);
end
end


function row = make_model_row_step06k_local()
row = struct('model_name', "", 'target_eo', NaN, ...
    'fit_point_count', NaN, 'fit_region_count', NaN, ...
    'K_mm_per_microstrain', NaN, 'K_um_per_microstrain', NaN, ...
    'K_ratio_to_fixed', NaN, 'phase_lag_rad', NaN, ...
    'phase_lag_deg', NaN, 'fit_rmse_mm', NaN, ...
    'fit_weighted_rmse_mm', NaN, ...
    'fit_phase_coverage_fraction', NaN, ...
    'unconstrained_BTT_K_um_per_microstrain', NaN, ...
    'FE_band_K_um_per_microstrain_min', NaN, ...
    'FE_band_K_um_per_microstrain_max', NaN, ...
    'constraint_status', "", ...
    'frequency_hz', NaN, 'frequency_prior_center_hz', NaN, ...
    'frequency_prior_half_width_hz', NaN, ...
    'objective', NaN, 'sensor_layering_score_mm', NaN, ...
    'target_u_rms_mm', NaN, 'target_A_peak_equiv_mm', NaN, ...
    'target_A_half_range_mm', NaN);
end


function Summary = build_summary_step06k_local(caseName, targetEO, ...
    RegionTable, RegionFits, EOFits, AsyncRegionFits, Model, Boot, TS, TargetPoints, ...
    SensorResiduals, FEReference, FEWidthSensitivity, S06I, cfg)
main = Model(Model.model_name == "primary_FE_prior_frequency_grid_BTT_reference", :);
sameMode = Model(Model.model_name == "diagnostic_same_mode_all_regions", :);
targetEOAll = Model(Model.model_name == "diagnostic_target_EO_all_windows", :);
targetOnly = Model(Model.model_name == "diagnostic_target_window_only", :);
feMode = Model(Model.model_name == "diagnostic_FE_mode_shape_K_main_phase", :);
bttOnly = Model(Model.model_name == "diagnostic_BTT_only_target_window", :);
asyncAmp = Model(Model.model_name == "diagnostic_async_amp_ratio_BTT_phase", :);
okBoot = strcmp(string(Boot.status), "ok") & isfinite(Boot.K_um_per_microstrain);
Kci = [NaN NaN];
Aci = [NaN NaN];
if any(okBoot)
    Kci = prctile(Boot.K_um_per_microstrain(okBoot), [2.5 97.5]);
    Aci = prctile(Boot.target_A_peak_equiv_mm(okBoot), [2.5 97.5]);
end

targetEORow = EOFits(EOFits.dominant_order == targetEO, :);
regionK = RegionFits.K_um_per_microstrain(strcmp(string(RegionFits.status), "ok"));
eoK = EOFits.K_um_per_microstrain(strcmp(string(EOFits.status), "ok"));

status = "accepted_FE_prior_frequency_grid_BTT_reference";
if main.fit_phase_coverage_fraction < 0.50
    status = "diagnostic_low_phase_coverage";
elseif main.target_A_peak_equiv_mm <= 0
    status = "diagnostic_invalid_amplitude";
elseif height(SensorResiduals) > 0 && ...
        max(abs(SensorResiduals.median_residual_mm), [], 'omitnan') > ...
        0.25 * sameMode.target_A_peak_equiv_mm
    status = "accepted_primary_but_same_mode_diagnostic_sensor_layering";
end

Summary = table();
Summary.case_name = string(caseName);
Summary.reference_definition = ...
    "FE_amplitude_prior_and_strain_frequency_grid_constrained_by_BTT";
Summary.reference_status = status;
Summary.reference_uses_step05 = false;
Summary.step05_usage = "comparison_only";
Summary.target_time_start_s = cfg.step06i_target_time_range_s(1);
Summary.target_time_end_s = cfg.step06i_target_time_range_s(2);
Summary.target_blade = cfg.step06i_target_blade;
Summary.target_eo = targetEO;
Summary.target_region_ids = string(mat2str(unique(TargetPoints.region_id(:)).'));
Summary.target_eo_fit_region_ids = string(targetEORow.region_ids(1));
Summary.main_fit_region_ids = string(targetEORow.region_ids(1));
Summary.same_frequency_region_count = height(RegionTable);
Summary.same_frequency_eo_count = height(EOFits);
Summary.target_window_point_count = height(TargetPoints);
Summary.target_eo_all_window_point_count = targetEOAll.fit_point_count;
Summary.same_mode_all_region_point_count = sameMode.fit_point_count;
Summary.same_mode_phase_coverage_all_regions = sameMode.fit_phase_coverage_fraction;
Summary.target_eo_phase_coverage_all_windows = targetEOAll.fit_phase_coverage_fraction;
Summary.target_window_only_phase_coverage = ...
    targetOnly.fit_phase_coverage_fraction;
Summary.K_um_per_microstrain = main.K_um_per_microstrain;
Summary.primary_unconstrained_EO14_BTT_K_um_per_microstrain = ...
    main.unconstrained_BTT_K_um_per_microstrain;
Summary.primary_constraint_status = string(main.constraint_status(1));
Summary.primary_frequency_hz = main.frequency_hz;
Summary.primary_frequency_prior_center_hz = main.frequency_prior_center_hz;
Summary.primary_frequency_prior_half_width_hz = ...
    main.frequency_prior_half_width_hz;
Summary.primary_objective = main.objective;
Summary.primary_sensor_layering_score_mm = main.sensor_layering_score_mm;
Summary.K_um_per_microstrain_ci_low = Kci(1);
Summary.K_um_per_microstrain_ci_high = Kci(2);
Summary.K_ratio_to_fixed = main.K_ratio_to_fixed;
Summary.phase_lag_rad = main.phase_lag_rad;
Summary.phase_lag_deg = main.phase_lag_deg;
Summary.target_u_rms_mm = main.target_u_rms_mm;
Summary.target_A_peak_equiv_mm = main.target_A_peak_equiv_mm;
Summary.target_A_peak_equiv_ci_low_mm = Aci(1);
Summary.target_A_peak_equiv_ci_high_mm = Aci(2);
Summary.target_A_half_range_mm = main.target_A_half_range_mm;
Summary.fit_rmse_mm = main.fit_rmse_mm;
Summary.target_only_K_um_per_microstrain = targetOnly.K_um_per_microstrain;
Summary.target_only_A_peak_equiv_mm = targetOnly.target_A_peak_equiv_mm;
Summary.target_EO_all_windows_K_um_per_microstrain = targetEOAll.K_um_per_microstrain;
Summary.target_EO_all_windows_A_peak_equiv_mm = targetEOAll.target_A_peak_equiv_mm;
Summary.same_mode_all_region_K_um_per_microstrain = sameMode.K_um_per_microstrain;
Summary.same_mode_all_region_A_peak_equiv_mm = sameMode.target_A_peak_equiv_mm;
Summary.global_allEO_K_um_per_microstrain = sameMode.K_um_per_microstrain;
Summary.global_allEO_A_peak_equiv_mm = sameMode.target_A_peak_equiv_mm;
Summary.FE_K_um_per_microstrain = feMode.K_um_per_microstrain;
Summary.FE_transmissibility_strain_per_m = cfg.step06k_fe_transmissibility_per_m;
Summary.FE_A_peak_equiv_mm = feMode.target_A_peak_equiv_mm;
if ~isempty(FEReference)
    Summary.FE_nominal_z_mm = FEReference.nominal_z_mm(1);
    Summary.FE_z_uncertainty_mm = FEReference.z_uncertainty_mm(1);
    Summary.FE_band_K_um_per_microstrain_min = ...
        FEReference.band_K_um_per_microstrain_min(1);
    Summary.FE_band_K_um_per_microstrain_median = ...
        FEReference.band_K_um_per_microstrain_median(1);
    Summary.FE_band_K_um_per_microstrain_max = ...
        FEReference.band_K_um_per_microstrain_max(1);
    Summary.FE_band_A_peak_equiv_mm_min = ...
        FEReference.band_A_peak_equiv_mm_min(1);
    Summary.FE_band_A_peak_equiv_mm_median = ...
        FEReference.band_A_peak_equiv_mm_median(1);
    Summary.FE_band_A_peak_equiv_mm_max = ...
        FEReference.band_A_peak_equiv_mm_max(1);
else
    Summary.FE_nominal_z_mm = NaN;
    Summary.FE_z_uncertainty_mm = NaN;
    Summary.FE_band_K_um_per_microstrain_min = NaN;
    Summary.FE_band_K_um_per_microstrain_median = NaN;
    Summary.FE_band_K_um_per_microstrain_max = NaN;
    Summary.FE_band_A_peak_equiv_mm_min = NaN;
    Summary.FE_band_A_peak_equiv_mm_median = NaN;
    Summary.FE_band_A_peak_equiv_mm_max = NaN;
end
if ~isempty(FEWidthSensitivity)
    Summary.FE_width_sweep_K_um_per_microstrain_min = ...
        min(FEWidthSensitivity.K_um_per_microstrain, [], 'omitnan');
    Summary.FE_width_sweep_K_um_per_microstrain_max = ...
        max(FEWidthSensitivity.K_um_per_microstrain, [], 'omitnan');
else
    Summary.FE_width_sweep_K_um_per_microstrain_min = NaN;
    Summary.FE_width_sweep_K_um_per_microstrain_max = NaN;
end
Summary.btt_only_target_A_peak_equiv_mm = bttOnly.target_A_peak_equiv_mm;
Summary.async_amp_ratio_K_um_per_microstrain = asyncAmp.K_um_per_microstrain;
Summary.async_amp_ratio_A_peak_equiv_mm = asyncAmp.target_A_peak_equiv_mm;
asyncOK = strcmp(string(AsyncRegionFits.status), "ok");
if any(asyncOK)
    Summary.async_region_K_amp_median_um_per_microstrain = ...
        median(AsyncRegionFits.K_amp_um_per_microstrain(asyncOK), 'omitnan');
    Summary.async_region_K_amp_mad_um_per_microstrain = ...
        robust_scale_step06k_local(AsyncRegionFits.K_amp_um_per_microstrain(asyncOK));
else
    Summary.async_region_K_amp_median_um_per_microstrain = NaN;
    Summary.async_region_K_amp_mad_um_per_microstrain = NaN;
end
    Summary.region_K_median_um_per_microstrain = median(regionK, 'omitnan');
    Summary.region_K_mad_um_per_microstrain = robust_scale_step06k_local(regionK);
    Summary.eo_K_median_um_per_microstrain = median(eoK, 'omitnan');
    Summary.eo_K_mad_um_per_microstrain = robust_scale_step06k_local(eoK);
targetRegionRows = RegionFits(RegionFits.dominant_order == targetEO & ...
    strcmp(string(RegionFits.status), "ok"), :);
if ~isempty(targetRegionRows)
    Summary.target_eo_region_K_min_um_per_microstrain = ...
        min(targetRegionRows.K_um_per_microstrain, [], 'omitnan');
    Summary.target_eo_region_K_max_um_per_microstrain = ...
        max(targetRegionRows.K_um_per_microstrain, [], 'omitnan');
    Summary.target_eo_region_projected_A_min_mm = ...
        min(targetRegionRows.target_projection_A_peak_equiv_mm, [], 'omitnan');
    Summary.target_eo_region_projected_A_max_mm = ...
        max(targetRegionRows.target_projection_A_peak_equiv_mm, [], 'omitnan');
else
    Summary.target_eo_region_K_min_um_per_microstrain = NaN;
    Summary.target_eo_region_K_max_um_per_microstrain = NaN;
    Summary.target_eo_region_projected_A_min_mm = NaN;
    Summary.target_eo_region_projected_A_max_mm = NaN;
end
    Summary.max_abs_sensor_median_residual_mm = ...
        max(abs(SensorResiduals.median_residual_mm), [], 'omitnan');
if ismember('u_step05_mm', TS.Properties.VariableNames)
    mask = isfinite(TS.u_step05_mm) & isfinite(TS.u_ref_main_same_mode_mm);
    Summary.step05_time_rmse_mm = sqrt(mean((TS.u_step05_mm(mask) - ...
        TS.u_ref_main_same_mode_mm(mask)).^2, 'omitnan'));
    Summary.step05_time_corr = corr(TS.u_step05_mm(mask), ...
        TS.u_ref_main_same_mode_mm(mask), 'Rows', 'complete');
    Summary.step05_A_peak_equiv_mm = sqrt(2) * ...
        std(TS.u_step05_mm(mask), 'omitnan');
    Summary.step05_minus_reference_A_percent = 100 * ...
        (Summary.step05_A_peak_equiv_mm - Summary.target_A_peak_equiv_mm) / ...
        Summary.target_A_peak_equiv_mm;
else
    Summary.step05_time_rmse_mm = NaN;
    Summary.step05_time_corr = NaN;
    Summary.step05_A_peak_equiv_mm = NaN;
    Summary.step05_minus_reference_A_percent = NaN;
end
if isfield(S06I, 'Summary') && istable(S06I.Summary)
    Summary.step06i_target_EO_K_um_per_microstrain = ...
        1000 * S06I.Summary.target_EO_K_mm_per_microstrain(1);
    Summary.step06i_target_EO_A_peak_equiv_mm = ...
        S06I.Summary.target_EO_A_rms_mm(1);
else
    Summary.step06i_target_EO_K_um_per_microstrain = NaN;
    Summary.step06i_target_EO_A_peak_equiv_mm = NaN;
end
end


function [figFile, feFigFile] = plot_step06k_reference_local(TS, TargetPoints, SameFreqPoints, ...
    RegionTable, RegionFits, EOFits, AsyncRegionFits, Model, SensorResiduals, ...
    Boot, FEWidthSensitivity, Summary, figDir, cfg)
fig = figure('Visible', cfg.step06k_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 32, 24], ...
    'Name', 'Step06K same-frequency robust reference');
tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
green = [0.466 0.674 0.188];
purple = [0.494 0.184 0.556];
gray = [0.520 0.520 0.520];
lightBlue = [0.80 0.90 1.00];
red = [0.635 0.078 0.184];

targetEO = Summary.target_eo(1);
targetRange = [Summary.target_time_start_s(1), Summary.target_time_end_s(1)];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
eoVals = unique(RegionTable.dominant_order(:));
yLim = [min(eoVals)-0.8, max(eoVals)+0.8];
patch(ax1, [targetRange(1) targetRange(2) targetRange(2) targetRange(1)], ...
    [yLim(1) yLim(1) yLim(2) yLim(2)], lightBlue, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.75, 'HandleVisibility', 'off');
for i = 1:height(RegionTable)
    t0 = RegionTable.time_start_s(i);
    t1 = RegionTable.time_end_s(i);
    eo = RegionTable.dominant_order(i);
    fc = 0.72 * [1 1 1];
    ec = 0.45 * [1 1 1];
    if eo == targetEO
        fc = blue;
        ec = blue;
    end
    patch(ax1, [t0 t1 t1 t0], [eo-0.28 eo-0.28 eo+0.28 eo+0.28], ...
        fc, 'EdgeColor', ec, 'FaceAlpha', 0.72, 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
    text(ax1, mean([t0 t1]), eo + 0.34, sprintf('R%d', ...
        RegionTable.region_id(i)), 'HorizontalAlignment', 'center', ...
        'FontSize', 7.4, 'Color', [0.15 0.15 0.15]);
end
xlim(ax1, [min(RegionTable.time_start_s)-0.6, ...
    max(RegionTable.time_end_s)+0.6]);
ylim(ax1, yLim);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'EO');
title(ax1, '1 Same-frequency resonance windows');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
okR = strcmp(string(RegionFits.status), "ok");
xR = 1:height(RegionFits);
b = bar(ax2, xR, RegionFits.K_um_per_microstrain, 0.72);
b.DisplayName = 'per-region K';
b.FaceColor = 'flat';
b.CData = repmat(0.68*[1 1 1], height(RegionFits), 1);
b.CData(RegionFits.dominant_order == targetEO, :) = repmat(blue, ...
    nnz(RegionFits.dominant_order == targetEO), 1);
yline(ax2, Summary.FE_band_K_um_per_microstrain_min, ':', 'Color', red, ...
    'LineWidth', 0.8, 'DisplayName', 'FE K band');
yline(ax2, Summary.FE_band_K_um_per_microstrain_max, ':', 'Color', red, ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax2, Summary.FE_K_um_per_microstrain, ':', 'Color', red, ...
    'LineWidth', 1.0, 'DisplayName', 'FE K');
yline(ax2, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'primary K');
set(ax2, 'XTick', xR, 'XTickLabel', string(RegionFits.region_id));
xlabel(ax2, 'Region id');
ylabel(ax2, 'K (um/microstrain)');
title(ax2, '2 Per-region no-bias gain');
legend(ax2, 'Location', 'best');

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
xE = 1:height(EOFits);
b = bar(ax3, xE, EOFits.K_um_per_microstrain, 0.68);
b.DisplayName = 'per-EO K';
b.FaceColor = 'flat';
b.CData = repmat(0.68*[1 1 1], height(EOFits), 1);
iTarget = find(EOFits.dominant_order == targetEO, 1);
if ~isempty(iTarget)
    b.CData(iTarget, :) = purple;
end
yline(ax3, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'primary K');
yline(ax3, Summary.FE_K_um_per_microstrain, ':', 'Color', red, ...
    'LineWidth', 1.0, 'DisplayName', 'FE K');
set(ax3, 'XTick', xE, 'XTickLabel', string(EOFits.dominant_order));
xlabel(ax3, 'EO');
ylabel(ax3, 'K (um/microstrain)');
title(ax3, '3 Per-EO same-frequency gain');
legend(ax3, 'Location', 'best');

ax4 = nexttile([1 2]); hold(ax4, 'on'); box(ax4, 'on');
if all(isfinite(TS.u_ref_boot_low_mm))
    fill(ax4, [TS.time_s; flipud(TS.time_s)], ...
        [TS.u_ref_boot_low_mm; flipud(TS.u_ref_boot_high_mm)], ...
        [0.88 0.88 0.88], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.75, 'DisplayName', 'bootstrap 95% band');
end
plot_downsampled_step06k_local(ax4, TS.time_s, TS.u_ref_main_EOall_mm, ...
    blue, 'primary: FE-band constrained EO14 BTT', 1.15);
plot_downsampled_step06k_local(ax4, TS.time_s, ...
    TS.u_ref_target_EO_all_windows_mm, purple, 'EO14 BTT unconstrained', 0.85);
plot_downsampled_step06k_local(ax4, TS.time_s, ...
    TS.u_ref_same_mode_all_regions_mm, gray, 'same-mode all-region diagnostic', 0.75);
plot_downsampled_step06k_local(ax4, TS.time_s, TS.u_ref_target_only_mm, ...
    orange, 'target window only', 0.85);
plot_downsampled_step06k_local(ax4, TS.time_s, TS.u_ref_FE_mode_shape_mm, ...
    red, 'FE band median', 0.75);
xlim(ax4, targetRange);
xlabel(ax4, 'Time (s)');
ylabel(ax4, 'u_{tip} (mm)');
title(ax4, '4 Complete target-window tip displacement');
legend(ax4, 'Location', 'best');

axStep05 = nexttile; hold(axStep05, 'on'); box(axStep05, 'on');
zoomRange = cfg.step06k_zoom_time_range_s(:).';
if zoomRange(1) < targetRange(1) || zoomRange(2) > targetRange(2) || ...
        zoomRange(2) <= zoomRange(1)
    zoomRange = [targetRange(1), min(targetRange(2), targetRange(1) + 0.06)];
end
zoomMask = TS.time_s >= zoomRange(1) & TS.time_s <= zoomRange(2);
plot_downsampled_step06k_local(axStep05, TS.time_s(zoomMask), ...
    TS.u_ref_main_EOall_mm(zoomMask), ...
    blue, 'Step06K reference', 1.05);
if ismember('u_step05_mm', TS.Properties.VariableNames)
    plot_downsampled_step06k_local(axStep05, TS.time_s(zoomMask), ...
        TS.u_step05_mm(zoomMask), ...
        gray, 'Step05 comparison', 0.85);
end
if ismember('u_step06i_with_sensor_bias_mm', TS.Properties.VariableNames)
    plot_downsampled_step06k_local(axStep05, TS.time_s(zoomMask), ...
        TS.u_step06i_with_sensor_bias_mm(zoomMask), orange, ...
        'Step06I sensor-bias fit', 0.75);
end
plot_downsampled_step06k_local(axStep05, TS.time_s(zoomMask), ...
    TS.u_ref_FE_mode_shape_mm(zoomMask), purple, 'FE mode-shape K', 0.70);
xlim(axStep05, zoomRange);
xlabel(axStep05, 'Time (s)');
ylabel(axStep05, 'u_{tip} (mm)');
title(axStep05, '5 Zoomed waveform comparison');
legend(axStep05, 'Location', 'best');

ax5 = nexttile; hold(ax5, 'on'); box(ax5, 'on');
okA = strcmp(string(AsyncRegionFits.status), "ok");
bar(ax5, 1:height(AsyncRegionFits), ...
    AsyncRegionFits.K_amp_um_per_microstrain, 0.72, ...
    'FaceColor', [0.70 0.70 0.70], 'DisplayName', 'BTT amp / strain amp');
targetRows = AsyncRegionFits.dominant_order == targetEO;
if any(targetRows)
    scatter(ax5, find(targetRows), ...
        AsyncRegionFits.K_amp_um_per_microstrain(targetRows), 42, ...
        blue, 'filled', 'DisplayName', 'target EO windows');
end
yline(ax5, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'primary K');
yline(ax5, Summary.FE_K_um_per_microstrain, ':', 'Color', red, ...
    'LineWidth', 1.0, 'DisplayName', 'FE K');
set(ax5, 'XTick', 1:height(AsyncRegionFits), ...
    'XTickLabel', string(AsyncRegionFits.region_id));
xlabel(ax5, 'Region id');
ylabel(ax5, 'K_{amp} (um/microstrain)');
title(ax5, '6 Asynchronous amplitude-ratio check');
legend(ax5, 'Location', 'best');

axBoot = nexttile; hold(axBoot, 'on'); box(axBoot, 'on');
okBoot = strcmp(string(Boot.status), "ok");
if any(okBoot)
    histogram(axBoot, Boot.K_um_per_microstrain(okBoot), 18, ...
        'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    xline(axBoot, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
        'LineWidth', 1.1);
    xline(axBoot, Summary.FE_K_um_per_microstrain, ':', 'Color', red, ...
        'LineWidth', 1.0);
end
xlabel(axBoot, 'K (um/microstrain)');
ylabel(axBoot, 'Count');
title(axBoot, '7 Main EO block bootstrap');

ax6 = nexttile; hold(ax6, 'on'); box(ax6, 'on');
sids = unique(TargetPoints.sensor_id(:)).';
colors = lines(max(numel(sids), 3));
for i = 1:numel(sids)
    sid = sids(i);
    m = TargetPoints.sensor_id == sid;
    scatter(ax6, TargetPoints.arrival_time_s(m), ...
        TargetPoints.dynamic_BTT_main_offset_removed_mm(m), 10, ...
        colors(i, :), 'filled', 'MarkerFaceAlpha', 0.55, ...
        'DisplayName', sprintf('S%d BTT', sid));
end
plot_downsampled_step06k_local(ax6, TS.time_s, TS.u_ref_main_EOall_mm, ...
    blue, 'main reference', 1.2);
xlim(ax6, targetRange);
xlabel(ax6, 'Time (s)');
ylabel(ax6, 'Dynamic BTT / u_{ref} (mm)');
title(ax6, '8 BTT anchors after common offset removal');
legend(ax6, 'Location', 'best');

ax7 = nexttile; hold(ax7, 'on'); box(ax7, 'on');
mp = isfinite(TargetPoints.dynamic_BTT_main_offset_removed_mm) & ...
    isfinite(TargetPoints.u_ref_main_EOall_at_btt_mm);
scatter(ax7, TargetPoints.u_ref_main_EOall_at_btt_mm(mp), ...
    TargetPoints.dynamic_BTT_main_offset_removed_mm(mp), 13, blue, ...
    'filled', 'MarkerFaceAlpha', 0.42);
if any(mp)
    xy = [TargetPoints.u_ref_main_EOall_at_btt_mm(mp); ...
        TargetPoints.dynamic_BTT_main_offset_removed_mm(mp)];
    lim = [min(xy, [], 'omitnan'), max(xy, [], 'omitnan')];
    pad = 0.05 * max(diff(lim), eps);
    lim = lim + [-pad pad];
    plot(ax7, lim, lim, 'k--', 'LineWidth', 0.8);
    xlim(ax7, lim);
    ylim(ax7, lim);
    axis(ax7, 'square');
end
xlabel(ax7, 'u_{ref}(t_i) (mm)');
ylabel(ax7, 'BTT minus B_r (mm)');
title(ax7, '9 Target-window point consistency');

ax8 = nexttile; hold(ax8, 'on'); box(ax8, 'on');
bar(ax8, SensorResiduals.sensor_id, ...
    [SensorResiduals.median_residual_mm, SensorResiduals.mean_residual_mm], ...
    0.72);
yline(ax8, 0, 'k:', 'LineWidth', 0.8);
xlabel(ax8, 'Sensor id');
ylabel(ax8, 'Residual (mm)');
title(ax8, '10 Sensor residual diagnostic');
legend(ax8, {'median', 'mean'}, 'Location', 'best');

ax9 = nexttile; hold(ax9, 'on'); box(ax9, 'on');
modelNames = categorical(Model.model_name);
modelNames = reordercats(modelNames, string(Model.model_name));
bar(ax9, modelNames, Model.target_A_peak_equiv_mm, 0.70);
ylabel(ax9, 'Target A_{peak,eq} (mm)');
title(ax9, '11 Target amplitude sensitivity');
xtickangle(ax9, 28);

statusText = strrep(string(Summary.reference_status(1)), '_', ' ');
sg = sgtitle(fig, sprintf(['Step06K: same-frequency windows diagnose K, ', ...
    'FE band constrains EO14 BTT gain to define u_{ref}\n', ...
    'EO%d, %.2f-%.2f s, K=%.3f um/microstrain, A=%.3f mm, status=%s'], ...
    targetEO, targetRange(1), targetRange(2), ...
    Summary.K_um_per_microstrain(1), ...
    Summary.target_A_peak_equiv_mm(1), ...
    statusText));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.5);

style_step06k_figure_local(fig);
figFile = fullfile(figDir, ...
    'Step06K_SameFrequencyRobustReference_20250527.png');
save_figure_step06k_local(fig, figFile);
feFigFile = plot_step06k_fe_width_sensitivity_local(FEWidthSensitivity, ...
    RegionFits, EOFits, Model, Summary, figDir, cfg);
end


function feFigFile = plot_step06k_fe_width_sensitivity_local(FEWidthSensitivity, ...
    RegionFits, EOFits, Model, Summary, figDir, cfg)
feFigFile = "";
if isempty(FEWidthSensitivity)
    return;
end
fig = figure('Visible', cfg.step06k_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [3, 3, 24, 16], ...
    'Name', 'Step06K FE width sensitivity');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
green = [0.466 0.674 0.188];
purple = [0.494 0.184 0.556];
gray = [0.520 0.520 0.520];
red = [0.635 0.078 0.184];
lightRed = [1.00 0.88 0.86];

z = FEWidthSensitivity.z_mm;
kFe = FEWidthSensitivity.K_um_per_microstrain;
inBand = FEWidthSensitivity.in_nominal_uncertainty_band;
bandZ = [Summary.FE_nominal_z_mm(1) - Summary.FE_z_uncertainty_mm(1), ...
    Summary.FE_nominal_z_mm(1) + Summary.FE_z_uncertainty_mm(1)];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
yLim1 = [min(kFe, [], 'omitnan'), max(kFe, [], 'omitnan')];
yPad1 = max(0.08 * diff(yLim1), 0.02);
yLim1 = yLim1 + [-yPad1 yPad1];
patch(ax1, [bandZ(1) bandZ(2) bandZ(2) bandZ(1)], ...
    [yLim1(1) yLim1(1) yLim1(2) yLim1(2)], lightRed, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.75, 'DisplayName', 'Z uncertainty');
plot(ax1, z, kFe, '-o', 'Color', red, 'MarkerFaceColor', red, ...
    'LineWidth', 1.1, 'MarkerSize', 4.0, 'DisplayName', 'FE width sweep');
scatter(ax1, z(inBand), kFe(inBand), 42, blue, 'filled', ...
    'DisplayName', 'used band');
yline(ax1, Summary.FE_K_um_per_microstrain, ':', 'Color', red, ...
    'LineWidth', 1.1, 'DisplayName', 'FE band median');
yline(ax1, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'main BTT-strain K');
yline(ax1, Summary.target_EO_all_windows_K_um_per_microstrain, '--', ...
    'Color', blue, 'LineWidth', 1.0, 'DisplayName', 'EO14 all-window K');
xlabel(ax1, 'Gauge width position Z (mm)');
ylabel(ax1, 'K (um/microstrain)');
title(ax1, '1 FE width-position sensitivity');
ylim(ax1, yLim1);
legend(ax1, 'Location', 'best');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
plot(ax2, z, FEWidthSensitivity.target_A_peak_equiv_mm, '-o', ...
    'Color', red, 'MarkerFaceColor', red, 'LineWidth', 1.1, ...
    'MarkerSize', 4.0, 'DisplayName', 'FE-implied target A');
scatter(ax2, z(inBand), FEWidthSensitivity.target_A_peak_equiv_mm(inBand), ...
    42, blue, 'filled', 'DisplayName', 'used band');
yline(ax2, Summary.target_A_peak_equiv_mm, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'main reference A');
yline(ax2, Summary.target_EO_all_windows_A_peak_equiv_mm, '--', ...
    'Color', blue, 'LineWidth', 1.0, 'DisplayName', 'EO14 all-window A');
yline(ax2, Summary.step05_A_peak_equiv_mm, ':', 'Color', gray, ...
    'LineWidth', 1.0, 'DisplayName', 'Step05 comparison A');
xlabel(ax2, 'Gauge width position Z (mm)');
ylabel(ax2, 'A_{peak,eq} (mm)');
title(ax2, '2 Target amplitude implied by FE K');
legend(ax2, 'Location', 'best');

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
okR = strcmp(string(RegionFits.status), "ok");
scatter(ax3, RegionFits.dominant_order(okR), ...
    RegionFits.K_um_per_microstrain(okR), 52, ...
    RegionFits.region_id(okR), 'filled', 'MarkerFaceAlpha', 0.78, ...
    'DisplayName', 'per-region BTT-strain K');
yline(ax3, Summary.FE_band_K_um_per_microstrain_min, ':', ...
    'Color', red, 'LineWidth', 0.8, 'DisplayName', 'FE Z-band min/max');
yline(ax3, Summary.FE_band_K_um_per_microstrain_max, ':', ...
    'Color', red, 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax3, Summary.FE_K_um_per_microstrain, '-', 'Color', red, ...
    'LineWidth', 1.0, 'DisplayName', 'FE Z-band median');
yline(ax3, Summary.K_um_per_microstrain, '-', 'Color', purple, ...
    'LineWidth', 1.0, 'DisplayName', 'main same-mode K');
cb = colorbar(ax3);
cb.Label.String = 'Region id';
xlabel(ax3, 'EO');
ylabel(ax3, 'K (um/microstrain)');
title(ax3, '3 Experimental region gains vs FE band');
legend(ax3, 'Location', 'best');

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
names = categorical(Model.model_name);
names = reordercats(names, string(Model.model_name));
bar(ax4, names, Model.K_um_per_microstrain, 0.70, ...
    'FaceColor', [0.68 0.68 0.68]);
yline(ax4, Summary.FE_band_K_um_per_microstrain_min, ':', ...
    'Color', red, 'LineWidth', 0.8);
yline(ax4, Summary.FE_band_K_um_per_microstrain_max, ':', ...
    'Color', red, 'LineWidth', 0.8);
yline(ax4, Summary.FE_K_um_per_microstrain, '-', 'Color', red, ...
    'LineWidth', 1.0);
ylabel(ax4, 'K (um/microstrain)');
title(ax4, '4 K estimates used for benchmark diagnosis');
xtickangle(ax4, 28);

sg = sgtitle(fig, sprintf(['FE prior is a sensitivity band, not a single truth: ', ...
    'Z=%.2f +/- %.2f mm, K=%.3f-%.3f um/microstrain'], ...
    Summary.FE_nominal_z_mm(1), Summary.FE_z_uncertainty_mm(1), ...
    Summary.FE_band_K_um_per_microstrain_min(1), ...
    Summary.FE_band_K_um_per_microstrain_max(1)));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.5);
style_step06k_figure_local(fig);
feFigFile = fullfile(figDir, ...
    'Step06K_FEWidthSensitivity_20250527.png');
save_figure_step06k_local(fig, feFigFile);
end


function plot_downsampled_step06k_local(ax, x, y, color, name, lw)
if isempty(x) || isempty(y)
    return;
end
step = max(1, ceil(numel(x) / 30000));
plot(ax, x(1:step:end), y(1:step:end), '-', ...
    'Color', color, 'LineWidth', lw, 'DisplayName', name);
end


function style_step06k_figure_local(fig)
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    set(axs(i), 'FontName', 'Times New Roman', 'FontSize', 8.2, ...
        'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.75, ...
        'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12);
end
txt = findall(fig, 'Type', 'text');
set(txt, 'FontName', 'Times New Roman');
leg = findall(fig, 'Type', 'legend');
for i = 1:numel(leg)
    set(leg(i), 'FontName', 'Times New Roman', 'FontSize', 7.1, ...
        'Box', 'off');
end
end


function save_figure_step06k_local(fig, pngFile)
[folder, base, ~] = fileparts(pngFile);
ensure_dir_step06k_local(folder);
pngFile = fullfile(folder, [base '.png']);
figFile = fullfile(folder, [base '.fig']);
try
    exportgraphics(fig, pngFile, 'Resolution', 220);
catch
    print(fig, pngFile, '-dpng', '-r220');
end
try
    savefig(fig, figFile);
catch
end
end


function x = min_finite_step06k_local(a, b)
x = a;
m = isfinite(b) & (~isfinite(a) | b < a);
x(m) = b(m);
end


function y = clamp_step06k_local(x, lo, hi)
y = min(max(x, lo), hi);
end


function p = wrap_to_pi_step06k_local(p)
p = mod(p + pi, 2*pi) - pi;
end


function s = robust_scale_step06k_local(r)
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


function c = phase_coverage_step06k_local(phi)
phi = phi(isfinite(phi));
if numel(phi) < 2
    c = 0;
    return;
end
p = sort(mod(phi(:), 2*pi));
gaps = diff([p; p(1) + 2*pi]);
c = 2*pi - max(gaps);
end


function oprTimes = load_opr_times_step06k_local(file)
S = load(file, 'jiluOPR');
if ~isfield(S, 'jiluOPR') || isempty(S.jiluOPR)
    error('Invalid jiluOPR file:\n  %s', file);
end
oprTimes = S.jiluOPR(:, 1);
oprTimes = oprTimes(:);
oprTimes = oprTimes(isfinite(oprTimes));
if any(diff(oprTimes) <= 0)
    oprTimes = sort(oprTimes);
end
end


function thetaRot = map_time_to_rotor_phase_step06k_local( ...
    oprTimes, sampleTimes, oprEventsPerRev)
sampleTimes = sampleTimes(:);
thetaRot = nan(size(sampleTimes));
if numel(oprTimes) <= oprEventsPerRev
    return;
end
revAnchorTimes = oprTimes(1:oprEventsPerRev:end);
revAnchorPhase = 2*pi*(0:numel(revAnchorTimes)-1).';
thetaRot(:) = interp1(revAnchorTimes(:), revAnchorPhase(:), ...
    sampleTimes, 'linear', 'extrap');
thetaRot(sampleTimes < oprTimes(1) | sampleTimes > oprTimes(end)) = NaN;
end


function ensure_dir_step06k_local(d)
if ~isfolder(d)
    mkdir(d);
end
end
