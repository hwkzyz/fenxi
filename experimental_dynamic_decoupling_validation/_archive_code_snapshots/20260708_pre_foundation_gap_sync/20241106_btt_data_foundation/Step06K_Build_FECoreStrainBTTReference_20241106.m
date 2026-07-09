%% Step06K_Build_FECoreStrainBTTReference_20241106.m
% Build a reviewer-facing reference displacement for the 20241106 route.
%
% Branch 1:
%   Use the Step06A resonance region, trim the two sides, extract the
%   synchronous strain complex envelope, and convert its amplitude to
%   blade-tip displacement with an FE transfer-ratio band.
%
% Branch 2:
%   Keep the FE-strain amplitude band fixed and use sparse BTT points only
%   as a residual/phase-consistency correction with one common static
%   offset. No sensor-specific physical offset is used in the reference.

clc; clear; close all;

cfg = BTTDataConfig_20241106();
cfg = apply_step06k_defaults_20241106_local(cfg);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06K FE-core strain/BTT reference: %s ===\n', caseName);

    files = resolve_step06k_files_20241106_local(cfg, caseName);
    assert_step06k_files_20241106_local(files);

    S06A = load(files.step06a_mat, 'Evidence');
    Evidence = S06A.Evidence;
    Bundle = load(files.step03_bundle_mat, 'bundle');
    B = Bundle.bundle;
    Obs = B.Observation_Table;
    oprTimes = B.OPR_Times_s(:);

    Region = select_reference_region_20241106_local(Evidence.RegionTable, cfg);
    StrainBasis = build_region_strain_basis_20241106_local( ...
        Evidence, Region, oprTimes, cfg);
    FitPointsAll = build_btt_fit_points_20241106_local( ...
        Obs, StrainBasis, oprTimes, cfg);
    FESpectrumReference = build_fe_spectrum_reference_20241106_local( ...
        StrainBasis.Tref, cfg);
    BttCorrectionPoints = select_btt_correction_points_20241106_local( ...
        FitPointsAll, FESpectrumReference, cfg);
    [PrimaryFit, PrimaryGrid] = build_fe_btt_reference_20241106_local( ...
        BttCorrectionPoints, StrainBasis.Tref, FESpectrumReference, cfg);
    TargetTimeSeries = build_time_series_20241106_local( ...
        StrainBasis.Tref, PrimaryFit, FESpectrumReference);
    BttCorrectionPoints = attach_btt_predictions_20241106_local( ...
        BttCorrectionPoints, PrimaryFit);
    SensorResiduals = build_sensor_residuals_20241106_local( ...
        BttCorrectionPoints);
    Step05Comparison = load_step05_comparison_20241106_local(files.step05_csv);
    Summary = build_summary_20241106_local(caseName, Region, ...
        FESpectrumReference, PrimaryFit, BttCorrectionPoints, ...
        SensorResiduals, Step05Comparison, cfg);

    outDir = fullfile(cfg.step06k_output_dir, caseName);
    figDir = fullfile(cfg.step06k_figure_dir, caseName);
    ensure_dir_20241106_local(outDir);
    ensure_dir_20241106_local(figDir);

    summaryCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_Summary_20241106.csv');
    feCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_FESpectrumReference_20241106.csv');
    pointsCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_BTTCorrectionPoints_20241106.csv');
    tsCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_TargetTimeSeries_20241106.csv');
    gridCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_PrimaryGrid_20241106.csv');
    sensorCsv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_SensorResiduals_20241106.csv');
    step05Csv = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_Step05Comparison_20241106.csv');
    matFile = fullfile(outDir, ...
        'Step06K_FECoreStrainBTTReference_20241106.mat');

    writetable(Summary, summaryCsv);
    writetable(FESpectrumReference, feCsv);
    writetable(BttCorrectionPoints, pointsCsv);
    writetable(TargetTimeSeries, tsCsv);
    writetable(PrimaryGrid, gridCsv);
    writetable(SensorResiduals, sensorCsv);
    writetable(Step05Comparison, step05Csv);
    save(matFile, 'Summary', 'FESpectrumReference', 'BttCorrectionPoints', ...
        'TargetTimeSeries', 'PrimaryGrid', 'SensorResiduals', ...
        'Step05Comparison', 'PrimaryFit', 'StrainBasis', 'Region', ...
        'files', '-v7.3');

    figFile = plot_reference_20241106_local(TargetTimeSeries, ...
        BttCorrectionPoints, FESpectrumReference, SensorResiduals, ...
        Step05Comparison, Summary, figDir, cfg);

    fprintf('\nStep06K 20241106 reference summary:\n');
    disp(Summary);
    fprintf(['Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        '  %s\n  %s\n  %s\n  %s\n'], summaryCsv, feCsv, pointsCsv, ...
        tsCsv, gridCsv, sensorCsv, step05Csv, matFile, figFile);
end


function cfg = apply_step06k_defaults_20241106_local(cfg)
if ~isfield(cfg, 'step06k_target_blade') || isempty(cfg.step06k_target_blade)
    cfg.step06k_target_blade = cfg.step05_target_blades(1);
end
if ~isfield(cfg, 'step06k_sensor_ids') || isempty(cfg.step06k_sensor_ids)
    cfg.step06k_sensor_ids = cfg.step05_analysis_sensors(:).';
end
if ~isfield(cfg, 'opr_events_per_revolution') || isempty(cfg.opr_events_per_revolution)
    if isfield(cfg, 'opr_pulses_per_rev') && ~isempty(cfg.opr_pulses_per_rev)
        cfg.opr_events_per_revolution = cfg.opr_pulses_per_rev;
    else
        cfg.opr_events_per_revolution = 1;
    end
end
if ~isfield(cfg, 'step06k_output_dir') || isempty(cfg.step06k_output_dir)
    cfg.step06k_output_dir = fullfile(cfg.output_root, ...
        'step06k_fe_core_strain_btt_reference');
end
if ~isfield(cfg, 'step06k_figure_dir') || isempty(cfg.step06k_figure_dir)
    cfg.step06k_figure_dir = fullfile(cfg.figure_root, ...
        'step06k_fe_core_strain_btt_reference');
end
if ~isfield(cfg, 'step06k_fig_visible') || isempty(cfg.step06k_fig_visible)
    if usejava('desktop')
        cfg.step06k_fig_visible = 'on';
    else
        cfg.step06k_fig_visible = 'off';
    end
end
if ~isfield(cfg, 'step06k_core_trim_fraction') || isempty(cfg.step06k_core_trim_fraction)
    cfg.step06k_core_trim_fraction = 0.15;
end
if ~isfield(cfg, 'step06k_core_trim_min_s') || isempty(cfg.step06k_core_trim_min_s)
    cfg.step06k_core_trim_min_s = 0.5;
end
if ~isfield(cfg, 'step06k_core_trim_max_s') || isempty(cfg.step06k_core_trim_max_s)
    cfg.step06k_core_trim_max_s = 2.0;
end
if ~isfield(cfg, 'step06k_demod_margin_s') || isempty(cfg.step06k_demod_margin_s)
    cfg.step06k_demod_margin_s = 2.0;
end
if ~isfield(cfg, 'step06k_envelope_lowpass_hz') || isempty(cfg.step06k_envelope_lowpass_hz)
    cfg.step06k_envelope_lowpass_hz = 8.0;
end
if ~isfield(cfg, 'step06k_fe_prior_summary_csv') || ...
        isempty(cfg.step06k_fe_prior_summary_csv)
    cfg.step06k_fe_prior_summary_csv = fullfile(cfg.output_root, ...
        'step06j_fe632_transfer_prior', ...
        'Step06J_FE632_TransferPrior_Summary_20241106.csv');
end
if ~isfield(cfg, 'step06k_fe_K_um_per_microstrain_band') || ...
        isempty(cfg.step06k_fe_K_um_per_microstrain_band)
    if isfile(cfg.step06k_fe_prior_summary_csv)
        FE = readtable(cfg.step06k_fe_prior_summary_csv, ...
            'VariableNamingRule', 'preserve');
        cfg.step06k_fe_K_um_per_microstrain_band = [ ...
            FE.FE_K_min_um_per_microstrain(1), ...
            FE.FE_K_nominal_um_per_microstrain(1), ...
            FE.FE_K_max_um_per_microstrain(1)];
        if ismember('FE_source', FE.Properties.VariableNames)
            cfg.step06k_fe_source = sprintf('%s; loaded from %s', ...
                string(FE.FE_source(1)), cfg.step06k_fe_prior_summary_csv);
        end
    else
        % 632 Hz Workbench mode-shape prior from the blade model folder:
        %   f1 = 632.388 Hz; gauge x=1.4 mm, y=7.05 mm.
        % Values are displacement-per-strain gains in um/microstrain. The
        % middle value is the measured gauge-position value near Z=22.18 mm
        % (equivalent to Z=2.82 mm by symmetry); the lower/upper bounds
        % cover the inspected width-position uncertainty.
        cfg.step06k_fe_K_um_per_microstrain_band = [0.88608573343401, ...
            0.930000652601896, 0.963390968577232];
    end
end
if ~isfield(cfg, 'step06k_fe_source') || isempty(cfg.step06k_fe_source)
    cfg.step06k_fe_source = ...
        "632.388 Hz Workbench FE mode-shape prior, gauge x=1.4 mm y=7.05 mm, width-position sensitivity Z=0-25 mm; central value from measured Z=2.82/22.18 mm";
end
if ~isfield(cfg, 'step06k_grid_K_count') || isempty(cfg.step06k_grid_K_count)
    cfg.step06k_grid_K_count = 31;
end
if ~isfield(cfg, 'step06k_grid_phase_count') || isempty(cfg.step06k_grid_phase_count)
    cfg.step06k_grid_phase_count = 181;
end
if ~isfield(cfg, 'step06k_btt_correction_use_core_only') || ...
        isempty(cfg.step06k_btt_correction_use_core_only)
    cfg.step06k_btt_correction_use_core_only = true;
end
if ~isfield(cfg, 'step06k_step05_csv') || isempty(cfg.step06k_step05_csv)
    cfg.step06k_step05_csv = fullfile(cfg.step05_output_dir, ...
        cfg.dynamic_cases{1}, sprintf( ...
        'Trend_Step05_SingleSyncDirectTemplate_B%d_S%s_20241106.csv', ...
        cfg.step05_target_blades(1), sprintf('%d', cfg.step05_analysis_sensors)));
end
end


function files = resolve_step06k_files_20241106_local(cfg, caseName)
files = struct();
files.step06a_mat = fullfile(cfg.step06a_output_dir, caseName, ...
    'Step06A_StrainRPM_ResonanceEvidence_20241106.mat');
files.step03_bundle_mat = fullfile(cfg.step03_output_dir, caseName, ...
    'BTT_Observation_Bundle_20241106.mat');
files.step05_csv = cfg.step06k_step05_csv;
end


function assert_step06k_files_20241106_local(files)
mustExist = {'step06a_mat', 'step03_bundle_mat'};
for i = 1:numel(mustExist)
    f = files.(mustExist{i});
    if ~isfile(f)
        error('Missing input %s:\n  %s', mustExist{i}, f);
    end
end
end


function Region = select_reference_region_20241106_local(RegionTable, cfg)
if isempty(RegionTable)
    error('Step06A RegionTable is empty.');
end
R = RegionTable;
if isfield(cfg, 'step06k_reference_region_id') && ...
        ~isempty(cfg.step06k_reference_region_id)
    m = R.region_id == cfg.step06k_reference_region_id;
else
    m = true(height(R), 1);
end
R = R(m, :);
if isempty(R)
    error('No Step06A region matches the requested reference region.');
end
[~, ix] = max(R.peak_amp_microstrain);
Region = R(ix, :);
end


function Basis = build_region_strain_basis_20241106_local(Evidence, Region, oprTimes, cfg)
t0 = Region.time_start_s(1);
t1 = Region.time_end_s(1);
eo = Region.dominant_order(1);
t = Evidence.StrainTimeBtt(:);
y = Evidence.StrainValue(:);
wide = t >= t0 - cfg.step06k_demod_margin_s & ...
    t <= t1 + cfg.step06k_demod_margin_s & isfinite(t) & isfinite(y);
tw = t(wide);
yw = y(wide);
coreWide = tw >= t0 & tw <= t1;
yw = yw - median(yw(coreWide), 'omitnan');
theta = map_time_to_rotor_phase_20241106_local(oprTimes, tw, ...
    cfg.opr_events_per_revolution);
psi = eo .* theta;
dt = median(diff(tw), 'omitnan');
fs = 1 / dt;
z = yw .* exp(-1i * psi);
q = fft_lowpass_complex_20241106_local(z, fs, cfg.step06k_envelope_lowpass_hz);
h = 2 .* q .* exp(1i * psi);
core = tw >= t0 & tw <= t1 & isfinite(real(h)) & isfinite(imag(h));
Tref = table();
Tref.time_s = tw(core);
Tref.region_id = repmat(Region.region_id(1), nnz(core), 1);
Tref.dominant_order = repmat(eo, nnz(core), 1);
Tref.strain_basis_real_microstrain = real(h(core));
Tref.strain_basis_imag_microstrain = imag(h(core));
Tref.strain_basis_abs_microstrain = abs(h(core));
Tref.strain_raw_microstrain = yw(core);
Tref.theta_rot_rad = theta(core);
Tref.psi_rad = psi(core);
Basis = struct();
Basis.Tref = sortrows(Tref, 'time_s');
Basis.region = Region;
Basis.fs_hz = fs;
end


function FitPoints = build_btt_fit_points_20241106_local(Obs, Basis, oprTimes, cfg)
R = Basis.region;
Tref = Basis.Tref;
m = Obs.is_valid & Obs.blade_id == cfg.step06k_target_blade & ...
    ismember(Obs.sensor_id, cfg.step06k_sensor_ids) & ...
    Obs.arrival_time_s >= R.time_start_s(1) & ...
    Obs.arrival_time_s <= R.time_end_s(1) & ...
    isfinite(Obs.arrival_time_s) & isfinite(Obs.displacement_mm);
T = Obs(m, {'arrival_time_s', 'displacement_mm', 'sensor_id', ...
    'blade_id', 'rpm'});
T.region_id = repmat(R.region_id(1), height(T), 1);
T.dominant_order = repmat(R.dominant_order(1), height(T), 1);
T.theta_rot_rad = map_time_to_rotor_phase_20241106_local(oprTimes, ...
    T.arrival_time_s, cfg.opr_events_per_revolution);
T.strain_basis_real_microstrain = interp1(Tref.time_s, ...
    Tref.strain_basis_real_microstrain, T.arrival_time_s, 'linear', NaN);
T.strain_basis_imag_microstrain = interp1(Tref.time_s, ...
    Tref.strain_basis_imag_microstrain, T.arrival_time_s, 'linear', NaN);
T.strain_basis_abs_microstrain = hypot(T.strain_basis_real_microstrain, ...
    T.strain_basis_imag_microstrain);
ok = isfinite(T.strain_basis_real_microstrain) & ...
    isfinite(T.strain_basis_imag_microstrain);
FitPoints = sortrows(T(ok, :), {'arrival_time_s', 'sensor_id'});
if height(FitPoints) < 20
    error('Too few BTT correction points in the selected resonance region.');
end
end


function F = build_fe_spectrum_reference_20241106_local(Tref, cfg)
coreMask = select_core_mask_20241106_local(Tref.time_s, cfg);
hAbs = Tref.strain_basis_abs_microstrain;
strainAmp = median(hAbs(coreMask), 'omitnan');
if ~isfinite(strainAmp) || strainAmp <= 0
    strainAmp = sqrt(2) * std(Tref.strain_raw_microstrain(coreMask), ...
        'omitnan');
end
K = sort(cfg.step06k_fe_K_um_per_microstrain_band(:));
if numel(K) == 1
    K = [K; K; K];
elseif numel(K) == 2
    K = [K(1); mean(K); K(2)];
else
    K = [K(1); median(K, 'omitnan'); K(end)];
end
A = (K / 1000) .* strainAmp;
F = table();
F.reference_name = "FE_strain_spectrum_only_core_region";
F.region_id = mode(Tref.region_id);
F.dominant_order = mode(Tref.dominant_order);
F.time_start_s = min(Tref.time_s, [], 'omitnan');
F.time_end_s = max(Tref.time_s, [], 'omitnan');
F.core_time_start_s = min(Tref.time_s(coreMask), [], 'omitnan');
F.core_time_end_s = max(Tref.time_s(coreMask), [], 'omitnan');
F.trim_fraction = cfg.step06k_core_trim_fraction;
F.strain_spectrum_amp_microstrain = strainAmp;
F.FE_K_min_um_per_microstrain = K(1);
F.FE_K_median_um_per_microstrain = K(2);
F.FE_K_max_um_per_microstrain = K(3);
F.A_tip_min_mm = A(1);
F.A_tip_median_mm = A(2);
F.A_tip_max_mm = A(3);
F.core_sample_count = nnz(coreMask);
F.total_sample_count = height(Tref);
F.FE_source = string(cfg.step06k_fe_source);
end


function Tcorr = select_btt_correction_points_20241106_local(T, F, cfg)
Tcorr = T;
if cfg.step06k_btt_correction_use_core_only
    m = T.arrival_time_s >= F.core_time_start_s(1) & ...
        T.arrival_time_s <= F.core_time_end_s(1);
    if nnz(m) >= max(30, ceil(0.20 * height(T)))
        Tcorr = T(m, :);
    end
end
Tcorr.is_used_for_primary_BTT_correction = true(height(Tcorr), 1);
end


function [Fit, Grid] = build_fe_btt_reference_20241106_local(T, Tref, F, cfg)
hBtt = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
hRef = Tref.strain_basis_real_microstrain + 1i*Tref.strain_basis_imag_microstrain;
y = T.displacement_mm(:);
kGrid = linspace(F.FE_K_min_um_per_microstrain(1) / 1000, ...
    F.FE_K_max_um_per_microstrain(1) / 1000, cfg.step06k_grid_K_count);
phaseGrid = linspace(-pi, pi, cfg.step06k_grid_phase_count);
rows = repmat(make_grid_row_20241106_local(), numel(kGrid)*numel(phaseGrid), 1);
rowIdx = 0;
for iK = 1:numel(kGrid)
    for iP = 1:numel(phaseGrid)
        rowIdx = rowIdx + 1;
        C = kGrid(iK) * exp(1i*phaseGrid(iP));
        u = real(C .* hBtt);
        B = mean(y - u, 'omitnan');
        yFit = B + u;
        resid = y - yFit;
        uRef = real(C .* hRef);
        rows(rowIdx).K_mm_per_microstrain = kGrid(iK);
        rows(rowIdx).K_um_per_microstrain = 1000*kGrid(iK);
        rows(rowIdx).phase_lag_rad = wrap_to_pi_20241106_local(phaseGrid(iP));
        rows(rowIdx).common_offset_mm = B;
        rows(rowIdx).rmse_mm = sqrt(mean(resid.^2, 'omitnan'));
        rows(rowIdx).median_abs_residual_mm = median(abs(resid), 'omitnan');
        rows(rowIdx).sensor_layering_score_mm = ...
            sensor_layering_score_20241106_local(T, resid);
        rows(rowIdx).objective = rows(rowIdx).rmse_mm;
        rows(rowIdx).target_A_peak_equiv_mm = sqrt(2) * std(uRef, 'omitnan');
        rows(rowIdx).target_A_half_range_mm = 0.5 * ...
            (max(uRef, [], 'omitnan') - min(uRef, [], 'omitnan'));
    end
end
Grid = sortrows(struct2table(rows), 'objective');
best = Grid(1, :);
Cbest = best.K_mm_per_microstrain(1) * exp(1i*best.phase_lag_rad(1));
uBtt = real(Cbest .* hBtt);
Bbest = mean(y - uBtt, 'omitnan');
uRef = real(Cbest .* hRef);
Fit = struct();
Fit.model_name = "primary_FE_strain_BTT_residual_corrected";
Fit.C_complex_mm_per_microstrain = Cbest;
Fit.K_mm_per_microstrain = abs(Cbest);
Fit.K_um_per_microstrain = 1000 * abs(Cbest);
Fit.phase_lag_rad = angle(Cbest);
Fit.phase_lag_deg = angle(Cbest) * 180/pi;
Fit.common_offset_mm = Bbest;
Fit.u_ref_mm = uRef;
Fit.u_at_points_mm = uBtt;
Fit.total_fit_at_points_mm = Bbest + uBtt;
Fit.dynamic_BTT_at_points_mm = y - Bbest;
Fit.residual_at_points_mm = y - Fit.total_fit_at_points_mm;
Fit.rmse_mm = sqrt(mean(Fit.residual_at_points_mm.^2, 'omitnan'));
Fit.sensor_layering_score_mm = sensor_layering_score_20241106_local( ...
    T, Fit.residual_at_points_mm);
Fit.point_count = height(T);
Fit.sensor_ids = unique(T.sensor_id(:)).';
Fit.A_peak_equiv_mm = sqrt(2) * std(uRef, 'omitnan');
Fit.A_half_range_mm = 0.5 * (max(uRef, [], 'omitnan') - ...
    min(uRef, [], 'omitnan'));
tolK = 0.5 * (kGrid(end) - kGrid(1)) / max(numel(kGrid) - 1, 1);
if abs(Fit.K_mm_per_microstrain - kGrid(1)) <= tolK
    Fit.constraint_status = "at_FE_lower_bound";
elseif abs(Fit.K_mm_per_microstrain - kGrid(end)) <= tolK
    Fit.constraint_status = "at_FE_upper_bound";
else
    Fit.constraint_status = "inside_FE_band";
end
Fit.time_s = Tref.time_s;
end


function row = make_grid_row_20241106_local()
row = struct('K_mm_per_microstrain', NaN, ...
    'K_um_per_microstrain', NaN, ...
    'phase_lag_rad', NaN, ...
    'common_offset_mm', NaN, ...
    'rmse_mm', NaN, ...
    'median_abs_residual_mm', NaN, ...
    'sensor_layering_score_mm', NaN, ...
    'objective', NaN, ...
    'target_A_peak_equiv_mm', NaN, ...
    'target_A_half_range_mm', NaN);
end


function TS = build_time_series_20241106_local(Tref, Fit, F)
TS = Tref;
h = Tref.strain_basis_real_microstrain + 1i*Tref.strain_basis_imag_microstrain;
Cfe = (F.FE_K_median_um_per_microstrain(1) / 1000) * ...
    exp(1i * angle(Fit.C_complex_mm_per_microstrain));
TS.u_ref_FE_strain_only_median_mm = real(Cfe .* h);
TS.u_ref_FE_BTT_corrected_mm = Fit.u_ref_mm(:);
TS.u_ref_main_mm = TS.u_ref_FE_BTT_corrected_mm;
end


function P = attach_btt_predictions_20241106_local(P, Fit)
P.u_ref_FE_BTT_corrected_at_btt_mm = Fit.u_at_points_mm(:);
P.common_offset_main_mm = repmat(Fit.common_offset_mm, height(P), 1);
P.dynamic_BTT_common_offset_removed_mm = ...
    P.displacement_mm - Fit.common_offset_mm;
P.total_fit_main_mm = Fit.total_fit_at_points_mm(:);
P.residual_main_mm = Fit.residual_at_points_mm(:);
end


function S = build_sensor_residuals_20241106_local(P)
sids = unique(P.sensor_id(:)).';
rows = repmat(struct('sensor_id', NaN, 'point_count', NaN, ...
    'mean_residual_mm', NaN, 'median_residual_mm', NaN, ...
    'std_residual_mm', NaN, 'mean_abs_residual_mm', NaN), numel(sids), 1);
for i = 1:numel(sids)
    sid = sids(i);
    m = P.sensor_id == sid;
    r = P.residual_main_mm(m);
    rows(i).sensor_id = sid;
    rows(i).point_count = nnz(m);
    rows(i).mean_residual_mm = mean(r, 'omitnan');
    rows(i).median_residual_mm = median(r, 'omitnan');
    rows(i).std_residual_mm = std(r, 'omitnan');
    rows(i).mean_abs_residual_mm = mean(abs(r), 'omitnan');
end
S = struct2table(rows);
end


function T = load_step05_comparison_20241106_local(file)
if ~isfile(file)
    T = table();
    T.note = "Step05 CSV not found";
    return;
end
Raw = readtable(file, 'VariableNamingRule', 'preserve');
ok = true(height(Raw), 1);
if ismember('status', Raw.Properties.VariableNames)
    ok = strcmp(string(Raw.status), "ok");
end
if ismember('A_id', Raw.Properties.VariableNames)
    A = Raw.A_id(ok);
else
    A = NaN;
end
T = table();
T.source_file = string(file);
T.window_count = nnz(ok);
T.A_id_median_mm = median(A, 'omitnan');
T.A_id_max_mm = max(A, [], 'omitnan');
T.A_id_at_max_mm = T.A_id_max_mm;
if ismember('window_center_time_s', Raw.Properties.VariableNames)
    tt = Raw.window_center_time_s(ok);
    [~, ix] = max(A);
    if ~isempty(ix) && isfinite(ix)
        T.time_at_A_max_s = tt(ix);
    else
        T.time_at_A_max_s = NaN;
    end
else
    T.time_at_A_max_s = NaN;
end
end


function Summary = build_summary_20241106_local(caseName, Region, F, Fit, P, ...
    SensorResiduals, Step05Comparison, cfg)
Summary = table();
Summary.case_name = string(caseName);
Summary.reference_definition = ...
    "Branch1_FE_strain_core_range_plus_Branch2_FE_strain_BTT_residual_correction";
Summary.reference_status = ...
    "632Hz_FE_prior_with_width_position_uncertainty_plus_BTT_residual_diagnostic";
Summary.target_blade = cfg.step06k_target_blade;
Summary.sensor_ids = string(mat2str(cfg.step06k_sensor_ids));
Summary.region_id = Region.region_id(1);
Summary.dominant_order = Region.dominant_order(1);
Summary.peak_freq_hz = Region.peak_freq_hz(1);
Summary.time_start_s = Region.time_start_s(1);
Summary.time_end_s = Region.time_end_s(1);
Summary.core_time_start_s = F.core_time_start_s(1);
Summary.core_time_end_s = F.core_time_end_s(1);
Summary.strain_amp_microstrain = F.strain_spectrum_amp_microstrain(1);
Summary.FE_K_min_um_per_microstrain = F.FE_K_min_um_per_microstrain(1);
Summary.FE_K_median_um_per_microstrain = F.FE_K_median_um_per_microstrain(1);
Summary.FE_K_max_um_per_microstrain = F.FE_K_max_um_per_microstrain(1);
Summary.FE_strain_A_min_mm = F.A_tip_min_mm(1);
Summary.FE_strain_A_median_mm = F.A_tip_median_mm(1);
Summary.FE_strain_A_max_mm = F.A_tip_max_mm(1);
Summary.FE_BTT_corrected_K_um_per_microstrain = Fit.K_um_per_microstrain;
Summary.FE_BTT_corrected_core_A_mm = ...
    Fit.K_mm_per_microstrain * F.strain_spectrum_amp_microstrain(1);
Summary.FE_BTT_corrected_waveform_A_peak_equiv_mm = Fit.A_peak_equiv_mm;
Summary.FE_BTT_corrected_A_mm = Summary.FE_BTT_corrected_core_A_mm;
Summary.FE_BTT_corrected_phase_lag_deg = Fit.phase_lag_deg;
Summary.FE_BTT_constraint_status = string(Fit.constraint_status);
Summary.FE_BTT_common_offset_mm = Fit.common_offset_mm;
Summary.BTT_correction_point_count = height(P);
Summary.BTT_correction_rmse_mm = Fit.rmse_mm;
Summary.BTT_sensor_layering_score_mm = Fit.sensor_layering_score_mm;
Summary.FE_BTT_corrected_within_FE_strain_range = ...
    Summary.FE_BTT_corrected_core_A_mm >= F.A_tip_min_mm(1) & ...
    Summary.FE_BTT_corrected_core_A_mm <= F.A_tip_max_mm(1);
Summary.max_abs_sensor_median_residual_mm = ...
    max(abs(SensorResiduals.median_residual_mm), [], 'omitnan');
if ~isempty(Step05Comparison) && ismember('A_id_max_mm', ...
        Step05Comparison.Properties.VariableNames)
    Summary.step05_A_id_median_mm = Step05Comparison.A_id_median_mm(1);
    Summary.step05_A_id_max_mm = Step05Comparison.A_id_max_mm(1);
    Summary.step05_time_at_A_max_s = Step05Comparison.time_at_A_max_s(1);
else
    Summary.step05_A_id_median_mm = NaN;
    Summary.step05_A_id_max_mm = NaN;
    Summary.step05_time_at_A_max_s = NaN;
end
Summary.FE_source = string(cfg.step06k_fe_source);
end


function figFile = plot_reference_20241106_local(TS, P, F, SensorResiduals, ...
    Step05Comparison, Summary, figDir, cfg)
fig = figure('Visible', cfg.step06k_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 28, 18], ...
    'Name', 'Step06K 20241106 FE-core strain/BTT reference');
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
red = [0.635 0.078 0.184];
gray = [0.500 0.500 0.500];

targetRange = [Summary.time_start_s(1), Summary.time_end_s(1)];
coreRange = [Summary.core_time_start_s(1), Summary.core_time_end_s(1)];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
plot_downsampled_20241106_local(ax1, TS.time_s, ...
    TS.strain_basis_abs_microstrain, blue, 'strain envelope', 1.0);
xline(ax1, coreRange(1), ':', 'Color', red, 'LineWidth', 0.9);
xline(ax1, coreRange(2), ':', 'Color', red, 'LineWidth', 0.9);
xlim(ax1, targetRange);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Strain amp. (microstrain)');
title(ax1, '1 Trimmed strain core');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
bar(ax2, categorical(["FE+strain", "FE+strain+BTT"]), ...
    [Summary.FE_strain_A_median_mm, Summary.FE_BTT_corrected_A_mm], ...
    0.55, 'FaceColor', [0.70 0.70 0.70], ...
    'HandleVisibility', 'off');
errorbar(ax2, categorical("FE+strain"), Summary.FE_strain_A_median_mm, ...
    Summary.FE_strain_A_median_mm - Summary.FE_strain_A_min_mm, ...
    Summary.FE_strain_A_max_mm - Summary.FE_strain_A_median_mm, ...
    'k.', 'LineWidth', 1.0, 'CapSize', 9, ...
    'DisplayName', 'FE K band');
if ~isempty(Step05Comparison) && ismember('A_id_max_mm', ...
        Step05Comparison.Properties.VariableNames)
    yline(ax2, Step05Comparison.A_id_max_mm(1), '--', ...
        'Color', orange, 'DisplayName', 'Step05 max A');
end
ylabel(ax2, 'A_{peak,eq} (mm)');
title(ax2, '2 Amplitude branches');
legend(ax2, 'Location', 'best');

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
plot_downsampled_20241106_local(ax3, TS.time_s, ...
    TS.u_ref_FE_strain_only_median_mm, red, 'FE+strain median', 0.85);
plot_downsampled_20241106_local(ax3, TS.time_s, ...
    TS.u_ref_FE_BTT_corrected_mm, blue, 'FE+strain+BTT', 1.05);
yline(ax3, Summary.FE_strain_A_min_mm, ':', 'Color', red, ...
    'HandleVisibility', 'off');
yline(ax3, -Summary.FE_strain_A_min_mm, ':', 'Color', red, ...
    'HandleVisibility', 'off');
yline(ax3, Summary.FE_strain_A_max_mm, ':', 'Color', red, ...
    'HandleVisibility', 'off');
yline(ax3, -Summary.FE_strain_A_max_mm, ':', 'Color', red, ...
    'HandleVisibility', 'off');
xlim(ax3, targetRange);
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'u_{tip} (mm)');
title(ax3, '3 Complete time-domain reference');
legend(ax3, 'Location', 'best');

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
sids = unique(P.sensor_id(:)).';
colors = lines(max(numel(sids), 3));
for i = 1:numel(sids)
    m = P.sensor_id == sids(i);
    scatter(ax4, P.arrival_time_s(m), ...
        P.dynamic_BTT_common_offset_removed_mm(m), 10, colors(i, :), ...
        'filled', 'MarkerFaceAlpha', 0.45, ...
        'DisplayName', sprintf('S%d BTT', sids(i)));
end
plot_downsampled_20241106_local(ax4, TS.time_s, ...
    TS.u_ref_FE_BTT_corrected_mm, blue, 'reference', 1.1);
xlim(ax4, targetRange);
xlabel(ax4, 'Time (s)');
ylabel(ax4, 'BTT-B / u_{ref} (mm)');
title(ax4, '4 BTT residual correction anchors');
legend(ax4, 'Location', 'best');

ax5 = nexttile; hold(ax5, 'on'); box(ax5, 'on');
mp = isfinite(P.dynamic_BTT_common_offset_removed_mm) & ...
    isfinite(P.u_ref_FE_BTT_corrected_at_btt_mm);
scatter(ax5, P.u_ref_FE_BTT_corrected_at_btt_mm(mp), ...
    P.dynamic_BTT_common_offset_removed_mm(mp), 12, blue, ...
    'filled', 'MarkerFaceAlpha', 0.40);
if any(mp)
    xy = [P.u_ref_FE_BTT_corrected_at_btt_mm(mp); ...
        P.dynamic_BTT_common_offset_removed_mm(mp)];
    lim = [min(xy, [], 'omitnan'), max(xy, [], 'omitnan')];
    pad = 0.05 * max(diff(lim), eps);
    lim = lim + [-pad pad];
    plot(ax5, lim, lim, 'k--', 'LineWidth', 0.8);
    xlim(ax5, lim);
    ylim(ax5, lim);
    axis(ax5, 'square');
end
xlabel(ax5, 'u_{ref}(t_i) (mm)');
ylabel(ax5, 'BTT-B (mm)');
title(ax5, '5 Point consistency');

ax6 = nexttile; hold(ax6, 'on'); box(ax6, 'on');
bar(ax6, SensorResiduals.sensor_id, ...
    [SensorResiduals.median_residual_mm, SensorResiduals.mean_residual_mm], ...
    0.72);
yline(ax6, 0, 'k:', 'LineWidth', 0.8);
xlabel(ax6, 'Sensor id');
ylabel(ax6, 'Residual (mm)');
title(ax6, '6 Sensor residual diagnostic');
legend(ax6, {'median', 'mean'}, 'Location', 'best');

sg = sgtitle(fig, sprintf(['20241106 Step06K: EO%d, %.2f-%.2f s, ', ...
    'core %.2f-%.2f s, A=%.3f mm'], Summary.dominant_order(1), ...
    targetRange(1), targetRange(2), coreRange(1), coreRange(2), ...
    Summary.FE_BTT_corrected_A_mm(1)));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.5);
style_figure_20241106_local(fig);
figFile = fullfile(figDir, ...
    'Step06K_FECoreStrainBTTReference_20241106.png');
save_figure_20241106_local(fig, figFile);
end


function coreMask = select_core_mask_20241106_local(t, cfg)
t = t(:);
valid = isfinite(t);
t0 = min(t(valid));
t1 = max(t(valid));
dur = t1 - t0;
trim = cfg.step06k_core_trim_fraction * dur;
trim = max(trim, cfg.step06k_core_trim_min_s);
trim = min(trim, cfg.step06k_core_trim_max_s);
if 2*trim >= 0.8*dur
    trim = 0.1*dur;
end
coreMask = valid & t >= t0 + trim & t <= t1 - trim;
if nnz(coreMask) < 3
    coreMask = valid;
end
end


function thetaRot = map_time_to_rotor_phase_20241106_local( ...
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


function y = fft_lowpass_complex_20241106_local(x, fs, cutoffHz)
x = x(:);
valid = isfinite(real(x)) & isfinite(imag(x));
if nnz(valid) < 8
    y = nan(size(x));
    return;
end
fill = x;
idx = (1:numel(x)).';
fill(~valid) = interp1(idx(valid), x(valid), idx(~valid), ...
    'linear', 'extrap');
n = numel(fill);
X = fft(fill);
f = (0:n-1).' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
mask = abs(f) <= cutoffHz;
y = ifft(X .* mask);
y(~valid) = NaN;
end


function score = sensor_layering_score_20241106_local(T, resid)
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


function p = wrap_to_pi_20241106_local(p)
p = mod(p + pi, 2*pi) - pi;
end


function plot_downsampled_20241106_local(ax, x, y, color, name, lw)
if isempty(x) || isempty(y)
    return;
end
step = max(1, ceil(numel(x) / 30000));
plot(ax, x(1:step:end), y(1:step:end), '-', ...
    'Color', color, 'LineWidth', lw, 'DisplayName', name);
end


function style_figure_20241106_local(fig)
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


function save_figure_20241106_local(fig, pngFile)
[folder, base, ~] = fileparts(pngFile);
ensure_dir_20241106_local(folder);
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


function ensure_dir_20241106_local(d)
if ~isfolder(d)
    mkdir(d);
end
end
