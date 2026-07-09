%% Step06I_Build_EODependentJointReference_20250527.m
% Build an EO-dependent independent strain-BTT reference.
%
% Earlier diagnostic scripts used one complex gain C for all same-frequency
% resonance regions.  That is too strong when different engine orders may
% excite different effective mode/trajectory components.  This standalone
% step reads Step03 sparse BTT, Step06A strain evidence, and Step02 OPR
% phase directly, then fits one complex strain-to-tip coefficient per EO:
%
%   d_i = B_r + delta_s + Re{ C_EO(i) h_r(t_i) } + e_i
%
% The reference is independent of Step05.  Step05 is loaded only after the
% EO-dependent reference has been built, for comparison.

clc; clear; close all;

cfg = BTTDataConfig_20250527();
cfg = apply_step06i_defaults_local(cfg);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06I EO-dependent joint reference: %s ===\n', caseName);

    files = resolve_input_files_local(cfg, caseName);
    assert_input_files_local(files);

    EvidenceS = load(files.step06a_evidence_mat, 'Evidence');
    Evidence = EvidenceS.Evidence;
    RegionTableAll = get_evidence_table_local(Evidence, {'regionTable', 'RegionTable'});
    Obs = load_step03_observation_table_local(files.step03_bundle_mat);
    oprTimes = load_opr_times_local(files.jiluopr_mat);
    [tStrain, yStrain, fsStrain] = get_strain_arrays_local(Evidence);

    RegionTable = select_reference_regions_local(RegionTableAll, Obs, cfg);
    StrainBasis = build_region_strain_basis_local( ...
        RegionTable, tStrain, yStrain, fsStrain, oprTimes, cfg);
    FitPoints = build_btt_fit_points_local(Obs, StrainBasis, oprTimes, cfg);
    GlobalFit = fit_global_reference_model_local(FitPoints, cfg);
    EOFit = fit_eo_dependent_reference_local(FitPoints, cfg);
    FitPointsEO = attach_eo_fit_to_points_local(FitPoints, EOFit);
    EOGainTable = build_eo_gain_table_local(FitPointsEO, EOFit, cfg);
    ReferenceTimeSeriesBase = build_base_reference_time_series_local( ...
        StrainBasis, GlobalFit, cfg);
    ReferenceTimeSeries = build_eo_reference_time_series_local( ...
        ReferenceTimeSeriesBase, EOFit, GlobalFit, cfg);
    SensorResiduals = build_sensor_residual_table_local(FitPointsEO);
    CVTable = build_eo_cross_validation_local(FitPointsEO, cfg);
    BootstrapTable = bootstrap_eo_gain_local(FitPointsEO, cfg);

    Step05 = load_step05_for_comparison_local(files, cfg);
    ReferenceTimeSeries = attach_step05_time_series_local( ...
        ReferenceTimeSeries, Step05, oprTimes, cfg);
    WindowComparison = build_window_comparison_local( ...
        ReferenceTimeSeries, Step05, cfg);
    Summary = build_summary_table_local(caseName, GlobalFit, ...
        EOGainTable, ReferenceTimeSeries, FitPointsEO, WindowComparison, ...
        CVTable, BootstrapTable, cfg);

    outDir = fullfile(cfg.step06i_output_dir, caseName);
    figDir = fullfile(cfg.step06i_figure_dir, caseName);
    ensure_dir_local(outDir);
    ensure_dir_local(figDir);

    summaryCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_Summary_20250527.csv');
    regionsCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_SelectedRegions_20250527.csv');
    eoCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_EOGainTable_20250527.csv');
    pointsCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_FitPoints_20250527.csv');
    tsCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_TimeSeries_20250527.csv');
    sensorCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_SensorResiduals_20250527.csv');
    cvCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_CrossValidation_20250527.csv');
    bootCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_Bootstrap_20250527.csv');
    windowCsv = fullfile(outDir, ...
        'Step06I_EODependentReference_Step05WindowComparison_20250527.csv');
    matFile = fullfile(outDir, ...
        'Step06I_EODependentReference_20250527.mat');

    writetable(Summary, summaryCsv);
    writetable(RegionTable, regionsCsv);
    writetable(EOGainTable, eoCsv);
    writetable(FitPointsEO, pointsCsv);
    writetable(ReferenceTimeSeries, tsCsv);
    writetable(SensorResiduals, sensorCsv);
    writetable(CVTable, cvCsv);
    writetable(BootstrapTable, bootCsv);
    writetable(WindowComparison, windowCsv);
    save(matFile, 'Summary', 'EOGainTable', 'FitPointsEO', ...
        'ReferenceTimeSeries', 'SensorResiduals', 'CVTable', ...
        'BootstrapTable', 'WindowComparison', 'EOFit', 'Step05', ...
        'GlobalFit', 'RegionTable', 'StrainBasis', 'files', '-v7.3');

    figFile = plot_eo_reference_local(ReferenceTimeSeries, FitPointsEO, ...
        Summary, EOGainTable, RegionTable, SensorResiduals, CVTable, ...
        BootstrapTable, figDir, cfg);

    fprintf('\nStep06I EO-dependent reference summary:\n');
    disp(Summary);
    fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        summaryCsv, regionsCsv, eoCsv, pointsCsv, tsCsv, sensorCsv, ...
        cvCsv, bootCsv, windowCsv, matFile, figFile);
end


function cfg = apply_step06i_defaults_local(cfg)
if ~isfield(cfg, 'step06i_target_blade') || isempty(cfg.step06i_target_blade)
    cfg.step06i_target_blade = 1;
end
if ~isfield(cfg, 'step06i_sensor_ids') || isempty(cfg.step06i_sensor_ids)
    cfg.step06i_sensor_ids = cfg.sensor_ids(:).';
end
if ~isfield(cfg, 'step06i_sensor_group') || isempty(cfg.step06i_sensor_group)
    cfg.step06i_sensor_group = ['S', sprintf('%d', cfg.step06i_sensor_ids)];
end
if ~isfield(cfg, 'step06i_reference_sensor_id') || ...
        isempty(cfg.step06i_reference_sensor_id)
    cfg.step06i_reference_sensor_id = cfg.step06i_sensor_ids(1);
end
if ~isfield(cfg, 'step06i_target_time_range_s') || ...
        isempty(cfg.step06i_target_time_range_s)
    cfg.step06i_target_time_range_s = [1.0 5.0];
end
if ~isfield(cfg, 'step06i_zoom_time_range_s') || ...
        isempty(cfg.step06i_zoom_time_range_s)
    cfg.step06i_zoom_time_range_s = [1.50 1.56];
end
if ~isfield(cfg, 'step06i_fixed_strain_to_mm') || ...
        isempty(cfg.step06i_fixed_strain_to_mm)
    cfg.step06i_fixed_strain_to_mm = 1 / 1.266 / 1000;
end
if ~isfield(cfg, 'step06i_target_freq_hz') || isempty(cfg.step06i_target_freq_hz)
    cfg.step06i_target_freq_hz = NaN; % default: first Step06A region freq
end
if ~isfield(cfg, 'step06i_freq_tolerance_hz') || ...
        isempty(cfg.step06i_freq_tolerance_hz)
    cfg.step06i_freq_tolerance_hz = 5.0;
end
if ~isfield(cfg, 'step06i_min_region_btt_points') || ...
        isempty(cfg.step06i_min_region_btt_points)
    cfg.step06i_min_region_btt_points = 40;
end
if ~isfield(cfg, 'step06i_min_region_duration_s') || ...
        isempty(cfg.step06i_min_region_duration_s)
    cfg.step06i_min_region_duration_s = 0.20;
end
if ~isfield(cfg, 'step06i_demod_margin_s') || ...
        isempty(cfg.step06i_demod_margin_s)
    cfg.step06i_demod_margin_s = 0.80;
end
if ~isfield(cfg, 'step06i_envelope_lowpass_hz') || ...
        isempty(cfg.step06i_envelope_lowpass_hz)
    cfg.step06i_envelope_lowpass_hz = 8.0;
end
if ~isfield(cfg, 'step06i_min_eo_points') || isempty(cfg.step06i_min_eo_points)
    cfg.step06i_min_eo_points = 80;
end
if ~isfield(cfg, 'step06i_min_phase_coverage') || ...
        isempty(cfg.step06i_min_phase_coverage)
    cfg.step06i_min_phase_coverage = 0.50;
end
if ~isfield(cfg, 'step06i_huber_k') || isempty(cfg.step06i_huber_k)
    cfg.step06i_huber_k = 1.345;
end
if ~isfield(cfg, 'step06i_irls_iterations') || ...
        isempty(cfg.step06i_irls_iterations)
    cfg.step06i_irls_iterations = 12;
end
if ~isfield(cfg, 'step06i_bootstrap_count') || ...
        isempty(cfg.step06i_bootstrap_count)
    cfg.step06i_bootstrap_count = 120;
end
if ~isfield(cfg, 'step06i_random_seed') || isempty(cfg.step06i_random_seed)
    cfg.step06i_random_seed = 20250527;
end
if ~isfield(cfg, 'step06i_output_dir') || isempty(cfg.step06i_output_dir)
    cfg.step06i_output_dir = fullfile(cfg.output_root, ...
        'step06i_eo_dependent_joint_reference');
end
if ~isfield(cfg, 'step06i_figure_dir') || isempty(cfg.step06i_figure_dir)
    cfg.step06i_figure_dir = fullfile(cfg.figure_root, ...
        'step06i_eo_dependent_joint_reference');
end
if ~isfield(cfg, 'step06i_fig_visible') || isempty(cfg.step06i_fig_visible)
    if usejava('desktop')
        cfg.step06i_fig_visible = 'on';
    else
        cfg.step06i_fig_visible = 'off';
    end
end
end


function files = resolve_input_files_local(cfg, caseName)
step05Root = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template', caseName);

files = struct();
files.step03_bundle_mat = fullfile(cfg.step03_output_dir, caseName, ...
    'BTT_Observation_Bundle_20250527.mat');
files.step06a_evidence_mat = fullfile(cfg.output_root, ...
    'step06a_strain_rpm_resonance_evidence', caseName, ...
    'Step06A_StrainRPM_ResonanceEvidence_20250527.mat');
files.step05_trend_csv = fullfile(step05Root, sprintf( ...
    'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
    cfg.step06i_target_blade, cfg.step06i_sensor_group, cfg.dataset));
files.jiluopr_mat = fullfile(cfg.step02_output_dir, caseName, 'jiluOPR.mat');
end


function assert_input_files_local(files)
if ~isfile(files.step03_bundle_mat)
    error('Missing Step03 observation bundle:\n  %s', files.step03_bundle_mat);
end
if ~isfile(files.step06a_evidence_mat)
    error('Missing Step06A resonance evidence MAT:\n  %s', files.step06a_evidence_mat);
end
if ~isfile(files.jiluopr_mat)
    error('Missing Step02 OPR MAT:\n  %s', files.jiluopr_mat);
end
if ~isfile(files.step05_trend_csv)
    warning('Missing Step05 trend CSV; Step06I will still build reference:\n  %s', ...
        files.step05_trend_csv);
end
end


function T = get_evidence_table_local(Evidence, names)
T = table();
for i = 1:numel(names)
    if isfield(Evidence, names{i}) && istable(Evidence.(names{i}))
        T = Evidence.(names{i});
        return;
    end
end
error('Step06A Evidence does not contain a region table.');
end


function Obs = load_step03_observation_table_local(file)
S = load(file, 'bundle');
if ~isfield(S, 'bundle') || ~isfield(S.bundle, 'Observation_Table')
    error('Step03 bundle has no Observation_Table:\n  %s', file);
end
Obs = S.bundle.Observation_Table;
end


function [t, y, fs] = get_strain_arrays_local(Evidence)
if isfield(Evidence, 'StrainRaw') && isstruct(Evidence.StrainRaw)
    t = Evidence.StrainRaw.timeBtt(:);
    if isfield(Evidence.StrainRaw, 'detrendedValue')
        y = Evidence.StrainRaw.detrendedValue(:);
    else
        y = Evidence.StrainRaw.rawValue(:);
    end
    if isfield(Evidence.StrainRaw, 'sampleRateHz')
        fs = Evidence.StrainRaw.sampleRateHz;
    else
        fs = NaN;
    end
else
    t = Evidence.StrainTimeBtt(:);
    y = Evidence.StrainValue(:);
    fs = NaN;
end
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
[t, order] = sort(t);
y = y(order);
if ~isfinite(fs) || fs <= 0
    fs = 1 / median(diff(t), 'omitnan');
end
end


function RegionTable = select_reference_regions_local(R, Obs, cfg)
if height(R) == 0
    error('Step06A region table is empty.');
end

freq = get_table_numeric_local(R, {'dominantFreqHz', 'peak_freq_hz', ...
    'regionPeakFreqHz'});
if ~isfinite(cfg.step06i_target_freq_hz)
    targetFreq = freq(1);
else
    targetFreq = cfg.step06i_target_freq_hz;
end
duration = get_table_numeric_local(R, {'regionDuration', 'duration_s'});
startT = get_table_numeric_local(R, {'bttStartSec', 'time_start_s', 'regionStart'});
endT = get_table_numeric_local(R, {'bttEndSec', 'time_end_s', 'regionEnd'});
eo = get_table_numeric_local(R, {'dominantOrder', 'dominant_order', ...
    'regionPeakOrder'});
amp = get_table_numeric_local(R, {'peakFftAmpMicrostrain', ...
    'peak_amp_microstrain'});
rid = get_table_numeric_local(R, {'regionId', 'region_id'});

bttCount = zeros(height(R), 1);
for i = 1:height(R)
    m = Obs.is_valid & Obs.blade_id == cfg.step06i_target_blade & ...
        ismember(Obs.sensor_id, cfg.step06i_sensor_ids) & ...
        Obs.arrival_time_s >= startT(i) & Obs.arrival_time_s <= endT(i);
    bttCount(i) = nnz(m);
end

keep = abs(freq - targetFreq) <= cfg.step06i_freq_tolerance_hz & ...
    duration >= cfg.step06i_min_region_duration_s & ...
    bttCount >= cfg.step06i_min_region_btt_points & ...
    isfinite(eo) & eo > 0;

RegionTable = table();
RegionTable.region_id = rid(keep);
RegionTable.time_start_s = startT(keep);
RegionTable.time_end_s = endT(keep);
RegionTable.duration_s = duration(keep);
RegionTable.dominant_order = eo(keep);
RegionTable.peak_freq_hz = freq(keep);
RegionTable.peak_amp_microstrain = amp(keep);
RegionTable.btt_point_count = bttCount(keep);
RegionTable.target_freq_hz = repmat(targetFreq, nnz(keep), 1);
RegionTable.freq_delta_hz = RegionTable.peak_freq_hz - targetFreq;

if isempty(RegionTable)
    error('No Step06I reference regions selected.');
end
RegionTable = sortrows(RegionTable, 'time_start_s');
end


function x = get_table_numeric_local(T, names)
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


function Basis = build_region_strain_basis_local( ...
    R, tStrain, yStrain, fs, oprTimes, cfg)
rows = struct([]);
for i = 1:height(R)
    t0 = R.time_start_s(i);
    t1 = R.time_end_s(i);
    eo = R.dominant_order(i);
    wide = tStrain >= t0 - cfg.step06i_demod_margin_s & ...
        tStrain <= t1 + cfg.step06i_demod_margin_s;
    tw = tStrain(wide);
    yw = yStrain(wide);
    core = tw >= t0 & tw <= t1;
    yw = yw - median(yw(core), 'omitnan');
    theta = map_time_to_rotor_phase_local(oprTimes, tw, ...
        cfg.opr_events_per_revolution);
    psi = eo .* theta;
    z = yw .* exp(-1i * psi);
    q = fft_lowpass_complex_local(z, fs, cfg.step06i_envelope_lowpass_hz);
    h = 2 .* q .* exp(1i * psi);

    tc = tw(core);
    hc = h(core);
    qc = q(core);
    pc = psi(core);
    finite = isfinite(tc) & isfinite(real(hc)) & isfinite(imag(hc));
    rows(i).region_id = R.region_id(i); %#ok<AGROW>
    rows(i).time_start_s = t0;
    rows(i).time_end_s = t1;
    rows(i).dominant_order = eo;
    rows(i).time_s = tc(finite);
    rows(i).h_complex_microstrain = hc(finite);
    rows(i).q_complex_microstrain = qc(finite);
    rows(i).psi_rad = pc(finite);
end
Basis = rows;
end


function y = fft_lowpass_complex_local(x, fs, cutoffHz)
x = x(:);
valid = isfinite(real(x)) & isfinite(imag(x));
if nnz(valid) < 8
    y = nan(size(x));
    return;
end
fill = x;
fill(~valid) = interp1(find(valid), x(valid), find(~valid), 'linear', 'extrap');
n = numel(fill);
X = fft(fill);
f = (0:n-1).' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
mask = abs(f) <= cutoffHz;
y = ifft(X .* mask);
y(~valid) = NaN;
end


function FitPoints = build_btt_fit_points_local(Obs, Basis, oprTimes, cfg)
All = table();
for i = 1:numel(Basis)
    B = Basis(i);
    mask = Obs.is_valid & ...
        Obs.blade_id == cfg.step06i_target_blade & ...
        ismember(Obs.sensor_id, cfg.step06i_sensor_ids) & ...
        Obs.arrival_time_s >= B.time_start_s & ...
        Obs.arrival_time_s <= B.time_end_s & ...
        isfinite(Obs.arrival_time_s) & isfinite(Obs.displacement_mm);
    T = Obs(mask, {'arrival_time_s', 'displacement_mm', ...
        'sensor_id', 'blade_id', 'rpm'});
    if isempty(T)
        continue;
    end
    hr = interp1(B.time_s, real(B.h_complex_microstrain), ...
        T.arrival_time_s, 'linear', NaN);
    hi = interp1(B.time_s, imag(B.h_complex_microstrain), ...
        T.arrival_time_s, 'linear', NaN);
    T.region_id = repmat(B.region_id, height(T), 1);
    T.dominant_order = repmat(B.dominant_order, height(T), 1);
    T.theta_rot_rad = map_time_to_rotor_phase_local(oprTimes, ...
        T.arrival_time_s, cfg.opr_events_per_revolution);
    T.strain_basis_real_microstrain = hr;
    T.strain_basis_imag_microstrain = hi;
    T.strain_basis_abs_microstrain = hypot(hr, hi);
    T.strain_basis_phase_rad = atan2(hi, hr);
    T = T(isfinite(T.strain_basis_real_microstrain) & ...
        isfinite(T.strain_basis_imag_microstrain), :);
    All = [All; T]; %#ok<AGROW>
end
FitPoints = sortrows(All, 'arrival_time_s');
if height(FitPoints) < 20
    error('Too few joint strain-BTT fit points.');
end
end


function Fit = fit_global_reference_model_local(T, cfg)
[X, meta] = build_global_design_matrix_local(T, cfg, true, true);
y = T.displacement_mm(:);
[beta, weights] = robust_irls_local(X, y, cfg);
yFit = X * beta;
h = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
C = beta(1) + 1i*beta(2);

Fit = struct();
Fit.beta = beta;
Fit.weights = weights;
Fit.C_complex_mm_per_microstrain = C;
Fit.K_mm_per_microstrain = abs(C);
Fit.phase_lag_rad = angle(C);
Fit.K_ratio_to_fixed = abs(C) / cfg.step06i_fixed_strain_to_mm;
Fit.y_fit_mm = yFit;
Fit.u_global_at_points_mm = real(C .* h);
Fit.residual_mm = y - yFit;
Fit.rmse_mm = rmse_local(y, yFit);
Fit.r2 = compute_r2_local(y, yFit);
Fit.weighted_rmse_mm = sqrt(sum(weights .* Fit.residual_mm.^2) / sum(weights));
Fit.phase_coverage_fraction = phase_coverage_local(angle(C .* h)) / (2*pi);
Fit.region_ids = meta.regionIds;
Fit.sensor_ids = meta.sensorIds;
Fit.rank_X = rank(X);
Fit.cond_X = cond(X);
end


function [X, meta] = build_global_design_matrix_local(T, cfg, useRegion, useSensor)
h = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
X = [real(h), -imag(h)];
meta = struct();
meta.regionIds = unique(T.region_id(:)).';
if useRegion
    for rid = meta.regionIds
        X(:, end+1) = double(T.region_id == rid); %#ok<AGROW>
    end
end
meta.sensorIds = unique(T.sensor_id(:)).';
if useSensor
    for sid = meta.sensorIds
        if sid == cfg.step06i_reference_sensor_id
            continue;
        end
        X(:, end+1) = double(T.sensor_id == sid); %#ok<AGROW>
    end
end
end


function Tref = build_base_reference_time_series_local(Basis, GlobalFit, cfg)
Tref = table();
for i = 1:numel(Basis)
    B = Basis(i);
    h = B.h_complex_microstrain(:);
    t = B.time_s(:);
    Tr = table();
    Tr.time_s = t;
    Tr.region_id = repmat(B.region_id, numel(t), 1);
    Tr.dominant_order = repmat(B.dominant_order, numel(t), 1);
    Tr.strain_basis_real_microstrain = real(h);
    Tr.strain_basis_imag_microstrain = imag(h);
    Tr.strain_basis_abs_microstrain = abs(h);
    Tr.u_fixedK_mm = cfg.step06i_fixed_strain_to_mm .* real(h);
    Tr.u_ref_globalK_allEO_mm = real( ...
        GlobalFit.C_complex_mm_per_microstrain .* h);
    Tref = [Tref; Tr]; %#ok<AGROW>
end
Tref = sortrows(Tref, 'time_s');
end


function Fit = fit_eo_dependent_reference_local(T, cfg)
[X, meta] = build_eo_design_matrix_local(T, cfg, true, true);
y = T.displacement_mm(:);
[beta, weights] = robust_irls_local(X, y, cfg);
yFit = X * beta;
h = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
u = nan(height(T), 1);

Fit = struct();
Fit.beta = beta;
Fit.weights = weights;
Fit.eo_ids = meta.eoIds;
Fit.C_eo = nan(numel(meta.eoIds), 1);
for i = 1:numel(meta.eoIds)
    cols = meta.eoCols(i, :);
    C = beta(cols(1)) + 1i*beta(cols(2));
    Fit.C_eo(i) = C;
    m = T.dominant_order == meta.eoIds(i);
    u(m) = real(C .* h(m));
end
Fit.region_ids = meta.regionIds;
Fit.region_offsets_mm = beta(meta.regionCols);
Fit.sensor_ids = meta.sensorIds;
Fit.sensor_bias_mm = zeros(numel(meta.sensorIds), 1);
for i = 1:numel(meta.sensorIds)
    sid = meta.sensorIds(i);
    if sid == cfg.step06i_reference_sensor_id
        Fit.sensor_bias_mm(i) = 0;
    else
        col = meta.sensorCols(meta.sensorColIds == sid);
        Fit.sensor_bias_mm(i) = beta(col);
    end
end
Fit.y_fit_mm = yFit;
Fit.u_eo_at_points_mm = u;
Fit.residual_mm = y - yFit;
Fit.rmse_mm = rmse_local(y, yFit);
Fit.r2 = compute_r2_local(y, yFit);
Fit.weighted_rmse_mm = sqrt(sum(weights .* Fit.residual_mm.^2) / sum(weights));
Fit.rank_X = rank(X);
Fit.cond_X = cond(X);
end


function [X, meta] = build_eo_design_matrix_local(T, cfg, useRegion, useSensor)
h = T.strain_basis_real_microstrain + 1i*T.strain_basis_imag_microstrain;
X = [];

meta = struct();
meta.eoIds = unique(T.dominant_order(:)).';
meta.eoCols = zeros(numel(meta.eoIds), 2);
for i = 1:numel(meta.eoIds)
    eo = meta.eoIds(i);
    m = double(T.dominant_order == eo);
    X(:, end+1) = m .* real(h); %#ok<AGROW>
    X(:, end+1) = -m .* imag(h); %#ok<AGROW>
    meta.eoCols(i, :) = [size(X, 2)-1, size(X, 2)];
end

meta.regionIds = unique(T.region_id(:)).';
meta.regionCols = [];
if useRegion
    for rid = meta.regionIds
        X(:, end+1) = double(T.region_id == rid); %#ok<AGROW>
        meta.regionCols(end+1) = size(X, 2); %#ok<AGROW>
    end
end

meta.sensorIds = unique(T.sensor_id(:)).';
meta.sensorCols = [];
meta.sensorColIds = [];
if useSensor
    for sid = meta.sensorIds
        if sid == cfg.step06i_reference_sensor_id
            continue;
        end
        X(:, end+1) = double(T.sensor_id == sid); %#ok<AGROW>
        meta.sensorCols(end+1) = size(X, 2); %#ok<AGROW>
        meta.sensorColIds(end+1) = sid; %#ok<AGROW>
    end
end
end


function [beta, weights] = robust_irls_local(X, y, cfg)
valid = all(isfinite(X), 2) & isfinite(y);
Xv = X(valid, :);
yv = y(valid);
weightsV = ones(size(yv));
beta = zeros(size(X, 2), 1);
for iter = 1:cfg.step06i_irls_iterations
    sw = sqrt(weightsV);
    betaV = (Xv .* sw) \ (yv .* sw);
    r = yv - Xv * betaV;
    s = robust_scale_local(r);
    if ~isfinite(s) || s <= eps
        beta = betaV;
        break;
    end
    cutoff = cfg.step06i_huber_k * s;
    weightsV = min(1, cutoff ./ max(abs(r), eps));
    beta = betaV;
end
weights = zeros(size(y));
weights(valid) = weightsV;
end


function s = robust_scale_local(r)
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


function T = attach_eo_fit_to_points_local(T, Fit)
n = height(T);
T.u_eo_ref_mm = Fit.u_eo_at_points_mm;
T.eo_fit_total_mm = Fit.y_fit_mm;
T.eo_fit_residual_mm = Fit.residual_mm;
T.eo_fit_weight = Fit.weights;
T.eo_K_mm_per_microstrain = nan(n, 1);
T.eo_phase_lag_rad = nan(n, 1);
for i = 1:numel(Fit.eo_ids)
    m = T.dominant_order == Fit.eo_ids(i);
    T.eo_K_mm_per_microstrain(m) = abs(Fit.C_eo(i));
    T.eo_phase_lag_rad(m) = angle(Fit.C_eo(i));
end
T.region_offset_eo_mm = nan(n, 1);
for i = 1:numel(Fit.region_ids)
    m = T.region_id == Fit.region_ids(i);
    T.region_offset_eo_mm(m) = Fit.region_offsets_mm(i);
end
T.sensor_bias_eo_mm = nan(n, 1);
for i = 1:numel(Fit.sensor_ids)
    m = T.sensor_id == Fit.sensor_ids(i);
    T.sensor_bias_eo_mm(m) = Fit.sensor_bias_mm(i);
end
T.displacement_dynamic_eo_mm = T.displacement_mm - ...
    T.region_offset_eo_mm - T.sensor_bias_eo_mm;
end


function EOT = build_eo_gain_table_local(T, Fit, cfg)
rows = struct([]);
for i = 1:numel(Fit.eo_ids)
    eo = Fit.eo_ids(i);
    C = Fit.C_eo(i);
    m = T.dominant_order == eo;
    r = T.eo_fit_residual_mm(m);
    phaseCoverage = phase_coverage_local(angle( ...
        T.strain_basis_real_microstrain(m) + ...
        1i*T.strain_basis_imag_microstrain(m))) / (2*pi);
    status = "ok";
    if nnz(m) < cfg.step06i_min_eo_points
        status = "diagnostic_too_few_points";
    elseif phaseCoverage < cfg.step06i_min_phase_coverage
        status = "diagnostic_low_phase_coverage";
    end
    rows(i).dominant_order = eo; %#ok<AGROW>
    rows(i).region_ids = string(mat2str(unique(T.region_id(m)).'));
    rows(i).point_count = nnz(m);
    rows(i).phase_coverage_fraction = phaseCoverage;
    rows(i).K_mm_per_microstrain = abs(C);
    rows(i).phase_lag_rad = angle(C);
    rows(i).phase_lag_deg = angle(C) * 180 / pi;
    rows(i).K_ratio_to_fixed = abs(C) / cfg.step06i_fixed_strain_to_mm;
    rows(i).rmse_mm = sqrt(mean(r.^2, 'omitnan'));
    rows(i).mean_abs_residual_mm = mean(abs(r), 'omitnan');
    rows(i).r2_within_eo = compute_r2_local( ...
        T.displacement_mm(m), T.eo_fit_total_mm(m));
    rows(i).status = status;
end
EOT = struct2table(rows);
end


function Tref = build_eo_reference_time_series_local(Tref, Fit, GlobalFit, cfg)
h = Tref.strain_basis_real_microstrain + 1i*Tref.strain_basis_imag_microstrain;
Tref.u_ref_globalK_allEO_mm = real( ...
    GlobalFit.C_complex_mm_per_microstrain .* h);
Tref.u_ref_eo_mm = nan(height(Tref), 1);
Tref.u_fixedK_eo_phase_mm = nan(height(Tref), 1);
Tref.eo_K_mm_per_microstrain = nan(height(Tref), 1);
Tref.eo_phase_lag_rad = nan(height(Tref), 1);
for i = 1:numel(Fit.eo_ids)
    eo = Fit.eo_ids(i);
    C = Fit.C_eo(i);
    m = Tref.dominant_order == eo;
    Tref.u_ref_eo_mm(m) = real(C .* h(m));
    Cfixed = cfg.step06i_fixed_strain_to_mm * exp(1i*angle(C));
    Tref.u_fixedK_eo_phase_mm(m) = real(Cfixed .* h(m));
    Tref.eo_K_mm_per_microstrain(m) = abs(C);
    Tref.eo_phase_lag_rad(m) = angle(C);
end
Tref.global_K_mm_per_microstrain = repmat( ...
    GlobalFit.K_mm_per_microstrain, height(Tref), 1);
Tref.global_phase_lag_rad = repmat(GlobalFit.phase_lag_rad, height(Tref), 1);
end


function S = build_sensor_residual_table_local(T)
sids = unique(T.sensor_id(:)).';
rows = struct([]);
for i = 1:numel(sids)
    sid = sids(i);
    m = T.sensor_id == sid;
    r = T.eo_fit_residual_mm(m);
    rows(i).sensor_id = sid; %#ok<AGROW>
    rows(i).point_count = nnz(m);
    rows(i).sensor_bias_eo_mm = median(T.sensor_bias_eo_mm(m), 'omitnan');
    rows(i).mean_residual_mm = mean(r, 'omitnan');
    rows(i).median_residual_mm = median(r, 'omitnan');
    rows(i).std_residual_mm = std(r, 'omitnan');
    rows(i).mean_abs_residual_mm = mean(abs(r), 'omitnan');
end
S = struct2table(rows);
end


function CV = build_eo_cross_validation_local(T, cfg)
eoIds = unique(T.dominant_order(:)).';
rows = repmat(make_cv_row_local(), numel(eoIds), 1);
for i = 1:numel(eoIds)
    eo = eoIds(i);
    train = T(T.dominant_order ~= eo, :);
    test = T(T.dominant_order == eo, :);
    rows(i) = evaluate_leave_one_eo_local(train, test, eo, cfg);
end
CV = struct2table(rows);
end


function row = evaluate_leave_one_eo_local(train, test, eo, cfg)
row = make_cv_row_local();
row.held_out_eo = eo;
row.train_count = height(train);
row.test_count = height(test);
if height(train) < 50 || height(test) < cfg.step06i_min_eo_points
    row.status = "too_few_points";
    return;
end
FitTrain = fit_eo_dependent_reference_local(train, cfg);
globalC = weighted_global_c_from_eo_local(FitTrain, train);
hTest = test.strain_basis_real_microstrain + 1i*test.strain_basis_imag_microstrain;
u = real(globalC .* hTest);
[Xoff, ~] = build_validation_offset_matrix_local(test, cfg);
betaOff = Xoff \ (test.displacement_mm - u);
yPred = u + Xoff * betaOff;
r = test.displacement_mm - yPred;
row.K_train_global_mm_per_microstrain = abs(globalC);
row.phase_train_global_rad = angle(globalC);
row.rmse_test_mm = rmse_local(test.displacement_mm, yPred);
row.mean_abs_test_mm = mean(abs(r), 'omitnan');
row.r2_test = compute_r2_local(test.displacement_mm, yPred);
row.status = "ok_global_fallback";
end


function C = weighted_global_c_from_eo_local(Fit, T)
num = 0;
den = 0;
for i = 1:numel(Fit.eo_ids)
    n = nnz(T.dominant_order == Fit.eo_ids(i));
    num = num + n * Fit.C_eo(i);
    den = den + n;
end
C = num / max(den, 1);
end


function [Xoff, meta] = build_validation_offset_matrix_local(T, cfg)
meta = struct();
meta.regionIds = unique(T.region_id(:)).';
Xoff = zeros(height(T), 0);
for rid = meta.regionIds
    Xoff(:, end+1) = double(T.region_id == rid); %#ok<AGROW>
end
if isempty(meta.regionIds)
    Xoff = ones(height(T), 1);
end
sensorIds = unique(T.sensor_id(:)).';
for sid = sensorIds
    if sid == cfg.step06i_reference_sensor_id
        continue;
    end
    Xoff(:, end+1) = double(T.sensor_id == sid); %#ok<AGROW>
end
end


function row = make_cv_row_local()
row = struct('held_out_eo', NaN, 'train_count', NaN, 'test_count', NaN, ...
    'K_train_global_mm_per_microstrain', NaN, ...
    'phase_train_global_rad', NaN, 'rmse_test_mm', NaN, ...
    'mean_abs_test_mm', NaN, 'r2_test', NaN, 'status', "");
end


function Boot = bootstrap_eo_gain_local(T, cfg)
rng(cfg.step06i_random_seed);
targetEO = infer_target_eo_from_points_local(T, cfg);
nBoot = cfg.step06i_bootstrap_count;
rows = struct('bootstrap_id', num2cell((1:nBoot).'), ...
    'target_eo', targetEO, 'K_target_eo_mm_per_microstrain', NaN, ...
    'phase_target_eo_rad', NaN, 'rmse_mm', NaN, 'status', "");
n = height(T);
for b = 1:nBoot
    idx = randi(n, n, 1);
    Tb = T(idx, :);
    try
        Fit = fit_eo_dependent_reference_local(Tb, cfg);
        k = find(Fit.eo_ids == targetEO, 1);
        if isempty(k)
            rows(b).status = "missing_target_eo";
        else
            rows(b).K_target_eo_mm_per_microstrain = abs(Fit.C_eo(k));
            rows(b).phase_target_eo_rad = angle(Fit.C_eo(k));
            rows(b).rmse_mm = Fit.rmse_mm;
            rows(b).status = "ok";
        end
    catch
        rows(b).status = "failed";
    end
end
Boot = struct2table(rows);
end


function targetEO = infer_target_eo_from_points_local(T, cfg)
m = T.arrival_time_s >= cfg.step06i_target_time_range_s(1) & ...
    T.arrival_time_s <= cfg.step06i_target_time_range_s(2);
if any(m)
    targetEO = mode(T.dominant_order(m));
else
    targetEO = mode(T.dominant_order);
end
end


function Step05 = load_step05_for_comparison_local(files, cfg)
Step05 = struct('available', false, 'Trend', table());
if ~isfile(files.step05_trend_csv)
    return;
end
T = readtable(files.step05_trend_csv);
ok = strcmpi(string(T.status), "ok") & ...
    isfinite(T.A_id) & isfinite(T.EO_id) & isfinite(T.phi_id_wrapped) & ...
    T.window_center_time_s >= cfg.step06i_target_time_range_s(1) & ...
    T.window_center_time_s <= cfg.step06i_target_time_range_s(2);
Step05.available = any(ok);
Step05.Trend = T;
Step05.ok = ok;
end


function oprTimes = load_opr_times_local(file)
S = load(file, 'jiluOPR');
oprTimes = S.jiluOPR(:, 1);
oprTimes = oprTimes(:);
oprTimes = oprTimes(isfinite(oprTimes));
if any(diff(oprTimes) <= 0)
    oprTimes = sort(oprTimes);
end
end


function Tref = attach_step05_time_series_local(Tref, Step05, oprTimes, cfg)
Tref.u_step05_mm = nan(height(Tref), 1);
Tref.step05_A_mm = nan(height(Tref), 1);
Tref.step05_EO = nan(height(Tref), 1);
Tref.step05_fn_hz = nan(height(Tref), 1);
if ~Step05.available
    return;
end
T = Step05.Trend;
okIdx = find(Step05.ok);
centers = T.window_center_time_s(okIdx);
theta = map_time_to_rotor_phase_local(oprTimes, Tref.time_s, ...
    cfg.opr_events_per_revolution);
for i = 1:height(Tref)
    [~, j] = min(abs(centers - Tref.time_s(i)));
    row = okIdx(j);
    A = T.A_id(row);
    EO = T.EO_id(row);
    phi = T.phi_id_wrapped(row);
    Tref.u_step05_mm(i) = A .* sin(EO .* theta(i) + phi);
    Tref.step05_A_mm(i) = A;
    Tref.step05_EO(i) = EO;
    Tref.step05_fn_hz(i) = T.fn_id(row);
end
end


function W = build_window_comparison_local(Tref, Step05, cfg)
if ~Step05.available
    W = table();
    return;
end
T = Step05.Trend;
idx = find(Step05.ok);
rows = struct([]);
dt = median(diff(T.window_center_time_s(Step05.ok)), 'omitnan');
if ~isfinite(dt) || dt <= 0
    dt = 0.05;
end
for k = 1:numel(idx)
    i = idx(k);
    tr = T.window_center_time_s(i) + [-0.5 0.5] * dt;
    m = Tref.time_s >= tr(1) & Tref.time_s <= tr(2);
    uEO = Tref.u_ref_eo_mm(m);
    uGlobal = Tref.u_ref_globalK_allEO_mm(m);
    u05 = Tref.u_step05_mm(m);
    rows(k).window_id = T.window_id(i); %#ok<AGROW>
    rows(k).window_center_time_s = T.window_center_time_s(i);
    rows(k).Step05_A_mm = T.A_id(i);
    rows(k).Step05_EO = T.EO_id(i);
    rows(k).Step05_fn_hz = T.fn_id(i);
    rows(k).EO_reference_A_rms_mm = sqrt(2) * std(uEO, 'omitnan');
    rows(k).Global_reference_A_rms_mm = sqrt(2) * std(uGlobal, 'omitnan');
    rows(k).A_error_Step05_minus_EO_ref_mm = ...
        rows(k).Step05_A_mm - rows(k).EO_reference_A_rms_mm;
    rows(k).A_error_percent_of_EO_ref = ...
        100 * rows(k).A_error_Step05_minus_EO_ref_mm / ...
        rows(k).EO_reference_A_rms_mm;
    rows(k).time_rmse_EO_ref_mm = rmse_local(u05, uEO);
    if nnz(isfinite(u05) & isfinite(uEO)) >= 3
        rows(k).time_corr_EO_ref = corr(u05, uEO, 'Rows', 'complete');
    else
        rows(k).time_corr_EO_ref = NaN;
    end
end
W = struct2table(rows);
end


function Summary = build_summary_table_local(caseName, GlobalFit, EOT, ...
    Tref, FitPoints, Wcmp, CV, Boot, cfg)
targetMask = Tref.time_s >= cfg.step06i_target_time_range_s(1) & ...
    Tref.time_s <= cfg.step06i_target_time_range_s(2);
targetEO = mode(Tref.dominant_order(targetMask));
eoRow = EOT(EOT.dominant_order == targetEO, :);
okBoot = strcmp(string(Boot.status), "ok") & ...
    Boot.target_eo == targetEO & isfinite(Boot.K_target_eo_mm_per_microstrain);
Kci = [NaN NaN];
if any(okBoot)
    Kci = prctile(Boot.K_target_eo_mm_per_microstrain(okBoot), [2.5 97.5]);
end
status = string(eoRow.status(1));
if status == "ok"
    referenceStatus = "accepted_EO_dependent_reference";
else
    referenceStatus = "diagnostic_EO_dependent_reference";
end

Summary = table();
Summary.case_name = string(caseName);
Summary.reference_uses_step05 = false;
Summary.step05_usage = "comparison_only";
Summary.target_time_start_s = cfg.step06i_target_time_range_s(1);
Summary.target_time_end_s = cfg.step06i_target_time_range_s(2);
Summary.target_eo = targetEO;
Summary.target_eo_status = status;
Summary.reference_status = referenceStatus;
Summary.global_K_allEO_mm_per_microstrain = GlobalFit.K_mm_per_microstrain;
Summary.global_K_ratio_to_fixed = GlobalFit.K_ratio_to_fixed;
Summary.global_fit_rmse_mm = GlobalFit.rmse_mm;
Summary.global_phase_coverage_fraction = GlobalFit.phase_coverage_fraction;
Summary.target_EO_K_mm_per_microstrain = eoRow.K_mm_per_microstrain(1);
Summary.target_EO_K_ratio_to_fixed = eoRow.K_ratio_to_fixed(1);
Summary.target_EO_phase_lag_rad = eoRow.phase_lag_rad(1);
Summary.target_EO_phase_coverage_fraction = eoRow.phase_coverage_fraction(1);
Summary.target_EO_point_count = eoRow.point_count(1);
Summary.target_EO_boot_K_ci_low = Kci(1);
Summary.target_EO_boot_K_ci_high = Kci(2);
Summary.target_EO_rmse_mm = eoRow.rmse_mm(1);
Summary.fit_point_count_all = height(FitPoints);
Summary.eo_group_count = height(EOT);
Summary.target_EO_A_rms_mm = sqrt(2) * std(Tref.u_ref_eo_mm(targetMask), 'omitnan');
Summary.global_allEO_A_rms_mm = sqrt(2) * ...
    std(Tref.u_ref_globalK_allEO_mm(targetMask), 'omitnan');
Summary.fixedK_EOphase_A_rms_mm = sqrt(2) * ...
    std(Tref.u_fixedK_eo_phase_mm(targetMask), 'omitnan');
if ~isempty(Wcmp)
    Summary.step05_median_A_mm = median(Wcmp.Step05_A_mm, 'omitnan');
    Summary.step05_median_fn_hz = median(Wcmp.Step05_fn_hz, 'omitnan');
    Summary.step05_minus_EO_ref_A_percent = ...
        median(Wcmp.A_error_percent_of_EO_ref, 'omitnan');
    Summary.step05_vs_EO_ref_time_rmse_mm = ...
        rmse_local(Tref.u_step05_mm(targetMask), Tref.u_ref_eo_mm(targetMask));
    Summary.step05_vs_EO_ref_time_corr = corr( ...
        Tref.u_step05_mm(targetMask), Tref.u_ref_eo_mm(targetMask), ...
        'Rows', 'complete');
else
    Summary.step05_median_A_mm = NaN;
    Summary.step05_median_fn_hz = NaN;
    Summary.step05_minus_EO_ref_A_percent = NaN;
    Summary.step05_vs_EO_ref_time_rmse_mm = NaN;
    Summary.step05_vs_EO_ref_time_corr = NaN;
end
Summary.leave_one_eo_cv_max_rmse_mm = max(CV.rmse_test_mm, [], 'omitnan');
end


function figFile = plot_eo_reference_local(Tref, FitPoints, Summary, EOT, ...
    RegionTable, SensorResiduals, CV, Boot, figDir, cfg)
fig = figure('Visible', cfg.step06i_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 31, 25], ...
    'Name', 'Step06I EO-dependent joint reference');
tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

targetRange = cfg.step06i_target_time_range_s;
targetEO = Summary.target_eo(1);
targetMask = Tref.time_s >= targetRange(1) & ...
    Tref.time_s <= targetRange(2);
targetEOMask = targetMask & Tref.dominant_order == targetEO;
targetPointMask = FitPoints.arrival_time_s >= targetRange(1) & ...
    FitPoints.arrival_time_s <= targetRange(2) & ...
    FitPoints.dominant_order == targetEO;
zoomRange = cfg.step06i_zoom_time_range_s;
if ~any(FitPoints.arrival_time_s >= zoomRange(1) & ...
        FitPoints.arrival_time_s <= zoomRange(2) & ...
        FitPoints.dominant_order == targetEO)
    zoomRange = [targetRange(1), min(targetRange(2), targetRange(1) + 0.06)];
end

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
green = [0.466 0.674 0.188];
purple = [0.494 0.184 0.556];
lightGray = [0.640 0.640 0.640];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
eoVals = unique(RegionTable.dominant_order(:));
if isempty(eoVals)
    eoVals = targetEO;
end
yLim = [min(eoVals)-0.8, max(eoVals)+0.8];
patch(ax1, [targetRange(1) targetRange(2) targetRange(2) targetRange(1)], ...
    [yLim(1) yLim(1) yLim(2) yLim(2)], [0.90 0.95 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
for i = 1:height(RegionTable)
    t0 = RegionTable.time_start_s(i);
    t1 = RegionTable.time_end_s(i);
    eo = RegionTable.dominant_order(i);
    isTarget = eo == targetEO && t1 >= targetRange(1) && ...
        t0 <= targetRange(2);
    fc = 0.72*[1 1 1];
    ec = 0.45*[1 1 1];
    fa = 0.65;
    if isTarget
        fc = blue;
        ec = blue;
        fa = 0.78;
    end
    patch(ax1, [t0 t1 t1 t0], [eo-0.28 eo-0.28 eo+0.28 eo+0.28], ...
        fc, 'EdgeColor', ec, 'LineWidth', 0.7, 'FaceAlpha', fa, ...
        'HandleVisibility', 'off');
end
xline(ax1, targetRange(1), 'k:', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
xline(ax1, targetRange(2), 'k:', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
ylim(ax1, yLim);
xlim(ax1, [min(RegionTable.time_start_s)-0.5, ...
    max(RegionTable.time_end_s)+0.5]);
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'EO');
title(ax1, '1 Step06A selected EO regions');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
strainZoomMask = Tref.time_s >= zoomRange(1) & ...
    Tref.time_s <= zoomRange(2) & Tref.dominant_order == targetEO;
if any(strainZoomMask)
    hRe = Tref.strain_basis_real_microstrain(strainZoomMask);
    hAbs = Tref.strain_basis_abs_microstrain(strainZoomMask);
    plot_downsampled_local(ax2, Tref.time_s(strainZoomMask), hRe, blue, ...
        'Re\{h_{EO}(t)\}', 0.75);
    plot_downsampled_local(ax2, Tref.time_s(strainZoomMask), hAbs, ...
        orange, '|h_{EO}(t)|', 0.75);
    plot_downsampled_local(ax2, Tref.time_s(strainZoomMask), -hAbs, ...
        orange, '-|h_{EO}(t)|', 0.75);
end
xlim(ax2, zoomRange);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Strain basis (microstrain)');
title(ax2, '2 Strain waveform basis h_{EO}');
legend(ax2, 'Location', 'best');

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
tr = zoomRange;
mZoom = Tref.time_s >= tr(1) & Tref.time_s <= tr(2) & ...
    Tref.dominant_order == targetEO;
plot(ax3, Tref.time_s(mZoom), Tref.u_ref_eo_mm(mZoom), ...
    '-', 'Color', blue, 'LineWidth', 1.15, ...
    'DisplayName', 'EO reference');
mpZoom = FitPoints.arrival_time_s >= tr(1) & ...
    FitPoints.arrival_time_s <= tr(2) & ...
    FitPoints.dominant_order == targetEO;
sids = unique(FitPoints.sensor_id(mpZoom)).';
colors = lines(max(numel(sids), 3));
for i = 1:numel(sids)
    sid = sids(i);
    ms = mpZoom & FitPoints.sensor_id == sid;
    scatter(ax3, FitPoints.arrival_time_s(ms), ...
        FitPoints.displacement_dynamic_eo_mm(ms), 20, colors(i, :), ...
        'filled', 'MarkerFaceAlpha', 0.72, ...
        'DisplayName', sprintf('BTT S%d', sid));
end
xlim(ax3, tr);
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Dynamic displacement (mm)');
title(ax3, '3 Sparse BTT anchors scale h_{EO}');
legend(ax3, 'Location', 'best');

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
mp = targetPointMask & isfinite(FitPoints.displacement_dynamic_eo_mm) & ...
    isfinite(FitPoints.u_eo_ref_mm);
scatter(ax4, FitPoints.u_eo_ref_mm(mp), ...
    FitPoints.displacement_dynamic_eo_mm(mp), 12, blue, 'filled', ...
    'MarkerFaceAlpha', 0.32, 'DisplayName', 'BTT points');
if any(mp)
    xy = [FitPoints.u_eo_ref_mm(mp); ...
        FitPoints.displacement_dynamic_eo_mm(mp)];
    lim = [min(xy, [], 'omitnan'), max(xy, [], 'omitnan')];
    pad = 0.06 * max(diff(lim), eps);
    lim = lim + [-pad pad];
    plot(ax4, lim, lim, 'k--', 'LineWidth', 0.9, ...
        'DisplayName', 'y=x');
    xlim(ax4, lim);
    ylim(ax4, lim);
    axis(ax4, 'square');
    rFit = corr(FitPoints.u_eo_ref_mm(mp), ...
        FitPoints.displacement_dynamic_eo_mm(mp), 'Rows', 'complete');
    text(ax4, 0.04, 0.94, sprintf('N=%d, RMSE=%.3f mm, r=%.3f', ...
        nnz(mp), Summary.target_EO_rmse_mm(1), rFit), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 8.2, 'BackgroundColor', 'w', 'Margin', 1.5);
end
xlabel(ax4, 'Re\{C_{EO}h(t_i)\} (mm)');
ylabel(ax4, 'BTT, offset removed (mm)');
title(ax4, '4 Fit C_{EO} using BTT points');
legend(ax4, 'Location', 'southeast');

ax5 = nexttile; hold(ax5, 'on'); box(ax5, 'on');
x = 1:height(EOT);
b = bar(ax5, x, 1000 * EOT.K_mm_per_microstrain, 0.68, ...
    'DisplayName', 'K_{EO}');
b.FaceColor = 'flat';
b.CData = repmat(0.70*[1 1 1], height(EOT), 1);
iTarget = find(EOT.dominant_order == targetEO, 1);
if ~isempty(iTarget)
    b.CData(iTarget, :) = purple;
    k0 = 1000 * Summary.target_EO_K_mm_per_microstrain(1);
    kLow = 1000 * Summary.target_EO_boot_K_ci_low(1);
    kHigh = 1000 * Summary.target_EO_boot_K_ci_high(1);
    if all(isfinite([k0, kLow, kHigh]))
        errorbar(ax5, iTarget, k0, k0-kLow, kHigh-k0, 'k.', ...
            'LineWidth', 1.1, 'CapSize', 8, ...
            'DisplayName', 'target 95% CI');
    end
end
yline(ax5, 1000 * Summary.global_K_allEO_mm_per_microstrain(1), ...
    'b-', 'LineWidth', 1.0, 'DisplayName', 'global K');
yline(ax5, 1000 * cfg.step06i_fixed_strain_to_mm, ...
    'g--', 'LineWidth', 1.0, 'DisplayName', 'fixed K');
set(ax5, 'XTick', x, 'XTickLabel', string(EOT.dominant_order));
xlabel(ax5, 'Engine order');
ylabel(ax5, 'K (\mum/microstrain)');
title(ax5, '5 EO-dependent strain-to-tip gain');
legend(ax5, 'Location', 'best');

ax6 = nexttile; hold(ax6, 'on'); box(ax6, 'on');
okBoot = strcmp(string(Boot.status), "ok") & Boot.target_eo == targetEO & ...
    isfinite(Boot.K_target_eo_mm_per_microstrain);
if any(okBoot)
    histogram(ax6, 1000 * Boot.K_target_eo_mm_per_microstrain(okBoot), 18, ...
        'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none', ...
        'DisplayName', 'bootstrap K');
    xline(ax6, 1000 * Summary.target_EO_K_mm_per_microstrain(1), ...
        'b-', 'LineWidth', 1.1, 'DisplayName', 'target K');
    xline(ax6, 1000 * Summary.target_EO_boot_K_ci_low(1), ...
        'b:', 'LineWidth', 0.9, 'HandleVisibility', 'off');
    xline(ax6, 1000 * Summary.target_EO_boot_K_ci_high(1), ...
        'b:', 'LineWidth', 0.9, 'HandleVisibility', 'off');
    xline(ax6, 1000 * cfg.step06i_fixed_strain_to_mm, ...
        'g--', 'LineWidth', 1.0, 'DisplayName', 'fixed K');
    xlabel(ax6, 'K (\mum/microstrain)');
    ylabel(ax6, 'Count');
    text(ax6, 0.04, 0.94, sprintf('max CV RMSE = %.3f mm', ...
        Summary.leave_one_eo_cv_max_rmse_mm(1)), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 8.0, 'BackgroundColor', 'w', 'Margin', 1.5);
    legend(ax6, 'Location', 'best');
else
    text(ax6, 0.5, 0.5, 'No valid bootstrap samples', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax6, 'off');
end
title(ax6, '6 Bootstrap stability of target K');

ax7 = nexttile; hold(ax7, 'on'); box(ax7, 'on');
if any(isfinite(Tref.u_step05_mm(targetMask)))
    plot_downsampled_local(ax7, Tref.time_s(targetMask), ...
        Tref.u_step05_mm(targetMask), lightGray, ...
        'Step05 full waveform', 0.35);
end
plot_downsampled_local(ax7, Tref.time_s(targetMask), ...
    Tref.u_ref_eo_mm(targetMask), blue, 'EO reference', 0.95);
plot_downsampled_local(ax7, Tref.time_s(targetMask), ...
    Tref.u_ref_globalK_allEO_mm(targetMask), orange, ...
    'global K all EO', 0.70);
plot_downsampled_local(ax7, Tref.time_s(targetMask), ...
    Tref.u_fixedK_eo_phase_mm(targetMask), green, ...
    'fixed K, EO phase', 0.70);
xlim(ax7, targetRange);
xlabel(ax7, 'Time (s)');
ylabel(ax7, 'u_{tip} (mm)');
title(ax7, '7 Complete target-window displacement');
legend(ax7, 'Location', 'best');

ax8 = nexttile; hold(ax8, 'on'); box(ax8, 'on');
bar(ax8, SensorResiduals.sensor_id, ...
    [SensorResiduals.sensor_bias_eo_mm, SensorResiduals.mean_residual_mm], ...
    0.72);
yline(ax8, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
text(ax8, 0.04, 0.94, sprintf('max leave-one-EO RMSE = %.3f mm', ...
    Summary.leave_one_eo_cv_max_rmse_mm(1)), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8.0, 'BackgroundColor', 'w', 'Margin', 1.5);
xlabel(ax8, 'Sensor id');
ylabel(ax8, 'Residual quantity (mm)');
title(ax8, '8 Sensor-offset residual check');
legend(ax8, {'fitted sensor bias', 'mean residual'}, 'Location', 'best');

sg = sgtitle(fig, sprintf(['Step06I: regions -> h_{EO}(t) -> ', ...
    'BTT-fit C_{EO} -> u_{ref}(t)\n', ...
    'Target %.2f-%.2f s, EO%d, K=%.3f \\mum/microstrain, A=%.3f mm'], ...
    Summary.target_time_start_s(1), Summary.target_time_end_s(1), ...
    targetEO, 1000 * Summary.target_EO_K_mm_per_microstrain(1), ...
    Summary.target_EO_A_rms_mm(1)));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.8);

style_step06i_figure_local(fig);
figFile = fullfile(figDir, ...
    'Step06I_EODependentReference_20250527.png');
save_figure_local(fig, figFile);
end


function plot_downsampled_local(ax, x, y, color, name, lw)
n = numel(x);
if n == 0
    return;
end
step = max(1, ceil(n / 30000));
plot(ax, x(1:step:end), y(1:step:end), '-', ...
    'Color', color, 'LineWidth', lw, 'DisplayName', name);
end


function style_step06i_figure_local(fig)
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    set(axs(i), 'FontName', 'Times New Roman', 'FontSize', 8.4, ...
        'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.75, ...
        'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.13);
end
txt = findall(fig, 'Type', 'text');
set(txt, 'FontName', 'Times New Roman');
leg = findall(fig, 'Type', 'legend');
for i = 1:numel(leg)
    set(leg(i), 'FontName', 'Times New Roman', 'FontSize', 7.3, ...
        'Box', 'off');
end
end


function thetaRot = map_time_to_rotor_phase_local(oprTimes, sampleTimes, oprEventsPerRev)
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


function c = phase_coverage_local(phi)
phi = phi(isfinite(phi));
if numel(phi) < 2
    c = 0;
    return;
end
p = sort(mod(phi(:), 2*pi));
gaps = diff([p; p(1) + 2*pi]);
c = 2*pi - max(gaps);
end


function v = rmse_local(y, yFit)
valid = isfinite(y) & isfinite(yFit);
if ~any(valid)
    v = NaN;
else
    v = sqrt(mean((y(valid) - yFit(valid)).^2, 'omitnan'));
end
end


function r2 = compute_r2_local(y, yFit)
valid = isfinite(y) & isfinite(yFit);
if nnz(valid) < 2
    r2 = NaN;
    return;
end
y = y(valid);
yFit = yFit(valid);
ssRes = sum((y - yFit).^2);
ssTot = sum((y - mean(y, 'omitnan')).^2);
if ssTot <= eps
    r2 = NaN;
else
    r2 = 1 - ssRes / ssTot;
end
end


function ensure_dir_local(p)
if ~isfolder(p)
    mkdir(p);
end
end


function save_figure_local(fig, pngFile)
[folder, base, ~] = fileparts(pngFile);
ensure_dir_local(folder);
pngFile = fullfile(folder, [base '.png']);
pdfFile = fullfile(folder, [base '.pdf']);
figFile = fullfile(folder, [base '.fig']);
try
    exportgraphics(fig, pngFile, 'Resolution', 220);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
catch
    print(fig, pngFile, '-dpng', '-r220');
    print(fig, pdfFile, '-dpdf', '-painters');
end
savefig(fig, figFile);
end
