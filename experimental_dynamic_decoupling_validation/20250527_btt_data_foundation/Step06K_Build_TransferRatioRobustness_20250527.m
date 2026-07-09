%% Step06K_Build_TransferRatioRobustness_20250527.m
% Step06K: constrain K to the FE/APDL transfer-ratio interval,
% then use BTT only to choose the in-band phase/gain correction.

clc; clear; close all;

dataset = '20250527';
cfg = BTTDataConfig_20250527();

%RUN_STEP06K_TRANSFERRATIOROBUSTNESS
% Robustness layer for the Step06I strain-transfer reference.  The main
% Step06K route constrains K to the FE/APDL band and uses BTT only to choose
% the in-band phase/gain correction.

dataset = char(dataset);
cfg = apply_step06k_transfer_defaults_local(cfg);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06K transfer-ratio robustness: %s / %s ===\n', ...
        dataset, caseName);

    files = resolve_step06k_transfer_files_local(cfg, dataset, caseName);
    assert_step06k_transfer_files_local(files);
    D = load(files.step06i_mat, 'Summary', 'TransferPrior', ...
        'RegionSummary', 'ReferenceTimeSeries', 'WindowComparison');

    RobustRegion = build_robust_region_table_local(D.RegionSummary, ...
        D.TransferPrior, cfg);
    BandTimeSeries = build_band_time_series_local(D.ReferenceTimeSeries, ...
        D.TransferPrior);
    BTTCorrectionPoints = build_btt_correction_points_local( ...
        BandTimeSeries, D.RegionSummary, files, cfg);
    [KConstrainedFit, KGridSearch] = fit_fe_band_btt_reference_local( ...
        BTTCorrectionPoints, BandTimeSeries, RobustRegion, cfg);
    BTTCorrectedTimeSeries = build_btt_corrected_time_series_local( ...
        BandTimeSeries, KConstrainedFit);
    Step05BandCheck = build_step05_band_check_local(D.WindowComparison, ...
        RobustRegion, BTTCorrectedTimeSeries, KConstrainedFit);
    Summary = build_step06k_transfer_summary_local(caseName, dataset, ...
        RobustRegion, Step05BandCheck, KConstrainedFit, ...
        BTTCorrectionPoints, files);

    outDir = fullfile(cfg.step06k_transfer_output_dir, caseName);
    figDir = fullfile(cfg.step06k_transfer_figure_dir, caseName);
    ensure_dir_local(outDir);
    ensure_dir_local(figDir);

    summaryCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_Summary_%s.csv', dataset));
    regionCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_RegionBand_%s.csv', dataset));
    tsCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_TimeSeriesBand_%s.csv', dataset));
    bandCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_Step05BandCheck_%s.csv', dataset));
    bttCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_BTTCorrectionPoints_%s.csv', dataset));
    fitCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_KConstrainedFit_%s.csv', dataset));
    gridCsv = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_KGridSearch_%s.csv', dataset));
    matFile = fullfile(outDir, sprintf( ...
        'Step06K_TransferRatioRobustness_%s.mat', dataset));

    writetable(Summary, summaryCsv);
    writetable(RobustRegion, regionCsv);
    writetable(BTTCorrectedTimeSeries, tsCsv);
    writetable(Step05BandCheck, bandCsv);
    writetable(BTTCorrectionPoints, bttCsv);
    writetable(KConstrainedFit, fitCsv);
    writetable(KGridSearch, gridCsv);
    save(matFile, 'Summary', 'RobustRegion', 'BandTimeSeries', ...
        'BTTCorrectedTimeSeries', 'BTTCorrectionPoints', ...
        'KConstrainedFit', 'KGridSearch', 'Step05BandCheck', ...
        'files', 'cfg', '-v7.3');

    figFile = plot_step06k_transfer_robustness_local(BTTCorrectedTimeSeries, ...
        RobustRegion, Step05BandCheck, KConstrainedFit, ...
        BTTCorrectionPoints, Summary, figDir, dataset, cfg);

    fprintf('\nStep06K transfer-ratio robustness summary:\n');
    disp(Summary);
    fprintf(['Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        '  %s\n  %s\n  %s\n'], ...
        summaryCsv, regionCsv, tsCsv, bandCsv, bttCsv, fitCsv, ...
        gridCsv, matFile, figFile);
end

function cfg = apply_step06k_transfer_defaults_local(cfg)
if ~isfield(cfg, 'output_root') || isempty(cfg.output_root)
    cfg.output_root = fullfile(pwd, 'output');
end
if ~isfield(cfg, 'figure_root') || isempty(cfg.figure_root)
    cfg.figure_root = fullfile(cfg.output_root, 'figures');
end
if ~isfield(cfg, 'step06i_transfer_output_dir') || ...
        isempty(cfg.step06i_transfer_output_dir)
    cfg.step06i_transfer_output_dir = fullfile(cfg.output_root, ...
        'step06i_strain_transfer_displacement_reference');
end
if ~isfield(cfg, 'step06k_transfer_output_dir') || ...
        isempty(cfg.step06k_transfer_output_dir)
    cfg.step06k_transfer_output_dir = fullfile(cfg.output_root, ...
        'step06k_transfer_ratio_robustness');
end
if ~isfield(cfg, 'step06k_transfer_figure_dir') || ...
        isempty(cfg.step06k_transfer_figure_dir)
    cfg.step06k_transfer_figure_dir = fullfile(cfg.figure_root, ...
        'step06k_transfer_ratio_robustness');
end
if ~isfield(cfg, 'step06k_transfer_fig_visible') || ...
        isempty(cfg.step06k_transfer_fig_visible)
    if usejava('desktop')
        cfg.step06k_transfer_fig_visible = 'on';
    else
        cfg.step06k_transfer_fig_visible = 'off';
    end
end
if ~isfield(cfg, 'step06k_target_blade') || isempty(cfg.step06k_target_blade)
    if isfield(cfg, 'step06i_target_blade') && ~isempty(cfg.step06i_target_blade)
        cfg.step06k_target_blade = cfg.step06i_target_blade;
    elseif isfield(cfg, 'step05_target_blades') && ~isempty(cfg.step05_target_blades)
        cfg.step06k_target_blade = cfg.step05_target_blades(1);
    else
        cfg.step06k_target_blade = 1;
    end
end
if ~isfield(cfg, 'step06k_sensor_ids') || isempty(cfg.step06k_sensor_ids)
    if isfield(cfg, 'optical_fiber_ids') && ~isempty(cfg.optical_fiber_ids)
        cfg.step06k_sensor_ids = cfg.optical_fiber_ids(:).';
    elseif isfield(cfg, 'sensor_ids') && ~isempty(cfg.sensor_ids)
        cfg.step06k_sensor_ids = cfg.sensor_ids(:).';
    else
        cfg.step06k_sensor_ids = [];
    end
end
if ~isfield(cfg, 'step06k_gauge_position_uncertainty_mm') || ...
        isempty(cfg.step06k_gauge_position_uncertainty_mm)
    cfg.step06k_gauge_position_uncertainty_mm = 2.5;
end
if ~isfield(cfg, 'step06k_K_grid_count') || isempty(cfg.step06k_K_grid_count)
    cfg.step06k_K_grid_count = 31;
end
if ~isfield(cfg, 'step06k_phase_grid_count') || ...
        isempty(cfg.step06k_phase_grid_count)
    cfg.step06k_phase_grid_count = 73;
end
if ~isfield(cfg, 'step06k_btt_correction_use_core_only') || ...
        isempty(cfg.step06k_btt_correction_use_core_only)
    cfg.step06k_btt_correction_use_core_only = true;
end
if ~isfield(cfg, 'step06k_min_btt_points_per_region') || ...
        isempty(cfg.step06k_min_btt_points_per_region)
    cfg.step06k_min_btt_points_per_region = 20;
end
end


function files = resolve_step06k_transfer_files_local(cfg, dataset, caseName)
files = struct();
files.step06i_mat = fullfile(cfg.step06i_transfer_output_dir, caseName, ...
    sprintf('Step06I_StrainTransferDisplacementReference_%s.mat', dataset));
files.step03_bundle_mat = fullfile(cfg.step03_output_dir, caseName, ...
    sprintf('BTT_Observation_Bundle_%s.mat', dataset));
end


function assert_step06k_transfer_files_local(files)
if ~isfile(files.step06i_mat)
    error(['Missing Step06I strain-transfer MAT:\n  %s\n', ...
        'Run Step06I_Build_StrainTransferDisplacementReference first.'], ...
        files.step06i_mat);
end
if ~isfile(files.step03_bundle_mat)
    error('Missing Step03 BTT observation bundle:\n  %s', files.step03_bundle_mat);
end
end


function Rb = build_robust_region_table_local(RS, TP, cfg)
Rb = RS;
Rb.K_nominal_um_per_microstrain = RS.K_nominal_um_per_microstrain;
Rb.K_low_um_per_microstrain = RS.K_band_low_um_per_microstrain;
Rb.K_high_um_per_microstrain = RS.K_band_high_um_per_microstrain;
if isfield(cfg, 'step06k_K_low_um_per_microstrain') && ...
        ~isempty(cfg.step06k_K_low_um_per_microstrain)
    Rb.K_low_um_per_microstrain(:) = cfg.step06k_K_low_um_per_microstrain;
end
if isfield(cfg, 'step06k_K_high_um_per_microstrain') && ...
        ~isempty(cfg.step06k_K_high_um_per_microstrain)
    Rb.K_high_um_per_microstrain(:) = cfg.step06k_K_high_um_per_microstrain;
end
scaleLow = Rb.K_low_um_per_microstrain ./ Rb.K_nominal_um_per_microstrain;
scaleHigh = Rb.K_high_um_per_microstrain ./ Rb.K_nominal_um_per_microstrain;
Rb.A_peak_equiv_low_mm = RS.u_nominal_A_peak_equiv_mm .* scaleLow;
Rb.A_peak_equiv_high_mm = RS.u_nominal_A_peak_equiv_mm .* scaleHigh;
Rb.envelope_median_low_mm = RS.u_abs_envelope_median_mm .* scaleLow;
Rb.envelope_median_high_mm = RS.u_abs_envelope_median_mm .* scaleHigh;
Rb.transfer_uncertainty_definition = repmat(string(sprintf( ...
    ['FE K constrained band; gauge neighborhood uncertainty ', ...
    'about %.2f mm, using Step06I K_low/K_high unless explicitly overridden'], ...
    cfg.step06k_gauge_position_uncertainty_mm)), height(Rb), 1);
Rb.gauge_position_uncertainty_mm = repmat( ...
    cfg.step06k_gauge_position_uncertainty_mm, height(Rb), 1);
if ~isempty(TP) && ismember('source_note', TP.Properties.VariableNames)
    Rb.transfer_source_note = TP.source_note(1:min(height(TP), height(Rb)));
end
end


function TSb = build_band_time_series_local(TS, TP)
TSb = TS(:, {'time_s', 'region_id', 'dominant_order', 'peak_freq_hz', ...
    'strain_basis_real_microstrain', 'strain_basis_imag_microstrain', ...
    'strain_basis_abs_microstrain', ...
    'u_transfer_nominal_mm', 'u_transfer_abs_envelope_mm', ...
    'u_step05_mm', 'step05_A_mm', 'step05_EO'});
TSb.K_low_mm_per_microstrain = nan(height(TSb), 1);
TSb.K_high_mm_per_microstrain = nan(height(TSb), 1);
for i = 1:height(TP)
    m = TSb.region_id == TP.region_id(i);
    TSb.K_low_mm_per_microstrain(m) = TP.K_band_low_um_per_microstrain(i) / 1000;
    TSb.K_high_mm_per_microstrain(m) = TP.K_band_high_um_per_microstrain(i) / 1000;
end
TSb.u_signed_low_mm = TSb.K_low_mm_per_microstrain .* ...
    TSb.strain_basis_real_microstrain;
TSb.u_signed_high_mm = TSb.K_high_mm_per_microstrain .* ...
    TSb.strain_basis_real_microstrain;
TSb.u_abs_low_mm = TSb.K_low_mm_per_microstrain .* ...
    TSb.strain_basis_abs_microstrain;
TSb.u_abs_high_mm = TSb.K_high_mm_per_microstrain .* ...
    TSb.strain_basis_abs_microstrain;
end


function P = build_btt_correction_points_local(TS, RS, files, cfg)
Obs = load_step03_observation_table_local(files.step03_bundle_mat);
P = table();
if isempty(Obs)
    return;
end
need = {'sensor_id', 'blade_id', 'arrival_time_s', 'displacement_mm'};
for i = 1:numel(need)
    if ~ismember(need{i}, Obs.Properties.VariableNames)
        error('Step03 observation table is missing %s:\n  %s', ...
            need{i}, files.step03_bundle_mat);
    end
end
valid = Obs.blade_id == cfg.step06k_target_blade & ...
    ismember(Obs.sensor_id, cfg.step06k_sensor_ids(:).') & ...
    isfinite(Obs.arrival_time_s) & isfinite(Obs.displacement_mm);
if ismember('is_valid', Obs.Properties.VariableNames)
    valid = valid & logical(Obs.is_valid);
end
Obs = Obs(valid, :);
if isempty(Obs)
    return;
end
rows = cell(height(RS), 1);
for i = 1:height(RS)
    rid = RS.region_id(i);
    if cfg.step06k_btt_correction_use_core_only && ...
            ismember('core_time_start_s', RS.Properties.VariableNames) && ...
            isfinite(RS.core_time_start_s(i)) && isfinite(RS.core_time_end_s(i))
        t0 = RS.core_time_start_s(i);
        t1 = RS.core_time_end_s(i);
    else
        t0 = RS.time_start_s(i);
        t1 = RS.time_end_s(i);
    end
    mObs = Obs.arrival_time_s >= t0 & Obs.arrival_time_s <= t1;
    mTS = TS.region_id == rid & TS.time_s >= t0 & TS.time_s <= t1;
    if nnz(mObs) < 1 || nnz(mTS) < 2
        continue;
    end
    O = Obs(mObs, :);
    [tu, ia] = unique(TS.time_s(mTS));
    hReal = TS.strain_basis_real_microstrain(mTS);
    hImag = zeros(size(hReal));
    if ismember('strain_basis_imag_microstrain', TS.Properties.VariableNames)
        hImag = TS.strain_basis_imag_microstrain(mTS);
    end
    hAbs = TS.strain_basis_abs_microstrain(mTS);
    hReal = hReal(ia);
    hImag = hImag(ia);
    hAbs = hAbs(ia);
    O.region_id = repmat(rid, height(O), 1);
    O.dominant_order = repmat(RS.dominant_order(i), height(O), 1);
    O.strain_basis_real_microstrain = interp1(tu, hReal, ...
        O.arrival_time_s, 'linear', NaN);
    O.strain_basis_imag_microstrain = interp1(tu, hImag, ...
        O.arrival_time_s, 'linear', NaN);
    O.strain_basis_abs_microstrain = interp1(tu, hAbs, ...
        O.arrival_time_s, 'linear', NaN);
    O.core_time_start_s = repmat(t0, height(O), 1);
    O.core_time_end_s = repmat(t1, height(O), 1);
    rows{i} = O;
end
rows = rows(~cellfun('isempty', rows));
if isempty(rows)
    return;
end
P = vertcat(rows{:});
ok = isfinite(P.strain_basis_real_microstrain) & ...
    isfinite(P.strain_basis_imag_microstrain);
P = sortrows(P(ok, :), {'region_id', 'arrival_time_s', 'sensor_id'});
end


function Obs = load_step03_observation_table_local(file)
S = load(file);
if isfield(S, 'bundle')
    B = S.bundle;
elseif isfield(S, 'BTT_Observation_Bundle')
    B = S.BTT_Observation_Bundle;
else
    B = struct();
end
if isfield(B, 'Observation_Table')
    Obs = B.Observation_Table;
elseif isfield(B, 'ObservationTable')
    Obs = B.ObservationTable;
elseif isfield(B, 'observation_table')
    Obs = B.observation_table;
else
    Obs = table();
end
end


function [Fit, Grid] = fit_fe_band_btt_reference_local(P, TS, Rb, cfg)
rows = repmat(empty_k_fit_row_local(), height(Rb), 1);
gridCells = cell(height(Rb), 1);
for i = 1:height(Rb)
    rid = Rb.region_id(i);
    Tp = P(P.region_id == rid, :);
    rows(i).region_id = rid;
    rows(i).dominant_order = Rb.dominant_order(i);
    rows(i).K_low_um_per_microstrain = Rb.K_low_um_per_microstrain(i);
    rows(i).K_high_um_per_microstrain = Rb.K_high_um_per_microstrain(i);
    rows(i).K_nominal_um_per_microstrain = Rb.K_nominal_um_per_microstrain(i);
    rows(i).point_count = height(Tp);
    if height(Tp) < cfg.step06k_min_btt_points_per_region
        rows(i).status = "too_few_BTT_points";
        continue;
    end
    kGrid = linspace(Rb.K_low_um_per_microstrain(i), ...
        Rb.K_high_um_per_microstrain(i), cfg.step06k_K_grid_count) ./ 1000;
    phaseGrid = linspace(-pi, pi, cfg.step06k_phase_grid_count);
    phaseGrid(end) = [];
    [best, G] = scan_k_phase_grid_local(Tp, kGrid, phaseGrid);
    rows(i).K_um_per_microstrain = 1000 * best.K_mm_per_microstrain;
    rows(i).K_mm_per_microstrain = best.K_mm_per_microstrain;
    rows(i).phase_lag_rad = best.phase_lag_rad;
    rows(i).phase_lag_deg = 180 * best.phase_lag_rad / pi;
    rows(i).offset_removed_rmse_mm = best.rmse_mm;
    rows(i).mean_abs_residual_mm = best.mean_abs_residual_mm;
    rows(i).sensor_offset_span_mm = best.sensor_offset_span_mm;
    rows(i).constraint_status = best.constraint_status;
    rows(i).C_real_mm_per_microstrain = best.C_real_mm_per_microstrain;
    rows(i).C_imag_mm_per_microstrain = best.C_imag_mm_per_microstrain;
    rows(i).status = "ok";
    mt = TS.region_id == rid;
    u = best.K_mm_per_microstrain .* ...
        (cos(best.phase_lag_rad) .* TS.strain_basis_real_microstrain(mt) - ...
        sin(best.phase_lag_rad) .* TS.strain_basis_imag_microstrain(mt));
    rows(i).corrected_A_peak_equiv_mm = sqrt(2) * std(u, 'omitnan');
    rows(i).corrected_abs_envelope_median_mm = best.K_mm_per_microstrain .* ...
        median(TS.strain_basis_abs_microstrain(mt), 'omitnan');
    gridCells{i} = G;
end
Fit = struct2table(rows);
Grid = vertcat(gridCells{~cellfun('isempty', gridCells)});
if isempty(Grid)
    Grid = table();
end
end


function [best, G] = scan_k_phase_grid_local(T, kGrid, phaseGrid)
row = 0;
gridRows = repmat(empty_grid_row_local(), numel(kGrid) * numel(phaseGrid), 1);
best = empty_grid_row_local();
best.rmse_mm = inf;
y = T.displacement_mm(:);
hRe = T.strain_basis_real_microstrain(:);
hIm = T.strain_basis_imag_microstrain(:);
sensorIds = unique(T.sensor_id(:)).';
for iK = 1:numel(kGrid)
    kval = kGrid(iK);
    for iP = 1:numel(phaseGrid)
        phi = phaseGrid(iP);
        u = kval .* (cos(phi) .* hRe - sin(phi) .* hIm);
        offsets = nan(size(sensorIds));
        yFit = u;
        for is = 1:numel(sensorIds)
            m = T.sensor_id == sensorIds(is);
            offsets(is) = median(y(m) - u(m), 'omitnan');
            yFit(m) = yFit(m) + offsets(is);
        end
        r = y - yFit;
        rmse = sqrt(mean(r.^2, 'omitnan'));
        row = row + 1;
        gridRows(row).region_id = mode(T.region_id);
        gridRows(row).K_um_per_microstrain = 1000 * kval;
        gridRows(row).phase_lag_rad = phi;
        gridRows(row).rmse_mm = rmse;
        gridRows(row).mean_abs_residual_mm = mean(abs(r), 'omitnan');
        if rmse < best.rmse_mm
            best = gridRows(row);
            best.K_mm_per_microstrain = kval;
            best.sensor_offset_span_mm = max(offsets, [], 'omitnan') - ...
                min(offsets, [], 'omitnan');
            best.C_real_mm_per_microstrain = kval * cos(phi);
            best.C_imag_mm_per_microstrain = kval * sin(phi);
        end
    end
end
gridRows = gridRows(1:row);
G = struct2table(gridRows);
tol = 0.5 * min(abs(diff(kGrid)), [], 'omitnan');
if ~isfinite(tol)
    tol = 0;
end
if abs(best.K_mm_per_microstrain - min(kGrid)) <= tol
    best.constraint_status = "at_lower_K_bound";
elseif abs(best.K_mm_per_microstrain - max(kGrid)) <= tol
    best.constraint_status = "at_upper_K_bound";
else
    best.constraint_status = "inside_K_band";
end
end


function TS = build_btt_corrected_time_series_local(TS, Fit)
TS.u_fe_btt_corrected_mm = nan(height(TS), 1);
TS.u_fe_btt_corrected_abs_envelope_mm = nan(height(TS), 1);
TS.K_fe_btt_um_per_microstrain = nan(height(TS), 1);
TS.phase_fe_btt_rad = nan(height(TS), 1);
for i = 1:height(Fit)
    if ~strcmp(string(Fit.status(i)), "ok")
        continue;
    end
    m = TS.region_id == Fit.region_id(i);
    k = Fit.K_mm_per_microstrain(i);
    phi = Fit.phase_lag_rad(i);
    TS.u_fe_btt_corrected_mm(m) = k .* ...
        (cos(phi) .* TS.strain_basis_real_microstrain(m) - ...
        sin(phi) .* TS.strain_basis_imag_microstrain(m));
    TS.u_fe_btt_corrected_abs_envelope_mm(m) = k .* ...
        TS.strain_basis_abs_microstrain(m);
    TS.K_fe_btt_um_per_microstrain(m) = Fit.K_um_per_microstrain(i);
    TS.phase_fe_btt_rad(m) = phi;
end
end


function row = empty_k_fit_row_local()
row = struct('region_id', NaN, 'dominant_order', NaN, ...
    'point_count', NaN, 'K_low_um_per_microstrain', NaN, ...
    'K_high_um_per_microstrain', NaN, ...
    'K_nominal_um_per_microstrain', NaN, ...
    'K_um_per_microstrain', NaN, 'K_mm_per_microstrain', NaN, ...
    'phase_lag_rad', NaN, 'phase_lag_deg', NaN, ...
    'C_real_mm_per_microstrain', NaN, ...
    'C_imag_mm_per_microstrain', NaN, ...
    'corrected_A_peak_equiv_mm', NaN, ...
    'corrected_abs_envelope_median_mm', NaN, ...
    'offset_removed_rmse_mm', NaN, ...
    'mean_abs_residual_mm', NaN, ...
    'sensor_offset_span_mm', NaN, ...
    'constraint_status', "", 'status', "");
end


function row = empty_grid_row_local()
row = struct('region_id', NaN, 'K_um_per_microstrain', NaN, ...
    'K_mm_per_microstrain', NaN, 'phase_lag_rad', NaN, ...
    'rmse_mm', NaN, 'mean_abs_residual_mm', NaN, ...
    'sensor_offset_span_mm', NaN, ...
    'C_real_mm_per_microstrain', NaN, ...
    'C_imag_mm_per_microstrain', NaN, ...
    'constraint_status', "");
end


function C = build_step05_band_check_local(W, Rb, TS, Fit)
if isempty(W)
    C = table();
    return;
end
C = W;
C.reference_A_low_mm = nan(height(W), 1);
C.reference_A_high_mm = nan(height(W), 1);
C.reference_envelope_low_mm = nan(height(W), 1);
C.reference_envelope_high_mm = nan(height(W), 1);
C.FE_BTT_corrected_A_mm = nan(height(W), 1);
C.FE_BTT_corrected_envelope_median_mm = nan(height(W), 1);
C.Step05_vs_FE_BTT_rmse_mm = nan(height(W), 1);
C.Step05_vs_FE_BTT_corr = nan(height(W), 1);
for i = 1:height(W)
    rid = W.reference_region_id(i);
    j = find(Rb.region_id == rid, 1);
    if isempty(j)
        continue;
    end
    C.reference_A_low_mm(i) = Rb.A_peak_equiv_low_mm(j);
    C.reference_A_high_mm(i) = Rb.A_peak_equiv_high_mm(j);
    C.reference_envelope_low_mm(i) = Rb.envelope_median_low_mm(j);
    C.reference_envelope_high_mm(i) = Rb.envelope_median_high_mm(j);
    jf = find(Fit.region_id == rid & strcmp(string(Fit.status), "ok"), 1);
    if ~isempty(jf)
        C.FE_BTT_corrected_A_mm(i) = Fit.corrected_A_peak_equiv_mm(jf);
        C.FE_BTT_corrected_envelope_median_mm(i) = ...
            Fit.corrected_abs_envelope_median_mm(jf);
    end
    if ismember('time_start_s', W.Properties.VariableNames) && ...
            isfinite(W.time_start_s(i)) && isfinite(W.time_end_s(i))
        m = TS.time_s >= W.time_start_s(i) & TS.time_s <= W.time_end_s(i);
        C.Step05_vs_FE_BTT_rmse_mm(i) = rmse_local( ...
            TS.u_step05_mm(m), TS.u_fe_btt_corrected_mm(m));
        C.Step05_vs_FE_BTT_corr(i) = corr_local( ...
            TS.u_step05_mm(m), TS.u_fe_btt_corrected_mm(m));
    end
end
C.Step05_A_inside_A_band = C.Step05_A_mm >= C.reference_A_low_mm & ...
    C.Step05_A_mm <= C.reference_A_high_mm;
C.Step05_A_inside_envelope_band = C.Step05_A_mm >= C.reference_envelope_low_mm & ...
    C.Step05_A_mm <= C.reference_envelope_high_mm;
C.Step05_to_reference_A_ratio = C.Step05_A_mm ./ C.reference_A_peak_equiv_mm;
C.Step05_to_FE_BTT_A_ratio = C.Step05_A_mm ./ C.FE_BTT_corrected_A_mm;
end


function Summary = build_step06k_transfer_summary_local(caseName, dataset, Rb, C, ...
    Fit, P, files)
Summary = table();
Summary.case_name = string(caseName);
Summary.dataset = string(dataset);
Summary.robustness_definition = ...
    "FE_K_band_constrained_strain_reference_with_BTT_phase_gain_correction";
Summary.input_step06i_mat = string(files.step06i_mat);
Summary.input_step03_bundle_mat = string(files.step03_bundle_mat);
Summary.region_count = height(Rb);
Summary.K_nominal_median_um_per_microstrain = median( ...
    Rb.K_nominal_um_per_microstrain, 'omitnan');
Summary.K_low_min_um_per_microstrain = min(Rb.K_low_um_per_microstrain, [], 'omitnan');
Summary.K_high_max_um_per_microstrain = max(Rb.K_high_um_per_microstrain, [], 'omitnan');
okFit = ~isempty(Fit) & strcmp(string(Fit.status), "ok");
Summary.btt_correction_point_count = height(P);
if any(okFit)
    Summary.FE_BTT_corrected_K_median_um_per_microstrain = median( ...
        Fit.K_um_per_microstrain(okFit), 'omitnan');
    Summary.FE_BTT_corrected_A_median_mm = median( ...
        Fit.corrected_A_peak_equiv_mm(okFit), 'omitnan');
    Summary.FE_BTT_fit_rmse_median_mm = median( ...
        Fit.offset_removed_rmse_mm(okFit), 'omitnan');
    Summary.FE_BTT_constraint_status = strjoin( ...
        unique(string(Fit.constraint_status(okFit))), ';');
else
    Summary.FE_BTT_corrected_K_median_um_per_microstrain = NaN;
    Summary.FE_BTT_corrected_A_median_mm = NaN;
    Summary.FE_BTT_fit_rmse_median_mm = NaN;
    Summary.FE_BTT_constraint_status = "";
end
Summary.reference_A_low_median_mm = median(Rb.A_peak_equiv_low_mm, 'omitnan');
Summary.reference_A_nominal_median_mm = median(Rb.u_nominal_A_peak_equiv_mm, 'omitnan');
Summary.reference_A_high_median_mm = median(Rb.A_peak_equiv_high_mm, 'omitnan');
if isempty(C)
    Summary.step05_window_count = 0;
    Summary.step05_A_inside_A_band_fraction = NaN;
    Summary.step05_A_inside_envelope_band_fraction = NaN;
    Summary.step05_to_reference_A_ratio_median = NaN;
    Summary.step05_to_FE_BTT_A_ratio_median = NaN;
    Summary.step05_vs_FE_BTT_rmse_median_mm = NaN;
else
    okA = isfinite(C.Step05_A_mm) & isfinite(C.reference_A_low_mm) & ...
        isfinite(C.reference_A_high_mm);
    okE = isfinite(C.Step05_A_mm) & isfinite(C.reference_envelope_low_mm) & ...
        isfinite(C.reference_envelope_high_mm);
    Summary.step05_window_count = height(C);
    Summary.step05_A_inside_A_band_fraction = mean(C.Step05_A_inside_A_band(okA));
    Summary.step05_A_inside_envelope_band_fraction = mean( ...
        C.Step05_A_inside_envelope_band(okE));
    Summary.step05_to_reference_A_ratio_median = median( ...
        C.Step05_to_reference_A_ratio, 'omitnan');
    Summary.step05_to_FE_BTT_A_ratio_median = median( ...
        C.Step05_to_FE_BTT_A_ratio, 'omitnan');
    Summary.step05_vs_FE_BTT_rmse_median_mm = median( ...
        C.Step05_vs_FE_BTT_rmse_mm, 'omitnan');
end
end


function figFile = plot_step06k_transfer_robustness_local(TS, Rb, C, Fit, P, ...
    Summary, figDir, dataset, cfg)
fig = figure('Visible', cfg.step06k_transfer_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 28, 18], ...
    'Name', 'Step06K transfer-ratio robustness');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
gray = [0.55 0.55 0.55];
purple = [0.494 0.184 0.556];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
fill_band_local(ax1, TS.time_s, -TS.u_abs_high_mm, TS.u_abs_high_mm, ...
    [0.90 0.92 0.96], 'transfer high envelope');
fill_band_local(ax1, TS.time_s, -TS.u_abs_low_mm, TS.u_abs_low_mm, ...
    [0.82 0.88 0.98], 'transfer low envelope');
plot_downsampled_local(ax1, TS.time_s, TS.u_transfer_nominal_mm, blue, ...
    'nominal signed reference', 0.8);
if ismember('u_fe_btt_corrected_mm', TS.Properties.VariableNames) && ...
        any(isfinite(TS.u_fe_btt_corrected_mm))
    plot_downsampled_local(ax1, TS.time_s, TS.u_fe_btt_corrected_mm, purple, ...
        'FE-band+BTT corrected', 0.95);
end
if any(isfinite(TS.u_step05_mm))
    plot_downsampled_local(ax1, TS.time_s, TS.u_step05_mm, gray, ...
        'Step05 waveform', 0.45);
end
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Displacement (mm)');
title(ax1, 'Time-domain transfer-ratio band');
legend(ax1, 'Location', 'best');
grid(ax1, 'on');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
errorbar(ax2, Rb.region_id, Rb.u_nominal_A_peak_equiv_mm, ...
    Rb.u_nominal_A_peak_equiv_mm - Rb.A_peak_equiv_low_mm, ...
    Rb.A_peak_equiv_high_mm - Rb.u_nominal_A_peak_equiv_mm, ...
    'o-', 'Color', purple, 'MarkerFaceColor', purple, 'CapSize', 7);
hold(ax2, 'on');
okFit = ~isempty(Fit) & strcmp(string(Fit.status), "ok");
if any(okFit)
    plot(ax2, Fit.region_id(okFit), Fit.corrected_A_peak_equiv_mm(okFit), ...
        's-', 'Color', blue, 'MarkerFaceColor', blue, ...
        'DisplayName', 'FE-band+BTT corrected A');
end
xlabel(ax2, 'Region ID');
ylabel(ax2, 'A peak equiv. (mm)');
title(ax2, 'Region amplitude robustness');
grid(ax2, 'on');

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
if ~isempty(C)
    plot(ax3, C.window_center_time_s, C.Step05_A_mm, 'o-', ...
        'Color', gray, 'MarkerFaceColor', gray, 'DisplayName', 'Step05 A');
    plot(ax3, C.window_center_time_s, C.reference_A_low_mm, ':', ...
        'Color', orange, 'DisplayName', 'reference A band');
    plot(ax3, C.window_center_time_s, C.reference_A_high_mm, ':', ...
        'Color', orange, 'HandleVisibility', 'off');
    plot(ax3, C.window_center_time_s, C.reference_A_peak_equiv_mm, '-', ...
        'Color', blue, 'DisplayName', 'nominal reference A');
    if ismember('FE_BTT_corrected_A_mm', C.Properties.VariableNames)
        plot(ax3, C.window_center_time_s, C.FE_BTT_corrected_A_mm, 's-', ...
            'Color', purple, 'MarkerFaceColor', purple, ...
            'DisplayName', 'FE-band+BTT corrected A');
    end
    legend(ax3, 'Location', 'best');
else
    text(ax3, 0.5, 0.5, 'No Step05 comparison', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Amplitude (mm)');
title(ax3, 'Step05 amplitude inside transfer band');
grid(ax3, 'on');

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
bar(ax4, categorical({'A band', 'envelope band'}), ...
    [Summary.step05_A_inside_A_band_fraction(1), ...
    Summary.step05_A_inside_envelope_band_fraction(1)]);
ylim(ax4, [0 1]);
ylabel(ax4, 'Fraction');
title(ax4, 'Band-coverage summary');
grid(ax4, 'on');

if ~isempty(P) && any(isfinite(P.displacement_mm))
    yyaxis(ax4, 'right');
    scatter(ax4, categorical(repmat("BTT pts", height(P), 1)), ...
        P.displacement_mm - median(P.displacement_mm, 'omitnan'), ...
        8, 'filled', 'MarkerFaceAlpha', 0.20, ...
        'DisplayName', 'BTT offset-removed points');
    ylabel(ax4, 'BTT centered displacement (mm)');
end

sg = sgtitle(fig, sprintf(['Step06K %s: transfer-ratio robustness ', ...
    '+ FE-band BTT correction\nK %.3f [%.3f, %.3f] um/microstrain; corrected K %.3f'], ...
    dataset, Summary.K_nominal_median_um_per_microstrain(1), ...
    Summary.K_low_min_um_per_microstrain(1), Summary.K_high_max_um_per_microstrain(1), ...
    Summary.FE_BTT_corrected_K_median_um_per_microstrain(1)));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.5);
style_figure_local(fig);

ensure_dir_local(figDir);
figFile = fullfile(figDir, sprintf( ...
    'Step06K_TransferRatioRobustness_%s.png', dataset));
save_figure_local(fig, figFile);
end


function fill_band_local(ax, x, ylo, yhi, color, name)
step = max(1, ceil(numel(x) / 30000));
x = x(1:step:end);
ylo = ylo(1:step:end);
yhi = yhi(1:step:end);
valid = isfinite(x) & isfinite(ylo) & isfinite(yhi);
x = x(valid); ylo = ylo(valid); yhi = yhi(valid);
if isempty(x)
    return;
end
patch(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.60, 'DisplayName', name);
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


function value = get_optional_field_local(S, name, defaultValue)
value = defaultValue;
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
end
end


function r = corr_local(x, y)
valid = isfinite(x) & isfinite(y);
if nnz(valid) < 3
    r = NaN;
else
    r = corr(x(valid), y(valid));
end
end


function v = rmse_local(y, yFit)
valid = isfinite(y) & isfinite(yFit);
if ~any(valid)
    v = NaN;
else
    v = sqrt(mean((y(valid) - yFit(valid)).^2, 'omitnan'));
end
end


function style_figure_local(fig)
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
    set(leg(i), 'FontName', 'Times New Roman', 'FontSize', 7.2, ...
        'Box', 'off');
end
end


function save_figure_local(fig, pngFile)
[folder, base, ~] = fileparts(pngFile);
ensure_dir_local(folder);
pngFile = fullfile(folder, [base '.png']);
pdfFile = fullfile(folder, [base '.pdf']);
figFile = fullfile(folder, [base '.fig']);
try
    exportgraphics(fig, pngFile, 'Resolution', 240);
    exportgraphics(fig, pdfFile, 'ContentType', 'image', 'Resolution', 240);
catch
    print(fig, pngFile, '-dpng', '-r240');
    print(fig, pdfFile, '-dpdf', '-painters');
end
try
    savefig(fig, figFile);
catch ME
    warning('savefig failed: %s', ME.message);
end
end


function ensure_dir_local(p)
if ~isfolder(p)
    mkdir(p);
end
end

