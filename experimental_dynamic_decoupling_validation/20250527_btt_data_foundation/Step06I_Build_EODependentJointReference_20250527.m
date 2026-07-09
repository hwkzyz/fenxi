%% Step06I_Build_EODependentJointReference_20250527.m
% Step06I: strain + FE transfer ratio only. Build an independent
% displacement reference first; Step05 is used only for comparison.

clc; clear; close all;

dataset = '20250527';
cfg = BTTDataConfig_20250527();
cfg.step06i_target_blade = 1;
cfg.step06i_step05_sensor_ids = cfg.sensor_ids(:).';
cfg.step06i_reference_region_ids = 1;
cfg.step06i_target_time_range_s = [1.0 5.0];
cfg.step06i_zoom_time_range_s = [1.50 1.56];
cfg.step06i_simplified_rectangular_apdl_root = ...
    'D:\博士-国科\试验台数据\大台子\叶片模型\580Hz传递比分析\simplified_rectangular_APDL_saved';

%RUN_STEP06I_STRAINTRANSFERDISPLACEMENTREFERENCE
% Build an independent displacement reference from Step06A strain evidence
% and an FE/APDL strain-to-tip transfer ratio. Step05 is used only after the
% reference is built, for comparison.

dataset = char(dataset);
cfg = apply_step06i_transfer_defaults_local(cfg, dataset);

targetCases = cfg.dynamic_cases;
if ischar(targetCases) || isstring(targetCases)
    targetCases = cellstr(targetCases);
end

for iCase = 1:numel(targetCases)
    caseName = char(targetCases{iCase});
    fprintf('\n=== Step06I strain-transfer displacement reference: %s / %s ===\n', ...
        dataset, caseName);

    files = resolve_step06i_transfer_files_local(cfg, dataset, caseName);
    assert_step06i_transfer_files_local(files);

    S = load(files.step06a_mat, 'Evidence');
    Evidence = S.Evidence;
    RegionTable = normalize_region_table_local(get_evidence_table_local(Evidence));
    [tStrain, yStrain, fsStrain] = get_strain_arrays_local(Evidence);
    oprTimes = load_opr_times_local(files.jiluopr_mat);

    SelectedRegions = select_step06i_transfer_regions_local(RegionTable, cfg);
    TransferPrior = build_transfer_prior_table_local(SelectedRegions, cfg, dataset);
    ReferenceTimeSeries = build_transfer_reference_time_series_local( ...
        SelectedRegions, TransferPrior, tStrain, yStrain, fsStrain, oprTimes, cfg);
    RegionSummary = build_region_summary_local(ReferenceTimeSeries, SelectedRegions, ...
        TransferPrior, cfg);

    Step05 = load_step05_for_comparison_local(files.step05_trend_csv, cfg);
    ReferenceTimeSeries = attach_nearest_step05_waveform_local( ...
        ReferenceTimeSeries, Step05, oprTimes, cfg);
    WindowComparison = build_step05_window_comparison_local( ...
        ReferenceTimeSeries, RegionSummary, Step05, cfg);
    Summary = build_step06i_transfer_summary_local(caseName, dataset, ...
        TransferPrior, RegionSummary, WindowComparison, Step05, files, cfg);

    outDir = fullfile(cfg.step06i_transfer_output_dir, caseName);
    figDir = fullfile(cfg.step06i_transfer_figure_dir, caseName);
    ensure_dir_local(outDir);
    ensure_dir_local(figDir);

    summaryCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_Summary_%s.csv', dataset));
    regionCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_SelectedRegions_%s.csv', dataset));
    transferCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_TransferPrior_%s.csv', dataset));
    regionSummaryCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_RegionSummary_%s.csv', dataset));
    tsCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_TimeSeries_%s.csv', dataset));
    windowCsv = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_Step05WindowComparison_%s.csv', dataset));
    matFile = fullfile(outDir, sprintf( ...
        'Step06I_StrainTransferDisplacementReference_%s.mat', dataset));

    writetable(Summary, summaryCsv);
    writetable(SelectedRegions, regionCsv);
    writetable(TransferPrior, transferCsv);
    writetable(RegionSummary, regionSummaryCsv);
    writetable(ReferenceTimeSeries, tsCsv);
    writetable(WindowComparison, windowCsv);
    save(matFile, 'Summary', 'SelectedRegions', 'TransferPrior', ...
        'RegionSummary', 'ReferenceTimeSeries', 'WindowComparison', ...
        'Step05', 'files', 'cfg', '-v7.3');

    figFile = plot_step06i_transfer_reference_local(ReferenceTimeSeries, ...
        SelectedRegions, TransferPrior, RegionSummary, WindowComparison, ...
        Summary, figDir, dataset, cfg);

    fprintf('\nStep06I transfer-reference summary:\n');
    disp(Summary);
    fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
        summaryCsv, regionCsv, transferCsv, regionSummaryCsv, tsCsv, ...
        windowCsv, matFile, figFile);
end

function cfg = apply_step06i_transfer_defaults_local(cfg, dataset)
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
if ~isfield(cfg, 'step06i_transfer_figure_dir') || ...
        isempty(cfg.step06i_transfer_figure_dir)
    cfg.step06i_transfer_figure_dir = fullfile(cfg.figure_root, ...
        'step06i_strain_transfer_displacement_reference');
end
if ~isfield(cfg, 'step06i_transfer_fig_visible') || ...
        isempty(cfg.step06i_transfer_fig_visible)
    if usejava('desktop')
        cfg.step06i_transfer_fig_visible = 'on';
    else
        cfg.step06i_transfer_fig_visible = 'off';
    end
end
if ~isfield(cfg, 'step06i_target_blade') || isempty(cfg.step06i_target_blade)
    if isfield(cfg, 'step05_target_blades') && ~isempty(cfg.step05_target_blades)
        cfg.step06i_target_blade = cfg.step05_target_blades(1);
    else
        cfg.step06i_target_blade = 1;
    end
end
if ~isfield(cfg, 'step06i_step05_sensor_ids') || ...
        isempty(cfg.step06i_step05_sensor_ids)
    if isfield(cfg, 'step05_analysis_sensors') && ~isempty(cfg.step05_analysis_sensors)
        cfg.step06i_step05_sensor_ids = cfg.step05_analysis_sensors(:).';
    elseif isfield(cfg, 'sensor_ids') && ~isempty(cfg.sensor_ids)
        cfg.step06i_step05_sensor_ids = cfg.sensor_ids(:).';
    else
        cfg.step06i_step05_sensor_ids = [];
    end
end
if ~isfield(cfg, 'step06i_target_time_range_s') || ...
        isempty(cfg.step06i_target_time_range_s)
    if strcmp(dataset, '20241106')
        cfg.step06i_target_time_range_s = [];
    else
        cfg.step06i_target_time_range_s = [1.0 5.0];
    end
end
if ~isfield(cfg, 'step06i_zoom_time_range_s') || ...
        isempty(cfg.step06i_zoom_time_range_s)
    if strcmp(dataset, '20241106')
        cfg.step06i_zoom_time_range_s = [];
    else
        cfg.step06i_zoom_time_range_s = [1.50 1.56];
    end
end
if ~isfield(cfg, 'step06i_transfer_demod_margin_s') || ...
        isempty(cfg.step06i_transfer_demod_margin_s)
    cfg.step06i_transfer_demod_margin_s = 0.80;
end
if ~isfield(cfg, 'step06i_transfer_envelope_lowpass_hz') || ...
        isempty(cfg.step06i_transfer_envelope_lowpass_hz)
    cfg.step06i_transfer_envelope_lowpass_hz = 8.0;
end
if ~isfield(cfg, 'step06i_stable_trim_fraction') || ...
        isempty(cfg.step06i_stable_trim_fraction)
    cfg.step06i_stable_trim_fraction = 0.15;
end
if ~isfield(cfg, 'step06i_stable_trim_min_s') || ...
        isempty(cfg.step06i_stable_trim_min_s)
    cfg.step06i_stable_trim_min_s = 0.5;
end
if ~isfield(cfg, 'step06i_stable_trim_max_s') || ...
        isempty(cfg.step06i_stable_trim_max_s)
    cfg.step06i_stable_trim_max_s = 2.0;
end
if ~isfield(cfg, 'step06i_transfer_frequency_switch_hz') || ...
        isempty(cfg.step06i_transfer_frequency_switch_hz)
    cfg.step06i_transfer_frequency_switch_hz = 605.0;
end
if ~isfield(cfg, 'step06i_transfer_surface_max_as_uncertainty') || ...
        isempty(cfg.step06i_transfer_surface_max_as_uncertainty)
    cfg.step06i_transfer_surface_max_as_uncertainty = true;
end
if ~isfield(cfg, 'step06i_transfer_gauge_neighborhood_mm') || ...
        isempty(cfg.step06i_transfer_gauge_neighborhood_mm)
    cfg.step06i_transfer_gauge_neighborhood_mm = 2.5;
end
if ~isfield(cfg, 'step06i_simplified_rectangular_apdl_root') || ...
        isempty(cfg.step06i_simplified_rectangular_apdl_root)
    cfg.step06i_simplified_rectangular_apdl_root = ...
        'D:\鍗氬＋-鍥界\璇曢獙鍙版暟鎹甛澶у彴瀛怽鍙剁墖妯″瀷\580Hz浼犻€掓瘮鍒嗘瀽\simplified_rectangular_APDL_saved';
end
if ~isfield(cfg, 'opr_events_per_revolution') || isempty(cfg.opr_events_per_revolution)
    if isfield(cfg, 'step04_opr_events_per_revolution') && ...
            ~isempty(cfg.step04_opr_events_per_revolution)
        cfg.opr_events_per_revolution = cfg.step04_opr_events_per_revolution;
    elseif isfield(cfg, 'opr_pulses_per_rev') && ~isempty(cfg.opr_pulses_per_rev)
        cfg.opr_events_per_revolution = cfg.opr_pulses_per_rev;
    else
        cfg.opr_events_per_revolution = 1;
    end
end
end


function files = resolve_step06i_transfer_files_local(cfg, dataset, caseName)
files = struct();
files.step06a_mat = fullfile(cfg.output_root, ...
    'step06a_strain_rpm_resonance_evidence', caseName, ...
    sprintf('Step06A_StrainRPM_ResonanceEvidence_%s.mat', dataset));
files.jiluopr_mat = fullfile(cfg.step02_output_dir, caseName, 'jiluOPR.mat');
files.step05_trend_csv = resolve_step05_trend_csv_local(cfg, dataset, caseName);
end


function file = resolve_step05_trend_csv_local(cfg, dataset, caseName)
sensorGroup = ['S', sprintf('%d', cfg.step06i_step05_sensor_ids)];
bladeId = cfg.step06i_target_blade;
step05Root = fullfile(cfg.output_root, 'step05_single_sync_direct_template', caseName);
if isfield(cfg, 'step05_output_dir') && ~isempty(cfg.step05_output_dir)
    step05Root = fullfile(cfg.step05_output_dir, caseName);
elseif isfield(cfg, 'step05_direct_template_output_dir') && ...
        ~isempty(cfg.step05_direct_template_output_dir)
    step05Root = fullfile(cfg.step05_direct_template_output_dir, caseName);
end

exactNames = { ...
    sprintf('Trend_Step05_FixedEtaGradDispVP_B%d_%s_%s.csv', ...
        bladeId, sensorGroup, dataset), ...
    sprintf('Trend_Step05_FoundationMainPulseAdaptiveFixedJointEta_B%d_%s_%s.csv', ...
        bladeId, sensorGroup, dataset), ...
    sprintf('Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
        bladeId, sensorGroup, dataset)};
for i = 1:numel(exactNames)
    c = fullfile(step05Root, exactNames{i});
    if isfile(c) && is_step05_trend_candidate_usable_local(c, step05Root, ...
            caseName, bladeId, sensorGroup)
        file = c;
        return;
    end
end

patterns = { ...
    sprintf('Trend_Step05_*Fixed*Eta*_B%d_%s_%s.csv', bladeId, sensorGroup, dataset), ...
    sprintf('Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', bladeId, sensorGroup, dataset)};
for i = 1:numel(patterns)
    dd = dir(fullfile(step05Root, patterns{i}));
    dd = filter_noeta_files_local(dd);
    if ~isempty(dd)
        [~, ix] = sort([dd.datenum], 'descend');
        dd = dd(ix);
        for j = 1:numel(dd)
            c = fullfile(dd(j).folder, dd(j).name);
            if is_step05_trend_candidate_usable_local(c, step05Root, ...
                    caseName, bladeId, sensorGroup)
                file = c;
                return;
            end
        end
    end
end

file = '';
end


function ok = is_step05_trend_candidate_usable_local(file, step05Root, ...
    caseName, bladeId, sensorGroup)
ok = true;
try
    T = readtable(file, 'TextType', 'string');
catch
    return;
end
if height(T) >= 3
    return;
end

summaryFile = fullfile(step05Root, '..', 'Step05_TargetBlade_Summary.csv');
datasetSummary = dir(fullfile(step05Root, '..', 'Step05_TargetBlade_Summary_*.csv'));
if isempty(datasetSummary)
    return;
end
summaryFile = fullfile(datasetSummary(1).folder, datasetSummary(1).name);
try
    S = readtable(summaryFile, 'TextType', 'string');
catch
    return;
end
if ~all(ismember(["case_name", "blade_id", "sensor_tag"], string(S.Properties.VariableNames)))
    return;
end
row = strcmp(string(S.case_name), string(caseName)) & ...
    S.blade_id == bladeId & strcmp(string(S.sensor_tag), string(sensorGroup));
if ~any(row)
    return;
end
expected = NaN;
if ismember("valid_window_count", string(S.Properties.VariableNames))
    expected = S.valid_window_count(find(row, 1));
elseif ismember("total_window_count", string(S.Properties.VariableNames))
    expected = S.total_window_count(find(row, 1));
end
if isfinite(expected) && expected >= 3 && height(T) < expected
    warning(['Ignoring stale Step05 trend with only %d rows; summary expects ' ...
        '%d windows:\n  %s'], height(T), expected, file);
    ok = false;
end
end


function dd = filter_noeta_files_local(dd)
keep = true(size(dd));
for i = 1:numel(dd)
    nm = lower(dd(i).name);
    if contains(nm, 'noeta') || contains(nm, 'eta0')
        keep(i) = false;
    end
end
dd = dd(keep);
end


function assert_step06i_transfer_files_local(files)
if ~isfile(files.step06a_mat)
    error('Missing Step06A evidence MAT:\n  %s', files.step06a_mat);
end
if ~isfile(files.jiluopr_mat)
    error('Missing Step02 OPR MAT:\n  %s', files.jiluopr_mat);
end
if isempty(files.step05_trend_csv) || ~isfile(files.step05_trend_csv)
    warning('Step05 trend CSV not found. Reference will be built without Step05 comparison.');
end
end


function T = get_evidence_table_local(Evidence)
if isfield(Evidence, 'RegionTable') && istable(Evidence.RegionTable)
    T = Evidence.RegionTable;
elseif isfield(Evidence, 'regionTable') && istable(Evidence.regionTable)
    T = Evidence.regionTable;
else
    error('Step06A Evidence does not contain RegionTable/regionTable.');
end
end


function R = normalize_region_table_local(T)
R = table();
R.region_id = get_table_numeric_local(T, {'region_id', 'regionId'});
R.time_start_s = get_table_numeric_local(T, {'time_start_s', 'bttStartSec', 'regionStart'});
R.time_end_s = get_table_numeric_local(T, {'time_end_s', 'bttEndSec', 'regionEnd'});
R.peak_time_s = get_table_numeric_local(T, {'peak_time_s', 'peakBttTimeSec', 'regionPeakTime'});
R.dominant_order = get_table_numeric_local(T, {'dominant_order', 'dominantOrder', 'regionPeakOrder'});
R.peak_freq_hz = get_table_numeric_local(T, {'peak_freq_hz', 'dominantFreqHz', 'regionPeakFreqHz'});
R.peak_amp_microstrain = get_table_numeric_local(T, ...
    {'peak_amp_microstrain', 'peakFftAmpMicrostrain', 'fftPeakAmpMicrostrain'});
R.duration_s = R.time_end_s - R.time_start_s;
ok = isfinite(R.time_start_s) & isfinite(R.time_end_s) & ...
    R.time_end_s > R.time_start_s & isfinite(R.dominant_order) & ...
    isfinite(R.peak_freq_hz);
R = R(ok, :);
R = sortrows(R, 'time_start_s');
if isempty(R)
    error('No usable Step06A resonance regions.');
end
end


function Selected = select_step06i_transfer_regions_local(R, cfg)
if isfield(cfg, 'step06i_reference_region_ids') && ...
        ~isempty(cfg.step06i_reference_region_ids)
    keep = ismember(R.region_id, cfg.step06i_reference_region_ids(:));
    Selected = R(keep, :);
    if isempty(Selected)
        error('No Step06A regions match cfg.step06i_reference_region_ids.');
    end
elseif isempty(cfg.step06i_target_time_range_s)
    Selected = R;
else
    tr = cfg.step06i_target_time_range_s(:).';
    keep = R.time_end_s >= tr(1) & R.time_start_s <= tr(2);
    Selected = R(keep, :);
    if isempty(Selected)
        Selected = R;
        warning('No Step06A regions overlap target_time_range; using all regions.');
    end
end
Selected = sortrows(Selected, 'time_start_s');
end


function TP = build_transfer_prior_table_local(R, cfg, dataset)
rows = repmat(empty_transfer_row_local(), height(R), 1);
priorPackage = load_simplified_rectangular_apdl_package_local(cfg);
for i = 1:height(R)
    freq = R.peak_freq_hz(i);
    p = select_simplified_rectangular_apdl_prior_local(priorPackage, freq, cfg);
    if isempty(p)
        if freq < cfg.step06i_transfer_frequency_switch_hz
            p = same_apdl_580_prior_local();
        else
            p = same_apdl_630_prior_local();
        end
    end
    if freq < cfg.step06i_transfer_frequency_switch_hz
        p.case_name = strrep(p.case_name, '630', '580');
    else
        p.case_name = strrep(p.case_name, '580', '630');
    end
    rows(i).region_id = R.region_id(i);
    rows(i).dataset = string(dataset);
    rows(i).peak_freq_hz = freq;
    rows(i).prior_case = string(p.case_name);
    rows(i).actual_mode_frequency_hz = p.actual_mode_frequency_hz;
    rows(i).transfer_nominal_strain_per_m = p.g1_nominal_strain_per_m;
    rows(i).transfer_front_strain_per_m = p.g1_front_strain_per_m;
    rows(i).transfer_back_strain_per_m = p.g1_back_strain_per_m;
    rows(i).surface_max_front_strain_per_m = p.surface_max_front_strain_per_m;
    rows(i).surface_max_back_strain_per_m = p.surface_max_back_strain_per_m;
    rows(i).K_nominal_mm_per_microstrain = 1e-3 / p.g1_nominal_strain_per_m;
    rows(i).K_nominal_um_per_microstrain = 1 / p.g1_nominal_strain_per_m;
    if isfield(p, 'local_min_strain_per_m') && ...
            isfinite(p.local_min_strain_per_m) && ...
            isfinite(p.local_max_strain_per_m) && p.local_min_strain_per_m > 0
        tMin = p.local_min_strain_per_m;
        tMax = p.local_max_strain_per_m;
    elseif cfg.step06i_transfer_surface_max_as_uncertainty
        tMin = p.g1_nominal_strain_per_m;
        tMax = max([p.surface_max_front_strain_per_m, ...
            p.surface_max_back_strain_per_m, p.g1_nominal_strain_per_m]);
    else
        tMin = min([p.g1_front_strain_per_m, p.g1_back_strain_per_m]);
        tMax = max([p.g1_front_strain_per_m, p.g1_back_strain_per_m]);
    end
    rows(i).K_band_low_um_per_microstrain = 1 / tMax;
    rows(i).K_band_high_um_per_microstrain = 1 / tMin;
    rows(i).source_note = string(p.source_note);
end
TP = struct2table(rows);
end


function pkg = load_simplified_rectangular_apdl_package_local(cfg)
pkg = struct('available', false, 'summary', table(), 'faceSummary', table(), ...
    'surface', table(), 'g1', table(), 'root', string(''));
root = string(get_optional_field_local(cfg, ...
    'step06i_simplified_rectangular_apdl_root', ''));
if strlength(root) == 0
    return;
end
summaryCsv = fullfile(root, 'simplified_rectangular_apdl_summary.csv');
faceCsv = fullfile(root, 'surface_transfer_dataset_same_apdl', ...
    'surface_transfer_same_apdl_face_summary.csv');
surfaceCsv = fullfile(root, 'surface_transfer_dataset_same_apdl', ...
    'surface_transfer_same_apdl_all.csv');
g1Csv = fullfile(root, 'surface_transfer_dataset_same_apdl', ...
    'surface_transfer_same_apdl_at_G1.csv');
if ~isfile(summaryCsv) || ~isfile(faceCsv)
    return;
end
try
    pkg.summary = readtable(summaryCsv, 'TextType', 'string');
    pkg.faceSummary = readtable(faceCsv, 'TextType', 'string');
    if isfile(surfaceCsv)
        pkg.surface = readtable(surfaceCsv, 'TextType', 'string');
    end
    if isfile(g1Csv)
        pkg.g1 = readtable(g1Csv, 'TextType', 'string');
    end
    pkg.root = root;
    pkg.available = true;
catch
    pkg.available = false;
end
end


function p = select_simplified_rectangular_apdl_prior_local(pkg, freq, cfg)
p = [];
if ~isstruct(pkg) || ~isfield(pkg, 'available') || ~pkg.available
    return;
end
if freq < cfg.step06i_transfer_frequency_switch_hz
    label = "580Hz";
    faceLabel = "580 Hz";
else
    label = "630Hz";
    faceLabel = "630 Hz";
end
S = pkg.summary;
F = pkg.faceSummary;
ms = strcmp(string(S.frequency_group), label);
mf = strcmp(string(F.case_label), faceLabel);
if ~any(ms) || ~any(mf)
    return;
end
s = S(find(ms, 1), :);
ff = F(mf, :);
frontRow = ff(strcmpi(string(ff.face), "front"), :);
backRow = ff(strcmpi(string(ff.face), "back"), :);
if isempty(frontRow) || isempty(backRow)
    return;
end
p = struct();
p.case_name = sprintf('simplified_rectangular_APDL_%s_G1', char(label));
p.actual_mode_frequency_hz = s.actual_frequency_hz(1);
p.g1_front_strain_per_m = s.front_G1_transfer_1_per_m(1);
p.g1_back_strain_per_m = s.back_G1_transfer_1_per_m(1);
p.g1_nominal_strain_per_m = mean([p.g1_front_strain_per_m, ...
    p.g1_back_strain_per_m], 'omitnan');
p.surface_max_front_strain_per_m = frontRow.max_transfer_ratio_1_per_m(1);
p.surface_max_back_strain_per_m = backRow.max_transfer_ratio_1_per_m(1);
p.local_min_strain_per_m = NaN;
p.local_max_strain_per_m = NaN;
p.local_neighborhood_mm = cfg.step06i_transfer_gauge_neighborhood_mm;
if ~isempty(pkg.surface) && ~isempty(pkg.g1)
    [p.local_min_strain_per_m, p.local_max_strain_per_m] = ...
        local_transfer_band_around_gauge_local(pkg.surface, pkg.g1, ...
        faceLabel, cfg.step06i_transfer_gauge_neighborhood_mm);
end
p.source_note = sprintf(['Simplified rectangular APDL saved package: ', ...
    'transfer_ratio_1_per_m = abs(EPELX / max_abs_tip_UY), ', ...
    'G1 nominal; K band from +/-%.2f mm gauge neighborhood when available; root=%s'], ...
    cfg.step06i_transfer_gauge_neighborhood_mm, char(pkg.root));
end


function [tMin, tMax] = local_transfer_band_around_gauge_local(S, G1, ...
    caseLabel, neighborhoodMm)
tMin = NaN;
tMax = NaN;
g = G1(strcmp(string(G1.case_name), string(caseLabel)) & ...
    strcmpi(string(G1.face), "front"), :);
if isempty(g)
    g = G1(strcmp(string(G1.case_name), string(caseLabel)), :);
end
if isempty(g) || ~all(ismember({'length_mm', 'width_mm'}, ...
        G1.Properties.VariableNames))
    return;
end
x0 = g.length_mm(1);
z0 = g.width_mm(1);
m = strcmp(string(S.case_label), string(caseLabel)) & ...
    abs(S.x_length_mm - x0) <= neighborhoodMm & ...
    abs(S.z_width_mm - z0) <= neighborhoodMm & ...
    isfinite(S.transfer_ratio_1_per_m) & S.transfer_ratio_1_per_m > 0;
if ~any(m)
    return;
end
vals = S.transfer_ratio_1_per_m(m);
tMin = min(vals, [], 'omitnan');
tMax = max(vals, [], 'omitnan');
end


function row = empty_transfer_row_local()
row = struct('region_id', NaN, 'dataset', "", 'peak_freq_hz', NaN, ...
    'prior_case', "", 'actual_mode_frequency_hz', NaN, ...
    'transfer_nominal_strain_per_m', NaN, ...
    'transfer_front_strain_per_m', NaN, ...
    'transfer_back_strain_per_m', NaN, ...
    'surface_max_front_strain_per_m', NaN, ...
    'surface_max_back_strain_per_m', NaN, ...
    'K_nominal_mm_per_microstrain', NaN, ...
    'K_nominal_um_per_microstrain', NaN, ...
    'K_band_low_um_per_microstrain', NaN, ...
    'K_band_high_um_per_microstrain', NaN, ...
    'source_note', "");
end


function p = same_apdl_580_prior_local()
p = struct();
p.case_name = 'same_apdl_580_Hz_G1';
p.actual_mode_frequency_hz = 579.993;
p.g1_front_strain_per_m = 0.9466;
p.g1_back_strain_per_m = 0.9465;
p.g1_nominal_strain_per_m = mean([p.g1_front_strain_per_m, p.g1_back_strain_per_m]);
p.surface_max_front_strain_per_m = 1.2349;
p.surface_max_back_strain_per_m = 1.2257;
p.source_note = ['Same APDL parameterized SOLID186 model: ', ...
    'transfer = abs(EPELX / max_abs_tip_UY), G1 nominal.'];
end


function p = same_apdl_630_prior_local()
p = struct();
p.case_name = 'same_apdl_630_Hz_G1';
p.actual_mode_frequency_hz = 629.893;
p.g1_front_strain_per_m = 1.0289;
p.g1_back_strain_per_m = 1.0289;
p.g1_nominal_strain_per_m = mean([p.g1_front_strain_per_m, p.g1_back_strain_per_m]);
p.surface_max_front_strain_per_m = 1.3296;
p.surface_max_back_strain_per_m = 1.3290;
p.source_note = ['Same APDL parameterized SOLID186 model: ', ...
    'transfer = abs(EPELX / max_abs_tip_UY), G1 nominal.'];
end


function TS = build_transfer_reference_time_series_local(R, TP, tStrain, yStrain, ...
    fs, oprTimes, cfg)
TS = table();
for i = 1:height(R)
    t0 = R.time_start_s(i);
    t1 = R.time_end_s(i);
    eo = R.dominant_order(i);
    rid = R.region_id(i);
    margin = cfg.step06i_transfer_demod_margin_s;
    wide = tStrain >= t0 - margin & tStrain <= t1 + margin;
    tw = tStrain(wide);
    yw = yStrain(wide);
    core = tw >= t0 & tw <= t1;
    if nnz(core) < 16
        continue;
    end
    yw = yw - median(yw(core), 'omitnan');
    theta = map_time_to_rotor_phase_local(oprTimes, tw, cfg.opr_events_per_revolution);
    psi = eo .* theta;
    z = yw .* exp(-1i .* psi);
    q = fft_lowpass_complex_local(z, fs, cfg.step06i_transfer_envelope_lowpass_hz);
    h = 2 .* q .* exp(1i .* psi);

    tc = tw(core);
    hc = h(core);
    thetaC = theta(core);
    finite = isfinite(tc) & isfinite(real(hc)) & isfinite(imag(hc)) & ...
        isfinite(thetaC);
    tc = tc(finite);
    hc = hc(finite);
    thetaC = thetaC(finite);
    if isempty(tc)
        continue;
    end

    tp = TP(TP.region_id == rid, :);
    K = tp.K_nominal_mm_per_microstrain(1);
    kLow = tp.K_band_low_um_per_microstrain(1) / 1000;
    kHigh = tp.K_band_high_um_per_microstrain(1) / 1000;

    T = table();
    T.time_s = tc(:);
    T.region_id = repmat(rid, numel(tc), 1);
    T.dominant_order = repmat(eo, numel(tc), 1);
    T.peak_freq_hz = repmat(R.peak_freq_hz(i), numel(tc), 1);
    T.theta_rot_rad = thetaC(:);
    T.strain_basis_real_microstrain = real(hc(:));
    T.strain_basis_imag_microstrain = imag(hc(:));
    T.strain_basis_abs_microstrain = abs(hc(:));
    T.u_transfer_nominal_mm = K .* real(hc(:));
    T.u_transfer_negative_sign_mm = -T.u_transfer_nominal_mm;
    T.u_transfer_abs_envelope_mm = K .* abs(hc(:));
    T.u_transfer_band_low_abs_mm = kLow .* abs(hc(:));
    T.u_transfer_band_high_abs_mm = kHigh .* abs(hc(:));
    TS = [TS; T]; %#ok<AGROW>
end
if isempty(TS)
    error('No Step06I transfer reference samples were built.');
end
TS = sortrows(TS, 'time_s');
end


function coreMask = select_stable_core_mask_local(t, regionMask, cfg)
coreMask = false(size(t));
idx = find(regionMask & isfinite(t));
if isempty(idx)
    return;
end
t0 = min(t(idx), [], 'omitnan');
t1 = max(t(idx), [], 'omitnan');
duration = t1 - t0;
if ~isfinite(duration) || duration <= 0
    coreMask(idx) = true;
    return;
end
trim = cfg.step06i_stable_trim_fraction .* duration;
trim = min(max(trim, cfg.step06i_stable_trim_min_s), ...
    cfg.step06i_stable_trim_max_s);
if 2 .* trim >= 0.80 .* duration
    trim = 0.10 .* duration;
end
tc0 = t0 + trim;
tc1 = t1 - trim;
if ~isfinite(tc0) || ~isfinite(tc1) || tc0 >= tc1
    coreMask(idx) = true;
else
    coreMask = regionMask & t >= tc0 & t <= tc1 & isfinite(t);
end
end


function RS = build_region_summary_local(TS, R, TP, cfg)
rows = repmat(empty_region_summary_row_local(), height(R), 1);
for i = 1:height(R)
    rid = R.region_id(i);
    m = TS.region_id == rid;
    tp = TP(TP.region_id == rid, :);
    coreMask = m & select_stable_core_mask_local(TS.time_s, m, cfg);
    if nnz(coreMask) < max(16, ceil(0.20 * nnz(m)))
        coreMask = m;
    end
    strainAmp = median(TS.strain_basis_abs_microstrain(coreMask), 'omitnan');
    if ~isfinite(strainAmp) || strainAmp <= 0
        strainAmp = median(TS.strain_basis_abs_microstrain(m), 'omitnan');
    end
    KNom = tp.K_nominal_um_per_microstrain(1) / 1000;
    KLow = tp.K_band_low_um_per_microstrain(1) / 1000;
    KHigh = tp.K_band_high_um_per_microstrain(1) / 1000;
    rows(i).region_id = rid;
    rows(i).time_start_s = R.time_start_s(i);
    rows(i).time_end_s = R.time_end_s(i);
    rows(i).core_time_start_s = min(TS.time_s(coreMask), [], 'omitnan');
    rows(i).core_time_end_s = max(TS.time_s(coreMask), [], 'omitnan');
    rows(i).dominant_order = R.dominant_order(i);
    rows(i).peak_freq_hz = R.peak_freq_hz(i);
    rows(i).transfer_prior_case = string(tp.prior_case(1));
    rows(i).K_nominal_um_per_microstrain = tp.K_nominal_um_per_microstrain(1);
    rows(i).K_band_low_um_per_microstrain = tp.K_band_low_um_per_microstrain(1);
    rows(i).K_band_high_um_per_microstrain = tp.K_band_high_um_per_microstrain(1);
    rows(i).amplitude_definition = ...
        "stable_core_EO_frequency_amplitude_median_abs_demod";
    rows(i).stable_frequency_amp_microstrain = strainAmp;
    rows(i).strain_abs_median_microstrain = median( ...
        TS.strain_basis_abs_microstrain(coreMask), 'omitnan');
    rows(i).strain_abs_p95_microstrain = percentile_local( ...
        TS.strain_basis_abs_microstrain(coreMask), 95);
    rows(i).u_signed_waveform_std_A_mm = sqrt(2) * ...
        std(TS.u_transfer_nominal_mm(coreMask), 'omitnan');
    rows(i).u_frequency_amp_nominal_mm = KNom .* strainAmp;
    rows(i).u_frequency_amp_band_low_mm = KLow .* strainAmp;
    rows(i).u_frequency_amp_band_high_mm = KHigh .* strainAmp;
    rows(i).u_nominal_A_peak_equiv_mm = rows(i).u_frequency_amp_nominal_mm;
    rows(i).u_abs_envelope_median_mm = median( ...
        TS.u_transfer_abs_envelope_mm(coreMask), 'omitnan');
    rows(i).u_abs_envelope_p95_mm = percentile_local( ...
        TS.u_transfer_abs_envelope_mm(coreMask), 95);
    rows(i).u_band_low_abs_median_mm = median( ...
        TS.u_transfer_band_low_abs_mm(coreMask), 'omitnan');
    rows(i).u_band_high_abs_median_mm = median( ...
        TS.u_transfer_band_high_abs_mm(coreMask), 'omitnan');
end
RS = struct2table(rows);
end


function row = empty_region_summary_row_local()
row = struct('region_id', NaN, 'time_start_s', NaN, 'time_end_s', NaN, ...
    'core_time_start_s', NaN, 'core_time_end_s', NaN, ...
    'dominant_order', NaN, 'peak_freq_hz', NaN, ...
    'transfer_prior_case', "", 'K_nominal_um_per_microstrain', NaN, ...
    'K_band_low_um_per_microstrain', NaN, ...
    'K_band_high_um_per_microstrain', NaN, ...
    'amplitude_definition', "", ...
    'stable_frequency_amp_microstrain', NaN, ...
    'strain_abs_median_microstrain', NaN, ...
    'strain_abs_p95_microstrain', NaN, ...
    'u_signed_waveform_std_A_mm', NaN, ...
    'u_frequency_amp_nominal_mm', NaN, ...
    'u_frequency_amp_band_low_mm', NaN, ...
    'u_frequency_amp_band_high_mm', NaN, ...
    'u_nominal_A_peak_equiv_mm', NaN, ...
    'u_abs_envelope_median_mm', NaN, ...
    'u_abs_envelope_p95_mm', NaN, ...
    'u_band_low_abs_median_mm', NaN, ...
    'u_band_high_abs_median_mm', NaN);
end


function Step05 = load_step05_for_comparison_local(file, cfg)
Step05 = struct('available', false, 'Trend', table(), 'ok', false(0,1), ...
    'WindowTable', table(), 'source_file', string(file));
if isempty(file) || ~isfile(file)
    return;
end
T = readtable(file);
if ~ismember('status', T.Properties.VariableNames)
    statusOk = true(height(T), 1);
else
    statusOk = strcmpi(string(T.status), "ok");
end
required = {'window_center_time_s', 'A_id', 'EO_id', 'phi_id_wrapped'};
for i = 1:numel(required)
    if ~ismember(required{i}, T.Properties.VariableNames)
        warning('Step05 trend is missing column %s:\n  %s', required{i}, file);
        return;
    end
end
ok = statusOk & isfinite(T.window_center_time_s) & ...
    isfinite(T.A_id) & isfinite(T.EO_id) & isfinite(T.phi_id_wrapped);
if ~isempty(cfg.step06i_target_time_range_s)
    tr = cfg.step06i_target_time_range_s(:).';
    ok = ok & T.window_center_time_s >= tr(1) & T.window_center_time_s <= tr(2);
end
Step05.available = any(ok);
Step05.Trend = T;
Step05.ok = ok;
Step05.WindowTable = build_step05_window_table_local(T, ok);
end


function W = build_step05_window_table_local(T, ok)
W = T(ok, :);
if isempty(W)
    W = table();
    return;
end
if ~ismember('window_id', W.Properties.VariableNames)
    W.window_id = (1:height(W)).';
end
if ~ismember('fn_id', W.Properties.VariableNames)
    W.fn_id = nan(height(W), 1);
end
rotHz = get_table_numeric_local(W, {'rot_freq_mean_hz'});
validLapCount = get_table_numeric_local(W, {'valid_lap_count'});
lapStart = get_table_numeric_local(W, {'lap_start'});
lapEnd = get_table_numeric_local(W, {'lap_end'});
lapCount = validLapCount;
fromLapSpan = isfinite(lapStart) & isfinite(lapEnd) & lapEnd >= lapStart;
lapCount(~isfinite(lapCount) & fromLapSpan) = ...
    lapEnd(~isfinite(lapCount) & fromLapSpan) - ...
    lapStart(~isfinite(lapCount) & fromLapSpan) + 1;
duration = lapCount ./ rotHz;
centers = W.window_center_time_s;
fallbackStep = median(diff(centers), 'omitnan');
if ~isfinite(fallbackStep) || fallbackStep <= 0
    fallbackStep = 0.05;
end
bad = ~isfinite(duration) | duration <= 0;
duration(bad) = fallbackStep;
W.time_start_s = centers - 0.5 .* duration;
W.time_end_s = centers + 0.5 .* duration;
W.window_duration_s = duration;
W = sortrows(W, 'window_center_time_s');
end


function TS = attach_nearest_step05_waveform_local(TS, Step05, oprTimes, cfg)
TS.u_step05_mm = nan(height(TS), 1);
TS.step05_A_mm = nan(height(TS), 1);
TS.step05_EO = nan(height(TS), 1);
TS.step05_phi_rad = nan(height(TS), 1);
TS.step05_fn_hz = nan(height(TS), 1);
if ~Step05.available
    return;
end
T = Step05.WindowTable;
centers = T.window_center_time_s;
theta = map_time_to_rotor_phase_local(oprTimes, TS.time_s, cfg.opr_events_per_revolution);
for i = 1:height(TS)
    inWin = TS.time_s(i) >= T.time_start_s & TS.time_s(i) <= T.time_end_s;
    if any(inWin)
        jj = find(inWin, 1, 'first');
    else
        [~, jj] = min(abs(centers - TS.time_s(i)));
    end
    A = T.A_id(jj);
    EO = T.EO_id(jj);
    phi = T.phi_id_wrapped(jj);
    TS.u_step05_mm(i) = A .* sin(EO .* theta(i) + phi);
    TS.step05_A_mm(i) = A;
    TS.step05_EO(i) = EO;
    TS.step05_phi_rad(i) = phi;
    TS.step05_fn_hz(i) = T.fn_id(jj);
end
end


function W = build_step05_window_comparison_local(TS, RS, Step05, cfg)
if ~Step05.available
    W = table();
    return;
end
T = Step05.WindowTable;
rows = repmat(empty_window_row_local(), height(T), 1);
for k = 1:height(T)
    tr = [T.time_start_s(k), T.time_end_s(k)];
    m = TS.time_s >= tr(1) & TS.time_s <= tr(2);
    uRef = TS.u_transfer_nominal_mm(m);
    uRefNeg = TS.u_transfer_negative_sign_mm(m);
    u05 = TS.u_step05_mm(m);
    env = TS.u_transfer_abs_envelope_mm(m);
    hAbs = TS.strain_basis_abs_microstrain(m);
    kNom = median(env ./ max(hAbs, eps), 'omitnan');
    rows(k).window_id = T.window_id(k);
    rows(k).time_start_s = tr(1);
    rows(k).time_end_s = tr(2);
    if any(m)
        rows(k).time_overlap_s = max(0, max(TS.time_s(m), [], 'omitnan') - ...
            min(TS.time_s(m), [], 'omitnan'));
    else
        rows(k).time_overlap_s = 0;
    end
    rows(k).window_center_time_s = T.window_center_time_s(k);
    rows(k).Step05_A_mm = T.A_id(k);
    rows(k).Step05_EO = T.EO_id(k);
    rows(k).Step05_fn_hz = T.fn_id(k);
    rows(k).reference_region_id = mode_or_nan_local(TS.region_id(m));
    rows(k).reference_EO = mode_or_nan_local(TS.dominant_order(m));
    rsIdx = find(RS.region_id == rows(k).reference_region_id, 1);
    if nnz(m) >= 8
        rows(k).reference_A_peak_equiv_mm = kNom .* median(hAbs, 'omitnan');
        rows(k).reference_abs_envelope_median_mm = median(env, 'omitnan');
        rows(k).reference_abs_envelope_p95_mm = percentile_local(env, 95);
        rows(k).reference_A_definition = ...
            "window_overlap_frequency_amplitude_median_abs_demod";
    elseif ~isempty(rsIdx)
        rows(k).reference_A_peak_equiv_mm = RS.u_nominal_A_peak_equiv_mm(rsIdx);
        rows(k).reference_abs_envelope_median_mm = ...
            RS.u_abs_envelope_median_mm(rsIdx);
        rows(k).reference_abs_envelope_p95_mm = ...
            RS.u_abs_envelope_p95_mm(rsIdx);
        rows(k).reference_A_definition = "region_stable_core_fallback";
    else
        rows(k).reference_A_peak_equiv_mm = median(env, 'omitnan');
        rows(k).reference_abs_envelope_median_mm = median(env, 'omitnan');
        rows(k).reference_abs_envelope_p95_mm = percentile_local(env, 95);
        rows(k).reference_A_definition = "window_envelope_fallback";
    end
    rows(k).Step05_minus_reference_A_mm = ...
        rows(k).Step05_A_mm - rows(k).reference_A_peak_equiv_mm;
    rows(k).Step05_minus_reference_A_percent = ...
        100 .* rows(k).Step05_minus_reference_A_mm ./ ...
        rows(k).reference_A_peak_equiv_mm;
    rows(k).rmse_plus_sign_mm = rmse_local(u05, uRef);
    rows(k).rmse_negative_sign_mm = rmse_local(u05, uRefNeg);
    rows(k).rmse_best_sign_mm = min(rows(k).rmse_plus_sign_mm, ...
        rows(k).rmse_negative_sign_mm);
    rows(k).best_sign = string(ternary_local( ...
        rows(k).rmse_plus_sign_mm <= rows(k).rmse_negative_sign_mm, ...
        '+', '-'));
    rows(k).corr_plus_sign = corr_local(u05, uRef);
    rows(k).corr_negative_sign = corr_local(u05, uRefNeg);
end
W = struct2table(rows);
end


function row = empty_window_row_local()
row = struct('window_id', NaN, 'time_start_s', NaN, 'time_end_s', NaN, ...
    'time_overlap_s', NaN, 'window_center_time_s', NaN, ...
    'Step05_A_mm', NaN, 'Step05_EO', NaN, 'Step05_fn_hz', NaN, ...
    'reference_region_id', NaN, 'reference_EO', NaN, ...
    'reference_A_definition', "", ...
    'reference_A_peak_equiv_mm', NaN, ...
    'reference_abs_envelope_median_mm', NaN, ...
    'reference_abs_envelope_p95_mm', NaN, ...
    'Step05_minus_reference_A_mm', NaN, ...
    'Step05_minus_reference_A_percent', NaN, ...
    'rmse_plus_sign_mm', NaN, 'rmse_negative_sign_mm', NaN, ...
    'rmse_best_sign_mm', NaN, 'best_sign', "", ...
    'corr_plus_sign', NaN, 'corr_negative_sign', NaN);
end


function Summary = build_step06i_transfer_summary_local(caseName, dataset, TP, RS, ...
    W, Step05, files, cfg)
Summary = table();
Summary.case_name = string(caseName);
Summary.dataset = string(dataset);
Summary.reference_definition = ...
    "Step06A_stable_core_frequency_amplitude_plus_same_APDL_G1_transfer_ratio";
Summary.reference_uses_BTT = false;
Summary.reference_uses_Step05 = false;
Summary.step05_usage = "comparison_only";
Summary.target_blade = cfg.step06i_target_blade;
Summary.target_time_start_s = min(RS.time_start_s, [], 'omitnan');
Summary.target_time_end_s = max(RS.time_end_s, [], 'omitnan');
Summary.region_count = height(RS);
Summary.region_ids = string(mat2str(RS.region_id(:).'));
Summary.dominant_orders = string(mat2str(unique(RS.dominant_order(:)).'));
Summary.transfer_cases = strjoin(unique(string(TP.prior_case)), ';');
Summary.K_nominal_um_per_microstrain_median = median( ...
    TP.K_nominal_um_per_microstrain, 'omitnan');
Summary.K_band_low_um_per_microstrain_min = min( ...
    TP.K_band_low_um_per_microstrain, [], 'omitnan');
Summary.K_band_high_um_per_microstrain_max = max( ...
    TP.K_band_high_um_per_microstrain, [], 'omitnan');
Summary.reference_A_peak_equiv_mm_median = median( ...
    RS.u_nominal_A_peak_equiv_mm, 'omitnan');
Summary.reference_abs_envelope_median_mm = median( ...
    RS.u_abs_envelope_median_mm, 'omitnan');
Summary.reference_abs_envelope_p95_mm = median( ...
    RS.u_abs_envelope_p95_mm, 'omitnan');
Summary.step05_available = Step05.available;
Summary.step05_trend_csv = string(files.step05_trend_csv);
if ~isempty(W)
    Summary.step05_window_count = height(W);
    Summary.step05_median_A_mm = median(W.Step05_A_mm, 'omitnan');
    Summary.step05_median_EO = median(W.Step05_EO, 'omitnan');
    Summary.step05_minus_reference_A_percent_median = median( ...
        W.Step05_minus_reference_A_percent, 'omitnan');
    Summary.step05_vs_reference_best_sign_rmse_median_mm = median( ...
        W.rmse_best_sign_mm, 'omitnan');
else
    Summary.step05_window_count = 0;
    Summary.step05_median_A_mm = NaN;
    Summary.step05_median_EO = NaN;
    Summary.step05_minus_reference_A_percent_median = NaN;
    Summary.step05_vs_reference_best_sign_rmse_median_mm = NaN;
end
end


function figFile = plot_step06i_transfer_reference_local(TS, R, TP, RS, W, Summary, ...
    figDir, dataset, cfg)
fig = figure('Visible', cfg.step06i_transfer_fig_visible, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 30, 22], ...
    'Name', 'Step06I strain transfer displacement reference');
tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

blue = [0.000 0.447 0.741];
orange = [0.850 0.325 0.098];
green = [0.466 0.674 0.188];
purple = [0.494 0.184 0.556];
gray = [0.55 0.55 0.55];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
for i = 1:height(R)
    patch(ax1, [R.time_start_s(i), R.time_end_s(i), R.time_end_s(i), R.time_start_s(i)], ...
        [R.dominant_order(i)-0.28, R.dominant_order(i)-0.28, ...
        R.dominant_order(i)+0.28, R.dominant_order(i)+0.28], ...
        blue, 'FaceAlpha', 0.30, 'EdgeColor', blue, 'HandleVisibility', 'off');
    text(ax1, mean([R.time_start_s(i), R.time_end_s(i)]), ...
        R.dominant_order(i), sprintf('R%d', R.region_id(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'EO');
title(ax1, 'Step06A regions used by transfer reference');
grid(ax1, 'on');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
plot(ax2, TP.region_id, TP.K_nominal_um_per_microstrain, 'o-', ...
    'Color', purple, 'MarkerFaceColor', purple, 'DisplayName', 'G1 nominal');
plot(ax2, TP.region_id, TP.K_band_low_um_per_microstrain, ':', ...
    'Color', gray, 'DisplayName', 'surface-max lower K');
plot(ax2, TP.region_id, TP.K_band_high_um_per_microstrain, ':', ...
    'Color', gray, 'HandleVisibility', 'off');
xlabel(ax2, 'Region ID');
ylabel(ax2, 'K (um/microstrain)');
title(ax2, 'Same-APDL transfer ratio prior');
legend(ax2, 'Location', 'best');
grid(ax2, 'on');

ax3 = nexttile([1 2]); hold(ax3, 'on'); box(ax3, 'on');
plot_downsampled_local(ax3, TS.time_s, TS.u_transfer_nominal_mm, blue, ...
    'strain-transfer reference', 0.9);
plot_downsampled_local(ax3, TS.time_s, TS.u_transfer_abs_envelope_mm, ...
    orange, 'abs envelope', 0.75);
plot_downsampled_local(ax3, TS.time_s, -TS.u_transfer_abs_envelope_mm, ...
    orange, '-abs envelope', 0.75);
if any(isfinite(TS.u_step05_mm))
    plot_downsampled_local(ax3, TS.time_s, TS.u_step05_mm, gray, ...
        'Step05 nearest-window waveform', 0.55);
end
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Displacement (mm)');
title(ax3, 'Transfer-ratio displacement reference and Step05 comparison');
legend(ax3, 'Location', 'best');
grid(ax3, 'on');

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
bar(ax4, RS.region_id, [RS.u_nominal_A_peak_equiv_mm, ...
    RS.u_abs_envelope_median_mm], 0.72);
errorbar(ax4, RS.region_id, RS.u_abs_envelope_median_mm, ...
    RS.u_abs_envelope_median_mm - RS.u_band_low_abs_median_mm, ...
    RS.u_band_high_abs_median_mm - RS.u_abs_envelope_median_mm, ...
    'k.', 'LineWidth', 0.9, 'CapSize', 7, 'DisplayName', 'transfer band');
xlabel(ax4, 'Region ID');
ylabel(ax4, 'A or envelope (mm)');
title(ax4, 'Stable-core frequency amplitude by region');
legend(ax4, {'stable frequency A', 'median abs envelope', 'transfer band'}, ...
    'Location', 'best');
grid(ax4, 'on');

ax5 = nexttile; hold(ax5, 'on'); box(ax5, 'on');
if ~isempty(W)
    plot(ax5, W.window_center_time_s, W.Step05_A_mm, 'o-', ...
        'Color', gray, 'MarkerFaceColor', gray, 'DisplayName', 'Step05 A');
    plot(ax5, W.window_center_time_s, W.reference_A_peak_equiv_mm, 's-', ...
        'Color', blue, 'MarkerFaceColor', blue, 'DisplayName', 'transfer A');
    plot(ax5, W.window_center_time_s, W.reference_abs_envelope_median_mm, '^-', ...
        'Color', orange, 'MarkerFaceColor', orange, 'DisplayName', 'transfer envelope');
    yyaxis(ax5, 'right');
    plot(ax5, W.window_center_time_s, W.Step05_minus_reference_A_percent, '.-', ...
        'Color', green, 'DisplayName', 'A error %');
    ylabel(ax5, 'Step05 - reference (%)');
    yyaxis(ax5, 'left');
    legend(ax5, 'Location', 'best');
else
    text(ax5, 0.5, 0.5, 'No Step05 comparison file', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel(ax5, 'Time (s)');
ylabel(ax5, 'Amplitude (mm)');
title(ax5, 'Step05 window comparison');
grid(ax5, 'on');

sg = sgtitle(fig, sprintf(['Step06I %s: strain -> FE/APDL transfer -> ', ...
    'tip displacement (no BTT fit)\nK=%.3f um/microstrain, reference A=%.3f mm, Step05 A=%.3f mm'], ...
    dataset, Summary.K_nominal_um_per_microstrain_median(1), ...
    Summary.reference_A_peak_equiv_mm_median(1), Summary.step05_median_A_mm(1)));
set(sg, 'FontName', 'Times New Roman', 'FontSize', 10.5);
style_figure_local(fig);

ensure_dir_local(figDir);
figFile = fullfile(figDir, sprintf( ...
    'Step06I_StrainTransferDisplacementReference_%s.png', dataset));
save_figure_local(fig, figFile);
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


function oprTimes = load_opr_times_local(file)
S = load(file, 'jiluOPR');
if ~isfield(S, 'jiluOPR')
    error('OPR MAT does not contain jiluOPR:\n  %s', file);
end
oprTimes = S.jiluOPR(:, 1);
oprTimes = oprTimes(:);
oprTimes = oprTimes(isfinite(oprTimes));
if any(diff(oprTimes) <= 0)
    oprTimes = sort(oprTimes);
end
end


function thetaRot = map_time_to_rotor_phase_local(oprTimes, sampleTimes, eventsPerRev)
sampleTimes = sampleTimes(:);
thetaRot = nan(size(sampleTimes));
if isempty(eventsPerRev) || ~isfinite(eventsPerRev) || eventsPerRev < 1
    eventsPerRev = 1;
end
eventsPerRev = max(1, round(eventsPerRev));
if numel(oprTimes) <= eventsPerRev
    return;
end
revAnchorTimes = oprTimes(1:eventsPerRev:end);
revAnchorPhase = 2*pi*(0:numel(revAnchorTimes)-1).';
thetaRot(:) = interp1(revAnchorTimes(:), revAnchorPhase(:), ...
    sampleTimes, 'linear', 'extrap');
thetaRot(sampleTimes < oprTimes(1) | sampleTimes > oprTimes(end)) = NaN;
end


function y = fft_lowpass_complex_local(x, fs, cutoffHz)
x = x(:);
valid = isfinite(real(x)) & isfinite(imag(x));
if nnz(valid) < 8
    y = nan(size(x));
    return;
end
fill = x;
bad = ~valid;
if any(bad)
    fill(bad) = interp1(find(valid), x(valid), find(bad), 'linear', 'extrap');
end
n = numel(fill);
X = fft(fill);
f = (0:n-1).' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
mask = abs(f) <= cutoffHz;
y = ifft(X .* mask);
y(~valid) = NaN;
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


function value = get_table_value_or_default_local(T, row, name, defaultValue)
value = defaultValue;
if ismember(name, T.Properties.VariableNames)
    value = T.(name)(row);
end
end


function m = mode_or_nan_local(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mode(x);
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


function p = percentile_local(x, q)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    p = NaN;
    return;
end
x = sort(x);
q = max(0, min(100, q));
pos = 1 + (numel(x)-1) * q / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    p = x(lo);
else
    w = pos - lo;
    p = (1-w) * x(lo) + w * x(hi);
end
end


function value = get_optional_field_local(S, name, defaultValue)
value = defaultValue;
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
end
end


function out = ternary_local(cond, a, b)
if cond
    out = a;
else
    out = b;
end
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

