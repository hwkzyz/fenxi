%% Step06B_Validate_Step05_With_Strain_TARC_20241106.m
% Validate Step05 direct-template identification using Step06A strain/RPM evidence.
%
% Position in workflow:
%   Step06A:
%       method-neutral strain/RPM/resonance evidence.
%
%   Step06B:
%       reads Step05 direct-template identification results;
%       reads Step06A external evidence;
%       computes frequency consistency, strain synchronous amplitude,
%       TARC/TRAC, response curve comparison, and representative waveform.
%
% This script does NOT rerun Step05 and does NOT redetect resonance regions.

clc; clear; close all;

%% 0. Config
if exist('BTTDataConfig_20241106', 'file') == 2
    cfg = BTTDataConfig_20241106();
elseif exist('Get_20241106_BTT_Config', 'file') == 2
    cfg = Get_20241106_BTT_Config();
else
    error('Cannot find BTTDataConfig_20241106.m or Get_20241106_BTT_Config.m.');
end

script_dir = fileparts(mfilename('fullpath'));

if ~isfield(cfg, 'dataset') || isempty(cfg.dataset)
    cfg.dataset = '20241106';
end

if ~isfield(cfg, 'route_dir') || isempty(cfg.route_dir)
    cfg.route_dir = script_dir;
end

if ~isfield(cfg, 'output_root') || isempty(cfg.output_root)
    cfg.output_root = fullfile(script_dir, 'output');
end

if ~isfield(cfg, 'figure_root') || isempty(cfg.figure_root)
    cfg.figure_root = fullfile(cfg.output_root, 'figures');
end

if ~isfield(cfg, 'save_figures') || isempty(cfg.save_figures)
    cfg.save_figures = true;
end

if ~isfield(cfg, 'dynamic_cases') || isempty(cfg.dynamic_cases)
    if isfield(cfg, 'default_dynamic_cases') && ~isempty(cfg.default_dynamic_cases)
        cfg.dynamic_cases = cfg.default_dynamic_cases;
    else
        cfg.dynamic_cases = {'3000_3150'};
    end
end

if ischar(cfg.dynamic_cases) || isstring(cfg.dynamic_cases)
    cfg.dynamic_cases = cellstr(cfg.dynamic_cases);
end

if ~isfield(cfg, 'sample_rate_hz') || isempty(cfg.sample_rate_hz)
    if isfield(cfg, 'pinlv') && ~isempty(cfg.pinlv)
        cfg.sample_rate_hz = cfg.pinlv;
    else
        cfg.sample_rate_hz = 5e6;
    end
end

%% 1. Run control
target_cases = cfg.dynamic_cases;

target_blades = cfg.step05_target_blades;
analysis_sensors = cfg.step05_analysis_sensors;
sensor_tag = ['S', sprintf('%d', analysis_sensors)];

show_plots = true;

show_plots_env = strtrim(getenv('STEP06B_SHOW_PLOTS'));
if ~isempty(show_plots_env)
    show_plots = any(strcmpi(show_plots_env, {'1', 'true', 'yes', 'on'}));
end

if show_plots
    set(groot, 'DefaultFigureVisible', 'on');
end

% Make force_rebuild configurable instead of hardcoded constant
if ~isfield(cfg, 'force_rebuild') || isempty(cfg.force_rebuild)
    force_rebuild = false;  % Default to use cache
else
    force_rebuild = cfg.force_rebuild;
end

force_rebuild_env = strtrim(getenv('STEP06B_FORCE_REBUILD'));
if ~isempty(force_rebuild_env)
    force_rebuild = any(strcmpi(force_rebuild_env, {'1', 'true', 'yes', 'on'}));
end

%% 2. Input/output folders
step06a_output_root = cfg.step06a_output_dir;
step06b_output_root = cfg.step06b_output_dir;
step06b_figure_root = cfg.step06b_figure_dir;

step06b_output_tag_env = strtrim(getenv('STEP06B_OUTPUT_TAG'));
if ~isempty(step06b_output_tag_env)
    step06b_output_tag_env = regexprep(step06b_output_tag_env, '[^\w\-]', '_');
    step06b_output_root = fullfile(cfg.output_root, ...
        ['step06b_step05_strain_tarc_validation_', step06b_output_tag_env]);
    step06b_figure_root = fullfile(cfg.figure_root, ...
        ['step06b_step05_strain_tarc_validation_', step06b_output_tag_env]);
end

% Step05 result search roots.
step05_search_roots = { ...
    cfg.step05_output_dir, ...
    fullfile(cfg.output_root, 'step05_single_sync_direct_template'), ...
    fullfile(cfg.output_root, 'step05_single_sync_direct_template_identification'), ...
    fullfile(cfg.output_root, 'step05_single_sync_directtemplate_identification'), ...
    fullfile(cfg.output_root, 'step05'), ...
    cfg.output_root};

step05_root_env = strtrim(getenv('STEP06B_STEP05_ROOT'));
if ~isempty(step05_root_env)
    step05_search_roots = [{step05_root_env}, step05_search_roots];
end

%% 3. Validation settings
ValidationSetting = struct();

ValidationSetting.target_cases = target_cases;
ValidationSetting.target_blades = target_blades;
ValidationSetting.analysis_sensors = analysis_sensors;
ValidationSetting.sensor_tag = sensor_tag;

% Step12-style strain FFT validation: use the selected resonance region as
% the strain spectral reference, not each short Step05 fitting window.
ValidationSetting.strain_peak_search_hz = cfg.step06_raw_fft_peak_band_hz;
ValidationSetting.strain_fft_nfft_min = 4096;
ValidationSetting.use_step12_style_reference = true;
ValidationSetting.reference_region_rank = 1;
ValidationSetting.validate_only_reference_region = true;
ValidationSetting.min_reference_overlap_ratio = 0.50;

ValidationSetting.use_equivalent_displacement = true;
% 1.266 microstrain / micrometer, therefore:
% 1 microstrain = 1 / 1.266 micrometer = 1 / 1.266 / 1000 mm.
ValidationSetting.strain_to_eq_disp_mm_per_microstrain = 1 / 1.266 / 1000;
ValidationSetting.eq_disp_calibration_note = ...
    'Equivalent displacement conversion: 1.266 microstrain per micrometer.';

% Step05 frequency tolerance for synchronous amplitude extraction.
ValidationSetting.sync_amp_freq_half_width_hz = 3.0;

% TARC/TRAC settings.
ValidationSetting.time_shift_search_s = -0.003:0.00005:0.003;
ValidationSetting.max_tarc_points = 20000;
ValidationSetting.min_tarc_points = 16;

% Representative window selection.
% Options:
%   'max_strain_sync_amp'
%   'max_best_TARC'
%   'min_freq_error'
ValidationSetting.representative_window_rule = 'max_strain_sync_amp';

% Channel sensitivity.
ValidationSetting.do_channel_sensitivity = true;
ValidationSetting.channel_sensitivity_mode = 'representative_window';  % 'representative_window' or 'all_valid_windows'
ValidationSetting.channel_file_patterns = [cfg.step06_strain_files, {'AI1-*.mat'}];
ValidationSetting.channel_max_files = inf;

% Plot.
ValidationSetting.save_figures = cfg.save_figures;
ValidationSetting.fig_visible = ternary_local(show_plots, 'on', 'off');

fprintf('\n=== Step06B: validate Step05 with strain/TARC ===\n');
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Target cases: %s\n', strjoin(cellstr(target_cases), ', '));
fprintf('Target blades: %s\n', mat2str(target_blades));
fprintf('Sensor tag: %s\n', sensor_tag);

%% 4. Main loop
for iCase = 1:numel(target_cases)
    caseName = char(target_cases{iCase});

    fprintf('\n--- Step06B case: %s ---\n', caseName);

    %% 4.1 Load Step06A evidence with error handling
    try
        [Evidence, evidenceFile] = load_step06a_evidence_local(step06a_output_root, caseName);
        fprintf('Loaded Step06A evidence:\n  %s\n', evidenceFile);

        % Validate Evidence schema
        if ~validate_evidence_schema_local(Evidence)
            warning('Evidence schema validation failed for case %s. Skipping.', caseName);
            continue;
        end
    catch ME
        warning('Failed to load Step06A evidence for case %s: %s', caseName, ME.message);
        continue;
    end

    %% 4.2 Validate each blade
    for ib = 1:numel(target_blades)
        bladeId = target_blades(ib);

        fprintf('\nValidate blade %d, sensor group %s.\n', bladeId, sensor_tag);

        outDir = fullfile(step06b_output_root, caseName);
        figDir = fullfile(step06b_figure_root, caseName);

        ensure_dir_local(outDir);
        ensure_dir_local(figDir);

        outMat = fullfile(outDir, sprintf( ...
            'Step06B_Validate_Step05_With_Strain_TARC_B%d_%s_20241106.mat', ...
            bladeId, sensor_tag));

        if isfile(outMat) && ~force_rebuild
            fprintf(['Existing Step06B result found, but figures will be refreshed ', ...
                'instead of skipping:\n  %s\n'], outMat);
        end

        %% 4.2.1 Find and load Step05 result with error handling
        try
            step05File = find_step05_result_file_local( ...
                step05_search_roots, caseName, bladeId, sensor_tag);
            fprintf('Loaded Step05 result:\n  %s\n', step05File);

            S05 = load_step05_result_local(step05File);
            loadedBlade = get_step05_target_blade_local(S05);
            if isfinite(loadedBlade) && loadedBlade ~= bladeId
                error('Loaded Step05 result blade mismatch: requested B%d but file contains B%d: %s', ...
                    bladeId, loadedBlade, step05File);
            end

            Step05WindowTable = extract_step05_window_table_local( ...
                S05, bladeId, sensor_tag, Evidence);

            if isempty(Step05WindowTable) || height(Step05WindowTable) == 0
                warning('No Step05 window result extracted for blade %d %s.', bladeId, sensor_tag);
                continue;
            end

            %% 4.2.2 Build window-level validation table
            WindowValidationTable = validate_windows_with_strain_local( ...
                Step05WindowTable, Evidence, ValidationSetting);

            %% 4.2.3 Build tables
            TARC_Table = build_tarc_table_local(WindowValidationTable);
            ResponseAmpComparisonTable = build_response_amp_comparison_local( ...
                WindowValidationTable, ValidationSetting);
            WindowStrainStatsTable = build_window_strain_stats_table_local(WindowValidationTable);
            ValidationSummary = build_validation_summary_local( ...
                WindowValidationTable, ResponseAmpComparisonTable, bladeId, sensor_tag);

            %% 4.2.4 Representative window
            RepresentativeWindow = choose_representative_window_local( ...
                WindowValidationTable, ValidationSetting.representative_window_rule);

            %% 4.2.5 Channel sensitivity
            if ValidationSetting.do_channel_sensitivity
                ChannelSensitivityTable = build_channel_sensitivity_local( ...
                    Evidence, cfg, ValidationSetting, RepresentativeWindow);
            else
                ChannelSensitivityTable = table();
            end

        catch ME
            warning('Validation failed for case %s blade %d: %s', caseName, bladeId, ME.message);
            continue;
        end

        %% 4.2.6 Save
        Step06B_Result = struct();

        Step06B_Result.Dataset = cfg.dataset;
        Step06B_Result.CaseName = caseName;
        Step06B_Result.BladeId = bladeId;
        Step06B_Result.SensorTag = sensor_tag;

        Step06B_Result.CreatedBy = mfilename;
        Step06B_Result.CreatedOn = datestr(now, 31);

        Step06B_Result.Step05File = step05File;
        Step06B_Result.Step06AEvidenceFile = evidenceFile;

        Step06B_Result.Step05WindowTable = Step05WindowTable;
        Step06B_Result.WindowValidationTable = WindowValidationTable;
        Step06B_Result.TARC_Table = TARC_Table;
        Step06B_Result.ResponseAmpComparisonTable = ResponseAmpComparisonTable;
        Step06B_Result.WindowStrainStatsTable = WindowStrainStatsTable;
        Step06B_Result.ChannelSensitivityTable = ChannelSensitivityTable;
        Step06B_Result.ValidationSummary = ValidationSummary;
        Step06B_Result.RepresentativeWindow = RepresentativeWindow;
        Step06B_Result.ValidationSetting = ValidationSetting;

        save(outMat, ...
            'Step06B_Result', ...
            'Step05WindowTable', ...
            'WindowValidationTable', ...
            'TARC_Table', ...
            'ResponseAmpComparisonTable', ...
            'WindowStrainStatsTable', ...
            'ChannelSensitivityTable', ...
            'ValidationSummary', ...
            'RepresentativeWindow', ...
            'ValidationSetting', ...
            '-v7.3');

        writetable(WindowStrainStatsTable, fullfile(outDir, sprintf( ...
            'Step06B_WindowStrainStats_B%d_%s_20241106.csv', bladeId, sensor_tag)));

        writetable(TARC_Table, fullfile(outDir, sprintf( ...
            'Step06B_TARC_Table_B%d_%s_20241106.csv', bladeId, sensor_tag)));

        writetable(ResponseAmpComparisonTable, fullfile(outDir, sprintf( ...
            'Step06B_ResponseAmp_Comparison_B%d_%s_20241106.csv', bladeId, sensor_tag)));

        writetable(ChannelSensitivityTable, fullfile(outDir, sprintf( ...
            'Step06B_ChannelSensitivity_B%d_%s_20241106.csv', bladeId, sensor_tag)));

        writetable(ValidationSummary, fullfile(outDir, sprintf( ...
            'Step06B_ValidationSummary_B%d_%s_20241106.csv', bladeId, sensor_tag)));

        fprintf('Saved Step06B result:\n  %s\n', outMat);

        %% 4.2.7 Plot
        if show_plots || ValidationSetting.save_figures
            plot_step06b_main_validation_local( ...
                Evidence, ...
                WindowValidationTable, ...
                ResponseAmpComparisonTable, ...
                RepresentativeWindow, ...
                ValidationSetting, ...
                figDir, ...
                ValidationSetting.fig_visible, ...
                ValidationSetting.save_figures, ...
                bladeId, ...
                sensor_tag);

            plot_step06b_tarc_trend_local( ...
                WindowValidationTable, ...
                figDir, ...
                ValidationSetting.fig_visible, ...
                ValidationSetting.save_figures, ...
                bladeId, ...
                sensor_tag);

            plot_step06b_channel_sensitivity_local( ...
                ChannelSensitivityTable, ...
                figDir, ...
                ValidationSetting.fig_visible, ...
                ValidationSetting.save_figures, ...
                bladeId, ...
                sensor_tag);

            drawnow;
            fprintf('Step06B figures refreshed:\n  %s\n', figDir);
        end
    end
end

fprintf('\n=== Step06B finished ===\n');

%% ========================================================================
% Local functions
% ========================================================================

function [Evidence, evidenceFile] = load_step06a_evidence_local(step06aRoot, caseName)

evidenceFile = fullfile(step06aRoot, caseName, ...
    'Step06A_StrainRPM_ResonanceEvidence_20241106.mat');

if ~isfile(evidenceFile)
    dd = dir(fullfile(step06aRoot, '**', 'Step06A_StrainRPM_ResonanceEvidence_20241106.mat'));

    evidenceFile = '';

    for k = 1:numel(dd)
        f = fullfile(dd(k).folder, dd(k).name);
        if contains(f, caseName)
            evidenceFile = f;
            break;
        end
    end
end

if isempty(evidenceFile) || ~isfile(evidenceFile)
    error('Missing Step06A evidence for case %s. Run Step06A first.', caseName);
end

D = load(evidenceFile);

if isfield(D, 'Evidence')
    Evidence = D.Evidence;
else
    Evidence = struct();

    if isfield(D, 'rawFft')
        Evidence.rawFft = D.rawFft;
        Evidence.RawFft = D.rawFft;
    end

    if isfield(D, 'regionTable')
        Evidence.regionTable = D.regionTable;
        Evidence.RegionTable = D.regionTable;
    end

    if isfield(D, 'scoreTable')
        Evidence.scoreTable = D.scoreTable;
        Evidence.ScoreTable = D.scoreTable;
    end

    if isfield(D, 'localFftTable')
        Evidence.localFftTable = D.localFftTable;
        Evidence.LocalFftTable = D.localFftTable;
    end

    if isfield(D, 'strainSource')
        Evidence.StrainSource = D.strainSource;
    end
end

Evidence = normalize_evidence_fields_local(Evidence);

end


function Evidence = normalize_evidence_fields_local(Evidence)

if ~isfield(Evidence, 'RegionTable')
    if isfield(Evidence, 'regionTable')
        Evidence.RegionTable = Evidence.regionTable;
    else
        Evidence.RegionTable = table();
    end
end

if ~isfield(Evidence, 'regionTable')
    Evidence.regionTable = Evidence.RegionTable;
end

if ~isfield(Evidence, 'RawFft')
    if isfield(Evidence, 'rawFft')
        Evidence.RawFft = Evidence.rawFft;
    else
        Evidence.RawFft = struct();
    end
end

if ~isfield(Evidence, 'rawFft')
    Evidence.rawFft = Evidence.RawFft;
end

if ~isfield(Evidence, 'StrainRaw')
    if isfield(Evidence, 'StrainRef')
        Evidence.StrainRaw.timeBtt = Evidence.StrainRef.time_s;
        Evidence.StrainRaw.rawValue = Evidence.StrainRef.raw_microstrain;
        Evidence.StrainRaw.detrendedValue = Evidence.StrainRef.detrended_microstrain;
        Evidence.StrainRaw.sampleRateHz = Evidence.StrainRef.sample_rate_hz;
        Evidence.StrainRaw.file = Evidence.StrainRef.strain_file;
        Evidence.StrainRaw.timeOffsetSec = Evidence.StrainRef.time_offset_sec;
    elseif isfield(Evidence, 'StrainTimeBtt') && isfield(Evidence, 'StrainValue')
        Evidence.StrainRaw.timeBtt = Evidence.StrainTimeBtt;
        Evidence.StrainRaw.rawValue = Evidence.StrainValue;
        Evidence.StrainRaw.detrendedValue = Evidence.StrainValue;
        Evidence.StrainRaw.sampleRateHz = estimate_sample_rate_local(Evidence.StrainTimeBtt);
        Evidence.StrainRaw.file = '';
        Evidence.StrainRaw.timeOffsetSec = 0;
    else
        error('Evidence does not contain StrainRaw or compatible strain fields.');
    end
end

if ~isfield(Evidence, 'StrainTimeBtt')
    Evidence.StrainTimeBtt = Evidence.StrainRaw.timeBtt;
end

if ~isfield(Evidence, 'StrainValue')
    Evidence.StrainValue = Evidence.StrainRaw.detrendedValue;
end

if ~isfield(Evidence, 'RPM')
    error('Evidence does not contain RPM.');
end

if ~isfield(Evidence.RPM, 'timeSec')
    if isfield(Evidence.RPM, 'time_s')
        Evidence.RPM.timeSec = Evidence.RPM.time_s;
    else
        error('Evidence.RPM lacks timeSec/time_s.');
    end
end

if ~isfield(Evidence.RPM, 'time_s')
    Evidence.RPM.time_s = Evidence.RPM.timeSec;
end

if ~isfield(Evidence.RPM, 'value')
    error('Evidence.RPM lacks value.');
end

if ~isfield(Evidence, 'StrainSource')
    Evidence.StrainSource = struct();
    Evidence.StrainSource.strain_file = get_optional_field_local(Evidence.StrainRaw, 'file', '');
    Evidence.StrainSource.strain_time_offset_total_sec = get_optional_field_local(Evidence.StrainRaw, 'timeOffsetSec', 0);
    Evidence.StrainSource.source_mode = 'from_Evidence_StrainRaw';
end

end


function step05File = find_step05_result_file_local(searchRoots, caseName, bladeId, sensorTag)

patterns = { ...
    sprintf('Result_Step05_SingleSyncDirectTemplate_B%d_%s_20241106.mat', bladeId, sensorTag), ...
    sprintf('*Step05*B%d*%s*.mat', bladeId, sensorTag), ...
    sprintf('*SingleSync*B%d*%s*.mat', bladeId, sensorTag), ...
    sprintf('*DirectTemplate*B%d*%s*.mat', bladeId, sensorTag), ...
    sprintf('*B%d*%s*20241106*.mat', bladeId, sensorTag)};

files = {};

for ir = 1:numel(searchRoots)
    root = searchRoots{ir};

    if ~isfolder(root)
        continue;
    end

    for ip = 1:numel(patterns)
        dd = dir(fullfile(root, '**', patterns{ip}));

        for k = 1:numel(dd)
            f = fullfile(dd(k).folder, dd(k).name);
            nameNoExt = erase(dd(k).name, '.mat');

            if contains(lower(f), 'step06')
                continue;
            end

            if ~contains(nameNoExt, sprintf('B%d_', bladeId)) && ...
                    ~contains(nameNoExt, sprintf('B%d-', bladeId)) && ...
                    ~endsWith(nameNoExt, sprintf('B%d', bladeId))
                continue;
            end

            if ~contains(f, caseName) && ~contains(dd(k).folder, caseName)
                continue;
            end

            files{end+1} = f; %#ok<AGROW>
        end
    end
end

files = unique(files, 'stable');

if isempty(files)
    error('Cannot find Step05 result for case %s, blade %d, sensor %s.', ...
        caseName, bladeId, sensorTag);
end

% Prefer newest file.
datenumList = nan(numel(files), 1);

for i = 1:numel(files)
    d = dir(files{i});
    datenumList(i) = d.datenum;
end

[~, im] = max(datenumList);
step05File = files{im};

end


function S05 = load_step05_result_local(step05File)

D = load(step05File);

preferred = { ...
    'ResultSet', ...
    'Step05_Result', ...
    'Step05Result', ...
    'Result', ...
    'result', ...
    'Results'};

S05 = [];

for i = 1:numel(preferred)
    if isfield(D, preferred{i})
        S05 = D.(preferred{i});
        break;
    end
end

if isempty(S05)
    fn = fieldnames(D);
    S05 = D.(fn{1});
end

end


function bladeId = get_step05_target_blade_local(S05)

bladeId = NaN;

if isstruct(S05)
    candidates = {'TargetBlade', 'target_blade', 'blade_id', 'BladeID'};
    for i = 1:numel(candidates)
        if isfield(S05, candidates{i}) && isnumeric(S05.(candidates{i})) && ...
                isscalar(S05.(candidates{i}))
            bladeId = double(S05.(candidates{i}));
            return;
        end
    end

    if isfield(S05, 'Cases') && isstruct(S05.Cases) && ...
            isfield(S05.Cases, 'BladeResult') && ~isempty(S05.Cases(1).BladeResult)
        bladeId = get_step05_target_blade_local(S05.Cases(1).BladeResult(1).Result);
    end
end

end


function T = extract_step05_window_table_local(S05, bladeId, sensorTag, Evidence)

WR = find_window_result_array_local(S05);

if isempty(WR)
    error('Cannot find WindowResult/windowResult in Step05 result.');
end

if istable(WR)
    T = extract_step05_window_table_from_table_local(WR, bladeId, sensorTag, Evidence);
    return;
end

if ~isstruct(WR)
    error('Step05 WindowResult is neither struct nor table.');
end

n = numel(WR);

window_id = nan(n,1);
time_start_s = nan(n,1);
time_end_s = nan(n,1);
time_center_s = nan(n,1);

blade_id = bladeId * ones(n,1);
sensor_group = repmat(string(sensorTag), n, 1);

EO_id = nan(n,1);
fn_id_hz = nan(n,1);
A_id_mm = nan(n,1);
phi_rad = nan(n,1);
dx_c_mm = nan(n,1);
eta_s_mm = nan(n,1);

weighted_voltage_rmse_mv = nan(n,1);
fit_stage = strings(n,1);

valid_point_count = nan(n,1);
core_point_count = nan(n,1);
expanded_point_count = nan(n,1);

rpm_mean = nan(n,1);

for i = 1:n
    wr = WR(i);

    window_id(i) = get_numeric_field_local(wr, {'window_id','windowId','id','WindowID'}, i);

    tw = get_numeric_vector_field_local(wr, { ...
        'time_window_s', ...
        'timeWindow', ...
        'time_window', ...
        'timeWindowSec', ...
        'window_time_s'});

    if numel(tw) >= 2
        time_start_s(i) = tw(1);
        time_end_s(i) = tw(2);
    else
        time_start_s(i) = get_numeric_field_local(wr, {'time_start_s','t_start_s','start_time_s','window_start_s'}, NaN);
        time_end_s(i) = get_numeric_field_local(wr, {'time_end_s','t_end_s','end_time_s','window_end_s'}, NaN);
    end

    if isfinite(time_start_s(i)) && isfinite(time_end_s(i))
        time_center_s(i) = 0.5 * (time_start_s(i) + time_end_s(i));
    else
        time_center_s(i) = get_numeric_field_local(wr, {'time_center_s','t_center_s','center_time_s'}, NaN);
    end

    R = get_struct_field_local(wr, {'Result','result','Fit','fit','BestResult','bestResult'}, struct());

    EO_id(i) = get_numeric_field_multi_local(wr, R, {'EO_id','EO','eo','EO_best','best_EO'}, NaN);
    fn_id_hz(i) = get_numeric_field_multi_local(wr, R, {'fn_id','freqHz','freq_hz','frequency_hz','f_id_hz','fn_hz'}, NaN);
    A_id_mm(i) = get_numeric_field_multi_local(wr, R, {'A_id','amplitudeMm','amplitude_mm','A_mm','A'}, NaN);
    phi_rad(i) = get_numeric_field_multi_local(wr, R, {'phi_id','phi_rad','phaseRad','phase_rad','phi'}, NaN);
    dx_c_mm(i) = get_numeric_field_multi_local(wr, R, {'dx_c','dx_c_mm','dx0','dx0_mm'}, NaN);
    eta_s_mm(i) = get_numeric_field_multi_local(wr, R, {'eta_s','eta_s_mm','sensor_eta_mm','eta'}, NaN);

    rmseMv = get_numeric_field_multi_local(wr, R, {'weighted_voltage_rmse_mv','rmse_mv','weighted_rmse_mv'}, NaN);

    if ~isfinite(rmseMv)
        rmseV = get_numeric_field_multi_local(wr, R, {'weighted_voltage_rmse','rmse_v','weighted_rmse'}, NaN);
        if isfinite(rmseV)
            rmseMv = rmseV * 1000;
        end
    end

    weighted_voltage_rmse_mv(i) = rmseMv;

    fit_stage(i) = string(get_char_field_multi_local(wr, R, {'fit_stage','stage','fitStage'}, ''));

    valid_point_count(i) = get_numeric_field_multi_local(wr, R, {'valid_point_count','n_valid','valid_count'}, NaN);
    core_point_count(i) = get_numeric_field_multi_local(wr, R, {'core_point_count','n_core','core_count'}, NaN);
    expanded_point_count(i) = get_numeric_field_multi_local(wr, R, {'expanded_point_count','n_expanded','expanded_count'}, NaN);

    rpm_mean(i) = get_numeric_field_multi_local(wr, R, {'rpm_mean','mean_rpm','RPM_mean'}, NaN);

    if ~isfinite(rpm_mean(i)) && isfinite(time_center_s(i))
        rpm_mean(i) = interp_rpm_local(Evidence, time_center_s(i));
    end

    if ~isfinite(fn_id_hz(i)) && isfinite(EO_id(i)) && isfinite(rpm_mean(i))
        fn_id_hz(i) = EO_id(i) * rpm_mean(i) / 60;
    end
end

T = table( ...
    window_id, ...
    time_start_s, ...
    time_end_s, ...
    time_center_s, ...
    blade_id, ...
    sensor_group, ...
    EO_id, ...
    fn_id_hz, ...
    A_id_mm, ...
    phi_rad, ...
    dx_c_mm, ...
    eta_s_mm, ...
    weighted_voltage_rmse_mv, ...
    fit_stage, ...
    valid_point_count, ...
    core_point_count, ...
    expanded_point_count, ...
    rpm_mean, ...
    'VariableNames', { ...
    'window_id', ...
    'time_start_s', ...
    'time_end_s', ...
    'time_center_s', ...
    'blade_id', ...
    'sensor_group', ...
    'EO_id', ...
    'fn_id_hz', ...
    'A_id_mm', ...
    'phi_rad', ...
    'dx_c_mm', ...
    'eta_s_mm', ...
    'weighted_voltage_rmse_mv', ...
    'fit_stage', ...
    'valid_point_count', ...
    'core_point_count', ...
    'expanded_point_count', ...
    'rpm_mean'});

validTime = isfinite(T.time_start_s) & isfinite(T.time_end_s) & T.time_end_s > T.time_start_s;
T = T(validTime, :);

end


function T = extract_step05_window_table_from_table_local(W, bladeId, sensorTag, Evidence)

n = height(W);

T = table();
T.window_id = get_table_numeric_local(W, {'window_id','windowId','id'}, (1:n).');
T.time_start_s = get_table_numeric_local(W, {'time_start_s','t_start_s','start_time_s','window_start_s'}, nan(n,1));
T.time_end_s = get_table_numeric_local(W, {'time_end_s','t_end_s','end_time_s','window_end_s'}, nan(n,1));

if ismember('time_window_s', W.Properties.VariableNames)
    tw = W.time_window_s;
    if size(tw,2) >= 2
        T.time_start_s = tw(:,1);
        T.time_end_s = tw(:,2);
    end
end

T.time_center_s = 0.5 * (T.time_start_s + T.time_end_s);

T.blade_id = bladeId * ones(n,1);
T.sensor_group = repmat(string(sensorTag), n, 1);

T.EO_id = get_table_numeric_local(W, {'EO_id','EO','eo'}, nan(n,1));
T.fn_id_hz = get_table_numeric_local(W, {'fn_id','freqHz','freq_hz','frequency_hz'}, nan(n,1));
T.A_id_mm = get_table_numeric_local(W, {'A_id','A_mm','amplitudeMm','amplitude_mm'}, nan(n,1));
T.phi_rad = get_table_numeric_local(W, {'phi_id','phi_rad','phaseRad','phase_rad'}, nan(n,1));
T.dx_c_mm = get_table_numeric_local(W, {'dx_c','dx_c_mm'}, nan(n,1));
T.eta_s_mm = get_table_numeric_local(W, {'eta_s','eta_s_mm'}, nan(n,1));

T.weighted_voltage_rmse_mv = get_table_numeric_local(W, {'weighted_voltage_rmse_mv','rmse_mv'}, nan(n,1));

if all(~isfinite(T.weighted_voltage_rmse_mv))
    rmseV = get_table_numeric_local(W, {'weighted_voltage_rmse','rmse_v'}, nan(n,1));
    T.weighted_voltage_rmse_mv = rmseV * 1000;
end

if ismember('fit_stage', W.Properties.VariableNames)
    T.fit_stage = string(W.fit_stage);
else
    T.fit_stage = strings(n,1);
end

T.valid_point_count = get_table_numeric_local(W, {'valid_point_count','n_valid'}, nan(n,1));
T.core_point_count = get_table_numeric_local(W, {'core_point_count','n_core'}, nan(n,1));
T.expanded_point_count = get_table_numeric_local(W, {'expanded_point_count','n_expanded'}, nan(n,1));
T.rpm_mean = get_table_numeric_local(W, {'rpm_mean','mean_rpm'}, nan(n,1));

for i = 1:n
    if ~isfinite(T.rpm_mean(i)) && isfinite(T.time_center_s(i))
        T.rpm_mean(i) = interp_rpm_local(Evidence, T.time_center_s(i));
    end

    if ~isfinite(T.fn_id_hz(i)) && isfinite(T.EO_id(i)) && isfinite(T.rpm_mean(i))
        T.fn_id_hz(i) = T.EO_id(i) * T.rpm_mean(i) / 60;
    end
end

validTime = isfinite(T.time_start_s) & isfinite(T.time_end_s) & T.time_end_s > T.time_start_s;
T = T(validTime, :);

end


function WR = find_window_result_array_local(S)

WR = [];

if istable(S)
    WR = S;
    return;
end

if ~isstruct(S)
    return;
end

directNames = { ...
    'WindowResult', ...
    'windowResult', ...
    'WindowResults', ...
    'windowResults', ...
    'Window_Table', ...
    'WindowTable', ...
    'window_table', ...
    'Results', ...
    'results'};

for i = 1:numel(directNames)
    name = directNames{i};

    if isfield(S, name)
        candidate = S.(name);

        if istable(candidate)
            WR = candidate;
            return;
        end

        if isstruct(candidate)
            if numel(candidate) > 1 || looks_like_window_struct_local(candidate)
                WR = candidate;
                return;
            end
        end
    end
end

if numel(S) > 1 && looks_like_window_struct_local(S)
    WR = S;
    return;
end

nestedNames = {'ResultSet','Step05_Result','Step05Result','Result','result','Output','output'};

for i = 1:numel(nestedNames)
    name = nestedNames{i};

    if isfield(S, name) && isstruct(S.(name))
        WR = find_window_result_array_local(S.(name));

        if ~isempty(WR)
            return;
        end
    end
end

end


function tf = looks_like_window_struct_local(S)

tf = false;

if isempty(S) || ~isstruct(S)
    return;
end

names = fieldnames(S);

keyNames = { ...
    'time_window_s', ...
    'timeWindow', ...
    'time_start_s', ...
    'Result', ...
    'result', ...
    'EO_id', ...
    'fn_id', ...
    'A_id'};

for i = 1:numel(keyNames)
    if any(strcmp(names, keyNames{i}))
        tf = true;
        return;
    end
end

end


function WindowValidationTable = validate_windows_with_strain_local(W, Evidence, P)

tStrain = Evidence.StrainRaw.timeBtt(:);
yStrain = Evidence.StrainRaw.detrendedValue(:);

validStrain = isfinite(tStrain) & isfinite(yStrain);

tStrain = tStrain(validStrain);
yStrain = yStrain(validStrain);

[tStrain, order] = sort(tStrain);
yStrain = yStrain(order);

R = Evidence.RegionTable;
StrainReference = build_step12_style_strain_reference_local(tStrain, yStrain, R, P);
k_eq = get_eq_disp_scale_local(P);

n = height(W);

WindowValidationTable = W;

WindowValidationTable.validated = false(n,1);
WindowValidationTable.region_id = nan(n,1);

WindowValidationTable.strain_fft_peak_hz = nan(n,1);
WindowValidationTable.strain_fft_peak_amp_microstrain = nan(n,1);
WindowValidationTable.strain_fft_peak_amp_eq_mm = nan(n,1);
WindowValidationTable.strain_reference_region_id = nan(n,1);
WindowValidationTable.strain_reference_time_start_s = nan(n,1);
WindowValidationTable.strain_reference_time_end_s = nan(n,1);
WindowValidationTable.local_strain_fft_peak_hz = nan(n,1);
WindowValidationTable.local_strain_fft_peak_amp_microstrain = nan(n,1);
WindowValidationTable.local_strain_fft_peak_amp_eq_mm = nan(n,1);
WindowValidationTable.strain_sync_amp_at_step05_fn_microstrain = nan(n,1);
WindowValidationTable.strain_sync_amp_at_step05_fn_eq_mm = nan(n,1);
WindowValidationTable.strain_sync_phase_at_step05_fn_rad = nan(n,1);

WindowValidationTable.freq_error_vs_strain_fft_peak_hz = nan(n,1);
WindowValidationTable.freq_error_vs_strain_fft_peak_percent = nan(n,1);
WindowValidationTable.freq_error_vs_local_strain_peak_hz = nan(n,1);

WindowValidationTable.strain_order_freq_hz = nan(n,1);
WindowValidationTable.freq_error_vs_order_line_hz = nan(n,1);

WindowValidationTable.strain_rms_microstrain = nan(n,1);
WindowValidationTable.strain_rms_eq_mm = nan(n,1);
WindowValidationTable.strain_half_peak_to_peak_microstrain = nan(n,1);
WindowValidationTable.strain_half_peak_to_peak_eq_mm = nan(n,1);

WindowValidationTable.raw_TARC = nan(n,1);
WindowValidationTable.best_TARC = nan(n,1);
WindowValidationTable.raw_TRAC = nan(n,1);
WindowValidationTable.best_TRAC = nan(n,1);
WindowValidationTable.best_time_shift_s = nan(n,1);
WindowValidationTable.phase_opt_TARC = nan(n,1);

WindowValidationTable.tarc_point_count = nan(n,1);

for i = 1:n
    t0 = W.time_start_s(i);
    t1 = W.time_end_s(i);
    tc = W.time_center_s(i);
    fn = W.fn_id_hz(i);

    if ~isfinite(t0) || ~isfinite(t1) || t1 <= t0 || ~isfinite(fn) || fn <= 0
        continue;
    end

    if isfield(P, 'validate_only_reference_region') && P.validate_only_reference_region && ...
            StrainReference.isValid
        overlapRatio = overlap_ratio_local([t0, t1], StrainReference.timeRange);

        minOverlap = 0.50;
        if isfield(P, 'min_reference_overlap_ratio') && isfinite(P.min_reference_overlap_ratio)
            minOverlap = P.min_reference_overlap_ratio;
        end

        if overlapRatio < minOverlap
            continue;
        end
    end

    mask = tStrain >= t0 & tStrain <= t1;

    if nnz(mask) < P.min_tarc_points
        continue;
    end

    ts = tStrain(mask);
    ys = yStrain(mask);

    ys = ys(:);
    ts = ts(:);

    finite = isfinite(ts) & isfinite(ys);
    ts = ts(finite);
    ys = ys(finite);

    if numel(ys) < P.min_tarc_points
        continue;
    end

    ys = detrend(ys, 'linear');
    ys = ys - mean(ys, 'omitnan');

    WindowValidationTable.strain_rms_microstrain(i) = rms_local(ys);
    WindowValidationTable.strain_half_peak_to_peak_microstrain(i) = ...
        0.5 * (max(ys, [], 'omitnan') - min(ys, [], 'omitnan'));

    if isfinite(k_eq)
        WindowValidationTable.strain_rms_eq_mm(i) = ...
            WindowValidationTable.strain_rms_microstrain(i) * k_eq;
        WindowValidationTable.strain_half_peak_to_peak_eq_mm(i) = ...
            WindowValidationTable.strain_half_peak_to_peak_microstrain(i) * k_eq;
    end

    if StrainReference.isValid
        pf = StrainReference.peakFreqHz;
        pamp = StrainReference.peakAmpMicrostrain;

        WindowValidationTable.strain_fft_peak_hz(i) = pf;
        WindowValidationTable.strain_fft_peak_amp_microstrain(i) = pamp;
        if isfinite(k_eq)
            WindowValidationTable.strain_fft_peak_amp_eq_mm(i) = pamp * k_eq;
        end
        WindowValidationTable.strain_reference_region_id(i) = StrainReference.regionId;
        WindowValidationTable.strain_reference_time_start_s(i) = StrainReference.timeRange(1);
        WindowValidationTable.strain_reference_time_end_s(i) = StrainReference.timeRange(2);

        WindowValidationTable.freq_error_vs_strain_fft_peak_hz(i) = W.fn_id_hz(i) - pf;

        if isfinite(pf) && pf > 0
            WindowValidationTable.freq_error_vs_strain_fft_peak_percent(i) = ...
                100 * (W.fn_id_hz(i) - pf) / pf;
        end
    end

    [localF, localA] = single_fft_amp_local(ts, ys, P.strain_peak_search_hz, P.strain_fft_nfft_min);

    if ~isempty(localF) && ~isempty(localA)
        localOk = isfinite(localF) & isfinite(localA);

        if any(localOk)
            localFok = localF(localOk);
            localAok = localA(localOk);
            [localPeakAmp, localIm] = max(localAok);
            localPeakFreq = localFok(localIm);

            WindowValidationTable.local_strain_fft_peak_hz(i) = localPeakFreq;
            WindowValidationTable.local_strain_fft_peak_amp_microstrain(i) = localPeakAmp;
            WindowValidationTable.freq_error_vs_local_strain_peak_hz(i) = ...
                W.fn_id_hz(i) - localPeakFreq;

            if isfinite(k_eq)
                WindowValidationTable.local_strain_fft_peak_amp_eq_mm(i) = localPeakAmp * k_eq;
            end
        end
    end

    [syncAmp, syncPhase, syncFit] = sine_projection_local(ts, ys, fn);

    WindowValidationTable.strain_sync_amp_at_step05_fn_microstrain(i) = syncAmp;
    WindowValidationTable.strain_sync_phase_at_step05_fn_rad(i) = syncPhase;

    if isfinite(k_eq)
        WindowValidationTable.strain_sync_amp_at_step05_fn_eq_mm(i) = syncAmp * k_eq;
    end

    if isfinite(W.EO_id(i)) && isfinite(W.rpm_mean(i))
        orderFreq = W.EO_id(i) * W.rpm_mean(i) / 60;
        WindowValidationTable.strain_order_freq_hz(i) = orderFreq;
        WindowValidationTable.freq_error_vs_order_line_hz(i) = W.fn_id_hz(i) - orderFreq;
    end

    WindowValidationTable.region_id(i) = find_region_id_for_time_local(R, tc);

    [ts2, ys2] = thin_signal_for_tarc_local(ts, ys, P.max_tarc_points);

    A = W.A_id_mm(i);
    if ~isfinite(A) || abs(A) < eps
        A = 1;
    end

    phi = W.phi_rad(i);
    if ~isfinite(phi)
        phi = 0;
    end

    rawTarc = calc_tarc_with_sine_local(ts2, ys2, fn, A, phi, 0);

    [bestTarc, bestShift] = calc_best_tarc_shift_local( ...
        ts2, ys2, fn, A, phi, P.time_shift_search_s);

    phaseOptTarc = calc_tarc_direct_local(ys, syncFit);

    WindowValidationTable.raw_TARC(i) = rawTarc;
    WindowValidationTable.best_TARC(i) = bestTarc;
    WindowValidationTable.raw_TRAC(i) = rawTarc;
    WindowValidationTable.best_TRAC(i) = bestTarc;
    WindowValidationTable.best_time_shift_s(i) = bestShift;
    WindowValidationTable.phase_opt_TARC(i) = phaseOptTarc;
    WindowValidationTable.tarc_point_count(i) = numel(ts2);

    WindowValidationTable.validated(i) = true;
end

end


function Ref = build_step12_style_strain_reference_local(tStrain, yStrain, R, P)

Ref = struct();
Ref.isValid = false;
Ref.regionId = NaN;
Ref.timeRange = [NaN NaN];
Ref.freqHz = [];
Ref.ampMicrostrain = [];
Ref.peakFreqHz = NaN;
Ref.peakAmpMicrostrain = NaN;

if isempty(R) || height(R) == 0
    return;
end

rank = 1;

if isfield(P, 'reference_region_rank') && isfinite(P.reference_region_rank)
    rank = max(1, round(P.reference_region_rank));
end

rank = min(rank, height(R));

row = R(rank, :);

if ismember('region_id', R.Properties.VariableNames)
    Ref.regionId = row.region_id;
elseif ismember('regionId', R.Properties.VariableNames)
    Ref.regionId = row.regionId;
else
    Ref.regionId = rank;
end

if ismember('time_start_s', R.Properties.VariableNames)
    t0 = row.time_start_s;
    t1 = row.time_end_s;
elseif ismember('bttStartSec', R.Properties.VariableNames)
    t0 = row.bttStartSec;
    t1 = row.bttEndSec;
else
    return;
end

Ref.timeRange = [t0, t1];

mask = tStrain >= t0 & tStrain <= t1 & isfinite(tStrain) & isfinite(yStrain);

if nnz(mask) < 16
    return;
end

ts = tStrain(mask);
ys = yStrain(mask);

ys = detrend(ys(:), 'linear');
ys = ys - mean(ys, 'omitnan');

[ff, aa] = single_fft_amp_local(ts, ys, P.strain_peak_search_hz, P.strain_fft_nfft_min);

if isempty(ff)
    return;
end

[pamp, im] = max(aa);

Ref.freqHz = ff;
Ref.ampMicrostrain = aa;
Ref.peakFreqHz = ff(im);
Ref.peakAmpMicrostrain = pamp;
Ref.isValid = isfinite(Ref.peakFreqHz) && isfinite(Ref.peakAmpMicrostrain);

end


function TARC_Table = build_tarc_table_local(W)

vars = { ...
    'window_id', ...
    'time_start_s', ...
    'time_end_s', ...
    'time_center_s', ...
    'rpm_mean', ...
    'region_id', ...
    'EO_id', ...
    'fn_id_hz', ...
    'A_id_mm', ...
    'weighted_voltage_rmse_mv', ...
    'raw_TARC', ...
    'best_TARC', ...
    'raw_TRAC', ...
    'best_TRAC', ...
    'best_time_shift_s', ...
    'phase_opt_TARC', ...
    'tarc_point_count', ...
    'validated'};

TARC_Table = W(:, select_existing_vars_local(W, vars));

end


function ResponseAmpComparisonTable = build_response_amp_comparison_local(W, P)

if isempty(W)
    ResponseAmpComparisonTable = table();
    return;
end

k_eq = get_eq_disp_scale_local(P);

if ismember('strain_sync_amp_at_step05_fn_eq_mm', W.Properties.VariableNames)
    strainEq = W.strain_sync_amp_at_step05_fn_eq_mm;
elseif isfinite(k_eq) && ismember('strain_sync_amp_at_step05_fn_microstrain', W.Properties.VariableNames)
    strainEq = W.strain_sync_amp_at_step05_fn_microstrain * k_eq;
else
    strainEq = nan(height(W), 1);
end

valid = W.validated & isfinite(W.A_id_mm) & isfinite(strainEq);

T = W(valid, :);
strainEq = strainEq(valid);

ResponseAmpComparisonTable = table();

if height(T) == 0
    return;
end

ResponseAmpComparisonTable.window_id = T.window_id;
ResponseAmpComparisonTable.time_center_s = T.time_center_s;
ResponseAmpComparisonTable.rpm_mean = T.rpm_mean;
ResponseAmpComparisonTable.region_id = T.region_id;

ResponseAmpComparisonTable.Step05_A_mm = abs(T.A_id_mm);
ResponseAmpComparisonTable.Strain_EqDispAmp_mm = abs(strainEq);

ResponseAmpComparisonTable.EqDisp_Error_mm = ...
    ResponseAmpComparisonTable.Step05_A_mm - ResponseAmpComparisonTable.Strain_EqDispAmp_mm;

ResponseAmpComparisonTable.EqDisp_AbsError_mm = abs(ResponseAmpComparisonTable.EqDisp_Error_mm);

ResponseAmpComparisonTable.EqDisp_Ratio_StrainOverStep05 = ...
    safe_divide_local(ResponseAmpComparisonTable.Strain_EqDispAmp_mm, ...
                      ResponseAmpComparisonTable.Step05_A_mm);

ResponseAmpComparisonTable.EqDisp_RelError_percent = ...
    100 * safe_divide_local(ResponseAmpComparisonTable.EqDisp_Error_mm, ...
                            ResponseAmpComparisonTable.Step05_A_mm);

ResponseAmpComparisonTable.A_id_mm = T.A_id_mm;

if ismember('strain_sync_amp_at_step05_fn_microstrain', T.Properties.VariableNames)
    ResponseAmpComparisonTable.strain_sync_amp_at_step05_fn_microstrain = ...
        T.strain_sync_amp_at_step05_fn_microstrain;
end

ResponseAmpComparisonTable.strain_sync_amp_at_step05_fn_eq_mm = strainEq;

ResponseAmpComparisonTable.Step05_A_norm = normalize_response_curve_local( ...
    ResponseAmpComparisonTable.Step05_A_mm, 'max_abs');

ResponseAmpComparisonTable.Strain_EqDispAmp_norm = normalize_response_curve_local( ...
    ResponseAmpComparisonTable.Strain_EqDispAmp_mm, 'max_abs');

if ismember('best_TARC', T.Properties.VariableNames)
    ResponseAmpComparisonTable.best_TARC = T.best_TARC;
end

if ismember('phase_opt_TARC', T.Properties.VariableNames)
    ResponseAmpComparisonTable.phase_opt_TARC = T.phase_opt_TARC;
end

end


function k = get_eq_disp_scale_local(P)

k = NaN;

if isstruct(P) && isfield(P, 'use_equivalent_displacement') && ...
        ~P.use_equivalent_displacement
    return;
end

if isstruct(P) && ...
        isfield(P, 'strain_to_eq_disp_mm_per_microstrain') && ...
        isfinite(P.strain_to_eq_disp_mm_per_microstrain)

    k = P.strain_to_eq_disp_mm_per_microstrain;
    return;
end

% Backward compatibility with older Step06B parameter files.
if isstruct(P) && isfield(P, 'strain_to_mm') && isfinite(P.strain_to_mm)
    k = P.strain_to_mm;
end

end


function y = safe_divide_local(a, b)

y = nan(size(a));
ok = isfinite(a) & isfinite(b) & abs(b) > eps;
y(ok) = a(ok) ./ b(ok);

end


function y = normalize_response_curve_local(x, mode)

if nargin < 2 || isempty(mode)
    mode = 'max_abs';
end

x = x(:);
y = nan(size(x));

switch lower(mode)
    case 'max_abs'
        denom = max(abs(x), [], 'omitnan');

        if isfinite(denom) && denom > 0
            y = x ./ denom;
        end

    otherwise
        y = x;
end

end


function WindowStrainStatsTable = build_window_strain_stats_table_local(W)

vars = { ...
    'window_id', ...
    'time_start_s', ...
    'time_end_s', ...
    'time_center_s', ...
    'rpm_mean', ...
    'region_id', ...
    'EO_id', ...
    'fn_id_hz', ...
    'strain_fft_peak_hz', ...
    'strain_fft_peak_amp_microstrain', ...
    'strain_fft_peak_amp_eq_mm', ...
    'local_strain_fft_peak_hz', ...
    'local_strain_fft_peak_amp_microstrain', ...
    'local_strain_fft_peak_amp_eq_mm', ...
    'strain_sync_amp_at_step05_fn_microstrain', ...
    'strain_sync_amp_at_step05_fn_eq_mm', ...
    'freq_error_vs_strain_fft_peak_hz', ...
    'freq_error_vs_strain_fft_peak_percent', ...
    'freq_error_vs_local_strain_peak_hz', ...
    'strain_order_freq_hz', ...
    'freq_error_vs_order_line_hz', ...
    'strain_rms_microstrain', ...
    'strain_rms_eq_mm', ...
    'strain_half_peak_to_peak_microstrain', ...
    'strain_half_peak_to_peak_eq_mm', ...
    'validated'};

WindowStrainStatsTable = W(:, select_existing_vars_local(W, vars));

end


function Summary = build_validation_summary_local(W, Resp, bladeId, sensorTag)

valid = W.validated;

numWindows = height(W);
numValidated = nnz(valid);

medianFreqErrorHz = median(abs(W.freq_error_vs_strain_fft_peak_hz(valid)), 'omitnan');
meanFreqErrorHz = mean(abs(W.freq_error_vs_strain_fft_peak_hz(valid)), 'omitnan');

medianBestTARC = median(W.best_TARC(valid), 'omitnan');
meanBestTARC = mean(W.best_TARC(valid), 'omitnan');
maxBestTARC = max(W.best_TARC(valid), [], 'omitnan');

medianPhaseOptTARC = median(W.phase_opt_TARC(valid), 'omitnan');
meanPhaseOptTARC = mean(W.phase_opt_TARC(valid), 'omitnan');

if ~isempty(Resp) && height(Resp) >= 3
    x = Resp.Step05_A_norm;
    y = Resp.Strain_EqDispAmp_norm;

    ok = isfinite(x) & isfinite(y);

    if nnz(ok) >= 3
        C = corrcoef(x(ok), y(ok));
        corrAmpCurve = C(1,2);
    else
        corrAmpCurve = NaN;
    end
else
    corrAmpCurve = NaN;
end

if ~isempty(Resp) && height(Resp) >= 3
    e = Resp.EqDisp_AbsError_mm;
    ratio = Resp.EqDisp_Ratio_StrainOverStep05;

    medianAbsEqError = median(e, 'omitnan');
    meanAbsEqError = mean(e, 'omitnan');
    medianEqRatio = median(ratio, 'omitnan');
    meanEqRatio = mean(ratio, 'omitnan');

    ok = isfinite(Resp.Step05_A_mm) & isfinite(Resp.Strain_EqDispAmp_mm);

    if nnz(ok) >= 3
        C = corrcoef(Resp.Step05_A_mm(ok), Resp.Strain_EqDispAmp_mm(ok));
        corrEqDisp = C(1,2);
    else
        corrEqDisp = NaN;
    end
else
    medianAbsEqError = NaN;
    meanAbsEqError = NaN;
    medianEqRatio = NaN;
    meanEqRatio = NaN;
    corrEqDisp = NaN;
end

bestRegionId = most_frequent_region_local(W.region_id(valid));

representativeWindow = choose_representative_window_local(W, 'max_strain_sync_amp');

if isempty(representativeWindow)
    representativeWindowId = NaN;
else
    representativeWindowId = representativeWindow.window_id;
end

Summary = table( ...
    bladeId, ...
    string(sensorTag), ...
    numWindows, ...
    numValidated, ...
    medianFreqErrorHz, ...
    meanFreqErrorHz, ...
    medianBestTARC, ...
    meanBestTARC, ...
    maxBestTARC, ...
    medianPhaseOptTARC, ...
    meanPhaseOptTARC, ...
    corrAmpCurve, ...
    medianAbsEqError, ...
    meanAbsEqError, ...
    medianEqRatio, ...
    meanEqRatio, ...
    corrEqDisp, ...
    bestRegionId, ...
    representativeWindowId, ...
    'VariableNames', { ...
    'blade_id', ...
    'sensor_group', ...
    'num_windows', ...
    'num_validated_windows', ...
    'median_abs_freq_error_hz', ...
    'mean_abs_freq_error_hz', ...
    'median_best_TARC', ...
    'mean_best_TARC', ...
    'max_best_TARC', ...
    'median_phase_opt_TARC', ...
    'mean_phase_opt_TARC', ...
    'corr_response_amp_curve', ...
    'median_abs_eq_disp_error_mm', ...
    'mean_abs_eq_disp_error_mm', ...
    'median_eq_disp_ratio', ...
    'mean_eq_disp_ratio', ...
    'corr_eq_disp_amp_mm', ...
    'dominant_region_id', ...
    'representative_window_id'});

end


function RepresentativeWindow = choose_representative_window_local(W, rule)

RepresentativeWindow = table();

if isempty(W) || height(W) == 0 || ~ismember('validated', W.Properties.VariableNames)
    return;
end

valid = W.validated;

if ~any(valid)
    return;
end

Wv = W(valid, :);

switch lower(char(rule))
    case 'max_best_tarc'
        score = Wv.best_TARC;

    case 'min_freq_error'
        score = -abs(Wv.freq_error_vs_strain_fft_peak_hz);

    otherwise
        score = Wv.strain_sync_amp_at_step05_fn_microstrain;
end

score(~isfinite(score)) = -inf;

[~, im] = max(score);

if isempty(im) || ~isfinite(score(im))
    RepresentativeWindow = Wv(1,:);
else
    RepresentativeWindow = Wv(im,:);
end

end


function ChannelSensitivityTable = build_channel_sensitivity_local(Evidence, cfg, P, RepresentativeWindow)

ChannelSensitivityTable = table();

if isempty(RepresentativeWindow) || height(RepresentativeWindow) == 0
    return;
end

primaryFile = '';

if isfield(Evidence, 'StrainSource') && isfield(Evidence.StrainSource, 'strain_file')
    primaryFile = Evidence.StrainSource.strain_file;
elseif isfield(Evidence, 'StrainRaw') && isfield(Evidence.StrainRaw, 'file')
    primaryFile = Evidence.StrainRaw.file;
end

if isempty(primaryFile) || ~isfile(resolve_mat_file_local(primaryFile))
    return;
end

primaryFile = resolve_mat_file_local(primaryFile);
primaryDir = fileparts(primaryFile);

searchDirs = {primaryDir};

if isfield(cfg, 'strain_root') && ~isempty(cfg.strain_root)
    searchDirs{end+1} = cfg.strain_root;
end

files = find_candidate_files_local(searchDirs, P.channel_file_patterns);
files = unique(files, 'stable');

if isfinite(P.channel_max_files)
    files = files(1:min(numel(files), P.channel_max_files));
end

if isempty(files)
    return;
end

offset = get_optional_field_local(Evidence.StrainSource, 'strain_time_offset_total_sec', ...
    get_optional_field_local(Evidence.StrainRaw, 'timeOffsetSec', 0));

windows = RepresentativeWindow;
k_eq = get_eq_disp_scale_local(P);

rows = [];

for iFile = 1:numel(files)
    file = resolve_mat_file_local(files{iFile});

    try
        [tRaw, yRaw] = read_strain_file_generic_local(file, cfg);
    catch
        continue;
    end

    t = tRaw + offset;
    y = preprocess_strain_for_validation_local(yRaw, t);

    for iw = 1:height(windows)
        W = windows(iw,:);

        t0 = W.time_start_s;
        t1 = W.time_end_s;
        fn = W.fn_id_hz;

        if ~isfinite(t0) || ~isfinite(t1) || ~isfinite(fn)
            continue;
        end

        mask = t >= t0 & t <= t1 & isfinite(y);

        if nnz(mask) < P.min_tarc_points
            continue;
        end

        ts = t(mask);
        ys = y(mask);
        ys = detrend(ys, 'linear');
        ys = ys - mean(ys, 'omitnan');

        [ts2, ys2] = thin_signal_for_tarc_local(ts, ys, P.max_tarc_points);

        A = W.A_id_mm;
        if ~isfinite(A) || abs(A) < eps
            A = 1;
        end

        phi = W.phi_rad;
        if ~isfinite(phi)
            phi = 0;
        end

        rawTarc = calc_tarc_with_sine_local(ts2, ys2, fn, A, phi, 0);
        [bestTarc, bestShift] = calc_best_tarc_shift_local(ts2, ys2, fn, A, phi, P.time_shift_search_s);

        [syncAmp, syncPhase, syncFit] = sine_projection_local(ts, ys, fn);
        phaseOptTarc = calc_tarc_direct_local(ys, syncFit);

        if isfinite(k_eq)
            syncAmpEqMm = syncAmp * k_eq;
        else
            syncAmpEqMm = NaN;
        end

        rows = [rows; { ...
            W.window_id, ...
            W.time_center_s, ...
            string(get_strain_channel_label_local(file)), ...
            string(file), ...
            syncAmp, ...
            syncAmpEqMm, ...
            syncPhase, ...
            rawTarc, ...
            bestTarc, ...
            rawTarc, ...
            bestTarc, ...
            bestShift, ...
            phaseOptTarc, ...
            rms_local(ys), ...
            numel(ts2)}]; %#ok<AGROW>
    end
end

if isempty(rows)
    return;
end

ChannelSensitivityTable = cell2table(rows, 'VariableNames', { ...
    'window_id', ...
    'time_center_s', ...
    'strain_channel', ...
    'strain_file', ...
    'strain_sync_amp_microstrain', ...
    'strain_sync_amp_eq_mm', ...
    'strain_sync_phase_rad', ...
    'raw_TARC', ...
    'best_TARC', ...
    'raw_TRAC', ...
    'best_TRAC', ...
    'best_time_shift_s', ...
    'phase_opt_TARC', ...
    'strain_rms_microstrain', ...
    'point_count'});

end


function plot_step06b_main_validation_local(Evidence, W, Resp, Rep, P, figDir, visible, saveFigures, bladeId, sensorTag)

% Enhanced figure with paper-quality styling
fig = figure( ...
    'Name', sprintf('Step06B Main Validation B%d %s', bladeId, sensorTag), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [3, 3, 17, 11], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0, 0, 17, 11], ...
    'Visible', visible);

tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% (a) Step12-style selected-region strain FFT reference
ax1 = nexttile(tl);

tStrain = Evidence.StrainRaw.timeBtt(:);
yStrain = Evidence.StrainRaw.detrendedValue(:);
Ref = build_step12_style_strain_reference_local(tStrain, yStrain, Evidence.RegionTable, P);

hold(ax1, 'on');

if Ref.isValid
    hFft = plot(ax1, Ref.freqHz, Ref.ampMicrostrain, 'k-', 'LineWidth', 1.2);
    hPeak = xline(ax1, Ref.peakFreqHz, '-', 'Color', [0.8 0 0], 'LineWidth', 1.2);

    validFn = W.validated & isfinite(W.fn_id_hz);
    fnMed = median(W.fn_id_hz(validFn), 'omitnan');

    if isfinite(fnMed)
        hStep05 = xline(ax1, fnMed, '--', 'Color', [0 0 0.8], 'LineWidth', 1.2);
        legend(ax1, [hFft hPeak hStep05], ...
            {'Strain FFT', sprintf('Strain peak (%.1f Hz)', Ref.peakFreqHz), ...
            sprintf('Step05 median (%.1f Hz)', fnMed)}, ...
            'Location', 'northeast', 'FontSize', 7);
    else
        legend(ax1, [hFft hPeak], ...
            {'Strain FFT', sprintf('Strain peak (%.1f Hz)', Ref.peakFreqHz)}, ...
            'Location', 'northeast', 'FontSize', 7);
    end

    xlim(ax1, P.strain_peak_search_hz);
    inBand = Ref.freqHz >= P.strain_peak_search_hz(1) & Ref.freqHz <= P.strain_peak_search_hz(2);
    yMax = max(Ref.ampMicrostrain(inBand), [], 'omitnan');
    if isfinite(yMax) && yMax > 0
        ylim(ax1, [0, 1.18 * yMax]);
    end
    title(ax1, '(a) Strain FFT reference', 'FontWeight', 'normal');
else
    title(ax1, '(a) Strain FFT reference', 'FontWeight', 'normal');
end

xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'FFT amplitude (\mu蔚)');
apply_paper_axes_style_local(ax1);

%% (b) Frequency error
ax2 = nexttile(tl);

valid = W.validated & isfinite(W.freq_error_vs_strain_fft_peak_hz);

plot(ax2, W.window_id(valid), W.freq_error_vs_strain_fft_peak_hz(valid), 'o-', ...
    'LineWidth', 1.2, 'MarkerSize', 4, 'Color', [0 0 0.8]);
hold(ax2, 'on');
yline(ax2, 0, 'k--', 'LineWidth', 0.8);

xlabel(ax2, 'Step05 window');
ylabel(ax2, 'Frequency error (Hz)');
title(ax2, '(b) Frequency error', 'FontWeight', 'normal');
apply_paper_axes_style_local(ax2);

%% (c) Response amplitude curve
ax3 = nexttile(tl);

if ~isempty(Resp) && height(Resp) > 0
    x = Resp.window_id;

    plot(ax3, x, Resp.Step05_A_mm, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'Color', [0 0 0.8], 'DisplayName', 'Step05');
    hold(ax3, 'on');
    plot(ax3, x, Resp.Strain_EqDispAmp_mm, 's-', 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'Color', [0.8 0 0], 'DisplayName', 'Strain equiv.');
    legend(ax3, 'Location', 'best', 'FontSize', 7);
    xlabel(ax3, 'Step05 window');
end

ylabel(ax3, 'Amplitude (mm)');
title(ax3, '(c) Amplitude comparison', 'FontWeight', 'normal');
apply_paper_axes_style_local(ax3);

%% (d) Representative window waveform
ax4 = nexttile(tl);

if ~isempty(Rep) && height(Rep) > 0
    [tRep, yRep, ~] = get_representative_waveform_local(Evidence, Rep);

    if ~isempty(tRep)
        tr = 1000 * (tRep - min(tRep));

        k_eq = get_eq_disp_scale_local(P);

        if isfinite(k_eq)
            yRepEqMm = yRep * k_eq;
        else
            yRepEqMm = yRep;
        end

        [syncAmpEqMm, ~, ySyncEqMm] = sine_projection_local(tRep, yRepEqMm, Rep.fn_id_hz);

        if isfinite(syncAmpEqMm) && syncAmpEqMm > eps && ...
                isfinite(Rep.A_id_mm) && abs(Rep.A_id_mm) > eps
            yStep05Aligned = ySyncEqMm - mean(ySyncEqMm, 'omitnan');
            yStep05Aligned = abs(Rep.A_id_mm) * yStep05Aligned / syncAmpEqMm;
        else
            yStep05Aligned = nan(size(ySyncEqMm));
        end

        plot(ax4, tr, yRepEqMm, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, ...
            'DisplayName', 'Raw strain equiv.');
        hold(ax4, 'on');
        plot(ax4, tr, ySyncEqMm, 'Color', [0.8 0 0], 'LineWidth', 1.2, ...
            'DisplayName', 'Strain sync.');
        plot(ax4, tr, yStep05Aligned, '--', 'Color', [0 0 0.8], 'LineWidth', 1.2, ...
            'DisplayName', 'Step05 scaled to strain phase');

        legend(ax4, 'Location', 'best', 'FontSize', 7);

        title(ax4, sprintf('(d) Representative window (TARC=%.2f)', Rep.phase_opt_TARC), ...
            'FontWeight', 'normal');
    else
        title(ax4, '(d) Representative window', 'FontWeight', 'normal');
    end
else
    title(ax4, '(d) Representative window', 'FontWeight', 'normal');
end

xlabel(ax4, 'Time (ms)');
ylabel(ax4, 'Displacement (mm)');
apply_paper_axes_style_local(ax4);

if saveFigures
    save_figure_local(fig, figDir, sprintf('Step06B_Fig01_MainValidation_2x2_B%d_%s', bladeId, sensorTag));
end

end


function plot_step06b_tarc_trend_local(W, figDir, visible, saveFigures, bladeId, sensorTag)

fig = figure( ...
    'Name', sprintf('Step06B TARC Trend B%d %s', bladeId, sensorTag), ...
    'Color', 'w', ...
    'Position', [120, 80, 1200, 620], ...
    'Visible', visible);

tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);

valid = W.validated;
x = W.window_id(valid);

plot(ax1, x, W.raw_TARC(valid), 'o-', 'LineWidth', 1.0);
hold(ax1, 'on');
plot(ax1, x, W.best_TARC(valid), 's-', 'LineWidth', 1.0);
plot(ax1, x, W.phase_opt_TARC(valid), '.-', 'LineWidth', 1.0);

legend(ax1, {'raw TARC', 'best-shift TARC', 'phase-opt TARC'}, 'Location', 'best');
grid(ax1, 'on');
ylabel(ax1, 'TARC');
title(ax1, '(a) TARC / TRAC trend');

ax2 = nexttile(tl);

if ismember('strain_sync_amp_at_step05_fn_eq_mm', W.Properties.VariableNames)
    ySyncTrend = W.strain_sync_amp_at_step05_fn_eq_mm(valid);
    yLabelText = 'Sync strain equivalent displacement (mm)';
else
    ySyncTrend = W.strain_sync_amp_at_step05_fn_microstrain(valid);
    yLabelText = 'Sync strain amp. (\mue)';
end

plot(ax2, x, ySyncTrend, 'o-', 'LineWidth', 1.0);
grid(ax2, 'on');
xlabel(ax2, 'Step05 window');
ylabel(ax2, yLabelText);
title(ax2, '(b) Strain synchronous amplitude at Step05 fn');

if saveFigures
    save_figure_local(fig, figDir, sprintf('Step06B_Fig02_TARC_Trend_B%d_%s', bladeId, sensorTag));
end

end


function plot_step06b_channel_sensitivity_local(T, figDir, visible, saveFigures, bladeId, sensorTag)

if isempty(T) || height(T) == 0
    return;
end

fig = figure( ...
    'Name', sprintf('Step06B Channel Sensitivity B%d %s', bladeId, sensorTag), ...
    'Color', 'w', ...
    'Position', [140, 90, 1100, 720], ...
    'Visible', visible);

channels = string(T.strain_channel);
[uniqueChannels, ~, groupId] = unique(channels, 'stable');
phaseTarc = splitapply(@(x) median(x, 'omitnan'), T.phase_opt_TARC, groupId);

if ismember('strain_sync_amp_eq_mm', T.Properties.VariableNames)
    syncAmp = splitapply(@(x) median(x, 'omitnan'), T.strain_sync_amp_eq_mm, groupId);
    syncLabel = 'Sync equivalent displacement (mm)';
else
    syncAmp = splitapply(@(x) median(x, 'omitnan'), T.strain_sync_amp_microstrain, groupId);
    syncLabel = 'Sync strain amp. (\mue)';
end

[~, order] = sort(phaseTarc, 'descend');
uniqueChannels = uniqueChannels(order);
phaseTarc = phaseTarc(order);
syncAmp = syncAmp(order);

tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
bar(ax1, categorical(uniqueChannels), phaseTarc);
grid(ax1, 'on');
ylabel(ax1, 'Phase-opt TARC');
title(ax1, sprintf('Step06B channel sensitivity, B%d %s', bladeId, sensorTag));
xtickangle(ax1, 35);

ax2 = nexttile(tl);

bar(ax2, categorical(uniqueChannels), syncAmp);
ylabel(ax2, syncLabel);

grid(ax2, 'on');
xlabel(ax2, 'Strain channel');
xtickangle(ax2, 35);

if saveFigures
    save_figure_local(fig, figDir, sprintf('Step06B_Fig05_ChannelSensitivity_B%d_%s', bladeId, sensorTag));
end

end


function [tRep, yRep, yRecon] = get_representative_waveform_local(Evidence, Rep)

tRep = [];
yRep = [];
yRecon = [];

t = Evidence.StrainRaw.timeBtt(:);
y = Evidence.StrainRaw.detrendedValue(:);

mask = t >= Rep.time_start_s & t <= Rep.time_end_s & isfinite(t) & isfinite(y);

if nnz(mask) < 16
    return;
end

tRep = t(mask);
yRep = y(mask);

yRep = detrend(yRep, 'linear');
yRep = yRep - mean(yRep, 'omitnan');

fn = Rep.fn_id_hz;
A = Rep.A_id_mm;
phi = Rep.phi_rad;
shift = Rep.best_time_shift_s;

if ~isfinite(A) || abs(A) < eps
    A = 1;
end

if ~isfinite(phi)
    phi = 0;
end

if ~isfinite(shift)
    shift = 0;
end

yRecon = A * sin(2*pi*fn*(tRep + shift) + phi);

end


function [ff, aa] = single_fft_amp_local(t, y, freqRange, nfftMin)

ff = [];
aa = [];

t = t(:);
y = y(:);

valid = isfinite(t) & isfinite(y);

t = t(valid);
y = y(valid);

if numel(y) < 16
    return;
end

Fs = estimate_sample_rate_local(t);

if ~isfinite(Fs) || Fs <= 0
    return;
end

y = y - mean(y, 'omitnan');

nfft = max(nfftMin, 2^nextpow2(numel(y)));

win = hann_local(numel(y));

Y = fft(y .* win, nfft);

f = (0:floor(nfft/2)).' * Fs / nfft;

amp = 2 * abs(Y(1:numel(f))) / max(sum(win), eps);

mask = f >= freqRange(1) & f <= freqRange(2);

ff = f(mask);
aa = amp(mask);

end


function ratio = overlap_ratio_local(a, b)

ratio = 0;

if numel(a) < 2 || numel(b) < 2
    return;
end

if ~all(isfinite(a)) || ~all(isfinite(b))
    return;
end

overlap = max(0, min(a(2), b(2)) - max(a(1), b(1)));
duration = max(eps, a(2) - a(1));

ratio = overlap / duration;

end


function [amp, phase, fitSig] = sine_projection_local(t, y, f)

t = t(:);
y = y(:);

valid = isfinite(t) & isfinite(y);

t = t(valid);
y = y(valid);

if numel(y) < 3 || ~isfinite(f) || f <= 0
    amp = NaN;
    phase = NaN;
    fitSig = nan(size(y));
    return;
end

tr = t - mean(t, 'omitnan');

X = [sin(2*pi*f*tr), cos(2*pi*f*tr), ones(size(tr))];

coef = X \ y;

bSin = coef(1);
bCos = coef(2);

amp = hypot(bSin, bCos);
phase = atan2(bCos, bSin);

fitSig = X * coef;

end


function TARC = calc_tarc_with_sine_local(t, y, f, A, phi, shift)

t = t(:);
y = y(:);

valid = isfinite(t) & isfinite(y);

t = t(valid);
y = y(valid);

if numel(y) < 3 || ~isfinite(f) || f <= 0
    TARC = NaN;
    return;
end

if ~isfinite(A) || abs(A) < eps
    A = 1;
end

if ~isfinite(phi)
    phi = 0;
end

if ~isfinite(shift)
    shift = 0;
end

y = y - mean(y, 'omitnan');

r = A * sin(2*pi*f*(t + shift) + phi);
r = r - mean(r, 'omitnan');

TARC = calc_tarc_direct_local(y, r);

end


function [bestTarc, bestShift] = calc_best_tarc_shift_local(t, y, f, A, phi, shifts)

bestTarc = NaN;
bestShift = NaN;

if isempty(shifts)
    bestTarc = calc_tarc_with_sine_local(t, y, f, A, phi, 0);
    bestShift = 0;
    return;
end

vals = nan(numel(shifts), 1);

for i = 1:numel(shifts)
    vals(i) = calc_tarc_with_sine_local(t, y, f, A, phi, shifts(i));
end

[bestTarc, im] = max(vals, [], 'omitnan');

if isempty(im) || ~isfinite(bestTarc)
    bestShift = NaN;
else
    bestShift = shifts(im);
end

end


function [rho, trac] = calc_corr_trac_direct_local(y, r)
% Calculate correlation and TRAC (Time-domain Assurance Criterion)
%
% rho: signed correlation coefficient (-1 to 1)
% trac: squared normalized correlation (0 to 1), same as Modal Assurance Criterion

y = y(:);
r = r(:);

valid = isfinite(y) & isfinite(r);

y = y(valid);
r = r(valid);

if numel(y) < 3
    rho = NaN;
    trac = NaN;
    return;
end

y = y - mean(y, 'omitnan');
r = r - mean(r, 'omitnan');

den = sqrt((y' * y) * (r' * r));

if den <= 0 || ~isfinite(den)
    rho = NaN;
    trac = NaN;
else
    rho = (y' * r) / den;
    trac = rho^2;  % TRAC = squared correlation, range [0, 1]
end

end


function TARC = calc_tarc_direct_local(y, r)
% Legacy wrapper for backward compatibility
% Returns TRAC (squared correlation)

[~, TARC] = calc_corr_trac_direct_local(y, r);

end


function [t2, y2] = thin_signal_for_tarc_local(t, y, maxPts)

t = t(:);
y = y(:);

if numel(t) <= maxPts
    t2 = t;
    y2 = y;
    return;
end

idx = unique(round(linspace(1, numel(t), maxPts)));

t2 = t(idx);
y2 = y(idx);

end


function regionId = find_region_id_for_time_local(R, t)

regionId = NaN;

if isempty(R) || ~istable(R) || height(R) == 0 || ~isfinite(t)
    return;
end

% Robust column name handling with fallback for region ID
if ismember('time_start_s', R.Properties.VariableNames)
    t0 = R.time_start_s;
    t1 = R.time_end_s;

    % Handle region_id column with fallback
    if ismember('region_id', R.Properties.VariableNames)
        rid = R.region_id;
    elseif ismember('regionId', R.Properties.VariableNames)
        rid = R.regionId;
    else
        % Fallback: use row indices as region IDs
        rid = (1:height(R)).';
    end

elseif ismember('bttStartSec', R.Properties.VariableNames)
    t0 = R.bttStartSec;
    t1 = R.bttEndSec;

    % Handle region_id column with fallback
    if ismember('regionId', R.Properties.VariableNames)
        rid = R.regionId;
    elseif ismember('region_id', R.Properties.VariableNames)
        rid = R.region_id;
    else
        % Fallback: use row indices as region IDs
        rid = (1:height(R)).';
    end
else
    return;
end

hit = find(t >= t0 & t <= t1, 1, 'first');

if ~isempty(hit)
    regionId = rid(hit);
end

end


function rpmVal = interp_rpm_local(Evidence, t)

rpmVal = NaN;

if ~isfield(Evidence, 'RPM')
    return;
end

rpm = Evidence.RPM;

if ~isfield(rpm, 'timeSec')
    if isfield(rpm, 'time_s')
        tr = rpm.time_s;
    else
        return;
    end
else
    tr = rpm.timeSec;
end

if ~isfield(rpm, 'value')
    return;
end

vr = rpm.value;

valid = isfinite(tr) & isfinite(vr);

tr = tr(valid);
vr = vr(valid);

if numel(tr) < 2
    return;
end

rpmVal = interp1(tr, vr, t, 'linear', NaN);

end


function [t, y] = read_strain_file_generic_local(file, cfg)

file = resolve_mat_file_local(file);

D = load(file);

if isfield(D, 'Datas') && isnumeric(D.Datas) && size(D.Datas,2) >= 2
    M = D.Datas;
else
    M = [];
    fn = fieldnames(D);
    bestScore = -inf;

    for i = 1:numel(fn)
        v = D.(fn{i});

        if isnumeric(v) && ismatrix(v) && size(v,1) >= 16
            score = size(v,1) + 1000 * double(size(v,2) >= 2);

            if score > bestScore
                M = v;
                bestScore = score;
            end
        end
    end

    if isempty(M)
        error('No numeric matrix found in strain file: %s', file);
    end
end

M = double(M);

if size(M,2) < 2
    M = [(0:size(M,1)-1).', M(:)];
end

x = M(:,1);
y = M(:,2);

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);

[x, order] = sort(x);
y = y(order);

[x, ia] = unique(x, 'stable');
y = y(ia);

dx = median(diff(x), 'omitnan');

if isfinite(dx) && dx > 0 && dx < 0.1 && max(x) < 1e6
    t = x;
else
    t = x ./ cfg.sample_rate_hz;
end

t = t(:);
y = y(:);

end


function y = preprocess_strain_for_validation_local(yRaw, t)

yRaw = yRaw(:);
t = t(:);

valid = isfinite(yRaw);

if any(valid)
    yRaw(~valid) = median(yRaw(valid), 'omitnan');
else
    yRaw(:) = 0;
end

if numel(yRaw) >= 3
    y = detrend(yRaw, 'linear');
else
    y = yRaw - mean(yRaw, 'omitnan');
end

y = y - median(y, 'omitnan');

end


function files = find_candidate_files_local(dirs0, patterns)

files = {};

for id = 1:numel(dirs0)
    root = dirs0{id};

    if ~isfolder(root)
        continue;
    end

    for ip = 1:numel(patterns)
        dd = dir(fullfile(root, '**', patterns{ip}));

        for k = 1:numel(dd)
            files{end+1} = fullfile(dd(k).folder, dd(k).name); %#ok<AGROW>
        end
    end
end

files = unique(files, 'stable');

end


function [t, y] = normalize_xy_local(t, y)

t = t(:);
y = y(:);

valid = isfinite(t) & isfinite(y);

t = t(valid);
y = y(valid);

end


function y = normalize_waveform_local(y)

y = y(:);

y = y - mean(y, 'omitnan');

s = max(abs(y), [], 'omitnan');

if ~isfinite(s) || s <= 0
    s = 1;
end

y = y ./ s;

end


function shade_regions_local(ax, R)

if isempty(R) || ~istable(R) || height(R) == 0
    return;
end

if ismember('time_start_s', R.Properties.VariableNames)
    t0 = R.time_start_s;
    t1 = R.time_end_s;
elseif ismember('bttStartSec', R.Properties.VariableNames)
    t0 = R.bttStartSec;
    t1 = R.bttEndSec;
else
    return;
end

yl = ylim(ax);
hold(ax, 'on');

for i = 1:numel(t0)
    patch(ax, [t0(i), t1(i), t1(i), t0(i)], [yl(1), yl(1), yl(2), yl(2)], ...
        [0.85 0.85 0.85], ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

ylim(ax, yl);

end


function vars = select_existing_vars_local(T, vars0)

vars = {};

for i = 1:numel(vars0)
    if ismember(vars0{i}, T.Properties.VariableNames)
        vars{end+1} = vars0{i}; %#ok<AGROW>
    end
end

end


function v = get_numeric_field_multi_local(S1, S2, names, defaultValue)

v = get_numeric_field_local(S2, names, NaN);

if ~isfinite(v)
    v = get_numeric_field_local(S1, names, defaultValue);
end

end


function s = get_char_field_multi_local(S1, S2, names, defaultValue)

s = get_char_field_local(S2, names, '');

if isempty(s)
    s = get_char_field_local(S1, names, defaultValue);
end

end


function v = get_numeric_field_local(S, names, defaultValue)

v = defaultValue;

if ~isstruct(S)
    return;
end

for i = 1:numel(names)
    name = names{i};

    if isfield(S, name)
        x = S.(name);

        if isnumeric(x) && isscalar(x)
            v = double(x);
            return;
        end

        if islogical(x) && isscalar(x)
            v = double(x);
            return;
        end
    end
end

end


function v = get_numeric_vector_field_local(S, names)

v = [];

if ~isstruct(S)
    return;
end

for i = 1:numel(names)
    name = names{i};

    if isfield(S, name)
        x = S.(name);

        if isnumeric(x)
            v = double(x(:));
            return;
        end
    end
end

end


function S2 = get_struct_field_local(S, names, defaultValue)

S2 = defaultValue;

if ~isstruct(S)
    return;
end

for i = 1:numel(names)
    name = names{i};

    if isfield(S, name) && isstruct(S.(name))
        S2 = S.(name);
        return;
    end
end

end


function s = get_char_field_local(S, names, defaultValue)

s = defaultValue;

if ~isstruct(S)
    return;
end

for i = 1:numel(names)
    name = names{i};

    if isfield(S, name)
        x = S.(name);

        if ischar(x)
            s = x;
            return;
        elseif isstring(x) && isscalar(x)
            s = char(x);
            return;
        elseif isnumeric(x) && isscalar(x)
            s = num2str(x);
            return;
        end
    end
end

end


function x = get_table_numeric_local(T, candidates, defaultValue)

x = defaultValue;

for i = 1:numel(candidates)
    name = candidates{i};

    if ismember(name, T.Properties.VariableNames)
        v = T.(name);

        if isnumeric(v) || islogical(v)
            x = double(v(:));
            return;
        end
    end
end

end


function name = get_file_name_local(file)

[~, name, ext] = fileparts(file);
name = [name ext];

end


function label = get_strain_channel_label_local(file)

name = get_file_name_local(file);
tok = regexp(name, '^(AI\d+-\d+)', 'tokens', 'once');

if ~isempty(tok)
    label = tok{1};
else
    [~, label] = fileparts(file);
end

end


function file = resolve_mat_file_local(file0)

file = char(file0);

if isfile(file)
    return;
end

[folder, name, ext] = fileparts(file);

if isempty(ext)
    file2 = fullfile(folder, [name '.mat']);

    if isfile(file2)
        file = file2;
        return;
    end
end

end


function regionId = most_frequent_region_local(x)

x = x(:);
x = x(isfinite(x));

if isempty(x)
    regionId = NaN;
    return;
end

u = unique(x);
cnt = zeros(size(u));

for i = 1:numel(u)
    cnt(i) = nnz(x == u(i));
end

[~, im] = max(cnt);
regionId = u(im);

end


function Fs = estimate_sample_rate_local(t)

t = t(:);
t = t(isfinite(t));

if numel(t) < 2
    Fs = NaN;
    return;
end

dt = median(diff(t), 'omitnan');

if ~isfinite(dt) || dt <= 0
    Fs = NaN;
else
    Fs = 1 / dt;
end

end


function r = rms_local(x)

x = x(:);
x = x(isfinite(x));

if isempty(x)
    r = NaN;
else
    r = sqrt(mean(x.^2));
end

end


function w = hann_local(n)

if exist('hann', 'file') == 2
    w = hann(n);
else
    w = 0.5 - 0.5 * cos(2*pi*(0:n-1).'/(n-1));
end

end


function cmap = high_contrast_colormap_local(n)

if nargin < 1
    n = 256;
end

try
    cmap = turbo(n);
catch
    cmap = jet(n);
end

end


function save_figure_local(fig, figDir, name)

ensure_dir_local(figDir);

figPath = fullfile(figDir, [name '.fig']);
pngPath = fullfile(figDir, [name '.png']);
pdfPath = fullfile(figDir, [name '.pdf']);

try
    oldDir = pwd;
    cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(figDir);

    savefig(fig, [name '.fig']);

catch ME
    warning('savefig failed for:\n  %s\nReason: %s\nPNG export will continue.', ...
        figPath, ME.message);
end

% Export high-resolution PNG
try
    exportgraphics(fig, pngPath, 'Resolution', 600);
catch ME1
    try
        print(fig, pngPath, '-dpng', '-r600');
    catch ME2
        warning('PNG export failed for:\n  %s\nexportgraphics: %s\nprint: %s', ...
            pngPath, ME1.message, ME2.message);
    end
end

% Export vector PDF for publication
try
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
catch ME
    warning('PDF export failed for:\n  %s\nReason: %s', pdfPath, ME.message);
end

end


function apply_paper_axes_style_local(ax)
% Apply consistent paper-quality styling to axes

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 8.5, ...
    'LineWidth', 0.8, ...
    'TickDir', 'in', ...
    'Box', 'on', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

grid(ax, 'off');

end


function ensure_dir_local(d)

if exist(d, 'dir') ~= 7
    mkdir(d);
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


function ok = validate_evidence_schema_local(Evidence)
% Validate that Evidence structure has required fields and consistent dimensions

ok = false;

if ~isstruct(Evidence)
    warning('Evidence is not a struct.');
    return;
end

% Check StrainRaw fields
if ~isfield(Evidence, 'StrainRaw')
    warning('Evidence missing field: StrainRaw');
    return;
end

if ~isfield(Evidence.StrainRaw, 'timeBtt') || ~isfield(Evidence.StrainRaw, 'detrendedValue')
    warning('Evidence.StrainRaw missing timeBtt or detrendedValue');
    return;
end

if numel(Evidence.StrainRaw.timeBtt) ~= numel(Evidence.StrainRaw.detrendedValue)
    warning('Evidence.StrainRaw: timeBtt and detrendedValue length mismatch');
    return;
end

if isempty(Evidence.StrainRaw.timeBtt)
    warning('Evidence.StrainRaw.timeBtt is empty');
    return;
end

% Check RegionTable
if ~isfield(Evidence, 'RegionTable')
    warning('Evidence missing field: RegionTable');
    return;
end

if ~istable(Evidence.RegionTable)
    warning('Evidence.RegionTable is not a table');
    return;
end

% Check RPM fields (optional but recommended)
if isfield(Evidence, 'RPM')
    rpm = Evidence.RPM;
    hasTime = isfield(rpm, 'timeSec') || isfield(rpm, 'time_s');
    hasValue = isfield(rpm, 'value');

    if ~hasTime || ~hasValue
        warning('Evidence.RPM exists but missing timeSec/time_s or value fields');
    end
end

ok = true;

end
