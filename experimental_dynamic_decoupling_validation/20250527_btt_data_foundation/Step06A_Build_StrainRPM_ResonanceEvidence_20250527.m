%% Step06A_Build_StrainRPM_ResonanceEvidence_20250527.m
% Build method-neutral strain/RPM/resonance evidence for the 20250527 BTT route.
%
% This script follows the old Step04_Detect_Resonance_By_BTT_STE_20250527 logic:
%   1) use an explicitly selected strain gauge channel;
%   2) align strain time to BTT time;
%   3) crop strain to valid BTT/RPM time range;
%   4) use raw segmented strain FFT surface amplitude to detect resonance regions;
%   5) save old-style and new-style evidence fields for Step06B.
%
% This script does NOT read Step05 and does NOT validate Step05.

clc; clear; close all;

%% 0. Config
if exist('BTTDataConfig_20250527', 'file') == 2
    cfg = BTTDataConfig_20250527();
elseif exist('Get_20250527_BTT_Config', 'file') == 2
    cfg = Get_20250527_BTT_Config();
else
    error('Cannot find BTTDataConfig_20250527.m or Get_20250527_BTT_Config.m.');
end

script_dir = fileparts(mfilename('fullpath'));
Step06A_Revision = '2026-07-02-step04-legacy-reference-v4-AI1-04';

route_parent_dir = fileparts(script_dir);

if ~isfield(cfg, 'step04_reference_legacy_output_root') || ...
        isempty(cfg.step04_reference_legacy_output_root)
    cfg.step04_reference_legacy_output_root = fullfile(route_parent_dir, ...
        '20250527_low_speed_gap_prior_decoupling', 'legacy', 'output');
end

if ~isfield(cfg, 'dataset') || isempty(cfg.dataset)
    cfg.dataset = '20250527';
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

if ~isfield(cfg, 'sample_rate_hz') || isempty(cfg.sample_rate_hz)
    if isfield(cfg, 'pinlv') && ~isempty(cfg.pinlv)
        cfg.sample_rate_hz = cfg.pinlv;
    else
        cfg.sample_rate_hz = 5e6;
    end
end

if ~isfield(cfg, 'strain_root') || isempty(cfg.strain_root)
    cfg.strain_root = 'E:\试验数据\20250527\应变片数据20250527';
end

if ~isfield(cfg, 'dynamic_cases') || isempty(cfg.dynamic_cases)
    if isfield(cfg, 'default_dynamic_cases') && ~isempty(cfg.default_dynamic_cases)
        cfg.dynamic_cases = cfg.default_dynamic_cases;
    else
        cfg.dynamic_cases = {'20250526_2500-3500_t400'};
    end
end

if ischar(cfg.dynamic_cases) || isstring(cfg.dynamic_cases)
    cfg.dynamic_cases = cellstr(cfg.dynamic_cases);
end

if ~isfield(cfg, 'step03_output_dir') || isempty(cfg.step03_output_dir)
    cfg.step03_output_dir = fullfile(cfg.output_root, 'step03_dynamic_observation_bundle');
end

if ~isfield(cfg, 'step3_default_btt_sensor_ids') || isempty(cfg.step3_default_btt_sensor_ids)
    if isfield(cfg, 'capacitance_ids') && ~isempty(cfg.capacitance_ids)
        cfg.step3_default_btt_sensor_ids = cfg.capacitance_ids;
    elseif isfield(cfg, 'sensor_ids') && ~isempty(cfg.sensor_ids)
        cfg.step3_default_btt_sensor_ids = cfg.sensor_ids;
    else
        cfg.step3_default_btt_sensor_ids = [1 3 6];
    end
end

if ~isfield(cfg, 'step3_default_strain_channel') || isempty(cfg.step3_default_strain_channel)
    cfg.step3_default_strain_channel = 'AI1-04';
end

if ~isfield(cfg, 'step3_default_alignment_offset_sec') || isempty(cfg.step3_default_alignment_offset_sec)
    cfg.step3_default_alignment_offset_sec = 0;
end

if ~isfield(cfg, 'step3_align_order_candidates') || isempty(cfg.step3_align_order_candidates)
    cfg.step3_align_order_candidates = 8:20;
end

%% 1. Run control
target_cases = cfg.dynamic_cases;
show_plots = true;
force_rebuild = true;

%% 2. Output folders
step06a_output_root = fullfile(cfg.output_root, 'step06a_strain_rpm_resonance_evidence');
step06a_figure_root = fullfile(cfg.figure_root, 'step06a_strain_rpm_resonance_evidence');

%% 3. Strain/RPM/BTT settings
targetBladeId = 1;

% Only used for BTT displacement display, not for resonance detection.
btt_sensor_ids_for_plot = cfg.step3_default_btt_sensor_ids;

baselineTimeSec = 3.0;

% Explicit strain file override.
% This is the strain gauge channel used by Step06A.
strain_file_override = 'E:\试验数据\20250527\应变片数据20250527\20250527_2500_3500_time400\AI1-04_20250527171419.mat';
strain_file_override = resolve_existing_mat_or_empty_local(strain_file_override, ...
    'configured strain_file_override');

% Official strain root and default channel from old config.
strain_root = cfg.strain_root;
strain_primary_keyword = cfg.step3_default_strain_channel;

% Fallback search patterns if override is empty.
strain_file_patterns = { ...
    [strain_primary_keyword, '*.mat'], ...
    [strain_primary_keyword, '_*.mat'], ...
    'AI1-*.mat', ...
    '*strain*.mat', ...
    '*Strain*.mat'};

strain_primary_index_if_no_keyword = 1;

% If the strain file Datas(:,2) is already microstrain, keep scale = 1.
strain_value_scale_to_microstrain = 1;
strain_value_offset_microstrain = 0;

% If legacy Step3_Result exists, its strain_time_offset_total_sec is used.
% Otherwise this manual offset is used.
manual_strain_time_offset_total_sec = cfg.step3_default_alignment_offset_sec;

% Old Step04 style.
detrendMethod = 'linear';       % 'linear', 'movingmedian', 'none'

%% 4. Resonance detection settings
orders = cfg.step3_align_order_candidates(:).';

% Old Step04 core rule:
% raw segmented strain FFT surface amplitude >= threshold.
resonanceFftAmpThresholdMicrostrain = 100.0;

% Only reported as check, not primary detection rule.
resonanceTimeAmpThresholdMicrostrain = 100.0;
resonanceLowThresholdRatio = 0.55;
resonanceCombineMode = "raw_fft_surface_only";

% Region post-processing.
regionMergeGapSec = 0.00;
minRegionDurationSec = 0.00;
maxRegionWidthSec = inf;
preferredRegionTimeSec = [];

% Old-style equal segmented FFT.
rawFftNumSegments = 200;
rawFftPeakBandHz = [20, 800];

% Display frequency band.
dashboardFreqBandHz = [20, 800];

make3DFigure = true;

%% 5. Valid-time and plot settings
crop_strain_to_valid_time_range = true;
rpm_min_valid_rpm = 100;
rpm_plot_y_margin_ratio = 0.08;
max_plot_points = 60000;

save_figures = cfg.save_figures;
fig_visible = ternary_local(show_plots, 'on', 'off');

%% 6. Pack settings
Step06A_Setting = struct();

Step06A_Setting.target_cases = target_cases;
Step06A_Setting.targetBladeId = targetBladeId;
Step06A_Setting.btt_sensor_ids_for_plot = btt_sensor_ids_for_plot;
Step06A_Setting.baselineTimeSec = baselineTimeSec;

Step06A_Setting.strain_root = strain_root;
Step06A_Setting.strain_file_override = strain_file_override;
Step06A_Setting.strain_file_patterns = strain_file_patterns;
Step06A_Setting.strain_primary_keyword = strain_primary_keyword;
Step06A_Setting.strain_primary_index_if_no_keyword = strain_primary_index_if_no_keyword;
Step06A_Setting.strain_value_scale_to_microstrain = strain_value_scale_to_microstrain;
Step06A_Setting.strain_value_offset_microstrain = strain_value_offset_microstrain;
Step06A_Setting.manual_strain_time_offset_total_sec = manual_strain_time_offset_total_sec;
Step06A_Setting.detrendMethod = detrendMethod;

Step06A_Setting.orders = orders;
Step06A_Setting.rawFftNumSegments = rawFftNumSegments;
Step06A_Setting.rawFftPeakBandHz = rawFftPeakBandHz;
Step06A_Setting.dashboardFreqBandHz = dashboardFreqBandHz;

Step06A_Setting.resonanceFftAmpThresholdMicrostrain = resonanceFftAmpThresholdMicrostrain;
Step06A_Setting.resonanceTimeAmpThresholdMicrostrain = resonanceTimeAmpThresholdMicrostrain;
Step06A_Setting.resonanceLowThresholdRatio = resonanceLowThresholdRatio;
Step06A_Setting.resonanceCombineMode = resonanceCombineMode;
Step06A_Setting.regionMergeGapSec = regionMergeGapSec;
Step06A_Setting.minRegionDurationSec = minRegionDurationSec;
Step06A_Setting.maxRegionWidthSec = maxRegionWidthSec;
Step06A_Setting.preferredRegionTimeSec = preferredRegionTimeSec;

Step06A_Setting.crop_strain_to_valid_time_range = crop_strain_to_valid_time_range;
Step06A_Setting.rpm_min_valid_rpm = rpm_min_valid_rpm;
Step06A_Setting.rpm_plot_y_margin_ratio = rpm_plot_y_margin_ratio;
Step06A_Setting.max_plot_points = max_plot_points;
Step06A_Setting.make3DFigure = make3DFigure;
Step06A_Setting.save_figures = save_figures;
Step06A_Setting.revision = Step06A_Revision;
Step06A_Setting.script_fullpath = mfilename('fullpath');

fprintf('\n=== Step06A: strain/RPM/resonance evidence, old-Step04-style ===\n');
fprintf('Step06A revision: %s\n', Step06A_Revision);
fprintf('Running script file:\n  %s\n', mfilename('fullpath'));
fprintf('Dataset: %s\n', cfg.dataset);
fprintf('Target cases: %s\n', strjoin(cellstr(target_cases), ', '));
fprintf('Strain override:\n  %s\n', strain_file_override);
fprintf('Rule: raw segmented strain FFT peak amplitude >= %.1f microstrain in %.1f-%.1f Hz.\n', ...
    resonanceFftAmpThresholdMicrostrain, rawFftPeakBandHz(1), rawFftPeakBandHz(2));

%% 7. Main loop
for iCase = 1:numel(target_cases)
    caseName = char(target_cases{iCase});

    fprintf('\n--- Step06A case: %s ---\n', caseName);

    outDir = fullfile(step06a_output_root, caseName);
    figDir = fullfile(step06a_figure_root, caseName);

    ensure_dir_local(outDir);
    ensure_dir_local(figDir);

    evidenceMat = fullfile(outDir, 'Step06A_StrainRPM_ResonanceEvidence_20250527.mat');

    if isfile(evidenceMat) && ~force_rebuild
        fprintf('Existing Step06A evidence found. Skip because force_rebuild=false:\n  %s\n', evidenceMat);
        continue;
    end

    %% 7.1 Load Step03 bundle and optional legacy results
    [bundle, bundleFile] = load_step03_bundle_local(cfg, caseName);
    rpm = extract_rpm_from_bundle_local(bundle);

    legacy = load_legacy_alignment_local(cfg, caseName);

    if legacy.hasAlignResult
        alignResult = legacy.alignResult;

        rpmLegacy = extract_rpm_from_legacy_align_local(alignResult);
        if ~isempty(rpmLegacy.timeSec)
            rpm = rpmLegacy;
        end

        validTimeRange = compute_valid_time_range_from_legacy_local( ...
            alignResult, rpm, rpm_min_valid_rpm);
    else
        alignResult = struct();
        validTimeRange = compute_valid_time_range_from_rpm_local(rpm, rpm_min_valid_rpm);
    end

    fprintf('Valid BTT/RPM time range for evidence: %.6f to %.6f s.\n', ...
        validTimeRange(1), validTimeRange(2));

    %% 7.2 Load strain
    [strainTimeBtt, strainValue, strainSource, strainRaw] = load_strain_evidence_local( ...
        cfg, caseName, Step06A_Setting, legacy, validTimeRange);

    if isempty(strainTimeBtt)
        error('No valid strain evidence found for case %s.', caseName);
    end

    fprintf('\nStep06A strain source:\n');
    fprintf('  mode   : %s\n', strainSource.source_mode);
    fprintf('  file   : %s\n', strainSource.strain_file);
    fprintf('  offset : %.9g s\n', strainSource.strain_time_offset_total_sec);

    %% 7.3 Crop strain to valid BTT/RPM time range
    if crop_strain_to_valid_time_range
        keep = strainRaw.timeBtt >= validTimeRange(1) & ...
               strainRaw.timeBtt <= validTimeRange(2) & ...
               isfinite(strainRaw.detrendedValue);

        strainRaw.timeBtt = strainRaw.timeBtt(keep);
        strainRaw.rawValue = strainRaw.rawValue(keep);
        strainRaw.detrendedValue = strainRaw.detrendedValue(keep);

        strainTimeBtt = strainRaw.timeBtt;
        strainValue = strainRaw.detrendedValue;

        fprintf('Cropped strain to valid BTT/RPM time range: %d points remain.\n', ...
            numel(strainTimeBtt));
    end

    if numel(strainTimeBtt) < 32
        error('Too few strain points after valid-time cropping.');
    end

    %% 7.4 Load target-blade BTT displacement for display only
    [bttTime, bttVib, bttSource] = load_btt_vibration_for_display_local( ...
        cfg, caseName, legacy, bundle, btt_sensor_ids_for_plot, ...
        targetBladeId, baselineTimeSec);

    %% 7.5 Raw segmented FFT and resonance regions
    rawFft = raw_segmented_fft_local( ...
        strainRaw, validTimeRange, rawFftNumSegments, rawFftPeakBandHz);

    resonanceParams = struct();
    resonanceParams.freqBandHz = dashboardFreqBandHz;
    resonanceParams.peakBandHz = rawFftPeakBandHz;
    resonanceParams.fftAmpThresholdMicrostrain = resonanceFftAmpThresholdMicrostrain;
    resonanceParams.timeAmpThresholdMicrostrain = resonanceTimeAmpThresholdMicrostrain;
    resonanceParams.lowThresholdRatio = resonanceLowThresholdRatio;
    resonanceParams.combineMode = resonanceCombineMode;
    resonanceParams.regionMergeGapSec = regionMergeGapSec;
    resonanceParams.minRegionDurationSec = minRegionDurationSec;
    resonanceParams.maxRegionWidthSec = maxRegionWidthSec;
    resonanceParams.preferredRegionTimeSec = preferredRegionTimeSec;

    [scoreTable, peakTable, regionTable, localFftTable, localFftResult] = ...
        detect_resonance_regions_from_raw_fft_local(rawFft, rpm, orders, resonanceParams);

    %% 7.6 Build Evidence
    Evidence = struct();

    Evidence.Dataset = cfg.dataset;
    Evidence.CaseName = caseName;
    Evidence.CreatedBy = mfilename;
    Evidence.CreatedOn = datestr(now, 31);

    Evidence.BundleFile = bundleFile;
    Evidence.LegacyAlignment = legacy;
    Evidence.AlignResult = alignResult;
    Evidence.ValidTimeRange = validTimeRange;

    Evidence.StrainTimeBtt = strainTimeBtt;
    Evidence.StrainValue = strainValue;
    Evidence.StrainSource = strainSource;
    Evidence.StrainRaw = strainRaw;

    Evidence.BTT_Source = bttSource;
    Evidence.BTT_TargetBlade.time_s = bttTime;
    Evidence.BTT_TargetBlade.vibration_mm = bttVib;
    Evidence.BTT_TargetBlade.blade_id = targetBladeId;
    Evidence.BTT_TargetBlade.sensor_ids = btt_sensor_ids_for_plot;

    Evidence.RPM = rpm;

    Evidence.rawFft = rawFft;
    Evidence.RawFft = rawFft;

    Evidence.scoreTable = scoreTable;
    Evidence.peakTable = peakTable;
    Evidence.regionTable = regionTable;
    Evidence.localFftTable = localFftTable;
    Evidence.localFftResult = localFftResult;

    Evidence.ScoreTable = scoreTable;
    Evidence.PeakTable = peakTable;
    Evidence.RegionTable = regionTable;
    Evidence.LocalFftTable = localFftTable;
    Evidence.LocalFftResult = localFftResult;

    Evidence.resonanceParams = resonanceParams;
    Evidence.ResonanceParams = resonanceParams;
    Evidence.Step06A_Setting = Step06A_Setting;

    %% 7.7 Save
    save(evidenceMat, ...
        'Evidence', ...
        'scoreTable', ...
        'peakTable', ...
        'regionTable', ...
        'localFftTable', ...
        'localFftResult', ...
        'rawFft', ...
        'strainSource', ...
        'resonanceParams', ...
        'alignResult', ...
        'caseName', ...
        'Step06A_Setting', ...
        '-v7.3');

    writetable(regionTable, fullfile(outDir, 'Step06A_RegionTable_20250527.csv'));
    writetable(scoreTable, fullfile(outDir, 'Step06A_ScoreTable_20250527.csv'));
    writetable(localFftTable, fullfile(outDir, 'Step06A_LocalStrainFFT_20250527.csv'));

    fprintf('\nRegion summary:\n');
    disp(regionTable);

    fprintf('Saved Step06A evidence:\n  %s\n', evidenceMat);

    %% 7.8 Plot
    if show_plots || save_figures
        plot_dashboard_local( ...
            strainTimeBtt, ...
            strainValue, ...
            bttTime, ...
            bttVib, ...
            targetBladeId, ...
            rpm, ...
            rawFft, ...
            peakTable, ...
            regionTable, ...
            validTimeRange, ...
            dashboardFreqBandHz, ...
            fig_visible, ...
            figDir, ...
            save_figures, ...
            Step06A_Setting);

        if make3DFigure
            plot_3d_auxiliary_local( ...
                rawFft, ...
                validTimeRange, ...
                dashboardFreqBandHz, ...
                regionTable, ...
                fig_visible, ...
                figDir, ...
                save_figures);
        end
    end
end

fprintf('\n=== Step06A finished ===\n');

%% ========================================================================
% Local functions
% ========================================================================

function [bundle, bundleFile] = load_step03_bundle_local(cfg, caseName)

candidateFiles = {};

if isfield(cfg, 'step03_output_dir') && ~isempty(cfg.step03_output_dir)
    candidateFiles{end+1} = fullfile(cfg.step03_output_dir, caseName, 'BTT_Observation_Bundle_20250527.mat');
end

candidateFiles{end+1} = fullfile(cfg.output_root, 'step03_dynamic_observation_bundle', caseName, 'BTT_Observation_Bundle_20250527.mat');
candidateFiles{end+1} = fullfile(cfg.output_root, caseName, 'BTT_Observation_Bundle_20250527.mat');

candidateFiles = unique(candidateFiles, 'stable');

bundleFile = '';

for i = 1:numel(candidateFiles)
    if isfile(candidateFiles{i})
        bundleFile = candidateFiles{i};
        break;
    end
end

if isempty(bundleFile)
    dd = dir(fullfile(cfg.output_root, '**', 'BTT_Observation_Bundle_20250527.mat'));
    for k = 1:numel(dd)
        fullName = fullfile(dd(k).folder, dd(k).name);
        if contains(fullName, caseName)
            bundleFile = fullName;
            break;
        end
    end
end

if isempty(bundleFile)
    error('Missing Step03 BTT_Observation_Bundle_20250527.mat for case %s.', caseName);
end

D = load(bundleFile);

if isfield(D, 'bundle')
    bundle = D.bundle;
elseif isfield(D, 'BTT_Observation_Bundle')
    bundle = D.BTT_Observation_Bundle;
else
    fn = fieldnames(D);
    bundle = D.(fn{1});
end

fprintf('Loaded Step03 bundle:\n  %s\n', bundleFile);

end


function rpm = extract_rpm_from_bundle_local(bundle)

rpm = struct();
rpm.timeSec = [];
rpm.time_s = [];
rpm.value = [];
rpm.source = 'Step03 bundle';

if isfield(bundle, 'Omega_Time_s') && isfield(bundle, 'Omega_RPM')
    rpm.timeSec = bundle.Omega_Time_s(:);
    rpm.value = bundle.Omega_RPM(:);
    rpm.source = get_optional_field_local(bundle, 'Speed_Source', 'Step03_Omega_RPM');

elseif isfield(bundle, 'Omega_Diagnostic_Time_s') && isfield(bundle, 'Omega_Diagnostic_RPM')
    rpm.timeSec = bundle.Omega_Diagnostic_Time_s(:);
    rpm.value = bundle.Omega_Diagnostic_RPM(:);
    rpm.source = 'Step03_Omega_Diagnostic_RPM';

elseif isfield(bundle, 'Observation_Table') && istable(bundle.Observation_Table)
    T = bundle.Observation_Table;

    timeVar = find_table_var_local(T, {'rpm_time_s', 'time_s', 'arrival_time_s', 'btt_time_s'});
    rpmVar = find_table_var_local(T, {'rpm', 'RPM', 'Omega_RPM'});

    if ~isempty(timeVar) && ~isempty(rpmVar)
        rpm.timeSec = T.(timeVar)(:);
        rpm.value = T.(rpmVar)(:);
        rpm.source = 'Step03_Observation_Table';
    end
end

valid = isfinite(rpm.timeSec) & isfinite(rpm.value);

rpm.timeSec = rpm.timeSec(valid);
rpm.value = rpm.value(valid);

[rpm.timeSec, order] = sort(rpm.timeSec);
rpm.value = rpm.value(order);

[rpm.timeSec, ia] = unique(rpm.timeSec, 'stable');
rpm.value = rpm.value(ia);

rpm.time_s = rpm.timeSec;

if isempty(rpm.timeSec)
    error('No valid RPM data found in Step03 bundle.');
end

end


function legacy = load_legacy_alignment_local(cfg, caseName)

legacy = struct();
legacy.hasAlignResult = false;
legacy.hasStep3Result = false;
legacy.alignFile = '';
legacy.step3ResultFile = '';
legacy.alignResult = struct();
legacy.step3Result = struct();

candidateDirs = {};

if isfield(cfg, 'step04_reference_legacy_output_root') && ...
        ~isempty(cfg.step04_reference_legacy_output_root)
    candidateDirs{end+1} = fullfile(cfg.step04_reference_legacy_output_root, caseName);
    candidateDirs{end+1} = cfg.step04_reference_legacy_output_root;
end

candidateDirs{end+1} = fullfile(cfg.output_root, caseName);
candidateDirs{end+1} = cfg.output_root;

if isfield(cfg, 'dataset_root') && ~isempty(cfg.dataset_root)
    candidateDirs{end+1} = fullfile(cfg.dataset_root, caseName);
end

if isfield(cfg, 'route_dir') && ~isempty(cfg.route_dir)
    candidateDirs{end+1} = fullfile(cfg.route_dir, 'output', caseName);
    candidateDirs{end+1} = fullfile(cfg.route_dir, 'outputs', caseName);
    candidateDirs{end+1} = cfg.route_dir;
end

alignPatterns = { ...
    'Step3_Spectrum_RPM_Alignment_20250527.mat', ...
    '*Spectrum_RPM_Alignment*.mat'};

resultPatterns = { ...
    'Step3_Result_20250527.mat'};

alignFiles = find_candidate_files_local(candidateDirs, alignPatterns);

for i = 1:numel(alignFiles)
    try
        A = load(alignFiles{i});

        if isfield(A, 'result')
            legacy.alignResult = A.result;
            legacy.hasAlignResult = true;
            legacy.alignFile = alignFiles{i};
            fprintf('Found legacy alignment file:\n  %s\n', alignFiles{i});
            break;

        elseif isfield(A, 'alignResult')
            legacy.alignResult = A.alignResult;
            legacy.hasAlignResult = true;
            legacy.alignFile = alignFiles{i};
            fprintf('Found legacy alignment file:\n  %s\n', alignFiles{i});
            break;
        end
    catch
    end
end

resultFiles = find_candidate_files_local(candidateDirs, resultPatterns);

for i = 1:numel(resultFiles)
    try
        R = load(resultFiles{i});

        if isfield(R, 'result')
            legacy.step3Result = R.result;
            legacy.hasStep3Result = true;
            legacy.step3ResultFile = resultFiles{i};
            fprintf('Found legacy Step3 result file:\n  %s\n', resultFiles{i});
            break;
        end
    catch
    end
end

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


function rpm = extract_rpm_from_legacy_align_local(alignResult)

rpm = struct();
rpm.timeSec = [];
rpm.time_s = [];
rpm.value = [];
rpm.source = 'legacy alignment';

if isfield(alignResult, 'rpm_time_s') && isfield(alignResult, 'rpm_values')
    rpm.timeSec = alignResult.rpm_time_s(:);
    rpm.value = alignResult.rpm_values(:);
elseif isfield(alignResult, 'rpmTimeSec') && isfield(alignResult, 'rpmValue')
    rpm.timeSec = alignResult.rpmTimeSec(:);
    rpm.value = alignResult.rpmValue(:);
end

valid = isfinite(rpm.timeSec) & isfinite(rpm.value);

rpm.timeSec = rpm.timeSec(valid);
rpm.value = rpm.value(valid);

[rpm.timeSec, order] = sort(rpm.timeSec);
rpm.value = rpm.value(order);

rpm.time_s = rpm.timeSec;

end


function validTimeRange = compute_valid_time_range_from_legacy_local(alignResult, rpm, rpmMin)

if isfield(alignResult, 'stft_time_s') && isfield(alignResult, 'best_tau_sec')
    tBtt = alignResult.stft_time_s(:) - alignResult.best_tau_sec;
else
    tBtt = [];
end

validRpm = isfinite(rpm.timeSec) & isfinite(rpm.value) & rpm.value > rpmMin;

if isempty(tBtt) || ~any(validRpm)
    validTimeRange = compute_valid_time_range_from_rpm_local(rpm, rpmMin);
else
    validTimeRange = [ ...
        max(min(tBtt), min(rpm.timeSec(validRpm))), ...
        min(max(tBtt), max(rpm.timeSec(validRpm)))];
end

if ~all(isfinite(validTimeRange)) || validTimeRange(2) <= validTimeRange(1)
    validTimeRange = compute_valid_time_range_from_rpm_local(rpm, rpmMin);
end

end


function validTimeRange = compute_valid_time_range_from_rpm_local(rpm, rpmMin)

valid = isfinite(rpm.timeSec) & isfinite(rpm.value) & rpm.value > rpmMin;

if ~any(valid)
    error('No valid RPM points above %.3f rpm.', rpmMin);
end

validTimeRange = [min(rpm.timeSec(valid)), max(rpm.timeSec(valid))];

end


function [strainTimeBtt, strainValue, strainSource, strainRaw] = load_strain_evidence_local( ...
    cfg, caseName, setting, legacy, validTimeRange) %#ok<INUSD>

strainTimeBtt = [];
strainValue = [];

strainSource = struct();
strainSource.source_mode = '';
strainSource.step3_file = '';
strainSource.strain_file = '';
strainSource.strain_time_offset_total_sec = NaN;

strainRaw = struct();
strainRaw.timeBtt = [];
strainRaw.rawValue = [];
strainRaw.detrendedValue = [];
strainRaw.sampleRateHz = NaN;
strainRaw.file = '';
strainRaw.timeOffsetSec = NaN;

% Highest priority: explicit override file.
if isfield(setting, 'strain_file_override') && ~isempty(setting.strain_file_override)

    file = resolve_mat_file_local(setting.strain_file_override);
    assert_strain_mat_file_candidate_local(file, 'strain_file_override');

    if ~isfile(file)
        error('strain_file_override does not exist:\n  %s', file);
    end

    [offset, offsetSource] = get_strain_offset_local(setting, legacy);

    fprintf('Use explicit strain_file_override:\n  %s\n', file);
    fprintf('Use strain time offset: %.9g s (%s)\n', offset, offsetSource);

    if isstruct(legacy) && isfield(legacy, 'hasStep3Result') && legacy.hasStep3Result && ...
            isfield(legacy.step3Result, 'strain_file') && ...
            ~isempty(legacy.step3Result.strain_file) && ...
            ~strcmpi(char(resolve_mat_file_local(legacy.step3Result.strain_file)), char(file))
        fprintf(['Legacy Step3_Result recorded a different strain file:\n  %s\n', ...
            'Step06A keeps the explicit override above and reuses only the legacy time offset.\n'], ...
            char(resolve_mat_file_local(legacy.step3Result.strain_file)));
    end

    [tRaw, yRaw] = read_strain_file_local(file, cfg, setting);

    strainTimeBtt = tRaw + offset;
    strainValue = preprocess_strain_local(yRaw, strainTimeBtt, setting.detrendMethod);

    valid = isfinite(strainTimeBtt) & isfinite(strainValue) & isfinite(yRaw);

    strainTimeBtt = strainTimeBtt(valid);
    strainValue = strainValue(valid);
    yRaw = yRaw(valid);

    [strainTimeBtt, order] = sort(strainTimeBtt);
    strainValue = strainValue(order);
    yRaw = yRaw(order);

    strainSource.source_mode = ['explicit_override_', offsetSource];
    strainSource.step3_file = get_optional_field_local(legacy, 'step3ResultFile', '');
    strainSource.strain_file = file;
    strainSource.legacy_step3_strain_file = get_legacy_step3_strain_file_local(legacy);
    strainSource.strain_time_offset_total_sec = offset;

    strainRaw.timeBtt = strainTimeBtt;
    strainRaw.rawValue = yRaw;
    strainRaw.detrendedValue = strainValue;
    strainRaw.sampleRateHz = estimate_sample_rate_local(strainTimeBtt);
    strainRaw.file = file;
    strainRaw.timeOffsetSec = offset;

    return;
end

% Second priority: legacy Step3_Result.
legacyStrainFile = '';

if legacy.hasStep3Result && isfield(legacy.step3Result, 'strain_file') && ...
        ~isempty(legacy.step3Result.strain_file)

    candidate = resolve_mat_file_local(legacy.step3Result.strain_file);

    if is_valid_strain_mat_file_candidate_local(candidate) && isfile(candidate)
        legacyStrainFile = candidate;
    elseif isfile(candidate)
        warning('Ignoring invalid legacy Step3 strain_file candidate:\n  %s', candidate);
    end
end

if ~isempty(legacyStrainFile)

    file = legacyStrainFile;
    offset = get_optional_field_local(legacy.step3Result, ...
        'strain_time_offset_total_sec', setting.manual_strain_time_offset_total_sec);

    fprintf('Use strain file from legacy Step3_Result:\n  %s\n', file);

    [tRaw, yRaw] = read_strain_file_local(file, cfg, setting);

    strainTimeBtt = tRaw + offset;
    strainValue = preprocess_strain_local(yRaw, strainTimeBtt, setting.detrendMethod);

    valid = isfinite(strainTimeBtt) & isfinite(strainValue) & isfinite(yRaw);

    strainTimeBtt = strainTimeBtt(valid);
    strainValue = strainValue(valid);
    yRaw = yRaw(valid);

    [strainTimeBtt, order] = sort(strainTimeBtt);
    strainValue = strainValue(order);
    yRaw = yRaw(order);

    strainSource.source_mode = 'legacy_Step3_Result';
    strainSource.step3_file = legacy.step3ResultFile;
    strainSource.strain_file = file;
    strainSource.strain_time_offset_total_sec = offset;

    strainRaw.timeBtt = strainTimeBtt;
    strainRaw.rawValue = yRaw;
    strainRaw.detrendedValue = strainValue;
    strainRaw.sampleRateHz = estimate_sample_rate_local(strainTimeBtt);
    strainRaw.file = file;
    strainRaw.timeOffsetSec = offset;

    return;
end

% Third priority: search strain_root.
files = find_strain_files_local(cfg, caseName, setting);
files = filter_strain_mat_file_candidates_local(files);

if isempty(files)
    error('No strain file found for case %s.', caseName);
end

idx = choose_primary_file_local(files, setting.strain_primary_keyword, ...
    setting.strain_primary_index_if_no_keyword);

file = files{idx};

[offset, offsetSource] = get_strain_offset_local(setting, legacy);

fprintf('Use searched strain file:\n  %s\n', file);
fprintf('Use strain time offset: %.9g s (%s)\n', offset, offsetSource);

[tRaw, yRaw] = read_strain_file_local(file, cfg, setting);

strainTimeBtt = tRaw + offset;
strainValue = preprocess_strain_local(yRaw, strainTimeBtt, setting.detrendMethod);

valid = isfinite(strainTimeBtt) & isfinite(strainValue) & isfinite(yRaw);

strainTimeBtt = strainTimeBtt(valid);
strainValue = strainValue(valid);
yRaw = yRaw(valid);

[strainTimeBtt, order] = sort(strainTimeBtt);
strainValue = strainValue(order);
yRaw = yRaw(order);

strainSource.source_mode = ['searched_strain_file_', offsetSource];
strainSource.step3_file = get_optional_field_local(legacy, 'step3ResultFile', '');
strainSource.strain_file = file;
strainSource.strain_time_offset_total_sec = offset;

strainRaw.timeBtt = strainTimeBtt;
strainRaw.rawValue = yRaw;
strainRaw.detrendedValue = strainValue;
strainRaw.sampleRateHz = estimate_sample_rate_local(strainTimeBtt);
strainRaw.file = file;
strainRaw.timeOffsetSec = offset;

end


function [offset, offsetSource] = get_strain_offset_local(setting, legacy)

offset = setting.manual_strain_time_offset_total_sec;
offsetSource = 'manual_offset';

if isstruct(legacy) && ...
        isfield(legacy, 'hasStep3Result') && legacy.hasStep3Result && ...
        isfield(legacy, 'step3Result') && isstruct(legacy.step3Result) && ...
        isfield(legacy.step3Result, 'strain_time_offset_total_sec') && ...
        isfinite(legacy.step3Result.strain_time_offset_total_sec)

    offset = legacy.step3Result.strain_time_offset_total_sec;
    offsetSource = 'legacy_Step3_Result_offset';
end

end


function file = get_legacy_step3_strain_file_local(legacy)

file = '';

if isstruct(legacy) && ...
        isfield(legacy, 'hasStep3Result') && legacy.hasStep3Result && ...
        isfield(legacy, 'step3Result') && isstruct(legacy.step3Result) && ...
        isfield(legacy.step3Result, 'strain_file') && ...
        ~isempty(legacy.step3Result.strain_file)
    file = char(resolve_mat_file_local(legacy.step3Result.strain_file));
end

end


function files = find_strain_files_local(cfg, caseName, setting)

dirs0 = {};

if isfield(setting, 'strain_root') && ~isempty(setting.strain_root)
    dirs0{end+1} = setting.strain_root;
end

if isfield(cfg, 'strain_root') && ~isempty(cfg.strain_root)
    dirs0{end+1} = cfg.strain_root;
end

if isfield(cfg, 'output_root') && ~isempty(cfg.output_root)
    dirs0{end+1} = fullfile(cfg.output_root, caseName);
    dirs0{end+1} = cfg.output_root;
end

if isfield(cfg, 'dataset_root') && ~isempty(cfg.dataset_root)
    dirs0{end+1} = fullfile(cfg.dataset_root, caseName);
    dirs0{end+1} = cfg.dataset_root;
end

if isfield(cfg, 'route_dir') && ~isempty(cfg.route_dir)
    dirs0{end+1} = cfg.route_dir;
end

dirs0 = unique(dirs0, 'stable');

files = find_candidate_files_local(dirs0, setting.strain_file_patterns);

if isfield(setting, 'strain_primary_keyword') && ~isempty(setting.strain_primary_keyword) && ~isempty(files)
    key = setting.strain_primary_keyword;
    hit = contains(files, key, 'IgnoreCase', true);
    files = [files(hit), files(~hit)];
end

files = unique(files, 'stable');

end


function filesOut = filter_strain_mat_file_candidates_local(files)

filesOut = {};
skipped = 0;

for i = 1:numel(files)
    candidate = resolve_mat_file_local(files{i});

    if is_valid_strain_mat_file_candidate_local(candidate) && isfile(candidate)
        filesOut{end+1} = candidate; %#ok<AGROW>
    else
        skipped = skipped + 1;
    end
end

filesOut = unique(filesOut, 'stable');

if skipped > 0
    fprintf('Filtered %d invalid strain candidate(s); only existing non-attachment .mat files are allowed.\n', ...
        skipped);
end

end


function idx = choose_primary_file_local(files, keyword, fallbackIndex)

idx = [];

if ~isempty(keyword)
    for i = 1:numel(files)
        if contains(files{i}, keyword, 'IgnoreCase', true)
            idx = i;
            break;
        end
    end
end

if isempty(idx)
    idx = max(1, min(fallbackIndex, numel(files)));
end

end


function [t, y] = read_strain_file_local(file, cfg, setting)

file = resolve_mat_file_local(file);
assert_strain_mat_file_candidate_local(file, 'strain data file');

D = load(file);

if isfield(D, 'Datas') && isnumeric(D.Datas) && size(D.Datas, 2) >= 2
    M = D.Datas;
else
    M = [];
    fn = fieldnames(D);
    bestScore = -inf;

    for i = 1:numel(fn)
        v = D.(fn{i});

        if isnumeric(v) && ismatrix(v) && size(v, 1) >= 16
            score = size(v, 1) + 1000 * double(size(v, 2) >= 2);

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

if size(M, 2) < 2
    M = [(0:size(M,1)-1).', M(:)];
end

x = M(:,1);
y = M(:,2) * setting.strain_value_scale_to_microstrain + ...
    setting.strain_value_offset_microstrain;

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


function y = preprocess_strain_local(x, t, detrendMethod)

x = x(:);
t = t(:);

finiteX = isfinite(x);

if any(finiteX)
    x(~finiteX) = median(x(finiteX), 'omitnan');
else
    x(:) = 0;
end

switch lower(char(detrendMethod))
    case 'none'
        y = x - median(x, 'omitnan');

    case 'linear'
        y = detrend(x, 'linear');

    case 'movingmedian'
        Fs = estimate_sample_rate_local(t);
        win = max(5, round(2.0 * Fs));

        if mod(win, 2) == 0
            win = win + 1;
        end

        y = x - movmedian(x, win, 'omitnan');

    otherwise
        error('Unknown detrendMethod: %s', detrendMethod);
end

y = y - median(y, 'omitnan');

end


function [allTime, allVib, bttSource] = load_btt_vibration_for_display_local( ...
    cfg, caseName, legacy, bundle, sensorIds, targetBladeId, baselineTimeSec)

[allTime, allVib, bttSource] = load_btt_vibration_from_legacy_step04_local( ...
    cfg, caseName, legacy, sensorIds, targetBladeId, baselineTimeSec);

if ~isempty(allTime)
    return;
end

[allTime, allVib, bttSource] = load_btt_vibration_from_bundle_local( ...
    bundle, sensorIds, targetBladeId, baselineTimeSec);

end


function [allTime, allVib, bttSource] = load_btt_vibration_from_legacy_step04_local( ...
    cfg, caseName, legacy, sensorIds, targetBladeId, baselineTimeSec)

allTime = [];
allVib = [];

bttSource = struct();
bttSource.source = 'legacy Step04 jilublade displacement';
bttSource.targetBladeId = targetBladeId;
bttSource.sensorIds = sensorIds;
bttSource.caseOutputDir = '';

candidateDirs = {};

if isstruct(legacy) && isfield(legacy, 'alignFile') && ~isempty(legacy.alignFile)
    candidateDirs{end+1} = fileparts(legacy.alignFile);
end

if isstruct(legacy) && isfield(legacy, 'step3ResultFile') && ~isempty(legacy.step3ResultFile)
    candidateDirs{end+1} = fileparts(legacy.step3ResultFile);
end

if isfield(cfg, 'step04_reference_legacy_output_root') && ...
        ~isempty(cfg.step04_reference_legacy_output_root)
    candidateDirs{end+1} = fullfile(cfg.step04_reference_legacy_output_root, caseName);
end

candidateDirs = unique(candidateDirs, 'stable');

caseOutputDir = '';

for i = 1:numel(candidateDirs)
    if isfolder(candidateDirs{i}) && ...
            isfile(fullfile(candidateDirs{i}, sprintf('jilublade_probe%d_vib_final.mat', sensorIds(1))))
        caseOutputDir = candidateDirs{i};
        break;
    end
end

if isempty(caseOutputDir)
    return;
end

bttSource.caseOutputDir = caseOutputDir;

for sid = sensorIds(:).'
    vibFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d_vib_final.mat', sid));

    if ~isfile(vibFile)
        warning('Missing legacy BTT displacement file: %s', vibFile);
        continue;
    end

    S = load(vibFile, 'jilublade');

    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 6
        warning('Invalid legacy BTT displacement file: %s', vibFile);
        continue;
    end

    bladeMask = S.jilublade(:, 4) == targetBladeId;
    t = S.jilublade(bladeMask, 3);
    y = S.jilublade(bladeMask, 6);

    valid = isfinite(t) & isfinite(y);
    t = t(valid);
    y = y(valid);

    if numel(t) < 10
        continue;
    end

    [t, orderIdx] = sort(t);
    y = y(orderIdx);

    dtBlade = median(diff(t), 'omitnan');

    if ~isfinite(dtBlade) || dtBlade <= 0
        dtBlade = 0.05;
    end

    baselinePts = max(5, round(baselineTimeSec / dtBlade));

    if mod(baselinePts, 2) == 0
        baselinePts = baselinePts + 1;
    end

    yVib = y - movmedian(y, baselinePts, 'omitnan');

    allTime = [allTime; t(:)]; %#ok<AGROW>
    allVib = [allVib; yVib(:)]; %#ok<AGROW>
end

[allTime, order] = sort(allTime);
allVib = allVib(order);

end


function [allTime, allVib, bttSource] = load_btt_vibration_from_bundle_local( ...
    bundle, sensorIds, targetBladeId, baselineTimeSec)

allTime = [];
allVib = [];

bttSource = struct();
bttSource.source = 'Step03 bundle Observation_Table';
bttSource.targetBladeId = targetBladeId;
bttSource.sensorIds = sensorIds;

if ~isfield(bundle, 'Observation_Table') || ~istable(bundle.Observation_Table)
    warning('No Observation_Table in Step03 bundle. BTT display trace will be empty.');
    return;
end

T = bundle.Observation_Table;

timeVar = find_table_var_local(T, {'arrival_time_s', 'time_s', 'btt_time_s', 'Time_s'});
dispVar = find_table_var_local(T, {'displacement_mm', 'vibration_mm', 'Vibration_mm', 'disp_mm'});
sensorVar = find_table_var_local(T, {'sensor_id', 'Sensor_ID', 'sensorId'});
bladeVar = find_table_var_local(T, {'blade_id', 'Blade_ID', 'bladeId'});
validVar = find_table_var_local(T, {'is_valid', 'valid', 'IsValid'});

if isempty(timeVar) || isempty(dispVar) || isempty(sensorVar) || isempty(bladeVar)
    warning('Observation_Table lacks required variables for BTT display.');
    return;
end

for sid = sensorIds(:).'
    mask = T.(sensorVar) == sid & ...
           T.(bladeVar) == targetBladeId & ...
           isfinite(T.(timeVar)) & ...
           isfinite(T.(dispVar));

    if ~isempty(validVar)
        mask = mask & logical(T.(validVar));
    end

    tt = T.(timeVar)(mask);
    yy = T.(dispVar)(mask);

    if isempty(tt)
        continue;
    end

    baseMask = tt <= min(tt) + baselineTimeSec;

    if nnz(baseMask) >= 3
        yy = yy - median(yy(baseMask), 'omitnan');
    else
        yy = yy - median(yy, 'omitnan');
    end

    allTime = [allTime; tt(:)]; %#ok<AGROW>
    allVib = [allVib; yy(:)]; %#ok<AGROW>
end

[allTime, order] = sort(allTime);
allVib = allVib(order);

end


function rawFft = raw_segmented_fft_local(strainRaw, timeRange, numSegments, peakBandHz)

tRaw = strainRaw.timeBtt(:);
xRaw = strainRaw.rawValue(:);
xDetrended = strainRaw.detrendedValue(:);
Fs = strainRaw.sampleRateHz;

validRange = isfinite(tRaw) & ...
             isfinite(xRaw) & ...
             tRaw >= timeRange(1) & ...
             tRaw <= timeRange(2);

idxAll = find(validRange);

if numel(idxAll) < 16
    error('Not enough raw strain samples inside BTT time range %.3f-%.3f s.', ...
        timeRange(1), timeRange(2));
end

i1 = idxAll(1);
i2 = idxAll(end);

x = xRaw(i1:i2);
xd = xDetrended(i1:i2);
t = tRaw(i1:i2);

N = numel(x);
numSegments = min(numSegments, floor(N / 16));
segmentLength = floor(N / numSegments);

if segmentLength < 2
    error('Raw FFT segment length is too short: %d samples.', segmentLength);
end

Nused = segmentLength * numSegments;

x = x(1:Nused);
xd = xd(1:Nused);
t = t(1:Nused);

freq = (0:floor(segmentLength / 2)).' * (Fs / segmentLength);

amplitudes = zeros(numSegments, numel(freq));
segmentCenterTimeSec = zeros(numSegments, 1);
segmentPeakFreqHz = zeros(numSegments, 1);
segmentPeakAmpMicrostrain = zeros(numSegments, 1);
segmentDetrendedMaxAbsMicrostrain = zeros(numSegments, 1);
segmentDetrendedHalfPeakToPeakMicrostrain = zeros(numSegments, 1);
segmentDetrendedP95AbsMicrostrain = zeros(numSegments, 1);

peakMask = freq >= peakBandHz(1) & freq <= peakBandHz(2);
freqForPeak = freq(peakMask);

for i = 1:numSegments
    localIdx = (i - 1) * segmentLength + (1:segmentLength);

    segmentData = x(localIdx);

    Y = fft(segmentData);
    P2 = abs(Y / segmentLength);
    P1 = P2(1:(floor(segmentLength / 2) + 1));

    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end

    P1(1) = 0;

    amplitudes(i, :) = P1(:).';
    segmentCenterTimeSec(i) = mean(t(localIdx), 'omitnan');

    if any(peakMask)
        [segmentPeakAmpMicrostrain(i), peakIdx] = max(P1(peakMask));
        segmentPeakFreqHz(i) = freqForPeak(peakIdx);
    else
        segmentPeakAmpMicrostrain(i) = NaN;
        segmentPeakFreqHz(i) = NaN;
    end

    segmentDetrended = xd(localIdx);

    halfPeakToPeak = 0.5 * ( ...
        max(segmentDetrended, [], 'omitnan') - ...
        min(segmentDetrended, [], 'omitnan'));

    absDetrended = abs(segmentDetrended);

    segmentDetrendedHalfPeakToPeakMicrostrain(i) = halfPeakToPeak;
    segmentDetrendedMaxAbsMicrostrain(i) = halfPeakToPeak;
    segmentDetrendedP95AbsMicrostrain(i) = percentile_local(absDetrended, 95);
end

rawFft = struct();

rawFft.sampleStart = i1;
rawFft.sampleEnd = i1 + Nused - 1;
rawFft.segmentLength = segmentLength;
rawFft.numSegments = numSegments;
rawFft.Nused = Nused;
rawFft.sampleRateHz = Fs;

rawFft.freqHz = freq;
rawFft.segmentCenterTimeSec = segmentCenterTimeSec;
rawFft.amplitudesMicrostrain = amplitudes;

rawFft.peakBandHz = peakBandHz;
rawFft.segmentPeakFreqHz = segmentPeakFreqHz;
rawFft.segmentPeakAmpMicrostrain = segmentPeakAmpMicrostrain;
rawFft.segmentDetrendedMaxAbsMicrostrain = segmentDetrendedMaxAbsMicrostrain;
rawFft.segmentDetrendedHalfPeakToPeakMicrostrain = segmentDetrendedHalfPeakToPeakMicrostrain;
rawFft.segmentDetrendedP95AbsMicrostrain = segmentDetrendedP95AbsMicrostrain;

% Modern aliases.
rawFft.freq_hz = rawFft.freqHz;
rawFft.segment_center_time_s = rawFft.segmentCenterTimeSec;
rawFft.amp_microstrain = rawFft.amplitudesMicrostrain.';
rawFft.segment_peak_freq_hz = rawFft.segmentPeakFreqHz;
rawFft.segment_peak_amp_microstrain = rawFft.segmentPeakAmpMicrostrain;

end


function [scoreTable, peakTable, regionTable, localFftTable, localFftResult] = ...
    detect_resonance_regions_from_raw_fft_local(rawFft, rpm, orders, P)

t = rawFft.segmentCenterTimeSec(:);
fftAmp = rawFft.segmentPeakAmpMicrostrain(:);

if isfield(rawFft, 'segmentDetrendedHalfPeakToPeakMicrostrain')
    timeAmp = rawFft.segmentDetrendedHalfPeakToPeakMicrostrain(:);
else
    timeAmp = rawFft.segmentDetrendedMaxAbsMicrostrain(:);
end

timeP95 = rawFft.segmentDetrendedP95AbsMicrostrain(:);
freq = rawFft.segmentPeakFreqHz(:);

rpmAt = interp1(rpm.timeSec, rpm.value, t, 'linear', NaN);

valid = isfinite(t) & ...
        isfinite(fftAmp) & ...
        isfinite(timeAmp) & ...
        isfinite(freq) & ...
        isfinite(rpmAt);

surfaceAmp = fftAmp;
timeAmpOut = timeAmp;

highMask = valid & surfaceAmp >= P.fftAmpThresholdMicrostrain;
lowMask = highMask;
score = surfaceAmp;

shaftFreqHz = rpmAt / 60;
orderApprox = freq ./ shaftFreqHz;

nearestOrder = NaN(size(orderApprox));
syncFreqAtOrderHz = NaN(size(orderApprox));

for i = 1:numel(orderApprox)
    if ~isfinite(orderApprox(i))
        continue;
    end

    [~, idxOrder] = min(abs(orders - orderApprox(i)));

    nearestOrder(i) = orders(idxOrder);
    syncFreqAtOrderHz(i) = nearestOrder(i) * shaftFreqHz(i);
end

scoreTable = table( ...
    t, ...
    rpmAt, ...
    freq, ...
    surfaceAmp, ...
    timeAmpOut, ...
    timeP95, ...
    orderApprox, ...
    nearestOrder, ...
    syncFreqAtOrderHz, ...
    score, ...
    highMask, ...
    lowMask, ...
    'VariableNames', { ...
    'bttTimeSec', ...
    'rpm', ...
    'peakFreqHz', ...
    'peakFftAmpMicrostrain', ...
    'peakTimeAmpMicrostrain', ...
    'peakTimeP95AbsMicrostrain', ...
    'fftPeakOrderApprox', ...
    'nearestOrder', ...
    'syncFreqAtNearestOrderHz', ...
    'resonanceScore', ...
    'isHighMask', ...
    'isLowMask'});

scoreTable.time_center_s = scoreTable.bttTimeSec;
scoreTable.peak_freq_hz = scoreTable.peakFreqHz;
scoreTable.peak_amp_microstrain = scoreTable.peakFftAmpMicrostrain;
scoreTable.is_resonance = scoreTable.isHighMask;

regionTable = build_raw_region_table_local(scoreTable, rawFft, orders, P);
peakTable = build_peak_table_from_regions_local(regionTable);
localFftTable = build_local_fft_table_from_regions_local(regionTable);

localFftResult = struct();
localFftResult.method = 'raw_segmented_fft_surface_threshold';
localFftResult.description = ['Raw strain segmented FFT amplitude and continuous intervals above ', ...
    'the surface-amplitude threshold. Detrended half peak-to-peak amplitude is only reported as check.'];
localFftResult.scoreTable = scoreTable;
localFftResult.params = P;

end


function R = build_raw_region_table_local(S, rawFft, orders, P)

edge = diff([false; S.isHighMask(:); false]);

startIdx = find(edge == 1);
endIdx = find(edge == -1) - 1;

rows = [];

for k = 1:numel(startIdx)
    segIdx = startIdx(k):endIdx(k);

    [~, rel] = max(S.peakFftAmpMicrostrain(segIdx));
    ip = segIdx(rel);

    i0 = startIdx(k);
    i1 = endIdx(k);

    dt = median(diff(S.bttTimeSec(isfinite(S.bttTimeSec))), 'omitnan');

    if ~isfinite(dt) || dt <= 0
        dt = rawFft.segmentLength / rawFft.sampleRateHz;
    end

    t0 = S.bttTimeSec(i0) - 0.5 * dt;
    t1 = S.bttTimeSec(i1) + 0.5 * dt;

    if t1 - t0 < P.minRegionDurationSec
        continue;
    end

    peakOrder = S.nearestOrder(ip);

    if ~isfinite(peakOrder)
        peakOrder = nearest_finite_order_local(S.fftPeakOrderApprox(ip), orders);
    end

    peakAmpDb = 20 * log10(max(S.peakFftAmpMicrostrain(ip), eps));

    idx = find(S.bttTimeSec >= t0 & S.bttTimeSec <= t1);

    rows = [rows; ...
        k, ...
        t0, ...
        t1, ...
        t1 - t0, ...
        S.bttTimeSec(ip), ...
        S.peakFreqHz(ip), ...
        peakOrder, ...
        peakAmpDb, ...
        S.peakFftAmpMicrostrain(ip), ...
        S.peakTimeAmpMicrostrain(ip), ...
        S.rpm(ip), ...
        S.syncFreqAtNearestOrderHz(ip), ...
        S.fftPeakOrderApprox(ip), ...
        S.resonanceScore(ip), ...
        numel(idx), ...
        double(overlap_duration_local([t0, t1], P.preferredRegionTimeSec))]; %#ok<AGROW>
end

if isempty(rows)
    R = empty_region_table_local();
    return;
end

R = array2table(rows, 'VariableNames', { ...
    'regionId', ...
    'bttStartSec', ...
    'bttEndSec', ...
    'regionDuration', ...
    'peakBttTimeSec', ...
    'dominantFreqHz', ...
    'dominantOrder', ...
    'peakAmpDb', ...
    'peakFftAmpMicrostrain', ...
    'peakTimeAmpMicrostrain', ...
    'peakRpm', ...
    'syncFreqAtRegionOrderHz', ...
    'fftPeakOrderApprox', ...
    'regionPeakScore', ...
    'peakCount', ...
    'preferredOverlapSec'});

R = merge_raw_regions_local(R, P.regionMergeGapSec, P.maxRegionWidthSec);

R.ordersIncluded = strings(height(R), 1);

for i = 1:height(R)
    inRegion = S.bttTimeSec >= R.bttStartSec(i) & ...
               S.bttTimeSec <= R.bttEndSec(i) & ...
               isfinite(S.nearestOrder);

    ord = unique(S.nearestOrder(inRegion)).';

    if isempty(ord)
        R.ordersIncluded(i) = "";
    else
        R.ordersIncluded(i) = join(string(ord), ",");
    end
end

R.regionStart = R.bttStartSec;
R.regionEnd = R.bttEndSec;
R.regionPeakTime = R.peakBttTimeSec;
R.regionPeakFreqHz = R.dominantFreqHz;
R.regionPeakOrder = R.dominantOrder;

R = sortrows(R, {'preferredOverlapSec', 'regionPeakScore'}, {'descend', 'descend'});

R.regionScanRank = (1:height(R)).';
R.regionId = R.regionScanRank;

R.region_id = R.regionId;
R.time_start_s = R.bttStartSec;
R.time_end_s = R.bttEndSec;
R.time_center_s = 0.5 * (R.bttStartSec + R.bttEndSec);
R.duration_s = R.regionDuration;
R.peak_time_s = R.peakBttTimeSec;
R.peak_freq_hz = R.dominantFreqHz;
R.peak_amp_microstrain = R.peakFftAmpMicrostrain;
R.peak_rpm = R.peakRpm;
R.dominant_order = R.dominantOrder;

end


function R2 = merge_raw_regions_local(R, mergeGapSec, maxRegionWidthSec)

if isempty(R)
    R2 = R;
    return;
end

R = sortrows(R, 'bttStartSec', 'ascend');

R2 = R([], :);

newId = 0;
i = 1;

while i <= height(R)
    rows = i;
    tEnd = R.bttEndSec(i);
    j = i + 1;

    while j <= height(R)
        mergedStart = min(R.bttStartSec(rows(1)), R.bttStartSec(j));
        mergedEnd = max(tEnd, R.bttEndSec(j));

        if R.bttStartSec(j) - tEnd <= mergeGapSec && ...
                mergedEnd - mergedStart <= maxRegionWidthSec

            rows(end + 1) = j; %#ok<AGROW>
            tEnd = mergedEnd;
            j = j + 1;
        else
            break;
        end
    end

    sub = R(rows, :);

    [~, idxBest] = max(sub.regionPeakScore);

    newId = newId + 1;

    newRow = sub(idxBest, :);
    newRow.regionId = newId;
    newRow.bttStartSec = min(sub.bttStartSec);
    newRow.bttEndSec = max(sub.bttEndSec);
    newRow.regionDuration = newRow.bttEndSec - newRow.bttStartSec;
    newRow.peakCount = sum(sub.peakCount);
    newRow.preferredOverlapSec = max(sub.preferredOverlapSec);

    R2 = [R2; newRow]; %#ok<AGROW>

    i = j;
end

end


function T = build_peak_table_from_regions_local(R)

if isempty(R)
    T = empty_peak_table_local();
    return;
end

T = table( ...
    R.dominantOrder, ...
    R.regionScanRank, ...
    R.peakBttTimeSec, ...
    R.dominantFreqHz, ...
    R.peakFftAmpMicrostrain, ...
    R.peakAmpDb, ...
    R.peakRpm, ...
    R.peakTimeAmpMicrostrain, ...
    R.regionId, ...
    'VariableNames', { ...
    'order', ...
    'rankInOrder', ...
    'bttTimeSec', ...
    'orderFreqHz', ...
    'amp', ...
    'ampDb', ...
    'rpm', ...
    'timeAmpMicrostrain', ...
    'regionId'});

T.region_id = T.regionId;
T.time_s = T.bttTimeSec;
T.freq_hz = T.orderFreqHz;
T.amp_microstrain = T.amp;

end


function T = build_local_fft_table_from_regions_local(R)

if isempty(R)
    T = empty_local_fft_table_local();
    return;
end

T = table( ...
    R.regionId, ...
    R.regionStart, ...
    R.regionEnd, ...
    R.regionPeakTime, ...
    R.regionPeakOrder, ...
    R.regionPeakFreqHz, ...
    R.regionPeakFreqHz, ...
    R.peakFftAmpMicrostrain, ...
    R.fftPeakOrderApprox, ...
    R.syncFreqAtRegionOrderHz, ...
    R.peakTimeAmpMicrostrain, ...
    'VariableNames', { ...
    'regionId', ...
    'regionStart', ...
    'regionEnd', ...
    'regionPeakTime', ...
    'regionPeakOrder', ...
    'regionPeakFreqHz', ...
    'fftPeakFreqHz', ...
    'fftPeakAmpMicrostrain', ...
    'fftPeakOrderApprox', ...
    'syncFreqAtRegionOrderHz', ...
    'peakTimeAmpMicrostrain'});

T.region_id = T.regionId;
T.time_start_s = T.regionStart;
T.time_end_s = T.regionEnd;
T.peak_freq_hz = T.fftPeakFreqHz;
T.peak_amp_microstrain = T.fftPeakAmpMicrostrain;

end


function T = empty_region_table_local()

T = table();

numericNames = { ...
    'regionId', ...
    'bttStartSec', ...
    'bttEndSec', ...
    'regionDuration', ...
    'peakBttTimeSec', ...
    'dominantFreqHz', ...
    'dominantOrder', ...
    'peakAmpDb', ...
    'peakFftAmpMicrostrain', ...
    'peakTimeAmpMicrostrain', ...
    'peakRpm', ...
    'syncFreqAtRegionOrderHz', ...
    'fftPeakOrderApprox', ...
    'regionPeakScore', ...
    'peakCount', ...
    'preferredOverlapSec', ...
    'regionStart', ...
    'regionEnd', ...
    'regionPeakTime', ...
    'regionPeakFreqHz', ...
    'regionPeakOrder', ...
    'regionScanRank', ...
    'region_id', ...
    'time_start_s', ...
    'time_end_s', ...
    'time_center_s', ...
    'duration_s', ...
    'peak_time_s', ...
    'peak_freq_hz', ...
    'peak_amp_microstrain', ...
    'peak_rpm', ...
    'dominant_order'};

for i = 1:numel(numericNames)
    T.(numericNames{i}) = zeros(0, 1);
end

T.ordersIncluded = strings(0, 1);

end


function T = empty_peak_table_local()

T = table( ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
    'order', ...
    'rankInOrder', ...
    'bttTimeSec', ...
    'orderFreqHz', ...
    'amp', ...
    'ampDb', ...
    'rpm', ...
    'timeAmpMicrostrain', ...
    'regionId'});

T.region_id = zeros(0,1);
T.time_s = zeros(0,1);
T.freq_hz = zeros(0,1);
T.amp_microstrain = zeros(0,1);

end


function T = empty_local_fft_table_local()

T = table();

names = { ...
    'regionId', ...
    'regionStart', ...
    'regionEnd', ...
    'regionPeakTime', ...
    'regionPeakOrder', ...
    'regionPeakFreqHz', ...
    'fftPeakFreqHz', ...
    'fftPeakAmpMicrostrain', ...
    'fftPeakOrderApprox', ...
    'syncFreqAtRegionOrderHz', ...
    'peakTimeAmpMicrostrain', ...
    'region_id', ...
    'time_start_s', ...
    'time_end_s', ...
    'peak_freq_hz', ...
    'peak_amp_microstrain'};

for i = 1:numel(names)
    T.(names{i}) = zeros(0,1);
end

end


function plot_dashboard_local(strainTimeBtt, strainValue, bttTime, bttVib, targetBladeId, ...
    rpm, rawFft, peakTable, regionTable, timeRange, freqBand, visible, figDir, saveFigures, setting)

rawTime = rawFft.segmentCenterTimeSec(:).';

freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);
rawAmp = rawFft.amplitudesMicrostrain(:, freqMask).';

fig = figure( ...
    'Name', 'Step06A old-Step04-style resonance dashboard', ...
    'Color', 'w', ...
    'Position', [50, 35, 1420, 940], ...
    'NumberTitle', 'off', ...
    'Visible', visible);

tiledlayout(fig, 4, 1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

ax1 = nexttile;

plot_target_blade_btt_and_strain_local( ...
    ax1, ...
    strainTimeBtt, ...
    strainValue, ...
    bttTime, ...
    bttVib, ...
    targetBladeId, ...
    regionTable, ...
    timeRange, ...
    setting.max_plot_points);

ax2 = nexttile;

plot_speed_trace_local( ...
    ax2, ...
    rpm, ...
    peakTable, ...
    regionTable, ...
    timeRange, ...
    setting.rpm_plot_y_margin_ratio);

ax3 = nexttile;

imagesc(ax3, rawTime, rawFft.freqHz(freqMask), rawAmp);
axis(ax3, 'xy');

grid(ax3, 'on');
box(ax3, 'on');

xlim(ax3, timeRange);
ylim(ax3, freqBand);

colormap(ax3, high_contrast_colormap_local(256));
clim_local(ax3, raw_amplitude_color_limits_local(rawAmp, 99.5));

cb = colorbar(ax3);
cb.Label.String = 'Amplitude (\mue)';

ylabel(ax3, 'Frequency (Hz)');
title(ax3, sprintf('3) Raw segmented FFT %.0f-%.0f Hz with direct %.0f-\\mue regions', ...
    freqBand(1), freqBand(2), setting.resonanceFftAmpThresholdMicrostrain));

hold(ax3, 'on');

draw_resonance_ellipses_local(ax3, regionTable, freqBand);
plot_major_peak_points_local(ax3, peakTable, regionTable, freqBand);

ax4 = nexttile;

if ~isempty(regionTable)
    bh = bar(ax4, regionTable.regionId, regionTable.peakAmpDb, 0.65, 'FaceColor', 'flat');
    bh.CData = region_palette_local(height(regionTable));

    ylim(ax4, [0, max(regionTable.peakAmpDb) * 1.22 + eps]);

    for i = 1:height(regionTable)
        txt = sprintf('R%d: %dX\n%.1f-%.1fs', ...
            regionTable.regionId(i), ...
            regionTable.dominantOrder(i), ...
            regionTable.bttStartSec(i), ...
            regionTable.bttEndSec(i));

        text(ax4, ...
            regionTable.regionId(i), ...
            regionTable.peakAmpDb(i), ...
            txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 8, ...
            'Interpreter', 'none');
    end
end

grid(ax4, 'on');
box(ax4, 'on');

xlabel(ax4, 'Region ID');
ylabel(ax4, 'Peak raw FFT amplitude (dB)');
title(ax4, '4) Region summary: dominant order and peak level');

linkaxes([ax1, ax2, ax3], 'x');

if saveFigures
    save_figure_local(fig, figDir, 'Step06A_Fig01_Strain_RPM_FFT_Dashboard');
end

end


function plot_target_blade_btt_and_strain_local(ax, strainTimeBtt, strainValue, ...
    bttTime, bttVib, targetBladeId, regionTable, timeRange, maxPlotPoints)

validBtt = bttTime >= timeRange(1) & ...
           bttTime <= timeRange(2) & ...
           isfinite(bttVib);

tBtt = bttTime(validBtt);
yBtt = bttVib(validBtt);

idxThin = decimate_index_local(numel(tBtt), maxPlotPoints);

if isempty(tBtt)
    warning('No BTT samples for target blade %d inside display range.', targetBladeId);
    tBtt = NaN;
    yBtt = NaN;
    idxThin = 1;
end

yyaxis(ax, 'left');

plot(ax, tBtt(idxThin), yBtt(idxThin), '.', ...
    'Color', [0.18 0.18 0.18], ...
    'MarkerSize', 3);

hold(ax, 'on');

if numel(tBtt) > 10
    bttEnv = moving_rms_local(yBtt, 0.15, estimate_sample_rate_local(tBtt));
    idxEnv = decimate_index_local(numel(tBtt), 12000);

    plot(ax, tBtt(idxEnv), bttEnv(idxEnv), 'r-', 'LineWidth', 1.0);
    plot(ax, tBtt(idxEnv), -bttEnv(idxEnv), 'r-', 'LineWidth', 1.0);
end

ylabel(ax, sprintf('Blade %d BTT vib. (mm)', targetBladeId));
ax.YColor = [0.10 0.10 0.10];

yyaxis(ax, 'right');

validStrain = strainTimeBtt >= timeRange(1) & ...
              strainTimeBtt <= timeRange(2) & ...
              isfinite(strainValue);

tStrain = strainTimeBtt(validStrain);
yStrain = strainValue(validStrain);

if ~isempty(tStrain)
    yStrain = yStrain - median(yStrain, 'omitnan');

    strainScale = percentile_local(abs(yStrain), 98);

    if ~isfinite(strainScale) || strainScale <= 0
        strainScale = max(abs(yStrain), [], 'omitnan');
    end

    if ~isfinite(strainScale) || strainScale <= 0
        strainScale = 1;
    end

    yStrainNorm = yStrain / strainScale;

    idxStrain = decimate_index_local(numel(tStrain), maxPlotPoints);

    plot(ax, tStrain(idxStrain), yStrainNorm(idxStrain), ...
        'Color', [0.05 0.35 0.85], ...
        'LineWidth', 0.6);

    if numel(tStrain) > 10
        strainEnv = moving_rms_local(yStrainNorm, 0.15, estimate_sample_rate_local(tStrain));
        idxStrainEnv = decimate_index_local(numel(tStrain), 12000);

        plot(ax, tStrain(idxStrainEnv), strainEnv(idxStrainEnv), ...
            'Color', [0.05 0.35 0.85], ...
            'LineWidth', 1.0);
    end
end

ylabel(ax, 'Strain norm.');
ax.YColor = [0.05 0.35 0.85];

yyaxis(ax, 'left');

draw_region_patches_local(ax, regionTable);

grid(ax, 'on');
box(ax, 'on');

xlim(ax, timeRange);

title(ax, sprintf('1) Strain response and BTT vibration of target blade %d', targetBladeId));

end


function plot_speed_trace_local(ax, rpm, peakTable, regionTable, timeRange, marginRatio)

valid = rpm.timeSec >= timeRange(1) & ...
        rpm.timeSec <= timeRange(2) & ...
        isfinite(rpm.value);

t = rpm.timeSec(valid);
y = rpm.value(valid);

plot(ax, t, y, 'k-', 'LineWidth', 1.0);
hold(ax, 'on');

if numel(y) > 9
    smoothPts = max(5, round(numel(y) / 350));

    plot(ax, t, movmean(y, smoothPts, 'omitnan'), ...
        'Color', [0.85 0.20 0.10], ...
        'LineWidth', 1.3);
end

draw_region_patches_local(ax, regionTable);
plot_peak_rpm_points_local(ax, peakTable, regionTable);

grid(ax, 'on');
box(ax, 'on');

xlim(ax, timeRange);

apply_rpm_ylim_local(ax, rpm, timeRange, marginRatio);

ylabel(ax, 'RPM');
title(ax, '2) OPR speed variation and selected resonance peak times');

end


function plot_peak_rpm_points_local(ax, T, R)

if isempty(T)
    return;
end

colors = colors_for_peak_rows_local(T, R);
sizes = scale_marker_sizes_local(T.ampDb, 18, 70);

for k = 1:height(T)
    if ~isfinite(T.rpm(k))
        continue;
    end

    scatter(ax, T.bttTimeSec(k), T.rpm(k), sizes(k), colors(k, :), ...
        'filled', ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.7, ...
        'MarkerFaceAlpha', 0.85, ...
        'HandleVisibility', 'off');
end

end


function plot_major_peak_points_local(ax, T, R, freqBand)

if isempty(T)
    return;
end

colors = colors_for_peak_rows_local(T, R);
sizes = scale_marker_sizes_local(T.ampDb, 30, 100);

for k = 1:height(T)
    y = T.orderFreqHz(k);

    if y < freqBand(1) || y > freqBand(2)
        continue;
    end

    scatter(ax, T.bttTimeSec(k), y, sizes(k), colors(k, :), ...
        'filled', ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.8, ...
        'MarkerFaceAlpha', 0.82, ...
        'HandleVisibility', 'off');

    text(ax, T.bttTimeSec(k), min(freqBand(2)-5, y+7), sprintf('%dX', T.order(k)), ...
        'Color', colors(k, :), ...
        'FontWeight', 'bold', ...
        'FontSize', 7, ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0 0 0], ...
        'Margin', 0.5, ...
        'Interpreter', 'none');
end

end


function draw_resonance_ellipses_local(ax, R, freqBand)

if isempty(R)
    return;
end

colors = region_palette_local(height(R));

theta = linspace(0, 2*pi, 240);
freqSpan = diff(freqBand);

for i = 1:height(R)
    c = colors(i, :);

    cx = R.peakBttTimeSec(i);
    cy = R.dominantFreqHz(i);

    rx = max(2.0, 0.58 * (R.bttEndSec(i) - R.bttStartSec(i) + 1));
    ry = max(12.0, 0.070 * freqSpan);

    x = cx + rx * cos(theta);
    y = cy + ry * sin(theta);

    outside = y < freqBand(1) | y > freqBand(2);

    x(outside) = NaN;
    y(outside) = NaN;

    plot(ax, x, y, 'w-', 'LineWidth', 2.4, 'HandleVisibility', 'off');
    plot(ax, x, y, '--', 'Color', c, 'LineWidth', 1.8, 'HandleVisibility', 'off');

    plot(ax, cx, cy, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', c, ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    text(ax, cx, min(freqBand(2)-6, cy+ry+8), ...
        sprintf('R%d  %dX  %.0f Hz', R.regionId(i), R.dominantOrder(i), R.dominantFreqHz(i)), ...
        'Color', c, ...
        'FontWeight', 'bold', ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0 0 0], ...
        'Margin', 1, ...
        'Interpreter', 'none');
end

end


function draw_region_patches_local(ax, R)

if isempty(R)
    return;
end

yl = ylim(ax);
colors = region_palette_local(height(R));

for i = 1:height(R)
    c = colors(i, :);

    patch(ax, ...
        [R.bttStartSec(i), R.bttEndSec(i), R.bttEndSec(i), R.bttStartSec(i)], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        c, ...
        'FaceAlpha', 0.055, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    xline(ax, R.bttStartSec(i), ':', ...
        'Color', c, ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    xline(ax, R.bttEndSec(i), ':', ...
        'Color', c, ...
        'LineWidth', 0.8, ...
        'HandleVisibility', 'off');

    xline(ax, R.peakBttTimeSec(i), '--', sprintf('R%d', R.regionId(i)), ...
        'Color', c, ...
        'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
end

ylim(ax, yl);

end


function plot_3d_auxiliary_local(rawFft, timeRange, freqBand, regionTable, visible, figDir, saveFigures)

freqMask = rawFft.freqHz >= freqBand(1) & rawFft.freqHz <= freqBand(2);

t = rawFft.segmentCenterTimeSec(:);
f = rawFft.freqHz(freqMask);
A = rawFft.amplitudesMicrostrain(:, freqMask).';

[T,F] = meshgrid(t, f);

fig = figure( ...
    'Name', 'Step06A raw segmented FFT 3D auxiliary', ...
    'Color', 'w', ...
    'Position', [100, 80, 1320, 760], ...
    'NumberTitle', 'off', ...
    'Visible', visible);

surf(T, F, A, 'EdgeColor', 'none');

view(48, 32);

grid on;
box on;

xlim(timeRange);
ylim(freqBand);

xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
zlabel('Amplitude (\mue)');

title('Raw segmented strain FFT surface');

colormap(high_contrast_colormap_local(256));
colorbar;

hold on;

for i = 1:height(regionTable)
    plot3( ...
        regionTable.peakBttTimeSec(i), ...
        regionTable.dominantFreqHz(i), ...
        regionTable.peakFftAmpMicrostrain(i), ...
        'wo', ...
        'MarkerFaceColor', region_palette_local(height(regionTable), i), ...
        'MarkerSize', 7, ...
        'LineWidth', 1.2);
end

if saveFigures
    save_figure_local(fig, figDir, 'Step06A_Fig02_RawSegmentedFFT_3D');
end

end


function name = find_table_var_local(T, candidates)

name = '';

for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return;
    end
end

end


function file = resolve_existing_mat_or_empty_local(file0, label)

file = '';

if nargin < 2 || isempty(label)
    label = 'file';
end

if isempty(file0)
    return;
end

candidate = resolve_mat_file_local(file0);

if ~is_valid_strain_mat_file_candidate_local(candidate)
    warning('Ignoring invalid %s; expected a non-attachment .mat file:\n  %s', ...
        label, candidate);
    return;
end

if isfile(candidate)
    file = candidate;
else
    warning('Configured %s was not found, so Step06A will search strain_root instead:\n  %s', ...
        label, candidate);
end

end


function assert_strain_mat_file_candidate_local(file, label)

if nargin < 2 || isempty(label)
    label = 'strain file';
end

if is_forbidden_runtime_attachment_path_local(file)
    error(['Refusing to use %s from a Codex/runtime attachment path. ', ...
        'Run the project copy of Step06A and set strain_file_override or cfg.strain_root to the real data directory:\n  %s'], ...
        label, char(file));
end

[~, ~, ext] = fileparts(char(file));

if ~strcmpi(ext, '.mat')
    error('Refusing to use %s because it is not a .mat file:\n  %s', ...
        label, char(file));
end

end


function tf = is_valid_strain_mat_file_candidate_local(file)

if isempty(file)
    tf = false;
    return;
end

[~, ~, ext] = fileparts(char(file));

tf = strcmpi(ext, '.mat') && ...
    ~is_forbidden_runtime_attachment_path_local(file);

end


function tf = is_forbidden_runtime_attachment_path_local(file)

if isempty(file)
    tf = false;
    return;
end

bs = char(92);
p = lower(strrep(char(file), '/', bs));

tf = contains(p, [bs '.codex' bs 'attachments' bs]) || ...
     contains(p, [bs '.codex' bs 'attachments']) || ...
     contains(p, [bs 'codex' bs 'attachments' bs]) || ...
     contains(p, 'pasted-text.txt');

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


function order = nearest_finite_order_local(orderApprox, orders)

order = NaN;

if ~isfinite(orderApprox) || isempty(orders)
    return;
end

[~, idx] = min(abs(orders - orderApprox));
order = orders(idx);

end


function dt = overlap_duration_local(a, b)

if isempty(b) || numel(b) < 2
    dt = 0;
    return;
end

dt = max(0, min(a(2), b(2)) - max(a(1), b(1)));

end


function sizes = scale_marker_sizes_local(x, smallSize, largeSize)

x = x(:);

if isempty(x) || all(~isfinite(x))
    sizes = smallSize * ones(size(x));
    return;
end

lo = min(x(isfinite(x)));
hi = max(x(isfinite(x)));

if hi <= lo
    sizes = mean([smallSize, largeSize]) * ones(size(x));
else
    sizes = smallSize + (largeSize - smallSize) * (x - lo) / (hi - lo);
end

sizes(~isfinite(sizes)) = smallSize;

end


function colors = colors_for_peak_rows_local(T, R)

if isempty(T)
    colors = zeros(0,3);
    return;
end

palette = region_palette_local(max(height(R), 1));

colors = zeros(height(T), 3);

for i = 1:height(T)
    rid = T.regionId(i);

    idx = find(R.regionId == rid, 1, 'first');

    if isempty(idx)
        colors(i,:) = [0.3 0.3 0.3];
    else
        colors(i,:) = palette(idx,:);
    end
end

end


function colors = region_palette_local(n, idx)

if nargin < 2
    idx = [];
end

base = lines(max(n, 1));

if isempty(idx)
    colors = base(1:n, :);
else
    colors = base(idx, :);
end

end


function colorLimits = raw_amplitude_color_limits_local(x, highPct)

x = x(:);
x = x(isfinite(x));

if isempty(x)
    colorLimits = [0 1];
    return;
end

hi = percentile_local(x, highPct);

if ~isfinite(hi) || hi <= 0
    hi = max(x);
end

if ~isfinite(hi) || hi <= 0
    hi = 1;
end

colorLimits = [0 hi];

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


function idx = decimate_index_local(n, maxPoints)

if n <= 0
    idx = [];
elseif n <= maxPoints
    idx = 1:n;
else
    idx = unique(round(linspace(1, n, maxPoints)));
end

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


function y = moving_rms_local(x, windowSec, Fs)

x = x(:);

if isempty(x)
    y = x;
    return;
end

if ~isfinite(Fs) || Fs <= 0
    y = abs(x);
    return;
end

win = max(3, round(windowSec * Fs));

if mod(win, 2) == 0
    win = win + 1;
end

y = sqrt(movmean(x.^2, win, 'omitnan'));

end


function apply_rpm_ylim_local(ax, rpm, tRange, marginRatio)

if nargin < 3 || isempty(tRange)
    tRange = [-inf inf];
end

if nargin < 4 || isempty(marginRatio)
    marginRatio = 0.08;
end

t = rpm.timeSec(:);
v = rpm.value(:);

valid = isfinite(t) & ...
        isfinite(v) & ...
        v > 0 & ...
        t >= tRange(1) & ...
        t <= tRange(2);

if nnz(valid) < 2
    return;
end

vv = v(valid);

lo = percentile_local(vv, 1);
hi = percentile_local(vv, 99);

if ~isfinite(lo) || ~isfinite(hi)
    lo = min(vv);
    hi = max(vv);
end

if hi <= lo
    mid = mean(vv, 'omitnan');
    margin = max(10, 0.05 * abs(mid));

    lo = mid - margin;
    hi = mid + margin;
else
    margin = marginRatio * (hi - lo);

    lo = lo - margin;
    hi = hi + margin;
end

if isfinite(lo) && isfinite(hi) && hi > lo
    ylim(ax, [lo hi]);
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


function clim_local(ax, lims)

try
    clim(ax, lims);
catch
    caxis(ax, lims);
end

end


function save_figure_local(fig, figDir, name)

ensure_dir_local(figDir);

figPath = fullfile(figDir, [name '.fig']);
pngPath = fullfile(figDir, [name '.png']);

try
    oldDir = pwd;
    cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(figDir);

    savefig(fig, [name '.fig']);
catch ME
    warning('savefig failed for:\n  %s\nReason: %s\nPNG export will continue.', ...
        figPath, ME.message);
end

try
    exportgraphics(fig, pngPath, 'Resolution', 300);
catch ME1
    try
        print(fig, pngPath, '-dpng', '-r300');
    catch ME2
        warning('PNG export failed for:\n  %s\nexportgraphics: %s\nprint: %s', ...
            pngPath, ME1.message, ME2.message);
    end
end

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
