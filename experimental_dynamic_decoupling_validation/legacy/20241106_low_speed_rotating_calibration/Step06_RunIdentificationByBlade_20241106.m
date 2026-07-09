%% Step06_RunIdentificationByBlade_20241106
% Run the direct low-speed-template identification by physical blade.
% This script now extracts the needed high-speed waveform slices inline and
% directly applies the same VP top-K synchronous waveform identification idea
% used by the 20250527/20251222 low-speed routes:
%
%   V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi))
%
% First-order VP screens EO candidates; the final result comes from the full
% template waveform objective.

if ~exist('STEP06_SKIP_CLEAR_FOR_DRIVER', 'var') || ~STEP06_SKIP_CLEAR_FOR_DRIVER
    clear; close all; clc;
end

P = NewFlow_Config_20241106();

% Step06 local parameter block.
% Edit here first. This script is intentionally self-tunable so blade-wise
% debugging does not require jumping back to NewFlow_Config_20241106.
S06 = struct();
S06.analysisSensors = [2 5 7];
S06.startTimeSec = 75.0;
S06.targetLaps = P.region.targetLaps;
S06.blades = 4;                     % default: first debug blade 4
S06.inlineHighSpeedExtraction = true;
S06.autoBuildCurrentRegionNumbering = true;
S06.saveAutoBuiltRegionAndNumbering = true;
S06.preferPeakCache = true;
S06.autoBuildPeakCacheIfMissing = true;
S06.windowLaps = 4;
S06.slidingStepLaps = 1;
S06.pulsePadSec = 0;
S06.dynamicWindowMode = 'legacy_row_bounds';
S06.windowRule = 't_start=jilublade(row,1)-pulsePadSec; t_end=jilublade(row,2)+pulsePadSec';
S06.windowIds = [];                 % e.g. 1:6
S06.lapRange = [];                  % e.g. [4 12]
S06.windowCenterTimeRangeSec = [];  % e.g. [75.05 75.20]
S06.freqSearchHz = [0 1000];
S06.topKEO = 3;
S06.amplitudeLimitMM = 0.50;
S06.dxCLimitMM = 0.35;
S06.sensorEtaLimitMM = 0.20;       % global hard ceiling
S06.sensorEtaAdaptiveLimitEnable = true;
S06.sensorEtaAdaptiveQuantile = 85;
S06.sensorEtaAdaptiveSafetyFactor = 1.20;
S06.sensorEtaAdaptiveMinMM = 0.04;
S06.sensorEtaAdaptiveMaxMM = S06.sensorEtaLimitMM;
S06.sensorEtaRegWeightVPerMM = 0.02;
S06.overshootPenaltyWeight = 50;
S06.pulseSelectionMode = 'single';
S06.domainSelectionMode = 'hard';
S06.domainMarginMM = 0.02;
S06.domainSoftMarginMM = 0.10;
S06.queryGuardMode = 'adaptive';
S06.queryGuardMM = S06.amplitudeLimitMM + S06.dxCLimitMM + 0.05;
S06.queryGuardQuantile = 95;
S06.queryGuardSafetyMM = 0.05;
S06.queryGuardMinMM = 0.12;
S06.queryGuardMaxMM = S06.queryGuardMM;
S06.sensorDomainExpandMMTable = [];
S06.sensorQueryGuardScaleTable = [5 0.80; 7 0.75];
% sensorDomainExpandMMTable: positive = wider, negative = narrower, applied
% symmetrically to both low/high domains.
% sensorQueryGuardScaleTable: scale < 1 widens both sides by shrinking the
% adaptive query guard for that sensor only.
S06.dynamicEffectiveMode = 'gradient';
S06.dynamicTemplateGradientMinRatio = 0.08;
S06.dynamicTimeGradientMinRatio = 0.15;
S06.dynamicPeakQuantile = 85;
S06.defaultSensorThreshold = 0.5;
S06.weightFloor = 0.05;
S06.phaseSafeExpansionEnable = true;
S06.phaseSafeReferenceMode = 'current_window'; % 'prev_window' | 'current_window'
S06.phaseSafeFallbackMode = 'final_core';      % current-window phase-safe reference uses fitCore
S06.phaseSafeRefreshEvery = 0;                 % only used when phaseSafeReferenceMode='prev_window'
S06.phaseSafeMarginMM = 0.02;
S06.phaseSafeMinPointCount = 20;
S06.phaseSafeMinPointCountPerSensor = 6;
S06.phaseSafePrevMaxClampFraction = 0.10;
S06.phaseSafePrevMinFinalGapRatio = 0.01;
S06.chooseBetterPass = true;        % true: keep lower-RMSE pass between core/expanded
S06.betterPassRmseToleranceV = 0;   % expanded must beat core by this margin
S06.modelMode = 'direct_only';      % direct_only | direct_and_gap | gap_only
S06.gapLibrarySensors = [5 7];
S06.gapLibraryTopKEO = 3;
S06.gapLibraryAmplitudeLimitMM = S06.amplitudeLimitMM;
S06.gapLibraryDxLimitMM = 0.60;
S06.gapLibraryDeltaGapLimitMM = 0.80;
S06.gapLibraryStaticRegWeight = 0.50;
S06.gapLibrarySensorEtaFallbackLimitMM = 0.06;
S06.gapLibrarySensorEtaRegWeightMV = 0.35;
S06.gapLibraryEOCandidateMode = 'direct_plus_vp';     % vp | direct | direct_plus_vp
S06.gapLibraryPrevWindowCandidateMode = 'soft';       % off | soft
S06.gapLibraryPrevWindowRefreshEvery = 5;
S06.gapLibraryBundleSourceMode = 'raw_gap_observation'; % direct_bundle_reuse | raw_gap_observation
S06.gapLibraryModelNames = {'fixed','gap_only','gap_tilt'};
S06.gapLibraryDeltaMuLimit = 0.030;
S06.gapLibraryCorrectedLibFile = '';
S06.outputLabel = 'B4_only';
S06.viewEnable = true;

analysisSensorsOverride = parse_step06_int_env_local('STEP06_ANALYSIS_SENSORS');
if ~isempty(analysisSensorsOverride)
    S06.analysisSensors = analysisSensorsOverride;
end
bladesOverride = parse_step06_int_env_local('STEP06_TARGET_BLADES');
if ~isempty(bladesOverride)
    S06.blades = bladesOverride;
end
S06.startTimeSec = parse_step06_scalar_env_local('STEP06_START_TIME_SEC', S06.startTimeSec);
S06.targetLaps = round(parse_step06_scalar_env_local('STEP06_TARGET_LAPS', S06.targetLaps));
S06.windowLaps = round(parse_step06_scalar_env_local('STEP06_WINDOW_LAPS', S06.windowLaps));
S06.slidingStepLaps = round(parse_step06_scalar_env_local('STEP06_SLIDING_STEP_LAPS', S06.slidingStepLaps));
S06.outputLabel = parse_step06_text_env_local('STEP06_OUTPUT_LABEL', S06.outputLabel);

windowIdsOverride = parse_step06_int_env_local('STEP06_WINDOW_IDS');
if ~isempty(windowIdsOverride)
    S06.windowIds = windowIdsOverride;
end
P = apply_step06_local_options_local(P, S06);

[WaveformLibrary, FilteredWaveformLibrary, waveformSourceInfo] = ...
    load_step06_waveform_inputs_local(P);

outCsv = P.files.identificationSummary;
outMat = P.files.identificationResult;
outDir = fileparts(outCsv);
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

cfg = build_identification_config_local(P);
method = build_method_info_local(cfg);

trendRows = [];
windowResults = repmat(empty_window_result_local(), 0, 1);
bestByBlade = repmat(struct('blade_id', NaN, 'weighted_voltage_rmse', inf), ...
    P.machine.bladeCount, 1);
gapComparisonRows = [];
gapComparisonResults = repmat(empty_gap_compare_window_result_local(), 0, 1);

[gapBranch, gapBranchInfo] = prepare_gap_library_branch_local(P, cfg, waveformSourceInfo);

fprintf('\n=== Step06: direct template identification by physical blade ===\n');
fprintf('Input filtered waveform library: %s\n', waveformSourceInfo.filtered_file);
if isfield(waveformSourceInfo, 'note') && ~isempty(waveformSourceInfo.note)
    fprintf('Waveform source note: %s\n', waveformSourceInfo.note);
end
fprintf('Model: V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi)).\n');
fprintf('EO candidates: %.1f-%.1f Hz, VP top-%d, amplitude <= %.3f mm, dx_c <= %.3f mm.\n', ...
    cfg.freq_search_hz(1), cfg.freq_search_hz(2), cfg.vp_top_k_eo, ...
    cfg.amplitude_limit_mm, cfg.dx_c_limit_mm);
fprintf('Model mode: %s\n', string(gapBranchInfo.modelMode));
if gapBranch.enabled
    fprintf('Gap-library branch: sensors %s, corrected library %s\n', ...
        mat2str(gapBranch.sensor_ids), gapBranch.corrected_lib_file);
else
    fprintf('Gap-library branch: disabled (%s)\n', gapBranchInfo.reason);
end
print_step06_scope_local(cfg);
flush_step06_console_local();

resultIndex = 0;
for bladeId = cfg.blades
    B = FilteredWaveformLibrary.Blade(bladeId);
    if isempty(B.Window)
        warning('B%d has no filtered windows. Skipping.', bladeId);
        continue;
    end
    prevPhaseResult = [];
    prevGapCompare = [];
    for w = 1:numel(B.Window)
        W = B.Window(w);
        if ~use_step06_window_local(W, cfg)
            continue;
        end
        fprintf('Step06 solving B%d W%02d laps %s ...\n', bladeId, W.window_id, mat2str(W.lap_range));
        flush_step06_console_local();
        if isempty(W.Bundle) || W.Bundle.point_count < cfg.min_point_count
            warning('B%d W%d has too few filtered points. Skipping.', bladeId, W.window_id);
            continue;
        end

        bundleCore = attach_template_and_speed_local(W.Bundle, WaveformLibrary, bladeId, W.lap_range, cfg);
        eoCandidates = build_eo_candidates_local(bundleCore.rot_freq_mean_hz, cfg.freq_search_hz, cfg.eo_pad);
        seedTableCore = solve_template_seed_eo_scan_local(bundleCore, eoCandidates);
        eoKeepCore = select_vp_topk_eo_local(seedTableCore, cfg.vp_top_k_eo);
        fitCore = refine_template_waveform_fit_local(bundleCore, seedTableCore, eoKeepCore, cfg);
        [phaseRef, phaseRefInfo] = build_phase_safe_reference_local( ...
            prevPhaseResult, fitCore, seedTableCore, W.window_id, cfg);
        [bundleExpanded, expandInfo] = expand_bundle_by_phase_safe_local( ...
            WaveformLibrary, bladeId, W.lap_range, bundleCore, phaseRef, phaseRefInfo, cfg);
        fitExpanded = [];
        seedTableExpanded = [];
        eoKeepExpanded = [];
        if expandInfo.applied
            seedTableExpanded = solve_template_seed_eo_scan_local(bundleExpanded, eoCandidates);
            eoKeepExpanded = select_vp_topk_eo_local(seedTableExpanded, cfg.vp_top_k_eo);
            fitExpanded = refine_template_waveform_fit_local(bundleExpanded, seedTableExpanded, eoKeepExpanded, cfg);
        end
        [bundleChosen, fitChosen, seedTableChosen, eoKeepChosen, expandInfo] = ...
            choose_better_identification_pass_local( ...
            bundleCore, fitCore, seedTableCore, eoKeepCore, ...
            bundleExpanded, fitExpanded, seedTableExpanded, eoKeepExpanded, expandInfo, cfg);

        resultIndex = resultIndex + 1;
        windowResults(resultIndex) = empty_window_result_local(); %#ok<SAGROW>
        windowResults(resultIndex).blade_id = bladeId;
        windowResults(resultIndex).window_id = W.window_id;
        windowResults(resultIndex).lap_range = W.lap_range;
        windowResults(resultIndex).Bundle = bundleChosen;
        windowResults(resultIndex).CoreBundle = bundleCore;
        windowResults(resultIndex).ExpandedBundle = bundleExpanded;
        windowResults(resultIndex).SeedTable = struct2table(rmfield(seedTableChosen, 'sensor_eta'));
        windowResults(resultIndex).SelectedEO = eoKeepChosen(:).';
        windowResults(resultIndex).Result = fitChosen;
        windowResults(resultIndex).CoreResult = fitCore;
        windowResults(resultIndex).ExpandedResult = fitExpanded;
        windowResults(resultIndex).ExpansionInfo = expandInfo;
        windowResults(resultIndex).PhaseSafeReferenceInfo = phaseRefInfo;
        windowResults(resultIndex).ChosenPass = expandInfo.chosen_pass;
        windowResults(resultIndex).SensorEtaInitMM = get_optional_field_local(bundleChosen, 'eta_initial_guess_by_sensor', []);
        windowResults(resultIndex).SensorEtaInitIQRMM = get_optional_field_local(bundleChosen, 'eta_initial_iqr_by_sensor', []);
        windowResults(resultIndex).GapLibraryCompare = [];

        trendRows = [trendRows; struct( ... %#ok<AGROW>
            'BladeID', bladeId, ...
            'WindowID', W.window_id, ...
            'LapStart', W.lap_range(1), ...
            'LapEnd', W.lap_range(end), ...
            'WindowCenterTimeSec', median(bundleChosen.T, 'omitnan'), ...
            'RotFreqHz', bundleChosen.rot_freq_mean_hz, ...
            'RotRPM', bundleChosen.rot_rpm_mean, ...
            'EO', fitChosen.EO_id, ...
            'FrequencyHz', fitChosen.fn_id, ...
            'AmplitudeMM', fitChosen.A_id, ...
            'PhiRad', fitChosen.phi_id_wrapped, ...
            'DxCMM', fitChosen.dx_c_id, ...
            'SensorEtaInitMaxAbsMM', max(abs(get_optional_field_local(bundleChosen, 'eta_initial_guess_by_sensor', 0)), [], 'omitnan'), ...
            'SensorEtaMaxAbsMM', fitChosen.sensor_eta_max_abs_mm, ...
            'CoreWeightedVoltageRMSE', fitCore.weighted_voltage_rmse, ...
            'ExpandedWeightedVoltageRMSE', get_refined_rmse_local(fitExpanded), ...
            'WeightedVoltageRMSE', fitChosen.weighted_voltage_rmse, ...
            'PlainVoltageRMSE', fitChosen.plain_voltage_rmse, ...
            'CorePointCount', fitCore.point_count, ...
            'PointCount', fitChosen.point_count, ...
            'ExpandedPointCount', expandInfo.final_point_count, ...
            'PhaseSafeApplied', expandInfo.applied, ...
            'PhaseSafePassCount', expandInfo.pass_count, ...
            'ChosenPass', expandInfo.chosen_pass, ...
            'PhaseSafeReferenceSource', string(phaseRefInfo.source), ...
            'PhaseSafeReferenceReason', string(phaseRefInfo.reason), ...
            'ValidSegmentCount', fitChosen.valid_segment_count, ...
            'ClampFraction', fitChosen.Coverage.clamp_fraction, ...
            'MaxQueryOvershootMM', fitChosen.Coverage.max_query_overshoot_mm, ...
            'TopKEO', mat2str(eoKeepChosen(:).'), ...
            'Status', expandInfo.status)];

        if gapBranch.enabled && any(strcmpi(gapBranch.model_mode, {'direct_and_gap', 'gap_only'}))
            gapCompare = run_gap_library_compare_local( ...
                bundleChosen, bladeId, W.window_id, W.lap_range, fitChosen, ...
                eoKeepChosen, gapBranch, prevGapCompare, WaveformLibrary);
            windowResults(resultIndex).GapLibraryCompare = gapCompare;
            gapComparisonResults(end + 1) = gapCompare; %#ok<AGROW>
            gapComparisonRows = [gapComparisonRows; gapCompare.trendRow]; %#ok<AGROW>
            prevGapCompare = gapCompare;
        end

        if fitChosen.weighted_voltage_rmse < bestByBlade(bladeId).weighted_voltage_rmse
            bestByBlade(bladeId).blade_id = bladeId;
            bestByBlade(bladeId).weighted_voltage_rmse = fitChosen.weighted_voltage_rmse;
            bestByBlade(bladeId).window_id = W.window_id;
            bestByBlade(bladeId).lap_range = W.lap_range;
            bestByBlade(bladeId).result = fitChosen;
        end
        prevPhaseResult = fitChosen;
        prevPhaseResult.window_id = W.window_id;
        prevPhaseResult.lap_range = W.lap_range;

        fprintf(['B%d W%02d laps %s: chosen=%s, EO=%d, f=%.3f Hz, A=%.4f mm, ' ...
            'dx_c=%.4f mm, RMSE core=%.5f V, expanded=%.5f V, final=%.5f V, points %d -> %d\n'], ...
            bladeId, W.window_id, mat2str(W.lap_range), expandInfo.chosen_pass, ...
            fitChosen.EO_id, fitChosen.fn_id, fitChosen.A_id, fitChosen.dx_c_id, ...
            fitCore.weighted_voltage_rmse, get_refined_rmse_local(fitExpanded), ...
            fitChosen.weighted_voltage_rmse, fitCore.point_count, fitChosen.point_count);
        flush_step06_console_local();
    end
end

if isempty(trendRows)
    error('No identification result was produced.');
end

IdentificationSummary = struct2table(trendRows);
IdentificationResult = struct();
IdentificationResult.dataset = '20241106';
IdentificationResult.method = method;
IdentificationResult.config = cfg;
IdentificationResult.source_filtered_waveform_library = string(waveformSourceInfo.filtered_file);
IdentificationResult.source_waveform_library = string(waveformSourceInfo.waveform_file);
IdentificationResult.source_waveform_note = string(waveformSourceInfo.note);
IdentificationResult.Trend = IdentificationSummary;
IdentificationResult.WindowResult = windowResults;
IdentificationResult.BestByBlade = bestByBlade;
IdentificationResult.ResonanceSummary = build_resonance_summary_local(IdentificationSummary);
IdentificationResult.GapLibraryBranch = gapBranchInfo;

if ~isempty(gapComparisonRows)
    GapLibraryComparisonSummary = struct2table(gapComparisonRows);
    IdentificationResult.GapLibraryComparisonSummary = GapLibraryComparisonSummary;
    IdentificationResult.GapLibraryWindowResult = gapComparisonResults;
else
    GapLibraryComparisonSummary = table();
    IdentificationResult.GapLibraryComparisonSummary = GapLibraryComparisonSummary;
    IdentificationResult.GapLibraryWindowResult = gapComparisonResults;
end

save(outMat, 'IdentificationResult', '-v7.3');
writetable(IdentificationSummary, outCsv);
if ~isempty(GapLibraryComparisonSummary)
    writetable(GapLibraryComparisonSummary, P.files.gapIdentificationSummary);
end

fprintf('\n=== Step06 identification trend ===\n');
disp(IdentificationSummary(:, {'BladeID','WindowID','ChosenPass','EO','FrequencyHz','AmplitudeMM','DxCMM','WeightedVoltageRMSE'}));
fprintf('Saved identification result: %s\n', outMat);
fprintf('Saved identification summary: %s\n', outCsv);
if ~isempty(GapLibraryComparisonSummary)
    fprintf('Saved gap-library comparison summary: %s\n', P.files.gapIdentificationSummary);
end
flush_step06_console_local();

if P.view.enable
    visualize_identification_trend_local(IdentificationSummary, P);
end

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s:\n  %s', label, pathText);
end
end

function P = apply_step06_local_options_local(P, S06)
P.step06 = struct();
P.step06.debug_enabled = true;
P.step06.windowIds = [];
P.step06.lapRange = [];
P.step06.windowCenterTimeRangeSec = [];
P.step06.outputLabel = '';
P.step06.inlineHighSpeedExtraction = true;
P.step06.autoBuildCurrentRegionNumbering = true;
P.step06.saveAutoBuiltRegionAndNumbering = true;
P.step06.preferPeakCache = true;
P.step06.autoBuildPeakCacheIfMissing = true;
P.step06.chooseBetterPass = true;
P.step06.betterPassRmseToleranceV = 0;
P.files.identificationFigure = fullfile(P.view.figureDir, '06_identification', ...
    'Step06_DirectTemplateIdentification_20241106.png');

P.sensors.analysis = S06.analysisSensors(:).';
P.region.startTimeSec = S06.startTimeSec;
P.region.targetLaps = S06.targetLaps;
P.waveform.windowLaps = S06.windowLaps;
P.waveform.slidingStepLaps = S06.slidingStepLaps;
P.waveform.pulsePadSec = S06.pulsePadSec;
P.waveform.dynamicWindowMode = char(S06.dynamicWindowMode);
P.waveform.windowRule = char(S06.windowRule);
P.identification.blades = S06.blades(:).';
P.identification.freqSearchHz = S06.freqSearchHz(:).';
P.identification.topKEO = S06.topKEO;
P.identification.amplitudeLimitMM = S06.amplitudeLimitMM;
P.identification.dxCLimitMM = S06.dxCLimitMM;
P.identification.sensorEtaLimitMM = S06.sensorEtaLimitMM;
P.identification.sensorEtaAdaptiveLimitEnable = S06.sensorEtaAdaptiveLimitEnable;
P.identification.sensorEtaAdaptiveQuantile = S06.sensorEtaAdaptiveQuantile;
P.identification.sensorEtaAdaptiveSafetyFactor = S06.sensorEtaAdaptiveSafetyFactor;
P.identification.sensorEtaAdaptiveMinMM = S06.sensorEtaAdaptiveMinMM;
P.identification.sensorEtaAdaptiveMaxMM = S06.sensorEtaAdaptiveMaxMM;
P.identification.sensorEtaRegWeightVPerMM = S06.sensorEtaRegWeightVPerMM;
P.identification.overshootPenaltyWeight = S06.overshootPenaltyWeight;
P.identification.pulseSelectionMode = S06.pulseSelectionMode;
P.identification.domainSelectionMode = S06.domainSelectionMode;
P.identification.domainMarginMM = S06.domainMarginMM;
P.identification.domainSoftMarginMM = S06.domainSoftMarginMM;
P.identification.queryGuardMode = S06.queryGuardMode;
P.identification.queryGuardMM = S06.queryGuardMM;
P.identification.queryGuardQuantile = S06.queryGuardQuantile;
P.identification.queryGuardSafetyMM = S06.queryGuardSafetyMM;
P.identification.queryGuardMinMM = S06.queryGuardMinMM;
P.identification.queryGuardMaxMM = S06.queryGuardMaxMM;
P.identification.sensorDomainExpandMMTable = S06.sensorDomainExpandMMTable;
P.identification.sensorQueryGuardScaleTable = S06.sensorQueryGuardScaleTable;
P.identification.dynamicEffectiveMode = S06.dynamicEffectiveMode;
P.identification.dynamicTemplateGradientMinRatio = S06.dynamicTemplateGradientMinRatio;
P.identification.dynamicTimeGradientMinRatio = S06.dynamicTimeGradientMinRatio;
P.identification.dynamicPeakQuantile = S06.dynamicPeakQuantile;
P.identification.defaultSensorThreshold = S06.defaultSensorThreshold;
P.identification.weightFloor = S06.weightFloor;
P.identification.phaseSafeExpansionEnable = S06.phaseSafeExpansionEnable;
P.identification.phaseSafeReferenceMode = char(S06.phaseSafeReferenceMode);
P.identification.phaseSafeFallbackMode = char(S06.phaseSafeFallbackMode);
P.identification.phaseSafeRefreshEvery = S06.phaseSafeRefreshEvery;
P.identification.phaseSafeMarginMM = S06.phaseSafeMarginMM;
P.identification.phaseSafeMinPointCount = S06.phaseSafeMinPointCount;
P.identification.phaseSafeMinPointCountPerSensor = S06.phaseSafeMinPointCountPerSensor;
P.identification.phaseSafePrevMaxClampFraction = S06.phaseSafePrevMaxClampFraction;
P.identification.phaseSafePrevMinFinalGapRatio = S06.phaseSafePrevMinFinalGapRatio;
P.step06.windowIds = S06.windowIds(:).';
P.step06.lapRange = S06.lapRange(:).';
P.step06.windowCenterTimeRangeSec = S06.windowCenterTimeRangeSec(:).';
P.step06.outputLabel = char(S06.outputLabel);
P.step06.inlineHighSpeedExtraction = logical(S06.inlineHighSpeedExtraction);
P.step06.autoBuildCurrentRegionNumbering = logical(S06.autoBuildCurrentRegionNumbering);
P.step06.saveAutoBuiltRegionAndNumbering = logical(S06.saveAutoBuiltRegionAndNumbering);
P.step06.preferPeakCache = logical(S06.preferPeakCache);
P.step06.autoBuildPeakCacheIfMissing = logical(S06.autoBuildPeakCacheIfMissing);
P.step06.chooseBetterPass = logical(S06.chooseBetterPass);
P.step06.betterPassRmseToleranceV = S06.betterPassRmseToleranceV;
P.step06.modelMode = char(S06.modelMode);
P.step06.gapLibrarySensors = S06.gapLibrarySensors(:).';
P.step06.gapLibraryTopKEO = round(S06.gapLibraryTopKEO);
P.step06.gapLibraryAmplitudeLimitMM = S06.gapLibraryAmplitudeLimitMM;
P.step06.gapLibraryDxLimitMM = S06.gapLibraryDxLimitMM;
P.step06.gapLibraryDeltaGapLimitMM = S06.gapLibraryDeltaGapLimitMM;
P.step06.gapLibraryStaticRegWeight = S06.gapLibraryStaticRegWeight;
P.step06.gapLibrarySensorEtaFallbackLimitMM = S06.gapLibrarySensorEtaFallbackLimitMM;
P.step06.gapLibrarySensorEtaRegWeightMV = S06.gapLibrarySensorEtaRegWeightMV;
P.step06.gapLibraryEOCandidateMode = char(S06.gapLibraryEOCandidateMode);
P.step06.gapLibraryPrevWindowCandidateMode = char(S06.gapLibraryPrevWindowCandidateMode);
P.step06.gapLibraryPrevWindowRefreshEvery = S06.gapLibraryPrevWindowRefreshEvery;
P.step06.gapLibraryBundleSourceMode = char(S06.gapLibraryBundleSourceMode);
P.step06.gapLibraryModelNames = S06.gapLibraryModelNames;
P.step06.gapLibraryDeltaMuLimit = S06.gapLibraryDeltaMuLimit;
P.step06.gapLibraryCorrectedLibFile = char(S06.gapLibraryCorrectedLibFile);
P.view.enable = logical(S06.viewEnable);

P = refresh_step06_artifact_paths_local(P);
end

function P = refresh_step06_artifact_paths_local(P)
sensorTag = step06_sensor_tag_local(P.sensors.analysis);
timeLabel = step06_time_label_local(P.region.startTimeSec);
suffix = step06_output_suffix_local(P.step06.outputLabel);
bladeTag = step06_blade_tag_local(P.identification.blades);

P.files.lowSpeedNumbering = fullfile(P.outputDir, '02_low_speed_numbering', sensorTag, ...
    'LowSpeedNumbering_20241106.mat');
P.files.regionSelection = fullfile(P.outputDir, '03_region_selection', ...
    sprintf('HighSpeedRegion_%s_20241106.mat', timeLabel));
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel));
P.files.waveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('WaveformLibrary_%s_20241106.mat', timeLabel));
P.files.filteredWaveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('FilteredWaveformLibrary_%s_20241106.mat', timeLabel));
P.files.waveformCompactCase = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('Step06_CompactWaveformCase_%s_%s_20241106.mat', bladeTag, timeLabel));
P.files.filteredWaveformSummary = fullfile(P.outputDir, '05_waveform_library', sensorTag, ...
    sprintf('FilteredWaveformSummary_%s_20241106.csv', timeLabel));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensorTag, ...
    'LowSpeedTemplateLibrary_20241106.mat');
P.files.highSpeedFileRanges = fullfile(P.outputDir, '04_high_speed_numbering', ...
    sprintf('HighSpeedFileRanges_OPR%d_SR%d_20241106.mat', ...
    P.machine.oprChannel, round(P.machine.sampleRateHz)));
P.files.identificationSummary = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationSummary_%s_20241106%s.csv', timeLabel, suffix));
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
P.files.gapIdentificationSummary = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('GapLibraryComparison_%s_20241106%s.csv', timeLabel, suffix));
P.files.gapCompareCacheDir = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('GapLibraryCompareCache_%s_20241106%s', timeLabel, suffix));
P.files.auditReport = fullfile(P.outputDir, '07_audit', sensorTag, ...
    sprintf('NumberingAudit_%s_20241106.csv', timeLabel));
P.files.identificationFigure = fullfile(P.view.figureDir, '06_identification', ...
    sprintf('Step06_DirectTemplateIdentification_20241106%s.png', suffix));
end

function tag = step06_sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = step06_time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
        label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

function suffix = step06_output_suffix_local(label)
suffix = '';
if isempty(label)
    return;
end
label = regexprep(char(label), '[^\w\d-]', '_');
label = regexprep(label, '_+', '_');
label = strtrim(label);
if ~isempty(label)
    suffix = ['_', label];
end
end

function tag = step06_blade_tag_local(bladeIds)
bladeIds = bladeIds(:).';
bladeIds = bladeIds(isfinite(bladeIds));
if isempty(bladeIds)
    tag = 'Ball';
    return;
end
bladeIds = unique(round(bladeIds), 'stable');
tag = ['B', sprintf('%d', bladeIds)];
end

function [WaveformLibrary, FilteredWaveformLibrary, sourceInfo] = load_step06_waveform_inputs_local(P)
sourceInfo = struct('waveform_file', "", 'filtered_file', "", 'note', "");

if isfield(P, 'step06') && isfield(P.step06, 'inlineHighSpeedExtraction') && P.step06.inlineHighSpeedExtraction
    [WaveformLibrary, FilteredWaveformLibrary, sourceInfo] = ...
        build_step06_waveform_inputs_inline_local(P);
    return;
end

error(['Step06 is configured to avoid external waveform-library dependence.\n' ...
    'Set S06.inlineHighSpeedExtraction = true so this script directly extracts\n' ...
    'the high-speed waveform slices needed for identification.']);
end

function [WaveformLibrary, FilteredWaveformLibrary, sourceInfo] = build_step06_waveform_inputs_inline_local(P)
require_file_local(P.files.lowSpeedFingerprint, 'Step01 low-speed fingerprint');
[HighSpeedNumbering, LowSpeedTemplateLibrary, inlineSourceInfo] = load_step06_inline_inputs_local(P);
loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
LowSpeedReference = loadedRef.LowSpeedReference;
requestedBlades = resolve_step06_requested_blades_local(P);
[WaveformLibrary, waveformCaseInfo] = load_or_build_step06_waveform_case_local( ...
    P, HighSpeedNumbering, LowSpeedTemplateLibrary, LowSpeedReference, inlineSourceInfo, requestedBlades);

FilteredWaveformLibrary = build_step06_filtered_waveform_library_inline_local(WaveformLibrary, P);
sourceInfo.waveform_file = waveformCaseInfo.waveform_file;
sourceInfo.filtered_file = "inline_step06_filtered_bundle_in_memory";
sourceInfo.note = sprintf(['%s; numbering source: %s; template source: %s'], ...
    waveformCaseInfo.note, inlineSourceInfo.numbering_note, inlineSourceInfo.template_note);
end

function [WaveformLibrary, info] = load_or_build_step06_waveform_case_local( ...
        P, HighSpeedNumbering, LowSpeedTemplateLibrary, LowSpeedReference, inlineSourceInfo, requestedBlades)
info = struct('waveform_file', "", 'note', "");
[hasCache, WaveformLibrary] = try_load_step06_waveform_case_cache_local(P, inlineSourceInfo, requestedBlades);
if hasCache
    info.waveform_file = string(P.files.waveformCompactCase);
    info.note = 'reused compact waveform case cache';
    return;
end

fileRanges = load_or_build_case_file_ranges_inline_local(P);
rawWindow = build_numbering_time_window_inline_local(HighSpeedNumbering, P);
raw = load_raw_subset_inline_local(P.data.highSpeedDir, fileRanges, P.sensors.analysis, P.machine.sampleRateHz, rawWindow);
oprTimes = load_opr_center_times_inline_local(P);
F_omega_deg = build_phase_speed_inline_local(P, oprTimes);

WaveformLibrary = struct();
WaveformLibrary.mode = 'step06_inline_high_speed_extraction';
WaveformLibrary.analysis_sensors = P.sensors.analysis;
WaveformLibrary.blade_count = P.machine.bladeCount;
WaveformLibrary.region_selection = HighSpeedNumbering.region_selection;
WaveformLibrary.high_speed_numbering_file = inlineSourceInfo.numbering_file;
WaveformLibrary.extraction_window_mode = P.waveform.dynamicWindowMode;
WaveformLibrary.extraction_window_rule = P.waveform.windowRule;
WaveformLibrary.requested_blades = requestedBlades(:).';
WaveformLibrary.note = ['Step06 builds high-speed waveform slices inline from raw probe files. ' ...
    'No Step06 artifact is required for identification reruns.'];
WaveformLibrary.Blade = repmat(struct('blade_id', NaN, 'Sensor', []), P.machine.bladeCount, 1);

for bladeId = requestedBlades
    Sensor = repmat(struct('sensor_id', NaN, 'Lap', []), numel(P.sensors.analysis), 1);
    for is = 1:numel(P.sensors.analysis)
        sid = P.sensors.analysis(is);
        sensorNumbering = HighSpeedNumbering.sensor([HighSpeedNumbering.sensor.sensor_id] == sid);
        probe = load_probe_inline_local(P, sid);
        templateSensor = load_template_sensor_inline_local(LowSpeedTemplateLibrary, bladeId, sid);
        thetaStd = read_opr_center_standard_angle_inline_local(LowSpeedReference, sid, bladeId);
        rows = sensorNumbering.selected_rows_by_physical(:, bladeId);
        Lap = repmat(struct( ...
            'lap_id', NaN, ...
            'row_id', NaN, ...
            't', [], ...
            'V', [], ...
            't_peak', NaN, ...
            't_start', NaN, ...
            't_end', NaN, ...
            'x_abs', [], ...
            'x_rel', [], ...
            'theta', []), numel(rows), 1);
        for k = 1:numel(rows)
            rowId = rows(k);
            tPeak = probe.jilublade(rowId, 3);
            [t0, t1] = build_dynamic_segment_window_inline_local(probe.jilublade, rowId, P);
            keep = raw(sid).T >= t0 & raw(sid).T <= t1;
            tSeg = raw(sid).T(keep);
            vSeg = raw(sid).V(keep);
            idxPrev = find(oprTimes < tPeak, 1, 'last');
            if isempty(idxPrev)
                xAbs = nan(size(tSeg));
                thetaRot = nan(size(tSeg));
            else
                thetaPointsDeg = map_segment_to_relative_angle_inline_local(oprTimes(idxPrev), tSeg, F_omega_deg);
                thetaDiffDeg = wrap_to_signed_period_inline_local(thetaPointsDeg - thetaStd, 360);
                xAbs = thetaDiffDeg * (pi / 180) * P.machine.tipRadiusMM;
                thetaRot = map_time_to_rotor_phase_inline_local(oprTimes, tSeg, P.machine.oprPulsesPerRev);
            end
            xRel = xAbs - templateSensor.xc;
            valid = isfinite(tSeg) & isfinite(vSeg) & isfinite(xRel);
            Lap(k).lap_id = k;
            Lap(k).row_id = rowId;
            Lap(k).t = tSeg(valid);
            Lap(k).V = vSeg(valid);
            Lap(k).t_peak = tPeak;
            Lap(k).t_start = t0;
            Lap(k).t_end = t1;
            Lap(k).x_abs = xAbs(valid);
            Lap(k).x_rel = xRel(valid);
            Lap(k).theta = thetaRot(valid);
        end
        Sensor(is).sensor_id = sid;
        Sensor(is).theta_std_deg = thetaStd;
        Sensor(is).template_x = templateSensor.x_grid(:);
        Sensor(is).template_v = templateSensor.v_grid(:);
        Sensor(is).template_dv_dx = templateSensor.dv_dx(:);
        Sensor(is).template_xc = templateSensor.xc;
        Sensor(is).template_domain = templateSensor.x_domain(:).';
        if isfield(templateSensor, 'eta_app_mm')
            Sensor(is).template_eta_app_mm = templateSensor.eta_app_mm(:);
        else
            Sensor(is).template_eta_app_mm = [];
        end
        if isfield(templateSensor, 'eta_median_mm')
            Sensor(is).template_eta_median_mm = templateSensor.eta_median_mm;
        else
            Sensor(is).template_eta_median_mm = NaN;
        end
        if isfield(templateSensor, 'eta_iqr_mm')
            Sensor(is).template_eta_iqr_mm = templateSensor.eta_iqr_mm;
        else
            Sensor(is).template_eta_iqr_mm = NaN;
        end
        if isfield(templateSensor, 'eta_limit_mm')
            Sensor(is).template_eta_limit_mm = templateSensor.eta_limit_mm;
        else
            Sensor(is).template_eta_limit_mm = NaN;
        end
        Sensor(is).Lap = Lap;
    end
    WaveformLibrary.Blade(bladeId).blade_id = bladeId;
    WaveformLibrary.Blade(bladeId).Sensor = Sensor;
end

save_step06_waveform_case_cache_local(P, WaveformLibrary, inlineSourceInfo, requestedBlades);
info.waveform_file = string(P.files.waveformCompactCase);
info.note = 'built compact waveform case cache from high-speed raw data';
end

function [hasCache, WaveformLibrary] = try_load_step06_waveform_case_cache_local(P, inlineSourceInfo, requestedBlades)
hasCache = false;
WaveformLibrary = struct();
if exist(P.files.waveformCompactCase, 'file') ~= 2
    return;
end
loaded = load(P.files.waveformCompactCase, 'CompactWaveformCase');
if ~isfield(loaded, 'CompactWaveformCase')
    return;
end
cache = loaded.CompactWaveformCase;
if ~isfield(cache, 'WaveformLibrary') || ~isfield(cache, 'meta')
    return;
end
meta = cache.meta;
if ~isequal(meta.requested_blades(:).', requestedBlades(:).')
    return;
end
if ~isequal(meta.analysis_sensors(:).', P.sensors.analysis(:).')
    return;
end
if ~isfinite(meta.pulse_pad_sec) || abs(meta.pulse_pad_sec - P.waveform.pulsePadSec) > eps
    return;
end
if ~strcmpi(char(meta.dynamic_window_mode), char(P.waveform.dynamicWindowMode))
    return;
end
if ~strcmp(char(meta.window_rule), char(P.waveform.windowRule))
    return;
end
if ~strcmp(char(meta.numbering_file), char(inlineSourceInfo.numbering_file))
    return;
end
if ~strcmp(char(meta.template_file), char(inlineSourceInfo.template_file))
    return;
end
if local_file_datenum_local(char(meta.numbering_file)) ~= meta.numbering_datenum
    return;
end
if local_file_datenum_local(char(meta.template_file)) ~= meta.template_datenum
    return;
end
WaveformLibrary = cache.WaveformLibrary;
hasCache = true;
end

function save_step06_waveform_case_cache_local(P, WaveformLibrary, inlineSourceInfo, requestedBlades)
cacheDir = fileparts(P.files.waveformCompactCase);
if exist(cacheDir, 'dir') ~= 7
    mkdir(cacheDir);
end
CompactWaveformCase = struct();
CompactWaveformCase.WaveformLibrary = WaveformLibrary;
CompactWaveformCase.meta = struct( ...
    'requested_blades', requestedBlades(:).', ...
    'analysis_sensors', P.sensors.analysis(:).', ...
    'pulse_pad_sec', P.waveform.pulsePadSec, ...
    'dynamic_window_mode', char(P.waveform.dynamicWindowMode), ...
    'window_rule', char(P.waveform.windowRule), ...
    'numbering_file', char(inlineSourceInfo.numbering_file), ...
    'numbering_datenum', local_file_datenum_local(char(inlineSourceInfo.numbering_file)), ...
    'template_file', char(inlineSourceInfo.template_file), ...
    'template_datenum', local_file_datenum_local(char(inlineSourceInfo.template_file)));
save(P.files.waveformCompactCase, 'CompactWaveformCase', '-v7.3');
end

function [HighSpeedNumbering, LowSpeedTemplateLibrary, sourceInfo] = load_step06_inline_inputs_local(P)
sourceInfo = struct('numbering_file', "", 'template_file', "", 'numbering_note', "", 'template_note', "");

if exist(P.files.highSpeedNumbering, 'file') == 2
    loadedNumbering = load(P.files.highSpeedNumbering, 'HighSpeedNumbering');
    HighSpeedNumbering = loadedNumbering.HighSpeedNumbering;
    sourceInfo.numbering_file = string(P.files.highSpeedNumbering);
    sourceInfo.numbering_note = "exact Step04 numbering artifact";
else
    [hasSubsetArtifact, subsetNumbering, subsetFile] = try_load_subset_high_speed_numbering_inline_local(P);
    if hasSubsetArtifact
        HighSpeedNumbering = subsetNumbering;
        sourceInfo.numbering_file = string(subsetFile);
        sourceInfo.numbering_note = "subset sensors cut from compatible Step04 superset artifact";
    elseif P.step06.autoBuildCurrentRegionNumbering
        require_file_local(P.files.lowSpeedFingerprint, 'Step01 low-speed fingerprint');
        loadedRef = load(P.files.lowSpeedFingerprint, 'LowSpeedReference');
        HighSpeedNumbering = build_step06_high_speed_numbering_inline_local(P, loadedRef.LowSpeedReference);
        if P.step06.saveAutoBuiltRegionAndNumbering
            sourceInfo.numbering_file = string(P.files.highSpeedNumbering);
            sourceInfo.numbering_note = "inline rebuilt Step04-style numbering and auto-saved for current startTimeSec";
        else
            sourceInfo.numbering_file = "inline_step06_high_speed_numbering_in_memory";
            sourceInfo.numbering_note = "inline rebuilt Step04-style numbering in memory for current startTimeSec";
        end
    else
        error('Missing Step04 numbering for sensors %s at %.3f s.', mat2str(P.sensors.analysis), P.region.startTimeSec);
    end
end

if exist(P.files.lowSpeedTemplateLibrary, 'file') == 2
    loadedTemplate = load(P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary');
    LowSpeedTemplateLibrary = loadedTemplate.LowSpeedTemplateLibrary;
    sourceInfo.template_file = string(P.files.lowSpeedTemplateLibrary);
    sourceInfo.template_note = "exact Step05 template artifact";
else
    [LowSpeedTemplateLibrary, sourceFile] = load_compatible_low_speed_template_library_inline_local(P);
    sourceInfo.template_file = string(sourceFile);
    sourceInfo.template_note = "compatible Step05 superset artifact";
end
end

function [tf, HighSpeedNumbering, sourceFile] = try_load_subset_high_speed_numbering_inline_local(P)
tf = false;
HighSpeedNumbering = struct();
sourceFile = '';
try
    [HighSpeedNumbering, sourceFile] = load_subset_high_speed_numbering_inline_local(P);
    tf = true;
catch ME
    if contains(ME.message, 'Missing Step04 numbering')
        tf = false;
        HighSpeedNumbering = struct();
        sourceFile = '';
    else
        rethrow(ME);
    end
end
end

function HighSpeedNumbering = build_step06_high_speed_numbering_inline_local(P, LowSpeedReference)
RegionSelection = build_step06_region_selection_inline_local(P);
validate_step06_numbering_inputs_local(LowSpeedReference, P.sensors.analysis);
fprintf('Step06 inline numbering: startTimeSec = %.6f s, sensors = %s\n', ...
    P.region.startTimeSec, mat2str(P.sensors.analysis));
flush_step06_console_local();

oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
require_file_local(oprFile, 'high-speed OPR pulse file');
loadedOpr = load(oprFile, 'jiluOPR');
oprTimes = loadedOpr.jiluOPR(:, 1);

probeCache = load_step06_probe_cache_for_numbering_local(P.data.highSpeedPulseDir, P.sensors.analysis, oprTimes);
selectedRevIds = select_step06_region_revolutions_local(probeCache, RegionSelection, P);
fprintf('Step06 inline numbering: selected %d revolutions starting from %.6f s\n', ...
    numel(selectedRevIds), RegionSelection.start_time_sec);
flush_step06_console_local();
analysisSensors = P.sensors.analysis;
bladeCount = P.machine.bladeCount;
matchPolyDegree = P.numbering.matchPolyDegree;
pulsePadSec = P.waveform.pulsePadSec;
timeWindow = build_step06_selected_region_time_window_local(probeCache, selectedRevIds, P.waveform.pulsePadSec);
peakCacheSet = [];
if P.step06.preferPeakCache
    peakCacheSet = ensure_step06_peak_cache_set_local( ...
        P, P.sensors.analysis, matchPolyDegree, P.machine.sampleRateHz, P.waveform.pulsePadSec, ...
        P.step06.autoBuildPeakCacheIfMissing);
end

raw = [];
if isempty(peakCacheSet)
    fileRanges = build_case_file_ranges_inline_local(P.data.highSpeedDir, P.machine.oprChannel, P.machine.sampleRateHz);
    raw = load_raw_subset_inline_local(P.data.highSpeedDir, fileRanges, P.sensors.analysis, P.machine.sampleRateHz, timeWindow);
end
sensorResults = repmat(empty_step06_sensor_numbering_result_local(bladeCount), numel(analysisSensors), 1);
summaryRows = repmat(struct( ...
    'SensorID', NaN, ...
    'DominantShift', NaN, ...
    'LocalSlotOfB1', NaN, ...
    'ValidRevolutionCount', NaN, ...
    'MeanBestScore', NaN, ...
    'MinBestScore', NaN, ...
    'MeanScoreMargin', NaN), numel(analysisSensors), 1);

for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    fprintf('Step06 inline numbering: matching CH%d (%d/%d)\n', sid, i, numel(analysisSensors));
    flush_step06_console_local();
    probe = probeCache(i);
    refFingerprint = read_step06_low_speed_fingerprint_local(LowSpeedReference, sid, bladeCount);
    if ~isempty(peakCacheSet)
        [revLib, numberingTable, summaryStruct] = build_step06_high_speed_revolution_library_from_cache_local( ...
            peakCacheSet{i}, probe, selectedRevIds, refFingerprint, bladeCount);
    else
        [revLib, numberingTable, summaryStruct] = build_step06_high_speed_revolution_library_local( ...
            raw(sid), probe, selectedRevIds, refFingerprint, bladeCount, pulsePadSec, matchPolyDegree);
    end

    sensorResults(i).sensor_id = sid;
    sensorResults(i).dominant_shift = summaryStruct.dominant_shift;
    sensorResults(i).local_to_physical_blade_ids = summaryStruct.local_to_physical;
    sensorResults(i).physical_to_local_blade_ids = summaryStruct.physical_to_local;
    sensorResults(i).best_score = summaryStruct.mean_best_score;
    sensorResults(i).score_margin = summaryStruct.mean_score_margin;
    sensorResults(i).selected_revolution_ids = revLib.revolution_ids(:);
    sensorResults(i).selected_rows_by_physical = revLib.selected_rows_by_physical;
    sensorResults(i).revolution_local_rows = revLib.local_rows;
    sensorResults(i).revolution_local_peak_values = revLib.local_peak_values;
    sensorResults(i).revolution_physical_peak_values = revLib.physical_peak_values;
    sensorResults(i).revolution_best_shifts = revLib.best_shifts(:);
    sensorResults(i).revolution_best_scores = revLib.best_scores(:);
    sensorResults(i).revolution_score_margins = revLib.score_margins(:);
    sensorResults(i).revolution_shift_scores = revLib.shift_scores;
    sensorResults(i).revolution_peak_times = revLib.peak_times;
    sensorResults(i).revolution_row_windows = revLib.row_windows;
    sensorResults(i).numbering_table = numberingTable;

    summaryRows(i).SensorID = sid;
    summaryRows(i).DominantShift = summaryStruct.dominant_shift;
    summaryRows(i).LocalSlotOfB1 = summaryStruct.physical_to_local(1);
    summaryRows(i).ValidRevolutionCount = numel(revLib.revolution_ids);
    summaryRows(i).MeanBestScore = summaryStruct.mean_best_score;
    summaryRows(i).MinBestScore = summaryStruct.min_best_score;
    summaryRows(i).MeanScoreMargin = summaryStruct.mean_score_margin;
end

HighSpeedNumbering = struct();
HighSpeedNumbering.mode = 'step06_inline_opr_defined_revolutions_with_six_peak_cyclic_matching';
HighSpeedNumbering.analysis_sensors = analysisSensors;
HighSpeedNumbering.blade_count = bladeCount;
HighSpeedNumbering.region_selection = RegionSelection;
HighSpeedNumbering.selected_revolution_ids = selectedRevIds(:);
HighSpeedNumbering.time_window = timeWindow(:).';
HighSpeedNumbering.low_speed_reference_file = P.files.lowSpeedFingerprint;
HighSpeedNumbering.match_poly_degree = matchPolyDegree;
HighSpeedNumbering.summary = struct2table(summaryRows);
HighSpeedNumbering.sensor = sensorResults;
HighSpeedNumbering.note = ['Step06 auto-built the Step04-style numbering for the current startTimeSec. ' ...
    'Each sensor contains a per-revolution six-peak table in physical blade order.'];

if P.step06.saveAutoBuiltRegionAndNumbering
    save_step06_region_and_numbering_artifacts_local(P, RegionSelection, HighSpeedNumbering);
    fprintf('Step06 inline numbering: saved numbering artifact for current startTimeSec\n');
    flush_step06_console_local();
end
end

function RegionSelection = build_step06_region_selection_inline_local(P)
RegionSelection = struct( ...
    'selection_source', 'step06_inline_manual_start', ...
    'region_id', NaN, ...
    'start_mode', 'manual', ...
    'analysis_sensors', P.sensors.analysis, ...
    'target_laps', P.region.targetLaps, ...
    'start_time_sec', P.region.startTimeSec, ...
    'region_start_sec', P.region.startTimeSec, ...
    'end_time_sec', NaN, ...
    'peak_time_sec', NaN, ...
    'dominant_order', NaN, ...
    'dominant_freq_hz', NaN, ...
    'region_plan_file', P.region.planFile);
end

function save_step06_region_and_numbering_artifacts_local(P, RegionSelection, HighSpeedNumbering)
numberingDir = fileparts(P.files.highSpeedNumbering);
if exist(numberingDir, 'dir') ~= 7
    mkdir(numberingDir);
end
save(P.files.highSpeedNumbering, 'HighSpeedNumbering', '-v7.3');
writetable(HighSpeedNumbering.summary, strrep(P.files.highSpeedNumbering, '.mat', '.csv'));
end

function validate_step06_numbering_inputs_local(LowSpeedReference, analysisSensors)
if isfield(LowSpeedReference, 'analysis_sensors')
    available = unique(LowSpeedReference.analysis_sensors(:).', 'stable');
    requested = unique(analysisSensors(:).', 'stable');
    if ~all(ismember(requested, available))
        error(['Step06 analysisSensors %s are not covered by Step01 artifact sensors %s.'], ...
            mat2str(analysisSensors), mat2str(LowSpeedReference.analysis_sensors));
    end
end
end

function probeCache = load_step06_probe_cache_for_numbering_local(highSpeedPulseDir, sensorIds, oprTimes)
probeCache = repmat(struct('sensor_id', NaN, 'jilublade', [], 'revolution_index', []), numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    probeFile = fullfile(highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
    require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
    loaded = load(probeFile, 'jilublade');
    probeCache(i).sensor_id = sid;
    probeCache(i).jilublade = loaded.jilublade;
    probeCache(i).revolution_index = map_step06_pulses_to_revolutions_local(loaded.jilublade(:, 3), oprTimes);
end
end

function revIndex = map_step06_pulses_to_revolutions_local(pulseTimes, oprTimes)
revIndex = nan(size(pulseTimes));
for i = 1:numel(pulseTimes)
    idx = find(oprTimes <= pulseTimes(i), 1, 'last');
    if ~isempty(idx)
        revIndex(i) = idx;
    end
end
end

function selectedRevIds = select_step06_region_revolutions_local(probeCache, RegionSelection, P)
anchorSensor = P.sensors.anchorPreference(1);
if ~ismember(anchorSensor, [probeCache.sensor_id])
    anchorSensor = probeCache(1).sensor_id;
end
anchorProbe = probeCache([probeCache.sensor_id] == anchorSensor);
rows = find(anchorProbe.jilublade(:, 3) >= RegionSelection.start_time_sec);
if isempty(rows)
    error('No anchor pulses after %.6f s.', RegionSelection.start_time_sec);
end
revIds = unique(anchorProbe.revolution_index(rows), 'stable');
revIds = revIds(isfinite(revIds));
if numel(revIds) < RegionSelection.target_laps
    error('Only %d revolutions found after %.6f s; need %d.', ...
        numel(revIds), RegionSelection.start_time_sec, RegionSelection.target_laps);
end
selectedRevIds = revIds(1:RegionSelection.target_laps);
end

function timeWindow = build_step06_selected_region_time_window_local(probeCache, selectedRevIds, pulsePadSec)
timeWindow = [inf, -inf];
for i = 1:numel(probeCache)
    probe = probeCache(i);
    keep = ismember(probe.revolution_index, selectedRevIds);
    rows = find(keep);
    if isempty(rows)
        continue;
    end
    timeWindow(1) = min(timeWindow(1), min(probe.jilublade(rows, 1)) - pulsePadSec);
    timeWindow(2) = max(timeWindow(2), max(probe.jilublade(rows, 2)) + pulsePadSec);
end
if ~all(isfinite(timeWindow))
    error('Could not determine selected high-speed region time window.');
end
end

function peakCacheSet = ensure_step06_peak_cache_set_local(P, sensorIds, matchPolyDegree, sampleRateHz, pulsePadSec, autoBuildIfMissing)
peakCacheSet = cell(numel(sensorIds), 1);
fileRanges = [];
oprTimes = [];
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    cacheFile = step06_peak_cache_file_local(P, sid);
    cache = struct();
    if exist(cacheFile, 'file') == 2
        loaded = load(cacheFile, 'HighSpeedPeakCacheSensor');
        if isfield(loaded, 'HighSpeedPeakCacheSensor') && ...
                step06_peak_cache_matches_request_local(loaded.HighSpeedPeakCacheSensor, sid, P, matchPolyDegree, sampleRateHz, pulsePadSec)
            cache = loaded.HighSpeedPeakCacheSensor;
        end
    end
    if isempty(fieldnames(cache))
        if ~autoBuildIfMissing
            peakCacheSet = [];
            return;
        end
        fprintf('Step06 inline numbering: building missing CH%d peak cache ...\n', sid);
        flush_step06_console_local();
        if isempty(oprTimes)
            oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
            require_file_local(oprFile, 'high-speed OPR timing file');
            loadedOpr = load(oprFile, 'jiluOPR');
            oprTimes = loadedOpr.jiluOPR(:, 1);
        end
        if isempty(fileRanges)
            fileRanges = build_case_file_ranges_inline_local(P.data.highSpeedDir, P.machine.oprChannel, sampleRateHz);
        end
        probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
        require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
        loadedProbe = load(probeFile, 'jilublade');
        jilublade = loadedProbe.jilublade;
        revolutionIndex = map_step06_pulses_to_revolutions_local(jilublade(:, 3), oprTimes);
        [tRaw, vRaw] = load_step06_full_raw_channel_local(P.data.highSpeedDir, fileRanges, sid, sampleRateHz);
        cache = build_step06_sensor_peak_cache_local( ...
            sid, jilublade, revolutionIndex, tRaw, vRaw, P, pulsePadSec, matchPolyDegree, sampleRateHz);
        cacheDir = fileparts(cacheFile);
        if exist(cacheDir, 'dir') ~= 7
            mkdir(cacheDir);
        end
        HighSpeedPeakCacheSensor = cache; %#ok<NASGU>
        save(cacheFile, 'HighSpeedPeakCacheSensor', '-v7.3');
    end
    peakCacheSet{i} = cache;
end
end

function cacheFile = step06_peak_cache_file_local(P, sid)
cacheFile = fullfile(P.outputDir, '04P_high_speed_peak_cache', ...
    sprintf('HighSpeedPeakCache_CH%d_20241106.mat', sid));
end

function tf = step06_peak_cache_matches_request_local(cache, sid, P, matchPolyDegree, sampleRateHz, pulsePadSec)
tf = isstruct(cache) && ...
    isfield(cache, 'sensor_id') && isequal(cache.sensor_id, sid) && ...
    isfield(cache, 'dataset') && strcmpi(string(cache.dataset), "20241106") && ...
    isfield(cache, 'high_speed_case') && strcmpi(string(cache.high_speed_case), string(P.data.highSpeedCase)) && ...
    isfield(cache, 'sample_rate_hz') && isequal(cache.sample_rate_hz, sampleRateHz) && ...
    isfield(cache, 'match_poly_degree') && isequal(cache.match_poly_degree, matchPolyDegree) && ...
    isfield(cache, 'pulse_pad_sec') && abs(cache.pulse_pad_sec - pulsePadSec) <= eps;
end

function cache = build_step06_sensor_peak_cache_local(sensorId, jilublade, revolutionIndex, tRaw, vRaw, P, pulsePadSec, matchPolyDegree, sampleRateHz)
nRows = size(jilublade, 1);
cache = struct();
cache.dataset = '20241106';
cache.high_speed_case = P.data.highSpeedCase;
cache.sensor_id = sensorId;
cache.sample_rate_hz = sampleRateHz;
cache.match_poly_degree = matchPolyDegree;
cache.pulse_pad_sec = pulsePadSec;
cache.row_count = nRows;
cache.row_id = (1:nRows).';
cache.t_start = jilublade(:, 1) - pulsePadSec;
cache.t_end = jilublade(:, 2) + pulsePadSec;
cache.t_peak_nominal = jilublade(:, 3);
cache.revolution_index = revolutionIndex(:);
cache.raw_peak_value = nan(nRows, 1);
cache.raw_peak_time = nan(nRows, 1);
cache.fit_peak_value = nan(nRows, 1);
cache.fit_peak_time = nan(nRows, 1);
cache.fit_succeeded = false(nRows, 1);
cache.degree_used = nan(nRows, 1);
cache.sample_count = zeros(nRows, 1);

[~, order] = sort(cache.t_peak_nominal, 'ascend');
sortedRowIds = cache.row_id(order);
sortedStarts = cache.t_start(order);
sortedEnds = cache.t_end(order);
idxStart = 1;
idxEnd = 1;
nRaw = numel(tRaw);
for k = 1:nRows
    rowId = sortedRowIds(k);
    while idxStart <= nRaw && tRaw(idxStart) < sortedStarts(k)
        idxStart = idxStart + 1;
    end
    if idxEnd < idxStart
        idxEnd = idxStart;
    end
    while idxEnd <= nRaw && tRaw(idxEnd) <= sortedEnds(k)
        idxEnd = idxEnd + 1;
    end
    i0 = idxStart;
    i1 = idxEnd - 1;
    if i0 > i1 || i0 > nRaw || i1 < 1
        continue;
    end
    tSeg = tRaw(i0:i1);
    vSeg = vRaw(i0:i1);
    cache.sample_count(rowId) = numel(tSeg);
    if numel(tSeg) < max(8, matchPolyDegree + 2)
        continue;
    end
    fitInfo = fit_step06_polynomial_peak_trace_local(tSeg, vSeg, matchPolyDegree);
    cache.raw_peak_value(rowId) = fitInfo.rawPeakValue;
    cache.raw_peak_time(rowId) = fitInfo.rawPeakTime;
    cache.fit_peak_value(rowId) = fitInfo.peakValue;
    cache.fit_peak_time(rowId) = fitInfo.peakTime;
    cache.fit_succeeded(rowId) = fitInfo.fitSucceeded;
    cache.degree_used(rowId) = fitInfo.degreeUsed;
end
end

function [tRaw, vRaw] = load_step06_full_raw_channel_local(caseDir, fileRanges, sensorId, sampleRateHz)
tRaw = [];
vRaw = [];
for k = 1:numel(fileRanges)
    fileId = fileRanges(k).file_id;
    offset = fileRanges(k).offset;
    [tLocal, vLocal] = load_raw_channel_inline_local(caseDir, sensorId, fileId, sampleRateHz);
    tRaw = [tRaw; tLocal(:) + offset]; %#ok<AGROW>
    vRaw = [vRaw; vLocal(:)]; %#ok<AGROW>
end
end

function [revLib, numberingTable, summaryStruct] = build_step06_high_speed_revolution_library_from_cache_local( ...
        peakCache, probe, selectedRevIds, refFingerprint, bladeCount)
nRev = numel(selectedRevIds);
localRows = nan(nRev, bladeCount);
peakTimes = nan(nRev, bladeCount);
localPeakValues = nan(nRev, bladeCount);
shiftScores = nan(nRev, bladeCount);
bestShifts = nan(nRev, 1);
bestScores = nan(nRev, 1);
scoreMargins = nan(nRev, 1);
rowWindows = nan(nRev, bladeCount, 2);

for k = 1:nRev
    revId = selectedRevIds(k);
    rows = find(probe.revolution_index == revId);
    [~, order] = sort(probe.jilublade(rows, 3), 'ascend');
    rows = rows(order);
    if numel(rows) ~= bladeCount
        continue;
    end
    localRows(k, :) = rows(:).';
    peakTimes(k, :) = probe.jilublade(rows, 3).';
    rowWindows(k, :, 1) = peakCache.t_start(rows);
    rowWindows(k, :, 2) = peakCache.t_end(rows);
    localPeakValues(k, :) = peakCache.fit_peak_value(rows);
    peakVec = normalize_step06_vector_max_local(localPeakValues(k, :));
    if ~all(isfinite(peakVec))
        continue;
    end
    for shift = 0:(bladeCount - 1)
        shifted = circshift(peakVec(:), -shift);
        shiftScores(k, shift + 1) = corr_step06_scalar_local(shifted, refFingerprint(:));
    end
    [bestScores(k), pos] = max(shiftScores(k, :));
    if ~isempty(pos) && isfinite(bestScores(k))
        bestShifts(k) = pos - 1;
        scoreMargins(k) = score_margin_step06_local(shiftScores(k, :));
    end
end

dominantShift = choose_step06_dominant_shift_local(bestShifts, bestScores, scoreMargins);
localToPhysical = mod((1:bladeCount) - 1 - dominantShift, bladeCount) + 1;
physicalToLocal = invert_step06_mapping_local(localToPhysical);
physicalPeakValues = reorder_step06_peak_matrix_local(localPeakValues, physicalToLocal);
selectedRowsByPhysical = reorder_step06_peak_matrix_local(localRows, physicalToLocal);

revRows = struct([]);
for k = 1:nRev
    row = struct();
    row.SensorID = probe.sensor_id;
    row.RevolutionID = selectedRevIds(k);
    row.DominantShift = dominantShift;
    row.BestShift = bestShifts(k);
    row.BestScore = bestScores(k);
    row.ScoreMargin = scoreMargins(k);
    row.LocalSlotOfB1 = physicalToLocal(1);
    for b = 1:bladeCount
        row.(sprintf('B%d', b)) = physicalPeakValues(k, b);
    end
    revRows = [revRows; row]; %#ok<AGROW>
end

summaryStruct = struct();
summaryStruct.dominant_shift = dominantShift;
summaryStruct.local_to_physical = localToPhysical;
summaryStruct.physical_to_local = physicalToLocal;
summaryStruct.mean_best_score = mean(bestScores, 'omitnan');
summaryStruct.min_best_score = min(bestScores, [], 'omitnan');
summaryStruct.mean_score_margin = mean(scoreMargins, 'omitnan');

revLib = struct();
revLib.revolution_ids = selectedRevIds(:).';
revLib.local_rows = localRows;
revLib.peak_times = peakTimes;
revLib.local_peak_values = localPeakValues;
revLib.physical_peak_values = physicalPeakValues;
revLib.best_shifts = bestShifts;
revLib.best_scores = bestScores;
revLib.score_margins = scoreMargins;
revLib.shift_scores = shiftScores;
revLib.selected_rows_by_physical = selectedRowsByPhysical;
revLib.row_windows = rowWindows;
numberingTable = struct2table(revRows);
end

function result = empty_step06_sensor_numbering_result_local(bladeCount)
result = struct( ...
    'sensor_id', NaN, ...
    'dominant_shift', NaN, ...
    'local_to_physical_blade_ids', nan(1, bladeCount), ...
    'physical_to_local_blade_ids', nan(1, bladeCount), ...
    'best_score', NaN, ...
    'score_margin', NaN, ...
    'selected_revolution_ids', [], ...
    'selected_rows_by_physical', [], ...
    'revolution_local_rows', [], ...
    'revolution_local_peak_values', [], ...
    'revolution_physical_peak_values', [], ...
    'revolution_best_shifts', [], ...
    'revolution_best_scores', [], ...
    'revolution_score_margins', [], ...
    'revolution_shift_scores', [], ...
    'revolution_peak_times', [], ...
    'revolution_row_windows', [], ...
    'numbering_table', table());
end

function fp = read_step06_low_speed_fingerprint_local(LowSpeedReference, sid, bladeCount)
if ~isKey(LowSpeedReference.fingerprints, sid)
    error('Low-speed fingerprint missing CH%d.', sid);
end
raw = LowSpeedReference.fingerprints(sid);
if numel(raw) < bladeCount
    error('Low-speed fingerprint for CH%d has only %d values.', sid, numel(raw));
end
fp = normalize_step06_vector_max_local(raw(1:bladeCount));
end

function [revLib, numberingTable, summaryStruct] = build_step06_high_speed_revolution_library_local( ...
        rawSensor, probe, selectedRevIds, refFingerprint, bladeCount, pulsePadSec, matchPolyDegree)
nRev = numel(selectedRevIds);
localRows = nan(nRev, bladeCount);
peakTimes = nan(nRev, bladeCount);
localPeakValues = nan(nRev, bladeCount);
shiftScores = nan(nRev, bladeCount);
bestShifts = nan(nRev, 1);
bestScores = nan(nRev, 1);
scoreMargins = nan(nRev, 1);
rowWindows = nan(nRev, bladeCount, 2);

for k = 1:nRev
    revId = selectedRevIds(k);
    rows = find(probe.revolution_index == revId);
    [~, order] = sort(probe.jilublade(rows, 3), 'ascend');
    rows = rows(order);
    if numel(rows) ~= bladeCount
        continue;
    end
    localRows(k, :) = rows(:).';
    peakTimes(k, :) = probe.jilublade(rows, 3).';
    for slot = 1:bladeCount
        rowId = rows(slot);
        t0 = probe.jilublade(rowId, 1) - pulsePadSec;
        t1 = probe.jilublade(rowId, 2) + pulsePadSec;
        rowWindows(k, slot, :) = [t0, t1];
        keep = rawSensor.T >= t0 & rawSensor.T <= t1;
        if nnz(keep) < max(8, matchPolyDegree + 2)
            continue;
        end
        fitInfo = fit_step06_polynomial_peak_trace_local(rawSensor.T(keep), rawSensor.V(keep), matchPolyDegree);
        if fitInfo.fitSucceeded
            localPeakValues(k, slot) = fitInfo.peakValue;
        end
    end
    peakVec = normalize_step06_vector_max_local(localPeakValues(k, :));
    if ~all(isfinite(peakVec))
        continue;
    end
    for shift = 0:(bladeCount - 1)
        shifted = circshift(peakVec(:), -shift);
        shiftScores(k, shift + 1) = corr_step06_scalar_local(shifted, refFingerprint(:));
    end
    [bestScores(k), pos] = max(shiftScores(k, :));
    if ~isempty(pos) && isfinite(bestScores(k))
        bestShifts(k) = pos - 1;
        scoreMargins(k) = score_margin_step06_local(shiftScores(k, :));
    end
end

dominantShift = choose_step06_dominant_shift_local(bestShifts, bestScores, scoreMargins);
localToPhysical = mod((1:bladeCount) - 1 - dominantShift, bladeCount) + 1;
physicalToLocal = invert_step06_mapping_local(localToPhysical);
physicalPeakValues = reorder_step06_peak_matrix_local(localPeakValues, physicalToLocal);
selectedRowsByPhysical = reorder_step06_peak_matrix_local(localRows, physicalToLocal);

revRows = struct([]);
for k = 1:nRev
    row = struct();
    row.SensorID = probe.sensor_id;
    row.RevolutionID = selectedRevIds(k);
    row.DominantShift = dominantShift;
    row.BestShift = bestShifts(k);
    row.BestScore = bestScores(k);
    row.ScoreMargin = scoreMargins(k);
    row.LocalSlotOfB1 = physicalToLocal(1);
    for b = 1:bladeCount
        row.(sprintf('B%d', b)) = physicalPeakValues(k, b);
    end
    revRows = [revRows; row]; %#ok<AGROW>
end

summaryStruct = struct();
summaryStruct.dominant_shift = dominantShift;
summaryStruct.local_to_physical = localToPhysical;
summaryStruct.physical_to_local = physicalToLocal;
summaryStruct.mean_best_score = mean(bestScores, 'omitnan');
summaryStruct.min_best_score = min(bestScores, [], 'omitnan');
summaryStruct.mean_score_margin = mean(scoreMargins, 'omitnan');

revLib = struct();
revLib.revolution_ids = selectedRevIds(:).';
revLib.local_rows = localRows;
revLib.peak_times = peakTimes;
revLib.local_peak_values = localPeakValues;
revLib.physical_peak_values = physicalPeakValues;
revLib.best_shifts = bestShifts;
revLib.best_scores = bestScores;
revLib.score_margins = scoreMargins;
revLib.shift_scores = shiftScores;
revLib.selected_rows_by_physical = selectedRowsByPhysical;
revLib.row_windows = rowWindows;
numberingTable = struct2table(revRows);
end

function physicalMat = reorder_step06_peak_matrix_local(localMat, physicalToLocal)
physicalMat = nan(size(localMat));
for bladeId = 1:numel(physicalToLocal)
    localSlot = physicalToLocal(bladeId);
    physicalMat(:, bladeId) = localMat(:, localSlot);
end
end

function dominantShift = choose_step06_dominant_shift_local(bestShiftByRev, bestScoreByRev, scoreMarginByRev)
valid = isfinite(bestShiftByRev) & isfinite(bestScoreByRev);
if ~any(valid)
    dominantShift = 0;
    return;
end
shiftValues = bestShiftByRev(valid);
scoreValues = bestScoreByRev(valid);
marginValues = scoreMarginByRev(valid);
uniqueShifts = unique(shiftValues(:).');
rankRows = nan(numel(uniqueShifts), 4);
for i = 1:numel(uniqueShifts)
    keep = shiftValues == uniqueShifts(i);
    rankRows(i, :) = [nnz(keep), mean(scoreValues(keep), 'omitnan'), ...
        mean(marginValues(keep), 'omitnan'), -uniqueShifts(i)];
end
[~, order] = sortrows(rankRows, [-1 -2 -3 -4]);
dominantShift = uniqueShifts(order(1));
end

function fitInfo = fit_step06_polynomial_peak_trace_local(t, v, degree)
fitInfo = struct( ...
    'peakValue', NaN, ...
    'peakTime', NaN, ...
    'rawPeakValue', NaN, ...
    'rawPeakTime', NaN, ...
    't', [], ...
    'v', [], ...
    'tFine', [], ...
    'vFine', [], ...
    'fitSucceeded', false, ...
    'degreeUsed', NaN);
t = t(:);
v = v(:);
valid = isfinite(t) & isfinite(v);
t = t(valid);
v = v(valid);
if numel(t) < max(3, degree + 1) || range(t) <= eps
    return;
end
[t, order] = sort(t);
v = v(order);
[rawPeakValue, rawIdx] = max(v);
fitInfo.rawPeakValue = rawPeakValue;
fitInfo.rawPeakTime = t(rawIdx);
fitInfo.t = t;
fitInfo.v = v;
realDegree = min(degree, numel(t) - 1);
fitInfo.degreeUsed = realDegree;
[coef, ~, mu] = polyfit(t, v, realDegree);
tFine = linspace(min(t), max(t), 200).';
vFine = polyval(coef, tFine, [], mu);
[peakValue, peakIdx] = max(vFine);
fitInfo.peakValue = peakValue;
fitInfo.peakTime = tFine(peakIdx);
fitInfo.tFine = tFine;
fitInfo.vFine = vFine;
fitInfo.fitSucceeded = true;
end

function y = normalize_step06_vector_max_local(x)
x = x(:);
mx = max(x);
if ~isfinite(mx) || mx <= eps
    y = nan(size(x));
else
    y = x / mx;
end
end

function value = corr_step06_scalar_local(a, b)
good = isfinite(a) & isfinite(b);
if nnz(good) < 3
    value = NaN;
    return;
end
a = a(good);
b = b(good);
if std(a) <= eps || std(b) <= eps
    value = NaN;
    return;
end
C = corrcoef(a(:), b(:));
value = C(1, 2);
end

function margin = score_margin_step06_local(scoreRow)
sortedScores = sort(scoreRow, 'descend', 'MissingPlacement', 'last');
if numel(sortedScores) >= 2 && isfinite(sortedScores(1)) && isfinite(sortedScores(2))
    margin = sortedScores(1) - sortedScores(2);
else
    margin = NaN;
end
end

function invMap = invert_step06_mapping_local(map)
invMap = nan(size(map));
for i = 1:numel(map)
    v = map(i);
    if isfinite(v) && v >= 1 && v <= numel(map)
        invMap(v) = i;
    end
end
end

function [HighSpeedNumbering, sourceFile] = load_subset_high_speed_numbering_inline_local(P)
sourceFile = '';
timeLabel = step06_time_label_local(P.region.startTimeSec);
pattern = sprintf('HighSpeedNumbering_%s_20241106.mat', timeLabel);
candidates = dir(fullfile(P.outputDir, '04_high_speed_numbering', 'S*', pattern));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'HighSpeedNumbering');
    if ~isfield(loaded, 'HighSpeedNumbering') || ~isfield(loaded.HighSpeedNumbering, 'sensor')
        continue;
    end
    available = [loaded.HighSpeedNumbering.sensor.sensor_id];
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing Step04 numbering for sensors %s at %.3f s.', mat2str(requested), P.region.startTimeSec);
end
sourceFile = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(sourceFile, 'HighSpeedNumbering');
HighSpeedNumbering = loaded.HighSpeedNumbering;
keep = ismember([HighSpeedNumbering.sensor.sensor_id], requested);
HighSpeedNumbering.sensor = HighSpeedNumbering.sensor(keep);
end

function [LowSpeedTemplateLibrary, sourceFile] = load_compatible_low_speed_template_library_inline_local(P)
sourceFile = '';
candidates = dir(fullfile(P.outputDir, '05A_low_speed_template_library', 'S*', 'LowSpeedTemplateLibrary_20241106.mat'));
requested = unique(P.sensors.analysis(:).', 'stable');
bestIdx = [];
bestCount = inf;
for i = 1:numel(candidates)
    fileNow = fullfile(candidates(i).folder, candidates(i).name);
    loaded = load(fileNow, 'LowSpeedTemplateLibrary');
    if ~isfield(loaded, 'LowSpeedTemplateLibrary')
        continue;
    end
    if isfield(loaded.LowSpeedTemplateLibrary, 'sensor_ids')
        available = unique(loaded.LowSpeedTemplateLibrary.sensor_ids(:).', 'stable');
    elseif isfield(loaded.LowSpeedTemplateLibrary, 'analysis_sensors')
        available = unique(loaded.LowSpeedTemplateLibrary.analysis_sensors(:).', 'stable');
    else
        continue;
    end
    if all(ismember(requested, available)) && numel(available) < bestCount
        bestIdx = i;
        bestCount = numel(available);
    end
end
if isempty(bestIdx)
    error('Missing Step05 template library for sensors %s.', mat2str(requested));
end
sourceFile = fullfile(candidates(bestIdx).folder, candidates(bestIdx).name);
loaded = load(sourceFile, 'LowSpeedTemplateLibrary');
LowSpeedTemplateLibrary = loaded.LowSpeedTemplateLibrary;
end

function probe = load_probe_inline_local(P, sid)
probeFile = fullfile(P.data.highSpeedPulseDir, sprintf('jilublade_probe%d.mat', sid));
require_file_local(probeFile, sprintf('high-speed CH%d blade pulse file', sid));
loaded = load(probeFile, 'jilublade');
probe = struct('sensor_id', sid, 'jilublade', loaded.jilublade);
end

function templateSensor = load_template_sensor_inline_local(LowSpeedTemplateLibrary, bladeId, sid)
entry = LowSpeedTemplateLibrary.entry(bladeId);
if ~entry.exists
    error('Missing low-speed template bundle for B%d.', bladeId);
end
loadedBundle = load(entry.template_file, 'TemplateBundle');
bundle = loadedBundle.TemplateBundle;
row = bundle.SensorFiles(bundle.SensorFiles.sensor_id == sid, :);
if height(row) ~= 1
    error('Template bundle for B%d does not contain CH%d.', bladeId, sid);
end
loadedSensor = load(char(row.template_file(1)));
if isfield(loadedSensor, 'Template')
    template = loadedSensor.Template;
elseif isfield(loadedSensor, 'TemplateSingle')
    template = loadedSensor.TemplateSingle;
elseif isfield(loadedSensor, 'sensor_template')
    template = loadedSensor.sensor_template;
else
    error('Unsupported template sensor file: %s', char(row.template_file(1)));
end
if ~isfield(template, 'Sensor') || isempty(template.Sensor)
    error('Template sensor file has no Sensor field: %s', char(row.template_file(1)));
end
templateSensor = template.Sensor;
if numel(templateSensor) > 1
    templateSensor = templateSensor([templateSensor.sensor_id] == sid);
end
if isempty(templateSensor)
    error('Template sensor file does not contain CH%d: %s', sid, char(row.template_file(1)));
end
end

function thetaStd = read_opr_center_standard_angle_inline_local(LowSpeedReference, sid, bladeId)
if isfield(LowSpeedReference, 'standard_angles_opr_center') && ...
        size(LowSpeedReference.standard_angles_opr_center, 1) >= sid && ...
        size(LowSpeedReference.standard_angles_opr_center, 2) >= bladeId
    thetaStd = LowSpeedReference.standard_angles_opr_center(sid, bladeId);
else
    error('Low-speed reference does not contain OPR-center angle for CH%d B%d.', sid, bladeId);
end
end

function oprTimes = load_opr_center_times_inline_local(P)
oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
require_file_local(oprFile, 'high-speed OPR timing file');
loaded = load(oprFile, 'jiluOPR');
oprTimes = loaded.jiluOPR(:, 1);
end

function F = build_phase_speed_inline_local(P, oprTimes)
omegaFile = fullfile(P.data.highSpeedPulseDir, 'omega.mat');
if exist(omegaFile, 'file') == 2
    omegaVars = whos('-file', omegaFile);
    if ismember('omega', {omegaVars.name})
        loadedOmega = load(omegaFile, 'omega');
        omega = loadedOmega.omega;
        if ~isempty(omega)
            F = griddedInterpolant(omega(:, 1), omega(:, 2), 'linear', 'nearest');
            return;
        end
    end
end
spdT = oprTimes(1:(end - P.machine.oprPulsesPerRev));
spdV = 360 ./ max(oprTimes((P.machine.oprPulsesPerRev + 1):end) - ...
    oprTimes(1:(end - P.machine.oprPulsesPerRev)), eps);
F = griddedInterpolant(spdT, spdV, 'linear', 'nearest');
end

function thetaPoints = map_segment_to_relative_angle_inline_local(tRef, tSeg, F_omega_deg)
if isempty(tSeg)
    thetaPoints = zeros(size(tSeg));
    return;
end
dtFirst = linspace(tRef, tSeg(1), 10);
thetaBase = trapz(dtFirst, F_omega_deg(dtFirst));
wSeg = F_omega_deg(tSeg);
thetaRel = cumtrapz(tSeg, wSeg);
thetaPoints = thetaBase + thetaRel;
end

function thetaRot = map_time_to_rotor_phase_inline_local(oprTimes, sampleTimes, pulsesPerRev)
thetaRot = nan(size(sampleTimes));
if isempty(sampleTimes)
    return;
end
for i = 1:numel(sampleTimes)
    prevIdx = find(oprTimes <= sampleTimes(i), 1, 'last');
    nextIdx = prevIdx + pulsesPerRev;
    if isempty(prevIdx) || nextIdx > numel(oprTimes)
        continue;
    end
    revDt = oprTimes(nextIdx) - oprTimes(prevIdx);
    if revDt <= eps
        continue;
    end
    thetaRot(i) = 2 * pi * (sampleTimes(i) - oprTimes(prevIdx)) / revDt;
end
end

function angle = wrap_to_signed_period_inline_local(angle, period)
angle = mod(angle + period / 2, period) - period / 2;
end

function timeWindow = build_numbering_time_window_inline_local(HighSpeedNumbering, P)
timeWindow = [inf, -inf];
requestedBlades = resolve_step06_requested_blades_local(P);
for is = 1:numel(HighSpeedNumbering.sensor)
    sid = HighSpeedNumbering.sensor(is).sensor_id;
    probe = load_probe_inline_local(P, sid);
    rows = HighSpeedNumbering.sensor(is).selected_rows_by_physical;
    rows = rows(:, requestedBlades);
    rows = rows(isfinite(rows));
    for rowId = rows(:).'
        [t0, t1] = build_dynamic_segment_window_inline_local(probe.jilublade, rowId, P);
        timeWindow(1) = min(timeWindow(1), t0);
        timeWindow(2) = max(timeWindow(2), t1);
    end
end
if ~all(isfinite(timeWindow))
    error('Could not build inline high-speed raw-data window.');
end
end

function [t0, t1] = build_dynamic_segment_window_inline_local(jilublade, row, P)
t0 = jilublade(row, 1) - P.waveform.pulsePadSec;
t1 = jilublade(row, 2) + P.waveform.pulsePadSec;
if ~all(isfinite([t0, t1])) || t1 <= t0
    error('Invalid jilublade row-bound window for row %d: [%.9f, %.9f].', row, t0, t1);
end
end

function FilteredWaveformLibrary = build_step06_filtered_waveform_library_inline_local(WaveformLibrary, P)
method = build_step06_inline_filter_cfg_local(P);
requestedBlades = resolve_step06_requested_blades_local(P);
FilteredWaveformLibrary = struct();
FilteredWaveformLibrary.mode = 'step06_inline_filtered_identification_input_waveforms';
FilteredWaveformLibrary.source_waveform_library = 'inline_memory';
FilteredWaveformLibrary.extraction_window_mode = P.waveform.dynamicWindowMode;
FilteredWaveformLibrary.extraction_window_rule = P.waveform.windowRule;
FilteredWaveformLibrary.note = ['Step06 filters its own inline-extracted high-speed waveforms. ' ...
    'Changing S06.startTimeSec/windowLaps/slidingStepLaps immediately changes identification input.'];
FilteredWaveformLibrary.analysis_sensors = P.sensors.analysis;
FilteredWaveformLibrary.window_laps = P.waveform.windowLaps;
FilteredWaveformLibrary.sliding_step_laps = P.waveform.slidingStepLaps;
FilteredWaveformLibrary.Blade = repmat(struct('blade_id', NaN, 'Window', []), P.machine.bladeCount, 1);

for bladeId = requestedBlades
    B = WaveformLibrary.Blade(bladeId);
    if isempty(B.Sensor)
        continue;
    end
    B = add_core_filter_lap_cache_to_blade_local(B);
    lapCount = min(arrayfun(@(s) numel(s.Lap), B.Sensor));
    numWindows = floor((lapCount - P.waveform.windowLaps) / P.waveform.slidingStepLaps) + 1;
    if numWindows < 1
        continue;
    end
    candidateWindowIds = resolve_step06_requested_window_ids_local(numWindows, P);
    if isempty(candidateWindowIds)
        FilteredWaveformLibrary.Blade(bladeId).blade_id = bladeId;
        FilteredWaveformLibrary.Blade(bladeId).Window = repmat(struct( ...
            'window_id', NaN, 'lap_range', [], 'Bundle', [], 'SensorMeta', []), 0, 1);
        continue;
    end
    Window = repmat(struct('window_id', NaN, 'lap_range', [], 'Bundle', [], 'SensorMeta', []), numel(candidateWindowIds), 1);
    for iw = 1:numel(candidateWindowIds)
        w = candidateWindowIds(iw);
        lapStart = 1 + (w - 1) * P.waveform.slidingStepLaps;
        lapRange = lapStart:(lapStart + P.waveform.windowLaps - 1);
        [bundle, sensorMeta] = build_filtered_bundle_for_blade_window_inline_local(B, lapRange, method);
        Window(iw).window_id = w;
        Window(iw).lap_range = lapRange;
        Window(iw).Bundle = bundle;
        Window(iw).SensorMeta = sensorMeta;
    end
    FilteredWaveformLibrary.Blade(bladeId).blade_id = bladeId;
    FilteredWaveformLibrary.Blade(bladeId).Window = Window;
end
end

function B = add_core_filter_lap_cache_to_blade_local(B)
for is = 1:numel(B.Sensor)
    Tpl = template_from_sensor_library_local(B.Sensor(is));
    for k = 1:numel(B.Sensor(is).Lap)
        xRaw = B.Sensor(is).Lap(k).x_rel(:);
        B.Sensor(is).Lap(k).core_f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw, 'pchip', NaN);
        B.Sensor(is).Lap(k).core_fx = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw, 'pchip', NaN);
    end
end
end

function [bundle, sensorMeta] = build_filtered_bundle_for_blade_window_inline_local(B, lapRange, method)
X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
F0 = [];
Fx = [];
sensorIndex = [];
sensorMeta = repmat(empty_sensor_meta_local(), numel(B.Sensor), 1);
etaInitBySensor = nan(1, numel(B.Sensor));
etaSpreadBySensor = nan(1, numel(B.Sensor));

for is = 1:numel(B.Sensor)
    Sdata = B.Sensor(is);
    Tpl = template_from_sensor_library_local(Sdata);
    [tRaw, xRaw, vRaw, thetaRaw, f0Raw, fxRaw] = concatenate_raw_laps_with_core_cache_local(Sdata.Lap(lapRange));
    [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        select_core_sensor_points_inline_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, Sdata.sensor_id, f0Raw, fxRaw);
    meta.sensor_id = Sdata.sensor_id;
    [etaGuess, etaSpread] = estimate_sensor_eta_initial_guess_local(Tpl, xSel, vSel);
    meta.eta_initial_guess_mm = etaGuess;
    meta.eta_initial_iqr_mm = etaSpread;
    sensorMeta(is) = meta;
    etaInitBySensor(is) = etaGuess;
    etaSpreadBySensor(is) = etaSpread;

    X = [X; xSel(:)]; %#ok<AGROW>
    T = [T; tSel(:)]; %#ok<AGROW>
    V = [V; vSel(:)]; %#ok<AGROW>
    S = [S; repmat(Sdata.sensor_id, numel(tSel), 1)]; %#ok<AGROW>
    W = [W; wSel(:)]; %#ok<AGROW>
    Theta = [Theta; thetaSel(:)]; %#ok<AGROW>
    F0 = [F0; f0Sel(:)]; %#ok<AGROW>
    Fx = [Fx; fxSel(:)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, numel(tSel), 1)]; %#ok<AGROW>
end

if isempty(T)
    error('No valid inline high-speed points were built for B%d laps %s.', B.blade_id, mat2str(lapRange));
end

bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.sensor_index = sensorIndex;
bundle.sensor_ids = [B.Sensor.sensor_id];
bundle.window_meta = sensorMeta;
bundle.WaveformSensor = B.Sensor;
bundle.valid_segment_count = sum([sensorMeta.pulse_segment_count]);
bundle.point_count = numel(T);
bundle.time_window = [min(T), max(T)];
bundle.selection_pass = 'step06_inline_core_filtering';
bundle.eta_initial_guess_by_sensor = anchor_sensor_eta_vector_local(etaInitBySensor);
bundle.eta_initial_iqr_by_sensor = etaSpreadBySensor;
end

function [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        select_core_sensor_points_inline_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, method, sensorId, f0, fx0)
meta = empty_sensor_meta_local();
meta.raw_points = numel(tRaw);
meta.domain_selection_mode = method.domain_selection_mode;
meta.domain_soft_margin_mm = method.domain_soft_margin_mm;
xDomainUse = expand_sensor_domain_local(Tpl.x_domain, method, sensorId);
meta.template_domain = xDomainUse;

if nargin < 8 || isempty(f0)
    f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw(:), 'pchip', NaN);
end
if nargin < 9 || isempty(fx0)
    fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw(:), 'pchip', NaN);
end
threshold = method.default_sensor_threshold;
if isfield(Tpl, 'threshold') && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
end

xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
if strcmpi(method.pulse_selection_mode, 'all')
    [maskPulse, pulseSegmentCount] = isolate_all_pulses_local(vRaw, threshold);
else
    [maskPulse, pulseSegmentCount] = isolate_main_pulse_local(vRaw, threshold);
end
maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, method);
maskDomain = xRaw >= xDomainUse(1) + method.domain_margin_mm & ...
             xRaw <= xDomainUse(2) - method.domain_margin_mm;
maskFinite = isfinite(f0) & isfinite(fx0) & isfinite(vRaw) & isfinite(thetaRaw);
maskBase = maskPulse & maskEffective & maskDomain & maskFinite;
queryGuardMM = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, maskBase, method, sensorId);
maskQuerySafe = xRaw >= xDomain(1) + queryGuardMM & xRaw <= xDomain(2) - queryGuardMM;
mask = maskBase;
if strcmpi(method.domain_selection_mode, 'hard')
    safeMask = mask & maskQuerySafe;
    if nnz(safeMask) >= 8
        mask = safeMask;
    end
end
if nnz(mask) < 8
    mask = maskEffective & maskDomain & maskFinite;
    if strcmpi(method.domain_selection_mode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
end

tSel = tRaw(mask);
xSel = xRaw(mask);
vSel = vRaw(mask);
thetaSel = thetaRaw(mask);
f0Sel = f0(mask);
fxSel = fx0(mask);
wSel = build_filtered_weight_local_refine(tSel, xSel, vSel, Tpl, queryGuardMM, method);

if isempty(tSel)
    meta.valid_points = 0;
else
    vStatic = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xSel(:), 'pchip', NaN);
    meta.valid_points = numel(tSel);
    meta.selected_x_range = [min(xSel), max(xSel)];
    meta.shape_rmse = sqrt(mean((vSel(:) - vStatic(:)).^2, 'omitnan'));
end
meta.pulse_segment_count = pulseSegmentCount;
meta.dynamic_effective_points = nnz(maskEffective);
meta.query_safe_points = nnz(maskQuerySafe(mask));
meta.query_guard_mm = queryGuardMM;
end

function method = build_step06_inline_filter_cfg_local(P)
method = struct();
method.pulse_selection_mode = lower(strtrim(P.identification.pulseSelectionMode));
method.domain_selection_mode = lower(strtrim(P.identification.domainSelectionMode));
method.domain_margin_mm = P.identification.domainMarginMM;
method.domain_soft_margin_mm = P.identification.domainSoftMarginMM;
method.query_guard_mode = lower(strtrim(P.identification.queryGuardMode));
method.query_guard_mm = P.identification.queryGuardMM;
method.query_guard_quantile = P.identification.queryGuardQuantile;
method.query_guard_safety_mm = P.identification.queryGuardSafetyMM;
method.query_guard_min_mm = P.identification.queryGuardMinMM;
method.query_guard_max_mm = P.identification.queryGuardMaxMM;
method.sensor_domain_expand_mm_table = P.identification.sensorDomainExpandMMTable;
method.sensor_query_guard_scale_table = P.identification.sensorQueryGuardScaleTable;
method.dynamic_effective_mode = lower(strtrim(P.identification.dynamicEffectiveMode));
method.dynamic_template_gradient_min_ratio = P.identification.dynamicTemplateGradientMinRatio;
method.dynamic_time_gradient_min_ratio = P.identification.dynamicTimeGradientMinRatio;
method.dynamic_peak_quantile = P.identification.dynamicPeakQuantile;
method.default_sensor_threshold = P.identification.defaultSensorThreshold;
method.weight_floor = P.identification.weightFloor;
end

function fileRanges = load_or_build_case_file_ranges_inline_local(P)
currentFileIds = list_case_file_ids_inline_local(P.data.highSpeedDir, P.machine.oprChannel);
if exist(P.files.highSpeedFileRanges, 'file') == 2
    loaded = load(P.files.highSpeedFileRanges, 'HighSpeedFileRangesCache');
    if isfield(loaded, 'HighSpeedFileRangesCache')
        cache = loaded.HighSpeedFileRangesCache;
        if isfield(cache, 'case_dir') && strcmp(cache.case_dir, P.data.highSpeedDir) && ...
                isfield(cache, 'opr_channel') && cache.opr_channel == P.machine.oprChannel && ...
                isfield(cache, 'sample_rate_hz') && cache.sample_rate_hz == P.machine.sampleRateHz && ...
                isfield(cache, 'file_ids') && isequal(cache.file_ids(:), currentFileIds(:)) && ...
                isfield(cache, 'fileRanges')
            fileRanges = cache.fileRanges;
            return;
        end
    end
end

fileRanges = build_case_file_ranges_inline_local(P.data.highSpeedDir, P.machine.oprChannel, P.machine.sampleRateHz, currentFileIds);
cacheDir = fileparts(P.files.highSpeedFileRanges);
if exist(cacheDir, 'dir') ~= 7
    mkdir(cacheDir);
end
HighSpeedFileRangesCache = struct( ...
    'case_dir', P.data.highSpeedDir, ...
    'opr_channel', P.machine.oprChannel, ...
    'sample_rate_hz', P.machine.sampleRateHz, ...
    'file_ids', currentFileIds(:), ...
    'fileRanges', fileRanges);
save(P.files.highSpeedFileRanges, 'HighSpeedFileRangesCache', '-v7.3');
end

function fileRanges = build_case_file_ranges_inline_local(caseDir, oprChannel, sampleRateHz, fileIds)
if nargin < 4 || isempty(fileIds)
    fileIds = list_case_file_ids_inline_local(caseDir, oprChannel);
end
if isempty(fileIds)
    error('No OPR raw files found in %s.', caseDir);
end
fileIds = fileIds(:);
fileIds = sort(fileIds);
fileIds = fileIds(isfinite(fileIds));
fileIds = unique(fileIds, 'stable');
fileIds = fileIds(:);
fileIds = fileIds.';
fileRanges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(fileIds), 1);
lastEnd = 0;
for i = 1:numel(fileIds)
    [tOpr, ~] = load_raw_channel_inline_local(caseDir, oprChannel, fileIds(i), sampleRateHz);
    if i == 1
        offset = 0;
    else
        offset = lastEnd + 1 / sampleRateHz - tOpr(1);
    end
    fileRanges(i).file_id = fileIds(i);
    fileRanges(i).offset = offset;
    fileRanges(i).t_start = tOpr(1) + offset;
    fileRanges(i).t_end = tOpr(end) + offset;
    lastEnd = fileRanges(i).t_end;
end
end

function fileIds = list_case_file_ids_inline_local(caseDir, oprChannel)
files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.dat', oprChannel)));
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('Probe%d_Data_*.txt', oprChannel)));
end
if isempty(files)
    files = dir(fullfile(caseDir, sprintf('4-%d-*.mat', oprChannel)));
end
if isempty(files)
    error('No OPR raw files found in %s.', caseDir);
end
fileIds = nan(numel(files), 1);
for i = 1:numel(files)
    token = regexp(files(i).name, 'Data_(\d+)', 'tokens', 'once');
    if isempty(token)
        token = regexp(files(i).name, '4-\d+-(\d+)\.mat', 'tokens', 'once');
    end
    fileIds(i) = str2double(token{1});
end
[fileIds, order] = sort(fileIds);
files = files(order); %#ok<NASGU>
fileIds = fileIds(:);
end

function raw = load_raw_subset_inline_local(caseDir, fileRanges, sensorIds, sampleRateHz, timeWindow)
raw(max(sensorIds)) = struct('T', [], 'V', []);
useFiles = fileRanges([fileRanges.t_end] >= timeWindow(1) & [fileRanges.t_start] <= timeWindow(2));
for k = 1:numel(useFiles)
    fileId = useFiles(k).file_id;
    offset = useFiles(k).offset;
    for sid = sensorIds
        [tLocal, vLocal] = load_raw_channel_inline_local(caseDir, sid, fileId, sampleRateHz);
        tGlobal = tLocal(:) + offset;
        keep = tGlobal >= timeWindow(1) & tGlobal <= timeWindow(2);
        raw(sid).T = [raw(sid).T; tGlobal(keep)]; %#ok<AGROW>
        raw(sid).V = [raw(sid).V; vLocal(keep)]; %#ok<AGROW>
    end
end
end

function [tSec, v] = load_raw_channel_inline_local(caseDir, channelId, fileId, sampleRateHz)
patterns = { ...
    sprintf('Probe%d_Data_%d.dat', channelId, fileId), ...
    sprintf('Probe%d_Data_%d.txt', channelId, fileId), ...
    sprintf('4-%d-%d.mat', channelId, fileId)};
filePath = '';
for i = 1:numel(patterns)
    candidate = fullfile(caseDir, patterns{i});
    if exist(candidate, 'file') == 2
        filePath = candidate;
        break;
    end
end
if isempty(filePath)
    error('Missing raw file for CH%d file %d in %s.', channelId, fileId, caseDir);
end
varName = sprintf('jilu%02d', channelId);
vars = whos('-file', filePath);
if ~ismember(varName, {vars.name})
    error('Raw MAT file does not contain %s: %s', varName, filePath);
end
loaded = load(filePath, varName);
raw = loaded.(varName);
raw(raw(:, 1) == 0, :) = [];
tSec = raw(:, 1) / sampleRateHz;
v = raw(:, 2);
end

function requestedBlades = resolve_step06_requested_blades_local(P)
requestedBlades = P.identification.blades(:).';
requestedBlades = requestedBlades(isfinite(requestedBlades));
requestedBlades = unique(round(requestedBlades), 'stable');
requestedBlades = requestedBlades(requestedBlades >= 1 & requestedBlades <= P.machine.bladeCount);
if isempty(requestedBlades)
    requestedBlades = 1:P.machine.bladeCount;
end
end

function dt = local_file_datenum_local(pathText)
dt = NaN;
if isempty(pathText)
    return;
end
if startsWith(pathText, 'inline_', 'IgnoreCase', true)
    return;
end
if exist(pathText, 'file') ~= 2
    return;
end
info = dir(pathText);
if isempty(info)
    return;
end
dt = info.datenum;
end

function candidateWindowIds = resolve_step06_requested_window_ids_local(numWindows, P)
candidateWindowIds = 1:numWindows;
if isfield(P, 'step06') && ~isempty(P.step06.windowIds)
    requestedWindowIds = unique(round(P.step06.windowIds(:).'), 'stable');
    candidateWindowIds = intersect(candidateWindowIds, requestedWindowIds, 'stable');
end
if isfield(P, 'step06') && ~isempty(P.step06.lapRange)
    keep = false(size(candidateWindowIds));
    for i = 1:numel(candidateWindowIds)
        w = candidateWindowIds(i);
        lapStart = 1 + (w - 1) * P.waveform.slidingStepLaps;
        lapRange = lapStart:(lapStart + P.waveform.windowLaps - 1);
        keep(i) = min(lapRange) >= P.step06.lapRange(1) && max(lapRange) <= P.step06.lapRange(end);
    end
    candidateWindowIds = candidateWindowIds(keep);
end
end

function cfg = build_identification_config_local(P)
cfg = struct();
cfg.blades = P.identification.blades;
cfg.analysis_sensors = P.sensors.analysis;
cfg.freq_search_hz = P.identification.freqSearchHz;
cfg.eo_pad = P.identification.eoPad;
cfg.vp_top_k_eo = P.identification.topKEO;
cfg.amplitude_limit_mm = P.identification.amplitudeLimitMM;
cfg.dx_c_limit_mm = P.identification.dxCLimitMM;
cfg.sensor_eta_limit_mm = P.identification.sensorEtaLimitMM;
cfg.sensor_eta_adaptive_limit_enable = P.identification.sensorEtaAdaptiveLimitEnable;
cfg.sensor_eta_adaptive_quantile = P.identification.sensorEtaAdaptiveQuantile;
cfg.sensor_eta_adaptive_safety_factor = P.identification.sensorEtaAdaptiveSafetyFactor;
cfg.sensor_eta_adaptive_min_mm = P.identification.sensorEtaAdaptiveMinMM;
cfg.sensor_eta_adaptive_max_mm = P.identification.sensorEtaAdaptiveMaxMM;
cfg.sensor_eta_reg_weight_v_per_mm = P.identification.sensorEtaRegWeightVPerMM;
cfg.overshoot_penalty_weight = P.identification.overshootPenaltyWeight;
cfg.pulse_selection_mode = lower(strtrim(P.identification.pulseSelectionMode));
cfg.domain_selection_mode = lower(strtrim(P.identification.domainSelectionMode));
cfg.domain_margin_mm = P.identification.domainMarginMM;
cfg.domain_soft_margin_mm = P.identification.domainSoftMarginMM;
cfg.query_guard_mode = lower(strtrim(P.identification.queryGuardMode));
cfg.query_guard_mm = P.identification.queryGuardMM;
cfg.query_guard_quantile = P.identification.queryGuardQuantile;
cfg.query_guard_safety_mm = P.identification.queryGuardSafetyMM;
cfg.query_guard_min_mm = P.identification.queryGuardMinMM;
cfg.query_guard_max_mm = P.identification.queryGuardMaxMM;
cfg.sensor_domain_expand_mm_table = P.identification.sensorDomainExpandMMTable;
cfg.sensor_query_guard_scale_table = P.identification.sensorQueryGuardScaleTable;
cfg.dynamic_effective_mode = lower(strtrim(P.identification.dynamicEffectiveMode));
cfg.dynamic_template_gradient_min_ratio = P.identification.dynamicTemplateGradientMinRatio;
cfg.dynamic_time_gradient_min_ratio = P.identification.dynamicTimeGradientMinRatio;
cfg.dynamic_peak_quantile = P.identification.dynamicPeakQuantile;
cfg.default_sensor_threshold = P.identification.defaultSensorThreshold;
cfg.weight_floor = P.identification.weightFloor;
cfg.phase_safe_expansion_enable = P.identification.phaseSafeExpansionEnable;
cfg.phase_safe_reference_mode = lower(strtrim(P.identification.phaseSafeReferenceMode));
cfg.phase_safe_fallback_mode = lower(strtrim(P.identification.phaseSafeFallbackMode));
cfg.phase_safe_refresh_every = P.identification.phaseSafeRefreshEvery;
cfg.phase_safe_margin_mm = P.identification.phaseSafeMarginMM;
cfg.phase_safe_min_point_count = P.identification.phaseSafeMinPointCount;
cfg.phase_safe_min_point_count_per_sensor = P.identification.phaseSafeMinPointCountPerSensor;
cfg.phase_safe_prev_max_clamp_fraction = P.identification.phaseSafePrevMaxClampFraction;
cfg.phase_safe_prev_min_final_gap_ratio = P.identification.phaseSafePrevMinFinalGapRatio;
cfg.inline_high_speed_extraction = isfield(P, 'step06') && P.step06.inlineHighSpeedExtraction;
cfg.auto_build_current_region_numbering = isfield(P, 'step06') && P.step06.autoBuildCurrentRegionNumbering;
cfg.prefer_peak_cache = isfield(P, 'step06') && P.step06.preferPeakCache;
cfg.window_laps = P.waveform.windowLaps;
cfg.sliding_step_laps = P.waveform.slidingStepLaps;
cfg.pulse_pad_sec = P.waveform.pulsePadSec;
cfg.min_point_count = 20;
cfg.optim_max_iter = 400;
cfg.optim_max_fun_evals = 1200;
cfg.debug_enabled = isfield(P, 'step06') && P.step06.debug_enabled;
cfg.window_ids = P.step06.windowIds;
cfg.lap_range = P.step06.lapRange;
cfg.window_center_time_range_sec = P.step06.windowCenterTimeRangeSec;
cfg.choose_better_pass = isfield(P, 'step06') && P.step06.chooseBetterPass;
cfg.better_pass_rmse_tolerance_v = P.step06.betterPassRmseToleranceV;
end

function print_step06_scope_local(cfg)
if ~cfg.debug_enabled
    return;
end
fprintf('Step06 local parameter block is ON.\n');
fprintf('  Sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('  Blades: %s\n', mat2str(cfg.blades));
fprintf('  Inline high-speed extraction: %d\n', cfg.inline_high_speed_extraction);
fprintf('  Auto-build current region numbering: %d\n', cfg.auto_build_current_region_numbering);
fprintf('  Prefer peak cache: %d\n', cfg.prefer_peak_cache);
fprintf('  Window laps / step: %d / %d\n', cfg.window_laps, cfg.sliding_step_laps);
if ~isempty(cfg.window_ids)
    fprintf('  Window IDs: %s\n', mat2str(cfg.window_ids));
end
if ~isempty(cfg.lap_range)
    fprintf('  Lap range filter: [%g %g]\n', cfg.lap_range(1), cfg.lap_range(end));
end
if ~isempty(cfg.window_center_time_range_sec)
    fprintf('  Window center time filter: [%.6f %.6f] s\n', ...
        cfg.window_center_time_range_sec(1), cfg.window_center_time_range_sec(end));
end
flush_step06_console_local();
end

function flush_step06_console_local()
drawnow;
end

function values = parse_step06_int_env_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    values = [];
    return;
end
tokens = regexp(raw, '(-?\d+)\s*:\s*(-?\d+)|-?\d+', 'match');
values = [];
for i = 1:numel(tokens)
    token = strtrim(tokens{i});
    if contains(token, ':')
        parts = sscanf(token, '%d:%d').';
        if numel(parts) == 2
            if parts(1) <= parts(2)
                values = [values, parts(1):parts(2)]; %#ok<AGROW>
            else
                values = [values, parts(1):-1:parts(2)]; %#ok<AGROW>
            end
        end
    else
        v = str2double(token);
        if isfinite(v)
            values = [values, v]; %#ok<AGROW>
        end
    end
end
values = unique(round(values), 'stable');
end

function value = parse_step06_scalar_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value)
    error('%s must be numeric.', name);
end
end

function value = parse_step06_text_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
else
    value = raw;
end
end

function tf = use_step06_window_local(W, cfg)
tf = true;
if ~isempty(cfg.window_ids) && ~ismember(W.window_id, cfg.window_ids)
    tf = false;
    return;
end
if ~isempty(cfg.lap_range)
    lr = W.lap_range;
    tf = min(lr) >= cfg.lap_range(1) && max(lr) <= cfg.lap_range(end);
    if ~tf
        return;
    end
end
if ~isempty(cfg.window_center_time_range_sec)
    tCenter = median(W.Bundle.T, 'omitnan');
    tf = isfinite(tCenter) && ...
        tCenter >= cfg.window_center_time_range_sec(1) && ...
        tCenter <= cfg.window_center_time_range_sec(end);
end
end

function method = build_method_info_local(cfg)
method = struct();
method.name = 'template_only_vp_topk_synchronous_waveform';
method.selection_rule = 'first_order_template_vp_topk_then_prev_window_phase_safe_expansion';
method.model = 'V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi))';
method.vp_seed_model = 'V_minus_Tx_equals_minus_Tprime_times_u';
method.final_waveform_objective = 'fixed_eo_direct_template_voltage_residual';
method.vp_top_k_eo = cfg.vp_top_k_eo;
method.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
method.sensor_eta_adaptive_limit_enable = cfg.sensor_eta_adaptive_limit_enable;
method.overshoot_penalty_weight = cfg.overshoot_penalty_weight;
method.phase_safe_expansion_enable = cfg.phase_safe_expansion_enable;
method.phase_safe_reference_mode = cfg.phase_safe_reference_mode;
method.phase_safe_fallback_mode = cfg.phase_safe_fallback_mode;
method.phase_safe_refresh_every = cfg.phase_safe_refresh_every;
method.phase_safe_margin_mm = cfg.phase_safe_margin_mm;
method.choose_better_pass = cfg.choose_better_pass;
end

function W = empty_window_result_local()
W = struct('blade_id', NaN, 'window_id', NaN, 'lap_range', [], ...
    'Bundle', [], 'CoreBundle', [], 'ExpandedBundle', [], ...
    'SeedTable', [], 'SelectedEO', [], 'Result', [], ...
    'CoreResult', [], 'ExpandedResult', [], 'ExpansionInfo', [], ...
    'PhaseSafeReferenceInfo', [], ...
    'SensorEtaInitMM', [], 'SensorEtaInitIQRMM', [], ...
    'ChosenPass', '', 'GapLibraryCompare', []);
end

function [bundleChosen, fitChosen, seedTableChosen, eoKeepChosen, expandInfo] = ...
        choose_better_identification_pass_local( ...
        bundleCore, fitCore, seedTableCore, eoKeepCore, ...
        bundleExpanded, fitExpanded, seedTableExpanded, eoKeepExpanded, expandInfo, cfg)
bundleChosen = bundleCore;
fitChosen = fitCore;
seedTableChosen = seedTableCore;
eoKeepChosen = eoKeepCore;
expandInfo.core_rmse = fitCore.weighted_voltage_rmse;
expandInfo.expanded_rmse = get_refined_rmse_local(fitExpanded);
expandInfo.chosen_pass = 'core';

if isempty(fitExpanded)
    expandInfo.status = 'identified_by_vp_topk_synchronous_template_core_only';
    return;
end

if cfg.choose_better_pass
    if fitExpanded.weighted_voltage_rmse < fitCore.weighted_voltage_rmse - cfg.better_pass_rmse_tolerance_v
        bundleChosen = bundleExpanded;
        fitChosen = fitExpanded;
        seedTableChosen = seedTableExpanded;
        eoKeepChosen = eoKeepExpanded;
        expandInfo.chosen_pass = 'expanded';
        expandInfo.status = 'identified_by_vp_topk_synchronous_template_phase_safe_expanded_kept';
    else
        expandInfo.status = 'identified_by_vp_topk_synchronous_template_core_kept_after_compare';
    end
else
    bundleChosen = bundleExpanded;
    fitChosen = fitExpanded;
    seedTableChosen = seedTableExpanded;
    eoKeepChosen = eoKeepExpanded;
    expandInfo.chosen_pass = 'expanded';
    expandInfo.status = 'identified_by_vp_topk_synchronous_template_phase_safe_expanded_forced';
end
end

function rmse = get_refined_rmse_local(fitRefined)
rmse = NaN;
if isempty(fitRefined)
    return;
end
if isfield(fitRefined, 'weighted_voltage_rmse')
    rmse = fitRefined.weighted_voltage_rmse;
end
end

function R = empty_gap_compare_result_local()
R = struct( ...
    'blade_id', NaN, ...
    'window_id', NaN, ...
    'lap_range', [], ...
    'direct_result', [], ...
    'gap_bundle', [], ...
    'seed_scan', table(), ...
    'eo_candidate_info', struct(), ...
    'gap_static_fit', [], ...
    'gap_model_fits', struct(), ...
    'gap_best_model_name', "", ...
    'gap_best_fit', [], ...
    'gap_eo_scan_table', table(), ...
    'summary', [], ...
    'trendRow', struct(), ...
    'cache_status', "");
end

function R = empty_gap_compare_window_result_local()
R = empty_gap_compare_result_local();
end

function [gapBranch, info] = prepare_gap_library_branch_local(P, cfg, waveformSourceInfo)
gapBranch = struct();
gapBranch.enabled = false;
gapBranch.model_mode = string(get_optional_field_local(P.step06, 'modelMode', 'direct_only'));
gapBranch.sensor_ids = get_optional_field_local(P.step06, 'gapLibrarySensors', [5 7]);
gapBranch.top_k_eo = max(1, round(get_optional_field_local(P.step06, 'gapLibraryTopKEO', 3)));
gapBranch.cfg = [];
gapBranch.corrected_lib_file = "";
gapBranch.CorrectedGapLibrary = [];
gapBranch.responseSurface = [];
gapBranch.cache_dir = "";

info = struct();
info.enabled = false;
info.modelMode = string(gapBranch.model_mode);
info.reason = 'disabled_by_mode';
info.correctedLibFile = "";
info.sensorIds = gapBranch.sensor_ids;

if strcmpi(gapBranch.model_mode, 'direct_only')
    return;
end

sensorIds = intersect(gapBranch.sensor_ids(:).', P.sensors.analysis(:).', 'stable');
if numel(sensorIds) < 2
    info.reason = 'gap_library_sensors_not_in_analysis_set';
    return;
end
gapBranch.sensor_ids = sensorIds;

cfgGap = struct();
cfgGap.targetBlade = NaN;
cfgGap.analysisSensors = sensorIds(:).';
cfgGap.eoCandidates = [];
cfgGap.amplitudeLimitMm = get_optional_field_local(P.step06, 'gapLibraryAmplitudeLimitMM', cfg.amplitude_limit_mm);
cfgGap.dxLimitMm = get_optional_field_local(P.step06, 'gapLibraryDxLimitMM', 0.60);
cfgGap.deltaGapLimitMm = get_optional_field_local(P.step06, 'gapLibraryDeltaGapLimitMM', 0.80);
cfgGap.deltaMuLimit = get_optional_field_local(P.step06, 'gapLibraryDeltaMuLimit', 0.030);
cfgGap.weightFloor = cfg.weight_floor;
cfgGap.overshootPenaltyWeight = max(cfg.overshoot_penalty_weight, 50);
cfgGap.staticRegWeight = get_optional_field_local(P.step06, 'gapLibraryStaticRegWeight', 0.50);
cfgGap.sensorEtaFallbackLimitMm = get_optional_field_local(P.step06, 'gapLibrarySensorEtaFallbackLimitMM', 0.06);
cfgGap.sensorEtaRegWeightMv = get_optional_field_local(P.step06, 'gapLibrarySensorEtaRegWeightMV', 0.35);
cfgGap.directVpTopK = gapBranch.top_k_eo;
cfgGap.eoCandidateMode = string(get_optional_field_local(P.step06, 'gapLibraryEOCandidateMode', 'direct_plus_vp'));
cfgGap.prevWindowCandidateMode = string(get_optional_field_local(P.step06, 'gapLibraryPrevWindowCandidateMode', 'soft'));
cfgGap.prevWindowRefreshEvery = max(0, round(get_optional_field_local(P.step06, 'gapLibraryPrevWindowRefreshEvery', 5)));
cfgGap.bundleSourceMode = string(get_optional_field_local(P.step06, 'gapLibraryBundleSourceMode', 'raw_gap_observation'));
cfgGap.modelNames = get_optional_field_local(P.step06, 'gapLibraryModelNames', {'fixed','gap_only','gap_tilt'});
cfgGap.phaseStarts = [0, pi/2, pi, 3*pi/2];
cfgGap.fminOptions = optimset('Display', 'off', 'MaxIter', 450, ...
    'MaxFunEvals', 1800, 'TolX', 1e-6, 'TolFun', 1e-6);

correctedLibOverride = strtrim(get_optional_field_local(P.step06, 'gapLibraryCorrectedLibFile', ''));
if ~isempty(correctedLibOverride)
    correctedLibFile = correctedLibOverride;
else
    correctedLibFile = find_gap_corrected_library_file_local(P, sensorIds, waveformSourceInfo);
end
if exist(correctedLibFile, 'file') ~= 2
    info.reason = 'corrected_gap_library_missing';
    info.correctedLibFile = string(correctedLibFile);
    return;
end

S = load(correctedLibFile, 'CorrectedGapLibrary');
if ~isfield(S, 'CorrectedGapLibrary')
    info.reason = 'corrected_gap_library_variable_missing';
    info.correctedLibFile = string(correctedLibFile);
    return;
end

CorrectedGapLibrary = S.CorrectedGapLibrary;
check_gap_library_sensor_coverage_local(CorrectedGapLibrary, sensorIds);

gapBranch.enabled = true;
gapBranch.corrected_lib_file = string(correctedLibFile);
gapBranch.CorrectedGapLibrary = CorrectedGapLibrary;
gapBranch.responseSurface = prepare_gap_response_surface_local(CorrectedGapLibrary.responseSurface);
gapBranch.cfg = cfgGap;
if isfield(P.files, 'gapCompareCacheDir')
    gapBranch.cache_dir = string(P.files.gapCompareCacheDir);
end

info.enabled = true;
info.reason = 'enabled';
info.correctedLibFile = string(correctedLibFile);
info.sensorIds = sensorIds;
end

function correctedLibFile = find_gap_corrected_library_file_local(P, sensorIds, waveformSourceInfo)
gapDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), '20241106_low_speed_gap_prior_decoupling', 'outputs');
bladeIds = P.identification.blades(:).';
sensorTag = ['S', sprintf('%d', sensorIds)];
correctedLibFile = '';
for i = 1:numel(bladeIds)
    candidate = fullfile(gapDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
        bladeIds(i), sensorTag));
    if exist(candidate, 'file') == 2
        correctedLibFile = candidate;
        return;
    end
end
if isfield(waveformSourceInfo, 'template_file')
    tplFile = char(waveformSourceInfo.template_file);
    tok = regexp(tplFile, '_B(\d+)_', 'tokens', 'once');
    if ~isempty(tok)
        bladeId = str2double(tok{1});
        candidate = fullfile(gapDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
            bladeId, sensorTag));
        if exist(candidate, 'file') == 2
            correctedLibFile = candidate;
            return;
        end
    end
end
correctedLibFile = fullfile(gapDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20241106_B%d_%s.mat', ...
    P.identification.blades(1), sensorTag));
end

function check_gap_library_sensor_coverage_local(CorrectedGapLibrary, sensorIds)
libIds = [CorrectedGapLibrary.sensor.sensorId];
for sid = sensorIds(:).'
    if ~ismember(sid, libIds)
        error('CorrectedGapLibrary does not contain gap-library sensor %d.', sid);
    end
end
end

function gapCompare = run_gap_library_compare_local(bundleChosen, bladeId, windowId, lapRange, fitChosen, eoKeepChosen, gapBranch, prevGapCompare, WaveformLibrary)
gapCompare = empty_gap_compare_result_local();
gapCompare.blade_id = bladeId;
gapCompare.window_id = windowId;
gapCompare.lap_range = lapRange;
gapCompare.direct_result = fitChosen;
directRmseMv = 1000 * fitChosen.weighted_voltage_rmse;

cacheMeta = build_gap_compare_cache_meta_local(bladeId, windowId, lapRange, fitChosen, eoKeepChosen, gapBranch);
[hasCache, cachedGapCompare] = try_load_gap_compare_cache_local(gapBranch, cacheMeta);
if hasCache
    gapCompare = cachedGapCompare;
    gapCompare.cache_status = "reused";
    return;
end

gapBundle = build_gap_bundle_for_compare_local(bundleChosen, bladeId, windowId, lapRange, gapBranch, WaveformLibrary);
eoAll = build_gap_eo_candidates_local(gapBundle, fitChosen, eoKeepChosen, gapBranch.cfg, prevGapCompare, windowId);
seedScan = solve_gap_direct_template_seed_scan_local(gapBundle, eoAll, gapBranch.cfg);
eoInfo = summarize_gap_eo_candidates_local(seedScan, fitChosen, eoKeepChosen, gapBranch.cfg, prevGapCompare, windowId);
directSelectedEO = eoInfo.selectedEO;
modelNames = string(gapBranch.cfg.modelNames);
modelFits = struct();
bestFit = [];
bestRmse = inf;
bestModelName = "";
staticFit = [];
eoScan = table();
for im = 1:numel(modelNames)
    modelName = char(modelNames(im));
    staticFitNow = fit_gap_static_only_local(gapBundle, gapBranch.cfg, modelName);
    modelRows = cell(numel(directSelectedEO), 1);
    modelBest = [];
    modelBestRmse = inf;
    for i = 1:numel(directSelectedEO)
        eo = directSelectedEO(i);
        fit = fit_gap_one_eo_local(gapBundle, eo, gapBranch.cfg, staticFitNow, modelName);
        modelRows{i} = table(eo, fit.freqHz, fit.amplitudeMm, fit.phaseRad, fit.dxMm, ...
            fit.weightedRmseMv, fit.unweightedRmseMv, ...
            'VariableNames', {'EO','frequencyHz','amplitudeMm','phaseRad','dxMm','weightedRmseMv','unweightedRmseMv'});
        if fit.weightedRmseMv < modelBestRmse
            modelBestRmse = fit.weightedRmseMv;
            modelBest = fit;
        end
    end
    if isempty(modelRows)
        eoScanNow = table();
    else
        eoScanNow = vertcat(modelRows{:});
    end
    modelFits.(modelName) = struct( ...
        'static_fit', staticFitNow, ...
        'best_fit', modelBest, ...
        'eo_scan_table', eoScanNow);
    if isempty(modelBest)
        modelBest = staticFitNow;
        modelBestRmse = staticFitNow.weightedRmseMv;
    end
    if modelBestRmse < bestRmse
        bestRmse = modelBestRmse;
        bestFit = modelBest;
        bestModelName = string(modelName);
        staticFit = staticFitNow;
        eoScan = eoScanNow;
    end
end
if isempty(bestFit)
    bestModelName = "fixed";
    staticFit = fit_gap_static_only_local(gapBundle, gapBranch.cfg, 'fixed');
    bestFit = staticFit;
end

gapCompare.gap_bundle = gapBundle;
gapCompare.seed_scan = seedScan;
gapCompare.eo_candidate_info = eoInfo;
gapCompare.gap_static_fit = staticFit;
gapCompare.gap_model_fits = modelFits;
gapCompare.gap_best_model_name = bestModelName;
gapCompare.gap_best_fit = bestFit;
gapCompare.gap_eo_scan_table = eoScan;
gapCompare.summary = struct( ...
    'blade_id', bladeId, ...
    'window_id', windowId, ...
    'direct_rmse_mv', directRmseMv, ...
    'gap_static_rmse_mv', staticFit.weightedRmseMv, ...
    'gap_best_rmse_mv', bestFit.weightedRmseMv, ...
    'gap_best_model', bestModelName, ...
    'gap_best_eo', bestFit.EO, ...
    'gap_best_freq_hz', bestFit.freqHz, ...
    'gap_best_amp_mm', bestFit.amplitudeMm, ...
    'gap_candidate_mode', string(eoInfo.mode));
gapCompare.trendRow = struct( ...
    'BladeID', bladeId, ...
    'WindowID', windowId, ...
    'LapStart', lapRange(1), ...
    'LapEnd', lapRange(end), ...
    'DirectEO', fitChosen.EO_id, ...
    'DirectFrequencyHz', fitChosen.fn_id, ...
    'DirectAmplitudeMM', fitChosen.A_id, ...
    'DirectDxMM', fitChosen.dx_c_id, ...
    'DirectWeightedVoltageRMSE', directRmseMv, ...
    'GapStaticWeightedVoltageRMSE', staticFit.weightedRmseMv, ...
    'GapBestModel', bestModelName, ...
    'GapEO', bestFit.EO, ...
    'GapFrequencyHz', bestFit.freqHz, ...
    'GapAmplitudeMM', bestFit.amplitudeMm, ...
    'GapDxMM', bestFit.dxMm, ...
    'GapWeightedVoltageRMSE', bestFit.weightedRmseMv, ...
    'DirectVsGapRmseDeltaMV', bestFit.weightedRmseMv - directRmseMv, ...
    'GapCandidateMode', string(eoInfo.mode), ...
    'GapVPTopK', mat2str(eoInfo.vpEO(:).'), ...
    'GapDirectEOSeed', eoInfo.directEO, ...
    'GapPrevWindowEOSeed', eoInfo.prevEO, ...
    'GapSelectedEO', mat2str(directSelectedEO(:).'), ...
    'GapBundleSourceMode', string(gapBranch.cfg.bundleSourceMode), ...
    'GapDeltaGapP5MM', local_gap_value_or_nan(bestFit, gapBranch.sensor_ids, 5), ...
    'GapDeltaGapP7MM', local_gap_value_or_nan(bestFit, gapBranch.sensor_ids, 7), ...
    'GapEtaP5MM', local_eta_value_or_nan(bestFit, gapBranch.sensor_ids, 5), ...
    'GapEtaP7MM', local_eta_value_or_nan(bestFit, gapBranch.sensor_ids, 7), ...
    'ChosenComparisonTag', string(local_choose_compare_tag(bestFit.weightedRmseMv, directRmseMv)));
gapCompare.cache_status = "built";
save_gap_compare_cache_local(gapBranch, cacheMeta, gapCompare);
end

function meta = build_gap_compare_cache_meta_local(bladeId, windowId, lapRange, fitChosen, eoKeepChosen, gapBranch)
cfg = gapBranch.cfg;
meta = struct();
meta.version = 1;
meta.blade_id = bladeId;
meta.window_id = windowId;
meta.lap_range = lapRange(:).';
meta.direct = [ ...
    get_optional_field_local(fitChosen, 'EO_id', NaN), ...
    get_optional_field_local(fitChosen, 'fn_id', NaN), ...
    get_optional_field_local(fitChosen, 'A_id', NaN), ...
    get_optional_field_local(fitChosen, 'phi_id_wrapped', NaN), ...
    get_optional_field_local(fitChosen, 'dx_c_id', NaN), ...
    get_optional_field_local(fitChosen, 'weighted_voltage_rmse', NaN), ...
    get_optional_field_local(fitChosen, 'point_count', NaN)];
meta.eo_keep = round(eoKeepChosen(:).');
meta.sensor_ids = gapBranch.sensor_ids(:).';
meta.corrected_lib_file = char(gapBranch.corrected_lib_file);
meta.corrected_lib_datenum = local_file_datenum_local(char(gapBranch.corrected_lib_file));
meta.cfg_numeric = [ ...
    cfg.amplitudeLimitMm, cfg.dxLimitMm, cfg.deltaGapLimitMm, cfg.deltaMuLimit, ...
    cfg.weightFloor, cfg.overshootPenaltyWeight, cfg.staticRegWeight, ...
    cfg.sensorEtaFallbackLimitMm, cfg.sensorEtaRegWeightMv, cfg.directVpTopK, ...
    cfg.prevWindowRefreshEvery];
meta.cfg_text = char(strjoin([ ...
    string(cfg.eoCandidateMode), string(cfg.prevWindowCandidateMode), ...
    string(cfg.bundleSourceMode), string(cfg.modelNames(:).')], '|'));
end

function [hasCache, gapCompare] = try_load_gap_compare_cache_local(gapBranch, meta)
hasCache = false;
gapCompare = [];
cacheFile = gap_compare_cache_file_local(gapBranch, meta);
if exist(cacheFile, 'file') ~= 2
    return;
end
loaded = load(cacheFile, 'GapCompareCache');
if ~isfield(loaded, 'GapCompareCache')
    return;
end
cache = loaded.GapCompareCache;
if ~isfield(cache, 'meta') || ~isfield(cache, 'gapCompare')
    return;
end
if ~isequaln(cache.meta, meta)
    return;
end
gapCompare = cache.gapCompare;
hasCache = true;
end

function save_gap_compare_cache_local(gapBranch, meta, gapCompare)
cacheFile = gap_compare_cache_file_local(gapBranch, meta);
cacheDir = fileparts(cacheFile);
if exist(cacheDir, 'dir') ~= 7
    mkdir(cacheDir);
end
GapCompareCache = struct('meta', meta, 'gapCompare', gapCompare);
save(cacheFile, 'GapCompareCache', '-v7.3');
end

function cacheFile = gap_compare_cache_file_local(gapBranch, meta)
cacheDir = gapBranch.cache_dir;
name = sprintf('GapCompare_B%d_W%03d_L%d_%d.mat', ...
    meta.blade_id, meta.window_id, meta.lap_range(1), meta.lap_range(end));
cacheFile = fullfile(cacheDir, name);
end

function gapBundle = build_gap_bundle_for_compare_local(bundleChosen, bladeId, windowId, lapRange, gapBranch, WaveformLibrary)
modeName = lower(strtrim(string(gapBranch.cfg.bundleSourceMode)));
if strcmpi(modeName, 'raw_gap_observation')
    try
        gapBundle = build_gap_bundle_from_raw_waveform_local(bladeId, windowId, lapRange, gapBranch, WaveformLibrary);
        return;
    catch ME
        warning('Gap raw observation bundle failed for B%d W%d. Falling back to direct bundle reuse. Reason: %s', ...
            bladeId, windowId, ME.message);
    end
end
gapBundle = build_gap_bundle_from_direct_bundle_local(bundleChosen, gapBranch);
end

function value = local_gap_value_or_nan(fit, sensorIds, sid)
value = NaN;
idx = find(sensorIds == sid, 1, 'first');
if isempty(idx) || ~isfield(fit, 'deltaGapMm') || numel(fit.deltaGapMm) < idx
    return;
end
value = fit.deltaGapMm(idx);
end

function value = local_eta_value_or_nan(fit, sensorIds, sid)
value = NaN;
idx = find(sensorIds == sid, 1, 'first');
if isempty(idx) || ~isfield(fit, 'sensorEtaMm') || numel(fit.sensorEtaMm) < idx
    return;
end
value = fit.sensorEtaMm(idx);
end

function tag = local_choose_compare_tag(gapRmse, directRmse)
if gapRmse < directRmse
    tag = 'gap_better';
elseif gapRmse > directRmse
    tag = 'direct_better';
else
    tag = 'tie';
end
end

function eoAll = build_gap_eo_candidates_local(bundle, fitChosen, eoKeepChosen, cfg, prevGapCompare, windowId)
eoAll = [];
if isfield(fitChosen, 'EO_id') && isfinite(fitChosen.EO_id)
    eoAll = round(fitChosen.EO_id);
end
eoAll = unique([eoAll(:).', round(eoKeepChosen(:).')], 'stable');
if isempty(eoAll)
    eoAll = round(build_eo_candidates_local(bundle.rot_freq_mean_hz, [0 1000], 0));
end
eoVpAll = round(build_eo_candidates_local(bundle.rot_freq_mean_hz, [0 1000], 0));
vpScan = solve_gap_direct_template_seed_scan_local(bundle, eoVpAll, cfg);
vpTop = vpScan.EO(1:min(cfg.directVpTopK, height(vpScan))).';
prevEO = get_prev_gap_window_eo_local(prevGapCompare, cfg, windowId);
mode = lower(strtrim(string(cfg.eoCandidateMode)));
switch mode
    case "vp"
        eoAll = vpTop;
    case "direct"
        eoAll = eoAll(:).';
    otherwise
        eoAll = unique([eoAll(:).', vpTop(:).'], 'stable');
end
if isfinite(prevEO)
    eoAll = unique([eoAll(:).', round(prevEO)], 'stable');
end
eoAll = eoAll(isfinite(eoAll) & eoAll > 0);
if isempty(eoAll)
    eoAll = vpTop;
end
if isempty(eoAll) && isfield(fitChosen, 'EO_id') && isfinite(fitChosen.EO_id)
    eoAll = round(fitChosen.EO_id);
end
end

function eo = get_prev_gap_window_eo_local(prevGapCompare, cfg, windowId)
eo = NaN;
if isempty(prevGapCompare) || ~isstruct(prevGapCompare)
    return;
end
if strcmpi(string(cfg.prevWindowCandidateMode), 'off')
    return;
end
if cfg.prevWindowRefreshEvery > 0 && isfinite(windowId) && mod(windowId - 1, cfg.prevWindowRefreshEvery) == 0
    return;
end
if ~isfield(prevGapCompare, 'gap_best_fit') || isempty(prevGapCompare.gap_best_fit)
    return;
end
fit = prevGapCompare.gap_best_fit;
if isfield(fit, 'EO') && isfinite(fit.EO) && fit.EO > 0
    eo = round(fit.EO);
end
end

function info = summarize_gap_eo_candidates_local(seedScan, fitChosen, eoKeepChosen, cfg, prevGapCompare, windowId)
info = struct();
info.mode = string(cfg.eoCandidateMode);
info.vpEO = [];
info.directEO = NaN;
info.prevEO = NaN;
info.prevReason = "not_used";
info.selectedEO = [];

if ~isempty(seedScan) && height(seedScan) > 0
    info.vpEO = seedScan.EO(1:min(cfg.directVpTopK, height(seedScan))).';
end
directEO = [];
if isfield(fitChosen, 'EO_id') && isfinite(fitChosen.EO_id)
    directEO = round(fitChosen.EO_id);
end
directEO = unique([directEO(:).', round(eoKeepChosen(:).')], 'stable');
if ~isempty(directEO)
    info.directEO = directEO(1);
end
prevEO = get_prev_gap_window_eo_local(prevGapCompare, cfg, windowId);
if isfinite(prevEO)
    info.prevEO = prevEO;
    info.prevReason = "prev_gap_window";
end

mode = lower(strtrim(string(cfg.eoCandidateMode)));
switch mode
    case "vp"
        selectedEO = info.vpEO;
    case "direct"
        selectedEO = directEO;
    otherwise
        selectedEO = unique([directEO(:).', info.vpEO(:).'], 'stable');
end
if isfinite(prevEO)
    selectedEO = unique([selectedEO(:).', prevEO], 'stable');
end
selectedEO = selectedEO(isfinite(selectedEO) & selectedEO > 0);
if isempty(selectedEO)
    selectedEO = info.vpEO;
end
if isempty(selectedEO) && ~isempty(directEO)
    selectedEO = directEO;
end
info.selectedEO = selectedEO(:).';
end

function bundle = build_gap_bundle_from_direct_bundle_local(bundleDirect, gapBranch)
sensorIds = gapBranch.sensor_ids(:).';
mask = ismember(bundleDirect.S, sensorIds) & isfinite(bundleDirect.X) & isfinite(bundleDirect.V) & ...
    isfinite(bundleDirect.T) & isfinite(bundleDirect.W) & isfinite(bundleDirect.Theta);

T = bundleDirect.T(mask);
X = bundleDirect.X(mask);
V = 1000 * bundleDirect.V(mask);
W = bundleDirect.W(mask);
S = bundleDirect.S(mask);
Theta = bundleDirect.Theta(mask);

[~, order] = sort(T);
T = T(order); X = X(order); V = V(order); W = W(order); S = S(order); Theta = Theta(order);

sensorIndex = NaN(size(S));
F0 = NaN(size(S));
Fx = NaN(size(S));
sensorInfo = repmat(struct('sensorId', NaN, 'corr', [], 'template', [], 'xDomain', [NaN, NaN]), 1, numel(sensorIds));
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    sensorMask = S == sid;
    sensorIndex(sensorMask) = is;
    sensorInfo(is).sensorId = sid;
    idxDirect = find(bundleDirect.sensor_ids == sid, 1, 'first');
    if isempty(idxDirect)
        error('Direct bundle does not contain gap-library sensor %d.', sid);
    end
    sensorInfo(is).template = bundleDirect.WaveformSensor(idxDirect);
    sensorInfo(is).corr = get_corrected_sensor_from_library_local(gapBranch.CorrectedGapLibrary, sid);
    sensorInfo(is).xDomain = bundleDirect.x_domain_by_sensor(idxDirect, :);
    [xCorr, vCorr] = build_gap_reference_curve_local(sensorInfo(is).corr);
    if numel(xCorr) >= 3
        dvCorr = gradient(vCorr, xCorr);
        F0(sensorMask) = interp1(xCorr, vCorr, X(sensorMask), 'pchip', NaN);
        Fx(sensorMask) = interp1(xCorr, dvCorr, X(sensorMask), 'pchip', NaN);
    end
end

bundle = struct();
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.X = X;
bundle.V = V;
bundle.W = W ./ max(W);
bundle.S = S;
bundle.sensor_index = sensorIndex;
bundle.sensor_ids = sensorIds;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.rot_freq_mean_hz = bundleDirect.rot_freq_mean_hz;
bundle.rot_rpm_mean = bundleDirect.rot_rpm_mean;
bundle.point_count = numel(T);
bundle.sensorInfo = sensorInfo;
bundle.responseSurface = gapBranch.responseSurface;
bundle.time_window = [min(T), max(T)];
bundle = finalize_gap_bundle_local(bundle, gapBranch.cfg);
end

function bundle = build_gap_bundle_from_raw_waveform_local(bladeId, windowId, lapRange, gapBranch, WaveformLibrary)
sensorIds = gapBranch.sensor_ids(:).';
X = [];
T = [];
V = [];
W = [];
S = [];
Theta = [];
F0 = [];
Fx = [];
sensorIndex = [];
sensorInfo = repmat(struct('sensorId', NaN, 'corr', [], 'template', [], 'xDomain', [NaN, NaN]), 1, numel(sensorIds));

for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    hit = find([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid, 1, 'first');
    if isempty(hit)
        error('WaveformLibrary B%d missing sensor %d for raw gap bundle.', bladeId, sid);
    end
    Sdata = WaveformLibrary.Blade(bladeId).Sensor(hit);
    Tpl = template_from_sensor_library_local(Sdata);
    corr = get_corrected_sensor_from_library_local(gapBranch.CorrectedGapLibrary, sid);
    [tRaw, xRaw, vRaw, thetaRaw] = concatenate_raw_laps_local(Sdata.Lap(lapRange));
    phaseRefDummy = struct('reference_source', 'template_inverse_apparent_query');
    [tSel, xSel, vSel, thetaSel, ~, ~, wSel] = ...
        select_phase_safe_sensor_points_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, phaseRefDummy, is, gap_branch_cfg_as_step06_local(gapBranch.cfg), sid);
    if numel(tSel) < 8
        continue;
    end
    [f0Sel, ~] = eval_gap_corrected_template_local(Tpl, gapBranch.responseSurface, corr, corr.g0Mm, xSel);
    fxSel = eval_gap_corrected_template_derivative_local(Tpl, gapBranch.responseSurface, corr, corr.g0Mm, xSel, 1e-3);
    keep = isfinite(tSel) & isfinite(xSel) & isfinite(vSel) & isfinite(thetaSel) & ...
        isfinite(wSel) & isfinite(f0Sel) & isfinite(fxSel);
    if nnz(keep) < 8
        continue;
    end
    X = [X; xSel(keep)]; %#ok<AGROW>
    T = [T; tSel(keep)]; %#ok<AGROW>
    V = [V; 1000 * vSel(keep)]; %#ok<AGROW>
    W = [W; wSel(keep)]; %#ok<AGROW>
    S = [S; repmat(sid, nnz(keep), 1)]; %#ok<AGROW>
    Theta = [Theta; thetaSel(keep)]; %#ok<AGROW>
    F0 = [F0; f0Sel(keep)]; %#ok<AGROW>
    Fx = [Fx; fxSel(keep)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, nnz(keep), 1)]; %#ok<AGROW>
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).corr = corr;
    sensorInfo(is).template = Sdata;
    sensorInfo(is).xDomain = Tpl.x_domain;
end

if isempty(T)
    error('Raw gap observation bundle is empty for B%d W%d.', bladeId, windowId);
end

[T, order] = sort(T);
X = X(order);
V = V(order);
W = W(order);
S = S(order);
Theta = Theta(order);
F0 = F0(order);
Fx = Fx(order);
sensorIndex = sensorIndex(order);

bundle = struct();
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.X = X;
bundle.V = V;
bundle.W = W ./ max(W);
bundle.S = S;
bundle.sensor_index = sensorIndex;
bundle.sensor_ids = sensorIds;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.rot_freq_mean_hz = estimate_gap_bundle_rotfreq_local(WaveformLibrary, bladeId, lapRange, sensorIds);
bundle.rot_rpm_mean = 60 * bundle.rot_freq_mean_hz;
bundle.point_count = numel(T);
bundle.sensorInfo = sensorInfo;
bundle.responseSurface = gapBranch.responseSurface;
bundle.time_window = [min(T), max(T)];
bundle = finalize_gap_bundle_local(bundle, gapBranch.cfg);
end

function responseSurface = prepare_gap_response_surface_local(responseSurface)
responseSurface.xGrid = responseSurface.xGrid(:);
responseSurface.xMin = min(responseSurface.xGrid);
responseSurface.xMax = max(responseSurface.xGrid);
responseSurface.gMin = min(responseSurface.gTrainMm(:));
responseSurface.gMax = max(responseSurface.gTrainMm(:));
dx = diff(responseSurface.xGrid);
responseSurface.hasUniformXGrid = ~isempty(dx) && max(abs(dx - dx(1))) <= max(1e-12, 1e-10 * abs(dx(1)));
if responseSurface.hasUniformXGrid
    responseSurface.xStep = dx(1);
    responseSurface.invXStep = 1 / responseSurface.xStep;
    responseSurface.coeffBase = responseSurface.coeff(1:end-1, :);
    responseSurface.coeffSlope = diff(responseSurface.coeff, 1, 1) .* responseSurface.invXStep;
else
    responseSurface.xStep = NaN;
    responseSurface.invXStep = NaN;
    responseSurface.coeffBase = [];
    responseSurface.coeffSlope = [];
end
responseSurface.interpSignature = [ ...
    numel(responseSurface.xGrid), ...
    responseSurface.xGrid(1), ...
    responseSurface.xGrid(end), ...
    responseSurface.g0Mm, ...
    sum(responseSurface.coeff(:, 1), 'omitnan'), ...
    sum(responseSurface.coeff(:, 2), 'omitnan'), ...
    sum(responseSurface.coeff(:, 3), 'omitnan')];
end

function bundle = finalize_gap_bundle_local(bundle, cfg)
nSensor = numel(bundle.sensor_ids);
bundle.valid_base = isfinite(bundle.V) & isfinite(bundle.W) & isfinite(bundle.X) & isfinite(bundle.Theta);
bundle.weight_eval = max(bundle.W, cfg.weightFloor);
bundle.weight_eval_sum = sum(bundle.weight_eval(bundle.valid_base));
bundle.sensor_point_idx = cell(1, nSensor);
bundle.tau_by_point = zeros(size(bundle.X));
bundle.x_scale_by_point = ones(size(bundle.X));
bundle.mu_gap_by_point = zeros(size(bundle.X));
bundle.g0_by_point = zeros(size(bundle.X));
bundle.voltage_gain_by_point = ones(size(bundle.X));
bundle.voltage_offset_by_point = zeros(size(bundle.X));
for is = 1:nSensor
    idx = find(bundle.sensor_index == is);
    bundle.sensor_point_idx{is} = idx;
    corr = bundle.sensorInfo(is).corr;
    bundle.tau_by_point(idx) = corr.tauMm;
    bundle.x_scale_by_point(idx) = corr.xScale;
    bundle.mu_gap_by_point(idx) = corr.muGapPerXMm;
    bundle.g0_by_point(idx) = corr.g0Mm;
    bundle.voltage_gain_by_point(idx) = corr.voltageGain;
    bundle.voltage_offset_by_point(idx) = corr.voltageOffsetMv;
end
bundle.x_minus_tau_by_point = bundle.X - bundle.tau_by_point;
bundle.sensor_eta_limit_by_sensor = collect_gap_sensor_eta_limit_local(bundle, cfg);
end

function cfgStep06 = gap_branch_cfg_as_step06_local(cfgGap)
cfgStep06 = struct();
cfgStep06.weight_floor = cfgGap.weightFloor;
cfgStep06.domain_selection_mode = 'hard';
cfgStep06.domain_soft_margin_mm = 0.10;
cfgStep06.dynamic_effective_mode = 'gradient';
cfgStep06.dynamic_template_gradient_min_ratio = 0.08;
cfgStep06.dynamic_time_gradient_min_ratio = 0.15;
cfgStep06.dynamic_peak_quantile = 85;
cfgStep06.default_sensor_threshold = 0.5;
cfgStep06.query_guard_mode = 'adaptive';
cfgStep06.query_guard_mm = cfgGap.amplitudeLimitMm + cfgGap.dxLimitMm + 0.05;
cfgStep06.query_guard_quantile = 95;
cfgStep06.query_guard_safety_mm = 0.05;
cfgStep06.query_guard_min_mm = 0.12;
cfgStep06.query_guard_max_mm = cfgGap.amplitudeLimitMm + cfgGap.dxLimitMm + 0.05;
cfgStep06.sensor_query_guard_scale_table = [5 0.80; 7 0.75];
cfgStep06.sensor_domain_expand_mm_table = [];
cfgStep06.phase_safe_min_point_count_per_sensor = 6;
cfgStep06.phase_safe_margin_mm = 0.02;
cfgStep06.amplitude_limit_mm = cfgGap.amplitudeLimitMm;
cfgStep06.dx_c_limit_mm = cfgGap.dxLimitMm;
cfgStep06.pulse_selection_mode = 'single';
cfgStep06.domain_margin_mm = 0.02;
end

function rotFreqHz = estimate_gap_bundle_rotfreq_local(WaveformLibrary, bladeId, lapRange, sensorIds)
periods = [];
for sid = sensorIds(:).'
    hit = find([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid, 1, 'first');
    if isempty(hit)
        continue;
    end
    Sdata = WaveformLibrary.Blade(bladeId).Sensor(hit);
    if numel(Sdata.Lap) < max(lapRange)
        continue;
    end
    tPeak = arrayfun(@(x) x.t_peak, Sdata.Lap(lapRange));
    dt = diff(tPeak(:));
    periods = [periods; dt(isfinite(dt) & dt > eps)]; %#ok<AGROW>
end
if isempty(periods)
    error('Could not estimate raw gap bundle rotation speed for B%d laps %s.', bladeId, mat2str(lapRange));
end
rotFreqHz = 1 / median(periods, 'omitnan');
end

function [model, overshoot] = eval_gap_corrected_template_local(Tpl, responseSurface, corr, g, xOpr)
xLowRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xLowRaw, xLo), xHi);
overshootLow = max(xLo - xLowRaw, 0) + max(xLowRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), 1000 * Tpl.v_grid(:), xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_gap_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, overBase] = eval_gap_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function dMdx = eval_gap_corrected_template_derivative_local(Tpl, responseSurface, corr, g, xOpr, h)
[mp, ~] = eval_gap_corrected_template_local(Tpl, responseSurface, corr, g, xOpr + h);
[mm, ~] = eval_gap_corrected_template_local(Tpl, responseSurface, corr, g, xOpr - h);
dMdx = (mp - mm) ./ (2 * h);
end

function [model, overshoot] = eval_gap_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr(:) - tau);
gEffRaw = g0 + mu .* (xOpr(:) - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
xLib = min(max(xLibRaw, xMin), xMax);
gEff = min(max(gEffRaw, gMin), gMax);
overshoot = hypot(max(xMin - xLibRaw, 0) + max(xLibRaw - xMax, 0), ...
    max(gMin - gEffRaw, 0) + max(gEffRaw - gMax, 0));
model = eval_gap_response_surface_local(responseSurface, gEff, xLib);
end

function corr = get_corrected_sensor_from_library_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('CorrectedGapLibrary is missing sensor %d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function [xRef, vRef] = build_gap_reference_curve_local(corr)
xRef = [];
vRef = [];
if isfield(corr, 'x') && isfield(corr, 'vFitMv')
    xRef = corr.x(:);
    vRef = corr.vFitMv(:);
elseif isfield(corr, 'x') && isfield(corr, 'vLowMv')
    xRef = corr.x(:);
    vRef = corr.vLowMv(:);
end
valid = isfinite(xRef) & isfinite(vRef);
xRef = xRef(valid);
vRef = vRef(valid);
if isempty(xRef)
    return;
end
[xRef, ia] = unique(xRef, 'stable');
vRef = vRef(ia);
end

function T = solve_gap_direct_template_seed_scan_local(bundle, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
rows = cell(numel(eoCandidates), 1);
useLinearized = isfield(bundle, 'F0') && isfield(bundle, 'Fx') && ...
    any(isfinite(bundle.F0)) && any(isfinite(bundle.Fx));
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    if useLinearized
        valid = isfinite(bundle.F0(:)) & isfinite(bundle.Fx(:)) & isfinite(bundle.V(:));
        s1 = sin(eo .* bundle.Theta(valid));
        c1 = cos(eo .* bundle.Theta(valid));
        basis = [-bundle.Fx(valid), -bundle.Fx(valid) .* s1, -bundle.Fx(valid) .* c1];
        y = bundle.V(valid) - bundle.F0(valid);
        wSqrt = sqrt(max(bundle.W(valid), cfg.weightFloor));
        coeff = (basis .* wSqrt) \ (y .* wSqrt);
        dx = coeff(1);
        amp = hypot(coeff(2), coeff(3));
        phi = atan2(coeff(3), coeff(2));
        residual = y - basis * coeff;
        weightedRmse = sqrt(sum(max(bundle.W(valid), cfg.weightFloor) .* residual.^2) / max(sum(max(bundle.W(valid), cfg.weightFloor)), eps));
    else
        s1 = sin(eo .* bundle.Theta(:));
        c1 = cos(eo .* bundle.Theta(:));
        basis = [s1, c1, ones(size(s1))];
        y = bundle.V(:);
        wSqrt = sqrt(max(bundle.W(:), cfg.weightFloor));
        coeff = (basis .* wSqrt) \ (y .* wSqrt);
        dx = 0;
        amp = hypot(coeff(1), coeff(2));
        phi = atan2(coeff(2), coeff(1));
        residual = y - basis * coeff;
        weightedRmse = sqrt(sum(max(bundle.W(:), cfg.weightFloor) .* residual.^2) / max(sum(max(bundle.W(:), cfg.weightFloor)), eps));
    end
    detail = struct();
    detail.EO = eo;
    detail.frequencyHz = eo * bundle.rot_freq_mean_hz;
    detail.amplitudeMm = min(abs(amp), cfg.amplitudeLimitMm);
    detail.phaseRad = phi;
    detail.dxMm = dx;
    detail.weightedRmseMv = weightedRmse;
    detail.unweightedRmseMv = sqrt(mean(residual.^2, 'omitnan'));
    rows{i} = struct2table(detail);
end
T = vertcat(rows{:});
T = sortrows(T, {'weightedRmseMv','EO'}, {'ascend','ascend'});
end

function fit = fit_gap_static_only_local(bundle, cfg, modeName)
if nargin < 3 || isempty(modeName)
    modeName = 'fixed';
end
nSensor = numel(bundle.sensor_ids);
theta0 = zeros(1, 3 + nSensor + gap_mode_extra_count_local(modeName, nSensor));
fun = @(thetaRaw) gap_static_objective_local(thetaRaw, bundle, cfg, modeName);
theta = fminsearch(fun, theta0, cfg.fminOptions);
theta = bound_gap_static_params_local(theta, cfg, nSensor, bundle, modeName);
detail = evaluate_gap_model_local(gap_static_to_dynamic_local(theta, modeName), 0, bundle, cfg, false, modeName);
fit = detail;
fit.EO = 0;
fit.freqHz = 0;
fit.amplitudeMm = 0;
fit.phaseRad = 0;
fit.mode = string(modeName);
end

function score = gap_static_objective_local(thetaRaw, bundle, cfg, modeName)
nSensor = numel(bundle.sensor_ids);
theta = bound_gap_static_params_local(thetaRaw, cfg, nSensor, bundle, modeName);
detail = evaluate_gap_model_local(gap_static_to_dynamic_local(theta, modeName), 0, bundle, cfg, false, modeName);
score = detail.weightedRmseMv + cfg.staticRegWeight * norm(detail.deltaGapMm) + detail.etaPenaltyMv + detail.tiltPenaltyMv;
end

function theta = bound_gap_static_params_local(theta, cfg, nSensor, bundle, modeName)
theta = theta(:).';
need = 3 + nSensor + gap_mode_extra_count_local(modeName, nSensor);
if numel(theta) < need
    theta(end+1:need) = 0;
end
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
theta(3) = min(max(theta(3), -cfg.dxLimitMm), cfg.dxLimitMm);
idx = 4:(3+nSensor);
etaLimit = collect_gap_sensor_eta_limit_local(bundle, cfg);
theta(idx) = min(max(theta(idx), -etaLimit), etaLimit);
if any(strcmpi(modeName, {'gap_only','gap_tilt'}))
    idxGap = 4+nSensor:(3+2*nSensor);
    theta(idxGap) = min(max(theta(idxGap), -cfg.deltaGapLimitMm), cfg.deltaGapLimitMm);
end
if strcmpi(modeName, 'gap_tilt')
    idxMu = 4+2*nSensor:(3+3*nSensor);
    theta(idxMu) = min(max(theta(idxMu), -cfg.deltaMuLimit), cfg.deltaMuLimit);
end
theta(4:3+nSensor) = min(max(theta(4:3+nSensor), -etaLimit), etaLimit);
end

function dynamicTheta = gap_static_to_dynamic_local(staticTheta, modeName)
if nargin < 2
    modeName = 'fixed';
end
modeName = char(modeName); %#ok<NASGU>
dynamicTheta = staticTheta(:).';
dynamicTheta = [0, 0, dynamicTheta];
end

function fit = fit_gap_one_eo_local(bundle, eo, cfg, staticFit, modeName)
if nargin < 5 || isempty(modeName)
    modeName = 'fixed';
end
nSensor = numel(bundle.sensor_ids);
eoContext = build_gap_eo_context_local(bundle, eo);
bestScore = inf;
bestDetail = [];
ampStarts = [0.005, 0.02, 0.05];
for A0 = ampStarts
    for phi0 = cfg.phaseStarts
        theta0 = [A0, phi0, staticFit.dxMm, staticFit.sensorEtaMm(:).'];
        if any(strcmpi(modeName, {'gap_only','gap_tilt'}))
            theta0 = [theta0, staticFit.deltaGapMm(:).'];
        end
        if strcmpi(modeName, 'gap_tilt')
            theta0 = [theta0, staticFit.deltaMuGapPerXMm(:).'];
        end
        fun = @(thetaRaw) gap_dynamic_objective_local(thetaRaw, eo, bundle, cfg, modeName, eoContext);
        theta = fminsearch(fun, theta0, cfg.fminOptions);
        theta = bound_gap_dynamic_params_local(theta, cfg, nSensor, bundle, modeName);
        detail = evaluate_gap_model_local(theta, eo, bundle, cfg, true, modeName, eoContext);
        if detail.weightedRmseMv < bestScore
            bestScore = detail.weightedRmseMv;
            bestDetail = detail;
        end
    end
end
fit = bestDetail;
fit.mode = string(modeName);
end

function eoContext = build_gap_eo_context_local(bundle, eo)
eoContext = struct();
eoContext.eo = eo;
eoContext.freqHz = eo * bundle.rot_freq_mean_hz;
eoContext.sinTheta = sin(eo .* bundle.Theta(:));
eoContext.cosTheta = cos(eo .* bundle.Theta(:));
end

function score = gap_dynamic_objective_local(thetaRaw, eo, bundle, cfg, modeName, eoContext)
nSensor = numel(bundle.sensor_ids);
theta = bound_gap_dynamic_params_local(thetaRaw, cfg, nSensor, bundle, modeName);
if nargin < 6
    eoContext = [];
end
detail = evaluate_gap_model_local(theta, eo, bundle, cfg, true, modeName, eoContext);
score = detail.weightedRmseMv + cfg.staticRegWeight * norm(detail.deltaGapMm) + detail.etaPenaltyMv + ...
    detail.tiltPenaltyMv + cfg.overshootPenaltyWeight * detail.overshootRmseMm;
end

function theta = bound_gap_dynamic_params_local(theta, cfg, nSensor, bundle, modeName)
theta = theta(:).';
need = 3 + nSensor + gap_mode_extra_count_local(modeName, nSensor);
if numel(theta) < need
    theta(end+1:need) = 0;
end
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
theta(3) = min(max(theta(3), -cfg.dxLimitMm), cfg.dxLimitMm);
etaLimit = collect_gap_sensor_eta_limit_local(bundle, cfg);
etaIdx = 4:3+nSensor;
theta(etaIdx) = min(max(theta(etaIdx), -etaLimit), etaLimit);
if any(strcmpi(modeName, {'gap_only','gap_tilt'}))
    idxGap = 4+nSensor:(3+2*nSensor);
    theta(idxGap) = min(max(theta(idxGap), -cfg.deltaGapLimitMm), cfg.deltaGapLimitMm);
end
if strcmpi(modeName, 'gap_tilt')
    idxMu = 4+2*nSensor:(3+3*nSensor);
    theta(idxMu) = min(max(theta(idxMu), -cfg.deltaMuLimit), cfg.deltaMuLimit);
end
end

function detail = evaluate_gap_model_local(theta, eo, bundle, cfg, includeDynamic, modeName, eoContext)
if nargin < 7
    eoContext = [];
end
if nargin < 6 || isempty(modeName)
    modeName = 'fixed';
end
nSensor = numel(bundle.sensor_ids);
theta = bound_gap_dynamic_params_local(theta, cfg, nSensor, bundle, modeName);
A = theta(1);
phi = theta(2);
dx0 = theta(3);
eta = theta(4:3+nSensor).';
dg = zeros(nSensor, 1);
dmu = zeros(nSensor, 1);
offset = 3 + nSensor;
if any(strcmpi(modeName, {'gap_only','gap_tilt'}))
    dg = theta(offset + (1:nSensor)).';
    offset = offset + nSensor;
end
if strcmpi(modeName, 'gap_tilt')
    dmu = theta(offset + (1:nSensor)).';
end
if includeDynamic && eo > 0
    if ~isempty(eoContext) && isfield(eoContext, 'eo') && eoContext.eo == eo
        freqHz = eoContext.freqHz;
        u = A .* (eoContext.sinTheta .* cos(phi) + eoContext.cosTheta .* sin(phi));
    else
        freqHz = eo * bundle.rot_freq_mean_hz;
        u = A .* sin(eo .* bundle.Theta(:) + phi);
    end
else
    freqHz = 0;
    u = zeros(size(bundle.T));
end

etaByPoint = eta(bundle.sensor_index);
dgByPoint = dg(bundle.sensor_index);
xLocal = bundle.x_minus_tau_by_point - dx0 - u - etaByPoint;
xLib = bundle.x_scale_by_point .* xLocal;
gEff = bundle.g0_by_point + dgByPoint + bundle.mu_gap_by_point .* xLocal;
[vPred, overshoot] = eval_gap_pretransformed_library_local( ...
    bundle.responseSurface, bundle.voltage_gain_by_point, bundle.voltage_offset_by_point, xLib, gEff);

valid = bundle.valid_base & isfinite(vPred);
residual = bundle.V - vPred;
w = bundle.weight_eval(valid);
weightedRmse = sqrt(sum(w .* residual(valid).^2) / max(bundle.weight_eval_sum, eps));
unweightedRmse = sqrt(mean(residual(valid).^2, 'omitnan'));
etaPenaltyMv = compute_gap_eta_penalty_local(eta, bundle, cfg);

detail = struct();
detail.EO = eo;
detail.freqHz = freqHz;
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx0;
detail.deltaGapMm = dg(:);
detail.deltaMuGapPerXMm = dmu(:);
detail.sensorEtaMm = eta(:);
detail.uMm = u;
detail.uStdMm = std(u(valid), 'omitnan');
detail.uPeakToPeakMm = max(u(valid), [], 'omitnan') - min(u(valid), [], 'omitnan');
detail.uMinMm = min(u(valid), [], 'omitnan');
detail.uMaxMm = max(u(valid), [], 'omitnan');
detail.amplitudeAtLimit = abs(A - cfg.amplitudeLimitMm) < 1e-4;
detail.Vpred = vPred;
detail.weightedRmseMv = weightedRmse;
detail.unweightedRmseMv = unweightedRmse;
detail.overshootRmseMm = sqrt(mean(overshoot(valid).^2, 'omitnan'));
detail.validPointCount = nnz(valid);
detail.etaPenaltyMv = etaPenaltyMv;
detail.tiltPenaltyMv = cfg.staticRegWeight * sqrt(mean((dmu ./ max(cfg.deltaMuLimit, eps)).^2, 'omitnan'));
end

function n = gap_mode_extra_count_local(modeName, nSensor)
n = double(any(strcmpi(modeName, {'gap_only','gap_tilt'}))) * nSensor + double(strcmpi(modeName, 'gap_tilt')) * nSensor;
end

function [model, overshoot] = eval_gap_corrected_library_local(responseSurface, corr, dg, xOpr)
xLocal = xOpr - corr.tauMm;
xLib = corr.xScale .* xLocal;
gEff = corr.g0Mm + dg + corr.muGapPerXMm .* xLocal;
[model, overshoot] = eval_gap_pretransformed_library_local( ...
    responseSurface, corr.voltageGain, corr.voltageOffsetMv, xLib, gEff);
end

function [model, overshoot] = eval_gap_pretransformed_library_local(responseSurface, voltageGain, voltageOffsetMv, xLib, gEff)
xMin = responseSurface.xMin;
xMax = responseSurface.xMax;
gMin = responseSurface.gMin;
gMax = responseSurface.gMax;
overX = max(0, xMin - xLib) + max(0, xLib - xMax);
overG = max(0, gMin - gEff) + max(0, gEff - gMax);
overshoot = hypot(overX, overG);
xEval = min(max(xLib, xMin), xMax);
gEval = min(max(gEff, gMin), gMax);
raw = eval_gap_response_surface_local(responseSurface, gEval, xEval);
model = voltageGain .* raw + voltageOffsetMv;
end

function F = eval_gap_response_surface_local(responseSurface, g, xq)
persistent cache
signature = responseSurface.interpSignature;
if isempty(cache) || ~isfield(cache, 'signature') || ~isequal(cache.signature, signature)
    xGrid = responseSurface.xGrid(:);
    cache = struct();
    cache.signature = signature;
    cache.B0 = griddedInterpolant(xGrid, responseSurface.coeff(:, 1), 'linear', 'none');
    cache.B1 = griddedInterpolant(xGrid, responseSurface.coeff(:, 2), 'linear', 'none');
    cache.B2 = griddedInterpolant(xGrid, responseSurface.coeff(:, 3), 'linear', 'none');
end
B0 = cache.B0(xq);
B1 = cache.B1(xq);
B2 = cache.B2(xq);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function etaLimit = collect_gap_sensor_eta_limit_local(bundle, cfg)
if isfield(bundle, 'sensor_eta_limit_by_sensor') && ...
        numel(bundle.sensor_eta_limit_by_sensor) == numel(bundle.sensor_ids)
    etaLimit = bundle.sensor_eta_limit_by_sensor(:).';
    return;
end
nSensor = numel(bundle.sensor_ids);
etaLimit = cfg.sensorEtaFallbackLimitMm * ones(1, nSensor);
for is = 1:nSensor
    corr = bundle.sensorInfo(is).corr;
    if isfield(corr, 'etaLimitMm') && isfinite(corr.etaLimitMm) && corr.etaLimitMm > 0
        etaLimit(is) = corr.etaLimitMm;
    end
end
etaLimit = max(etaLimit, 1e-3);
end

function penalty = compute_gap_eta_penalty_local(eta, bundle, cfg)
eta = eta(:);
etaLimit = collect_gap_sensor_eta_limit_local(bundle, cfg).';
penalty = cfg.sensorEtaRegWeightMv * sqrt(mean((eta ./ max(etaLimit, eps)).^2, 'omitnan'));
end

function bundle = attach_template_and_speed_local(bundle, WaveformLibrary, bladeId, lapRange, cfg)
bundle.interp_v = cell(numel(bundle.sensor_ids), 1);
bundle.x_domain_by_sensor = nan(numel(bundle.sensor_ids), 2);
for is = 1:numel(bundle.sensor_ids)
    sid = bundle.sensor_ids(is);
    Sdata = WaveformLibrary.Blade(bladeId).Sensor([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid);
    if isempty(Sdata)
        error('WaveformLibrary B%d is missing CH%d.', bladeId, sid);
    end
    x = Sdata.template_x(:);
    v = Sdata.template_v(:);
    [x, ia] = unique(x, 'stable');
    v = v(ia);
    bundle.interp_v{is} = griddedInterpolant(x, v, 'pchip', 'none');
    if isfield(Sdata, 'template_domain') && numel(Sdata.template_domain) >= 2
        bundle.x_domain_by_sensor(is, :) = Sdata.template_domain(:).';
    else
        bundle.x_domain_by_sensor(is, :) = [min(x), max(x)];
    end
end
[rotFreqHz, rotRpm] = estimate_window_speed_local(WaveformLibrary, bladeId, lapRange, bundle.sensor_ids);
bundle.rot_freq_mean_hz = rotFreqHz;
bundle.rot_rpm_mean = rotRpm;
bundle.overshoot_penalty_weight = cfg.overshoot_penalty_weight;
bundle.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
bundle.sensor_eta_limit_by_sensor = build_sensor_eta_limit_by_sensor_local(bundle, cfg);
bundle.sensor_eta_reg_weight_v_per_mm = cfg.sensor_eta_reg_weight_v_per_mm;
bundle.query_guard_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm + 0.05;
end

function [phaseRef, info] = build_phase_safe_reference_local(prevPhaseResult, fitCore, seedTableCore, windowId, cfg)
if ~cfg.phase_safe_expansion_enable
    phaseRef = [];
    info = make_phase_safe_reference_info_local('', 'phase_safe_expansion_disabled', NaN, false);
    return;
end

if strcmpi(cfg.phase_safe_reference_mode, 'prev_window')
    [phaseRef, info] = build_prev_window_phase_safe_reference_local(prevPhaseResult, windowId, cfg);
    if isempty(phaseRef)
        [phaseRef, info] = build_current_window_phase_safe_reference_local( ...
            fitCore, seedTableCore, windowId, cfg, info.reason);
    end
else
    [phaseRef, info] = build_current_window_phase_safe_reference_local( ...
        fitCore, seedTableCore, windowId, cfg, 'configured_current_window_reference');
end
end

function info = make_phase_safe_reference_info_local(source, reason, sourceWindowId, usedPreviousWindow)
info = struct();
info.source = source;
info.reason = reason;
info.source_window_id = sourceWindowId;
info.used_previous_window = usedPreviousWindow;
end

function [phaseRef, info] = build_current_window_phase_safe_reference_local( ...
        fitCore, seedTableCore, windowId, cfg, fallbackReason)
mode = cfg.phase_safe_reference_mode;
if strcmpi(mode, 'prev_window') || strcmpi(mode, 'current_window')
    mode = cfg.phase_safe_fallback_mode;
end

if strcmpi(mode, 'linear_vp')
    phaseRef = build_phase_safe_reference_from_seed_local(seedTableCore, cfg);
elseif strcmpi(mode, 'template_inverse')
    phaseRef = build_template_inverse_phase_safe_reference_local();
elseif strcmpi(mode, 'final_core')
    phaseRef = fitCore;
    phaseRef.reference_source = 'current_core_final_waveform';
else
    error('Unsupported Step06 phase-safe fallback mode: %s.', mode);
end

if isempty(fallbackReason)
    reason = mode;
else
    reason = sprintf('%s_then_%s', fallbackReason, mode);
end
info = make_phase_safe_reference_info_local(phaseRef.reference_source, reason, windowId, false);
end

function phaseRef = build_phase_safe_reference_from_seed_local(seedTable, cfg)
seed = seedTable(1);
phaseRef = struct();
phaseRef.EO_id = seed.EO;
phaseRef.A_id = min(abs(seed.A), abs(cfg.amplitude_limit_mm));
phaseRef.phi_id_wrapped = wrap_to_pi_local(seed.phi);
phaseRef.dx_c_id = max(min(seed.dx_c, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phaseRef.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
phaseRef.weighted_voltage_rmse = seed.weighted_voltage_rmse;
phaseRef.plain_voltage_rmse = seed.plain_voltage_rmse;
phaseRef.reference_source = 'linear_vp_seed';
phaseRef.CandidateTable = struct2table(rmfield(seedTable, {'sensor_eta'}));
end

function phaseRef = build_template_inverse_phase_safe_reference_local()
phaseRef = struct();
phaseRef.EO_id = NaN;
phaseRef.A_id = NaN;
phaseRef.phi_id_wrapped = NaN;
phaseRef.dx_c_id = NaN;
phaseRef.sensor_eta_id = [];
phaseRef.weighted_voltage_rmse = NaN;
phaseRef.plain_voltage_rmse = NaN;
phaseRef.reference_source = 'template_inverse_apparent_query';
phaseRef.CandidateTable = table();
phaseRef.Coverage = struct('clamp_fraction', NaN);
end

function [phaseRef, info] = build_prev_window_phase_safe_reference_local(prevResult, windowId, cfg)
phaseRef = [];
[ok, reason] = is_prev_window_phase_reference_usable_local(prevResult, windowId, cfg);
if ~ok
    info = make_phase_safe_reference_info_local('', reason, NaN, false);
    return;
end

phaseRef = struct();
phaseRef.EO_id = prevResult.EO_id;
phaseRef.A_id = min(abs(prevResult.A_id), abs(cfg.amplitude_limit_mm));
phaseRef.phi_id_wrapped = wrap_to_pi_local(prevResult.phi_id_wrapped);
phaseRef.dx_c_id = max(min(prevResult.dx_c_id, abs(cfg.dx_c_limit_mm)), -abs(cfg.dx_c_limit_mm));
phaseRef.sensor_eta_id = zeros(1, numel(cfg.analysis_sensors));
if isfield(prevResult, 'sensor_eta_id') && ~isempty(prevResult.sensor_eta_id)
    n = min(numel(cfg.analysis_sensors), numel(prevResult.sensor_eta_id));
    phaseRef.sensor_eta_id(1:n) = prevResult.sensor_eta_id(1:n);
end
phaseRef.weighted_voltage_rmse = prevResult.weighted_voltage_rmse;
phaseRef.plain_voltage_rmse = prevResult.plain_voltage_rmse;
phaseRef.reference_source = 'prev_window_final';
phaseRef.CandidateTable = get_optional_field_local(prevResult, 'CandidateTable', table());
phaseRef.Coverage = get_optional_field_local(prevResult, 'Coverage', struct('clamp_fraction', NaN));
info = make_phase_safe_reference_info_local(phaseRef.reference_source, ...
    'prev_window_quality_pass', prevResult.window_id, true);
end

function [ok, reason] = is_prev_window_phase_reference_usable_local(prevResult, windowId, cfg)
ok = false;
if isempty(prevResult) || ~isstruct(prevResult)
    reason = 'no_previous_window';
    return;
end
if cfg.phase_safe_refresh_every > 0 && mod(windowId - 1, cfg.phase_safe_refresh_every) == 0
    reason = 'scheduled_refresh';
    return;
end
required = {'window_id', 'EO_id', 'A_id', 'phi_id_wrapped', 'dx_c_id', 'weighted_voltage_rmse', 'Coverage'};
for i = 1:numel(required)
    if ~isfield(prevResult, required{i})
        reason = ['previous_missing_', required{i}];
        return;
    end
end
vals = [prevResult.window_id, prevResult.EO_id, prevResult.A_id, ...
    prevResult.phi_id_wrapped, prevResult.dx_c_id, prevResult.weighted_voltage_rmse];
if any(~isfinite(vals))
    reason = 'previous_nonfinite_reference';
    return;
end
if isfield(prevResult.Coverage, 'clamp_fraction') && ...
        prevResult.Coverage.clamp_fraction > cfg.phase_safe_prev_max_clamp_fraction
    reason = 'previous_clamped';
    return;
end
gapRatio = final_candidate_gap_ratio_local(prevResult);
if isfinite(cfg.phase_safe_prev_min_final_gap_ratio) && ...
        gapRatio < cfg.phase_safe_prev_min_final_gap_ratio
    reason = 'previous_final_gap_too_small';
    return;
end
ok = true;
reason = 'prev_window_quality_pass';
end

function gapRatio = final_candidate_gap_ratio_local(result)
gapRatio = inf;
if ~isfield(result, 'CandidateTable') || ~istable(result.CandidateTable) || ...
        height(result.CandidateTable) < 2 || ...
        ~ismember('weighted_voltage_rmse', result.CandidateTable.Properties.VariableNames)
    return;
end
rmse = result.CandidateTable.weighted_voltage_rmse;
rmse = rmse(isfinite(rmse));
if numel(rmse) < 2
    return;
end
rmse = sort(rmse(:), 'ascend');
gapRatio = (rmse(2) - rmse(1)) / max(rmse(1), eps);
end

function [rotFreqHz, rotRpm] = estimate_window_speed_local(WaveformLibrary, bladeId, lapRange, sensorIds)
periods = [];
for sid = sensorIds(:).'
    Sdata = WaveformLibrary.Blade(bladeId).Sensor([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid);
    if isempty(Sdata) || numel(Sdata.Lap) < max(lapRange)
        continue;
    end
    tPeak = arrayfun(@(x) x.t_peak, Sdata.Lap(lapRange));
    dt = diff(tPeak(:));
    periods = [periods; dt(isfinite(dt) & dt > eps)]; %#ok<AGROW>
end
if isempty(periods)
    error('Could not estimate rotation speed for B%d laps %s.', bladeId, mat2str(lapRange));
end
rotFreqHz = 1 / median(periods, 'omitnan');
rotRpm = 60 * rotFreqHz;
end

function eoCandidates = build_eo_candidates_local(rotFreqHz, freqSearchHz, eoPad)
rotFreqHz = max(rotFreqHz, eps);
freqLo = min(freqSearchHz);
freqHi = max(freqSearchHz);
eoMin = max(1, ceil(freqLo / rotFreqHz) - max(0, floor(eoPad)));
eoMax = max(eoMin, floor(freqHi / rotFreqHz) + max(0, floor(eoPad)));
eoAll = eoMin:eoMax;
freqAll = eoAll .* rotFreqHz;
eoCandidates = eoAll(freqAll >= freqLo & freqAll <= freqHi);
if isempty(eoCandidates)
    eoCenter = max(1, round(mean([freqLo, freqHi]) / rotFreqHz));
    eoCandidates = max(1, eoCenter - eoPad):max(1, eoCenter + eoPad);
end
eoCandidates = unique(round(eoCandidates(:).'));
end

function seedTable = solve_template_seed_eo_scan_local(bundle, eoCandidates)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'weighted_voltage_rmse', inf, 'plain_voltage_rmse', inf, ...
    'linear_vp_rmse', inf), numel(eoCandidates), 1);
nSensor = numel(bundle.sensor_ids);
finite = isfinite(bundle.V) & isfinite(bundle.F0) & isfinite(bundle.Fx) & ...
    isfinite(bundle.W) & isfinite(bundle.Theta);
if nnz(finite) < 3
    error('Too few finite points for VP seed scan.');
end
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo * bundle.Theta(finite));
    c1 = cos(eo * bundle.Theta(finite));
    basis = [-bundle.Fx(finite), -bundle.Fx(finite) .* s1, -bundle.Fx(finite) .* c1];
    wSqrt = sqrt(bundle.W(finite));
    y = bundle.V(finite) - bundle.F0(finite);
    coeff = (basis .* wSqrt) \ (y .* wSqrt);
    dxC = coeff(1);
    eta = zeros(1, nSensor);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y - basis * coeff;
    [obj, plainRmse] = template_synchronous_objective_local([A, phi, dxC, eta], eo, bundle);
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx_c = dxC;
    rows(i).sensor_eta = eta;
    rows(i).sensor_eta_max_abs = max(abs(eta), [], 'omitnan');
    rows(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    rows(i).plain_voltage_rmse = plainRmse;
    rows(i).linear_vp_rmse = sqrt(sum(bundle.W(finite) .* linRes.^2) / max(sum(bundle.W(finite)), eps));
end
score = [[rows.weighted_voltage_rmse].', [rows.linear_vp_rmse].', [rows.plain_voltage_rmse].', [rows.EO].'];
[~, order] = sortrows(score, [1, 2, 3, 4]);
seedTable = rows(order);
end

function eoKeep = select_vp_topk_eo_local(seedTable, topK)
if nargin < 2 || isempty(topK) || ~isfinite(topK)
    eoKeep = unique(round([seedTable.EO]), 'stable');
    return;
end
topK = min(max(1, floor(topK)), numel(seedTable));
eoKeep = unique(round([seedTable(1:topK).EO]), 'stable');
end

function result = refine_template_waveform_fit_local(bundle, seedTable, eoCandidates, cfg)
ampLimit = abs(cfg.amplitude_limit_mm);
dxCLimit = abs(cfg.dx_c_limit_mm);
sensorEtaLimit = resolve_sensor_eta_limit_vector_local(bundle, cfg);
nSensor = numel(bundle.sensor_ids);
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'fn_hz', NaN, 'objective', inf, 'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, 'clamp_fraction', NaN, ...
    'max_query_overshoot_mm', NaN, 'overshoot_penalty', NaN, ...
    'V_pred', [], 'u_est', []), numel(eoCandidates), 1);
best = struct('objective', inf);

for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find([seedTable.EO] == eo, 1, 'first');
    if isempty(seedIdx)
        seed = seedTable(1);
    else
        seed = seedTable(seedIdx);
    end
    eta0 = zeros(1, nSensor);
    if isfield(seed, 'sensor_eta') && ~isempty(seed.sensor_eta)
        eta0(1:min(nSensor, numel(seed.sensor_eta))) = seed.sensor_eta(1:min(nSensor, numel(seed.sensor_eta)));
    end
    eta0(:) = 0;
    seedParams = [seed.A, seed.phi, seed.dx_c, eta0];
    fun = @(p) template_synchronous_objective_local(p, eo, bundle);
    opts = optimset('Display', 'off', 'MaxIter', cfg.optim_max_iter, 'MaxFunEvals', cfg.optim_max_fun_evals);
    pOpt = fminsearch(@(p) bounded_synchronous_objective_local( ...
        p, fun, ampLimit, dxCLimit, sensorEtaLimit, nSensor), seedParams, opts);
    pOpt(1) = min(abs(pOpt(1)), ampLimit);
    pOpt(2) = wrap_to_pi_local(pOpt(2));
    pOpt(3) = max(min(pOpt(3), dxCLimit), -dxCLimit);
    pOpt = bound_sensor_eta_params_local(pOpt, sensorEtaLimit, nSensor);
    [obj, plainRmse, vPred, uEst, coverage] = template_synchronous_objective_local(pOpt, eo, bundle);

    candidate(i).EO = eo;
    candidate(i).A = pOpt(1);
    candidate(i).phi = pOpt(2);
    candidate(i).dx_c = pOpt(3);
    candidate(i).sensor_eta = pOpt(4:3+nSensor);
    candidate(i).sensor_eta_max_abs = max(abs(candidate(i).sensor_eta), [], 'omitnan');
    candidate(i).fn_hz = eo * bundle.rot_freq_mean_hz;
    candidate(i).objective = obj;
    candidate(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    candidate(i).plain_voltage_rmse = plainRmse;
    candidate(i).clamp_fraction = coverage.clamp_fraction;
    candidate(i).max_query_overshoot_mm = coverage.max_query_overshoot_mm;
    candidate(i).overshoot_penalty = coverage.overshoot_penalty;
    candidate(i).V_pred = vPred;
    candidate(i).u_est = uEst;
    if obj < best.objective
        best = candidate(i);
    end
end

candidateTable = struct2table(rmfield(candidate, {'V_pred', 'u_est', 'sensor_eta'}));
candidateTable = sortrows(candidateTable, {'weighted_voltage_rmse', 'plain_voltage_rmse', 'EO'}, ...
    {'ascend', 'ascend', 'ascend'});

result = struct();
result.A_id = best.A;
result.EO_id = best.EO;
result.phi_id_wrapped = best.phi;
result.dx_c_id = best.dx_c;
result.d0_id = best.dx_c;
result.sensor_eta_id = best.sensor_eta;
result.sensor_eta_max_abs_mm = best.sensor_eta_max_abs;
result.fn_id = best.fn_hz;
result.weighted_voltage_rmse = best.weighted_voltage_rmse;
result.plain_voltage_rmse = best.plain_voltage_rmse;
result.V_pred = best.V_pred;
result.u_est = best.u_est;
result.objective_score = best.objective;
result.valid_segment_count = bundle.valid_segment_count;
result.point_count = bundle.point_count;
result.CandidateTable = candidateTable;
result.Coverage = struct( ...
    'clamp_fraction', best.clamp_fraction, ...
    'max_query_overshoot_mm', best.max_query_overshoot_mm, ...
    'query_guard_mm', bundle.query_guard_mm, ...
    'overshoot_penalty', best.overshoot_penalty, ...
    'overshoot_penalty_weight', bundle.overshoot_penalty_weight);
result.metric_note = 'Template-only VP top-K synchronous waveform fit with shared dx_c; no EO prior.';
result.sensor_eta_init_mm = get_optional_field_local(bundle, 'eta_initial_guess_by_sensor', []);
result.sensor_eta_init_iqr_mm = get_optional_field_local(bundle, 'eta_initial_iqr_by_sensor', []);
end

function [bundleOut, info] = expand_bundle_by_phase_safe_local( ...
        WaveformLibrary, bladeId, lapRange, bundleIn, phaseRef, phaseRefInfo, cfg)
bundleOut = bundleIn;
info = struct( ...
    'applied', false, ...
    'status', 'identified_by_vp_topk_synchronous_template_core_only', ...
    'pass_count', 0, ...
    'initial_point_count', bundleIn.point_count, ...
    'final_point_count', bundleIn.point_count, ...
    'reason', phaseRefInfo.reason, ...
    'reference_source', phaseRefInfo.source, ...
    'chosen_pass', 'core');
if ~cfg.phase_safe_expansion_enable || isempty(phaseRef)
    return;
end

bundleTry = build_phase_safe_expanded_bundle_from_raw_waveform_local( ...
    WaveformLibrary, bladeId, lapRange, bundleIn, phaseRef, cfg);
if isempty(bundleTry) || ~isfield(bundleTry, 'point_count') || ...
        bundleTry.point_count < max(cfg.min_point_count, cfg.phase_safe_min_point_count)
    info.reason = 'phase_safe_bundle_too_small';
    info.status = 'identified_by_vp_topk_synchronous_template_core_only';
    return;
end
if bundleTry.point_count <= bundleIn.point_count
    info.reason = 'phase_safe_added_no_extra_points';
    info.status = 'identified_by_vp_topk_synchronous_template_core_only';
    return;
end

bundleOut = bundleTry;
info.applied = true;
info.status = 'identified_by_vp_topk_synchronous_template_phase_safe_expanded';
info.pass_count = 1;
info.final_point_count = bundleTry.point_count;
info.reason = 'phase_safe_previous_window_reference_used';
end

function bundle = build_phase_safe_expanded_bundle_from_raw_waveform_local( ...
        WaveformLibrary, bladeId, lapRange, bundleRef, phaseRef, cfg)
X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
F0 = [];
Fx = [];
sensorIndex = [];
sensorMeta = repmat(empty_sensor_meta_local(), numel(bundleRef.sensor_ids), 1);

for is = 1:numel(bundleRef.sensor_ids)
    sid = bundleRef.sensor_ids(is);
    Sdata = WaveformLibrary.Blade(bladeId).Sensor([WaveformLibrary.Blade(bladeId).Sensor.sensor_id] == sid);
    if isempty(Sdata)
        continue;
    end
    Tpl = template_from_sensor_library_local(Sdata);
    [tRaw, xRaw, vRaw, thetaRaw] = concatenate_raw_laps_local(Sdata.Lap(lapRange));
    [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        select_phase_safe_sensor_points_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, phaseRef, is, cfg, sid);
    if numel(tSel) < cfg.phase_safe_min_point_count_per_sensor
        continue;
    end
    meta.sensor_id = sid;
    sensorMeta(is) = meta;

    X = [X; xSel(:)]; %#ok<AGROW>
    T = [T; tSel(:)]; %#ok<AGROW>
    V = [V; vSel(:)]; %#ok<AGROW>
    S = [S; repmat(sid, numel(tSel), 1)]; %#ok<AGROW>
    W = [W; wSel(:)]; %#ok<AGROW>
    Theta = [Theta; thetaSel(:)]; %#ok<AGROW>
    F0 = [F0; f0Sel(:)]; %#ok<AGROW>
    Fx = [Fx; fxSel(:)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, numel(tSel), 1)]; %#ok<AGROW>
end

if isempty(T)
    bundle = struct();
    return;
end

bundle = bundleRef;
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.sensor_index = sensorIndex;
bundle.window_meta = sensorMeta;
bundle.WaveformSensor = WaveformLibrary.Blade(bladeId).Sensor;
bundle.valid_segment_count = sum([sensorMeta.pulse_segment_count]);
bundle.point_count = numel(T);
bundle.time_window = [min(T), max(T)];
bundle.selection_pass = 'step06_prev_window_phase_safe_expansion';
end

function Tpl = template_from_sensor_library_local(Sdata)
Tpl = struct();
Tpl.x_grid = Sdata.template_x(:);
Tpl.v_grid = Sdata.template_v(:);
Tpl.dv_dx = Sdata.template_dv_dx(:);
Tpl.x_domain = Sdata.template_domain(:).';
Tpl.threshold = NaN;
end

function [t, x, v, theta] = concatenate_raw_laps_local(Lap)
t = [];
x = [];
v = [];
theta = [];
for k = 1:numel(Lap)
    t = [t; Lap(k).t(:)]; %#ok<AGROW>
    x = [x; Lap(k).x_rel(:)]; %#ok<AGROW>
    v = [v; Lap(k).V(:)]; %#ok<AGROW>
    theta = [theta; Lap(k).theta(:)]; %#ok<AGROW>
end
[t, order] = sort(t);
x = x(order);
v = v(order);
theta = theta(order);
end

function [t, x, v, theta, f0, fx] = concatenate_raw_laps_with_core_cache_local(Lap)
[t, x, v, theta] = concatenate_raw_laps_local(Lap);
f0 = [];
fx = [];
tCat = [];
for k = 1:numel(Lap)
    if isfield(Lap(k), 'core_f0') && isfield(Lap(k), 'core_fx') && ...
            numel(Lap(k).core_f0) == numel(Lap(k).t) && numel(Lap(k).core_fx) == numel(Lap(k).t)
        f0 = [f0; Lap(k).core_f0(:)]; %#ok<AGROW>
        fx = [fx; Lap(k).core_fx(:)]; %#ok<AGROW>
        tCat = [tCat; Lap(k).t(:)]; %#ok<AGROW>
    else
        f0 = [];
        fx = [];
        return;
    end
end
[~, order] = sort(tCat);
f0 = f0(order);
fx = fx(order);
end

function [tSel, xSel, vSel, thetaSel, f0Sel, fxSel, wSel, meta] = ...
        select_phase_safe_sensor_points_local(tRaw, xRaw, vRaw, thetaRaw, Tpl, phaseRef, sensorLocalIndex, cfg, sensorId)
meta = empty_sensor_meta_local();
meta.raw_points = numel(tRaw);
meta.domain_selection_mode = cfg.domain_selection_mode;
meta.domain_soft_margin_mm = cfg.domain_soft_margin_mm;
xDomainUse = expand_sensor_domain_local(Tpl.x_domain, cfg, sensorId);
meta.template_domain = xDomainUse;

f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw(:), 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw(:), 'pchip', NaN);
threshold = cfg.default_sensor_threshold;
if isfield(Tpl, 'threshold') && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
end

if strcmpi(cfg.pulse_selection_mode, 'all')
    [maskPulse, pulseSegmentCount] = isolate_all_pulses_local(vRaw, threshold);
else
    [maskPulse, pulseSegmentCount] = isolate_main_pulse_local(vRaw, threshold);
end

maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, cfg);
maskDomain = xRaw >= xDomainUse(1) + cfg.domain_margin_mm & ...
             xRaw <= xDomainUse(2) - cfg.domain_margin_mm;
maskFinite = isfinite(f0) & isfinite(fx0) & isfinite(vRaw) & isfinite(thetaRaw);
maskBase = maskPulse & maskEffective & maskDomain & maskFinite;

queryGuardMM = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, maskBase, cfg, sensorId);
xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
maskQuerySafe = xRaw >= xDomain(1) + queryGuardMM & ...
                xRaw <= xDomain(2) - queryGuardMM;
maskPhaseSafe = build_phase_safe_query_mask_local(Tpl, xRaw, vRaw, thetaRaw, xDomain, sensorLocalIndex, phaseRef, cfg);
if nnz(maskBase & maskPhaseSafe) >= cfg.phase_safe_min_point_count_per_sensor
    maskQuerySafe = maskPhaseSafe;
end

mask = maskBase;
if strcmpi(cfg.domain_selection_mode, 'hard')
    safeMask = mask & maskQuerySafe;
    if nnz(safeMask) >= cfg.phase_safe_min_point_count_per_sensor
        mask = safeMask;
    end
end
if nnz(mask) < cfg.phase_safe_min_point_count_per_sensor
    mask = maskBase;
end
if nnz(mask) < cfg.phase_safe_min_point_count_per_sensor
    mask = maskEffective & maskDomain & maskFinite;
    if strcmpi(cfg.domain_selection_mode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= cfg.phase_safe_min_point_count_per_sensor
            mask = safeMask;
        end
    end
end
if nnz(mask) < cfg.phase_safe_min_point_count_per_sensor
    mask = maskPulse & maskFinite;
end

tSel = tRaw(mask);
xSel = xRaw(mask);
vSel = vRaw(mask);
thetaSel = thetaRaw(mask);
f0Sel = f0(mask);
fxSel = fx0(mask);
wSel = build_filtered_weight_local_refine(tSel, xSel, vSel, Tpl, queryGuardMM, cfg);

if isempty(tSel)
    meta.valid_points = 0;
else
    vStatic = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xSel, 'pchip', NaN);
    meta.valid_points = numel(tSel);
    meta.selected_x_range = [min(xSel), max(xSel)];
    meta.shape_rmse = sqrt(mean((vSel(:) - vStatic(:)).^2, 'omitnan'));
end
meta.pulse_segment_count = pulseSegmentCount;
meta.dynamic_effective_points = nnz(maskEffective);
meta.query_safe_points = nnz(maskQuerySafe(maskBase));
meta.phase_safe_points = nnz(maskBase & maskPhaseSafe);
meta.query_guard_mm = queryGuardMM;
end

function meta = empty_sensor_meta_local()
meta = struct( ...
    'sensor_id', NaN, ...
    'raw_points', 0, ...
    'valid_points', 0, ...
    'eta_initial_guess_mm', NaN, ...
    'eta_initial_iqr_mm', NaN, ...
    'pulse_segment_count', 0, ...
    'dynamic_effective_points', 0, ...
    'query_safe_points', 0, ...
    'phase_safe_points', 0, ...
    'query_guard_mm', NaN, ...
    'domain_selection_mode', '', ...
    'domain_soft_margin_mm', NaN, ...
    'template_domain', [NaN NaN], ...
    'selected_x_range', [NaN NaN], ...
    'shape_rmse', NaN);
end

function wTotal = build_filtered_weight_local_refine(t, xComp, v, Tpl, queryGuardMM, cfg)
wEdge = build_edge_weight_local_refine(t, v, cfg.weight_floor);
wDomain = build_domain_soft_weight_local_refine(xComp, Tpl.x_domain, cfg.domain_soft_margin_mm, cfg.weight_floor);
wQuery = build_query_guard_soft_weight_local_refine(xComp, [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))], ...
    queryGuardMM, cfg.domain_soft_margin_mm, cfg.weight_floor);
wGradient = build_template_gradient_weight_local_refine(Tpl, xComp, cfg.domain_selection_mode, cfg.weight_floor);
wTotal = max(cfg.weight_floor, wEdge(:) .* wDomain(:) .* wQuery(:) .* wGradient(:));
if max(wTotal) > 0
    wTotal = max(cfg.weight_floor, wTotal ./ max(wTotal));
end
end

function [mask, segmentCount] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segmentCount = numel(starts);
end

function [mask, segmentCount] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segmentCount = 1;
    return;
end
segmentGap = idx(starts(2:end)) - idx(ends(1:end-1));
largeGapThreshold = max(10, round(0.02 * numel(v)));
pulseBreaks = find(segmentGap > largeGapThreshold);
groupStarts = [1; pulseBreaks + 1];
groupEnds = [pulseBreaks; numel(starts)];
for ig = 1:numel(groupStarts)
    segIds = groupStarts(ig):groupEnds(ig);
    bestSeg = segIds(1);
    bestPeak = -inf;
    for iseg = segIds
        seg = idx(starts(iseg):ends(iseg));
        peakVal = max(v(seg));
        if peakVal > bestPeak
            bestPeak = peakVal;
            bestSeg = iseg;
        end
    end
    mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end
segmentCount = numel(groupStarts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, cfg)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(cfg.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end

gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gTplMax = max(gTpl(finite), [], 'omitnan');
if ~isfinite(gTplMax) || gTplMax <= 0
    maskTpl = finite;
else
    maskTpl = gTpl >= cfg.dynamic_template_gradient_min_ratio * gTplMax;
end

gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if ~isfinite(gTimeMax) || gTimeMax <= 0
    maskTime = false(size(v));
else
    maskTime = gTime >= cfg.dynamic_time_gradient_min_ratio * gTimeMax;
end
peakLevel = prctile(v(finite), min(max(cfg.dynamic_peak_quantile, 0), 100));
maskPeak = v >= max(threshold, peakLevel);
mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function queryGuardMM = resolve_query_guard_mm_local(Tpl, x, v, baseMask, cfg, sensorId)
queryGuardMM = cfg.query_guard_mm;
if ~strcmpi(cfg.query_guard_mode, 'adaptive') || nnz(baseMask) < 8
    queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, cfg, sensorId);
    return;
end
xSel = x(baseMask);
vSel = v(baseMask);
xStat = invert_template_voltage_local_refine(Tpl, vSel, xSel);
uApp = abs(xSel(:) - xStat(:));
uApp = uApp(isfinite(uApp));
if isempty(uApp)
    return;
end
adaptiveGuard = prctile(uApp, min(max(cfg.query_guard_quantile, 0), 100)) + ...
    cfg.query_guard_safety_mm;
queryGuardMM = min(max(adaptiveGuard, cfg.query_guard_min_mm), cfg.query_guard_max_mm);
queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, cfg, sensorId);
end

function queryGuardMM = apply_sensor_query_guard_scale_local(queryGuardMM, cfg, sensorId)
if ~isfield(cfg, 'sensor_query_guard_scale_table') || isempty(cfg.sensor_query_guard_scale_table)
    return;
end
tableNow = cfg.sensor_query_guard_scale_table;
if size(tableNow, 2) < 2
    return;
end
idx = find(tableNow(:, 1) == sensorId, 1);
if isempty(idx)
    return;
end
scale = tableNow(idx, 2);
if ~isfinite(scale) || scale <= 0
    return;
end
queryGuardMM = queryGuardMM * scale;
if isfield(cfg, 'query_guard_min_mm') && isfinite(cfg.query_guard_min_mm)
    queryGuardMM = max(queryGuardMM, 0.5 * cfg.query_guard_min_mm);
end
if isfield(cfg, 'query_guard_max_mm') && isfinite(cfg.query_guard_max_mm)
    queryGuardMM = min(queryGuardMM, cfg.query_guard_max_mm);
end
end

function xDomainUse = expand_sensor_domain_local(xDomain, cfg, sensorId)
xDomainUse = xDomain;
if ~isfield(cfg, 'sensor_domain_expand_mm_table') || isempty(cfg.sensor_domain_expand_mm_table)
    return;
end
tableNow = cfg.sensor_domain_expand_mm_table;
if size(tableNow, 2) < 2
    return;
end
idx = find(tableNow(:, 1) == sensorId, 1);
if isempty(idx)
    return;
end
expandMM = tableNow(idx, 2);
if ~isfinite(expandMM) || abs(expandMM) <= eps
    return;
end
xCandidate = [xDomain(1) - expandMM, xDomain(2) + expandMM];
if xCandidate(2) <= xCandidate(1)
    return;
end
xDomainUse = xCandidate;
end

function xStat = invert_template_voltage_local_refine(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = nan(size(v));
for i = 1:numel(v)
    vv = v(i);
    if ~isfinite(vv)
        continue;
    end
    diffV = vGrid - vv;
    crossingX = [];
    for k = 1:(numel(vGrid) - 1)
        if diffV(k) == 0
            crossingX(end + 1, 1) = xGrid(k); %#ok<AGROW>
        elseif diffV(k) * diffV(k + 1) <= 0
            denom = vGrid(k + 1) - vGrid(k);
            if abs(denom) < eps
                continue;
            end
            alpha = (vv - vGrid(k)) / denom;
            crossingX(end + 1, 1) = xGrid(k) + alpha * (xGrid(k + 1) - xGrid(k)); %#ok<AGROW>
        end
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function [etaGuess, etaIqr] = estimate_sensor_eta_initial_guess_local(Tpl, xSel, vSel)
etaGuess = NaN;
etaIqr = NaN;
if isempty(Tpl) || ~isfield(Tpl, 'x_grid') || ~isfield(Tpl, 'v_grid')
    return;
end
xSel = xSel(:);
vSel = vSel(:);
keep = isfinite(xSel) & isfinite(vSel);
if nnz(keep) < 5
    return;
end
xGrid = Tpl.x_grid(:);
xDomain = [min(xGrid), max(xGrid)];
keep = keep & xSel >= xDomain(1) & xSel <= xDomain(2);
if nnz(keep) < 5
    return;
end
xTpl = invert_template_voltage_local_refine(Tpl, vSel(keep), xSel(keep));
etaApp = xSel(keep) - xTpl;
etaApp = etaApp(isfinite(etaApp));
if isempty(etaApp)
    return;
end
etaGuess = median(etaApp, 'omitnan');
etaIqr = iqr(etaApp);
end

function etaLimitBySensor = build_sensor_eta_limit_by_sensor_local(bundle, cfg)
nSensor = numel(bundle.sensor_ids);
etaLimitBySensor = abs(cfg.sensor_eta_limit_mm) * ones(1, nSensor);
if ~isfield(cfg, 'sensor_eta_adaptive_limit_enable') || ~cfg.sensor_eta_adaptive_limit_enable
    return;
end
minMM = max(cfg.sensor_eta_adaptive_min_mm, 0);
maxMM = abs(cfg.sensor_eta_adaptive_max_mm);
for is = 1:nSensor
    etaLim = NaN;
    if isfield(bundle, 'WaveformSensor') && numel(bundle.WaveformSensor) >= is
        sensorNow = bundle.WaveformSensor(is);
        if isfield(sensorNow, 'template_eta_limit_mm') && isfinite(sensorNow.template_eta_limit_mm) && sensorNow.template_eta_limit_mm > 0
            etaLim = sensorNow.template_eta_limit_mm;
        end
    end
    if ~isfinite(etaLim) || etaLim <= 0
        etaLimitBySensor(is) = min(max(abs(cfg.sensor_eta_limit_mm), minMM), maxMM);
        continue;
    end
    etaLimitBySensor(is) = min(max(etaLim, minMM), maxMM);
end
if ~isempty(etaLimitBySensor)
    etaLimitBySensor(1) = min(max(etaLimitBySensor(1), minMM), maxMM);
end
end

function sensorEtaLimit = resolve_sensor_eta_limit_vector_local(bundle, cfg)
if isfield(bundle, 'sensor_eta_limit_by_sensor') && ~isempty(bundle.sensor_eta_limit_by_sensor)
    sensorEtaLimit = abs(bundle.sensor_eta_limit_by_sensor(:).');
else
    sensorEtaLimit = abs(cfg.sensor_eta_limit_mm) * ones(1, numel(bundle.sensor_ids));
end
end

function etaVec = anchor_sensor_eta_vector_local(etaVec)
etaVec = etaVec(:).';
if isempty(etaVec)
    etaVec = [];
    return;
end
finite = isfinite(etaVec);
if ~any(finite)
    etaVec = zeros(size(etaVec));
    return;
end
ref = NaN;
if isfinite(etaVec(1))
    ref = etaVec(1);
else
    refIdx = find(finite, 1, 'first');
    if ~isempty(refIdx)
        ref = etaVec(refIdx);
    end
end
if ~isfinite(ref)
    ref = 0;
end
etaVec(~finite) = 0;
etaVec = etaVec - ref;
etaVec(1) = 0;
end

function maskSafe = build_phase_safe_query_mask_local(Tpl, x, v, theta, xDomain, sensorLocalIndex, phaseRef, cfg)
if isfield(phaseRef, 'reference_source') && strcmpi(phaseRef.reference_source, 'template_inverse_apparent_query')
    maskSafe = build_template_inverse_safe_mask_local(Tpl, x, v, xDomain, cfg);
    return;
end

A = phaseRef.A_id;
eo = phaseRef.EO_id;
phi = phaseRef.phi_id_wrapped;
dxC = phaseRef.dx_c_id;
eta = 0;
if isfield(phaseRef, 'sensor_eta_id') && numel(phaseRef.sensor_eta_id) >= sensorLocalIndex
    eta = phaseRef.sensor_eta_id(sensorLocalIndex);
end
uEst = A .* sin(eo .* theta(:) + phi);
xQuery = x(:) - dxC - eta - uEst;
margin = cfg.phase_safe_margin_mm;
maskSafe = xQuery >= xDomain(1) + margin & xQuery <= xDomain(2) - margin;
maskSafe = reshape(maskSafe, size(x));
end

function maskSafe = build_template_inverse_safe_mask_local(Tpl, x, v, xDomain, cfg)
xStat = invert_template_voltage_local_refine(Tpl, v(:), x(:));
margin = cfg.phase_safe_margin_mm;
uApp = abs(x(:) - xStat(:));
maskSafe = isfinite(xStat) & ...
    xStat >= xDomain(1) + margin & xStat <= xDomain(2) - margin & ...
    uApp <= cfg.amplitude_limit_mm + abs(cfg.dx_c_limit_mm) + margin;
maskSafe = reshape(maskSafe, size(x));
end

function value = get_optional_field_local(S, fieldName, defaultValue)
if isfield(S, fieldName)
    value = S.(fieldName);
else
    value = defaultValue;
end
end

function wEdge = build_edge_weight_local_refine(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_local_refine(x, xDomain, marginMM, floorW)
if marginMM <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMM, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_local_refine(x, xDomain, queryGuardMM, marginMM, floorW)
if queryGuardMM <= 0 || marginMM <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMM - marginMM, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMM, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_local_refine(Tpl, x, domainSelectionMode, floorW)
if ~strcmpi(domainSelectionMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end

function obj = bounded_synchronous_objective_local(p, fun, ampLimit, dxCLimit, sensorEtaLimit, nSensor)
pUse = p;
pUse(1) = min(abs(pUse(1)), ampLimit);
pUse(2) = wrap_to_pi_local(pUse(2));
pUse(3) = max(min(pUse(3), dxCLimit), -dxCLimit);
pUse = bound_sensor_eta_params_local(pUse, sensorEtaLimit, nSensor);
obj = fun(pUse);
end

function p = bound_sensor_eta_params_local(p, sensorEtaLimit, nSensor)
need = 3 + nSensor;
if numel(p) < need
    p(numel(p)+1:need) = 0;
end
p = p(1:need);
eta = p(4:end);
sensorEtaLimit = abs(sensorEtaLimit(:).');
if isempty(sensorEtaLimit)
    sensorEtaLimit = inf(1, nSensor);
elseif numel(sensorEtaLimit) == 1
    sensorEtaLimit = repmat(sensorEtaLimit, 1, nSensor);
elseif numel(sensorEtaLimit) < nSensor
    sensorEtaLimit(numel(sensorEtaLimit)+1:nSensor) = sensorEtaLimit(end);
else
    sensorEtaLimit = sensorEtaLimit(1:nSensor);
end
eta = max(min(eta, sensorEtaLimit), -sensorEtaLimit);
if ~isempty(eta)
    eta = eta - eta(1);
    eta = max(min(eta, sensorEtaLimit), -sensorEtaLimit);
end
p(4:end) = eta;
end

function [obj, plainRmse, vPred, uEst, coverage] = template_synchronous_objective_local(params, eo, bundle)
A = params(1);
phi = params(2);
dxC = params(3);
sensorEtaLimit = resolve_sensor_eta_limit_vector_local(bundle, struct('sensor_eta_limit_mm', bundle.sensor_eta_limit_mm));
params = bound_sensor_eta_params_local(params, sensorEtaLimit, numel(bundle.sensor_ids));
eta = params(4:end).';
uEst = A .* sin(eo .* bundle.Theta + phi);
etaSample = eta(bundle.sensor_index(:));
xIn = bundle.X - dxC - etaSample - uEst;
vPred = nan(size(bundle.V));
clamped = false(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    xLo = bundle.x_domain_by_sensor(is, 1);
    xHi = bundle.x_domain_by_sensor(is, 2);
    xEvalRaw = xIn(mask);
    xEval = min(max(xEvalRaw, xLo), xHi);
    clamped(mask) = xEvalRaw < xLo | xEvalRaw > xHi;
    overshoot(mask) = max(xLo - xEvalRaw, 0) + max(xEvalRaw - xHi, 0);
    vPred(mask) = bundle.interp_v{is}(xEval);
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
res = bundle.V(valid) - vPred(valid);
overshootPenalty = bundle.overshoot_penalty_weight * sum(bundle.W(valid) .* (overshoot(valid) .^ 2));
etaRegPenalty = bundle.point_count * (bundle.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta .^ 2))) ^ 2;
obj = sum(bundle.W(valid) .* (res .^ 2)) + overshootPenalty + etaRegPenalty + 1e6 * nnz(~valid);
plainRmse = sqrt(mean(res .^ 2));
coverage = struct();
coverage.clamp_fraction = nnz(clamped) / max(numel(clamped), 1);
coverage.max_query_overshoot_mm = max(overshoot, [], 'omitnan');
coverage.overshoot_penalty = overshootPenalty;
coverage.sensor_eta_reg_penalty = etaRegPenalty;
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function summary = build_resonance_summary_local(T)
summary = struct();
summary.dominant_eo = mode(T.EO);
summary.mean_frequency_hz = mean(T.FrequencyHz, 'omitnan');
summary.median_frequency_hz = median(T.FrequencyHz, 'omitnan');
summary.mean_amplitude_mm = mean(T.AmplitudeMM, 'omitnan');
summary.median_amplitude_mm = median(T.AmplitudeMM, 'omitnan');
summary.mean_weighted_rmse = mean(T.WeightedVoltageRMSE, 'omitnan');
summary.median_weighted_rmse = median(T.WeightedVoltageRMSE, 'omitnan');
end

function visualize_identification_trend_local(T, P)
figDir = fullfile(P.view.figureDir, '06_identification');
if P.view.saveFigures && exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end
fig = figure('Name', 'Step06 direct-template identification trend', 'Color', 'w', ...
    'Position', [80, 80, 1250, 760]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for bladeId = unique(T.BladeID).'
    rows = T(T.BladeID == bladeId, :);
    plot(rows.WindowID, rows.EO, 'o-', 'DisplayName', sprintf('B%d', bladeId));
end
ylabel('EO');
title('Selected EO');
legend('Location', 'bestoutside');

nexttile; hold on; grid on; box on;
for bladeId = unique(T.BladeID).'
    rows = T(T.BladeID == bladeId, :);
    plot(rows.WindowID, rows.AmplitudeMM, 'o-', 'DisplayName', sprintf('B%d', bladeId));
end
ylabel('A (mm)');
title('Amplitude');

nexttile; hold on; grid on; box on;
for bladeId = unique(T.BladeID).'
    rows = T(T.BladeID == bladeId, :);
    plot(rows.WindowID, rows.WeightedVoltageRMSE, 'o-', 'DisplayName', sprintf('B%d', bladeId));
end
xlabel('Window ID');
ylabel('Weighted RMSE (V)');
title('Fit quality');

if P.view.saveFigures
    if isfield(P.files, 'identificationFigure') && ~isempty(P.files.identificationFigure)
        figPath = P.files.identificationFigure;
    else
        figPath = fullfile(figDir, 'Step06_DirectTemplateIdentification_20241106.png');
    end
    exportgraphics(fig, figPath, 'Resolution', 300);
end
end
