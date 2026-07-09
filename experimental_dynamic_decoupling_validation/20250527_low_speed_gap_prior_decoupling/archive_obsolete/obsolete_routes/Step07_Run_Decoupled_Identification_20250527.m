%% Step07: Main multi-sensor clearance-vibration decoupled identification
% This step identifies clearance and vibration with a hybrid low-speed/gap
% library model:
%
%   V_s(x,t) = a_s { T_low,s(x - x0_s - u(t))
%              + F(g_s, x - x0_s - u(t)) - F(g0_s, x - x0_s - u(t)) } + b_s
%
% Step06B first fits the low-speed no-vibration waveform to estimate g0_s.
% Dynamic fitting then combines the measured low-speed waveform T_low,s with
% the static gap library correction F(g_s)-F(g0_s).

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
sensorOverride = strtrim(getenv('STEP07_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
cfg.eoCandidates = 3:24;
cfg.eoScanOffsetsHz = [-0.5, 0, 0.5];
cfg.eoRelRmseTolerance = 1.02;
cfg.frequencyFineOffsetsHz = -0.50:0.10:0.50;
cfg.eoLockWindowCount = 3;

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
gapPriorFile = fullfile(outDir, sprintf('Step06B_LowSpeed_GapPrior_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(highMapFile)
    error('Run Step06 first. Missing file: %s', highMapFile);
end
if ~isfile(gapPriorFile)
    error('Run Step06B first. Missing file: %s', gapPriorFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
H = load(highMapFile, 'highMap');
highMapFull = H.highMap;
[highMapFull, coordinateSource] = attach_model_coordinate_local(highMapFull);
P = load(gapPriorFile, 'GapPrior', 'priorTable');
GapPrior = validate_gap_prior_local(P.GapPrior, cfg.analysisSensors, responseSurface);
LowTemplate = load_low_speed_template_for_sensors_local(GapPrior, cfg.analysisSensors);
windowSpecs = build_sliding_window_specs_local(highMapFull);
maxWindowOverride = str2double(strtrim(getenv('STEP07_MAX_WINDOWS')));
if isfinite(maxWindowOverride) && maxWindowOverride > 0
    windowSpecs = windowSpecs(1:min(numel(windowSpecs), floor(maxWindowOverride)));
end

fitSettings = struct();
fitSettings.maxFitPoints = 350;
fitSettings.x0LimitMm = 0.08;
fitSettings.ampLimitMm = 1.50;
fitSettings.apparentGapGridCount = 17;
fitSettings.finalGapGridCount = 3;
fitSettings.staticX0GridMm = -0.06:0.03:0.06;
fitSettings.x0GridMm = 0;
fitSettings.coordinateIterations = 1;
fitSettings.outerCoordinatePasses = 1;
fitSettings.fullCoordinateFrequencyKeep = 3;
fitSettings.x0ProjectionLimitMm = 0.06;
fitSettings.lowSpeedDeltaGapHalfWidthMm = 0.35;
fitSettings.lowSpeedDeltaGapGridCount = 17;
fitSettings.lowSpeedDeltaGapFinalHalfWidthMm = 0.08;
fitSettings.lowSpeedDeltaGapFinalGridCount = 5;
fitSettings.lowSpeedFullWaveTopK = numel(cfg.eoCandidates);
fitSettings.lowSpeedFullWaveAmpLimitMm = 0.60;
fitSettings.lowSpeedFullWaveDxLimitMm = 0.20;
fitSettings.lowSpeedFullWaveDeltaGapHalfWidthMm = 0.12;
fitSettings.lowSpeedFullWaveDeltaGapGridCount = 7;
fitSettings.lowSpeedFullWaveOvershootPenalty = 1e5;

fprintf('\n=== Step07: sliding-window clearance-vibration decoupling ===\n');
fprintf('HighMap source: %s\n', highMapFile);
fprintf('Gap-prior source: %s\n', gapPriorFile);
fprintf('Analysis sensors: %s, target blade B%d\n', mat2str(cfg.analysisSensors), cfg.targetBlade);
fprintf('Sliding windows: %d windows, %d laps/window, step %d lap(s)\n', ...
    numel(windowSpecs), highMapFull.analysisWinSize, highMapFull.slidingStep);
fprintf('Model coordinate source: %s\n', coordinateSource);
disp(P.priorTable);

[globalEO, eoSelectionTable] = initialize_global_eo_local( ...
    highMapFull, windowSpecs, responseSurface, GapPrior, LowTemplate, cfg, fitSettings);
fprintf('Data-driven global EO initialization: EO%d selected from candidates %s.\n', ...
    globalEO, mat2str(cfg.eoCandidates));
disp(eoSelectionTable);

WindowResult = struct([]);
trendRows = cell(numel(windowSpecs), 1);
bestScore = inf;
bestIdx = 1;
for iw = 1:numel(windowSpecs)
    highMapWinFull = subset_highmap_by_laps_local(highMapFull, windowSpecs(iw).lapRange);
    highMapWinFull.window_id = windowSpecs(iw).window_id;
    fprintf('\nWindow %02d/%02d, laps %s, time %.6f-%.6f s\n', ...
        iw, numel(windowSpecs), mat2str(windowSpecs(iw).lapRange), ...
        min(highMapWinFull.t_v), max(highMapWinFull.t_v));
    [windowResult, trendRows{iw}] = fit_one_highmap_window_local( ...
        highMapWinFull, responseSurface, GapPrior, LowTemplate, cfg, fitSettings, coordinateSource, globalEO);
    if iw == 1
        WindowResult = repmat(windowResult, numel(windowSpecs), 1);
    else
        WindowResult(iw) = windowResult;
    end
    fprintf('  joint: f %.6f Hz, EO %.4f, A %.6f mm, g_mean %.6f mm, RMSE %.6f mV\n', ...
        WindowResult(iw).multiSensorJoint.freqHz, ...
        WindowResult(iw).multiSensorJoint.freqHz / highMapWinFull.rotFreqHz, ...
        WindowResult(iw).multiSensorJoint.amplitudeMm, ...
        mean(WindowResult(iw).multiSensorJoint.gBySensorMm, 'omitnan'), ...
        WindowResult(iw).multiSensorJoint.weightedRmseMv);
    if WindowResult(iw).multiSensorJoint.weightedRmseMv < bestScore
        bestScore = WindowResult(iw).multiSensorJoint.weightedRmseMv;
        bestIdx = iw;
    end
end
trendTable = vertcat(trendRows{:});

result = struct();
result.dataset = '20250527';
result.caseName = highMapFull.caseName;
result.targetBlade = highMapFull.targetBlade;
result.analysisSensors = highMapFull.analysisSensors(:).';
result.referenceFreqHz = highMapFull.referenceFreqHz;
result.referenceOrder = globalEO;
result.dataDrivenGlobalEO = globalEO;
result.eoSelectionTable = eoSelectionTable;
result.coordinateUnit = 'mm circumferential displacement';
result.gMinMm = min(responseSurface.gTrainMm);
result.gMaxMm = max(responseSurface.gTrainMm);
result.gapPrior = GapPrior;
result.gapPriorFile = gapPriorFile;
result.x0LimitMm = fitSettings.x0LimitMm;
result.ampLimitMm = fitSettings.ampLimitMm;
result.modelCoordinateSource = coordinateSource;
result.fullHighMapPointCount = highMapFull.pointCount;
result.dynamicBaselineBySensor = highMapFull.dynamicBaselineBySensor;
result.analysisStartTime = highMapFull.analysisStartTime;
result.targetLaps = highMapFull.targetLaps;
result.analysisWinSize = highMapFull.analysisWinSize;
result.slidingStep = highMapFull.slidingStep;
result.WindowResult = WindowResult;
result.Trend = trendTable;
result.BestWindowIndex = bestIdx;
result.BestWindow = WindowResult(bestIdx);
result.multiSensorJoint = result.BestWindow.multiSensorJoint;
result.peakCoordinateDiagnostic = result.BestWindow.peakCoordinateDiagnostic;
result.note = ['Hybrid low-speed-template/gap-library clearance-vibration decoupling. Step06B fits the low-speed ' ...
    'no-vibration waveform to estimate per-sensor effective-clearance centers and bounds. Step07 identifies the dynamic ' ...
    'vibration with the measured low-speed waveform as the baseline template and the static gap response surface as ' ...
    'a clearance-change correction around the low-speed gap prior. The global engine order is selected by integer-EO ' ...
    'variable-projection screening on the first sliding windows; vibration frequency and amplitude are then estimated ' ...
    'from each dynamic window rather than imported from a template-only reference.'];

sensorList = result.analysisSensors;
summaryTable = table(result.BestWindowIndex, ...
    result.multiSensorJoint.freqHz, result.multiSensorJoint.freqHz / result.BestWindow.rotFreqHz, ...
    result.multiSensorJoint.amplitudeMm, ...
    mean(result.multiSensorJoint.gBySensorMm, 'omitnan'), ...
    result.multiSensorJoint.weightedRmseMv, result.peakCoordinateDiagnostic.amplitudeMm, ...
    'VariableNames', {'bestWindowIndex', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv', 'peakCoordinateAmplitudeMm'});
sensorTable = make_sensor_result_table_local(result.multiSensorJoint, sensorList);

matFile = fullfile(outDir, sprintf('Step07_Main_Decoupled_Identification_20250527_B%d_%s.mat', ...
    result.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step07_Main_Decoupled_Identification_Summary_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag));
sensorCsvFile = fullfile(outDir, sprintf('Step07_Main_MultiSensor_Gap_Summary_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag));
save(matFile, 'result', 'summaryTable', 'sensorTable', '-v7.3');
writetable(summaryTable, csvFile);
writetable(sensorTable, sensorCsvFile);
writetable(trendTable, fullfile(outDir, sprintf('Step07_Main_Decoupled_Identification_Trend_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag)));

fprintf('\nStep07 complete.\n');
disp(summaryTable);
disp(sensorTable);
disp(trendTable);
fprintf('Result saved to:\n  %s\n', matFile);

%% Local functions
function windowSpecs = build_sliding_window_specs_local(highMap)
laps = unique(highMap.rev_v(:).');
laps = sort(laps(isfinite(laps)));
if isfield(highMap, 'analysisWinSize') && isfinite(highMap.analysisWinSize)
    winSize = highMap.analysisWinSize;
else
    winSize = numel(laps);
end
if isfield(highMap, 'slidingStep') && isfinite(highMap.slidingStep)
    step = highMap.slidingStep;
else
    step = winSize;
end
numWindows = floor((numel(laps) - winSize) / step) + 1;
if numWindows < 1
    error('Not enough laps (%d) for a %d-lap sliding window.', numel(laps), winSize);
end
windowSpecs = repmat(struct('window_id', NaN, 'lapRange', []), numWindows, 1);
for iw = 1:numWindows
    i0 = 1 + (iw - 1) * step;
    windowSpecs(iw).window_id = iw;
    windowSpecs(iw).lapRange = laps(i0:(i0 + winSize - 1));
end
end

function highMapWin = subset_highmap_by_laps_local(highMap, lapRange)
highMapWin = highMap;
n = numel(highMap.t_v);
mask = ismember(highMap.rev_v, lapRange);
if nnz(mask) < 20
    error('Sliding window laps %s have too few highMap points.', mat2str(lapRange));
end
fields = fieldnames(highMap);
for i = 1:numel(fields)
    f = fields{i};
    v = highMap.(f);
    if isnumeric(v) || islogical(v)
        if isvector(v) && numel(v) == n
            highMapWin.(f) = v(mask);
        end
    end
end
highMapWin.pointCount = nnz(mask);
highMapWin.lapRange = lapRange;
highMapWin.window_id = NaN;
if isfield(highMap, 't_v')
    highMapWin.time_window = [min(highMapWin.t_v), max(highMapWin.t_v)];
end
localRotFreqHz = compute_local_rot_freq_from_highmap_local(highMap, highMapWin.time_window);
if isfinite(localRotFreqHz)
    highMapWin.rotFreqHz = localRotFreqHz;
    highMapWin.rotRpm = 60 * localRotFreqHz;
end
end

function rotFreqHz = compute_local_rot_freq_from_highmap_local(highMap, timeWindow)
rotFreqHz = NaN;
if ~isfield(highMap, 'experiment_case') || ~isfield(highMap.experiment_case, 'Global_OPR_Times')
    return;
end
oprTimes = highMap.experiment_case.Global_OPR_Times(:);
if isfield(highMap.experiment_case, 'Config') && isfield(highMap.experiment_case.Config, 'opr_pulses_per_rev')
    pulsesPerRev = highMap.experiment_case.Config.opr_pulses_per_rev;
elseif isfield(highMap.experiment_case, 'Diagnostics') && isfield(highMap.experiment_case.Diagnostics, 'opr_pulses_per_rev')
    pulsesPerRev = highMap.experiment_case.Diagnostics.opr_pulses_per_rev;
else
    return;
end
mask = oprTimes >= timeWindow(1) & oprTimes <= timeWindow(2);
t = oprTimes(mask);
if numel(t) > pulsesPerRev
    rotFreqHz = median(1 ./ max(t(1+pulsesPerRev:end) - t(1:end-pulsesPerRev), eps), 'omitnan');
end
end

function [windowResult, trendRow] = fit_one_highmap_window_local(highMapWinFull, responseSurface, GapPrior, LowTemplate, cfg, fitSettings, coordinateSource, globalEO)
highMapFullPointCount = highMapWinFull.pointCount;
highMap = decimate_highmap_for_fit_local(highMapWinFull, fitSettings.maxFitPoints);
sensorList = highMap.analysisSensors(:).';

gapBounds = get_gap_bounds_for_sensors_local(GapPrior, sensorList);
gMin = gapBounds.lower(:).';
gMax = gapBounds.upper(:).';
gCenter = gapBounds.center(:).';
x0LimitMm = fitSettings.x0LimitMm;
lowSpeedPriorState = estimate_low_speed_prior_static_state_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, fitSettings);
freqCandidates = build_frequency_candidates_local(highMap, cfg, globalEO);
apparentGap = estimate_apparent_gap_windows_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gMin, gMax, gCenter, ...
    fitSettings, fitSettings.staticX0GridMm);
x0Grid = make_local_x0_grid_from_static_alignment_local(apparentGap, x0LimitMm);
apparentGap = refine_gap_range_with_projected_diagnostic_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, apparentGap, ...
    gMin, gMax, freqCandidates);
gGridBySensor = make_final_gap_grid_from_apparent_local(apparentGap, gMin, gMax, fitSettings.finalGapGridCount);

bestJoint = fit_projection_coordinate_descent_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    gGridBySensor, x0Grid, freqCandidates, fitSettings, apparentGap);
bestJoint.diagInfo.gapGridBySensor = gGridBySensor;
bestJoint.diagInfo.gapAtGridBoundary = check_gap_grid_boundary_local(bestJoint.gBySensorMm, gGridBySensor);
bestJoint.diagInfo.frequencyCandidatesHz = freqCandidates(:).';
referenceDiagnostics = build_reference_model_diagnostics_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, bestJoint, fitSettings);
bestJoint.diagInfo.referenceModelDiagnostics = referenceDiagnostics;
lowSpeedPriorVpBranch = fit_low_speed_prior_branch_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, ...
    lowSpeedPriorState, cfg, fitSettings);
lowSpeedPriorBranch = fit_low_speed_prior_full_wave_branch_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, ...
    lowSpeedPriorState, cfg, fitSettings);

windowResult = struct();
windowResult.window_id = highMapWinFull.window_id;
windowResult.lapRange = highMapWinFull.lapRange;
windowResult.timeWindow = highMapWinFull.time_window;
windowResult.coordinateSource = coordinateSource;
windowResult.frequencyCandidatesHz = freqCandidates(:);
windowResult.fullPointCount = highMapFullPointCount;
windowResult.fitPointCount = highMap.pointCount;
windowResult.gapPriorCenterMm = gCenter(:);
windowResult.gapPriorLowerMm = gMin(:);
windowResult.gapPriorUpperMm = gMax(:);
windowResult.lowSpeedPriorStaticState = lowSpeedPriorState;
windowResult.apparentGapDiagnostic = apparentGap;
windowResult.multiSensorJoint = pack_result_multisensor_local(bestJoint.objective, ...
    bestJoint.params, bestJoint.freqHz, bestJoint.pred, bestJoint.affine, ...
    bestJoint.diagInfo, highMap, sensorList);
windowResult.referenceModelDiagnostics = referenceDiagnostics;
windowResult.lowSpeedPriorVpBranch = lowSpeedPriorVpBranch;
windowResult.lowSpeedPriorBranch = lowSpeedPriorBranch;
windowResult.sensorTable = make_sensor_result_table_local(windowResult.multiSensorJoint, sensorList);
windowResult.peakCoordinateDiagnostic = fit_peak_coordinate_diagnostic_local( ...
    highMapWinFull, windowResult.multiSensorJoint.freqHz, sensorList);
windowResult.rotFreqHz = highMapWinFull.rotFreqHz;
windowResult.rotRpm = highMapWinFull.rotRpm;

trendRow = table(windowResult.window_id, windowResult.lapRange(1), windowResult.lapRange(end), ...
    mean(windowResult.timeWindow), highMapWinFull.rotFreqHz, highMapWinFull.rotRpm, ...
    windowResult.multiSensorJoint.freqHz, windowResult.multiSensorJoint.freqHz / highMapWinFull.rotFreqHz, ...
    windowResult.multiSensorJoint.amplitudeMm, ...
    mean(windowResult.multiSensorJoint.gBySensorMm, 'omitnan'), ...
    windowResult.multiSensorJoint.weightedRmseMv, ...
    referenceDiagnostics.templateOnlyWeightedRmseMv, ...
    referenceDiagnostics.fixedGapHybridWeightedRmseMv, ...
    referenceDiagnostics.templateToFreeRmseDropMv, ...
    referenceDiagnostics.fixedGapToFreeRmseDropMv, ...
    lowSpeedPriorBranch.deltaGapMm, ...
    mean(lowSpeedPriorBranch.gBySensorMm, 'omitnan'), ...
    lowSpeedPriorBranch.freqHz, ...
    lowSpeedPriorBranch.freqHz / highMapWinFull.rotFreqHz, ...
    lowSpeedPriorBranch.amplitudeMm, ...
    lowSpeedPriorBranch.weightedRmseMv, ...
    windowResult.peakCoordinateDiagnostic.amplitudeMm, highMap.pointCount, ...
    'VariableNames', {'window_id', 'lap_start', 'lap_end', 'window_center_time', ...
    'rot_freq_hz', 'rot_rpm', 'frequency_hz', 'EO', 'amplitude_mm', ...
    'mean_gap_mm', 'weighted_rmse_mV', 'template_only_rmse_mV', ...
    'fixed_gap_hybrid_rmse_mV', 'template_to_free_rmse_drop_mV', ...
    'fixed_gap_to_free_rmse_drop_mV', 'low_prior_delta_gap_mm', 'low_prior_mean_gap_mm', ...
    'low_prior_frequency_hz', 'low_prior_EO', 'low_prior_amplitude_mm', ...
    'low_prior_weighted_rmse_mV', 'peak_coordinate_amplitude_mm', 'point_count'});
end

function [highMap, coordinateSource] = attach_model_coordinate_local(highMap)
if isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(highMap.t_v) && ...
        any(isfinite(highMap.x_fit_v))
    highMap.x_model_v = highMap.x_fit_v(:);
    coordinateSource = 'x_fit_v_opr_restored_not_peak_centered';
else
    highMap.x_model_v = highMap.x_v(:);
    coordinateSource = 'x_v_legacy_peak_centered_fallback';
    warning(['highMap.x_fit_v is missing; Step07 falls back to peak-centered x_v. ' ...
        'Vibration amplitude may be underestimated. Re-run Step06 to generate x_fit_v.']);
end
end

function freqCandidates = build_frequency_candidates_local(highMap, cfg, globalEO)
if ~isfield(highMap, 'rotFreqHz') || ~isfinite(highMap.rotFreqHz) || highMap.rotFreqHz <= 0
    error('Cannot build Step07 frequency candidates without a valid local rotating frequency.');
end
freqCenter = globalEO * highMap.rotFreqHz;
offsets = cfg.frequencyFineOffsetsHz(:);
if ~isfinite(freqCenter) || freqCenter <= 0
    error('Cannot build Step07 frequency candidates without a valid EO-derived frequency center.');
end
freqCandidates = freqCenter + offsets;
freqCandidates = unique(freqCandidates(isfinite(freqCandidates) & freqCandidates > 0));
freqCandidates = sort(freqCandidates(:).');
end

function GapPrior = validate_gap_prior_local(GapPrior, sensorList, responseSurface)
if ~isfield(GapPrior, 'sensorPrior') || isempty(GapPrior.sensorPrior)
    error('GapPrior.sensorPrior is missing or empty. Re-run Step06B.');
end
available = [GapPrior.sensorPrior.sensorId];
missing = setdiff(sensorList, available);
if ~isempty(missing)
    error('Gap prior is missing sensor(s): %s. Re-run Step06B with the same sensor set.', mat2str(missing));
end
gTrainMin = min(responseSurface.gTrainMm);
gTrainMax = max(responseSurface.gTrainMm);
for i = 1:numel(GapPrior.sensorPrior)
    sp = GapPrior.sensorPrior(i);
    if ~isfinite(sp.gCenterMm) || ~isfinite(sp.gLowerMm) || ~isfinite(sp.gUpperMm)
        error('Invalid gap prior for CH%d. Re-run Step06B.', sp.sensorId);
    end
    GapPrior.sensorPrior(i).gLowerMm = max(gTrainMin, sp.gLowerMm);
    GapPrior.sensorPrior(i).gUpperMm = min(gTrainMax, sp.gUpperMm);
    if GapPrior.sensorPrior(i).gLowerMm > GapPrior.sensorPrior(i).gUpperMm
        error('Invalid gap-prior bound for CH%d.', sp.sensorId);
    end
end
end

function gapBounds = get_gap_bounds_for_sensors_local(GapPrior, sensorList)
gapBounds = struct();
gapBounds.center = nan(numel(sensorList), 1);
gapBounds.lower = nan(numel(sensorList), 1);
gapBounds.upper = nan(numel(sensorList), 1);
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = find([GapPrior.sensorPrior.sensorId] == sid, 1, 'first');
    if isempty(idx)
        error('Gap prior is missing CH%d.', sid);
    end
    sp = GapPrior.sensorPrior(idx);
    gapBounds.center(is) = min(max(sp.gCenterMm, sp.gLowerMm), sp.gUpperMm);
    gapBounds.lower(is) = sp.gLowerMm;
    gapBounds.upper(is) = sp.gUpperMm;
end
end

function LowTemplate = load_low_speed_template_for_sensors_local(GapPrior, sensorList)
if ~isfield(GapPrior, 'templateFile') || ~isfile(GapPrior.templateFile)
    error('GapPrior.templateFile is missing or invalid. Re-run Step06B.');
end
S = load(GapPrior.templateFile, 'Template');
Template = S.Template;
LowTemplate = struct();
LowTemplate.templateFile = GapPrior.templateFile;
LowTemplate.sensor = repmat(struct('sensorId', NaN, 'x', [], 'vMv', [], ...
    'gCenterMm', NaN, 'responseX0Mm', NaN, 'responseGain', NaN), numel(sensorList), 1);
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
    if isempty(idx)
        error('Low-speed template does not contain CH%d.', sid);
    end
    Tsen = Template.Sensor(idx);
    priorIdx = find([GapPrior.sensorPrior.sensorId] == sid, 1, 'first');
    if isempty(priorIdx)
        error('Gap prior does not contain CH%d.', sid);
    end
    sp = GapPrior.sensorPrior(priorIdx);
    LowTemplate.sensor(is).sensorId = sid;
    LowTemplate.sensor(is).x = Tsen.x_grid(:);
    LowTemplate.sensor(is).vMv = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
    LowTemplate.sensor(is).gCenterMm = sp.gCenterMm;
    LowTemplate.sensor(is).responseX0Mm = sp.x0Mm;
    LowTemplate.sensor(is).responseGain = sp.affineGain;
end
end

function [globalEO, eoSelectionTable] = initialize_global_eo_local( ...
    highMapFull, windowSpecs, responseSurface, GapPrior, LowTemplate, cfg, fitSettings)
sensorList = highMapFull.analysisSensors(:).';
% EO is initialized with the measured low-speed waveform only. This follows
% Step03's first-order displacement projection and avoids letting the static
% gap surface choose a wrong EO before vibration frequency is locked.
unused = {responseSurface, GapPrior}; %#ok<NASGU>
nLock = min([cfg.eoLockWindowCount, numel(windowSpecs)]);
nEO = numel(cfg.eoCandidates);
scoreByEO = nan(nEO, nLock);
bestFreqByEO = nan(nEO, nLock);
rotFreqByWindow = nan(1, nLock);
for iw = 1:nLock
    highMapWin = subset_highmap_by_laps_local(highMapFull, windowSpecs(iw).lapRange);
    highMap = decimate_highmap_for_fit_local(highMapWin, fitSettings.maxFitPoints);
    rotFreqByWindow(iw) = highMap.rotFreqHz;
    for ie = 1:nEO
        eo = cfg.eoCandidates(ie);
        freqCandidates = eo * highMap.rotFreqHz + cfg.eoScanOffsetsHz(:).';
        bestScore = inf;
        bestFreq = NaN;
        for f = freqCandidates
            if ~isfinite(f) || f <= 0
                continue;
            end
            rmse = score_low_template_first_order_local( ...
                highMap, LowTemplate, sensorList, f, fitSettings);
            if rmse < bestScore
                bestScore = rmse;
                bestFreq = f;
            end
        end
        scoreByEO(ie, iw) = bestScore;
        bestFreqByEO(ie, iw) = bestFreq;
    end
end
relScoreByEO = nan(size(scoreByEO));
for iw = 1:nLock
    windowBest = min(scoreByEO(:, iw), [], 'omitnan');
    for ie = 1:nEO
        eo = cfg.eoCandidates(ie);
        fErr = abs(bestFreqByEO(ie, iw) - eo * rotFreqByWindow(iw));
        relScoreByEO(ie, iw) = scoreByEO(ie, iw) / max(windowBest, eps) + 1e-3 * fErr;
    end
end
supportCount = sum(relScoreByEO <= cfg.eoRelRmseTolerance, 2, 'omitnan');
medianRelScore = median(relScoreByEO, 2, 'omitnan');
meanRelScore = mean(relScoreByEO, 2, 'omitnan');
selectionScore = medianRelScore + 0.02 * (nLock - supportCount);
medianScore = median(scoreByEO, 2, 'omitnan');
meanBestFreq = mean(bestFreqByEO, 2, 'omitnan');
eoSelectionTable = table(cfg.eoCandidates(:), supportCount(:), medianRelScore(:), ...
    meanRelScore(:), selectionScore(:), medianScore(:), meanBestFreq(:), ...
    'VariableNames', {'EO', 'support_count', 'median_rel_rmse', 'mean_rel_rmse', ...
    'score', 'medianWeightedRmseMv', 'meanBestFrequencyHz'});
eoSelectionTable = sortrows(eoSelectionTable, {'support_count', 'median_rel_rmse'}, ...
    {'descend', 'ascend'});
globalEO = eoSelectionTable.EO(1);
end

function rmse = score_low_template_first_order_local(highMap, LowTemplate, sensorList, freqHz, fitSettings)
[tAll, xAll, vAll, wAll, sAll, f0All, fxAll] = pack_low_template_window_local( ...
    highMap, LowTemplate, sensorList);
if numel(tAll) < 30
    rmse = inf;
    return;
end
tRel = tAll - min(tAll);
wAll = wAll ./ max(wAll);
swAll = sqrt(max(wAll, eps));

vNorm = nan(size(vAll));
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = sAll == sid;
    if nnz(idx) < 3
        continue;
    end
    Haff = [ones(nnz(idx), 1), f0All(idx)];
    betaAff = (Haff .* swAll(idx)) \ (vAll(idx) .* swAll(idx));
    gain0 = betaAff(2);
    if abs(gain0) < 1e-6
        gain0 = 1;
    end
    vNorm(idx) = (vAll(idx) - betaAff(1)) ./ gain0;
end
dv = vNorm - f0All;
valid = isfinite(dv) & isfinite(fxAll) & isfinite(wAll);
if nnz(valid) < 30
    rmse = inf;
    return;
end

s1 = sin(2*pi*freqHz*tRel);
c1 = cos(2*pi*freqHz*tRel);
Hx0 = zeros(numel(tRel), numel(sensorList));
for is = 1:numel(sensorList)
    Hx0(:, is) = -fxAll .* (sAll == sensorList(is));
end
H = [Hx0, -fxAll .* s1, -fxAll .* c1];
Hw = H(valid, :) .* swAll(valid);
yw = dv(valid) .* swAll(valid);
G = Hw.' * Hw;
rhs = Hw.' * yw;
ridge = 1e-8 * max(trace(G), eps);
beta = (G + ridge * eye(size(G))) \ rhs;

x0 = min(max(beta(1:numel(sensorList)).', -fitSettings.x0LimitMm), fitSettings.x0LimitMm);
c = beta(end-1);
d = beta(end);
amp = hypot(c, d);
if amp > fitSettings.ampLimitMm
    scale = fitSettings.ampLimitMm / max(amp, eps);
    c = c * scale;
    d = d * scale;
end
u = c .* s1 + d .* c1;
rAll = [];
wRes = [];
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = sAll == sid;
    T = get_low_template_sensor_local(LowTemplate, sid);
    xShift = xAll(idx) - x0(is) - u(idx);
    fShift = interp1(T.x, T.vMv, xShift, 'pchip', NaN);
    keep = isfinite(fShift) & isfinite(vAll(idx)) & isfinite(wAll(idx));
    if nnz(keep) < 20
        continue;
    end
    Haff = [ones(nnz(keep), 1), fShift(keep)];
    ww = wAll(idx);
    vv = vAll(idx);
    betaAff = (Haff .* sqrt(max(ww(keep), eps))) \ (vv(keep) .* sqrt(max(ww(keep), eps)));
    rr = vv(keep) - Haff * betaAff;
    rAll = [rAll; rr]; %#ok<AGROW>
    wRes = [wRes; ww(keep)]; %#ok<AGROW>
end
if isempty(rAll)
    rmse = inf;
else
    rmse = sqrt(sum(wRes .* rAll.^2) / sum(wRes));
end
end

function diagnostics = build_reference_model_diagnostics_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, bestJoint, fitSettings)
diagnostics = struct();
diagnostics.description = ['Side-branch diagnostics only. These scores do not affect Step07 main identification; ' ...
    'they quantify how much the free-gap hybrid result improves over fixed low-speed-template references.'];
diagnostics.templateOnlyWeightedRmseMv = NaN;
diagnostics.templateOnlyAmplitudeMm = NaN;
diagnostics.templateOnlyFrequencyHz = bestJoint.freqHz;
diagnostics.fixedGapHybridWeightedRmseMv = NaN;
diagnostics.fixedGapHybridAmplitudeMm = NaN;
diagnostics.fixedGapHybridFrequencyHz = bestJoint.freqHz;
diagnostics.freeGapHybridWeightedRmseMv = bestJoint.objective;
diagnostics.templateToFreeRmseDropMv = NaN;
diagnostics.fixedGapToFreeRmseDropMv = NaN;
diagnostics.fixedGapBySensorMm = gCenter(:);
diagnostics.freeGapBySensorMm = bestJoint.gBySensorMm(:);
diagnostics.gapShiftFromLowSpeedMm = bestJoint.gBySensorMm(:) - gCenter(:);

try
    templateOnly = evaluate_template_only_reference_local( ...
        highMap, LowTemplate, sensorList, bestJoint.freqHz, fitSettings);
    diagnostics.templateOnlyWeightedRmseMv = templateOnly.weightedRmseMv;
    diagnostics.templateOnlyAmplitudeMm = templateOnly.amplitudeMm;
    diagnostics.templateOnlyX0BySensorMm = templateOnly.x0BySensorMm(:);
catch ME
    diagnostics.templateOnlyError = ME.message;
end

try
    fixedGap = evaluate_projection_state_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        gCenter(:).', bestJoint.x0BySensorMm, bestJoint.freqHz, fitSettings);
    diagnostics.fixedGapHybridWeightedRmseMv = fixedGap.objective;
    diagnostics.fixedGapHybridAmplitudeMm = fixedGap.diagInfo.amplitudeMm;
    diagnostics.fixedGapHybridX0BySensorMm = fixedGap.x0BySensorMm(:);
catch ME
    diagnostics.fixedGapHybridError = ME.message;
end

diagnostics.templateToFreeRmseDropMv = diagnostics.templateOnlyWeightedRmseMv - ...
    diagnostics.freeGapHybridWeightedRmseMv;
diagnostics.fixedGapToFreeRmseDropMv = diagnostics.fixedGapHybridWeightedRmseMv - ...
    diagnostics.freeGapHybridWeightedRmseMv;
end

function result = evaluate_template_only_reference_local(highMap, LowTemplate, sensorList, freqHz, fitSettings)
[tAll, xAll, vAll, wAll, sAll, f0All, fxAll] = pack_low_template_window_local( ...
    highMap, LowTemplate, sensorList);
if numel(tAll) < 30
    result = struct('weightedRmseMv', inf, 'amplitudeMm', NaN, ...
        'cMm', NaN, 'dMm', NaN, 'x0BySensorMm', nan(numel(sensorList), 1), ...
        'pointCount', numel(tAll));
    return;
end
tRel = tAll - min(tAll);
wAll = wAll ./ max(wAll);
swAll = sqrt(max(wAll, eps));

vNorm = nan(size(vAll));
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = sAll == sid;
    if nnz(idx) < 3
        continue;
    end
    Haff = [ones(nnz(idx), 1), f0All(idx)];
    betaAff = (Haff .* swAll(idx)) \ (vAll(idx) .* swAll(idx));
    gain0 = betaAff(2);
    if abs(gain0) < 1e-6
        gain0 = 1;
    end
    vNorm(idx) = (vAll(idx) - betaAff(1)) ./ gain0;
end
dv = vNorm - f0All;
valid = isfinite(dv) & isfinite(fxAll) & isfinite(wAll);
if nnz(valid) < 30
    result = struct('weightedRmseMv', inf, 'amplitudeMm', NaN, ...
        'cMm', NaN, 'dMm', NaN, 'x0BySensorMm', nan(numel(sensorList), 1), ...
        'pointCount', nnz(valid));
    return;
end

s1 = sin(2*pi*freqHz*tRel);
c1 = cos(2*pi*freqHz*tRel);
Hx0 = zeros(numel(tRel), numel(sensorList));
for is = 1:numel(sensorList)
    Hx0(:, is) = -fxAll .* (sAll == sensorList(is));
end
H = [Hx0, -fxAll .* s1, -fxAll .* c1];
Hw = H(valid, :) .* swAll(valid);
yw = dv(valid) .* swAll(valid);
G = Hw.' * Hw;
rhs = Hw.' * yw;
ridge = 1e-8 * max(trace(G), eps);
beta = (G + ridge * eye(size(G))) \ rhs;

x0 = min(max(beta(1:numel(sensorList)).', -fitSettings.x0LimitMm), fitSettings.x0LimitMm);
c = beta(end-1);
d = beta(end);
amp = hypot(c, d);
if amp > fitSettings.ampLimitMm
    scale = fitSettings.ampLimitMm / max(amp, eps);
    c = c * scale;
    d = d * scale;
end
u = c .* s1 + d .* c1;
[weightedRmse, plainRmse, pointCount] = score_template_only_residual_local( ...
    xAll, vAll, wAll, sAll, u, x0, LowTemplate, sensorList);
result = struct('weightedRmseMv', weightedRmse, 'plainRmseMv', plainRmse, ...
    'amplitudeMm', hypot(c, d), 'cMm', c, 'dMm', d, ...
    'x0BySensorMm', x0(:), 'pointCount', pointCount);
end

function [weightedRmse, plainRmse, pointCount] = score_template_only_residual_local( ...
    xAll, vAll, wAll, sAll, u, x0BySensor, LowTemplate, sensorList)
rAll = [];
wRes = [];
for is = 1:numel(sensorList)
    sid = sensorList(is);
    idx = sAll == sid;
    T = get_low_template_sensor_local(LowTemplate, sid);
    xShift = xAll(idx) - x0BySensor(is) - u(idx);
    fShift = interp1(T.x, T.vMv, xShift, 'pchip', NaN);
    keep = isfinite(fShift) & isfinite(vAll(idx)) & isfinite(wAll(idx));
    if nnz(keep) < 20
        continue;
    end
    Haff = [ones(nnz(keep), 1), fShift(keep)];
    ww = wAll(idx);
    vv = vAll(idx);
    sw = sqrt(max(ww(keep), eps));
    betaAff = (Haff .* sw) \ (vv(keep) .* sw);
    rr = vv(keep) - Haff * betaAff;
    rAll = [rAll; rr]; %#ok<AGROW>
    wRes = [wRes; ww(keep)]; %#ok<AGROW>
end
pointCount = numel(rAll);
if isempty(rAll)
    weightedRmse = inf;
    plainRmse = inf;
else
    weightedRmse = sqrt(sum(wRes .* rAll.^2) / sum(wRes));
    plainRmse = sqrt(mean(rAll.^2, 'omitnan'));
end
end

function [tAll, xAll, vAll, wAll, sAll, f0All, fxAll] = pack_low_template_window_local(highMap, LowTemplate, sensorList)
tAll = [];
xAll = [];
vAll = [];
wAll = [];
sAll = [];
f0All = [];
fxAll = [];
for is = 1:numel(sensorList)
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    if nnz(mask) < 3
        continue;
    end
    T = get_low_template_sensor_local(LowTemplate, sid);
    x = highMap.x_model_v(mask);
    f0 = interp1(T.x, T.vMv, x, 'pchip', NaN);
    fx = finite_difference_interp_local(T.x, T.vMv, x);
    keep = isfinite(f0) & isfinite(fx) & isfinite(highMap.V_a(mask)) & ...
        isfinite(highMap.t_v(mask)) & isfinite(highMap.W_v(mask));
    t = highMap.t_v(mask);
    v = highMap.V_a(mask);
    w = highMap.W_v(mask);
    tAll = [tAll; t(keep)]; %#ok<AGROW>
    xAll = [xAll; x(keep)]; %#ok<AGROW>
    vAll = [vAll; v(keep)]; %#ok<AGROW>
    wAll = [wAll; max(w(keep), 0.05)]; %#ok<AGROW>
    sAll = [sAll; repmat(sid, nnz(keep), 1)]; %#ok<AGROW>
    f0All = [f0All; f0(keep)]; %#ok<AGROW>
    fxAll = [fxAll; fx(keep)]; %#ok<AGROW>
end
end

function apparentGap = estimate_apparent_gap_windows_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gMin, gMax, gCenter, fitSettings, x0Grid)
nSensor = numel(sensorList);
apparentGap = repmat(struct('sensorId', NaN, 'gAppMm', NaN, 'gLowerMm', NaN, ...
    'gUpperMm', NaN, 'x0AppMm', NaN, 'rmseMinMv', NaN, 'rmseAtPriorMv', NaN, ...
    'projectedGapScaleMm', NaN, 'separability', NaN, 'gridMm', [], ...
    'rmseCurveMv', [], 'projectedDiagnostic', []), nSensor, 1);
for is = 1:nSensor
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    gGrid = linspace(gMin(is), gMax(is), fitSettings.apparentGapGridCount);
    rmseCurve = nan(size(gGrid));
    bestX0ByGap = nan(size(gGrid));
    for ig = 1:numel(gGrid)
        bestRmse = inf;
        bestX0 = 0;
        for ix0 = 1:numel(x0Grid)
            raw = eval_anchored_template_local( ...
                responseSurface, LowTemplate, sid, gGrid(ig), gCenter(is), ...
                highMap.x_model_v(mask) - x0Grid(ix0));
            [pred, ~] = apply_single_sensor_offset_local(raw, highMap.V_a(mask), highMap.W_v(mask));
            valid = isfinite(pred) & isfinite(highMap.V_a(mask)) & isfinite(highMap.W_v(mask));
            if nnz(valid) >= 20
                res = highMap.V_a(mask) - pred;
                wLocal = highMap.W_v(mask);
                rmseTry = sqrt(sum(wLocal(valid) .* res(valid).^2) / sum(wLocal(valid)));
                if rmseTry < bestRmse
                    bestRmse = rmseTry;
                    bestX0 = x0Grid(ix0);
                end
            end
        end
        rmseCurve(ig) = bestRmse;
        bestX0ByGap(ig) = bestX0;
    end
    [rmseMin, bestIdx] = min(rmseCurve, [], 'omitnan');
    if isempty(bestIdx) || ~isfinite(rmseMin)
        bestIdx = round((numel(gGrid) + 1) / 2);
        rmseMin = NaN;
    end
    priorIdx = find(abs(gGrid - gCenter(is)) == min(abs(gGrid - gCenter(is))), 1, 'first');
    threshold = rmseMin + max(rmseMin * 0.03, 0.5);
    good = isfinite(rmseCurve) & rmseCurve <= threshold;
    if nnz(good) < 2
        expand = max(1, round(0.12 * numel(gGrid)));
        good(max(1, bestIdx-expand):min(numel(gGrid), bestIdx+expand)) = true;
    end
    apparentGap(is).sensorId = sid;
    apparentGap(is).gAppMm = gGrid(bestIdx);
    apparentGap(is).x0AppMm = bestX0ByGap(bestIdx);
    apparentGap(is).gLowerMm = max(gMin(is), min(gGrid(good)));
    apparentGap(is).gUpperMm = min(gMax(is), max(gGrid(good)));
    apparentGap(is).rmseMinMv = rmseMin;
    apparentGap(is).rmseAtPriorMv = rmseCurve(priorIdx);
    apparentGap(is).projectedGapScaleMm = 0.5 * (apparentGap(is).gUpperMm - apparentGap(is).gLowerMm);
    apparentGap(is).separability = abs(apparentGap(is).rmseAtPriorMv - rmseMin) / max(rmseMin, eps);
    apparentGap(is).gridMm = gGrid;
    apparentGap(is).rmseCurveMv = rmseCurve;
end
end

function gGridBySensor = make_final_gap_grid_from_apparent_local(apparentGap, gMin, gMax, nGrid)
gGridBySensor = cell(numel(apparentGap), 1);
for i = 1:numel(apparentGap)
    lo = max(gMin(i), apparentGap(i).gLowerMm);
    hi = min(gMax(i), apparentGap(i).gUpperMm);
    if ~isfinite(lo) || ~isfinite(hi) || lo > hi
        lo = gMin(i);
        hi = gMax(i);
    end
    if nGrid <= 2
        if abs(apparentGap(i).gAppMm - lo) >= abs(apparentGap(i).gAppMm - hi)
            edge = lo;
        else
            edge = hi;
        end
        gGrid = [apparentGap(i).gAppMm, edge];
    else
        gGrid = linspace(lo, hi, nGrid);
    end
    gGridBySensor{i} = unique(sort(gGrid));
end
end

function x0Grid = make_local_x0_grid_from_static_alignment_local(apparentGap, x0LimitMm)
x0Center = [apparentGap.x0AppMm];
x0Center(~isfinite(x0Center)) = 0;
x0Grid = cell(numel(apparentGap), 1);
for i = 1:numel(apparentGap)
    x0Grid{i} = min(max(x0Center(i), -x0LimitMm), x0LimitMm);
end
end

function staticState = estimate_low_speed_prior_static_state_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, fitSettings)
dgGrid = linspace(-fitSettings.lowSpeedDeltaGapHalfWidthMm, ...
    fitSettings.lowSpeedDeltaGapHalfWidthMm, fitSettings.lowSpeedDeltaGapGridCount);
x0Grid = fitSettings.staticX0GridMm(:).';

scoreByDg = inf(numel(dgGrid), 1);
x0ByDg = nan(numel(dgGrid), numel(sensorList));
for idg = 1:numel(dgGrid)
    dg = dgGrid(idg);
    gTry = min(max(gCenter + dg, gMin), gMax);
    [scoreByDg(idg), x0ByDg(idg, :)] = score_static_shared_delta_gap_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, gTry, x0Grid);
end
[bestScore, bestIdx] = min(scoreByDg);
if isempty(bestIdx) || ~isfinite(bestScore)
    bestIdx = ceil(numel(dgGrid) / 2);
end
deltaGap = dgGrid(bestIdx);
gBySensor = min(max(gCenter + deltaGap, gMin), gMax);
x0BySensor = x0ByDg(bestIdx, :);
x0BySensor(~isfinite(x0BySensor)) = 0;

staticState = struct();
staticState.method = 'low_speed_gap_prior_shared_delta_gap_static_fit';
staticState.deltaGapMm = deltaGap;
staticState.gBySensorMm = gBySensor(:);
staticState.x0BySensorMm = x0BySensor(:);
staticState.weightedRmseMv = bestScore;
staticState.deltaGapGridMm = dgGrid(:);
staticState.rmseCurveMv = scoreByDg(:);
staticState.x0ByDeltaGapMm = x0ByDg;
end

function [score, x0BestBySensor] = score_static_shared_delta_gap_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, x0Grid)
rAll = [];
wAll = [];
x0BestBySensor = nan(1, numel(sensorList));
for is = 1:numel(sensorList)
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    bestScore = inf;
    bestX0 = 0;
    bestResidual = [];
    bestWeight = [];
    for ix0 = 1:numel(x0Grid)
        x0 = x0Grid(ix0);
        raw = eval_anchored_template_local( ...
            responseSurface, LowTemplate, sid, gBySensor(is), gCenter(is), ...
            highMap.x_model_v(mask) - x0);
        [pred, ~] = apply_single_sensor_offset_local(raw, highMap.V_a(mask), highMap.W_v(mask));
        valid = isfinite(pred) & isfinite(highMap.V_a(mask)) & isfinite(highMap.W_v(mask));
        if nnz(valid) < 20
            continue;
        end
        res = highMap.V_a(mask) - pred;
        w = highMap.W_v(mask);
        scoreTry = sqrt(sum(w(valid) .* res(valid).^2) / sum(w(valid)));
        if scoreTry < bestScore
            bestScore = scoreTry;
            bestX0 = x0;
            bestResidual = res(valid);
            bestWeight = w(valid);
        end
    end
    x0BestBySensor(is) = bestX0;
    rAll = [rAll; bestResidual(:)]; %#ok<AGROW>
    wAll = [wAll; bestWeight(:)]; %#ok<AGROW>
end
if isempty(rAll)
    score = inf;
else
    score = sqrt(sum(wAll .* rAll.^2) / sum(wAll));
end
end

function branch = fit_low_speed_prior_branch_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, ...
    staticState, cfg, fitSettings)
freqCandidates = build_frequency_candidates_from_eo_scan_local(highMap, cfg);
scoreByFreq = inf(numel(freqCandidates), 1);
candidateByFreq = cell(numel(freqCandidates), 1);
for ifr = 1:numel(freqCandidates)
    candidate = evaluate_projection_state_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        staticState.gBySensorMm(:).', staticState.x0BySensorMm(:).', ...
        freqCandidates(ifr), fitSettings);
    scoreByFreq(ifr) = candidate.objective;
    candidateByFreq{ifr} = candidate;
end
[~, order] = sort(scoreByFreq, 'ascend');
keepN = min(fitSettings.fullCoordinateFrequencyKeep, numel(order));
best = struct('objective', inf);
bestDg = staticState.deltaGapMm;
dgCenters = staticState.deltaGapMm;
for ii = 1:keepN
    ifr = order(ii);
    if isempty(candidateByFreq{ifr})
        continue;
    end
    freqHz = freqCandidates(ifr);
    localDgGrid = make_local_delta_gap_grid_local( ...
        dgCenters, fitSettings.lowSpeedDeltaGapFinalHalfWidthMm, ...
        fitSettings.lowSpeedDeltaGapFinalGridCount, gCenter, gMin, gMax);
    for idg = 1:numel(localDgGrid)
        dg = localDgGrid(idg);
        gTry = min(max(gCenter + dg, gMin), gMax);
        candidate = evaluate_projection_state_local( ...
            highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
            gTry, staticState.x0BySensorMm(:).', freqHz, fitSettings);
        if candidate.objective < best.objective
            best = candidate;
            bestDg = dg;
        end
    end
end
if ~isfield(best, 'objective') || ~isfinite(best.objective)
    best = candidateByFreq{order(1)};
    bestDg = staticState.deltaGapMm;
end
branch = pack_low_speed_prior_branch_local(best, bestDg, staticState, freqCandidates, scoreByFreq, highMap);
end

function freqCandidates = build_frequency_candidates_from_eo_scan_local(highMap, cfg)
freqList = [];
for eo = cfg.eoCandidates(:).'
    freqList = [freqList, eo * highMap.rotFreqHz + cfg.eoScanOffsetsHz(:).']; %#ok<AGROW>
end
freqCandidates = unique(freqList(isfinite(freqList) & freqList > 0));
end

function dgGrid = make_local_delta_gap_grid_local(centerDg, halfWidth, nGrid, gCenter, gMin, gMax)
lo = max(centerDg - halfWidth, max(gMin - gCenter));
hi = min(centerDg + halfWidth, min(gMax - gCenter));
if ~isfinite(lo) || ~isfinite(hi) || lo > hi
    dgGrid = centerDg;
else
    dgGrid = linspace(lo, hi, max(1, nGrid));
end
end

function branch = pack_low_speed_prior_branch_local(best, deltaGap, staticState, freqCandidates, scoreByFreq, highMap)
branch = struct();
branch.method = 'low_speed_prior_shared_delta_gap_fixed_gap_vp';
branch.deltaGapMm = deltaGap;
branch.gBySensorMm = best.gBySensorMm(:);
branch.x0BySensorMm = best.x0BySensorMm(:);
branch.freqHz = best.freqHz;
branch.EO = best.freqHz / highMap.rotFreqHz;
branch.cMm = best.diagInfo.cMm;
branch.dMm = best.diagInfo.dMm;
branch.amplitudeMm = best.diagInfo.amplitudeMm;
branch.weightedRmseMv = best.objective;
branch.plainRmseMv = best.diagInfo.plainRmse;
branch.VPred = best.pred;
branch.affine = best.affine;
branch.diagInfo = best.diagInfo;
branch.staticState = staticState;
branch.frequencyCandidatesHz = freqCandidates(:);
branch.frequencyScoresMv = scoreByFreq(:);
end

function branch = fit_low_speed_prior_full_wave_branch_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gMin, gMax, ...
    staticState, cfg, fitSettings)
dgGrid = make_local_delta_gap_grid_local( ...
    staticState.deltaGapMm, fitSettings.lowSpeedFullWaveDeltaGapHalfWidthMm, ...
    fitSettings.lowSpeedFullWaveDeltaGapGridCount, gCenter, gMin, gMax);
eoCandidates = unique(round(cfg.eoCandidates(:).'));

seedRows = repmat(struct('EO', NaN, 'deltaGapMm', NaN, 'A', NaN, ...
    'phi', NaN, 'dx_c', NaN, 'weightedRmseMv', inf, ...
    'linearVpRmseMv', inf, 'freqHz', NaN), numel(dgGrid) * numel(eoCandidates), 1);
row = 0;
for idg = 1:numel(dgGrid)
    dg = dgGrid(idg);
    gTry = min(max(gCenter + dg, gMin), gMax);
    for eo = eoCandidates
        row = row + 1;
        seed = solve_low_prior_full_wave_seed_local( ...
            highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
            gTry, staticState.x0BySensorMm(:).', eo, fitSettings);
        seedRows(row).EO = eo;
        seedRows(row).deltaGapMm = dg;
        seedRows(row).A = seed.A;
        seedRows(row).phi = seed.phi;
        seedRows(row).dx_c = seed.dx_c;
        seedRows(row).weightedRmseMv = seed.weightedRmseMv;
        seedRows(row).linearVpRmseMv = seed.linearVpRmseMv;
        seedRows(row).freqHz = eo * highMap.rotFreqHz;
    end
end
seedRows = seedRows(1:row);
seedTable = struct2table(seedRows);
seedTable = sortrows(seedTable, {'weightedRmseMv', 'linearVpRmseMv', 'EO'}, ...
    {'ascend', 'ascend', 'ascend'});
refineSeedRows = repmat(seedRows(1), numel(eoCandidates), 1);
keepN = 0;
for eo = eoCandidates
    idx = find(seedTable.EO == eo, 1, 'first');
    if isempty(idx)
        continue;
    end
    keepN = keepN + 1;
    refineSeedRows(keepN) = table2struct(seedTable(idx, :));
end
refineSeedRows = refineSeedRows(1:keepN);
refineSeedTable = struct2table(refineSeedRows);

best = struct('objective', inf);
candidateRows = repmat(struct('EO', NaN, 'deltaGapMm', NaN, 'A', NaN, ...
    'phi', NaN, 'dx_c', NaN, 'freqHz', NaN, 'weightedRmseMv', inf, ...
    'plainRmseMv', inf, 'clampFraction', NaN, 'maxOvershootMm', NaN), keepN, 1);
dgLo = max(gMin - gCenter);
dgHi = min(gMax - gCenter);
for i = 1:keepN
    eo = refineSeedTable.EO(i);
    p0 = [refineSeedTable.A(i), refineSeedTable.phi(i), ...
        refineSeedTable.dx_c(i), refineSeedTable.deltaGapMm(i)];
    fun = @(p) low_prior_full_wave_objective_scalar_local( ...
        p, eo, highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        staticState.x0BySensorMm(:).', dgLo, dgHi, fitSettings);
    opts = optimset('Display', 'off', 'MaxIter', 350, 'MaxFunEvals', 1200);
    pOpt = fminsearch(fun, p0, opts);
    [score, detail] = evaluate_low_prior_full_wave_objective_local( ...
        pOpt, eo, highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        staticState.x0BySensorMm(:).', dgLo, dgHi, fitSettings);
    candidateRows(i).EO = eo;
    candidateRows(i).deltaGapMm = detail.deltaGapMm;
    candidateRows(i).A = detail.amplitudeMm;
    candidateRows(i).phi = detail.phi;
    candidateRows(i).dx_c = detail.dx_c;
    candidateRows(i).freqHz = eo * highMap.rotFreqHz;
    candidateRows(i).weightedRmseMv = detail.weightedRmseMv;
    candidateRows(i).plainRmseMv = detail.plainRmseMv;
    candidateRows(i).clampFraction = detail.clampFraction;
    candidateRows(i).maxOvershootMm = detail.maxOvershootMm;
    if score < best.objective
        best = detail;
        best.objective = score;
        best.EO = eo;
        best.freqHz = eo * highMap.rotFreqHz;
    end
end

branch = struct();
branch.method = 'low_speed_prior_gap_library_full_wave_all_eo';
branch.deltaGapMm = best.deltaGapMm;
branch.gBySensorMm = min(max(gCenter(:) + best.deltaGapMm, gMin(:)), gMax(:));
branch.x0BySensorMm = staticState.x0BySensorMm(:) + best.dx_c;
branch.staticX0BySensorMm = staticState.x0BySensorMm(:);
branch.dxCMm = best.dx_c;
branch.freqHz = best.freqHz;
branch.EO = best.EO;
branch.phaseRad = best.phi;
branch.amplitudeMm = best.amplitudeMm;
branch.weightedRmseMv = best.weightedRmseMv;
branch.plainRmseMv = best.plainRmseMv;
branch.VPred = best.vPred;
branch.sensorAffine = best.sensorAffine;
branch.diagInfo = rmfield(best, {'vPred'});
branch.staticState = staticState;
branch.seedTable = seedTable;
branch.refineSeedTable = refineSeedTable;
branch.CandidateTable = struct2table(candidateRows);
end

function seed = solve_low_prior_full_wave_seed_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, ...
    staticX0BySensor, eo, fitSettings)
nSensor = numel(sensorList);
s1 = sin(eo .* highMap.theta_v(:));
c1 = cos(eo .* highMap.theta_v(:));
base = nan(size(highMap.V_a));
dMdx = nan(size(highMap.V_a));
for is = 1:nSensor
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    xq = highMap.x_model_v(mask) - staticX0BySensor(is);
    base(mask) = eval_anchored_template_local(responseSurface, LowTemplate, ...
        sid, gBySensor(is), gCenter(is), xq);
    dMdx(mask) = eval_anchored_template_dfdx_local(responseSurface, LowTemplate, ...
        sid, gBySensor(is), gCenter(is), xq);
end
valid = isfinite(base) & isfinite(dMdx) & isfinite(highMap.V_a) & ...
    isfinite(highMap.W_v) & isfinite(s1) & isfinite(c1);
H = zeros(nnz(valid), nSensor + 3);
y = highMap.V_a(valid) - base(valid);
sensorValid = highMap.S_v(valid);
dValid = dMdx(valid);
rowW = highMap.W_v(valid);
for is = 1:nSensor
    H(:, is) = sensorValid == sensorList(is);
end
H(:, nSensor + 1) = -dValid;
H(:, nSensor + 2) = -dValid .* s1(valid);
H(:, nSensor + 3) = -dValid .* c1(valid);
if size(H, 1) < size(H, 2)
    seed = struct('A', 0, 'phi', 0, 'dx_c', 0, ...
        'weightedRmseMv', inf, 'linearVpRmseMv', inf);
    return;
end
sw = sqrt(max(rowW, eps));
beta = (H .* sw) \ (y .* sw);
dx_c = beta(nSensor + 1);
bs = beta(nSensor + 2);
bc = beta(nSensor + 3);
A = min(hypot(bs, bc), fitSettings.lowSpeedFullWaveAmpLimitMm);
phi = atan2(bc, bs);
dx_c = max(min(dx_c, fitSettings.lowSpeedFullWaveDxLimitMm), -fitSettings.lowSpeedFullWaveDxLimitMm);
linRes = y - H * beta;
dgSeed = mean(gBySensor - gCenter, 'omitnan');
[~, detail] = evaluate_low_prior_full_wave_objective_local( ...
    [A, phi, dx_c, dgSeed], eo, highMap, responseSurface, ...
    LowTemplate, sensorList, gCenter, staticX0BySensor, dgSeed, ...
    dgSeed, fitSettings);
seed = struct('A', A, 'phi', phi, 'dx_c', dx_c, ...
    'weightedRmseMv', detail.weightedRmseMv, ...
    'linearVpRmseMv', sqrt(sum(rowW .* linRes.^2) / max(sum(rowW), eps)));
end

function obj = low_prior_full_wave_objective_scalar_local( ...
    p, eo, highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    staticX0BySensor, dgLo, dgHi, fitSettings)
[obj, ~] = evaluate_low_prior_full_wave_objective_local( ...
    p, eo, highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    staticX0BySensor, dgLo, dgHi, fitSettings);
end

function [score, detail] = evaluate_low_prior_full_wave_objective_local( ...
    pRaw, eo, highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    staticX0BySensor, dgLo, dgHi, fitSettings)
p = pRaw(:).';
ampRaw = p(1);
phiRaw = p(2);
dxRaw = p(3);
dgRaw = p(4);
amp = min(abs(ampRaw), fitSettings.lowSpeedFullWaveAmpLimitMm);
phi = wrap_to_pi_local(phiRaw);
dx_c = max(min(dxRaw, fitSettings.lowSpeedFullWaveDxLimitMm), -fitSettings.lowSpeedFullWaveDxLimitMm);
deltaGap = max(min(dgRaw, dgHi), dgLo);
penalty = 1e6 * ((ampRaw - amp)^2 + (dxRaw - dx_c)^2 + (dgRaw - deltaGap)^2);

gBySensor = gCenter + deltaGap;
u = amp .* sin(eo .* highMap.theta_v(:) + phi);
vPred = nan(size(highMap.V_a));
overshoot = zeros(size(highMap.V_a));
clamped = false(size(highMap.V_a));
sensorAffine = repmat(struct('sensorId', NaN, 'gain', NaN, 'offsetMv', NaN), numel(sensorList), 1);
for is = 1:numel(sensorList)
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    T = get_low_template_sensor_local(LowTemplate, sid);
    xDomain = anchored_template_domain_local(responseSurface, T);
    xRaw = highMap.x_model_v(mask) - staticX0BySensor(is) - dx_c - u(mask);
    xEval = min(max(xRaw, xDomain(1)), xDomain(2));
    clamped(mask) = xRaw < xDomain(1) | xRaw > xDomain(2);
    overshoot(mask) = max(xDomain(1) - xRaw, 0) + max(xRaw - xDomain(2), 0);
    model = eval_anchored_template_local(responseSurface, LowTemplate, ...
        sid, gBySensor(is), gCenter(is), xEval);
    validOffset = isfinite(model) & isfinite(highMap.V_a(mask)) & isfinite(highMap.W_v(mask));
    if nnz(validOffset) < 3
        continue;
    end
    ww = highMap.W_v(mask);
    vv = highMap.V_a(mask);
    Haff = [model(validOffset), ones(nnz(validOffset), 1)];
    sw = sqrt(max(ww(validOffset), eps));
    betaAff = (Haff .* sw) \ (vv(validOffset) .* sw);
    gain = betaAff(1);
    offset = betaAff(2);
    sensorAffine(is).sensorId = sid;
    sensorAffine(is).gain = gain;
    sensorAffine(is).offsetMv = offset;
    predLocal = nan(size(model));
    predLocal(validOffset) = Haff * betaAff;
    vPred(mask) = predLocal;
end
valid = isfinite(vPred) & isfinite(highMap.V_a) & isfinite(highMap.W_v);
res = highMap.V_a(valid) - vPred(valid);
w = highMap.W_v(valid);
if isempty(res)
    mse = inf;
    plainRmse = inf;
else
    mse = sum(w .* res.^2) / max(sum(w), eps);
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
end
overshootPenalty = fitSettings.lowSpeedFullWaveOvershootPenalty * ...
    sum(highMap.W_v(valid) .* overshoot(valid).^2) / max(sum(highMap.W_v(valid)), eps);
score = mse + overshootPenalty + penalty;
detail = struct('amplitudeMm', amp, 'phi', phi, 'dx_c', dx_c, ...
    'deltaGapMm', deltaGap, 'gBySensorMm', gBySensor(:), ...
    'weightedRmseMv', sqrt(mse), 'plainRmseMv', plainRmse, ...
    'score', score, 'overshootPenalty', overshootPenalty, ...
    'clampFraction', nnz(clamped) / max(numel(clamped), 1), ...
    'maxOvershootMm', max(overshoot, [], 'omitnan'), ...
    'sensorAffine', sensorAffine, 'vPred', vPred, ...
    'pointCount', nnz(valid));
end

function xDomain = anchored_template_domain_local(responseSurface, T)
xLow = [min(T.x(:)), max(T.x(:))];
xResp = [min(responseSurface.xGrid(:)), max(responseSurface.xGrid(:))] + T.responseX0Mm;
xDomain = [max(xLow(1), xResp(1)), min(xLow(2), xResp(2))];
if xDomain(1) >= xDomain(2)
    xDomain = xLow;
end
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function apparentGap = refine_gap_range_with_projected_diagnostic_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, apparentGap, ...
    gMin, gMax, freqCandidates)
for is = 1:numel(sensorList)
    sid = sensorList(is);
    freqHz = median(freqCandidates, 'omitnan');
    D = compute_projected_gap_diagnostic_sensor_local( ...
        highMap, responseSurface, LowTemplate, sid, gCenter(is), ...
        apparentGap(is).gAppMm, apparentGap(is).x0AppMm, freqHz);
    if isfinite(D.delta_g_proj_mm)
        halfWidth = min(max(1.5 * D.delta_g_proj_mm, 0.03), 0.25);
        lo = max(gMin(is), apparentGap(is).gAppMm - halfWidth);
        hi = min(gMax(is), apparentGap(is).gAppMm + halfWidth);
        if lo <= hi
            apparentGap(is).gLowerMm = lo;
            apparentGap(is).gUpperMm = hi;
            apparentGap(is).projectedGapScaleMm = halfWidth;
        end
    end
    apparentGap(is).projectedDiagnostic = D;
end
end

function D = compute_projected_gap_diagnostic_sensor_local( ...
    highMap, responseSurface, LowTemplate, sid, gCenter, g, x0, freqHz)
mask = highMap.S_v == sid;
xq = highMap.x_model_v(mask) - x0;
V = highMap.V_a(mask);
w = highMap.W_v(mask);
tRel = highMap.t_v(mask) - min(highMap.t_v(mask));
base = eval_anchored_template_local(responseSurface, LowTemplate, sid, g, gCenter, xq);
[pred, ~] = apply_single_sensor_offset_local(base, V, w);
r = V - pred;
dMdx = eval_anchored_template_dfdx_local(responseSurface, LowTemplate, sid, g, gCenter, xq);
jg = eval_anchored_template_dfdg_local(responseSurface, LowTemplate, sid, g, xq);
Jnu = [-dMdx(:), -dMdx(:) .* sin(2*pi*freqHz*tRel(:)), ...
    -dMdx(:) .* cos(2*pi*freqHz*tRel(:))];
valid = isfinite(r) & isfinite(w) & isfinite(jg) & all(isfinite(Jnu), 2);
rw = sqrt(max(w(valid), eps)) .* r(valid);
jgw = sqrt(max(w(valid), eps)) .* jg(valid);
Jnuw = Jnu(valid, :) .* sqrt(max(w(valid), eps));
[Q, rankNu] = orth_basis_local(Jnuw);
if rankNu > 0
    rPerp = rw - Q * (Q' * rw);
    jgPerp = jgw - Q * (Q' * jgw);
else
    rPerp = rw;
    jgPerp = jgw;
end
Sg = norm(jgPerp);
epsilonPerp = norm(rPerp);
den = dot(jgPerp, jgPerp);
if den > eps
    deltaLs = dot(jgPerp, rPerp) / den;
else
    deltaLs = NaN;
end
D = struct('sensorId', sid, 'gInitMm', g, 'x0InitMm', x0, 'freqHz', freqHz, ...
    'delta_g_ls_mm', deltaLs, 'abs_delta_g_ls_mm', abs(deltaLs), ...
    'delta_g_proj_mm', epsilonPerp / max(Sg, eps), 'epsilon_perp', epsilonPerp, ...
    'Sg', Sg, 'rank_nuisance', rankNu, 'sample_count', nnz(valid));
end

function dMdg = eval_anchored_template_dfdg_local(responseSurface, LowTemplate, sid, g, xq)
T = get_low_template_sensor_local(LowTemplate, sid);
xResp = xq - T.responseX0Mm;
gapStep = max(1e-4, 1e-3 * max(abs(g), 1));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
gapStep = min(gapStep, 0.45 * max(g - gMin, eps));
gapStep = min(gapStep, 0.45 * max(gMax - g, eps));
gapStep = max(gapStep, 1e-6);
dFdg = (eval_response_surface_local(responseSurface, g + gapStep, xResp) - ...
    eval_response_surface_local(responseSurface, g - gapStep, xResp)) ./ (2 * gapStep);
dMdg = T.responseGain .* dFdg;
end

function [Q, rankJ] = orth_basis_local(J)
[Q0, R] = qr(J, 0);
d = abs(diag(R));
tol = max(size(J)) * eps(max(d, [], 'omitnan'));
rankJ = sum(d > tol);
if rankJ == 0
    Q = zeros(size(J, 1), 0);
else
    Q = Q0(:, 1:rankJ);
end
end

function isBoundary = check_gap_grid_boundary_local(gBySensor, gGridBySensor)
isBoundary = false(numel(gBySensor), 1);
for i = 1:numel(gBySensor)
    grid = gGridBySensor{i};
    if isempty(grid) || ~isfinite(gBySensor(i))
        isBoundary(i) = true;
    else
        isBoundary(i) = abs(gBySensor(i) - min(grid)) < 1e-9 || ...
            abs(gBySensor(i) - max(grid)) < 1e-9;
    end
end
end

function diagInfo = fit_peak_coordinate_diagnostic_local(highMap, freqHz, sensorList)
diagInfo = struct('available', false, 'freqHz', freqHz, 'amplitudeMm', NaN, ...
    'cMm', NaN, 'dMm', NaN, 'plainRmseMm', NaN, 'pointCount', 0, ...
    'note', 'x_peak_mm_v was not available.');
if ~isfield(highMap, 'x_peak_mm_v') || isempty(highMap.x_peak_mm_v) || ...
        ~isfinite(freqHz)
    return;
end
keep = isfinite(highMap.S_v) & isfinite(highMap.pulsePeakTime_v) & ...
    isfinite(highMap.x_peak_mm_v);
if nnz(keep) < numel(sensorList) + 4
    diagInfo.note = 'Too few peak-coordinate samples for diagnostic fit.';
    return;
end
T = table(highMap.S_v(keep), highMap.pulsePeakTime_v(keep), ...
    highMap.x_peak_mm_v(keep), 'VariableNames', {'sensorId', 'peakTime', 'xPeakMm'});
T = unique(T, 'rows');
if height(T) < numel(sensorList) + 4
    diagInfo.note = 'Too few unique peak-coordinate samples for diagnostic fit.';
    return;
end
tRel = T.peakTime - min(T.peakTime);
H = zeros(height(T), numel(sensorList) + 2);
for is = 1:numel(sensorList)
    H(:, is) = T.sensorId == sensorList(is);
end
H(:, end-1) = sin(2*pi*freqHz*tRel);
H(:, end) = cos(2*pi*freqHz*tRel);
beta = H \ T.xPeakMm;
pred = H * beta;
res = T.xPeakMm - pred;
diagInfo.available = true;
diagInfo.freqHz = freqHz;
diagInfo.sensorOffsetsMm = beta(1:numel(sensorList));
diagInfo.cMm = beta(end-1);
diagInfo.dMm = beta(end);
diagInfo.amplitudeMm = hypot(beta(end-1), beta(end));
diagInfo.plainRmseMm = sqrt(mean(res.^2, 'omitnan'));
diagInfo.pointCount = height(T);
diagInfo.note = 'Linear diagnostic fit of x_peak_mm_v with per-sensor offsets plus shared sinusoid.';
end

function p = make_param_vector_local(gBySensor, x0BySensor, c, d)
p = [gBySensor(:).', x0BySensor(:).', c, d];
end

function [gBySensor, x0BySensor, c, d] = unpack_param_vector_local(p, nSensor)
gBySensor = p(1:nSensor);
x0BySensor = p(nSensor+1:2*nSensor);
c = p(2*nSensor+1);
d = p(2*nSensor+2);
end

function highMapSmall = decimate_highmap_for_fit_local(highMap, maxFitPoints)
highMapSmall = highMap;
n = numel(highMap.t_v);
if n <= maxFitPoints
    highMapSmall.pointCount = n;
    return;
end
sensorList = unique(highMap.S_v(:).');
keep = false(n, 1);
for sid = sensorList
    idxSensor = find(highMap.S_v == sid);
    nSensorKeep = max(50, round(maxFitPoints * numel(idxSensor) / n));
    stride = max(1, floor(numel(idxSensor) / nSensorKeep));
    keep(idxSensor(1:stride:end)) = true;
end
idx = find(keep);
if numel(idx) > maxFitPoints
    pick = round(linspace(1, numel(idx), maxFitPoints));
    idx = idx(pick);
end
idx = sort(idx);
fields = {'t_v', 'x_v', 'x_timebased_ref_v', 'x_mm_v', 'V_a', 'S_v', ...
    'rev_v', 'W_v', 'theta_v', 'pulsePeakTime_v', 'x_fit_v', 'x_model_v', ...
    'x_peak_mm_v'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(highMapSmall, f) && numel(highMapSmall.(f)) == n
        highMapSmall.(f) = highMapSmall.(f)(idx);
    end
end
highMapSmall.pointCount = numel(idx);
highMapSmall.decimation = struct('enabled', true, ...
    'originalPointCount', n, 'fitPointCount', numel(idx), ...
    'maxFitPoints', maxFitPoints);
end

function bestJoint = fit_projection_coordinate_descent_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    gGridBySensor, x0Grid, freqCandidates, fitSettings, apparentGap)
nSensor = numel(sensorList);
gInitial = nan(1, nSensor);
x0Initial = nan(1, nSensor);
for is = 1:nSensor
    gInitial(is) = apparentGap(is).gAppMm;
    x0Initial(is) = apparentGap(is).x0AppMm;
end
freqScore = inf(numel(freqCandidates), 1);
freqLightCandidate = cell(numel(freqCandidates), 1);
for ifr = 1:numel(freqCandidates)
    candidate = evaluate_projection_state_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        gInitial, x0Initial, freqCandidates(ifr), fitSettings);
    freqScore(ifr) = candidate.objective;
    freqLightCandidate{ifr} = candidate;
end
[~, orderFreq] = sort(freqScore, 'ascend');
keepN = min(fitSettings.fullCoordinateFrequencyKeep, numel(orderFreq));
keepFreqIdx = orderFreq(1:keepN);
freqCandidates = sort(freqCandidates(keepFreqIdx));

bestJoint = struct('objective', inf);
for f = freqCandidates(:).'
    lightFreq = nan(numel(freqLightCandidate), 1);
    for il = 1:numel(freqLightCandidate)
        if ~isempty(freqLightCandidate{il})
            lightFreq(il) = freqLightCandidate{il}.freqHz;
        end
    end
    idxLight = find(abs(lightFreq - f) < 1e-9, 1, 'first');
    if isempty(idxLight)
        gState = gInitial;
        x0State = x0Initial;
    else
        gState = freqLightCandidate{idxLight}.gBySensorMm;
        x0State = freqLightCandidate{idxLight}.x0BySensorMm;
    end
    current = evaluate_projection_state_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
        gState, x0State, f, fitSettings);
    for pass = 1:fitSettings.outerCoordinatePasses
        for is = 1:nSensor
            for gTryValue = gGridBySensor{is}(:).'
                gTry = current.gBySensorMm;
                gTry(is) = gTryValue;
                candidate = evaluate_projection_state_local( ...
                    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
                    gTry, current.x0BySensorMm, f, fitSettings);
                if candidate.objective < current.objective
                    current = candidate;
                end
            end
            for x0TryValue = x0Grid{is}(:).'
                x0Try = current.x0BySensorMm;
                x0Try(is) = x0TryValue;
                candidate = evaluate_projection_state_local( ...
                    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
                    current.gBySensorMm, x0Try, f, fitSettings);
                if candidate.objective < current.objective
                    current = candidate;
                end
            end
        end
    end
    if current.objective < bestJoint.objective
        bestJoint = current;
        bestJoint.diagInfo.fullCoordinateFrequencyCandidatesHz = freqCandidates(:).';
        bestJoint.diagInfo.lightFrequencyScoresMv = freqScore(:).';
    end
end
end

function candidate = evaluate_projection_state_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    gBySensor, x0BySensor, freqHz, fitSettings)
t0 = min(highMap.t_v);
sinTerm = sin(2*pi*freqHz*(highMap.t_v - t0));
cosTerm = cos(2*pi*freqHz*(highMap.t_v - t0));
[obj, pred, affine, diagInfo] = solve_projection_for_fixed_gap_x0_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, ...
    gBySensor, x0BySensor, sinTerm, cosTerm, ...
    fitSettings.coordinateIterations, fitSettings.ampLimitMm, fitSettings.x0ProjectionLimitMm);
diagInfo.gBySensor = gBySensor;
if isfield(diagInfo, 'x0BySensorMm')
    x0BySensor = diagInfo.x0BySensorMm;
else
    diagInfo.x0BySensorMm = x0BySensor;
end
diagInfo.frequencyHz = freqHz;
candidate = struct('objective', obj, ...
    'params', make_param_vector_local(gBySensor, x0BySensor, diagInfo.cMm, diagInfo.dMm), ...
    'freqHz', freqHz, 'pred', pred, 'affine', affine, ...
    'diagInfo', diagInfo, 'gBySensorMm', gBySensor, 'x0BySensorMm', x0BySensor);
end

function [obj, vPred, affine, diagInfo] = solve_projection_for_fixed_gap_x0_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, x0BySensor, ...
    sinTerm, cosTerm, nIter, ampLimitMm, x0ProjectionLimitMm)
c = 0;
d = 0;
u = zeros(size(highMap.t_v));
for iter = 1:max(1, nIter)
    [base, jX0, jSin, jCos, validBase] = build_linearized_projection_terms_local( ...
        highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, ...
        x0BySensor, u, sinTerm, cosTerm);
    [c, d, diagLinear] = solve_cd_variable_projection_local( ...
        highMap.V_a, highMap.W_v, highMap.S_v, base, jX0, jSin, jCos, validBase, ...
        ampLimitMm);
    if isfield(diagLinear, 'x0CorrectionBySensorMm')
        x0BySensor = x0BySensor + diagLinear.x0CorrectionBySensorMm(:).';
        x0BySensor = min(max(x0BySensor, -x0ProjectionLimitMm), x0ProjectionLimitMm);
    end
    u = c .* sinTerm + d .* cosTerm;
end
[baseFinal, ~, ~, ~, validFinal] = build_linearized_projection_terms_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, ...
    x0BySensor, u, sinTerm, cosTerm);
[vPred, affine] = apply_sensor_offset_projection_local(baseFinal, highMap.V_a, highMap.W_v, highMap.S_v);
valid = validFinal & isfinite(vPred) & isfinite(highMap.V_a) & isfinite(highMap.W_v);
res = highMap.V_a - vPred;
if nnz(valid) < 20
    obj = inf;
else
    obj = sqrt(sum(highMap.W_v(valid) .* res(valid).^2) / sum(highMap.W_v(valid)));
end
diagInfo = struct('cMm', c, 'dMm', d, 'amplitudeMm', hypot(c, d), ...
    'plainRmse', sqrt(mean(res(valid).^2, 'omitnan')), ...
    'weightedRmse', obj, 'voltageObjective', obj, 'u', u, ...
    'x0BySensorMm', x0BySensor, ...
    'sinTerm', sinTerm, 'cosTerm', cosTerm, ...
    'variableProjection', diagLinear);
end

function [base, jX0, jSin, jCos, valid] = build_linearized_projection_terms_local( ...
    highMap, responseSurface, LowTemplate, sensorList, gCenter, gBySensor, x0BySensor, ...
    u, sinTerm, cosTerm)
base = nan(size(highMap.V_a));
jX0 = nan(size(highMap.V_a));
jSin = nan(size(highMap.V_a));
jCos = nan(size(highMap.V_a));
for is = 1:numel(sensorList)
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    xq = highMap.x_model_v(mask) - x0BySensor(is) - u(mask);
    base(mask) = eval_anchored_template_local(responseSurface, LowTemplate, ...
        sid, gBySensor(is), gCenter(is), xq);
    dMdx = eval_anchored_template_dfdx_local(responseSurface, LowTemplate, ...
        sid, gBySensor(is), gCenter(is), xq);
    jX0(mask) = -dMdx;
    jSin(mask) = -dMdx .* sinTerm(mask);
    jCos(mask) = -dMdx .* cosTerm(mask);
end
valid = isfinite(base) & isfinite(jX0) & isfinite(jSin) & isfinite(jCos) & ...
    isfinite(highMap.V_a) & isfinite(highMap.W_v);
end

function [c, d, diagInfo] = solve_cd_variable_projection_local( ...
    vObs, w, sensorIds, base, jX0, jSin, jCos, valid, ampLimitMm)
sensorList = unique(sensorIds(:).');
nValid = nnz(valid);
H = zeros(nValid, 2 * numel(sensorList) + 2);
y = zeros(nValid, 1);
ww = zeros(nValid, 1);
row0 = 0;
for i = 1:numel(sensorList)
    sid = sensorList(i);
    mask = valid & sensorIds == sid;
    if nnz(mask) < 3
        continue;
    end
    nRows = nnz(mask);
    rows = row0 + (1:nRows);
    H(rows, i) = 1;
    H(rows, numel(sensorList) + i) = jX0(mask);
    H(rows, end-1) = jSin(mask);
    H(rows, end) = jCos(mask);
    y(rows) = vObs(mask) - base(mask);
    ww(rows) = w(mask);
    row0 = row0 + nRows;
end
H = H(1:row0, :);
y = y(1:row0);
ww = ww(1:row0);
if size(H, 1) < size(H, 2)
    c = 0;
    d = 0;
    diagInfo = struct('rank', rank(H), 'pointCount', size(H, 1), 'note', 'Too few points.');
    return;
end
ws = sqrt(max(ww, eps));
beta = (H .* ws) \ (y .* ws);
c = beta(end-1);
d = beta(end);
amp = hypot(c, d);
ampUpper = ampLimitMm;
if amp > ampUpper
    scale = ampUpper / max(amp, eps);
    c = c * scale;
    d = d * scale;
end
diagInfo = struct('rank', rank(H), 'pointCount', size(H, 1), ...
    'offsetBySensorMv', beta(1:numel(sensorList)), ...
    'x0CorrectionBySensorMm', beta(numel(sensorList)+1:2*numel(sensorList)), ...
    'unclippedCMm', beta(end-1), ...
    'unclippedDMm', beta(end), 'unclippedAmplitudeMm', amp, ...
    'amplitudeUpperMm', ampUpper, 'amplitudeClipped', amp > ampUpper);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function M = eval_anchored_template_local(responseSurface, LowTemplate, sid, g, gCenter, xq)
T = get_low_template_sensor_local(LowTemplate, sid);
lowRaw = interp1(T.x, T.vMv, xq, 'linear', NaN);
xResp = xq - T.responseX0Mm;
gapDynamic = eval_response_surface_local(responseSurface, g, xResp);
gapBase = eval_response_surface_local(responseSurface, gCenter, xResp);
M = lowRaw + T.responseGain .* (gapDynamic - gapBase);
end

function dMdx = eval_anchored_template_dfdx_local(responseSurface, LowTemplate, sid, g, gCenter, xq)
T = get_low_template_sensor_local(LowTemplate, sid);
dLow = finite_difference_interp_local(T.x, T.vMv, xq);
xResp = xq - T.responseX0Mm;
dGapDynamic = eval_response_surface_dfdx_local(responseSurface, g, xResp);
dGapBase = eval_response_surface_dfdx_local(responseSurface, gCenter, xResp);
dMdx = dLow + T.responseGain .* (dGapDynamic - dGapBase);
end

function dFdx = eval_response_surface_dfdx_local(responseSurface, g, xq)
if isfield(responseSurface, 'dFdxGrid') && isfield(responseSurface, 'gTrainMm')
    dGrid = interp1(responseSurface.gTrainMm(:), responseSurface.dFdxGrid.', g, 'linear', 'extrap').';
else
    Fgrid = eval_response_surface_local(responseSurface, g, responseSurface.xGrid);
    dGrid = gradient(Fgrid, responseSurface.xGrid);
end
dFdx = interp1(responseSurface.xGrid, dGrid, xq, 'linear', NaN);
end

function dy = finite_difference_interp_local(x, y, xq)
dyGrid = gradient(y, x);
dy = interp1(x, dyGrid, xq, 'linear', NaN);
end

function T = get_low_template_sensor_local(LowTemplate, sid)
idx = find([LowTemplate.sensor.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('Low-speed template cache does not contain CH%d.', sid);
end
T = LowTemplate.sensor(idx);
end

function [vPred, affine] = apply_sensor_offset_projection_local(modelRaw, vObs, w, sensorIds)
vPred = nan(size(vObs));
sensorList = unique(sensorIds(:).');
affine = repmat(struct('sensorId', NaN, 'gain', 1, 'offset', NaN), numel(sensorList), 1);
for i = 1:numel(sensorList)
    sid = sensorList(i);
    mask = sensorIds == sid & isfinite(modelRaw) & isfinite(vObs) & isfinite(w);
    if nnz(mask) < 3
        continue;
    end
    H = [modelRaw(mask), ones(nnz(mask), 1)];
    sw = sqrt(max(w(mask), eps));
    coeff = (H .* sw) \ (vObs(mask) .* sw);
    gain = coeff(1);
    offset = coeff(2);
    vPred(mask) = H * coeff;
    affine(i).sensorId = sid;
    affine(i).gain = gain;
    affine(i).offset = offset;
end
vPred(~isfinite(vPred)) = 0;
end

function [vPred, offset] = apply_single_sensor_offset_local(modelRaw, vObs, w)
valid = isfinite(modelRaw) & isfinite(vObs) & isfinite(w);
vPred = nan(size(vObs));
if nnz(valid) < 3
    offset = 0;
    vPred(valid) = modelRaw(valid);
    return;
end
H = [modelRaw(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
coeff = (H .* sw) \ (vObs(valid) .* sw);
offset = coeff(2);
vPred(valid) = H * coeff;
end

function out = pack_result_multisensor_local(obj, p, freqHz, pred, affine, diagInfo, highMap, sensorList)
out = struct();
nSensor = numel(sensorList);
[gBySensor, x0BySensor, c, d] = unpack_param_vector_local(p, nSensor);
if isfield(diagInfo, 'gBySensor')
    gBySensor = diagInfo.gBySensor;
end
if isfield(diagInfo, 'x0BySensorMm')
    x0BySensor = diagInfo.x0BySensorMm;
end
out.gBySensorMm = gBySensor(:);
out.x0BySensorMm = x0BySensor(:);
out.cMm = c;
out.dMm = d;
out.freqHz = freqHz;
out.amplitudeMm = hypot(c, d);
if isfield(diagInfo, 'amplitudeMm')
    out.amplitudeMm = diagInfo.amplitudeMm;
end
out.objective = obj;
out.VPred = pred;
out.affine = affine;
out.plainRmseMv = sqrt(mean((highMap.V_a - pred).^2, 'omitnan'));
out.weightedRmseMv = sqrt(sum(highMap.W_v .* (highMap.V_a - pred).^2) / sum(highMap.W_v));
out.diagInfo = diagInfo;
end

function sensorTable = make_sensor_result_table_local(modelResult, sensorList)
gain = nan(numel(sensorList), 1);
offset = nan(numel(sensorList), 1);
for i = 1:numel(sensorList)
    sid = sensorList(i);
    idx = find([modelResult.affine.sensorId] == sid, 1, 'first');
    if ~isempty(idx)
        gain(i) = modelResult.affine(idx).gain;
        offset(i) = modelResult.affine(idx).offset;
    end
end
sensorTable = table(sensorList(:), modelResult.gBySensorMm(:), ...
    modelResult.x0BySensorMm(:), gain, offset, ...
    'VariableNames', {'sensorId', 'gHatMm', 'x0Mm', 'affineGain', 'affineOffsetMv'});
end
