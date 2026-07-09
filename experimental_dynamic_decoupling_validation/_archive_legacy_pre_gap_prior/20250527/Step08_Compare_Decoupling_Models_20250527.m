%% Step08: Compare clearance-vibration model variants
% This step fits one static effective clearance for each sensor and one
% shared vibration displacement for the selected blade:
%
%   V_s(x,t) = a_s F(g_s, x - x0_s - u(t)) + c_s
%   u(t)     = C sin(2*pi*f*t) + D cos(2*pi*f*t)
%
% F(g,x) is the static response surface from Step05.  Step06 already
% subtracts sensor-specific continuous-waveform dynamic baselines.

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
cfg.frequencyHalfWidthHz = 35;
cfg.frequencyFineOffsetsHz = -0.18:0.02:0.18;

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(highMapFile)
    error('Run Step06 first. Missing file: %s', highMapFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
H = load(highMapFile, 'highMap');
highMapFull = H.highMap;
[highMapFull, coordinateSource] = attach_model_coordinate_local(highMapFull);
windowSpecs = build_sliding_window_specs_local(highMapFull);

fitSettings = struct();
fitSettings.maxFitPoints = 3500;
fitSettings.x0LimitMm = 0.60;
fitSettings.ampLimitMm = 0.80;
fitSettings.options = optimset('Display', 'off', 'MaxIter', 160, 'MaxFunEvals', 900, ...
    'TolX', 1e-5, 'TolFun', 1e-5);

fprintf('\n=== Step08: sliding-window model comparison ===\n');
fprintf('HighMap source: %s\n', highMapFile);
fprintf('Analysis sensors: %s, target blade B%d\n', mat2str(cfg.analysisSensors), cfg.targetBlade);
fprintf('Sliding windows: %d windows, %d laps/window, step %d lap(s)\n', ...
    numel(windowSpecs), highMapFull.analysisWinSize, highMapFull.slidingStep);
fprintf('Model coordinate source: %s\n', coordinateSource);

WindowResult = [];
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
        highMapWinFull, responseSurface, cfg, fitSettings, coordinateSource);
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
result.referenceOrder = highMapFull.referenceOrder;
result.coordinateUnit = 'mm circumferential displacement';
result.gMinMm = min(responseSurface.gTrainMm);
result.gMaxMm = max(responseSurface.gTrainMm);
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
result.staticTemplate = result.BestWindow.staticTemplate;
result.fixedGapVibration = result.BestWindow.fixedGapVibration;
result.multiSensorJoint = result.BestWindow.multiSensorJoint;
result.peakCoordinateDiagnostic = result.BestWindow.peakCoordinateDiagnostic;
result.note = ['Multi-sensor effective-clearance validation. Each sensor has an independent static effective gap ' ...
    'and coordinate offset; the selected blade shares one sinusoidal vibration displacement. ' ...
    'Step05 and Step06 both subtract continuous no-pulse baselines before matching. ' ...
    'Step07 now evaluates the same 1.5 s / 20-lap / 3-lap sliding-window route as the low-speed rotating-template validation.'];

sensorList = result.analysisSensors;
summaryTable = table( ...
    ["fixed_static"; "fixed_gap_vibration"; "multi_sensor_gap_vibration"], ...
    [result.staticTemplate.plainRmseMv; result.fixedGapVibration.plainRmseMv; result.multiSensorJoint.plainRmseMv], ...
    [result.staticTemplate.weightedRmseMv; result.fixedGapVibration.weightedRmseMv; result.multiSensorJoint.weightedRmseMv], ...
    [result.staticTemplate.freqHz; result.fixedGapVibration.freqHz; result.multiSensorJoint.freqHz], ...
    [result.staticTemplate.amplitudeMm; result.fixedGapVibration.amplitudeMm; result.multiSensorJoint.amplitudeMm], ...
    'VariableNames', {'method', 'plainRmseMv', 'weightedRmseMv', 'freqHz', 'amplitudeMm'});

sensorTable = make_sensor_result_table_local(result.multiSensorJoint, sensorList);

matFile = fullfile(outDir, sprintf('Step08_Model_Comparison_20250527_B%d_%s.mat', ...
    result.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step08_Model_Comparison_Summary_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag));
sensorCsvFile = fullfile(outDir, sprintf('Step08_Model_Comparison_Gap_Summary_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag));
save(matFile, 'result', 'summaryTable', 'sensorTable', '-v7.3');
writetable(summaryTable, csvFile);
writetable(sensorTable, sensorCsvFile);
writetable(trendTable, fullfile(outDir, sprintf('Step08_Model_Comparison_Trend_20250527_B%d_%s.csv', ...
    result.targetBlade, sensorTag)));

fprintf('\nStep08 complete.\n');
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
    if isfield(highMapWin, 'referenceOrder') && isfinite(highMapWin.referenceOrder)
        highMapWin.referenceFreqHz = highMapWin.referenceOrder * localRotFreqHz;
    end
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

function [windowResult, trendRow] = fit_one_highmap_window_local(highMapWinFull, responseSurface, cfg, fitSettings, coordinateSource)
highMapFullPointCount = highMapWinFull.pointCount;
highMap = decimate_highmap_for_fit_local(highMapWinFull, fitSettings.maxFitPoints);
sensorList = highMap.analysisSensors(:).';
nSensor = numel(sensorList);

gMin = min(responseSurface.gTrainMm);
gMax = max(responseSurface.gTrainMm);
gFixed = median(responseSurface.gTrainMm);
x0LimitMm = fitSettings.x0LimitMm;
ampLimitMm = fitSettings.ampLimitMm;
freqCandidates = build_frequency_candidates_local(highMap, cfg);
vibrationStarts = make_vibration_starts_local();
opts = fitSettings.options;

staticFun = @(p) objective_multisensor_local( ...
    make_param_vector_local(gFixed * ones(1, nSensor), p(1:nSensor), 0, 0), ...
    highMap, responseSurface, NaN, sensorList, gMin, gMax, x0LimitMm, ampLimitMm);
pStatic0 = zeros(1, nSensor);
pStatic = fminsearch(staticFun, pStatic0, opts);
[objStatic, predStatic, affineStatic, diagStatic] = staticFun(pStatic);

bestFixedGap = struct('objective', inf);
for f = freqCandidates
    vibFun = @(p) objective_multisensor_local( ...
        make_param_vector_local(gFixed * ones(1, nSensor), p(1:nSensor), p(nSensor+1), p(nSensor+2)), ...
        highMap, responseSurface, f, sensorList, gMin, gMax, x0LimitMm, ampLimitMm);
    for istart = 1:size(vibrationStarts, 1)
        p0 = [zeros(1, nSensor), vibrationStarts(istart, :)];
        pOpt = fminsearch(vibFun, p0, opts);
        [obj, pred, affine, diagInfo] = vibFun(pOpt);
        if obj < bestFixedGap.objective
            bestFixedGap = struct('objective', obj, 'params', pOpt, 'freqHz', f, ...
                'pred', pred, 'affine', affine, 'diagInfo', diagInfo);
        end
    end
end

bestJoint = struct('objective', inf);
for f = freqCandidates
    jointFun = @(p) objective_multisensor_local(p, highMap, responseSurface, f, ...
        sensorList, gMin, gMax, x0LimitMm, ampLimitMm);
    for istart = 1:size(vibrationStarts, 1)
        p0 = make_param_vector_local(gFixed * ones(1, nSensor), zeros(1, nSensor), ...
            vibrationStarts(istart, 1), vibrationStarts(istart, 2));
        pOpt = fminsearch(jointFun, p0, opts);
        [obj, pred, affine, diagInfo] = jointFun(pOpt);
        if obj < bestJoint.objective
            bestJoint = struct('objective', obj, 'params', pOpt, 'freqHz', f, ...
                'pred', pred, 'affine', affine, 'diagInfo', diagInfo);
        end
    end
end

windowResult = struct();
windowResult.window_id = highMapWinFull.window_id;
windowResult.lapRange = highMapWinFull.lapRange;
windowResult.timeWindow = highMapWinFull.time_window;
windowResult.coordinateSource = coordinateSource;
windowResult.frequencyCandidatesHz = freqCandidates(:);
windowResult.fullPointCount = highMapFullPointCount;
windowResult.fitPointCount = highMap.pointCount;
windowResult.staticTemplate = pack_result_multisensor_local(objStatic, ...
    make_param_vector_local(gFixed * ones(1, nSensor), pStatic, 0, 0), NaN, ...
    predStatic, affineStatic, diagStatic, highMap, sensorList);
windowResult.fixedGapVibration = pack_result_multisensor_local(bestFixedGap.objective, ...
    make_param_vector_local(gFixed * ones(1, nSensor), bestFixedGap.params(1:nSensor), ...
    bestFixedGap.params(nSensor+1), bestFixedGap.params(nSensor+2)), ...
    bestFixedGap.freqHz, bestFixedGap.pred, bestFixedGap.affine, bestFixedGap.diagInfo, highMap, sensorList);
windowResult.multiSensorJoint = pack_result_multisensor_local(bestJoint.objective, ...
    bestJoint.params, bestJoint.freqHz, bestJoint.pred, bestJoint.affine, ...
    bestJoint.diagInfo, highMap, sensorList);
windowResult.summaryTable = table( ...
    ["fixed_static"; "fixed_gap_vibration"; "multi_sensor_gap_vibration"], ...
    [windowResult.staticTemplate.plainRmseMv; windowResult.fixedGapVibration.plainRmseMv; windowResult.multiSensorJoint.plainRmseMv], ...
    [windowResult.staticTemplate.weightedRmseMv; windowResult.fixedGapVibration.weightedRmseMv; windowResult.multiSensorJoint.weightedRmseMv], ...
    [windowResult.staticTemplate.freqHz; windowResult.fixedGapVibration.freqHz; windowResult.multiSensorJoint.freqHz], ...
    [windowResult.staticTemplate.amplitudeMm; windowResult.fixedGapVibration.amplitudeMm; windowResult.multiSensorJoint.amplitudeMm], ...
    'VariableNames', {'method', 'plainRmseMv', 'weightedRmseMv', 'freqHz', 'amplitudeMm'});
windowResult.sensorTable = make_sensor_result_table_local(windowResult.multiSensorJoint, sensorList);
windowResult.peakCoordinateDiagnostic = fit_peak_coordinate_diagnostic_local( ...
    highMapWinFull, windowResult.multiSensorJoint.freqHz, sensorList);

trendRow = table(windowResult.window_id, windowResult.lapRange(1), windowResult.lapRange(end), ...
    mean(windowResult.timeWindow), highMapWinFull.rotFreqHz, highMapWinFull.rotRpm, ...
    windowResult.multiSensorJoint.freqHz, windowResult.multiSensorJoint.freqHz / highMapWinFull.rotFreqHz, ...
    windowResult.multiSensorJoint.amplitudeMm, ...
    mean(windowResult.multiSensorJoint.gBySensorMm, 'omitnan'), ...
    windowResult.multiSensorJoint.weightedRmseMv, ...
    windowResult.peakCoordinateDiagnostic.amplitudeMm, highMap.pointCount, ...
    'VariableNames', {'window_id', 'lap_start', 'lap_end', 'window_center_time', ...
    'rot_freq_hz', 'rot_rpm', 'frequency_hz', 'EO', 'amplitude_mm', ...
    'mean_gap_mm', 'weighted_rmse_mV', 'peak_coordinate_amplitude_mm', 'point_count'});
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

function freqCandidates = build_frequency_candidates_local(highMap, cfg)
if isfield(highMap, 'rotFreqHz') && isfinite(highMap.rotFreqHz) && highMap.rotFreqHz > 0 && ...
        isfield(highMap, 'referenceOrder') && isfinite(highMap.referenceOrder)
    freqCenter = highMap.referenceOrder * highMap.rotFreqHz;
else
    freqCenter = highMap.referenceFreqHz;
end
if ~isfinite(freqCenter) || freqCenter <= 0
    error('Cannot build Step08 frequency candidates without a valid EO reference or reference frequency.');
end
if isfield(cfg, 'frequencyFineOffsetsHz')
    offsets = cfg.frequencyFineOffsetsHz(:);
else
    offsets = cfg.eoScanOffsetsHz(:);
end
freqCandidates = freqCenter + offsets;
freqCandidates = unique(freqCandidates(isfinite(freqCandidates) & freqCandidates > 0));
freqCandidates = sort(freqCandidates(:).');
end

function starts = make_vibration_starts_local()
ampStarts = [0.20, 0.38];
phaseStarts = 0;
starts = zeros(numel(ampStarts) * numel(phaseStarts), 2);
row = 0;
for ia = 1:numel(ampStarts)
    for ip = 1:numel(phaseStarts)
        row = row + 1;
        starts(row, :) = [ampStarts(ia) * cos(phaseStarts(ip)), ...
            ampStarts(ia) * sin(phaseStarts(ip))];
    end
end
starts = unique(starts, 'rows', 'stable');
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

function [obj, vPredCal, affine, diagInfo] = objective_multisensor_local(p, highMap, responseSurface, freqHz, sensorList, gMin, gMax, x0LimitMm, ampLimitMm)
nSensor = numel(sensorList);
[gBySensor, x0BySensor, c, d] = unpack_param_vector_local(p, nSensor);
penalty = 0;
gBySensorRaw = gBySensor;
gBySensor = min(max(gBySensor, gMin), gMax);
penalty = penalty + 1e12 * sum((gBySensorRaw - gBySensor).^2);
x0Raw = x0BySensor;
x0BySensor = min(max(x0BySensor, -x0LimitMm), x0LimitMm);
penalty = penalty + 1e12 * sum((x0Raw - x0BySensor).^2);
amp = hypot(c, d);
if amp > ampLimitMm
    penalty = penalty + 1e12 * (amp - ampLimitMm)^2;
    scale = ampLimitMm / max(amp, eps);
    c = c * scale;
    d = d * scale;
end

if isnan(freqHz)
    sinTerm = zeros(size(highMap.t_v));
    cosTerm = zeros(size(highMap.t_v));
else
    t0 = min(highMap.t_v);
    sinTerm = sin(2*pi*freqHz*(highMap.t_v - t0));
    cosTerm = cos(2*pi*freqHz*(highMap.t_v - t0));
end
u = c .* sinTerm + d .* cosTerm;

modelRaw = nan(size(highMap.V_a));
for is = 1:nSensor
    sid = sensorList(is);
    mask = highMap.S_v == sid;
    xq = highMap.x_model_v(mask) - x0BySensor(is) - u(mask);
    modelRaw(mask) = eval_response_surface_local(responseSurface, gBySensor(is), xq);
end
modelRaw(~isfinite(modelRaw)) = 0;
[vPredCal, affine] = apply_sensor_affine_projection_local(modelRaw, highMap.V_a, highMap.W_v, highMap.S_v);
res = highMap.V_a - vPredCal;
obj = sum(highMap.W_v .* res.^2) + penalty;
diagInfo = struct('gBySensor', gBySensor, 'x0BySensorMm', x0BySensor, ...
    'cMm', c, 'dMm', d, 'amplitudeMm', hypot(c, d), ...
    'plainRmse', sqrt(mean(res.^2, 'omitnan')), ...
    'weightedRmse', sqrt(sum(highMap.W_v .* res.^2) / sum(highMap.W_v)), ...
    'u', u, 'sinTerm', sinTerm, 'cosTerm', cosTerm);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function [vPred, affine] = apply_sensor_affine_projection_local(modelRaw, vObs, w, sensorIds)
vPred = nan(size(vObs));
sensorList = unique(sensorIds(:).');
affine = repmat(struct('sensorId', NaN, 'gain', NaN, 'offset', NaN), numel(sensorList), 1);
for i = 1:numel(sensorList)
    sid = sensorList(i);
    mask = sensorIds == sid & isfinite(modelRaw) & isfinite(vObs) & isfinite(w);
    if nnz(mask) < 3
        continue;
    end
    A = [modelRaw(mask), ones(nnz(mask), 1)];
    ws = sqrt(max(w(mask), eps));
    coeff = (A .* ws) \ (vObs(mask) .* ws);
    vPred(mask) = A * coeff;
    affine(i).sensorId = sid;
    affine(i).gain = coeff(1);
    affine(i).offset = coeff(2);
end
vPred(~isfinite(vPred)) = 0;
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
