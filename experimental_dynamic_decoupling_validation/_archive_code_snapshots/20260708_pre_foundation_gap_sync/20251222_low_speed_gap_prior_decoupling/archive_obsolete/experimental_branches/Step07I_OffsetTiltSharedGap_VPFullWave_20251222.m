%% Step07I: Low-speed template baseline + offset-tilt shared gap increment
% The first-order VP step only screens EO candidates. Final scoring uses the
% full waveform model:
%
%   V_s(t) = M_s(g0_s + dg_s, x_s(t) - dx_c - u(t))
%
% where M_s is the Step06I low-speed-template anchored incremental gap
% library. When dg_s = 0, M_s reduces to the measured low-speed template.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07i_offset_tilt_shared_gap_vp_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 2, 3];
cfg.freqSearchHz = [100, 1000];
cfg.eoPad = 2;
cfg.vpTopK = parse_positive_integer_env_local('STEP07I_TOP_K_EO', 3);
cfg.maxWindows = parse_positive_integer_env_local('STEP07I_MAX_WINDOWS', inf);
cfg.maxPointsPerWindow = parse_positive_integer_env_local('STEP07I_MAX_POINTS_PER_WINDOW', 900);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP07I_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP07I_DX_LIMIT_MM', 0.20);
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP07I_DELTA_G_LIMIT_MM', 0.25);
cfg.derivativeStepMm = 1e-3;
cfg.weightFloor = 0.05;
cfg.defaultSensorThreshold = 0.5;
cfg.pulseSelectionMode = lower(strtrim(getenv('STEP07I_PULSE_MODE')));
if isempty(cfg.pulseSelectionMode)
    cfg.pulseSelectionMode = 'all';
end
if ~ismember(cfg.pulseSelectionMode, {'single', 'all'})
    error('STEP07I_PULSE_MODE must be "single" or "all".');
end
cfg.dynamicGradientMinRatio = parse_nonnegative_numeric_env_local('STEP07I_DYNAMIC_GRADIENT_MIN_RATIO', 0.08);
cfg.domainSoftMarginMm = parse_nonnegative_numeric_env_local('STEP07I_DOMAIN_SOFT_MARGIN_MM', 0.20);
cfg.overshootPenaltyWeight = parse_nonnegative_numeric_env_local('STEP07I_OVERSHOOT_PENALTY_WEIGHT', 100);
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 450, 'MaxFunEvals', 1800, ...
    'TolX', 1e-5, 'TolFun', 1e-5);

sensorOverride = strtrim(getenv('STEP07I_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07I_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

rotDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
rotTemplateDir = fullfile(rotDir, 'output', 'templates');
rotDynamicDir = fullfile(rotDir, 'output', 'dynamic_maps');
rotResultDir = fullfile(rotDir, 'output', 'identification');

correctedLibFile = fullfile(outDir, sprintf('Step06I_OffsetTiltShared_GapLibrary_20251222_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
templateFile = find_low_speed_template_file_local(rotTemplateDir, cfg.targetBlade, cfg.analysisSensors);
dynamicFile = find_dynamic_map_file_local(rotDynamicDir, cfg.targetBlade, cfg.analysisSensors);
directMainFile = find_direct_main_result_file_local(rotResultDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(correctedLibFile)
    error('Run Step06I first. Missing file: %s', correctedLibFile);
end

St = load(templateFile, 'Template');
Sd = load(dynamicFile, 'DynamicMap');
Sc = load(correctedLibFile, 'CorrectedGapLibrary');
Template = filter_template_sensors_local(St.Template, cfg.analysisSensors, sensorTag);
DynamicMap = filter_dynamic_map_sensors_local(Sd.DynamicMap, cfg.analysisSensors, sensorTag);
CorrectedGapLibrary = Sc.CorrectedGapLibrary;
responseSurface = CorrectedGapLibrary.responseSurface;
directRef = load_direct_reference_local(directMainFile);
directTrend = directRef.Trend;

numWindows = min(numel(DynamicMap.Window), cfg.maxWindows);
fprintf('\n=== Step07I: offset-tilt shared gap VP/full-wave ===\n');
fprintf('Template: %s\n', templateFile);
fprintf('DynamicMap: %s\n', dynamicFile);
fprintf('Corrected gap library: %s\n', correctedLibFile);
fprintf('Windows: %d, sensors: %s, VP Top-K: %d\n', numWindows, mat2str(cfg.analysisSensors), cfg.vpTopK);
fprintf('Model: M_s(g0_s+dg_s, x-dx-u); dg limit %.3f mm. No Step04/EO prior.\n', cfg.deltaGapLimitMm);

WindowResult = struct([]);
trendRows = cell(numWindows, 1);
bestIdx = 1;
bestRmse = inf;

for iw = 1:numWindows
    Wmap = DynamicMap.Window(iw);
    if ~isempty(directRef.WindowResult) && numel(directRef.WindowResult) >= iw && ...
            isfield(directRef.WindowResult(iw), 'bundle')
        bundleFull = build_gap_bundle_from_direct_bundle_local( ...
            directRef.WindowResult(iw).bundle, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
    else
        bundleFull = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg);
    end
    bundle = decimate_bundle_local(bundleFull, cfg.maxPointsPerWindow);
    eoCandidates = build_eo_candidates_local(bundle.rotFreqMeanHz, cfg.freqSearchHz, cfg.eoPad);
    seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg);
    selectedEO = select_topk_eo_local(seedTable, cfg.vpTopK);

    fixedFit = refine_gap_waveform_fit_local(bundle, seedTable, selectedEO, cfg, false);
    freeFit = refine_gap_waveform_fit_local(bundle, seedTable, selectedEO, cfg, true);
    wr = pack_window_result_local(iw, Wmap, bundleFull, seedTable, selectedEO, fixedFit, freeFit);
    if iw == 1
        WindowResult = repmat(wr, numWindows, 1);
    else
        WindowResult(iw) = wr;
    end
    trendRows{iw} = make_trend_row_local(wr, directTrend);
    if freeFit.weightedRmseMv < bestRmse
        bestRmse = freeFit.weightedRmseMv;
        bestIdx = iw;
    end
    fprintf('Window %02d/%02d laps %s: fixed EO%d f %.3f A %.4f RMSE %.2f mV | free EO%d f %.3f A %.4f RMSE %.2f mV, mean dg %.4f\n', ...
        iw, numWindows, mat2str(Wmap.lap_range), fixedFit.EO, fixedFit.freqHz, fixedFit.amplitudeMm, ...
        fixedFit.weightedRmseMv, freeFit.EO, freeFit.freqHz, freeFit.amplitudeMm, ...
        freeFit.weightedRmseMv, mean(freeFit.deltaGapMm, 'omitnan'));
end

Trend = vertcat(trendRows{:});
Summary = build_summary_table_local(Trend, bestIdx);

Result = struct();
Result.dataset = '20251222';
Result.method = 'low_template_anchored_offset_tilt_shared_gap_vp_full_wave';
Result.description = ['Step06I offset-tilt shared gap library is used only as a clearance-change increment. ' ...
    'VP screens EO candidates; final results use full waveform nonlinear refinement.'];
Result.cfg = cfg;
Result.templateFile = templateFile;
Result.dynamicFile = dynamicFile;
Result.correctedLibFile = correctedLibFile;
Result.directMainFile = directMainFile;
Result.WindowResult = WindowResult;
Result.Trend = Trend;
Result.Summary = Summary;
Result.BestWindowIndex = bestIdx;
Result.BestWindow = WindowResult(bestIdx);

matFile = fullfile(outDir, sprintf('Step07I_OffsetTiltSharedGap_VPFullWave_20251222_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
trendCsv = fullfile(outDir, sprintf('Step07I_OffsetTiltSharedGap_VPFullWave_Trend_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
summaryCsv = fullfile(outDir, sprintf('Step07I_OffsetTiltSharedGap_VPFullWave_Summary_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figTrend = fullfile(figDir, sprintf('Step07I_Trend_20251222_B%d_%s.png', cfg.targetBlade, sensorTag));
figBest = fullfile(figDir, sprintf('Step07I_BestWindow_Collapse_20251222_B%d_%s.png', cfg.targetBlade, sensorTag));
figGap = fullfile(figDir, sprintf('Step07I_GapIdentification_20251222_B%d_%s.png', cfg.targetBlade, sensorTag));

save(matFile, 'Result', 'Trend', 'Summary', '-v7.3');
writetable(Trend, trendCsv);
writetable(Summary, summaryCsv);
plot_trend_local(Trend, figTrend);
plot_best_window_collapse_local(Result.BestWindow, cfg, figBest);
plot_gap_identification_local(Result, Trend, cfg, figGap);

fprintf('\nStep07I complete.\n');
disp(Summary);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, trendCsv, figTrend);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', targetBlade, sensorTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
error('No low-speed template file found under %s for %s.', templateDir, sensorTag);
end

function dynamicFile = find_dynamic_map_file_local(dynamicDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
candidate = fullfile(dynamicDir, sprintf('DynamicMap_B%d_%s_SlidingWindows_20251222.mat', targetBlade, sensorTag));
if isfile(candidate)
    dynamicFile = candidate;
    return;
end
error('No DynamicMap file found under %s for %s.', dynamicDir, sensorTag);
end

function resultFile = find_direct_main_result_file_local(resultDir, targetBlade, sensorIds)
sensorTag = ['S', sprintf('%d', sensorIds)];
preferred = {
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_Main_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s_GradientXRange030_20251222.mat', targetBlade, sensorTag)
    };
for i = 1:numel(preferred)
    candidate = fullfile(resultDir, preferred{i});
    if isfile(candidate)
        resultFile = candidate;
        return;
    end
end
files = dir(fullfile(resultDir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*_20251222.mat', ...
    targetBlade, sensorTag)));
if isempty(files)
    resultFile = '';
else
    [~, idx] = max([files.datenum]);
    resultFile = fullfile(files(idx).folder, files(idx).name);
end
end

function directRef = load_direct_reference_local(resultFile)
directRef = struct('Trend', table(), 'WindowResult', struct([]));
if isempty(resultFile) || ~isfile(resultFile)
    return;
end
S = load(resultFile, 'Result');
if isfield(S, 'Result') && isfield(S.Result, 'Trend')
    directRef.Trend = S.Result.Trend;
end
if isfield(S, 'Result') && isfield(S.Result, 'WindowResult')
    directRef.WindowResult = S.Result.WindowResult;
end
end

function bundle = build_gap_bundle_from_direct_bundle_local(Bsrc, Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
if isfield(Bsrc, 'sensor_ids')
    sensorIds = Bsrc.sensor_ids(:).';
else
    sensorIds = cfg.analysisSensors(:).';
end
X = Bsrc.X(:);
T = Bsrc.T(:);
Theta = Bsrc.Theta(:);
W = max(Bsrc.W(:), cfg.weightFloor);
sensorIndex = Bsrc.sensor_index(:);
S = sensorIds(sensorIndex).';
V = nan(size(X));
F0 = nan(size(X));
Fx = nan(size(X));
g0BySensor = NaN(numel(sensorIds), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = sensorIndex == is;
    Tpl = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    V(mask) = (Bsrc.V(mask) - Tpl.baseline) * 1000;
    F0(mask) = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask));
    Fx(mask) = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, X(mask), cfg.derivativeStepMm);
    g0BySensor(is) = corr.g0Mm;
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    sensorInfo(is).pointCount = nnz(mask);
end
valid = isfinite(X) & isfinite(T) & isfinite(Theta) & isfinite(W) & ...
    isfinite(V) & isfinite(F0) & isfinite(Fx) & isfinite(sensorIndex);
bundle = struct();
bundle.X = X(valid);
bundle.T = T(valid);
bundle.TRel = T(valid) - min(T(valid));
bundle.V = V(valid);
bundle.W = W(valid);
bundle.Theta = Theta(valid);
bundle.S = S(valid);
bundle.sensorIndex = sensorIndex(valid);
bundle.sensorIds = sensorIds;
bundle.F0 = F0(valid);
bundle.Fx = Fx(valid);
bundle.Template = Template;
bundle.CorrectedGapLibrary = CorrectedGapLibrary;
bundle.responseSurface = responseSurface;
bundle.g0BySensor = g0BySensor;
bundle.sensorInfo = sensorInfo;
bundle.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
bundle.rotRpmMean = Wmap.rot_rpm_mean;
bundle.windowId = Wmap.window_id;
bundle.lapRange = Wmap.lap_range;
bundle.pointCount = numel(bundle.T);
end

function Template = filter_template_sensors_local(Template, sensorIds, sensorTag)
available = [Template.Sensor.sensor_id];
idx = zeros(size(sensorIds));
for i = 1:numel(sensorIds)
    hit = find(available == sensorIds(i), 1);
    if isempty(hit)
        error('Template does not contain CH%d.', sensorIds(i));
    end
    idx(i) = hit;
end
Template.Sensor = Template.Sensor(idx);
Template.SensorIDs = sensorIds;
Template.SensorTag = sensorTag;
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, sensorIds, sensorTag)
for iw = 1:numel(DynamicMap.Window)
    available = [DynamicMap.Window(iw).Sensor.sensor_id];
    idx = zeros(size(sensorIds));
    for i = 1:numel(sensorIds)
        hit = find(available == sensorIds(i), 1);
        if isempty(hit)
            error('DynamicMap window %d does not contain CH%d.', iw, sensorIds(i));
        end
        idx(i) = hit;
    end
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(idx);
end
DynamicMap.SensorIDs = sensorIds;
DynamicMap.SensorTag = sensorTag;
end

function bundle = build_gap_observation_bundle_local(Wmap, Template, CorrectedGapLibrary, responseSurface, cfg)
X = []; T = []; V = []; W = []; Theta = []; S = []; sensorIndex = [];
F0 = []; Fx = []; g0BySensor = NaN(numel(cfg.analysisSensors), 1);
sensorInfo = repmat(struct('sensorId', NaN, 'xDomain', [NaN, NaN], 'pointCount', 0), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tpl = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    xRaw = D.x_rel(:);
    tRaw = D.t(:);
    thetaRaw = D.theta(:);
    vRawMv = (D.V(:) - Tpl.baseline) * 1000;
    wRaw = max(D.W(:), cfg.weightFloor);
    thresholdMv = (resolve_threshold_local(Tpl, cfg) - Tpl.baseline) * 1000;

    f0Raw = eval_corrected_template_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw);
    fxRaw = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, corr.g0Mm, xRaw, cfg.derivativeStepMm);
    maskPulse = build_pulse_mask_local(vRawMv, thresholdMv, cfg.pulseSelectionMode);
    maskEffective = build_dynamic_effective_mask_local(xRaw, vRawMv, Tpl, thresholdMv, cfg);
    maskDomain = isfinite(f0Raw) & isfinite(fxRaw) & isfinite(vRawMv) & isfinite(thetaRaw);
    mask = maskPulse & maskEffective & maskDomain;
    if nnz(mask) < 8
        mask = maskEffective & maskDomain;
    end
    if nnz(mask) < 8
        mask = maskDomain;
    end
    if nnz(mask) < 8
        continue;
    end

    X = [X; xRaw(mask)]; %#ok<AGROW>
    T = [T; tRaw(mask)]; %#ok<AGROW>
    V = [V; vRawMv(mask)]; %#ok<AGROW>
    W = [W; normalize_weight_local(wRaw(mask), cfg.weightFloor)]; %#ok<AGROW>
    Theta = [Theta; thetaRaw(mask)]; %#ok<AGROW>
    S = [S; repmat(sid, nnz(mask), 1)]; %#ok<AGROW>
    sensorIndex = [sensorIndex; repmat(is, nnz(mask), 1)]; %#ok<AGROW>
    F0 = [F0; f0Raw(mask)]; %#ok<AGROW>
    Fx = [Fx; fxRaw(mask)]; %#ok<AGROW>
    g0BySensor(is) = corr.g0Mm;
    sensorInfo(is).sensorId = sid;
    sensorInfo(is).xDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    sensorInfo(is).pointCount = nnz(mask);
end
if isempty(T)
    error('No valid points were constructed for window %d.', Wmap.window_id);
end
bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.TRel = T - min(T);
bundle.V = V;
bundle.W = W;
bundle.Theta = Theta;
bundle.S = S;
bundle.sensorIndex = sensorIndex;
bundle.sensorIds = cfg.analysisSensors(:).';
bundle.F0 = F0;
bundle.Fx = Fx;
bundle.Template = Template;
bundle.CorrectedGapLibrary = CorrectedGapLibrary;
bundle.responseSurface = responseSurface;
bundle.g0BySensor = g0BySensor;
bundle.sensorInfo = sensorInfo;
bundle.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
bundle.rotRpmMean = Wmap.rot_rpm_mean;
bundle.windowId = Wmap.window_id;
bundle.lapRange = Wmap.lap_range;
bundle.pointCount = numel(T);
end

function bundle = decimate_bundle_local(bundle, maxPoints)
if numel(bundle.T) <= maxPoints
    return;
end
idx = unique(round(linspace(1, numel(bundle.T), maxPoints)));
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    bundle.(f) = bundle.(f)(idx);
end
bundle.pointCount = numel(idx);
end

function seedTable = solve_vp_seed_scan_local(bundle, eoCandidates, cfg)
eoCandidates = unique(round(eoCandidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, ...
    'linearRmseMv', inf, 'fallbackRmseMv', inf), numel(eoCandidates), 1);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo .* bundle.Theta);
    c1 = cos(eo .* bundle.Theta);
    basis = [-bundle.Fx(:), -bundle.Fx(:) .* s1, -bundle.Fx(:) .* c1];
    y = bundle.V(:) - bundle.F0(:);
    sw = sqrt(max(bundle.W(:), cfg.weightFloor));
    coeff = (basis .* sw) \ (y .* sw);
    dx = coeff(1);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    linRes = y - basis * coeff;
    fixedParams = pack_params_local(A, phi, dx, zeros(numel(bundle.sensorIds), 1), false);
    fixedEval = evaluate_gap_waveform_local(fixedParams, eo, bundle, cfg, false);
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx = dx;
    rows(i).linearRmseMv = sqrt(sum(bundle.W(:) .* linRes.^2) / max(numel(linRes), 1));
    rows(i).fallbackRmseMv = fixedEval.weightedRmseMv;
end
[~, order] = sortrows([[rows.fallbackRmseMv].', [rows.linearRmseMv].', [rows.EO].'], [1, 2, 3]);
seedTable = struct2table(rows(order));
end

function eoKeep = select_topk_eo_local(seedTable, topK)
topK = min(max(1, floor(topK)), height(seedTable));
eoKeep = unique(seedTable.EO(1:topK).', 'stable');
end

function fit = refine_gap_waveform_fit_local(bundle, seedTable, eoCandidates, cfg, useGapIncrement)
eoCandidates = unique(round(eoCandidates(:).'));
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'deltaGapMm', [], 'VPred', [], 'uMm', [], 'score', inf), numel(eoCandidates), 1);
best = struct('score', inf);
for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    seedIdx = find(seedTable.EO == eo, 1, 'first');
    if isempty(seedIdx)
        seedIdx = 1;
    end
    seed = seedTable(seedIdx, :);
    theta0 = [seed.A, seed.phi, seed.dx, zeros(1, numel(bundle.sensorIds))];
    if ~useGapIncrement
        theta0 = theta0(1:3);
    end
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfg, useGapIncrement);
    thetaOpt = fminsearch(fun, theta0, cfg.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfg, useGapIncrement, numel(bundle.sensorIds));
    detail = evaluate_gap_waveform_local(thetaOpt, eo, bundle, cfg, useGapIncrement);
    candidate(i).EO = eo;
    candidate(i).freqHz = eo * bundle.rotFreqMeanHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    candidate(i).deltaGapMm = detail.deltaGapMm;
    candidate(i).VPred = detail.VPred;
    candidate(i).uMm = detail.uMm;
    candidate(i).score = detail.score;
    if detail.score < best.score
        best = candidate(i);
    end
end
fit = best;
fit.mode = ternary_local(useGapIncrement, 'free_gap_increment', 'fixed_gap_fallback');
fit.CandidateTable = sortrows(struct2table(rmfield(candidate, {'VPred','uMm','deltaGapMm'})), ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
end

function score = bounded_objective_local(thetaRaw, eo, bundle, cfg, useGapIncrement)
theta = bound_params_local(thetaRaw, cfg, useGapIncrement, numel(bundle.sensorIds));
detail = evaluate_gap_waveform_local(theta, eo, bundle, cfg, useGapIncrement);
score = detail.score + 1e4 * sum((thetaRaw(:) - theta(:)).^2);
end

function theta = bound_params_local(thetaRaw, cfg, useGapIncrement, nSensor)
theta = thetaRaw(:).';
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
theta(3) = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
if useGapIncrement
    if numel(theta) < 3 + nSensor
        theta(4:3+nSensor) = 0;
    end
    theta(4:3+nSensor) = max(min(theta(4:3+nSensor), cfg.deltaGapLimitMm), -cfg.deltaGapLimitMm);
else
    theta = theta(1:3);
end
end

function theta = pack_params_local(A, phi, dx, dg, useGapIncrement)
if useGapIncrement
    theta = [A, phi, dx, dg(:).'];
else
    theta = [A, phi, dx];
end
end

function detail = evaluate_gap_waveform_local(theta, eo, bundle, cfg, useGapIncrement)
theta = bound_params_local(theta, cfg, useGapIncrement, numel(bundle.sensorIds));
A = theta(1);
phi = theta(2);
dx = theta(3);
if useGapIncrement
    dg = theta(4:3+numel(bundle.sensorIds)).';
else
    dg = zeros(numel(bundle.sensorIds), 1);
end
u = A .* sin(eo .* bundle.Theta + phi);
vPred = nan(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensorIds)
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex == is;
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    g = corr.g0Mm + dg(is);
    xEval = bundle.X(mask) - dx - u(mask);
    [vPred(mask), overshoot(mask)] = eval_corrected_template_local(Tpl, ...
        bundle.responseSurface, corr, g, xEval);
end
valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
if nnz(valid) < 8
    weightedRmse = inf;
    plainRmse = inf;
    overshootPenalty = inf;
else
    res = bundle.V(valid) - vPred(valid);
    ww = max(bundle.W(valid), cfg.weightFloor);
    weightedRmse = sqrt(sum(ww .* res.^2) / max(nnz(valid), 1));
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
    overshootPenalty = cfg.overshootPenaltyWeight * sum(ww .* overshoot(valid).^2) / max(nnz(valid), 1);
end
detail = struct();
detail.EO = eo;
detail.freqHz = eo * bundle.rotFreqMeanHz;
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx;
detail.deltaGapMm = dg(:);
detail.weightedRmseMv = weightedRmse;
detail.plainRmseMv = plainRmse;
detail.overshootPenalty = overshootPenalty;
detail.score = weightedRmse + overshootPenalty;
detail.VPred = vPred;
detail.uMm = u;
end

function [model, overshoot] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr)
xLowRaw = xOpr(:);
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xLowRaw, xLo), xHi);
overshootLow = max(xLo - xLowRaw, 0) + max(xLowRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function dMdx = eval_corrected_template_derivative_local(Tpl, responseSurface, corr, g, xOpr, h)
[mp, ~] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr + h);
[mm, ~] = eval_corrected_template_local(Tpl, responseSurface, corr, g, xOpr - h);
dMdx = (mp - mm) ./ (2 * h);
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
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
model = eval_response_surface_local(responseSurface, gEff, xLib);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1);
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function threshold = resolve_threshold_local(Tpl, cfg)
if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = cfg.defaultSensorThreshold;
end
end

function mask = build_pulse_mask_local(vMv, thresholdMv, modeName)
if strcmpi(modeName, 'single')
    mask = isolate_main_pulse_local(vMv, thresholdMv);
else
    mask = isolate_all_pulses_local(vMv, thresholdMv);
end
end

function mask = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
end

function mask = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
bestSeg = 1;
bestPeak = -inf;
for i = 1:numel(starts)
    seg = idx(starts(i):ends(i));
    peakVal = max(v(seg));
    if peakVal > bestPeak
        bestPeak = peakVal;
        bestSeg = i;
    end
end
mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end

function mask = build_dynamic_effective_mask_local(x, vMv, Tpl, thresholdMv, cfg)
finite = isfinite(x) & isfinite(vMv);
gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:) * 1000, x(:), 'linear', 0));
end
gMax = max(gTpl(finite), [], 'omitnan');
if isfinite(gMax) && gMax > 0
    maskGradient = gTpl >= cfg.dynamicGradientMinRatio * gMax;
else
    maskGradient = false(size(x));
end
maskPeak = vMv >= thresholdMv;
mask = finite & (maskGradient | maskPeak);
if nnz(mask) < 8
    mask = finite;
end
end

function w = normalize_weight_local(wRaw, floorW)
w = max(wRaw(:), floorW);
if max(w) > 0
    w = max(floorW, w ./ max(w));
end
end

function eoCandidates = build_eo_candidates_local(rotFreqHz, freqSearchHz, eoPad)
rotFreqHz = max(rotFreqHz, eps);
eoMin = max(1, ceil(min(freqSearchHz) / rotFreqHz));
eoMax = max(eoMin, floor(max(freqSearchHz) / rotFreqHz));
eoCandidates = eoMin:eoMax;
if isempty(eoCandidates)
    eoCenter = max(1, round(mean(freqSearchHz) / rotFreqHz));
    eoCandidates = max(1, eoCenter - eoPad):max(1, eoCenter + eoPad);
end
freq = eoCandidates .* rotFreqHz;
eoCandidates = eoCandidates(freq >= min(freqSearchHz) & freq <= max(freqSearchHz));
end

function wr = pack_window_result_local(iw, Wmap, bundle, seedTable, selectedEO, fixedFit, freeFit)
wr = struct();
wr.windowId = iw;
wr.lapRange = Wmap.lap_range;
wr.timeWindow = Wmap.time_window;
wr.rotFreqMeanHz = Wmap.rot_freq_mean_hz;
wr.rotRpmMean = Wmap.rot_rpm_mean;
wr.bundle = bundle;
wr.VPSeedTable = seedTable;
wr.VPSelectedEO = selectedEO;
wr.fixedGapFallback = fixedFit;
wr.freeGapIncrement = freeFit;
end

function row = make_trend_row_local(wr, directTrend)
directEO = NaN; directFreq = NaN; directAmp = NaN; directRmseMv = NaN;
if ~isempty(directTrend) && height(directTrend) >= wr.windowId
    directEO = directTrend.EO_id(wr.windowId);
    directFreq = directTrend.fn_id(wr.windowId);
    directAmp = directTrend.A_id(wr.windowId);
    directRmseMv = 1000 * directTrend.weighted_voltage_rmse(wr.windowId);
end
dg = wr.freeGapIncrement.deltaGapMm(:);
row = table(wr.windowId, wr.lapRange(1), wr.lapRange(end), wr.rotFreqMeanHz, ...
    directEO, directFreq, directAmp, directRmseMv, ...
    wr.fixedGapFallback.EO, wr.fixedGapFallback.freqHz, wr.fixedGapFallback.amplitudeMm, ...
    wr.fixedGapFallback.dxMm, wr.fixedGapFallback.weightedRmseMv, ...
    wr.freeGapIncrement.EO, wr.freeGapIncrement.freqHz, wr.freeGapIncrement.amplitudeMm, ...
    wr.freeGapIncrement.dxMm, wr.freeGapIncrement.weightedRmseMv, ...
    mean(dg, 'omitnan'), dg(1), dg(min(2, numel(dg))), dg(min(3, numel(dg))), ...
    'VariableNames', {'window_id','lap_start','lap_end','rot_freq_hz', ...
    'direct_EO','direct_frequency_hz','direct_amplitude_mm','direct_rmse_mV', ...
    'fixed_EO','fixed_frequency_hz','fixed_amplitude_mm','fixed_dx_mm','fixed_rmse_mV', ...
    'free_EO','free_frequency_hz','free_amplitude_mm','free_dx_mm','free_rmse_mV', ...
    'free_mean_delta_gap_mm','free_delta_gap_ch1_mm','free_delta_gap_ch3_mm','free_delta_gap_ch6_mm'});
end

function Summary = build_summary_table_local(Trend, bestIdx)
method = ["direct_low_template_main"; "fixed_gap_fallback"; "free_gap_increment"];
dominantEO = [mode_finite_local(Trend.direct_EO); mode_finite_local(Trend.fixed_EO); mode_finite_local(Trend.free_EO)];
meanFreq = [mean(Trend.direct_frequency_hz, 'omitnan'); mean(Trend.fixed_frequency_hz, 'omitnan'); mean(Trend.free_frequency_hz, 'omitnan')];
stdFreq = [std(Trend.direct_frequency_hz, 'omitnan'); std(Trend.fixed_frequency_hz, 'omitnan'); std(Trend.free_frequency_hz, 'omitnan')];
meanAmp = [mean(Trend.direct_amplitude_mm, 'omitnan'); mean(Trend.fixed_amplitude_mm, 'omitnan'); mean(Trend.free_amplitude_mm, 'omitnan')];
meanRmse = [mean(Trend.direct_rmse_mV, 'omitnan'); mean(Trend.fixed_rmse_mV, 'omitnan'); mean(Trend.free_rmse_mV, 'omitnan')];
medianRmse = [median(Trend.direct_rmse_mV, 'omitnan'); median(Trend.fixed_rmse_mV, 'omitnan'); median(Trend.free_rmse_mV, 'omitnan')];
bestWindowIndex = repmat(bestIdx, 3, 1);
Summary = table(method, dominantEO, meanFreq, stdFreq, meanAmp, meanRmse, medianRmse, bestWindowIndex, ...
    'VariableNames', {'method','dominant_EO','mean_frequency_hz','std_frequency_hz', ...
    'mean_amplitude_mm','mean_rmse_mV','median_rmse_mV','best_free_gap_window'});
end

function plot_trend_local(T, figFile)
style = paper_style_local();
fig = figure('Name', 'Step07I trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 13.0]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_EO, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_EO, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.free_EO, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Free gap');
ylabel('EO'); title('Order'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_amplitude_mm, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_amplitude_mm, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.free_amplitude_mm, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Free gap');
ylabel('A (mm)', 'Interpreter', 'tex'); title('Amplitude'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
nexttile; hold on; grid on; box on;
plot(T.window_id, T.direct_rmse_mV, '-', 'Color', style.black, 'LineWidth', 1.0, ...
    'Marker', '.', 'MarkerSize', 9, 'DisplayName', 'Direct');
plot(T.window_id, T.fixed_rmse_mV, '-', 'Color', style.blue, 'LineWidth', 1.0, ...
    'Marker', 'o', 'MarkerSize', 4.2, 'DisplayName', 'Fixed gap');
plot(T.window_id, T.free_rmse_mV, '-', 'Color', style.red, 'LineWidth', 1.0, ...
    'Marker', 's', 'MarkerSize', 4.2, 'DisplayName', 'Free gap');
xlabel('Window'); ylabel('RMSE (mV)'); title('Voltage residual'); legend('Location', 'northeast', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function plot_best_window_collapse_local(wr, cfg, figFile)
style = paper_style_local();
bundle = decimate_bundle_local(wr.bundle, cfg.maxPointsPerWindow);
fit = wr.freeGapIncrement;
fig = figure('Name', 'Step07I best-window dynamic collapse', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 4.8 * numel(cfg.analysisSensors)]);
tiledlayout(numel(cfg.analysisSensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = bundle.sensorIndex == is;
    xRaw = bundle.X(mask) - fit.dxMm;
    xComp = xRaw - fit.uMm(mask);
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    g = corr.g0Mm + fit.deltaGapMm(is);
    [modelNo, ~] = eval_corrected_template_local(Tpl, bundle.responseSurface, corr, g, xRaw);
    [modelDyn, ~] = eval_corrected_template_local(Tpl, bundle.responseSurface, corr, g, xComp);
    [xs1, o1] = sort(xRaw);
    [xs2, o2] = sort(xComp);
    vObs = bundle.V(mask);

    nexttile; hold on; grid on; box on;
    plot(xRaw, vObs, '.', 'Color', style.gray, 'MarkerSize', 4.0, 'DisplayName', 'Samples');
    plot(xs1, modelNo(o1), '-', 'Color', style.red, 'LineWidth', 1.2, 'DisplayName', 'Static template');
    title(sprintf('CH%d before', sid));
    xlabel('x_{OPR} (mm)', 'Interpreter', 'tex'); ylabel('mV'); legend('Location', 'southeast', 'Box', 'off');
    format_axes_local(gca, style);

    nexttile; hold on; grid on; box on;
    plot(xComp, vObs, '.', 'Color', style.gray, 'MarkerSize', 4.0, 'DisplayName', 'Samples');
    plot(xs2, modelDyn(o2), '-', 'Color', style.red, 'LineWidth', 1.2, 'DisplayName', 'Dynamic template');
    title(sprintf('CH%d after, dg = %.4f mm', sid, fit.deltaGapMm(is)));
    xlabel('x_{OPR}-u(t) (mm)', 'Interpreter', 'tex'); ylabel('mV'); legend('Location', 'southeast', 'Box', 'off');
    format_axes_local(gca, style);
end
export_paper_figure_local(fig, figFile);
end

function plot_gap_identification_local(Result, T, cfg, figFile)
style = paper_style_local();
fig = figure('Name', 'Step07I gap identification', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17.0, 10.5]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sensorIds = cfg.analysisSensors(:).';
g0 = Result.BestWindow.bundle.g0BySensor(:).';
dgMat = collect_delta_gap_matrix_local(Result.WindowResult, sensorIds);
absG = dgMat + g0;
colors = sensor_colors_local(numel(sensorIds), style);

nexttile; hold on; box on;
for is = 1:numel(sensorIds)
    plot(T.window_id, dgMat(:, is), '-', 'Color', colors(is, :), ...
        'LineWidth', 1.1, 'Marker', sensor_marker_local(is), ...
        'MarkerSize', 4.2, 'DisplayName', sprintf('CH%d', sensorIds(is)));
end
yline(0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', 'dg = 0');
ylabel('\Delta g (mm)', 'Interpreter', 'tex');
title('Gap increment');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);

nexttile; hold on; box on;
for is = 1:numel(sensorIds)
    plot(T.window_id, absG(:, is), '-', 'Color', colors(is, :), ...
        'LineWidth', 1.1, 'Marker', sensor_marker_local(is), ...
        'MarkerSize', 4.2, 'DisplayName', sprintf('CH%d', sensorIds(is)));
    yline(g0(is), ':', 'Color', colors(is, :), 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end
xlabel('Window');
ylabel('g_0+\Delta g (mm)', 'Interpreter', 'tex');
title('Effective gap');
legend('Location', 'eastoutside', 'Box', 'off');
format_axes_local(gca, style);
export_paper_figure_local(fig, figFile);
end

function dgMat = collect_delta_gap_matrix_local(WindowResult, sensorIds)
nWin = numel(WindowResult);
nSensor = numel(sensorIds);
dgMat = NaN(nWin, nSensor);
for iw = 1:nWin
    if ~isfield(WindowResult(iw), 'freeGapIncrement') || ...
            ~isfield(WindowResult(iw).freeGapIncrement, 'deltaGapMm')
        continue;
    end
    dg = WindowResult(iw).freeGapIncrement.deltaGapMm(:).';
    nUse = min(nSensor, numel(dg));
    dgMat(iw, 1:nUse) = dg(1:nUse);
end
end

function style = paper_style_local()
style.fontName = 'Times New Roman';
style.fontSize = 9;
style.tickFontSize = 8.2;
style.black = [0.08, 0.08, 0.08];
style.gray = [0.35, 0.35, 0.35];
style.blue = [0.00, 0.28, 0.70];
style.red = [0.82, 0.10, 0.10];
style.green = [0.00, 0.45, 0.28];
style.orange = [0.85, 0.37, 0.05];
end

function colors = sensor_colors_local(n, style)
base = [
    style.blue
    style.red
    style.green
    style.orange
    0.45, 0.20, 0.70
    0.15, 0.55, 0.70
    ];
if n <= size(base, 1)
    colors = base(1:n, :);
else
    colors = lines(n);
end
end

function marker = sensor_marker_local(idx)
markers = {'o', 's', '^', 'd', 'v', '>'};
marker = markers{mod(idx - 1, numel(markers)) + 1};
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on', ...
    'XMinorTick', 'off', 'YMinorTick', 'off');
grid(ax, 'off');
xlabel(ax, get(get(ax, 'XLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
ylabel(ax, get(get(ax, 'YLabel'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize);
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
lgd = legend(ax);
if ~isempty(lgd) && isvalid(lgd)
    set(lgd, 'FontName', style.fontName, 'FontSize', style.tickFontSize);
end
end

function export_paper_figure_local(fig, pngFile)
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
pdfFile = fullfile(folder, [name, '.pdf']);
emfFile = fullfile(folder, [name, '.emf']);
try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
catch
    warning('Could not export PDF figure: %s', pdfFile);
end
try
    print(fig, emfFile, '-dmeta', '-r300');
catch
    warning('Could not export EMF figure: %s', emfFile);
end
end

function v = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    v = defaultValue;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = defaultValue;
else
    v = max(1, floor(tmp));
end
end

function v = parse_nonnegative_numeric_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    v = defaultValue;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp < 0
    error('%s must be a nonnegative numeric value.', name);
end
v = tmp;
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function out = ternary_local(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function m = mode_finite_local(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mode(x);
end
end

