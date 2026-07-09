%% Step07A: Gap-library static-template full-wave identification
% This diagnostic/main candidate avoids VP linearization. For each sliding
% window it first builds a high-speed no-vibration template directly from
% the static gap response surface, then scans integer EO values with the
% complete nonlinear waveform:
%
%   V_s(t) = a_s F(g_s, x_s(t) - x0_s - dx_s - dx_c - A sin(EO*theta(t)+phi)) + b_s
%
% The low-speed template enters through Step06B gap-prior bounds and the
% OPR-to-library x-axis registration x0_s. It is not used as the high-speed
% waveform body. By default, EO is selected independently in each window
% from this window's BTT waveform objective.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.eoCandidates = 3:24;
cfg.maxFitPointsPerWindow = 420;
cfg.globalStaticMaxFitPoints = 1800;
cfg.maxFitPointsPerSensorStatic = 220;
cfg.staticGapGridCount = 41;
cfg.staticDxGridMm = -0.30:0.01:0.30;
cfg.staticDxCorrectionGridMm = -0.08:0.005:0.08;
cfg.useLowSpeedXRegistration = false;
cfg.dynamicAmpLimitMm = 0.80;
cfg.dynamicDxLimitMm = 0.20;
cfg.dynamicAmpSeedsMm = [0, 0.20, 0.40, 0.60];
cfg.dynamicPhiSeedsRad = linspace(0, 2*pi, 9);
cfg.dynamicPhiSeedsRad(end) = [];
cfg.overshootPenaltyMvPerMm = 800;
cfg.minValidPoints = 60;
cfg.useStep04OrderGate = false;

sensorOverride = strtrim(getenv('STEP07A_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07A_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
maxWindowOverride = str2double(strtrim(getenv('STEP07A_MAX_WINDOWS')));

sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];
responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
gapPriorFile = fullfile(outDir, sprintf('Step06B_LowSpeed_GapPrior_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
step04RegionFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20250527.csv');

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
highMapFull = attach_model_coordinate_local(H.highMap);
P = load(gapPriorFile, 'GapPrior');
GapPrior = P.GapPrior;
if cfg.useStep04OrderGate && isfile(step04RegionFile)
    step04RegionTable = readtable(step04RegionFile);
else
    step04RegionTable = table();
end

windowSpecs = build_sliding_window_specs_local(highMapFull);
if isfinite(maxWindowOverride) && maxWindowOverride > 0
    windowSpecs = windowSpecs(1:min(numel(windowSpecs), floor(maxWindowOverride)));
end

fprintf('\n=== Step07A: gap-library static-template full-wave identification ===\n');
fprintf('HighMap: %s\n', highMapFile);
fprintf('Response surface: %s\n', responseFile);
fprintf('Gap prior: %s\n', gapPriorFile);
fprintf('Windows: %d, sensors: %s\n', numel(windowSpecs), mat2str(cfg.analysisSensors));

staticMap = decimate_highmap_for_fit_local(highMapFull, cfg.globalStaticMaxFitPoints);
globalStaticState = fit_static_gap_templates_local(staticMap, responseSurface, GapPrior, cfg);
fprintf('Global static template RMSE %.3f mV; g = %s mm; dx = %s mm\n', ...
    mean([globalStaticState.weightedRmseMv], 'omitnan'), ...
    mat2str(round([globalStaticState.gMm], 4)), ...
    mat2str(round([globalStaticState.dxMm], 4)));

WindowResult = struct([]);
trendRows = cell(numel(windowSpecs), 1);
bestIdx = 1;
bestScore = inf;
for iw = 1:numel(windowSpecs)
    Wfull = subset_highmap_by_laps_local(highMapFull, windowSpecs(iw).lapRange);
    W = decimate_highmap_for_fit_local(Wfull, cfg.maxFitPointsPerWindow);
    allowedEO = get_step04_allowed_eo_local(step04RegionTable, mean(Wfull.t_v, 'omitnan'), cfg.eoCandidates);
    fprintf('\nWindow %02d/%02d, laps %s, points %d\n', ...
        iw, numel(windowSpecs), mat2str(windowSpecs(iw).lapRange), numel(W.t_v));

    staticState = globalStaticState;
    dyn = fit_dynamic_full_wave_local(W, responseSurface, staticState, cfg, allowedEO);

    wr = pack_window_result_local(iw, windowSpecs(iw).lapRange, Wfull, staticState, dyn);
    if iw == 1
        WindowResult = repmat(wr, numel(windowSpecs), 1);
    else
        WindowResult(iw) = wr;
    end
    trendRows{iw} = make_trend_row_local(wr);
    if dyn.weightedRmseMv < bestScore
        bestScore = dyn.weightedRmseMv;
        bestIdx = iw;
    end
    fprintf('  static RMSE %.3f mV | EO%d, f %.3f Hz, A %.4f mm, RMSE %.3f mV\n', ...
        mean([staticState.weightedRmseMv], 'omitnan'), dyn.EO, dyn.freqHz, ...
        dyn.amplitudeMm, dyn.weightedRmseMv);
end

trendTable = vertcat(trendRows{:});
best = WindowResult(bestIdx).dynamic;
summaryTable = table(bestIdx, best.freqHz, best.EO, best.amplitudeMm, ...
    mean(best.gBySensorMm, 'omitnan'), best.weightedRmseMv, ...
    'VariableNames', {'bestWindowIndex', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv'});
sensorTable = table(cfg.analysisSensors(:), best.gBySensorMm(:), best.dxBySensorMm(:), ...
    [best.sensorAffine.gain].', [best.sensorAffine.offsetMv].', ...
    'VariableNames', {'sensorId', 'gHatMm', 'dxMm', 'affineGain', 'affineOffsetMv'});

result = struct();
result.dataset = '20250527';
result.method = 'gap_library_static_template_full_wave';
result.description = ['Per-window static high-speed templates are built from the gap response surface; ' ...
    'dynamic vibration is fitted with the complete nonlinear waveform without VP linearization.'];
result.cfg = cfg;
result.responseFile = responseFile;
result.highMapFile = highMapFile;
result.gapPriorFile = gapPriorFile;
result.WindowResult = WindowResult;
result.Trend = trendTable;
result.BestWindowIndex = bestIdx;
result.BestWindow = WindowResult(bestIdx);
result.dynamic = best;

matFile = fullfile(outDir, sprintf('Step07A_GapLibrary_StaticTemplate_FullWave_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step07A_GapLibrary_StaticTemplate_FullWave_Summary_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
trendCsv = fullfile(outDir, sprintf('Step07A_GapLibrary_StaticTemplate_FullWave_Trend_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
sensorCsv = fullfile(outDir, sprintf('Step07A_GapLibrary_StaticTemplate_FullWave_Sensor_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'result', 'summaryTable', 'sensorTable', '-v7.3');
writetable(summaryTable, csvFile);
writetable(trendTable, trendCsv);
writetable(sensorTable, sensorCsv);

fprintf('\nStep07A complete.\n');
disp(summaryTable);
disp(sensorTable);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, csvFile, trendCsv);

%% Local functions
function highMap = attach_model_coordinate_local(highMap)
if isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(highMap.t_v) && any(isfinite(highMap.x_fit_v))
    highMap.x_model_v = highMap.x_fit_v(:);
else
    highMap.x_model_v = highMap.x_v(:);
end
end

function windowSpecs = build_sliding_window_specs_local(highMap)
laps = unique(highMap.rev_v(:).');
laps = sort(laps(isfinite(laps)));
winSize = highMap.analysisWinSize;
step = highMap.slidingStep;
nWin = floor((numel(laps) - winSize) / step) + 1;
windowSpecs = repmat(struct('window_id', NaN, 'lapRange', []), nWin, 1);
for iw = 1:nWin
    i0 = 1 + (iw - 1) * step;
    windowSpecs(iw).window_id = iw;
    windowSpecs(iw).lapRange = laps(i0:(i0 + winSize - 1));
end
end

function W = subset_highmap_by_laps_local(highMap, lapRange)
mask = ismember(highMap.rev_v, lapRange);
W = highMap;
fields = {'t_v','x_v','x_fit_v','x_model_v','x_timebased_ref_v','x_mm_v', ...
    'x_peak_mm_v','V_raw_mV','V_a','baseline_mV_v','S_v','rev_v','W_v', ...
    'theta_v','pulsePeakTime_v'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(W, f) && numel(W.(f)) == numel(mask)
        W.(f) = W.(f)(mask);
    end
end
W.pointCount = nnz(mask);
end

function W = decimate_highmap_for_fit_local(W, maxPoints)
n = numel(W.t_v);
if n <= maxPoints
    return;
end
idx = unique(round(linspace(1, n, maxPoints)));
fields = {'t_v','x_v','x_fit_v','x_model_v','x_timebased_ref_v','x_mm_v', ...
    'x_peak_mm_v','V_raw_mV','V_a','baseline_mV_v','S_v','rev_v','W_v', ...
    'theta_v','pulsePeakTime_v'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(W, f) && numel(W.(f)) == n
        W.(f) = W.(f)(idx);
    end
end
W.pointCount = numel(idx);
end

function staticState = fit_static_gap_templates_local(W, responseSurface, GapPrior, cfg)
staticState = repmat(struct('sensorId', NaN, 'gMm', NaN, 'dxMm', NaN, ...
    'x0PriorMm', NaN, 'dxCorrectionMm', NaN, ...
    'weightedRmseMv', inf, 'pointCount', 0), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    sp = get_sensor_prior_local(GapPrior, sid);
    gLo = max(min(responseSurface.gTrainMm), sp.gLowerMm - 0.15);
    gHi = min(max(responseSurface.gTrainMm), sp.gUpperMm + 0.15);
    gGrid = linspace(gLo, gHi, cfg.staticGapGridCount);
    mask = W.S_v == sid & isfinite(W.x_model_v) & isfinite(W.V_a) & isfinite(W.W_v);
    idx = find(mask);
    if numel(idx) > cfg.maxFitPointsPerSensorStatic
        idx = idx(round(linspace(1, numel(idx), cfg.maxFitPointsPerSensorStatic)));
    end
    D = struct('x', W.x_model_v(idx), 'v', W.V_a(idx), 'w', max(W.W_v(idx), 0.05), ...
        'x0PriorMm', sp.x0Mm, 'useLowSpeedXRegistration', cfg.useLowSpeedXRegistration);
    if cfg.useLowSpeedXRegistration
        dxGrid = cfg.staticDxCorrectionGridMm;
    else
        dxGrid = cfg.staticDxGridMm;
    end
    fit = fit_static_one_sensor_local(D, responseSurface, gGrid, dxGrid);
    staticState(is).sensorId = sid;
    staticState(is).gMm = fit.gMm;
    staticState(is).dxMm = fit.x0PriorMm + fit.dxCorrectionMm;
    staticState(is).x0PriorMm = fit.x0PriorMm;
    staticState(is).dxCorrectionMm = fit.dxCorrectionMm;
    staticState(is).weightedRmseMv = fit.weightedRmseMv;
    staticState(is).pointCount = fit.pointCount;
end
end

function fit = fit_static_one_sensor_local(D, responseSurface, gGrid, dxGrid)
fit = struct('gMm', NaN, 'x0PriorMm', D.x0PriorMm, ...
    'dxCorrectionMm', NaN, 'weightedRmseMv', inf, 'pointCount', 0);
for g = gGrid(:).'
    for dxCorr = dxGrid(:).'
        if D.useLowSpeedXRegistration
            xEval = D.x - D.x0PriorMm - dxCorr;
            dxTotal = D.x0PriorMm + dxCorr;
        else
            xEval = D.x - dxCorr;
            dxTotal = dxCorr;
        end
        model = eval_response_surface_local(responseSurface, g, xEval);
        [rmse, nValid] = score_with_sensor_affine_local(model, D.v, D.w);
        if rmse < fit.weightedRmseMv
            fit = struct('gMm', g, 'x0PriorMm', D.x0PriorMm, ...
                'dxCorrectionMm', dxTotal - D.x0PriorMm, ...
                'weightedRmseMv', rmse, 'pointCount', nValid);
        end
    end
end
end

function dyn = fit_dynamic_full_wave_local(W, responseSurface, staticState, cfg, allowedEO)
best = struct('score', inf);
bestAllowed = struct('score', inf);
candidateRows = repmat(struct('EO', NaN, 'frequencyHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxCMm', NaN, 'weightedRmseMv', inf, 'score', inf, ...
    'overshootRmseMm', NaN), numel(cfg.eoCandidates), 1);
row = 0;
for eo = cfg.eoCandidates(:).'
    seed = coarse_dynamic_seed_local(W, responseSurface, staticState, cfg, eo);
    fun = @(p) dynamic_objective_scalar_local(p, eo, W, responseSurface, staticState, cfg);
    opts = optimset('Display', 'off', 'MaxIter', 220, 'MaxFunEvals', 900, ...
        'TolX', 1e-5, 'TolFun', 1e-5);
    pOpt = fminsearch(fun, seed, opts);
    [score, detail] = evaluate_dynamic_objective_local(pOpt, eo, W, responseSurface, staticState, cfg);
    row = row + 1;
    candidateRows(row).EO = eo;
    candidateRows(row).frequencyHz = detail.freqHz;
    candidateRows(row).amplitudeMm = detail.amplitudeMm;
    candidateRows(row).phaseRad = detail.phaseRad;
    candidateRows(row).dxCMm = detail.dxCMm;
    candidateRows(row).weightedRmseMv = detail.weightedRmseMv;
    candidateRows(row).score = score;
    candidateRows(row).overshootRmseMm = detail.overshootRmseMm;
    if score < best.score
        best = detail;
        best.score = score;
        best.EO = eo;
        best.freqHz = eo * W.rotFreqHz;
    end
    if any(eo == allowedEO) && score < bestAllowed.score
        bestAllowed = detail;
        bestAllowed.score = score;
        bestAllowed.EO = eo;
        bestAllowed.freqHz = eo * W.rotFreqHz;
    end
end
candidateRows = candidateRows(1:row);
candidateTable = sortrows(struct2table(candidateRows), {'score', 'EO'}, {'ascend', 'ascend'});
if ~isempty(allowedEO) && isfield(bestAllowed, 'weightedRmseMv') && isfinite(bestAllowed.weightedRmseMv)
    dyn = bestAllowed;
    dyn.eoSelectionSource = 'step04_order_gate';
    dyn.allowedEO = allowedEO(:).';
    dyn.unrestrictedEO = best.EO;
    dyn.unrestrictedWeightedRmseMv = best.weightedRmseMv;
else
    dyn = best;
    dyn.eoSelectionSource = 'unrestricted_waveform_rmse';
    dyn.allowedEO = allowedEO(:).';
    dyn.unrestrictedEO = best.EO;
    dyn.unrestrictedWeightedRmseMv = best.weightedRmseMv;
end
dyn.CandidateTable = candidateTable;
end

function seed = coarse_dynamic_seed_local(W, responseSurface, staticState, cfg, eo)
bestScore = inf;
seed = [0, 0, 0];
for A = cfg.dynamicAmpSeedsMm
    for phi = cfg.dynamicPhiSeedsRad
        p = [A, phi, 0];
        [score, ~] = evaluate_dynamic_objective_local(p, eo, W, responseSurface, staticState, cfg);
        if score < bestScore
            bestScore = score;
            seed = p;
        end
    end
end
end

function score = dynamic_objective_scalar_local(p, eo, W, responseSurface, staticState, cfg)
[score, ~] = evaluate_dynamic_objective_local(p, eo, W, responseSurface, staticState, cfg);
end

function [score, detail] = evaluate_dynamic_objective_local(pRaw, eo, W, responseSurface, staticState, cfg)
p = pRaw(:).';
ampRaw = p(1);
phiRaw = p(2);
dxRaw = p(3);
amp = min(abs(ampRaw), cfg.dynamicAmpLimitMm);
phi = wrap_to_pi_local(phiRaw);
dxC = max(min(dxRaw, cfg.dynamicDxLimitMm), -cfg.dynamicDxLimitMm);
penalty = 1e4 * ((ampRaw - amp)^2 + (dxRaw - dxC)^2);

u = amp .* sin(eo .* W.theta_v(:) + phi);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
vPred = nan(size(W.V_a));
overshoot = zeros(size(W.V_a));
sensorAffine = repmat(struct('sensorId', NaN, 'gain', NaN, 'offsetMv', NaN), numel(staticState), 1);
for is = 1:numel(staticState)
    sid = staticState(is).sensorId;
    mask = W.S_v == sid;
    xRaw = W.x_model_v(mask) - staticState(is).dxMm - dxC - u(mask);
    xEval = min(max(xRaw, xMin), xMax);
    overshoot(mask) = max(xMin - xRaw, 0) + max(xRaw - xMax, 0);
    model = eval_response_surface_local(responseSurface, staticState(is).gMm, xEval);
    vv = W.V_a(mask);
    ww = max(W.W_v(mask), 0.05);
    [predLocal, aff] = apply_affine_local(model, vv, ww);
    vPred(mask) = predLocal;
    sensorAffine(is).sensorId = sid;
    sensorAffine(is).gain = aff.gain;
    sensorAffine(is).offsetMv = aff.offsetMv;
end
valid = isfinite(vPred) & isfinite(W.V_a) & isfinite(W.W_v);
if nnz(valid) < cfg.minValidPoints
    weightedRmse = inf;
    plainRmse = inf;
else
    res = W.V_a(valid) - vPred(valid);
    ww = max(W.W_v(valid), 0.05);
    weightedRmse = sqrt(sum(ww .* res.^2) / max(sum(ww), eps));
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
end
overshootRmse = sqrt(sum(max(W.W_v(valid), 0.05) .* overshoot(valid).^2) / ...
    max(sum(max(W.W_v(valid), 0.05)), eps));
score = weightedRmse + cfg.overshootPenaltyMvPerMm * overshootRmse + penalty;
detail = struct('EO', eo, 'freqHz', eo * W.rotFreqHz, 'amplitudeMm', amp, ...
    'phaseRad', phi, 'dxCMm', dxC, 'weightedRmseMv', weightedRmse, ...
    'plainRmseMv', plainRmse, 'overshootRmseMm', overshootRmse, ...
    'gBySensorMm', [staticState.gMm].', 'dxBySensorMm', [staticState.dxMm].' + dxC, ...
    'sensorAffine', sensorAffine, 'VPred', vPred, 'pointCount', nnz(valid));
end

function [rmse, nValid] = score_with_sensor_affine_local(model, v, w)
[pred, ~] = apply_affine_local(model, v, w);
valid = isfinite(pred) & isfinite(v) & isfinite(w);
nValid = nnz(valid);
if nValid < 3
    rmse = inf;
    return;
end
res = v(valid) - pred(valid);
ww = max(w(valid), 0.05);
rmse = sqrt(sum(ww .* res.^2) / max(sum(ww), eps));
end

function [pred, affine] = apply_affine_local(model, v, w)
pred = nan(size(v));
valid = isfinite(model) & isfinite(v) & isfinite(w);
if nnz(valid) < 3
    affine = struct('gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), 0.05));
beta = (H .* sw) \ (v(valid) .* sw);
pred(valid) = H * beta;
affine = struct('gain', beta(1), 'offsetMv', beta(2));
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function sp = get_sensor_prior_local(GapPrior, sid)
idx = find([GapPrior.sensorPrior.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('Gap prior does not contain CH%d.', sid);
end
sp = GapPrior.sensorPrior(idx);
end

function allowedEO = get_step04_allowed_eo_local(regionTable, tCenter, eoCandidates)
allowedEO = [];
if isempty(regionTable) || height(regionTable) == 0 || ...
        ~all(ismember({'regionStart','regionEnd','regionPeakOrder'}, regionTable.Properties.VariableNames))
    return;
end
cover = regionTable.regionStart <= tCenter & regionTable.regionEnd >= tCenter;
if any(cover)
    regionNow = regionTable(cover, :);
    [~, idx] = max(regionNow.regionPeakScore);
else
    [~, idx] = min(abs(regionTable.regionPeakTime - tCenter));
    regionNow = regionTable(idx, :);
    idx = 1;
end
eo = round(regionNow.regionPeakOrder(idx));
if isfinite(eo) && any(eoCandidates == eo)
    allowedEO = eo;
end
end

function wr = pack_window_result_local(iw, lapRange, Wfull, staticState, dyn)
wr = struct();
wr.window_id = iw;
wr.lapRange = lapRange(:).';
wr.windowCenterTime = mean(Wfull.t_v, 'omitnan');
wr.rotFreqHz = Wfull.rotFreqHz;
wr.rotRpm = Wfull.rotFreqHz * 60;
wr.staticState = staticState;
wr.staticMeanRmseMv = mean([staticState.weightedRmseMv], 'omitnan');
wr.dynamic = dyn;
end

function row = make_trend_row_local(wr)
row = table(wr.window_id, wr.lapRange(1), wr.lapRange(end), wr.windowCenterTime, ...
    wr.rotFreqHz, wr.rotRpm, wr.staticMeanRmseMv, wr.dynamic.freqHz, wr.dynamic.EO, ...
    wr.dynamic.amplitudeMm, mean(wr.dynamic.gBySensorMm, 'omitnan'), ...
    wr.dynamic.weightedRmseMv, wr.dynamic.overshootRmseMm, wr.dynamic.pointCount, ...
    'VariableNames', {'window_id', 'lap_start', 'lap_end', 'window_center_time', ...
    'rot_freq_hz', 'rot_rpm', 'static_mean_rmse_mV', 'frequency_hz', 'EO', ...
    'amplitude_mm', 'mean_gap_mm', 'weighted_rmse_mV', 'overshoot_rmse_mm', 'point_count'});
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
