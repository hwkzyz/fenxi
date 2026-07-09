%% Step07B: Full-wave identification with tilt-corrected gap library
% This step uses the corrected library from Step06F. Unlike the static
% checks, it compares the high-speed waveform after applying vibration:
%
%   V_s(t) ~= T_corr,s(g_s, x_s(t) - u(t))
%   u(t) = A sin(EO*theta(t) + phi)
%
% EO is selected independently for each window from that window's waveform
% objective. No Step04 order gate is used by default.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step07b_tilt_corrected_full_wave');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.eoCandidates = 3:24;
cfg.maxFitPointsPerWindow = 420;
cfg.gGridCount = 31;
cfg.ampLimitMm = 0.80;
cfg.dxLimitMm = 0.20;
cfg.ampSeedsMm = [0, 0.15, 0.30, 0.45, 0.60];
cfg.phiSeedsRad = linspace(0, 2*pi, 9);
cfg.phiSeedsRad(end) = [];
cfg.sensorDxGridMm = -0.08:0.04:0.08;
cfg.fixStaticSensorStateBeforeDynamic = true;
cfg.overshootPenaltyMvPerMm = 600;
cfg.minValidPoints = 60;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 180, 'MaxFunEvals', 700, ...
    'TolX', 1e-5, 'TolFun', 1e-5);

sensorOverride = strtrim(getenv('STEP07B_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP07B_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
maxWindowOverride = str2double(strtrim(getenv('STEP07B_MAX_WINDOWS')));
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
correctedLibFile = fullfile(outDir, sprintf('Step06F_TiltCorrected_GapLibrary_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
if ~isfile(highMapFile)
    error('Run Step06 first. Missing file: %s', highMapFile);
end
if ~isfile(correctedLibFile)
    error('Run Step06F first. Missing file: %s', correctedLibFile);
end

H = load(highMapFile, 'highMap');
highMapFull = attach_model_coordinate_local(H.highMap);
C = load(correctedLibFile, 'CorrectedGapLibrary');
CorrectedGapLibrary = C.CorrectedGapLibrary;
responseSurface = CorrectedGapLibrary.responseSurface;
Template = load(CorrectedGapLibrary.templateFile, 'Template');
Template = Template.Template;

windowSpecs = build_sliding_window_specs_local(highMapFull);
if isfinite(maxWindowOverride) && maxWindowOverride > 0
    windowSpecs = windowSpecs(1:min(numel(windowSpecs), floor(maxWindowOverride)));
end

fprintf('\n=== Step07B: tilt-corrected full-wave identification ===\n');
fprintf('HighMap: %s\n', highMapFile);
fprintf('Corrected library: %s\n', correctedLibFile);
fprintf('Windows: %d, sensors: %s\n', numel(windowSpecs), mat2str(cfg.analysisSensors));

WindowResult = struct([]);
trendRows = cell(numel(windowSpecs), 1);
bestIdx = 1;
bestScore = inf;
for iw = 1:numel(windowSpecs)
    Wfull = subset_highmap_by_laps_local(highMapFull, windowSpecs(iw).lapRange);
    W = decimate_highmap_for_fit_local(Wfull, cfg.maxFitPointsPerWindow);
    fprintf('\nWindow %02d/%02d, laps %s, points %d\n', ...
        iw, numel(windowSpecs), mat2str(windowSpecs(iw).lapRange), numel(W.t_v));
    dyn = fit_window_full_wave_local(W, Template, CorrectedGapLibrary, responseSurface, cfg);
    wr = pack_window_result_local(iw, windowSpecs(iw).lapRange, Wfull, dyn);
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
    fprintf('  EO%d, f %.3f Hz, A %.4f mm, RMSE %.3f mV, static-like RMSE %.3f mV\n', ...
        dyn.EO, dyn.freqHz, dyn.amplitudeMm, dyn.weightedRmseMv, dyn.noVibWeightedRmseMv);
end

trendTable = vertcat(trendRows{:});
best = WindowResult(bestIdx).dynamic;
summaryTable = table(bestIdx, best.freqHz, best.EO, best.amplitudeMm, ...
    mean(best.gBySensorMm, 'omitnan'), best.weightedRmseMv, best.noVibWeightedRmseMv, ...
    'VariableNames', {'bestWindowIndex', 'frequencyHz', 'EO', 'amplitudeMm', ...
    'meanGapMm', 'weightedRmseMv', 'noVibWeightedRmseMv'});
sensorTable = table(cfg.analysisSensors(:), best.gBySensorMm(:), best.sensorDxMm(:), ...
    [best.sensorAffine.gain].', [best.sensorAffine.offsetMv].', ...
    'VariableNames', {'sensorId','gHatMm','sensorDxMm','affineGain','affineOffsetMv'});

result = struct();
result.dataset = '20250527';
result.method = 'tilt_corrected_gap_library_full_wave';
result.description = 'Per-window full-wave identification using Step06F tilt-corrected incremental gap library.';
result.cfg = cfg;
result.highMapFile = highMapFile;
result.correctedLibFile = correctedLibFile;
result.WindowResult = WindowResult;
result.Trend = trendTable;
result.BestWindowIndex = bestIdx;
result.BestWindow = WindowResult(bestIdx);
result.dynamic = best;

matFile = fullfile(outDir, sprintf('Step07B_TiltCorrected_FullWave_Identification_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step07B_TiltCorrected_FullWave_Identification_Summary_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
trendCsv = fullfile(outDir, sprintf('Step07B_TiltCorrected_FullWave_Identification_Trend_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
sensorCsv = fullfile(outDir, sprintf('Step07B_TiltCorrected_FullWave_Identification_Sensor_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figTrend = fullfile(figDir, sprintf('Step07B_TiltCorrected_FullWave_Trend_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figBest = fullfile(figDir, sprintf('Step07B_BestWindow_DynamicCollapse_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'result', 'summaryTable', 'sensorTable', '-v7.3');
writetable(summaryTable, csvFile);
writetable(trendTable, trendCsv);
writetable(sensorTable, sensorCsv);
plot_trend_local(trendTable, figTrend);
plot_best_window_dynamic_collapse_local(result.BestWindow, Template, CorrectedGapLibrary, responseSurface, cfg, figBest);

fprintf('\nStep07B complete.\n');
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

function dyn = fit_window_full_wave_local(W, Template, CorrectedGapLibrary, responseSurface, cfg)
staticSensorState = estimate_static_sensor_state_local(W, Template, CorrectedGapLibrary, responseSurface, cfg);
best = struct('score', inf);
candidateRows = repmat(struct('EO', NaN, 'frequencyHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'weightedRmseMv', inf, 'noVibWeightedRmseMv', inf, ...
    'score', inf), numel(cfg.eoCandidates), 1);
row = 0;
for eo = cfg.eoCandidates(:).'
    seed = coarse_seed_local(W, Template, CorrectedGapLibrary, responseSurface, cfg, eo, staticSensorState);
    fun = @(p) dynamic_objective_scalar_local(p, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState);
    pOpt = fminsearch(fun, seed, cfg.fminOptions);
    [score, detail] = evaluate_dynamic_local(pOpt, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState);
    row = row + 1;
    candidateRows(row).EO = eo;
    candidateRows(row).frequencyHz = detail.freqHz;
    candidateRows(row).amplitudeMm = detail.amplitudeMm;
    candidateRows(row).phaseRad = detail.phaseRad;
    candidateRows(row).weightedRmseMv = detail.weightedRmseMv;
    candidateRows(row).noVibWeightedRmseMv = detail.noVibWeightedRmseMv;
    candidateRows(row).score = score;
    if score < best.score
        best = detail;
        best.score = score;
    end
end
best.CandidateTable = sortrows(struct2table(candidateRows(1:row)), {'score','EO'}, {'ascend','ascend'});
dyn = best;
end

function staticSensorState = estimate_static_sensor_state_local(W, Template, CorrectedGapLibrary, responseSurface, cfg)
staticSensorState = repmat(struct('sensorId', NaN, 'gMm', NaN, 'dxMm', NaN), numel(cfg.analysisSensors), 1);
zeroU = zeros(size(W.V_a));
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = W.S_v == sid;
    Tsen = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    xObs = W.x_model_v(mask);
    vObs = W.V_a(mask);
    wObs = max(W.W_v(mask), 0.05);
    best = struct('rmse', inf, 'gMm', NaN, 'dxMm', NaN);
    gGrid = linspace(min(responseSurface.gTrainMm), max(responseSurface.gTrainMm), cfg.gGridCount);
    for g = gGrid(:).'
        for dx = cfg.sensorDxGridMm(:).'
            [modelRaw, ~] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xObs - dx - zeroU(mask));
            [rmse, ~, ~] = affine_rmse_local(vObs, modelRaw, wObs);
            if rmse < best.rmse
                best = struct('rmse', rmse, 'gMm', g, 'dxMm', dx);
            end
        end
    end
    staticSensorState(is).sensorId = sid;
    staticSensorState(is).gMm = best.gMm;
    staticSensorState(is).dxMm = best.dxMm;
end
end

function seed = coarse_seed_local(W, Template, CorrectedGapLibrary, responseSurface, cfg, eo, staticSensorState)
bestScore = inf;
seed = [0, 0];
for A = cfg.ampSeedsMm
    for phi = cfg.phiSeedsRad
        p = [A, phi];
        [score, ~] = evaluate_dynamic_local(p, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState);
        if score < bestScore
            bestScore = score;
            seed = p;
        end
    end
end
end

function score = dynamic_objective_scalar_local(p, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState)
[score, ~] = evaluate_dynamic_local(p, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState);
end

function [score, detail] = evaluate_dynamic_local(pRaw, eo, W, Template, CorrectedGapLibrary, responseSurface, cfg, staticSensorState)
p = pRaw(:).';
ampRaw = p(1);
phiRaw = p(2);
amp = min(abs(ampRaw), cfg.ampLimitMm);
phi = wrap_to_pi_local(phiRaw);
paramPenalty = 1e4 * (ampRaw - amp)^2;
u = amp .* sin(eo .* W.theta_v(:) + phi);

[vPred, modelRaw, noVibPred, sensorFit, overshoot] = evaluate_all_sensors_local( ...
    W, Template, CorrectedGapLibrary, responseSurface, cfg, u, staticSensorState);
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
validNo = isfinite(noVibPred) & isfinite(W.V_a) & isfinite(W.W_v);
resNo = W.V_a(validNo) - noVibPred(validNo);
wNo = max(W.W_v(validNo), 0.05);
noVibRmse = sqrt(sum(wNo .* resNo.^2) / max(sum(wNo), eps));
overshootRmse = sqrt(mean(overshoot(valid).^2, 'omitnan'));
score = weightedRmse + cfg.overshootPenaltyMvPerMm * overshootRmse + paramPenalty;
detail = struct('EO', eo, 'freqHz', eo * W.rotFreqHz, ...
    'amplitudeMm', amp, 'phaseRad', phi, 'weightedRmseMv', weightedRmse, ...
    'plainRmseMv', plainRmse, 'noVibWeightedRmseMv', noVibRmse, ...
    'overshootRmseMm', overshootRmse, 'pointCount', nnz(valid), ...
    'VPred', vPred, 'modelRaw', modelRaw, 'noVibPred', noVibPred, ...
    'uMm', u, 'gBySensorMm', [sensorFit.gMm].', ...
    'sensorDxMm', [sensorFit.dxMm].', 'sensorAffine', [sensorFit.affine].');
end

function [vPred, modelRawAll, noVibPred, sensorFit, overshootAll] = evaluate_all_sensors_local( ...
    W, Template, CorrectedGapLibrary, responseSurface, cfg, u, staticSensorState)
vPred = nan(size(W.V_a));
modelRawAll = nan(size(W.V_a));
noVibPred = nan(size(W.V_a));
overshootAll = zeros(size(W.V_a));
sensorFit = repmat(struct('sensorId', NaN, 'gMm', NaN, 'dxMm', NaN, ...
    'affine', struct('sensorId', NaN, 'gain', NaN, 'offsetMv', NaN)), numel(cfg.analysisSensors), 1);
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = W.S_v == sid;
    Tsen = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    xObs = W.x_model_v(mask);
    vObs = W.V_a(mask);
    wObs = max(W.W_v(mask), 0.05);
    uObs = u(mask);
    stateIdx = find([staticSensorState.sensorId] == sid, 1, 'first');
    g = staticSensorState(stateIdx).gMm;
    dx = staticSensorState(stateIdx).dxMm;
    [modelRaw, overshoot] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xObs - dx - uObs);
    [rmse, pred, affine] = affine_rmse_local(vObs, modelRaw, wObs); %#ok<ASGLU>
    best = struct('rmse', rmse, 'gMm', g, 'dxMm', dx, ...
        'pred', pred, 'modelRaw', modelRaw, 'affine', affine, 'overshoot', overshoot);
    [noVibRaw, ~] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, best.gMm, xObs - best.dxMm);
    noVibPredLocal = apply_affine_coeff_local(noVibRaw, best.affine);
    vPred(mask) = best.pred;
    modelRawAll(mask) = best.modelRaw;
    noVibPred(mask) = noVibPredLocal;
    overshootAll(mask) = best.overshoot;
    sensorFit(is).sensorId = sid;
    sensorFit(is).gMm = best.gMm;
    sensorFit(is).dxMm = best.dxMm;
    best.affine.sensorId = sid;
    sensorFit(is).affine = best.affine;
end
end

function [model, overshoot] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xOpr)
vLow = interp1(Tsen.x_grid(:), (Tsen.v_grid(:) - Tsen.baseline) * 1000, xOpr, 'pchip', NaN);
[Fdyn, overshootDyn] = eval_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, overshootBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(overshootDyn, overshootBase);
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr - tau);
gEffRaw = g0 + mu .* (xOpr - tau);
xMin = min(responseSurface.xGrid(:)); xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:)); gMax = max(responseSurface.gTrainMm(:));
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

function [rmse, pred, affine] = affine_rmse_local(v, model, w)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model) & isfinite(w);
if nnz(valid) < 3
    rmse = inf;
    affine = struct('sensorId', NaN, 'gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
beta = (H .* sw) \ (v(valid) .* sw);
pred(valid) = H * beta;
res = v(valid) - pred(valid);
rmse = sqrt(sum(w(valid) .* res.^2) / max(sum(w(valid)), eps));
affine = struct('sensorId', NaN, 'gain', beta(1), 'offsetMv', beta(2));
end

function pred = apply_affine_coeff_local(model, affine)
pred = nan(size(model));
valid = isfinite(model);
pred(valid) = affine.gain .* model(valid) + affine.offsetMv;
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function wr = pack_window_result_local(iw, lapRange, Wfull, dyn)
wr = struct();
wr.window_id = iw;
wr.lapRange = lapRange(:).';
wr.windowCenterTime = mean(Wfull.t_v, 'omitnan');
wr.rotFreqHz = Wfull.rotFreqHz;
wr.rotRpm = Wfull.rotFreqHz * 60;
wr.highMap = Wfull;
wr.dynamic = dyn;
end

function row = make_trend_row_local(wr)
row = table(wr.window_id, wr.lapRange(1), wr.lapRange(end), wr.windowCenterTime, ...
    wr.rotFreqHz, wr.rotRpm, wr.dynamic.freqHz, wr.dynamic.EO, wr.dynamic.amplitudeMm, ...
    mean(wr.dynamic.gBySensorMm, 'omitnan'), wr.dynamic.weightedRmseMv, ...
    wr.dynamic.noVibWeightedRmseMv, wr.dynamic.overshootRmseMm, wr.dynamic.pointCount, ...
    'VariableNames', {'window_id','lap_start','lap_end','window_center_time', ...
    'rot_freq_hz','rot_rpm','frequency_hz','EO','amplitude_mm','mean_gap_mm', ...
    'weighted_rmse_mV','no_vib_weighted_rmse_mV','overshoot_rmse_mm','point_count'});
end

function plot_trend_local(T, figFile)
fig = figure('Name', 'Step07B tilt-corrected full-wave trend', 'Color', 'w', ...
    'Position', [80, 80, 1300, 820]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; plot(T.window_id, T.EO, 'o-', 'LineWidth', 1.3); grid on; ylabel('EO');
title('Per-window EO');
nexttile; plot(T.window_id, T.amplitude_mm, 'o-', 'LineWidth', 1.3); grid on; ylabel('A (mm)');
title('Amplitude');
nexttile; hold on; grid on;
plot(T.window_id, T.weighted_rmse_mV, 'o-', 'LineWidth', 1.3, 'DisplayName', 'with vibration');
plot(T.window_id, T.no_vib_weighted_rmse_mV, 's--', 'LineWidth', 1.1, 'DisplayName', 'no vibration');
xlabel('window'); ylabel('RMSE (mV)'); title('RMSE comparison'); legend('Location', 'best');
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_best_window_dynamic_collapse_local(bestWindow, Template, CorrectedGapLibrary, responseSurface, cfg, figFile)
W = decimate_highmap_for_fit_local(bestWindow.highMap, cfg.maxFitPointsPerWindow);
dyn = bestWindow.dynamic;
fig = figure('Name', 'Step07B best-window dynamic collapse', 'Color', 'w', ...
    'Position', [50, 50, 1500, 360 * numel(cfg.analysisSensors)]);
tiledlayout(numel(cfg.analysisSensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    mask = W.S_v == sid;
    Tsen = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    g = dyn.gBySensorMm(is);
    dx = dyn.sensorDxMm(is);
    xRaw = W.x_model_v(mask) - dx;
    xComp = xRaw - dyn.uMm(mask);
    [modelNo, ~] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xRaw);
    [modelDyn, ~] = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xComp);
    predNo = apply_affine_coeff_local(modelNo, dyn.sensorAffine(is));
    predDyn = apply_affine_coeff_local(modelDyn, dyn.sensorAffine(is));
    [xSortRaw, ordRaw] = sort(xRaw);
    [xSortComp, ordComp] = sort(xComp);

    nexttile; hold on; grid on; box on;
    plot(xRaw, W.V_a(mask), '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 5, 'DisplayName', 'high-speed samples');
    plot(xSortRaw, predNo(ordRaw), 'r-', 'LineWidth', 1.3, 'DisplayName', 'no-vibration template');
    title(sprintf('CH%d before vibration compensation', sid));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(xComp, W.V_a(mask), '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 5, 'DisplayName', 'samples vs x-u');
    plot(xSortComp, predDyn(ordComp), 'r-', 'LineWidth', 1.3, 'DisplayName', 'dynamic template');
    title(sprintf('CH%d after compensation | g %.3f | dx %.3f', sid, g, dx));
    xlabel('x_{OPR}-u(t) (mm)'); ylabel('mV'); legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
