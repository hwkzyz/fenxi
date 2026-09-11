%% Diagnostic_SingleSync_VoltageFunnel_20251222
% Compare the saved VP Top-3 against an all-EO short-joint funnel. This is
% isolated diagnostics code and does not call or modify the formal Main10.
clear; clc;

thisDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(thisDir));
addpath(fullfile(packageRoot, 'functions', 'gap_aware'));

defaultSavedResult = fullfile(packageRoot, 'diagnostics', ...
    'projection_trust_region_20260729', 'results', ...
    'Main_GapAware_VPTopK_FullWave_20251222_B1_S123_diag_projection_allW_fullpoints.mat');
savedResult = strtrim(getenv('STEP07J_DIAG_RESULT_FILE'));
if isempty(savedResult)
    savedResult = defaultSavedResult;
elseif ~isfile(savedResult)
    savedResult = fullfile(packageRoot, 'diagnostics', ...
        'projection_trust_region_20260729', 'results', savedResult);
end
assert(exist(savedResult, 'file') == 2, 'Saved diagnostic result not found: %s', savedResult);
regionTag = region_tag_local(savedResult);

S = load(savedResult, 'Result');
windows = S.Result.WindowResult;
maxWindows = min(numel(windows), parse_positive_integer_env_local( ...
    'STEP07J_DIAG_MAX_WINDOWS', numel(windows)));
windowIndices = parse_window_ids_env_local('STEP07J_DIAG_WINDOW_IDS', ...
    1:maxWindows, numel(windows));
maxPoints = parse_positive_integer_env_local('STEP07J_DIAG_MAX_POINTS', 1200);

cfg = struct();
cfg.amplitudeLimitMm = 0.50;
cfg.dxLimitMm = 0.35;
cfg.deltaGapLimitMm = 0.25;
cfg.deltaGapPriorScaleMm = 0.25;
cfg.staticRegWeightMv = 0.75;
cfg.continuousFrequencyRefine = true;
cfg.frequencyRefineHalfWidthHz = 2.0;
cfg.funnelStage1Iterations = 1;
cfg.funnelStage2Iterations = 4;
cfg.funnelStage2Count = 12;
cfg.finalCandidateCount = 3;

rows = repmat(empty_row_local(), numel(windowIndices), 1);
detail = repmat(struct('windowId', NaN, 'baselineTable', table(), ...
    'funnelStage1', table(), 'funnelStage2', table(), ...
    'baselineFits', {{}}, 'funnelFits', {{}}), numel(windowIndices), 1);

for ik = 1:numel(windowIndices)
    iw = windowIndices(ik);
    ticWindow = tic;
    bundle = decimate_bundle_local(windows(iw).bundle, maxPoints);
    model = make_model_local(bundle);
    vpTable = windows(iw).VPSeedTable;
    vpTable = vpTable(isfinite(vpTable.EO) & isfinite(vpTable.A) & ...
        isfinite(vpTable.phi) & isfinite(vpTable.dx), :);
    assert(~isempty(vpTable), 'Window %d has no finite VP candidates.', iw);

    nFinal = min(cfg.finalCandidateCount, height(vpTable));
    baselineTable = vpTable(1:nFinal, :);
    baselineFits = run_complete_fits_local(bundle, model, baselineTable, cfg);

    stage1 = run_short_funnel_local(bundle, model, vpTable, ...
        cfg.funnelStage1Iterations, cfg);
    nStage2 = min(cfg.funnelStage2Count, height(stage1));
    stage2 = run_short_funnel_local(bundle, model, stage1(1:nStage2, :), ...
        cfg.funnelStage2Iterations, cfg);
    funnelTable = stage2(1:min(nFinal, height(stage2)), :);
    funnelFits = run_complete_fits_local(bundle, model, funnelTable, cfg);

    [baselineBest, baselineMargin] = select_best_fit_local(baselineFits);
    [funnelBest, funnelMargin] = select_best_fit_local(funnelFits);

    row = empty_row_local();
    row.windowId = windows(iw).windowId;
    row.rotFreqHz = bundle.rotFreqMeanHz;
    row.pointCount = bundle.pointCount;
    row.allEoCount = height(vpTable);
    row.baselineTopEO = join_eo_local(baselineTable.EO);
    row.funnelTopEO = join_eo_local(funnelTable.EO);
    row.top3Overlap = numel(intersect(baselineTable.EO, funnelTable.EO));
    row.baselineBestEO = baselineBest.EO;
    row.funnelBestEO = funnelBest.EO;
    row.baselineFrequencyHz = baselineBest.frequencyHz;
    row.funnelFrequencyHz = funnelBest.frequencyHz;
    row.baselineAmplitudeMm = baselineBest.amplitudeMm;
    row.funnelAmplitudeMm = funnelBest.amplitudeMm;
    row.baselinePlainRmseMv = baselineBest.plainRmseMv;
    row.funnelPlainRmseMv = funnelBest.plainRmseMv;
    row.rmseDeltaFunnelMinusBaselineMv = ...
        funnelBest.plainRmseMv - baselineBest.plainRmseMv;
    row.baselineCandidateMarginMv = baselineMargin;
    row.funnelCandidateMarginMv = funnelMargin;
    row.selectedEoChanged = baselineBest.EO ~= funnelBest.EO;
    row.funnelHitBoundary = logical(funnelBest.hitAmplitudeBoundary || ...
        funnelBest.hitDxBoundary || funnelBest.hitGapBoundary);
    row.runtimeSec = toc(ticWindow);
    rows(ik) = row;

    detail(ik).windowId = row.windowId;
    detail(ik).baselineTable = baselineTable;
    detail(ik).funnelStage1 = stage1;
    detail(ik).funnelStage2 = stage2;
    detail(ik).baselineFits = baselineFits;
    detail(ik).funnelFits = funnelFits;

    fprintf(['W%02d VP[%s] -> EO%d %.4f mV | Funnel[%s] -> EO%d %.4f mV ' ...
        '| delta %.4g mV | %.1f s\n'], row.windowId, row.baselineTopEO, ...
        row.baselineBestEO, row.baselinePlainRmseMv, row.funnelTopEO, ...
        row.funnelBestEO, row.funnelPlainRmseMv, ...
        row.rmseDeltaFunnelMinusBaselineMv, row.runtimeSec);
end

summary = struct2table(rows);
outputTag = output_tag_local(windowIndices, numel(windows));
baseName = ['SingleSync_VoltageFunnel_vs_VP_20251222_', regionTag, outputTag];
csvFile = fullfile(thisDir, [baseName, '.csv']);
matFile = fullfile(thisDir, [baseName, '.mat']);
writetable(summary, csvFile);
save(matFile, 'summary', 'detail', 'cfg', 'savedResult', '-v7.3');
fprintf('Saved %s\nSaved %s\n', csvFile, matFile);

function fits = run_complete_fits_local(bundle, model, candidateTable, cfg)
fits = cell(height(candidateTable), 1);
for i = 1:height(candidateTable)
    candidate = struct('EO', candidateTable.EO(i), 'A', candidateTable.A(i), ...
        'phi', candidateTable.phi(i), 'dx', candidateTable.dx(i));
    fits{i} = step07jcore.optimize_experimental_full_waveform( ...
        bundle, candidate, model, cfg);
end
end

function result = run_short_funnel_local(bundle, model, candidateTable, maxIterations, cfg)
n = height(candidateTable);
out = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, ...
    'shortPlainRmseMv', inf, 'shortExitFlag', NaN, ...
    'shortDeltaGapMm', []), n, 1);
for i = 1:n
    seed = struct('EO', candidateTable.EO(i), 'A', candidateTable.A(i), ...
        'phi', candidateTable.phi(i), 'dx', candidateTable.dx(i));
    if ismember('shortDeltaGapMm', candidateTable.Properties.VariableNames)
        seed.deltaGapMm = candidateTable.shortDeltaGapMm{i};
    else
        seed.deltaGapMm = zeros(numel(bundle.sensorIds), 1);
    end
    fit = short_joint_fit_local(bundle, model, seed, maxIterations, cfg);
    out(i).EO = fit.EO;
    out(i).A = fit.amplitudeMm;
    out(i).phi = fit.phaseRad;
    out(i).dx = fit.dxMm;
    out(i).shortPlainRmseMv = fit.plainRmseMv;
    out(i).shortExitFlag = fit.exitflag;
    out(i).shortDeltaGapMm = fit.deltaGapMm;
end
result = struct2table(out);
result = sortrows(result, {'shortPlainRmseMv', 'EO'}, {'ascend', 'ascend'});
end

function fit = short_joint_fit_local(bundle, model, seed, maxIterations, cfg)
nSensor = numel(bundle.sensorIds);
ampLim = cfg.amplitudeLimitMm;
dxLim = cfg.dxLimitMm;
dgLim = cfg.deltaGapLimitMm;
lb = [0, -pi, -dxLim, -dgLim * ones(1, nSensor)];
ub = [ampLim, pi, dxLim, dgLim * ones(1, nSensor)];
dg0 = reshape(seed.deltaGapMm, 1, []);
if numel(dg0) ~= nSensor
    dg0 = zeros(1, nSensor);
end
z0 = min(max([seed.A, seed.phi, seed.dx, dg0], lb), ub);
opts = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxIterations', maxIterations, 'MaxFunctionEvaluations', 250, ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6);
[z, ~, ~, exitflag] = lsqnonlin(@(q) residual_local( ...
    q, seed.EO, bundle, model, cfg), z0, lb, ub, opts);
[plainRmse, ~] = plain_rmse_local(z, seed.EO, bundle, model);
fit = struct('EO', round(seed.EO), 'amplitudeMm', z(1), ...
    'phaseRad', z(2), 'dxMm', z(3), 'deltaGapMm', z(4:end).', ...
    'plainRmseMv', plainRmse, 'exitflag', exitflag);
end

function r = residual_local(z, eo, bundle, model, cfg)
[~, residual] = plain_rmse_local(z, eo, bundle, model);
dg = z(4:end);
if cfg.staticRegWeightMv > 0
    reg = sqrt(numel(bundle.V)) * cfg.staticRegWeightMv * ...
        (dg(:) / cfg.deltaGapPriorScaleMm) / sqrt(max(numel(dg), 1));
    residual = [residual; reg]; %#ok<AGROW>
end
r = residual;
if isempty(r)
    r = 1e9;
end
end

function [plainRmse, residual] = plain_rmse_local(z, eo, bundle, model)
A = z(1); phi = z(2); dx = z(3); dg = z(4:end);
u = A * sin(eo * bundle.Theta(:) + phi);
vPred = nan(size(bundle.V(:)));
for is = 1:numel(bundle.sensorIds)
    mask = bundle.sensorIndex(:) == is;
    [v, ~] = model(is).evaluate(dg(is), bundle.X(mask) - dx - u(mask));
    vPred(mask) = v(:);
end
valid = isfinite(vPred) & isfinite(bundle.V(:));
residual = vPred(valid) - bundle.V(valid);
if isempty(residual)
    plainRmse = inf;
else
    plainRmse = sqrt(mean(residual .^ 2));
end
end

function [best, margin] = select_best_fit_local(fits)
rmse = cellfun(@(x) x.plainRmseMv, fits);
[~, order] = sort(rmse, 'ascend');
best = fits{order(1)};
if numel(order) >= 2
    margin = rmse(order(2)) - rmse(order(1));
else
    margin = NaN;
end
end

function model = make_model_local(bundle)
nSensor = numel(bundle.sensorIds);
model = repmat(struct('xGrid', [], 'evaluate', []), nSensor, 1);
for is = 1:nSensor
    sid = bundle.sensorIds(is);
    it = find([bundle.Template.Sensor.sensor_id] == sid, 1);
    ic = find([bundle.CorrectedGapLibrary.sensor.sensorId] == sid, 1);
    Tpl = bundle.Template.Sensor(it);
    corr = bundle.CorrectedGapLibrary.sensor(ic);
    rs = bundle.responseSurface;
    model(is).xGrid = Tpl.x_grid(:);
    model(is).evaluate = @(dg, x) local_eval_sensor(Tpl, rs, corr, dg, x);
end
end

function [v, aux] = local_eval_sensor(Tpl, rs, corr, dg, x)
[v, aux0] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x, 0, 0);
h = 1e-4;
[vp, ~] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x + h, 0, 0);
[vm, ~] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x - h, 0, 0);
aux = aux0;
aux.dVoltageDx = (vp - vm) / (2 * h);
end

function b = decimate_bundle_local(b, maxPoints)
if numel(b.X) <= maxPoints
    return;
end
keep = false(numel(b.X), 1);
for is = 1:numel(b.sensorIds)
    ids = find(b.sensorIndex == is);
    nk = min(numel(ids), max(1, round(maxPoints * numel(ids) / numel(b.X))));
    pick = unique(round(linspace(1, numel(ids), nk)));
    keep(ids(pick)) = true;
end
fields = {'X', 'T', 'TRel', 'V', 'W', 'Theta', 'S', 'sensorIndex', 'F0', 'Fx'};
for i = 1:numel(fields)
    name = fields{i};
    if isfield(b, name) && numel(b.(name)) == numel(keep)
        b.(name) = b.(name)(keep);
    end
end
b.pointCount = nnz(keep);
end

function row = empty_row_local()
row = struct('windowId', NaN, 'rotFreqHz', NaN, 'pointCount', NaN, ...
    'allEoCount', NaN, 'baselineTopEO', "", 'funnelTopEO', "", ...
    'top3Overlap', NaN, 'baselineBestEO', NaN, 'funnelBestEO', NaN, ...
    'baselineFrequencyHz', NaN, 'funnelFrequencyHz', NaN, ...
    'baselineAmplitudeMm', NaN, 'funnelAmplitudeMm', NaN, ...
    'baselinePlainRmseMv', NaN, 'funnelPlainRmseMv', NaN, ...
    'rmseDeltaFunnelMinusBaselineMv', NaN, ...
    'baselineCandidateMarginMv', NaN, 'funnelCandidateMarginMv', NaN, ...
    'selectedEoChanged', false, 'funnelHitBoundary', false, 'runtimeSec', NaN);
end

function value = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
else
    value = str2double(raw);
    assert(isfinite(value) && value >= 1 && value == floor(value), ...
        '%s must be a positive integer.', name);
end
end

function ids = parse_window_ids_env_local(name, defaultIds, maxId)
raw = strtrim(getenv(name));
if isempty(raw)
    ids = defaultIds;
    return;
end
ids = unique(round(sscanf(raw, '%f').'), 'stable');
assert(~isempty(ids) && all(isfinite(ids)) && all(ids >= 1) && ...
    all(ids <= maxId), '%s must contain window IDs from 1 to %d.', name, maxId);
end

function tag = output_tag_local(ids, totalCount)
if isequal(ids(:).', 1:totalCount)
    tag = '';
else
    tag = ['_W', char(join(compose('%02d', ids), '_'))];
end
end

function tag = region_tag_local(fileName)
name = lower(string(fileName));
if contains(name, 'r07')
    tag = 'R07';
elseif contains(name, 'r04')
    tag = 'R04';
else
    tag = 'R01';
end
end

function text = join_eo_local(eo)
text = join(string(eo(:).'), ',');
end
