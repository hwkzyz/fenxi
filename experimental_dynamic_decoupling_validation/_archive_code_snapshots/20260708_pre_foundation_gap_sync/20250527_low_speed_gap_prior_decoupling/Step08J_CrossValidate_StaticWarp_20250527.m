%% Step08J: Cross-validate high-speed static warp models
% This diagnostic branch reads Step07J results and refits each static-warp
% model on held-out folds inside every high-speed window. The goal is not to
% add another free final model. It checks whether dmu_s improves validation
% waveform RMSE, and whether dmu_s is stable rather than boundary-seeking.
%
% It does not overwrite Step07J or earlier outputs.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step08j_cross_validate_static_warp');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = C0.bladeId;
cfg.analysisSensors = C0.sensorIds;
cfg.kFold = parse_positive_integer_env_local('STEP08J_KFOLD', 3);
cfg.maxWindows = parse_positive_integer_env_local('STEP08J_MAX_WINDOWS', inf);
cfg.maxTrainPointsPerFold = parse_positive_integer_env_local('STEP08J_MAX_TRAIN_POINTS', inf);
cfg.maxValidPointsPerFold = parse_positive_integer_env_local('STEP08J_MAX_VALID_POINTS', inf);
cfg.vpTopK = parse_positive_integer_env_local('STEP08J_TOP_K_EO', 3);
cfg.amplitudeLimitMm = parse_nonnegative_numeric_env_local('STEP08J_AMP_LIMIT_MM', 0.50);
cfg.dxLimitMm = parse_nonnegative_numeric_env_local('STEP08J_DX_LIMIT_MM', 0.35);
cfg.deltaGapLimitMm = parse_nonnegative_numeric_env_local('STEP08J_DELTA_G_LIMIT_MM', 0.25);
cfg.deltaMuLimit = parse_nonnegative_numeric_env_local('STEP08J_DELTA_MU_LIMIT', 0.030);
cfg.deltaTauLimitMm = parse_nonnegative_numeric_env_local('STEP08J_DELTA_TAU_LIMIT_MM', 0.16);
cfg.staticRegWeightMv = parse_nonnegative_numeric_env_local('STEP08J_STATIC_REG_WEIGHT_MV', 0.75);
cfg.overshootPenaltyWeight = parse_nonnegative_numeric_env_local('STEP08J_OVERSHOOT_PENALTY_WEIGHT', 100);
cfg.weightFloor = 0.05;
cfg.boundHitRatio = 0.95;
cfg.modelNames = {'fixed','gap_only','gap_tilt','gap_tilt_shift'};
cfg.fminOptions = optimset('Display', 'off', ...
    'MaxIter', parse_positive_integer_env_local('STEP08J_MAX_ITER', 220), ...
    'MaxFunEvals', parse_positive_integer_env_local('STEP08J_MAX_FUN_EVALS', 850), ...
    'TolX', 2e-5, 'TolFun', 2e-5);

sensorOverride = strtrim(getenv('STEP08J_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP08J_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = C0.sensorTag;

modelOverride = strtrim(getenv('STEP08J_MODELS'));
if ~isempty(modelOverride)
    cfg.modelNames = strsplit(strrep(modelOverride, ',', ' '));
    cfg.modelNames = cfg.modelNames(~cellfun('isempty', cfg.modelNames));
end

step07jFileOverride = strtrim(getenv('STEP08J_STEP07J_FILE'));
step07jSuffix = sanitize_suffix_local(strtrim(getenv('STEP08J_STEP07J_SUFFIX')));
if ~isempty(step07jFileOverride)
    step07jFile = step07jFileOverride;
    outputSuffix = sanitize_suffix_local(strtrim(getenv('STEP08J_OUTPUT_SUFFIX')));
else
    step07jFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s%s.mat', ...
        C0.dataset, C0.caseTag, step07jSuffix));
    if ~isfile(step07jFile)
        step07jFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s%s.mat', ...
            C0.dataset, cfg.targetBlade, sensorTag, step07jSuffix));
    end
    if ~isfile(step07jFile) && isempty(step07jSuffix)
        step07jFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', ...
            C0.dataset, C0.caseTag));
    end
    if ~isfile(step07jFile) && isempty(step07jSuffix)
        step07jFile = fullfile(outDir, sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaptilt.mat', ...
            C0.dataset, cfg.targetBlade, sensorTag));
    end
    outputSuffix = step07jSuffix;
end
if isempty(outputSuffix)
    outputSuffix = sanitize_suffix_local(strtrim(getenv('STEP08J_OUTPUT_SUFFIX')));
end
if ~isfile(step07jFile)
    error('Run Step07J first. Missing file: %s', step07jFile);
end

S = load(step07jFile, 'Result', 'Trend', 'Summary');
Step07J = S.Result;
if isfield(S, 'Trend')
    Step07Trend = S.Trend;
else
    Step07Trend = table();
end

numWindows = min(numel(Step07J.WindowResult), cfg.maxWindows);
fprintf('\n=== Step08J: cross-validate static warp models ===\n');
fprintf('Input Step07J result: %s\n', step07jFile);
fprintf('Windows: %d, folds: %d, sensors: %s\n', numWindows, cfg.kFold, mat2str(cfg.analysisSensors));
fprintf('Models: %s\n', strjoin(cfg.modelNames, ', '));
fprintf('Fit limits: A %.3f mm, dx %.3f mm, dg %.3f mm, dmu %.3f, dtau %.3f mm\n', ...
    cfg.amplitudeLimitMm, cfg.dxLimitMm, cfg.deltaGapLimitMm, cfg.deltaMuLimit, cfg.deltaTauLimitMm);

cvRows = {};
paramRows = {};
fitStore = struct([]);
fitCount = 0;

for iw = 1:numWindows
    wr = Step07J.WindowResult(iw);
    bundleFull = filter_bundle_sensors_local(wr.bundle, cfg.analysisSensors);
    foldId = make_stratified_fold_ids_local(bundleFull, cfg.kFold);
    fprintf('\nWindow %02d/%02d laps %s, points %d\n', ...
        iw, numWindows, mat2str(wr.lapRange), numel(bundleFull.V));

    for k = 1:cfg.kFold
        trainMask = foldId ~= k;
        validMask = foldId == k;
        trainBundle = decimate_bundle_local(select_bundle_points_local(bundleFull, trainMask), ...
            cfg.maxTrainPointsPerFold);
        validBundle = decimate_bundle_local(select_bundle_points_local(bundleFull, validMask), ...
            cfg.maxValidPointsPerFold);

        for im = 1:numel(cfg.modelNames)
            modeName = cfg.modelNames{im};
            fitTrain = refine_static_warp_fit_local(trainBundle, wr.VPSeedTable, cfg, modeName);
            trainEval = evaluate_static_warp_waveform_local(fitTrain.theta, fitTrain.EO, ...
                trainBundle, cfg, modeName);
            validEval = evaluate_static_warp_waveform_local(fitTrain.theta, fitTrain.EO, ...
                validBundle, cfg, modeName);
            fitTrain.trainEval = rmfield_if_exists_local(trainEval, {'VPred','uMm'});
            fitTrain.validEval = rmfield_if_exists_local(validEval, {'VPred','uMm'});

            fitCount = fitCount + 1;
            fitStore(fitCount).windowId = iw;
            fitStore(fitCount).fold = k;
            fitStore(fitCount).mode = modeName;
            fitStore(fitCount).fit = rmfield_if_exists_local(fitTrain, {'CandidateTable'});

            cvRows{end+1, 1} = make_cv_row_local(iw, k, wr, modeName, fitTrain, trainEval, validEval, cfg); %#ok<SAGROW>
            sensorParam = make_sensor_param_rows_local(iw, k, modeName, fitTrain, bundleFull.sensorIds, cfg);
            for ir = 1:numel(sensorParam)
                paramRows{end+1, 1} = sensorParam{ir}; %#ok<SAGROW>
            end
        end

        foldText = summarize_fold_text_local(cvRows, iw, k);
        fprintf('  fold %d/%d: %s\n', k, cfg.kFold, foldText);
    end
end

CV = vertcat(cvRows{:});
SensorParams = vertcat(paramRows{:});
WindowSummary = build_window_summary_table_local(CV);
Summary = build_summary_table_local(CV, SensorParams, cfg);
Decision = build_decision_table_local(Summary);

Result = struct();
Result.dataset = '20250527';
Result.method = 'cross_validate_static_warp_models';
Result.description = ['Refits Step07J nested static-warp models on deterministic held-out folds. ' ...
    'No global EO consistency constraint is used; EO candidates come from each window VP scan.'];
Result.cfg = cfg;
Result.inputStep07JFile = step07jFile;
Result.Step07Trend = Step07Trend;
Result.CV = CV;
Result.SensorParams = SensorParams;
Result.WindowSummary = WindowSummary;
Result.Summary = Summary;
Result.Decision = Decision;
Result.FitStore = fitStore;

outputTag = [C0.caseTag outputSuffix];
matFile = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_%s_%s.mat', ...
    C0.dataset, outputTag));
cvCsv = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_CV_%s_%s.csv', ...
    C0.dataset, outputTag));
paramCsv = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_Params_%s_%s.csv', ...
    C0.dataset, outputTag));
windowCsv = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_WindowSummary_%s_%s.csv', ...
    C0.dataset, outputTag));
summaryCsv = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_Summary_%s_%s.csv', ...
    C0.dataset, outputTag));
decisionCsv = fullfile(outDir, sprintf('Step08J_CrossValidate_StaticWarp_Decision_%s_%s.csv', ...
    C0.dataset, outputTag));

save(matFile, 'Result', 'CV', 'SensorParams', 'WindowSummary', 'Summary', 'Decision', '-v7.3');
writetable(CV, cvCsv);
writetable(SensorParams, paramCsv);
writetable(WindowSummary, windowCsv);
writetable(Summary, summaryCsv);
writetable(Decision, decisionCsv);

plot_cv_rmse_local(CV, Summary, fullfile(figDir, sprintf('Step08J_CV_RMSE_%s_%s.png', ...
    C0.dataset, outputTag)));
plot_validation_trend_local(WindowSummary, fullfile(figDir, sprintf('Step08J_Validation_Trend_%s_%s.png', ...
    C0.dataset, outputTag)));
plot_dmu_stability_local(SensorParams, cfg, fullfile(figDir, sprintf('Step08J_Dmu_Stability_%s_%s.png', ...
    C0.dataset, outputTag)));

fprintf('\nStep08J complete.\n');
disp(Summary);
disp(Decision);
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, cvCsv, summaryCsv);

%% Local functions
function bundle = filter_bundle_sensors_local(bundleIn, sensorIds)
keepSensor = ismember(bundleIn.S(:), sensorIds(:));
bundle = select_bundle_points_local(bundleIn, keepSensor);
bundle.sensorIds = sensorIds(:).';
end

function foldId = make_stratified_fold_ids_local(bundle, kFold)
n = numel(bundle.V);
foldId = zeros(n, 1);
for is = 1:numel(bundle.sensorIds)
    idx = find(bundle.sensorIndex(:) == is);
    if isempty(idx)
        continue;
    end
    [~, order] = sortrows([bundle.X(idx), bundle.TRel(idx)]);
    idxSorted = idx(order);
    foldId(idxSorted) = 1 + mod((0:numel(idxSorted)-1).', kFold);
end
if any(foldId == 0)
    idx = find(foldId == 0);
    foldId(idx) = 1 + mod((0:numel(idx)-1).', kFold);
end
end

function sub = select_bundle_points_local(bundle, mask)
mask = mask(:);
pointFields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
sub = bundle;
for i = 1:numel(pointFields)
    f = pointFields{i};
    if isfield(bundle, f) && numel(bundle.(f)) == numel(mask)
        sub.(f) = bundle.(f)(mask);
    end
end
sub.pointCount = nnz(mask);
end

function sub = decimate_bundle_local(bundle, maxPoints)
n = numel(bundle.V);
if ~isfinite(maxPoints) || n <= maxPoints
    sub = bundle;
    return;
end
[~, order] = sort(bundle.TRel(:));
pickLocal = unique(round(linspace(1, n, maxPoints))).';
mask = false(n, 1);
mask(order(pickLocal)) = true;
sub = select_bundle_points_local(bundle, mask);
end

function fit = refine_static_warp_fit_local(bundle, seedTable, cfg, modeName)
seedTable = prepare_seed_table_local(seedTable, cfg.vpTopK);
candidate = repmat(struct('EO', NaN, 'freqHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'weightedRmseMv', inf, 'plainRmseMv', inf, ...
    'deltaGapMm', [], 'deltaMuGapPerXMm', [], 'deltaTauMm', [], ...
    'theta', [], 'score', inf), height(seedTable), 1);
best = struct('score', inf);
for i = 1:height(seedTable)
    seed = seedTable(i, :);
    eo = seed.EO;
    theta0 = initial_theta_local(seed.A, seed.phi, seed.dx, modeName, numel(bundle.sensorIds));
    fun = @(theta) bounded_objective_local(theta, eo, bundle, cfg, modeName);
    thetaOpt = fminsearch(fun, theta0, cfg.fminOptions);
    thetaOpt = bound_params_local(thetaOpt, cfg, modeName, numel(bundle.sensorIds));
    detail = evaluate_static_warp_waveform_local(thetaOpt, eo, bundle, cfg, modeName);
    candidate(i).EO = eo;
    candidate(i).freqHz = eo * bundle.rotFreqMeanHz;
    candidate(i).amplitudeMm = detail.amplitudeMm;
    candidate(i).phaseRad = detail.phaseRad;
    candidate(i).dxMm = detail.dxMm;
    candidate(i).weightedRmseMv = detail.weightedRmseMv;
    candidate(i).plainRmseMv = detail.plainRmseMv;
    candidate(i).deltaGapMm = detail.deltaGapMm;
    candidate(i).deltaMuGapPerXMm = detail.deltaMuGapPerXMm;
    candidate(i).deltaTauMm = detail.deltaTauMm;
    candidate(i).theta = thetaOpt;
    candidate(i).score = detail.score;
    if detail.score < best.score
        best = candidate(i);
    end
end
fit = best;
fit.mode = modeName;
fit.CandidateTable = sortrows(struct2table(rmfield_if_exists_local(candidate, ...
    {'deltaGapMm','deltaMuGapPerXMm','deltaTauMm','theta'})), ...
    {'score','weightedRmseMv','EO'}, {'ascend','ascend','ascend'});
end

function seedTable = prepare_seed_table_local(seedTable, topK)
if isempty(seedTable)
    error('VPSeedTable is empty.');
end
names = seedTable.Properties.VariableNames;
if ismember('linearRmseMv', names)
    seedTable = sortrows(seedTable, {'linearRmseMv','EO'}, {'ascend','ascend'});
elseif ismember('fallbackRmseMv', names)
    seedTable = sortrows(seedTable, {'fallbackRmseMv','EO'}, {'ascend','ascend'});
else
    seedTable = sortrows(seedTable, 'EO');
end
seedTable = seedTable(1:min(height(seedTable), topK), :);
end

function score = bounded_objective_local(thetaRaw, eo, bundle, cfg, modeName)
theta = bound_params_local(thetaRaw, cfg, modeName, numel(bundle.sensorIds));
detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName);
score = detail.score + 1e4 * sum((thetaRaw(:) - theta(:)).^2);
end

function theta0 = initial_theta_local(A, phi, dx, modeName, nSensor)
theta0 = [A, phi, dx];
if mode_has_gap_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
if mode_has_tilt_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
if mode_has_shift_local(modeName)
    theta0 = [theta0, zeros(1, nSensor)];
end
end

function theta = bound_params_local(thetaRaw, cfg, modeName, nSensor)
theta = thetaRaw(:).';
theta(1) = min(abs(theta(1)), cfg.amplitudeLimitMm);
theta(2) = wrap_to_pi_local(theta(2));
theta(3) = max(min(theta(3), cfg.dxLimitMm), -cfg.dxLimitMm);
need = 3 + mode_extra_count_local(modeName) * nSensor;
if numel(theta) < need
    theta(numel(theta)+1:need) = 0;
end
theta = theta(1:need);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
if ~isempty(idxDg)
    theta(idxDg) = max(min(theta(idxDg), cfg.deltaGapLimitMm), -cfg.deltaGapLimitMm);
end
if ~isempty(idxDmu)
    theta(idxDmu) = max(min(theta(idxDmu), cfg.deltaMuLimit), -cfg.deltaMuLimit);
end
if ~isempty(idxDtau)
    theta(idxDtau) = max(min(theta(idxDtau), cfg.deltaTauLimitMm), -cfg.deltaTauLimitMm);
end
end

function detail = evaluate_static_warp_waveform_local(theta, eo, bundle, cfg, modeName)
theta = bound_params_local(theta, cfg, modeName, numel(bundle.sensorIds));
A = theta(1);
phi = theta(2);
dx = theta(3);
nSensor = numel(bundle.sensorIds);
[idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor);
dg = zeros(nSensor, 1);
dmu = zeros(nSensor, 1);
dtau = zeros(nSensor, 1);
if ~isempty(idxDg), dg = theta(idxDg).'; end
if ~isempty(idxDmu), dmu = theta(idxDmu).'; end
if ~isempty(idxDtau), dtau = theta(idxDtau).'; end

u = A .* sin(eo .* bundle.Theta + phi);
vPred = nan(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:nSensor
    sid = bundle.sensorIds(is);
    mask = bundle.sensorIndex == is;
    Tpl = get_template_sensor_local(bundle.Template, sid);
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
    xEval = bundle.X(mask) - dx - u(mask);
    [vPred(mask), overshoot(mask)] = eval_static_warp_template_local(Tpl, ...
        bundle.responseSurface, corr, dg(is), dmu(is), dtau(is), xEval);
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
regularizationPenalty = cfg.staticRegWeightMv * ( ...
    sum((dg ./ max(cfg.deltaGapLimitMm, eps)).^2) + ...
    sum((dmu ./ max(cfg.deltaMuLimit, eps)).^2) + ...
    sum((dtau ./ max(cfg.deltaTauLimitMm, eps)).^2)) / max(nSensor, 1);

detail = struct();
detail.EO = eo;
detail.freqHz = eo * bundle.rotFreqMeanHz;
detail.amplitudeMm = A;
detail.phaseRad = phi;
detail.dxMm = dx;
detail.deltaGapMm = dg(:);
detail.deltaMuGapPerXMm = dmu(:);
detail.deltaTauMm = dtau(:);
detail.weightedRmseMv = weightedRmse;
detail.plainRmseMv = plainRmse;
detail.overshootPenalty = overshootPenalty;
detail.staticRegPenalty = regularizationPenalty;
detail.score = weightedRmse + overshootPenalty + regularizationPenalty;
detail.VPred = vPred;
detail.uMm = u;
end

function n = mode_extra_count_local(modeName)
n = double(mode_has_gap_local(modeName)) + double(mode_has_tilt_local(modeName)) + ...
    double(mode_has_shift_local(modeName));
end

function tf = mode_has_gap_local(modeName)
tf = any(strcmpi(modeName, {'gap_only','gap_tilt','gap_tilt_shift'}));
end

function tf = mode_has_tilt_local(modeName)
tf = any(strcmpi(modeName, {'gap_tilt','gap_tilt_shift'}));
end

function tf = mode_has_shift_local(modeName)
tf = strcmpi(modeName, 'gap_tilt_shift');
end

function [idxDg, idxDmu, idxDtau] = mode_param_indices_local(modeName, nSensor)
p = 3;
idxDg = [];
idxDmu = [];
idxDtau = [];
if mode_has_gap_local(modeName)
    idxDg = p + (1:nSensor);
    p = p + nSensor;
end
if mode_has_tilt_local(modeName)
    idxDmu = p + (1:nSensor);
    p = p + nSensor;
end
if mode_has_shift_local(modeName)
    idxDtau = p + (1:nSensor);
end
end

function [model, overshoot] = eval_static_warp_template_local(Tpl, responseSurface, corr, dg, dmu, dtau, xOpr)
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, xLow, 'pchip', NaN);
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
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

function row = make_cv_row_local(iw, fold, wr, modeName, fit, trainEval, validEval, cfg)
row = table(iw, fold, string(modeName), wr.lapRange(1), wr.lapRange(end), wr.rotFreqMeanHz, ...
    fit.EO, fit.freqHz, fit.amplitudeMm, fit.phaseRad, fit.dxMm, ...
    trainEval.weightedRmseMv, validEval.weightedRmseMv, validEval.weightedRmseMv - trainEval.weightedRmseMv, ...
    mean_or_nan_local(fit.deltaGapMm), mean_or_nan_local(fit.deltaMuGapPerXMm), mean_or_nan_local(fit.deltaTauMm), ...
    max_abs_or_nan_local(fit.deltaMuGapPerXMm), max_abs_or_nan_local(fit.deltaTauMm), ...
    max_abs_or_nan_local(fit.deltaMuGapPerXMm) / max(cfg.deltaMuLimit, eps), ...
    'VariableNames', {'window_id','fold','model','lap_start','lap_end','rot_freq_hz', ...
    'EO','frequency_hz','amplitude_mm','phase_rad','dx_mm', ...
    'train_rmse_mV','valid_rmse_mV','generalization_gap_mV', ...
    'mean_delta_gap_mm','mean_delta_mu','mean_delta_tau_mm', ...
    'max_abs_delta_mu','max_abs_delta_tau_mm','delta_mu_bound_ratio'});
end

function rows = make_sensor_param_rows_local(iw, fold, modeName, fit, sensorIds, cfg)
nSensor = numel(sensorIds);
dg = fit.deltaGapMm(:);
dmu = fit.deltaMuGapPerXMm(:);
dtau = fit.deltaTauMm(:);
if isempty(dg), dg = NaN(nSensor, 1); end
if isempty(dmu), dmu = NaN(nSensor, 1); end
if isempty(dtau), dtau = NaN(nSensor, 1); end
dg(end+1:nSensor, 1) = NaN;
dmu(end+1:nSensor, 1) = NaN;
dtau(end+1:nSensor, 1) = NaN;
rows = cell(nSensor, 1);
for is = 1:nSensor
    dmuRatio = abs(dmu(is)) / max(cfg.deltaMuLimit, eps);
    rows{is} = table(iw, fold, string(modeName), sensorIds(is), dg(is), dmu(is), dtau(is), ...
        dmuRatio, dmuRatio >= cfg.boundHitRatio, ...
        'VariableNames', {'window_id','fold','model','sensor_id', ...
        'delta_gap_mm','delta_mu','delta_tau_mm','delta_mu_bound_ratio','delta_mu_near_bound'});
end
end

function textOut = summarize_fold_text_local(cvRows, iw, fold)
T = vertcat(cvRows{:});
T = T(T.window_id == iw & T.fold == fold, :);
parts = cell(height(T), 1);
for i = 1:height(T)
    parts{i} = sprintf('%s %.2f->%.2f mV', char(T.model(i)), T.train_rmse_mV(i), T.valid_rmse_mV(i));
end
textOut = strjoin(parts, ' | ');
end

function WindowSummary = build_window_summary_table_local(CV)
models = unique(CV.model, 'stable');
rows = {};
for iw = unique(CV.window_id).'
    for im = 1:numel(models)
        model = models(im);
        mask = CV.window_id == iw & CV.model == model;
        if ~any(mask)
            continue;
        end
        rows{end+1, 1} = table(iw, model, ...
            mean(CV.EO(mask), 'omitnan'), mode_finite_local(CV.EO(mask)), ...
            mean(CV.frequency_hz(mask), 'omitnan'), mean(CV.amplitude_mm(mask), 'omitnan'), ...
            mean(CV.train_rmse_mV(mask), 'omitnan'), mean(CV.valid_rmse_mV(mask), 'omitnan'), ...
            std(CV.valid_rmse_mV(mask), 'omitnan'), mean(CV.delta_mu_bound_ratio(mask), 'omitnan'), ...
            'VariableNames', {'window_id','model','mean_EO','dominant_EO', ...
            'mean_frequency_hz','mean_amplitude_mm','mean_train_rmse_mV', ...
            'mean_valid_rmse_mV','std_valid_rmse_mV','mean_delta_mu_bound_ratio'}); %#ok<AGROW>
    end
end
WindowSummary = vertcat(rows{:});
end

function Summary = build_summary_table_local(CV, SensorParams, cfg)
models = unique(CV.model, 'stable');
rows = {};
gapOnlyMeanValid = NaN;
gapMask = CV.model == "gap_only";
if any(gapMask)
    gapOnlyMeanValid = mean(CV.valid_rmse_mV(gapMask), 'omitnan');
end
for im = 1:numel(models)
    model = models(im);
    mask = CV.model == model;
    pmask = SensorParams.model == model;
    rows{end+1, 1} = table(model, mode_finite_local(CV.EO(mask)), ...
        mean(CV.frequency_hz(mask), 'omitnan'), std(CV.frequency_hz(mask), 'omitnan'), ...
        mean(CV.amplitude_mm(mask), 'omitnan'), std(CV.amplitude_mm(mask), 'omitnan'), ...
        mean(CV.train_rmse_mV(mask), 'omitnan'), mean(CV.valid_rmse_mV(mask), 'omitnan'), ...
        median(CV.valid_rmse_mV(mask), 'omitnan'), ...
        gapOnlyMeanValid - mean(CV.valid_rmse_mV(mask), 'omitnan'), ...
        mean(CV.generalization_gap_mV(mask), 'omitnan'), ...
        mean(abs(SensorParams.delta_mu(pmask)), 'omitnan'), ...
        std(SensorParams.delta_mu(pmask), 'omitnan'), ...
        mean(SensorParams.delta_mu_near_bound(pmask), 'omitnan'), ...
        mean(CV.delta_mu_bound_ratio(mask), 'omitnan'), ...
        cfg.deltaMuLimit, ...
        'VariableNames', {'model','dominant_EO','mean_frequency_hz','std_frequency_hz', ...
        'mean_amplitude_mm','std_amplitude_mm','mean_train_rmse_mV','mean_valid_rmse_mV', ...
        'median_valid_rmse_mV','valid_improvement_vs_gap_only_mV', ...
        'mean_generalization_gap_mV','mean_abs_delta_mu','std_delta_mu', ...
        'delta_mu_near_bound_rate','mean_delta_mu_bound_ratio','delta_mu_limit'}); %#ok<AGROW>
end
Summary = vertcat(rows{:});
end

function Decision = build_decision_table_local(Summary)
hasGap = any(Summary.model == "gap_only");
hasTilt = any(Summary.model == "gap_tilt");
hasShift = any(Summary.model == "gap_tilt_shift");
gapTiltValidGain = NaN;
gapTiltNearBound = NaN;
shiftExtraGain = NaN;
if hasGap && hasTilt
    gap = Summary(Summary.model == "gap_only", :);
    tilt = Summary(Summary.model == "gap_tilt", :);
    gapTiltValidGain = gap.mean_valid_rmse_mV - tilt.mean_valid_rmse_mV;
    gapTiltNearBound = tilt.delta_mu_near_bound_rate;
end
if hasTilt && hasShift
    tilt = Summary(Summary.model == "gap_tilt", :);
    shift = Summary(Summary.model == "gap_tilt_shift", :);
    shiftExtraGain = tilt.mean_valid_rmse_mV - shift.mean_valid_rmse_mV;
end
if isfinite(gapTiltValidGain) && gapTiltValidGain > 2.0 && gapTiltNearBound < 0.25
    verdict = "gap_tilt passes validation";
elseif isfinite(gapTiltValidGain) && gapTiltValidGain > 2.0
    verdict = "gap_tilt improves validation but dmu hits bounds";
elseif isfinite(gapTiltValidGain) && gapTiltValidGain > 0
    verdict = "gap_tilt gives weak validation gain";
else
    verdict = "gap_tilt not validated";
end
Decision = table(verdict, gapTiltValidGain, gapTiltNearBound, shiftExtraGain, ...
    'VariableNames', {'verdict','gap_tilt_valid_gain_mV', ...
    'gap_tilt_dmu_near_bound_rate','shift_extra_valid_gain_mV'});
end

function plot_cv_rmse_local(CV, Summary, figFile)
models = Summary.model;
style = paper_style_local();
fig = figure('Name', 'Step08J cross-validation RMSE', 'Color', 'w', ...
    'Position', [120, 120, 980, 520]);
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
x = 1:numel(models);
errorbar(x - 0.08, Summary.mean_train_rmse_mV, std_by_model_local(CV, models, 'train_rmse_mV'), ...
    'o-', 'Color', style.blue, 'LineWidth', 1.2, 'DisplayName', 'train');
hold on;
errorbar(x + 0.08, Summary.mean_valid_rmse_mV, std_by_model_local(CV, models, 'valid_rmse_mV'), ...
    's-', 'Color', style.red, 'LineWidth', 1.2, 'DisplayName', 'validation');
set(gca, 'XTick', x, 'XTickLabel', cellstr(models));
xlabel('model');
ylabel('weighted RMSE (mV)');
title('Train/validation RMSE');
legend('Location', 'best');
format_axes_local(gca, style);

nexttile;
bar(x, Summary.valid_improvement_vs_gap_only_mV, 0.55, 'FaceColor', style.green);
yline(0, '--', 'Color', [0.35 0.35 0.35]);
set(gca, 'XTick', x, 'XTickLabel', cellstr(models));
xlabel('model');
ylabel('validation gain vs gap only (mV)');
title('Positive means lower validation RMSE');
format_axes_local(gca, style);

export_figure_local(fig, figFile);
end

function plot_validation_trend_local(WindowSummary, figFile)
models = unique(WindowSummary.model, 'stable');
style = paper_style_local();
colors = model_colors_local(numel(models), style);
fig = figure('Name', 'Step08J validation trend', 'Color', 'w', ...
    'Position', [160, 160, 980, 560]);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on;
for im = 1:numel(models)
    mask = WindowSummary.model == models(im);
    plot(WindowSummary.window_id(mask), WindowSummary.mean_valid_rmse_mV(mask), ...
        'o-', 'Color', colors(im, :), 'LineWidth', 1.1, 'DisplayName', char(models(im)));
end
xlabel('window');
ylabel('validation RMSE (mV)');
title('Held-out waveform error');
legend('Location', 'best');
format_axes_local(gca, style);

nexttile;
hold on;
for im = 1:numel(models)
    mask = WindowSummary.model == models(im);
    plot(WindowSummary.window_id(mask), WindowSummary.mean_amplitude_mm(mask), ...
        'o-', 'Color', colors(im, :), 'LineWidth', 1.1, 'DisplayName', char(models(im)));
end
xlabel('window');
ylabel('amplitude (mm)');
title('Amplitude estimated from training folds');
legend('Location', 'best');
format_axes_local(gca, style);

export_figure_local(fig, figFile);
end

function plot_dmu_stability_local(SensorParams, cfg, figFile)
style = paper_style_local();
models = ["gap_tilt", "gap_tilt_shift"];
models = models(ismember(models, unique(SensorParams.model)));
fig = figure('Name', 'Step08J dmu stability', 'Color', 'w', ...
    'Position', [200, 160, 1020, 560]);
tiledlayout(fig, numel(models), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for im = 1:numel(models)
    nexttile;
    hold on;
    model = models(im);
    T = SensorParams(SensorParams.model == model, :);
    sensorIds = unique(T.sensor_id, 'stable');
    colors = lines(numel(sensorIds));
    for is = 1:numel(sensorIds)
        sid = sensorIds(is);
        rows = aggregate_param_by_window_sensor_local(T, sid);
        plot(rows.window_id, rows.mean_delta_mu, 'o-', 'Color', colors(is, :), ...
            'LineWidth', 1.0, 'DisplayName', sprintf('CH%d', sid));
    end
    yline(cfg.deltaMuLimit, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', '+limit');
    yline(-cfg.deltaMuLimit, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'DisplayName', '-limit');
    yline(0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('window');
    ylabel('delta mu');
    title(sprintf('%s dmu by sensor', char(model)));
    legend('Location', 'best');
    format_axes_local(gca, style);
end
export_figure_local(fig, figFile);
end

function rows = aggregate_param_by_window_sensor_local(T, sensorId)
windows = unique(T.window_id(T.sensor_id == sensorId)).';
rows = table();
for iw = windows
    mask = T.window_id == iw & T.sensor_id == sensorId;
    row = table(iw, mean(T.delta_mu(mask), 'omitnan'), std(T.delta_mu(mask), 'omitnan'), ...
        'VariableNames', {'window_id','mean_delta_mu','std_delta_mu'});
    rows = [rows; row]; %#ok<AGROW>
end
end

function s = std_by_model_local(CV, models, fieldName)
s = NaN(numel(models), 1);
for im = 1:numel(models)
    mask = CV.model == models(im);
    s(im) = std(CV.(fieldName)(mask), 'omitnan');
end
end

function style = paper_style_local()
style = struct();
style.fontName = 'Arial';
style.fontSize = 10;
style.tickFontSize = 9;
style.blue = [0.1216 0.4667 0.7059];
style.red = [0.8392 0.1529 0.1569];
style.green = [0.1725 0.6275 0.1725];
style.orange = [1.0000 0.4980 0.0549];
style.black = [0.10 0.10 0.10];
end

function colors = model_colors_local(n, style)
base = [style.black; style.red; style.green; style.orange; style.blue];
if n <= size(base, 1)
    colors = base(1:n, :);
else
    colors = lines(n);
end
end

function format_axes_local(ax, style)
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', 0.8, 'TickDir', 'in', 'Box', 'on', ...
    'TickLabelInterpreter', 'none');
xlabel(ax, get(get(ax, 'XLabel'), 'String'), 'FontName', style.fontName, 'FontSize', style.fontSize);
ylabel(ax, get(get(ax, 'YLabel'), 'String'), 'FontName', style.fontName, 'FontSize', style.fontSize);
title(ax, get(get(ax, 'Title'), 'String'), 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
lgd = findobj(ancestor(ax, 'figure'), 'Type', 'Legend');
for i = 1:numel(lgd)
    set(lgd(i), 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
        'Interpreter', 'none');
end
end

function export_figure_local(fig, pngFile)
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, pngFile, 'Resolution', 300);
[folder, name] = fileparts(pngFile);
pdfFile = fullfile(folder, [name, '.pdf']);
try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
catch
    warning('Could not export PDF figure: %s', pdfFile);
end
end

function out = rmfield_if_exists_local(in, names)
out = in;
for i = 1:numel(names)
    if isstruct(out) && isfield(out, names{i})
        out = rmfield(out, names{i});
    end
end
end

function y = mean_or_nan_local(x)
if isempty(x)
    y = NaN;
else
    y = mean(x, 'omitnan');
end
end

function y = max_abs_or_nan_local(x)
if isempty(x)
    y = NaN;
else
    y = max(abs(x), [], 'omitnan');
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

function suffix = sanitize_suffix_local(raw)
suffix = strtrim(raw);
if isempty(suffix)
    return;
end
suffix = regexprep(suffix, '[^A-Za-z0-9_\\-]', '_');
if ~startsWith(suffix, '_')
    suffix = ['_', suffix];
end
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
