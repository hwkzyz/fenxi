%% Step06L: estimate shared high-speed equivalent sensor static state.
%
% M0: dg_s = 0
% M1: one common dg shared by all sensors
% M2: CH1/CH3/CH6 have independent dg_s
%
% This step creates an experimental StaticPrior structure for later
% protected multi-frequency decoupling.  The reported dg is an equivalent
% static waveform state, not an independently measured physical clearance.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(thisDir));
gapmfDir = fullfile(rootDir, '20250527_synthetic_spectrum_gap_bias_study');
addpath(gapmfDir);

outDir = fullfile(thisDir, 'outputs');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

sourceFile = fullfile(outDir, ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_foundation_fullbundle_template_nested_compare_v2.mat');
if ~isfile(sourceFile)
    sourceFile = fullfile(outDir, ...
        'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaponly.mat');
end
if ~isfile(sourceFile)
    hits = dir(fullfile(outDir, ...
        'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_*main_gaponly.mat'));
    if ~isempty(hits)
        [~, latestIdx] = max([hits.datenum]);
        sourceFile = fullfile(hits(latestIdx).folder, hits(latestIdx).name);
    end
end
if ~isfile(sourceFile)
    error('Run Step07J first. Missing file: %s', sourceFile);
end

S = load(sourceFile, 'Result');
WindowResult = S.Result.WindowResult;
maxWindows = str2double(getenv('STEP06L_MAX_WINDOWS'));
if ~isfinite(maxWindows) || maxWindows <= 0
    maxWindows = min(4, numel(WindowResult));
end
WindowResult = WindowResult(1:min(maxWindows, numel(WindowResult)));
if numel(WindowResult) < 2
    error('Step06L requires at least two windows for leave-one-window validation.');
end

[Spectrum, spectrumPriorFile] = load_spectrum_prior_local(gapmfDir, 4);
sideBetaWeight = str2double(getenv('STEP06L_SIDE_BETA_WEIGHT'));
if ~isfinite(sideBetaWeight) || sideBetaWeight < 0
    sideBetaWeight = 0;
end
measuredHighLowDgMm = [-0.224400, 0.129150, -0.005296];
opts = struct('outerGridOnly', false, ...
    'outerGridValuesMm', [-0.24, -0.12, -0.06, 0, 0.06, 0.12, 0.24], ...
    'maxInnerIter', 10, 'betaLimitMm', 0.60, ...
    'maxOuterIter', 20, 'maxOuterFunEvals', 80, ...
    'dgLimitMm', 0.25, 'priorWeight', 1, 'dataSigmaMv', [], ...
    'sideBetaWeight', sideBetaWeight, ...
    'maxBackgroundAmplificationVsStrain', 5);

Bundles = cell(numel(WindowResult), 1);
for iw = 1:numel(WindowResult)
    Bundles{iw} = prepare_window_bundle_local(WindowResult(iw).bundle, iw, 450);
end
sensorIds = Bundles{1}.sensorIds(:).';
measuredBySensor = map_sensor_values_local([1 3 6], measuredHighLowDgMm, sensorIds);
opts.outerCandidateDgMm = candidate_dg_starts_local(measuredBySensor);
opts.dataSigmaMv = map_sensor_values_local([1 3 6], [45.67 38.44 47.50], sensorIds);
StaticPrior0 = struct('sensorIds', sensorIds, 'meanDgMm', zeros(1, numel(sensorIds)), ...
    'covDgMm2', eye(numel(sensorIds)) .* 0.05^2, ...
    'modelType', "zero", 'evidenceLabel', "equivalent_static_state");

modelNames = ["S0_fixed_single_eo", "M0_fixed_multifrequency", ...
    "M1_common_gap", "M2_per_sensor_gap"];
methodNames = ["fixed_single_eo", "fixed_multi_frequency", ...
    "common_gap_ablation", "per_sensor_gap_multi_frequency"];

comparisonRows = [];
for im = 1:numel(modelNames)
    trainBundle = concat_bundles_local(Bundles);
    fit = gapmf_fit_nested_model_20250527(trainBundle, trainBundle, Spectrum, ...
        StaticPrior0, methodNames(im), opts);
    comparisonRows = [comparisonRows; make_comparison_row_local(modelNames(im), ...
        "all_windows_train", NaN, fit, trainBundle)]; %#ok<AGROW>
    for ihold = 1:numel(Bundles)
        trainIdx = setdiff(1:numel(Bundles), ihold);
        fitCv = gapmf_fit_nested_model_20250527(concat_bundles_local(Bundles(trainIdx)), ...
            Bundles{ihold}, Spectrum, StaticPrior0, methodNames(im), opts);
        comparisonRows = [comparisonRows; make_comparison_row_local(modelNames(im), ...
            "leave_one_window", ihold, fitCv, Bundles{ihold})]; %#ok<AGROW>
    end
end
ModelComparison = struct2table(comparisonRows);

windowRows = [];
for iw = 1:numel(Bundles)
    fitW = gapmf_fit_nested_model_20250527(Bundles{iw}, Bundles{iw}, Spectrum, ...
        StaticPrior0, "per_sensor_gap_multi_frequency", opts);
    windowRows = [windowRows; make_window_row_local(iw, WindowResult(iw), fitW)]; %#ok<AGROW>
end
WindowTable = struct2table(windowRows);

dgMat = [WindowTable.dg1_mm, WindowTable.dg3_mm, WindowTable.dg6_mm];
diagnosis = diagnose_static_prior_support_local(WindowTable, ModelComparison, opts, Spectrum);
sharedM2DgMm = shared_m2_dg_local(ModelComparison);
if any(~isfinite(sharedM2DgMm))
    sharedM2DgMm = mean(dgMat, 1, 'omitnan');
end
StaticPrior = struct();
StaticPrior.sensorIds = sensorIds;
StaticPrior.diagnosticSharedDgMm = sharedM2DgMm;
StaticPrior.diagnosticWindowMeanDgMm = mean(dgMat, 1, 'omitnan');
StaticPrior.diagnosticMeanDgMm = sharedM2DgMm;
StaticPrior.diagnosticWindowCovDgMm2 = cov(dgMat, 'omitrows');
StaticPrior.diagnosticCovDgMm2 = diagnosis.covDgMm2;
if diagnosis.isUsablePrior
    StaticPrior.meanDgMm = sharedM2DgMm;
    StaticPrior.covDgMm2 = diagnosis.covDgMm2;
else
    StaticPrior.meanDgMm = zeros(1, numel(sensorIds));
    StaticPrior.covDgMm2 = eye(numel(sensorIds)) .* 0.10^2;
end
if any(~isfinite(StaticPrior.covDgMm2), 'all') || rank(StaticPrior.covDgMm2) < numel(sensorIds)
    StaticPrior.covDgMm2 = StaticPrior.covDgMm2 + eye(numel(sensorIds)) .* 1e-4;
end
StaticPrior.diagnosticCi95DgMm = diagnosis.ci95DgMm;
if diagnosis.isUsablePrior
    StaticPrior.ci95DgMm = StaticPrior.diagnosticCi95DgMm;
    StaticPrior.modelType = "per_sensor";
else
    StaticPrior.ci95DgMm = NaN(numel(sensorIds), 2);
    StaticPrior.modelType = "not_supported";
end
StaticPrior.evidenceLabel = "equivalent_static_state";
StaticPrior.evidenceLevel = diagnosis.evidenceLevel;
StaticPrior.smokeTestPassed = diagnosis.smokeTestPassed;
StaticPrior.ciMethod = diagnosis.ciMethod;
StaticPrior.sourceWindows = WindowTable.window_id(:).';
StaticPrior.validationRmseMv = mean(ModelComparison.eval_rmse_mv( ...
    ModelComparison.model == "M2_per_sensor_gap" & ...
    ModelComparison.validation_type == "leave_one_window"), 'omitnan');
StaticPrior.sourceFile = sourceFile;
StaticPrior.note = "Equivalent high-speed static waveform state; not independent physical clearance.";
StaticPrior.validationMeaning = "Leave-one-window validation fixes the static state from training windows and re-estimates the holdout window dynamic coefficients.";
StaticPrior.isUsablePrior = diagnosis.isUsablePrior;
StaticPrior.rejectionReason = diagnosis.reason;
StaticPrior.boundHitFraction = diagnosis.boundHitFraction;
StaticPrior.leaveOneWindowRmseMv = diagnosis.leaveOneWindowRmseMv;
StaticPrior.independentBlockCount = diagnosis.independentBlockCount;
StaticPrior.hasOverlappingWindows = diagnosis.hasOverlappingWindows;
StaticPrior.spectrumPriorFile = spectrumPriorFile;
StaticPrior.fixedBackgroundFreqHz = Spectrum.fixedFreqHz;
StaticPrior.sideBetaWeight = sideBetaWeight;
StaticPrior.independentBlockWindowIds = WindowTable.window_id(diagnosis.independentBlockRows).';
StaticPrior.normalizedObjectivePerPoint = diagnosis.normalizedObjectivePerPoint;
StaticPrior.dataScaleAuditPassed = diagnosis.dataScaleAuditPassed;
StaticPrior.m2HoldoutA14MeanMm = diagnosis.m2HoldoutA14MeanMm;
StaticPrior.m2HoldoutA14StdMm = diagnosis.m2HoldoutA14StdMm;
StaticPrior.m2HoldoutA14Cv = diagnosis.m2HoldoutA14Cv;
StaticPrior.m2HoldoutA14RangeMm = diagnosis.m2HoldoutA14RangeMm;
StaticPrior.m2BackgroundToA14Ratio = diagnosis.m2BackgroundToA14Ratio;
StaticPrior.strainBackgroundRmsRatio = diagnosis.strainBackgroundRmsRatio;
StaticPrior.backgroundAmplificationVsStrain = diagnosis.backgroundAmplificationVsStrain;
StaticPrior.backgroundAmplitudeAuditPassed = diagnosis.backgroundAmplitudeAuditPassed;

windowPath = fullfile(outDir, 'Step06L_SharedStaticState_Window.csv');
comparisonPath = fullfile(outDir, 'Step06L_SharedStaticState_ModelComparison.csv');
priorPath = fullfile(outDir, 'Step06L_SharedStaticState_Prior.mat');
summaryPath = fullfile(outDir, 'Step06L_SharedStaticState_Summary.md');
writetable(WindowTable, windowPath);
writetable(ModelComparison, comparisonPath);
save(priorPath, 'StaticPrior', 'WindowTable', 'ModelComparison', '-v7');
write_summary_local(summaryPath, StaticPrior, ModelComparison);

fprintf('Step06L shared high-speed sensor state completed.\n');
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', windowPath, comparisonPath, priorPath, summaryPath);

function B = prepare_window_bundle_local(B, groupId, maxPoints)
fields = {'X','T','V','W','Theta','S'};
for i = 1:numel(fields)
    B.(fields{i}) = B.(fields{i})(:);
end

B.sensorIds = B.sensorIds(:).';
n = numel(B.X);
if n > maxPoints
    keep = false(n, 1);
    perSensor = max(40, floor(maxPoints / numel(B.sensorIds)));
    for is = 1:numel(B.sensorIds)
        idx = find(B.S(:) == B.sensorIds(is));
        take = idx(unique(round(linspace(1, numel(idx), min(perSensor, numel(idx))))));
        keep(take) = true;
    end
    B = subset_bundle_local(B, keep);
end
B.GroupId = repmat(groupId, numel(B.X), 1);
end

function [Spectrum, priorPath] = load_spectrum_prior_local(gapmfDir, maxComponents)
priorPath = fullfile(gapmfDir, 'output', 'Step11_experimental_spectrum_prior.mat');
if ~isfile(priorPath)
    error('Run Step11 before Step06L. Missing spectrum prior: %s', priorPath);
end
S = load(priorPath, 'prior');
T = S.prior.component_table;
side = T(T.component_type == "side", :);
side = sortrows(side, 'amp_ratio_to_main', 'descend');
nKeep = min(maxComponents, height(side));
Spectrum = struct('targetEO', 14, 'fixedFreqHz', side.freq_hz(1:nKeep).', ...
    'componentPhaseRad', side.phase_rad_at_tref(1:nKeep).', ...
    'componentAmplitudeRatio', side.amp_ratio_to_main(1:nKeep).', ...
    'freqSource', "Step11_strongest_strain_components", ...
    'selectionRule', "largest_strain_amplitude_ratio", ...
    'usePhasePrior', false, 'maxComponents', nKeep);
end

function Bout = concat_bundles_local(Bundles)
Bout = Bundles{1};
fields = {'X','T','V','W','Theta','S','GroupId'};
for i = 1:numel(fields)
    f = fields{i};
    x = cell(numel(Bundles), 1);
    for ib = 1:numel(Bundles)
        x{ib} = Bundles{ib}.(f);
    end
    Bout.(f) = vertcat(x{:});
end
Bout.pointCount = numel(Bout.X);
end

function B = subset_bundle_local(B, mask)
fields = {'X','T','V','W','Theta','S','TRel','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(B, f) && numel(B.(f)) == numel(mask)
        B.(f) = B.(f)(mask);
    end
end
B.pointCount = nnz(mask);
end

function row = make_comparison_row_local(modelName, valType, holdoutWindow, fit, bundle)
row = struct();
row.model = string(modelName);
row.validation_type = string(valType);
row.holdout_window = holdoutWindow;
row.point_count = numel(bundle.X);
row.A14_mean_mm = fit.A14Mm;
row.dg1_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 1);
row.dg3_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 3);
row.dg6_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 6);
row.eval_rmse_mv = fit.evalRmseMv;
row.eval_objective = fit.evalObjective;
row.fit_objective_per_point = fit.fitNormalizedObjectivePerPoint;
[row.background_rms_mm, row.max_background_amp_mm, ...
    row.background_to_A14_ratio] = background_diagnostics_local(fit);
row.fallback_to_nested = logical(fit.fallbackToNested);
end

function [rmsBg, maxBg, ratioBg] = background_diagnostics_local(fit)
T = fit.dynamicCoefficients;
groups = unique(T.group_id, 'stable');
groupRms = nan(numel(groups), 1);
allAmp = [];
for ig = 1:numel(groups)
    freq = unique(T.frequency_hz(T.group_id == groups(ig) & ...
        T.type == "background"), 'stable');
    amp = nan(numel(freq), 1);
    for jf = 1:numel(freq)
        idx = find(T.group_id == groups(ig) & T.type == "background" & ...
            abs(T.frequency_hz - freq(jf)) < 1e-9);
        if numel(idx) >= 2
            amp(jf) = hypot(T.coefficient_mm(idx(1)), T.coefficient_mm(idx(2)));
        end
    end
    groupRms(ig) = sqrt(sum(amp.^2, 'omitnan'));
    allAmp = [allAmp; amp(:)]; %#ok<AGROW>
end
rmsBg = mean(groupRms, 'omitnan');
maxBg = max(allAmp, [], 'omitnan');
ratioBg = rmsBg ./ max(fit.A14EvalMm, eps);
if isempty(allAmp)
    rmsBg = 0;
    maxBg = 0;
    ratioBg = 0;
end
end

function M = candidate_dg_starts_local(measuredDgMm)
measuredDgMm = measuredDgMm(:).';
ns = numel(measuredDgMm);
M = [zeros(1, ns); 0.5 .* measuredDgMm; measuredDgMm; ...
    -0.5 .* measuredDgMm; -measuredDgMm; ...
    repmat([-0.12; -0.06; 0.06; 0.12], 1, ns); ...
    0.06 .* eye(ns); -0.06 .* eye(ns)];
M = unique(M, 'rows', 'stable');
end

function row = make_window_row_local(iw, W, fit)
row = struct();
row.window_index = iw;
row.window_id = get_field_default_local(W, 'windowId', iw);
row.lap_start = min(get_field_default_local(W, 'lapRange', [NaN NaN]));
row.lap_end = max(get_field_default_local(W, 'lapRange', [NaN NaN]));
row.dg1_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 1);
row.dg3_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 3);
row.dg6_mm = get_sensor_value_local(fit.sensorIds, fit.dgMm, 6);
row.A14_mm = fit.A14Mm;
row.eval_rmse_mv = fit.evalRmseMv;
end

function dg = shared_m2_dg_local(ModelComparison)
idx = ModelComparison.model == "M2_per_sensor_gap" & ...
    ModelComparison.validation_type == "all_windows_train";
if ~any(idx)
    dg = [NaN NaN NaN];
    return;
end
k = find(idx, 1, 'first');
dg = [ModelComparison.dg1_mm(k), ModelComparison.dg3_mm(k), ModelComparison.dg6_mm(k)];
end

function v = get_sensor_value_local(sensorIds, values, sid)
idx = find(sensorIds(:).' == sid, 1);
if isempty(idx)
    v = NaN;
else
    v = values(idx);
end
end

function values = map_sensor_values_local(srcIds, srcValues, dstIds)
values = nan(1, numel(dstIds));
for i = 1:numel(dstIds)
    idx = find(srcIds(:).' == dstIds(i), 1);
    if isempty(idx)
        values(i) = NaN;
    else
        values(i) = srcValues(idx);
    end
end
end

function value = get_field_default_local(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function diagnosis = diagnose_static_prior_support_local(WindowTable, ModelComparison, opts, Spectrum)
dgMat = [WindowTable.dg1_mm, WindowTable.dg3_mm, WindowTable.dg6_mm];
tol = 0.995 * opts.dgLimitMm;
boundHit = abs(dgMat) >= tol;
[independentBlockRows, hasOverlappingWindows] = select_independent_blocks_local(WindowTable);
independentBlockCount = numel(independentBlockRows);
dgIndependent = dgMat(independentBlockRows, :);
[covDgMm2, ci95DgMm, ciMethod] = resampling_uncertainty_local( ...
    dgIndependent, independentBlockCount);
models = ["M0_fixed_multifrequency", "M1_common_gap", "M2_per_sensor_gap"];
cv = nan(1, numel(models));
for i = 1:numel(models)
    idx = ModelComparison.model == models(i) & ...
        ModelComparison.validation_type == "leave_one_window";
    cv(i) = mean(ModelComparison.eval_rmse_mv(idx), 'omitnan');
end
m2Improves = cv(3) < min(cv(1:2)) - 0.02 * max(min(cv(1:2)), eps);
fewBoundHits = mean(boundHit(:), 'omitnan') <= 0.25;
sufficientBlocks = independentBlockCount >= 4;
idxM2Train = ModelComparison.model == "M2_per_sensor_gap" & ...
    ModelComparison.validation_type == "all_windows_train";
normalizedObjectivePerPoint = mean(ModelComparison.fit_objective_per_point(idxM2Train), 'omitnan');
dataScaleAuditPassed = isfinite(normalizedObjectivePerPoint) && ...
    normalizedObjectivePerPoint >= 0.5 && normalizedObjectivePerPoint <= 2.0;
idxM2Holdout = ModelComparison.model == "M2_per_sensor_gap" & ...
    ModelComparison.validation_type == "leave_one_window";
m2A14 = ModelComparison.A14_mean_mm(idxM2Holdout);
m2A14Mean = mean(m2A14, 'omitnan');
m2A14Std = std(m2A14, 0, 'omitnan');
m2BackgroundRatio = mean(ModelComparison.background_to_A14_ratio(idxM2Train), 'omitnan');
strainBackgroundRmsRatio = sqrt(sum(Spectrum.componentAmplitudeRatio.^2, 'omitnan'));
backgroundAmplification = m2BackgroundRatio ./ max(strainBackgroundRmsRatio, eps);
backgroundAmplitudeAuditPassed = isfinite(backgroundAmplification) && ...
    backgroundAmplification <= opts.maxBackgroundAmplificationVsStrain;
smokeTestPassed = logical(m2Improves && fewBoundHits);
diagnosis = struct();
diagnosis.leaveOneWindowRmseMv = cv;
diagnosis.boundHitFraction = mean(boundHit(:), 'omitnan');
diagnosis.independentBlockCount = independentBlockCount;
diagnosis.independentBlockRows = independentBlockRows;
diagnosis.hasOverlappingWindows = hasOverlappingWindows;
diagnosis.smokeTestPassed = smokeTestPassed;
diagnosis.covDgMm2 = covDgMm2;
diagnosis.ci95DgMm = ci95DgMm;
diagnosis.ciMethod = ciMethod;
diagnosis.normalizedObjectivePerPoint = normalizedObjectivePerPoint;
diagnosis.dataScaleAuditPassed = dataScaleAuditPassed;
diagnosis.m2HoldoutA14MeanMm = m2A14Mean;
diagnosis.m2HoldoutA14StdMm = m2A14Std;
diagnosis.m2HoldoutA14Cv = m2A14Std ./ max(abs(m2A14Mean), eps);
diagnosis.m2HoldoutA14RangeMm = max(m2A14, [], 'omitnan') - ...
    min(m2A14, [], 'omitnan');
diagnosis.m2BackgroundToA14Ratio = m2BackgroundRatio;
diagnosis.strainBackgroundRmsRatio = strainBackgroundRmsRatio;
diagnosis.backgroundAmplificationVsStrain = backgroundAmplification;
diagnosis.backgroundAmplitudeAuditPassed = backgroundAmplitudeAuditPassed;
diagnosis.isUsablePrior = logical(smokeTestPassed && sufficientBlocks && ...
    dataScaleAuditPassed && backgroundAmplitudeAuditPassed);
if hasOverlappingWindows && height(WindowTable) == 2
    diagnosis.evidenceLevel = "two_overlapping_windows";
elseif ~sufficientBlocks
    diagnosis.evidenceLevel = "insufficient_independent_blocks";
else
    diagnosis.evidenceLevel = "independent_blocks";
end
if diagnosis.isUsablePrior
    diagnosis.reason = "M2 improves validation error, dg estimates do not saturate bounds, and independent block count is sufficient.";
elseif smokeTestPassed && ~sufficientBlocks
    diagnosis.reason = "Smoke test passed, but rejected as usable prior because independent non-overlapping block count is insufficient.";
elseif smokeTestPassed && sufficientBlocks && ~dataScaleAuditPassed
    diagnosis.reason = "Smoke test passed, but rejected because the normalized dynamic residual scale is not calibrated.";
elseif smokeTestPassed && sufficientBlocks && dataScaleAuditPassed && ...
        ~backgroundAmplitudeAuditPassed
    diagnosis.reason = "Smoke test passed, but rejected because fitted background motion is inconsistent with the Step11 strain-spectrum scale audit.";
elseif ~fewBoundHits && ~m2Improves
    diagnosis.reason = "Rejected: dg estimates saturate bounds and M2 does not improve leave-one-window error over M0/M1.";
elseif ~fewBoundHits
    diagnosis.reason = "Rejected: dg estimates saturate bounds.";
else
    diagnosis.reason = "Rejected: M2 does not improve leave-one-window error over M0/M1.";
end
end

function [selectedRows, hasOverlap] = select_independent_blocks_local(WindowTable)
ranges = [WindowTable.lap_start, WindowTable.lap_end];
valid = all(isfinite(ranges), 2);
validRows = find(valid);
ranges = ranges(validRows, :);
hasOverlap = false;
for i = 1:size(ranges, 1)
    for j = i + 1:size(ranges, 1)
        if ranges(i, 1) <= ranges(j, 2) && ranges(j, 1) <= ranges(i, 2)
            hasOverlap = true;
        end
    end
end
if isempty(ranges)
    selectedRows = zeros(0, 1);
    return;
end
[ranges, order] = sortrows(ranges, [2 1]);
validRows = validRows(order);
selectedRows = zeros(0, 1);
lastEnd = -inf;
for i = 1:size(ranges, 1)
    if ranges(i, 1) > lastEnd
        selectedRows(end + 1, 1) = validRows(i); %#ok<AGROW>
        lastEnd = ranges(i, 2);
    end
end
end

function [C, ci, method] = resampling_uncertainty_local(dgMat, independentBlockCount)
ns = size(dgMat, 2);
if independentBlockCount < 4 || size(dgMat, 1) < 4
    C = eye(ns) .* 0.10^2;
    ci = NaN(ns, 2);
    method = "insufficient_independent_blocks";
    return;
end
n = size(dgMat, 1);
jack = nan(n, ns);
for i = 1:n
    jack(i, :) = mean(dgMat(setdiff(1:n, i), :), 1, 'omitnan');
end
center = mean(jack, 1, 'omitnan');
J = jack - center;
C = (n - 1) / n .* (J.' * J);
if any(~isfinite(C), 'all') || rank(C) < ns
    C = C + eye(ns) .* 1e-4;
end
halfWidth = 1.96 .* sqrt(max(diag(C), 0));
meanDg = mean(dgMat, 1, 'omitnan');
ci = [meanDg(:) - halfWidth(:), meanDg(:) + halfWidth(:)];
method = "leave_one_block_jackknife";
end

function write_summary_local(pathText, StaticPrior, ModelComparison)
fid = fopen(pathText, 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot open summary file: %s', pathText);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# Step06L shared high-speed equivalent static state\n\n');
fprintf(fid, 'Terminology: equivalent static waveform state, not independent physical clearance.\n\n');
fprintf(fid, 'Smoke test passed: %d\n\n', StaticPrior.smokeTestPassed);
fprintf(fid, 'Usable as StaticPrior: %d\n\n', StaticPrior.isUsablePrior);
fprintf(fid, 'Evidence level: %s\n\n', StaticPrior.evidenceLevel);
fprintf(fid, 'Independent non-overlapping block count: %d\n\n', StaticPrior.independentBlockCount);
fprintf(fid, 'Fixed background frequencies from Step11 (Hz):');
fprintf(fid, ' %.6f', StaticPrior.fixedBackgroundFreqHz);
fprintf(fid, '\n\n');
fprintf(fid, 'Background coefficient ridge weight: %.6g\n\n', StaticPrior.sideBetaWeight);
fprintf(fid, 'Normalized M2 training objective per valid point: %.6f\n\n', StaticPrior.normalizedObjectivePerPoint);
fprintf(fid, 'Dynamic residual-scale audit passed: %d\n\n', StaticPrior.dataScaleAuditPassed);
fprintf(fid, 'M2 holdout A14 mean/std/CV/range: %.6f / %.6f / %.6f / %.6f mm\n\n', ...
    StaticPrior.m2HoldoutA14MeanMm, StaticPrior.m2HoldoutA14StdMm, ...
    StaticPrior.m2HoldoutA14Cv, StaticPrior.m2HoldoutA14RangeMm);
fprintf(fid, 'M2 background/EO14 ratio: %.6f\n\n', StaticPrior.m2BackgroundToA14Ratio);
fprintf(fid, 'Step11 strain background RMS/main ratio: %.6f\n\n', StaticPrior.strainBackgroundRmsRatio);
fprintf(fid, 'Background amplification versus strain ratio: %.6f\n\n', ...
    StaticPrior.backgroundAmplificationVsStrain);
fprintf(fid, 'Background-amplitude audit passed: %d\n\n', ...
    StaticPrior.backgroundAmplitudeAuditPassed);
fprintf(fid, 'Decision reason: %s\n\n', StaticPrior.rejectionReason);
fprintf(fid, 'Diagnostic shared M2 dg by sensor:\n\n');
for i = 1:numel(StaticPrior.sensorIds)
    fprintf(fid, '- CH%d: %.6f mm\n', ...
        StaticPrior.sensorIds(i), StaticPrior.diagnosticSharedDgMm(i));
end
fprintf(fid, '\nDiagnostic window-mean dg by sensor:\n\n');
for i = 1:numel(StaticPrior.sensorIds)
    fprintf(fid, '- CH%d: %.6f mm\n', ...
        StaticPrior.sensorIds(i), StaticPrior.diagnosticWindowMeanDgMm(i));
end
fprintf(fid, '\nCI method: %s\n', StaticPrior.ciMethod);
if all(isfinite(StaticPrior.diagnosticCi95DgMm), 'all')
    for i = 1:numel(StaticPrior.sensorIds)
        fprintf(fid, '- CH%d CI95 [%.6f, %.6f] mm\n', ...
            StaticPrior.sensorIds(i), StaticPrior.diagnosticCi95DgMm(i, 1), ...
            StaticPrior.diagnosticCi95DgMm(i, 2));
    end
else
    fprintf(fid, '- CI95: NaN because independent validation blocks are insufficient.\n');
end
fprintf(fid, '\nModel comparison should be judged on leave-one-window error, not only training error.\n');
fprintf(fid, 'The leave-one-window error is conditional: the static state is estimated from training windows, while the holdout window dynamic coefficients are re-estimated.\n\n');
models = unique(ModelComparison.model, 'stable');
for i = 1:numel(models)
    idx = ModelComparison.model == models(i) & ModelComparison.validation_type == "leave_one_window";
    fprintf(fid, '- %s leave-one-window mean RMSE: %.6f mV\n', ...
        models(i), mean(ModelComparison.eval_rmse_mv(idx), 'omitnan'));
end
end
