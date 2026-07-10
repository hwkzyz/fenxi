%% Step02: inject controlled synthetic voltage into the experimental bundle.
%
% This script fits fixed, gap_only and gap_tilt-style nested models
% independently. gap_tilt does not use gap_only A/phi/d0 as a prior.

cfg = study_config_20250527();
priorPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_prior.mat');
if ~isfile(priorPath)
    error('Run Step01 first. Missing file: %s', priorPath);
end
if ~isfile(cfg.resultFileS136)
    error('Missing Step5 S136 result: %s', cfg.resultFileS136);
end

P = load(priorPath, 'prior');
prior = P.prior;
S = load(cfg.resultFileS136, 'Result_Struct');
Result = S.Result_Struct;
baseBundle = Result.BestWindow.bundle;

caseNames = string(fieldnames(cfg.ch6ShiftCases));
rows = table();
rng(cfg.syntheticVoltageSeed, 'twister');

for iSet = 1:numel(cfg.referenceSensorSets)
    sensorSet = cfg.referenceSensorSets{iSet};
    sensorTag = cfg.referenceSensorTags(iSet);
    bundleSet = filter_bundle_by_sensors_local(baseBundle, sensorSet);

    for iCase = 1:numel(caseNames)
        caseName = caseNames(iCase);
        caseCfg = cfg.ch6ShiftCases.(caseName);

        % Use the same noise realization across bias cases within one sensor
        % set, so S13 differences do not masquerade as CH6-bias effects.
        synSeed = cfg.syntheticVoltageSeed + 1000 * iSet;
        bundleSyn = make_synthetic_voltage_bundle_local(bundleSet, cfg, prior, caseCfg, synSeed);
        fitBundle = thin_bundle_local(bundleSyn, cfg.maxFitPoints);

        for iMethod = 1:numel(cfg.fitMethods)
            methodName = cfg.fitMethods(iMethod);
            fit = fit_nested_sg_model_local(fitBundle, bundleSyn, methodName, cfg, prior);

            row = table();
            row.case_name = caseName;
            row.sensors = sensorTag;
            row.method = methodName;
            row.A_true_mm = prior.main_amp_tip_mm;
            row.A_id_mm = fit.A_mm;
            row.amp_error_mm = fit.A_mm - prior.main_amp_tip_mm;
            row.amp_error_percent = 100 * (fit.A_mm - prior.main_amp_tip_mm) ./ ...
                max(prior.main_amp_tip_mm, eps);
            row.phi_id_rad = fit.phi_rad;
            row.d0_id_mm = fit.d0_mm;
            row.rmse_v = fit.rmse_v;
            row.weighted_rmse_v = fit.weighted_rmse_v;
            row.true_dg6_mm = caseCfg.dg6_mm;
            row.true_dmu6_mm_per_mm = caseCfg.dmu6_mm_per_mm;
            row.noise_std_v = bundleSyn.synthetic_noise_std_v;
            row.point_count_full = numel(bundleSyn.T);
            row.point_count_fit = numel(fitBundle.T);
            row.fit_dg6_mm = get_sensor_value_local(fit.sensor_ids, fit.dg_mm, 6);
            row.fit_dmu6_mm_per_mm = get_sensor_value_local(fit.sensor_ids, fit.dmu_mm_per_mm, 6);
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end

detailPath = fullfile(cfg.outputDir, 'Step02_voltage_bias_detail.csv');
writetable(rows, detailPath);

summary = summarize_rows_local(rows);
summaryPath = fullfile(cfg.outputDir, 'Step02_voltage_bias_summary.csv');
writetable(summary, summaryPath);

disp('Step02 voltage-bias detail:');
disp(rows(:, {'case_name','sensors','method','A_true_mm','A_id_mm', ...
    'amp_error_mm','amp_error_percent','fit_dg6_mm','fit_dmu6_mm_per_mm'}));
fprintf('Step02 completed.\n  %s\n  %s\n', detailPath, summaryPath);

function bundleOut = make_synthetic_voltage_bundle_local(bundleIn, cfg, prior, caseCfg, seed)
bundleOut = bundleIn;
theta = bundleIn.Theta(:);
x = bundleIn.X(:);
s = bundleIn.S(:);

uTrue = prior.main_amp_tip_mm .* sin(cfg.targetEO .* theta + cfg.truePhaseRad) + ...
    cfg.trueD0Mm;

dgSample = zeros(size(x));
dmuSample = zeros(size(x));
isCh6 = s == 6;
dgSample(isCh6) = caseCfg.dg6_mm;
dmuSample(isCh6) = caseCfg.dmu6_mm_per_mm;
staticShift = dgSample + dmuSample .* (x - bundleIn.xc_pt(:));

xIn = x - uTrue - staticShift;
vClean = eval_super_gaussian_local(xIn, bundleIn);

rangeV = max(vClean) - min(vClean);
noiseStd = caseCfg.noise_scale .* cfg.voltageNoiseStdFractionOfRange .* max(rangeV, eps);
rng(seed, 'twister');
vSyn = vClean + noiseStd .* randn(size(vClean));

bundleOut.V = vSyn;
bundleOut.V_clean_synthetic = vClean;
bundleOut.u_true_mm = uTrue;
bundleOut.static_shift_true_mm = staticShift;
bundleOut.synthetic_noise_std_v = noiseStd;
bundleOut.synthetic_case = caseCfg;
end

function fit = fit_nested_sg_model_local(fitBundle, evalBundle, methodName, cfg, prior)
sensorIds = unique(fitBundle.S(:)).';
freeSensorIds = sensorIds(2:end);
nFree = numel(freeSensorIds);
methodName = string(methodName);

seedAmp = unique([cfg.fitSeedAmpMm(:); prior.main_amp_tip_mm]);
seedPhi = cfg.fitSeedPhiRad(:);
best = struct('obj', inf, 'p', []);

opts = optimset('Display', 'off', 'MaxIter', cfg.fitMaxIter, ...
    'MaxFunEvals', cfg.fitMaxFunEvals, 'TolX', 1e-6, 'TolFun', 1e-8);

for ia = 1:numel(seedAmp)
    for ip = 1:numel(seedPhi)
        p0 = build_initial_parameter_local(seedAmp(ia), seedPhi(ip), nFree, methodName);
        objFun = @(p) nested_objective_local(p, fitBundle, methodName, ...
            sensorIds, freeSensorIds, cfg);
        pOpt = fminsearch(objFun, p0, opts);
        obj = objFun(pOpt);
        if obj < best.obj
            best.obj = obj;
            best.p = pOpt;
        end
    end
end

[~, evalDetail] = nested_objective_local(best.p, evalBundle, methodName, ...
    sensorIds, freeSensorIds, cfg);

fit = struct();
fit.method = methodName;
fit.sensor_ids = sensorIds;
fit.A_mm = evalDetail.A;
fit.phi_rad = evalDetail.phi;
fit.d0_mm = evalDetail.d0;
fit.dg_mm = evalDetail.dgBySensor;
fit.dmu_mm_per_mm = evalDetail.dmuBySensor;
fit.rmse_v = evalDetail.rmse;
fit.weighted_rmse_v = evalDetail.weightedRmse;
fit.objective = evalDetail.objective;
fit.p = best.p;
end

function p0 = build_initial_parameter_local(seedAmp, seedPhi, nFree, methodName)
p0 = [seedAmp, seedPhi, 0];
if methodName == "gap_only" || methodName == "gap_tilt"
    p0 = [p0, zeros(1, nFree)];
end
if methodName == "gap_tilt"
    p0 = [p0, zeros(1, nFree)];
end
end

function [obj, detail] = nested_objective_local(p, bundle, methodName, sensorIds, freeSensorIds, cfg)
[A, phi, d0, dgBySensor, dmuBySensor] = unpack_params_local( ...
    p, methodName, sensorIds, freeSensorIds);

u = A .* sin(cfg.targetEO .* bundle.Theta(:) + phi) + d0;
dgSample = zeros(size(bundle.X(:)));
dmuSample = zeros(size(bundle.X(:)));
for i = 1:numel(sensorIds)
    mask = bundle.S(:) == sensorIds(i);
    dgSample(mask) = dgBySensor(i);
    dmuSample(mask) = dmuBySensor(i);
end
staticShift = dgSample + dmuSample .* (bundle.X(:) - bundle.xc_pt(:));
xIn = bundle.X(:) - u - staticShift;
vPred = eval_super_gaussian_local(xIn, bundle);
res = bundle.V(:) - vPred;
w = max(bundle.W(:), 0);
if ~any(w > 0)
    w = ones(size(res));
end

dataObj = sum(w .* res.^2);
penalty = parameter_penalty_local(A, d0, dgBySensor, dmuBySensor, cfg, numel(res));
obj = dataObj + penalty;

if nargout > 1
    detail = struct();
    detail.A = min(max(A, 0), cfg.fitAmpLimitMm);
    detail.phi = wrap_pi_local(phi);
    detail.d0 = min(max(d0, -cfg.fitD0LimitMm), cfg.fitD0LimitMm);
    detail.dgBySensor = dgBySensor;
    detail.dmuBySensor = dmuBySensor;
    detail.VPred = vPred;
    detail.rmse = sqrt(mean(res.^2, 'omitnan'));
    detail.weightedRmse = sqrt(sum(w .* res.^2) ./ max(sum(w), eps));
    detail.objective = dataObj;
end
end

function [A, phi, d0, dgBySensor, dmuBySensor] = unpack_params_local( ...
    p, methodName, sensorIds, freeSensorIds)
A = p(1);
phi = wrap_pi_local(p(2));
d0 = p(3);
dgBySensor = zeros(1, numel(sensorIds));
dmuBySensor = zeros(1, numel(sensorIds));

idx = 4;
if methodName == "gap_only" || methodName == "gap_tilt"
    for i = 1:numel(freeSensorIds)
        sensorIdx = find(sensorIds == freeSensorIds(i), 1);
        dgBySensor(sensorIdx) = p(idx);
        idx = idx + 1;
    end
end
if methodName == "gap_tilt"
    for i = 1:numel(freeSensorIds)
        sensorIdx = find(sensorIds == freeSensorIds(i), 1);
        dmuBySensor(sensorIdx) = p(idx);
        idx = idx + 1;
    end
end
end

function penalty = parameter_penalty_local(A, d0, dgBySensor, dmuBySensor, cfg, n)
penalty = 0;
penalty = penalty + soft_bound_penalty_local(A, 0, cfg.fitAmpLimitMm, 1);
penalty = penalty + soft_bound_penalty_local(d0, -cfg.fitD0LimitMm, cfg.fitD0LimitMm, 1);
for i = 1:numel(dgBySensor)
    penalty = penalty + soft_bound_penalty_local(dgBySensor(i), ...
        -cfg.fitDgLimitMm, cfg.fitDgLimitMm, 1);
    penalty = penalty + soft_bound_penalty_local(dmuBySensor(i), ...
        -cfg.fitDmuLimitMmPerMm, cfg.fitDmuLimitMmPerMm, 1);
end
if cfg.fitDgRegWeight > 0
    penalty = penalty + n * cfg.fitDgRegWeight * sum((dgBySensor ./ max(cfg.fitDgLimitMm, eps)).^2);
end
if cfg.fitDmuRegWeight > 0
    penalty = penalty + n * cfg.fitDmuRegWeight * sum((dmuBySensor ./ max(cfg.fitDmuLimitMmPerMm, eps)).^2);
end
penalty = penalty * 1e6;
end

function p = soft_bound_penalty_local(x, lo, hi, scale)
if x < lo
    p = ((lo - x) ./ max(scale, eps)).^2;
elseif x > hi
    p = ((x - hi) ./ max(scale, eps)).^2;
else
    p = 0;
end
end

function v = eval_super_gaussian_local(xIn, bundle)
v = bundle.B_pt(:) .* exp(-abs((xIn(:) - bundle.xc_pt(:)) ./ bundle.w_pt(:)).^bundle.n_pt(:)) + ...
    bundle.base_pt(:);
end

function bundleOut = filter_bundle_by_sensors_local(bundleIn, sensorSet)
mask = ismember(bundleIn.S(:), sensorSet(:));
bundleOut = subset_bundle_local(bundleIn, mask);
bundleOut.target_sensors = sensorSet(:).';
end

function bundleOut = thin_bundle_local(bundleIn, maxPoints)
n = numel(bundleIn.T);
if n <= maxPoints
    bundleOut = bundleIn;
    return;
end
idx = unique(round(linspace(1, n, maxPoints))).';
mask = false(n, 1);
mask(idx) = true;
bundleOut = subset_bundle_local(bundleIn, mask);
end

function bundleOut = subset_bundle_local(bundleIn, mask)
bundleOut = bundleIn;
sampleFields = {'X','T','V','S','W','Theta','Y_obs','B_pt','w_pt','n_pt', ...
    'xc_pt','base_pt','V_clean_synthetic','u_true_mm','static_shift_true_mm'};
for i = 1:numel(sampleFields)
    f = sampleFields{i};
    if isfield(bundleIn, f) && numel(bundleIn.(f)) == numel(mask)
        bundleOut.(f) = bundleIn.(f)(mask);
    end
end
bundleOut.point_count = nnz(mask);
end

function val = get_sensor_value_local(sensorIds, values, sensorId)
idx = find(sensorIds == sensorId, 1);
if isempty(idx)
    val = NaN;
else
    val = values(idx);
end
end

function summary = summarize_rows_local(rows)
caseVals = unique(rows.case_name, 'stable');
summary = table();
for ic = 1:numel(caseVals)
    idxCase = rows.case_name == caseVals(ic);
    sensorsVals = unique(rows.sensors(idxCase), 'stable');
    for is = 1:numel(sensorsVals)
        idxSensors = idxCase & rows.sensors == sensorsVals(is);
        methodVals = unique(rows.method(idxSensors), 'stable');
        for im = 1:numel(methodVals)
            idx = idxSensors & rows.method == methodVals(im);
            row = table();
            row.case_name = caseVals(ic);
            row.sensors = sensorsVals(is);
            row.method = methodVals(im);
            row.n = nnz(idx);
            row.mean_A_id_mm = mean(rows.A_id_mm(idx), 'omitnan');
            row.mean_amp_error_mm = mean(rows.amp_error_mm(idx), 'omitnan');
            row.mean_amp_error_percent = mean(rows.amp_error_percent(idx), 'omitnan');
            row.mean_rmse_v = mean(rows.rmse_v(idx), 'omitnan');
            row.mean_fit_dg6_mm = mean(rows.fit_dg6_mm(idx), 'omitnan');
            row.mean_fit_dmu6_mm_per_mm = mean(rows.fit_dmu6_mm_per_mm(idx), 'omitnan');
            summary = [summary; row]; %#ok<AGROW>
        end
    end
end
end

function x = wrap_pi_local(x)
x = mod(x + pi, 2*pi) - pi;
end
