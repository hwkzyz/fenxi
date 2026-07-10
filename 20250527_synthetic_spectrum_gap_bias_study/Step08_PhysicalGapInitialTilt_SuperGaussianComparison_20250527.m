%% Step08: physical gap/initial-tilt synthetic comparison.
%
% This script replaces the equivalent-displacement toy voltage injection with
% the Step07J physical static-warp response surface:
%
%   V = T_low(x - u) + gain * [F(g0 + dg, mu, x - u) - F(g0, mu, x - u)]
%
% The low-speed template absorbs the initial tilt. Therefore initial tilt by
% itself must not create an amplitude bias when the calibration uses the same
% mu. Bias is tested only when the high-speed gap changes, or when the
% calibrated tilt slope used by the method differs from the truth.

clear; clc; close all;

cfg = study_config_20250527();
cfg = extend_physical_cfg_local(cfg);

priorPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_prior.mat');
if ~isfile(priorPath)
    error('Run Step01 first. Missing file: %s', priorPath);
end
if ~isfile(cfg.physicalStep07JFile)
    error('Missing physical Step07J source file: %s', cfg.physicalStep07JFile);
end

P = load(priorPath, 'prior');
prior = P.prior;
S = load(cfg.physicalStep07JFile, 'Result');
sourceBundle = normalize_physical_bundle_local(S.Result.BestWindow.bundle);
sourceBundle = attach_super_gaussian_kernel_local(sourceBundle, cfg);

caseList = build_physical_case_list_local(sourceBundle, cfg, prior);
rows = repmat(empty_result_row_local(), 0, 1);

for ic = 1:numel(caseList)
    caseSpec = caseList(ic);
    bundleSynFull = make_physical_synthetic_bundle_local( ...
        sourceBundle, cfg, prior, caseSpec);

    for iSet = 1:numel(cfg.referenceSensorSets)
        sensorSet = cfg.referenceSensorSets{iSet};
        sensorTag = cfg.referenceSensorTags(iSet);
        bundleSet = filter_physical_bundle_by_sensors_local(bundleSynFull, sensorSet);
        fitBundle = thin_physical_bundle_local(bundleSet, cfg.physicalFitMaxPoints);

        methodList = methods_for_case_local(caseSpec, cfg);
        for im = 1:numel(methodList)
            methodName = methodList(im);
            fprintf('Step08 %s %s %s\n', caseSpec.case_name, sensorTag, methodName);
            fit = fit_physical_method_local(fitBundle, bundleSet, ...
                methodName, cfg, prior, caseSpec);
            rows(end + 1, 1) = make_result_row_local( ... %#ok<SAGROW>
                caseSpec, sensorTag, methodName, fit, bundleSet, fitBundle, prior);
        end
    end
end

Detail = struct2table(rows);
Summary = summarize_physical_rows_local(Detail);
Qa = build_quality_report_local(Detail, cfg);

detailPath = fullfile(cfg.outputDir, 'Step08_physical_gap_initial_tilt_detail.csv');
summaryPath = fullfile(cfg.outputDir, 'Step08_physical_gap_initial_tilt_summary.csv');
qaPath = fullfile(cfg.outputDir, 'Step08_physical_gap_initial_tilt_qa.csv');
figPath = fullfile(cfg.outputDir, 'Step08_physical_gap_initial_tilt_amplitude_error.png');

writetable(Detail, detailPath);
writetable(Summary, summaryPath);
writetable(Qa, qaPath);
plot_step08_summary_local(Detail, figPath, cfg);

if any(~Qa.pass)
    disp(Qa);
    error('Step08 physical-simulation QA failed. See %s', qaPath);
end

disp('Step08 physical gap/initial-tilt comparison completed.');
disp(Summary(:, {'case_name','sensors','method','mean_A_id_mm', ...
    'mean_amp_error_percent','mean_fit_dg6_mm','mean_fit_dmu6_mm_per_mm'}));
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', ...
    detailPath, summaryPath, qaPath, figPath);

%% Local functions

function cfg = extend_physical_cfg_local(cfg)
rootDir = fileparts(cfg.studyDir);
cfg.physicalStep07JFile = fullfile(rootDir, ...
    'experimental_dynamic_decoupling_validation', ...
    '20250527_low_speed_gap_prior_decoupling', ...
    'outputs', ...
    'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_foundation_fullbundle_template_nested_compare_v2.mat');
if ~isfile(cfg.physicalStep07JFile)
    cfg.physicalStep07JFile = fullfile(rootDir, ...
        'experimental_dynamic_decoupling_validation', ...
        '20250527_low_speed_gap_prior_decoupling', ...
        'outputs', ...
        'Step07J_NestedStaticWarp_VPFullWave_20250527_B1_S136_main_gaptilt.mat');
end

cfg.physicalMethods = ["superGaussian", "fixed_template", "gap_fixed_tilt"];
cfg.physicalDiagnosticMethods = "gap_tilt_residual";
cfg.physicalFitMaxPoints = min(cfg.maxFitPoints, 900);
cfg.physicalFitMaxIter = 120;
cfg.physicalFitMaxFunEvals = 450;
cfg.physicalAmpLimitMm = cfg.fitAmpLimitMm;
cfg.physicalDxLimitMm = 0.30;
cfg.physicalDgLimitMm = cfg.fitDgLimitMm;
cfg.physicalDmuLimitMmPerMm = 0.025;
cfg.physicalWeightFloor = 0.05;
cfg.physicalNoiseStdFractionOfRange = cfg.voltageNoiseStdFractionOfRange;
cfg.physicalVoltageSeed = 20260711;
cfg.physicalReferenceSensorId = 1;
cfg.physicalGap6Mm = 0.050;
cfg.physicalTiltCalErrorMmPerMm = 0.005;
cfg.physicalTiltChangeMmPerMm = 0.005;
cfg.physicalSgFitMaxIter = 500;
cfg.physicalSgFitMaxFunEvals = 1600;
cfg.qaCleanAmpTolMm = 5e-3;
cfg.qaZeroAmpTolMm = 5e-3;
cfg.qaDgTolMm = 8e-3;
end

function methods = methods_for_case_local(caseSpec, cfg)
methods = cfg.physicalMethods;
diagnosticCases = ["stage2_gap6_initial_tilt_noise", ...
    "stage3_gap6_tilt_cal_error_plus", ...
    "stage3_gap6_tilt_cal_error_minus", ...
    "stage4_tilt_change6_noise"];
if ismember(caseSpec.case_name, diagnosticCases)
    methods = [methods, cfg.physicalDiagnosticMethods];
end
end

function bundle = normalize_physical_bundle_local(bundle)
bundle.X = bundle.X(:);
bundle.T = bundle.T(:);
bundle.V = bundle.V(:);
bundle.W = bundle.W(:);
bundle.Theta = bundle.Theta(:);
bundle.S = bundle.S(:);
bundle.sensorIds = bundle.sensorIds(:).';
[present, sensorIndex] = ismember(bundle.S(:), bundle.sensorIds);
if ~all(present)
    error('Physical bundle contains samples from sensors outside sensorIds.');
end
bundle.sensorIndex = sensorIndex(:);
bundle.pointCount = numel(bundle.X);
end

function bundle = attach_super_gaussian_kernel_local(bundle, cfg)
sensorIds = bundle.sensorIds(:).';
kernel = repmat(struct('sensorId', NaN, 'B', NaN, 'xc', NaN, ...
    'w', NaN, 'n', NaN, 'base', NaN, 'templateRmseMv', NaN), ...
    numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    Tpl = get_template_sensor_local(bundle.Template, sensorIds(is));
    kernel(is) = fit_super_gaussian_template_local(Tpl, cfg);
end
bundle.SuperGaussianKernel = kernel;
end

function K = fit_super_gaussian_template_local(Tpl, cfg)
x = Tpl.x_grid(:);
y = (Tpl.v_grid(:) - Tpl.baseline) .* 1000;
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 12
    error('Template CH%d has too few valid points for super-Gaussian fit.', ...
        Tpl.sensor_id);
end

base0 = prctile(y, 5);
[yMax, iMax] = max(y);
B0 = max(yMax - base0, eps);
xc0 = x(iMax);
halfLevel = base0 + 0.5 * B0;
idxHalf = find(y >= halfLevel);
if numel(idxHalf) >= 2
    w0 = max((max(x(idxHalf)) - min(x(idxHalf))) / 2, 0.02);
else
    w0 = max((max(x) - min(x)) / 8, 0.02);
end
n0 = 4;
p0 = [B0, xc0, log(w0), log(n0), base0];

opts = optimset('Display', 'off', 'MaxIter', cfg.physicalSgFitMaxIter, ...
    'MaxFunEvals', cfg.physicalSgFitMaxFunEvals, 'TolX', 1e-8, 'TolFun', 1e-8);
obj = @(pRaw) super_gaussian_template_objective_local(pRaw, x, y);
pOpt = fminsearch(obj, p0, opts);
[~, p] = super_gaussian_template_objective_local(pOpt, x, y);
yFit = eval_super_gaussian_params_local(x, p);

K = struct();
K.sensorId = Tpl.sensor_id;
K.B = p(1);
K.xc = p(2);
K.w = p(3);
K.n = p(4);
K.base = p(5);
K.templateRmseMv = sqrt(mean((y - yFit).^2, 'omitnan'));
end

function [obj, p] = super_gaussian_template_objective_local(pRaw, x, y)
xLo = min(x);
xHi = max(x);
p = pRaw(:).';
p(1) = max(abs(p(1)), eps);
p(2) = min(max(p(2), xLo), xHi);
p(3) = min(max(exp(p(3)), 0.005), max(xHi - xLo, 0.01));
p(4) = min(max(exp(p(4)), 1.05), 14);
p(5) = min(max(p(5), min(y) - range(y)), max(y) + range(y));
yPred = eval_super_gaussian_params_local(x, p);
obj = mean((y - yPred).^2, 'omitnan') + 1e6 * sum((pRaw(:).' - ...
    [p(1), p(2), log(p(3)), log(p(4)), p(5)]).^2);
end

function y = eval_super_gaussian_params_local(x, p)
y = p(1) .* exp(-abs((x(:) - p(2)) ./ p(3)).^p(4)) + p(5);
end

function cases = build_physical_case_list_local(bundle, cfg, prior)
sensorIds = bundle.sensorIds(:).';
baseMu = collect_mu_by_sensor_from_library_local(bundle, sensorIds);
zero = zeros(1, numel(sensorIds));
dg6 = zero;
dmu6 = zero;
calErrPlus = zero;
calErrMinus = zero;
i6 = find(sensorIds == 6, 1);
if ~isempty(i6)
    dg6(i6) = cfg.physicalGap6Mm;
    dmu6(i6) = cfg.physicalTiltChangeMmPerMm;
    calErrPlus(i6) = cfg.physicalTiltCalErrorMmPerMm;
    calErrMinus(i6) = -cfg.physicalTiltCalErrorMmPerMm;
end

cases = repmat(empty_case_local(sensorIds), 0, 1);
cases(end + 1) = make_case_local("stage0_clean", "stage0_baseline", ...
    prior.main_amp_tip_mm, 0, zero, zero, zero, zero, ...
    "No gap, no initial tilt, no noise.");
cases(end + 1) = make_case_local("stage0_noise_only", "stage0_baseline", ...
    prior.main_amp_tip_mm, 1, zero, zero, zero, zero, ...
    "No gap, no initial tilt, with noise.");
cases(end + 1) = make_case_local("qa_zero_amp_clean", "qa", ...
    0, 0, zero, zero, zero, zero, ...
    "Zero vibration check.");
cases(end + 1) = make_case_local("stage1_gap6_clean", "stage1_gap_change", ...
    prior.main_amp_tip_mm, 0, dg6, zero, zero, zero, ...
    "CH6 gap change only, no initial tilt, no noise.");
cases(end + 1) = make_case_local("stage1_gap6_noise", "stage1_gap_change", ...
    prior.main_amp_tip_mm, 1, dg6, zero, zero, zero, ...
    "CH6 gap change only, no initial tilt, with noise.");
cases(end + 1) = make_case_local("stage2_initial_tilt_only_clean", ...
    "stage2_initial_tilt", prior.main_amp_tip_mm, 0, zero, baseMu, baseMu, zero, ...
    "Initial tilt is calibrated; no gap change.");
cases(end + 1) = make_case_local("stage2_gap6_initial_tilt_clean", ...
    "stage2_gap_plus_initial_tilt", prior.main_amp_tip_mm, 0, dg6, ...
    baseMu, baseMu, zero, "CH6 gap plus calibrated initial tilt, no noise.");
cases(end + 1) = make_case_local("stage2_gap6_initial_tilt_noise", ...
    "stage2_gap_plus_initial_tilt", prior.main_amp_tip_mm, 1, dg6, ...
    baseMu, baseMu, zero, "CH6 gap plus calibrated initial tilt, with noise.");
cases(end + 1) = make_case_local("stage3_gap6_tilt_cal_error_plus", ...
    "stage3_tilt_cal_error", prior.main_amp_tip_mm, 1, dg6, ...
    baseMu, baseMu + calErrPlus, zero, "CH6 tilt calibration error +0.005.");
cases(end + 1) = make_case_local("stage3_gap6_tilt_cal_error_minus", ...
    "stage3_tilt_cal_error", prior.main_amp_tip_mm, 1, dg6, ...
    baseMu, baseMu + calErrMinus, zero, "CH6 tilt calibration error -0.005.");
cases(end + 1) = make_case_local("stage4_tilt_change6_noise", ...
    "stage4_tilt_change", prior.main_amp_tip_mm, 1, zero, ...
    baseMu, baseMu, dmu6, "CH6 high-speed residual tilt change only.");

    function c = make_case_local(name, stage, Atrue, noiseScale, dg, muTruth, muCal, dmuTruth, notes)
        c = empty_case_local(sensorIds);
        c.case_name = string(name);
        c.stage = string(stage);
        c.A_true_mm = Atrue;
        c.noise_scale = noiseScale;
        c.dg_by_sensor = dg(:).';
        c.mu_truth_by_sensor = muTruth(:).';
        c.mu_cal_by_sensor = muCal(:).';
        c.dmu_truth_by_sensor = dmuTruth(:).';
        c.notes = string(notes);
    end
end

function c = empty_case_local(sensorIds)
c = struct();
c.case_name = "";
c.stage = "";
c.sensor_ids = sensorIds(:).';
c.A_true_mm = NaN;
c.noise_scale = NaN;
c.dg_by_sensor = zeros(1, numel(sensorIds));
c.mu_truth_by_sensor = zeros(1, numel(sensorIds));
c.mu_cal_by_sensor = zeros(1, numel(sensorIds));
c.dmu_truth_by_sensor = zeros(1, numel(sensorIds));
c.notes = "";
end

function bundleOut = make_physical_synthetic_bundle_local(bundleIn, cfg, prior, caseSpec)
bundleOut = bundleIn;
theta = bundleIn.Theta(:);
uTrue = caseSpec.A_true_mm .* sin(cfg.targetEO .* theta + cfg.truePhaseRad) + ...
    cfg.trueD0Mm;

vClean = nan(size(bundleIn.X(:)));
overshoot = zeros(size(vClean));
for is = 1:numel(bundleIn.sensorIds)
    sid = bundleIn.sensorIds(is);
    mask = bundleIn.sensorIndex == is;
    Tpl = get_template_sensor_local(bundleIn.Template, sid);
    corr = get_corrected_sensor_local(bundleIn.CorrectedGapLibrary, sid);
    corr.muGapPerXMm = caseSpec.mu_truth_by_sensor(is);
    xEval = bundleIn.X(mask) - uTrue(mask);
    [vClean(mask), overshoot(mask)] = eval_static_warp_template_local( ...
        Tpl, bundleIn.responseSurface, corr, caseSpec.dg_by_sensor(is), ...
        caseSpec.dmu_truth_by_sensor(is), 0, xEval);
end

rangeV = max(vClean, [], 'omitnan') - min(vClean, [], 'omitnan');
noiseStd = caseSpec.noise_scale .* cfg.physicalNoiseStdFractionOfRange .* max(rangeV, eps);
rng(cfg.physicalVoltageSeed, 'twister');
noise = noiseStd .* randn(size(vClean));

bundleOut.V = vClean + noise;
bundleOut.V_clean_synthetic = vClean;
bundleOut.synthetic_noise_mv = noise;
bundleOut.synthetic_noise_std_mv = noiseStd;
bundleOut.synthetic_overshoot_mm = overshoot;
bundleOut.u_true_mm = uTrue;
bundleOut.synthetic_case = caseSpec;
end

function fit = fit_physical_method_local(fitBundle, evalBundle, methodName, cfg, prior, caseSpec)
methodName = string(methodName);
sensorIds = fitBundle.sensorIds(:).';
freeMask = sensorIds ~= cfg.physicalReferenceSensorId;
if ~any(freeMask) && numel(sensorIds) > 1
    freeMask(1) = false;
    freeMask(2:end) = true;
end
freeSensorIds = sensorIds(freeMask);
nFree = numel(freeSensorIds);

initialList = build_physical_initial_list_local(methodName, nFree, cfg, prior, caseSpec);
best = struct('obj', inf, 'p', []);
opts = optimset('Display', 'off', 'MaxIter', cfg.physicalFitMaxIter, ...
    'MaxFunEvals', cfg.physicalFitMaxFunEvals, 'TolX', 1e-6, 'TolFun', 1e-7);

for i0 = 1:numel(initialList)
    p0 = initialList{i0};
    objFun = @(pRaw) physical_objective_local(pRaw, fitBundle, methodName, ...
        sensorIds, freeSensorIds, cfg, caseSpec);
    pOpt = fminsearch(objFun, p0, opts);
    obj = objFun(pOpt);
    if obj < best.obj
        best.obj = obj;
        best.p = pOpt;
    end
end

[~, evalDetail] = physical_objective_local(best.p, evalBundle, methodName, ...
    sensorIds, freeSensorIds, cfg, caseSpec);
fit = evalDetail;
fit.method = methodName;
fit.sensor_ids = sensorIds;
fit.free_sensor_ids = freeSensorIds;
fit.p = best.p;
end

function initialList = build_physical_initial_list_local(methodName, nFree, cfg, prior, caseSpec)
methodName = string(methodName);
ampSeeds = unique([caseSpec.A_true_mm, 0.85 * prior.main_amp_tip_mm, ...
    prior.main_amp_tip_mm, 0]);
ampSeeds = ampSeeds(isfinite(ampSeeds) & ampSeeds >= 0);
phaseSeeds = unique(wrap_pi_local([cfg.truePhaseRad, 0]));
initialList = {};
for ia = 1:numel(ampSeeds)
    for ip = 1:numel(phaseSeeds)
        p0 = [ampSeeds(ia), phaseSeeds(ip), 0];
        if method_has_gap_local(methodName)
            p0 = [p0, zeros(1, nFree)]; %#ok<AGROW>
        end
        if method_has_residual_tilt_local(methodName)
            p0 = [p0, zeros(1, nFree)]; %#ok<AGROW>
        end
        initialList{end + 1} = p0; %#ok<AGROW>
    end
end

if method_has_gap_local(methodName) && nFree > 0
    gapSeeds = 0.05;
    for ig = 1:numel(gapSeeds)
        p0 = [prior.main_amp_tip_mm, cfg.truePhaseRad, 0, ...
            repmat(gapSeeds(ig), 1, nFree)];
        if method_has_residual_tilt_local(methodName)
            p0 = [p0, zeros(1, nFree)]; %#ok<AGROW>
        end
        initialList{end + 1} = p0; %#ok<AGROW>
    end
end
end

function [obj, detail] = physical_objective_local(pRaw, bundle, methodName, ...
    sensorIds, freeSensorIds, cfg, caseSpec)
[p, boundPenalty] = bound_physical_params_local(pRaw, methodName, nfree_local(freeSensorIds), cfg);
detail = evaluate_physical_model_local(p, bundle, methodName, ...
    sensorIds, freeSensorIds, cfg, caseSpec);
residualObj = detail.residual_objective_mv2;
obj = residualObj + boundPenalty;
detail.objective = residualObj;
detail.bound_penalty = boundPenalty;
end

function n = nfree_local(freeSensorIds)
n = numel(freeSensorIds);
end

function [p, penalty] = bound_physical_params_local(pRaw, methodName, nFree, cfg)
pRaw = pRaw(:).';
need = 3 + double(method_has_gap_local(methodName)) * nFree + ...
    double(method_has_residual_tilt_local(methodName)) * nFree;
if numel(pRaw) < need
    pRaw(numel(pRaw) + 1:need) = 0;
end
p = pRaw(1:need);
p(1) = min(max(abs(p(1)), 0), cfg.physicalAmpLimitMm);
p(2) = wrap_pi_local(p(2));
p(3) = min(max(p(3), -cfg.physicalDxLimitMm), cfg.physicalDxLimitMm);
idx = 4;
if method_has_gap_local(methodName)
    p(idx:(idx + nFree - 1)) = min(max(p(idx:(idx + nFree - 1)), ...
        -cfg.physicalDgLimitMm), cfg.physicalDgLimitMm);
    idx = idx + nFree;
end
if method_has_residual_tilt_local(methodName)
    p(idx:(idx + nFree - 1)) = min(max(p(idx:(idx + nFree - 1)), ...
        -cfg.physicalDmuLimitMmPerMm), cfg.physicalDmuLimitMmPerMm);
end
penalty = 1e10 .* sum((pRaw(1:need) - p).^2);
end

function detail = evaluate_physical_model_local(p, bundle, methodName, ...
    sensorIds, freeSensorIds, cfg, caseSpec)
methodName = string(methodName);
A = p(1);
phi = p(2);
dx = p(3);
u = A .* sin(cfg.targetEO .* bundle.Theta(:) + phi) + cfg.trueD0Mm;

dgBySensor = zeros(1, numel(sensorIds));
dmuBySensor = zeros(1, numel(sensorIds));
idx = 4;
if method_has_gap_local(methodName)
    for i = 1:numel(freeSensorIds)
        k = find(sensorIds == freeSensorIds(i), 1);
        dgBySensor(k) = p(idx);
        idx = idx + 1;
    end
end
if method_has_residual_tilt_local(methodName)
    for i = 1:numel(freeSensorIds)
        k = find(sensorIds == freeSensorIds(i), 1);
        dmuBySensor(k) = p(idx);
        idx = idx + 1;
    end
end

muCalBySensor = collect_mu_by_sensor_local(sensorIds, caseSpec);
vPred = nan(size(bundle.V(:)));
overshoot = zeros(size(vPred));
for is = 1:numel(sensorIds)
    sid = sensorIds(is);
    mask = bundle.sensorIndex == is;
    xEval = bundle.X(mask) - dx - u(mask);
    if methodName == "superGaussian"
        vPred(mask) = eval_super_gaussian_kernel_local(bundle, sid, xEval);
        overshoot(mask) = 0;
    else
        Tpl = get_template_sensor_local(bundle.Template, sid);
        corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sid);
        corr.muGapPerXMm = muCalBySensor(is);
        [vPred(mask), overshoot(mask)] = eval_static_warp_template_local( ...
            Tpl, bundle.responseSurface, corr, dgBySensor(is), ...
            dmuBySensor(is), 0, xEval);
    end
end

valid = isfinite(vPred) & isfinite(bundle.V) & isfinite(bundle.W);
if nnz(valid) < 8
    residualObj = inf;
    weightedRmse = inf;
    plainRmse = inf;
else
    res = bundle.V(valid) - vPred(valid);
    w = max(bundle.W(valid), cfg.physicalWeightFloor);
    residualObj = sum(w .* res.^2);
    weightedRmse = sqrt(sum(w .* res.^2) ./ max(sum(w), eps));
    plainRmse = sqrt(mean(res.^2, 'omitnan'));
end

detail = struct();
detail.A_mm = A;
detail.phi_rad = wrap_pi_local(phi);
detail.dx_mm = dx;
detail.dg_mm = dgBySensor;
detail.dmu_mm_per_mm = dmuBySensor;
detail.mu_cal_mm_per_mm = muCalBySensor;
detail.VPred = vPred;
detail.u_mm = u;
detail.overshoot_mm = overshoot;
detail.weighted_rmse_mv = weightedRmse;
detail.rmse_mv = plainRmse;
detail.residual_objective_mv2 = residualObj;
end

function y = eval_super_gaussian_kernel_local(bundle, sid, x)
idx = find([bundle.SuperGaussianKernel.sensorId] == sid, 1);
if isempty(idx)
    error('Super-Gaussian kernel does not contain CH%d.', sid);
end
K = bundle.SuperGaussianKernel(idx);
p = [K.B, K.xc, K.w, K.n, K.base];
y = eval_super_gaussian_params_local(x, p);
end

function tf = method_has_gap_local(methodName)
methodName = string(methodName);
tf = methodName == "gap_fixed_tilt" || methodName == "gap_tilt_residual";
end

function tf = method_has_residual_tilt_local(methodName)
methodName = string(methodName);
tf = methodName == "gap_tilt_residual";
end

function bundleOut = filter_physical_bundle_by_sensors_local(bundleIn, sensorSet)
sensorSet = sensorSet(:).';
mask = ismember(bundleIn.S(:), sensorSet);
bundleOut = subset_physical_bundle_local(bundleIn, mask);
bundleOut.sensorIds = sensorSet;
[present, sensorIndex] = ismember(bundleOut.S(:), sensorSet);
if ~all(present)
    error('Filtered bundle contains samples outside the selected sensor set.');
end
bundleOut.sensorIndex = sensorIndex(:);
bundleOut.target_sensors = sensorSet;
end

function bundleOut = thin_physical_bundle_local(bundleIn, maxPoints)
n = numel(bundleIn.X);
if n <= maxPoints
    bundleOut = bundleIn;
    return;
end
idx = unique(round(linspace(1, n, maxPoints))).';
mask = false(n, 1);
mask(idx) = true;
bundleOut = subset_physical_bundle_local(bundleIn, mask);
end

function bundleOut = subset_physical_bundle_local(bundleIn, mask)
bundleOut = bundleIn;
sampleFields = {'X','T','TRel','V','W','Theta','S','sensorIndex', ...
    'F0','Fx','V_clean_synthetic','synthetic_noise_mv', ...
    'synthetic_overshoot_mm','u_true_mm'};
for i = 1:numel(sampleFields)
    f = sampleFields{i};
    if isfield(bundleIn, f) && numel(bundleIn.(f)) == numel(mask)
        bundleOut.(f) = bundleIn.(f)(mask);
    end
end
bundleOut.pointCount = nnz(mask);
end

function row = empty_result_row_local()
row = struct();
row.case_name = "";
row.stage = "";
row.sensors = "";
row.method = "";
row.A_true_mm = NaN;
row.A_id_mm = NaN;
row.amp_error_mm = NaN;
row.amp_error_percent = NaN;
row.phi_id_rad = NaN;
row.dx_id_mm = NaN;
row.rmse_mv = NaN;
row.weighted_rmse_mv = NaN;
row.true_dg6_mm = NaN;
row.true_mu6_mm_per_mm = NaN;
row.cal_mu6_mm_per_mm = NaN;
row.true_dmu6_mm_per_mm = NaN;
row.fit_dg6_mm = NaN;
row.fit_dmu6_mm_per_mm = NaN;
row.noise_std_mv = NaN;
row.point_count_full = NaN;
row.point_count_fit = NaN;
row.sg_template_rmse_mean_mv = NaN;
row.notes = "";
end

function row = make_result_row_local(caseSpec, sensorTag, methodName, ...
    fit, bundleSet, fitBundle, prior)
row = empty_result_row_local();
sensorIds = fit.sensor_ids;
row.case_name = caseSpec.case_name;
row.stage = caseSpec.stage;
row.sensors = sensorTag;
row.method = string(methodName);
row.A_true_mm = caseSpec.A_true_mm;
row.A_id_mm = fit.A_mm;
row.amp_error_mm = fit.A_mm - caseSpec.A_true_mm;
if abs(caseSpec.A_true_mm) > eps
    row.amp_error_percent = 100 .* row.amp_error_mm ./ caseSpec.A_true_mm;
else
    row.amp_error_percent = NaN;
end
row.phi_id_rad = fit.phi_rad;
row.dx_id_mm = fit.dx_mm;
row.rmse_mv = fit.rmse_mv;
row.weighted_rmse_mv = fit.weighted_rmse_mv;
row.true_dg6_mm = get_sensor_value_local(caseSpec.sensor_ids, caseSpec.dg_by_sensor, 6);
row.true_mu6_mm_per_mm = get_sensor_value_local(caseSpec.sensor_ids, ...
    caseSpec.mu_truth_by_sensor, 6);
row.cal_mu6_mm_per_mm = get_sensor_value_local(caseSpec.sensor_ids, ...
    caseSpec.mu_cal_by_sensor, 6);
row.true_dmu6_mm_per_mm = get_sensor_value_local(caseSpec.sensor_ids, ...
    caseSpec.dmu_truth_by_sensor, 6);
row.fit_dg6_mm = get_sensor_value_local(sensorIds, fit.dg_mm, 6);
row.fit_dmu6_mm_per_mm = get_sensor_value_local(sensorIds, fit.dmu_mm_per_mm, 6);
row.noise_std_mv = bundleSet.synthetic_noise_std_mv;
row.point_count_full = numel(bundleSet.X);
row.point_count_fit = numel(fitBundle.X);
row.sg_template_rmse_mean_mv = mean_sg_template_rmse_local(bundleSet, sensorIds);
row.notes = caseSpec.notes;
if nargin < 7 || isempty(prior) %#ok<INUSD>
    return;
end
end

function mu = collect_mu_by_sensor_local(sensorIds, caseSpec)
mu = zeros(1, numel(sensorIds));
if isempty(caseSpec)
    return;
end
for i = 1:numel(sensorIds)
    idx = find(caseSpec.sensor_ids == sensorIds(i), 1);
    if ~isempty(idx)
        mu(i) = caseSpec.mu_cal_by_sensor(idx);
    end
end
end

function m = mean_sg_template_rmse_local(bundle, sensorIds)
vals = nan(1, numel(sensorIds));
for i = 1:numel(sensorIds)
    idx = find([bundle.SuperGaussianKernel.sensorId] == sensorIds(i), 1);
    if ~isempty(idx)
        vals(i) = bundle.SuperGaussianKernel(idx).templateRmseMv;
    end
end
m = mean(vals, 'omitnan');
end

function summary = summarize_physical_rows_local(T)
caseVals = unique(T.case_name, 'stable');
template = struct('case_name', "", 'stage', "", 'sensors', "", ...
    'method', "", 'n', NaN, 'mean_A_id_mm', NaN, ...
    'mean_amp_error_mm', NaN, 'mean_amp_error_percent', NaN, ...
    'mean_weighted_rmse_mv', NaN, 'mean_fit_dg6_mm', NaN, ...
    'mean_fit_dmu6_mm_per_mm', NaN, 'mean_sg_template_rmse_mv', NaN);
summaryRows = repmat(template, 0, 1);
for ic = 1:numel(caseVals)
    idxCase = T.case_name == caseVals(ic);
    sensorVals = unique(T.sensors(idxCase), 'stable');
    for is = 1:numel(sensorVals)
        idxSensor = idxCase & T.sensors == sensorVals(is);
        methodVals = unique(T.method(idxSensor), 'stable');
        for im = 1:numel(methodVals)
            idx = idxSensor & T.method == methodVals(im);
            r = template;
            r.case_name = caseVals(ic);
            r.stage = T.stage(find(idx, 1));
            r.sensors = sensorVals(is);
            r.method = methodVals(im);
            r.n = nnz(idx);
            r.mean_A_id_mm = mean(T.A_id_mm(idx), 'omitnan');
            r.mean_amp_error_mm = mean(T.amp_error_mm(idx), 'omitnan');
            r.mean_amp_error_percent = mean(T.amp_error_percent(idx), 'omitnan');
            r.mean_weighted_rmse_mv = mean(T.weighted_rmse_mv(idx), 'omitnan');
            r.mean_fit_dg6_mm = mean(T.fit_dg6_mm(idx), 'omitnan');
            r.mean_fit_dmu6_mm_per_mm = mean(T.fit_dmu6_mm_per_mm(idx), 'omitnan');
            r.mean_sg_template_rmse_mv = mean(T.sg_template_rmse_mean_mv(idx), 'omitnan');
            summaryRows(end + 1, 1) = r; %#ok<AGROW>
        end
    end
end
summary = struct2table(summaryRows);
end

function Qa = build_quality_report_local(T, cfg)
rows = repmat(struct('check_name', "", 'pass', false, 'value', NaN, ...
    'tolerance', NaN, 'message', ""), 0, 1);
exactMethods = ["fixed_template", "gap_fixed_tilt"];
for im = 1:numel(exactMethods)
    rows(end + 1) = qa_row_local("clean_A_" + exactMethods(im), ...
        max_abs_error_local(T, "stage0_clean", exactMethods(im), "A_id_mm", "A_true_mm"), ...
        cfg.qaCleanAmpTolMm, "<=", ...
        "No-gap clean exact physical methods must recover A."); %#ok<AGROW>
    rows(end + 1) = qa_row_local("zero_A_" + exactMethods(im), ...
        max_abs_value_local(T, "qa_zero_amp_clean", exactMethods(im), "A_id_mm"), ...
        cfg.qaZeroAmpTolMm, "<=", ...
        "Zero-amplitude clean exact physical methods must recover near zero."); %#ok<AGROW>
end
rows(end + 1) = qa_row_local("gap6_clean_dg_gap_fixed_tilt", ...
    max_abs_dg6_error_local(T, "stage1_gap6_clean", "gap_fixed_tilt"), ...
    cfg.qaDgTolMm, "<=", ...
    "Gap-fixed-tilt must estimate injected CH6 dg in the clean gap-only case."); %#ok<AGROW>
rows(end + 1) = qa_row_local("gap6_initial_tilt_clean_dg_gap_fixed_tilt", ...
    max_abs_dg6_error_local(T, "stage2_gap6_initial_tilt_clean", "gap_fixed_tilt"), ...
    cfg.qaDgTolMm, "<=", ...
    "Gap-fixed-tilt must estimate CH6 dg when initial tilt is correctly calibrated."); %#ok<AGROW>
Qa = struct2table(rows);
end

function r = qa_row_local(name, value, tol, relation, message)
r = struct();
r.check_name = string(name);
r.value = value;
r.tolerance = tol;
r.pass = isfinite(value) && value <= tol;
r.message = string(message);
if ~strcmp(string(relation), "<=")
    error('Unsupported QA relation.');
end
end

function v = max_abs_error_local(T, caseName, methodName, valueCol, truthCol)
idx = T.case_name == caseName & T.method == methodName;
if ~any(idx)
    v = NaN;
else
    v = max(abs(T.(valueCol)(idx) - T.(truthCol)(idx)), [], 'omitnan');
end
end

function v = max_abs_value_local(T, caseName, methodName, valueCol)
idx = T.case_name == caseName & T.method == methodName;
if ~any(idx)
    v = NaN;
else
    v = max(abs(T.(valueCol)(idx)), [], 'omitnan');
end
end

function v = max_abs_dg6_error_local(T, caseName, methodName)
idx = T.case_name == caseName & T.method == methodName & T.sensors == "S136";
if ~any(idx)
    v = NaN;
else
    v = max(abs(T.fit_dg6_mm(idx) - T.true_dg6_mm(idx)), [], 'omitnan');
end
end

function plot_step08_summary_local(T, figPath, cfg)
idx = T.sensors == "S136" & T.case_name ~= "qa_zero_amp_clean";
T = T(idx, :);
if isempty(T)
    return;
end
caseVals = unique(T.case_name, 'stable');
methodVals = cfg.physicalMethods;
E = nan(numel(caseVals), numel(methodVals));
for ic = 1:numel(caseVals)
    for im = 1:numel(methodVals)
        idx = T.case_name == caseVals(ic) & T.method == methodVals(im);
        E(ic, im) = mean(T.amp_error_percent(idx), 'omitnan');
    end
end
fig = figure('Name', 'Step08 physical synthetic amplitude error', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 22, 9], ...
    'Visible', 'off');
bar(E);
grid on; box on;
set(gca, 'XTick', 1:numel(caseVals), 'XTickLabel', cellstr(caseVals), ...
    'XTickLabelRotation', 35, 'TickLabelInterpreter', 'none', ...
    'FontName', 'Times New Roman', 'FontSize', 7);
ylabel('Amplitude error (%)');
legend(cellstr(methodVals), 'Location', 'eastoutside', 'Box', 'off', ...
    'Interpreter', 'none');
title('Physical gap/initial-tilt synthetic comparison, S136', ...
    'FontWeight', 'normal', 'Interpreter', 'none');
exportgraphics(fig, figPath, 'Resolution', 250);
close(fig);
end

function mu = collect_mu_by_sensor_from_library_local(bundle, sensorIds)
mu = zeros(1, numel(sensorIds));
for i = 1:numel(sensorIds)
    corr = get_corrected_sensor_local(bundle.CorrectedGapLibrary, sensorIds(i));
    mu(i) = corr.muGapPerXMm;
end
end

function [model, overshoot] = eval_static_warp_template_local( ...
    Tpl, responseSurface, corr, dg, dmu, dtau, xOpr)
xWarpRaw = xOpr(:) - dtau;
xLo = min(Tpl.x_grid(:));
xHi = max(Tpl.x_grid(:));
xLow = min(max(xWarpRaw, xLo), xHi);
overshootLow = max(xLo - xWarpRaw, 0) + max(xWarpRaw - xHi, 0);
vLow = interp1(Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, ...
    xLow, 'pchip', NaN);
if abs(dg) <= eps && abs(dmu) <= eps
    model = vLow;
    overshoot = overshootLow;
    return;
end
[Fdyn, overDyn] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm + dg, ...
    corr.tauMm, corr.xScale, corr.muGapPerXMm + dmu, xWarpRaw);
[Fbase, overBase] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, ...
    corr.tauMm, corr.xScale, corr.muGapPerXMm, xWarpRaw);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
overshoot = max(max(overDyn, overBase), overshootLow);
end

function [model, overshoot] = eval_raw_tilt_path_local( ...
    responseSurface, g0, tau, k, mu, xOpr)
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
B0 = interp_response_coeff_finite_local(responseSurface.xGrid, ...
    responseSurface.coeff(:, 1), xq);
B1 = interp_response_coeff_finite_local(responseSurface.xGrid, ...
    responseSurface.coeff(:, 2), xq);
B2 = interp_response_coeff_finite_local(responseSurface.xGrid, ...
    responseSurface.coeff(:, 3), xq);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function yq = interp_response_coeff_finite_local(xGrid, yGrid, xq)
x = xGrid(:);
y = yGrid(:);
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 2
    yq = NaN(size(xq));
    return;
end
xOk = x(ok);
yOk = y(ok);
xqClip = min(max(xq(:), min(xOk)), max(xOk));
yq = interp1(xOk, yOk, xqClip, 'linear', NaN);
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

function val = get_sensor_value_local(sensorIds, values, sensorId)
idx = find(sensorIds == sensorId, 1);
if isempty(idx)
    val = NaN;
else
    val = values(idx);
end
end

function x = wrap_pi_local(x)
x = mod(x + pi, 2*pi) - pi;
end
