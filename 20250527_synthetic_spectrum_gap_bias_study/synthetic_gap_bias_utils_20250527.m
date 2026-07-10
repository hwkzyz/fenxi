function varargout = synthetic_gap_bias_utils_20250527(action, varargin)
%synthetic_gap_bias_utils_20250527 Shared helpers for the synthetic study.

action = lower(strtrim(char(action)));
switch action
    case 'filter_bundle_by_sensors'
        varargout{1} = filter_bundle_by_sensors_local(varargin{:});
    case 'thin_bundle'
        varargout{1} = thin_bundle_local(varargin{:});
    case 'make_synthetic_voltage_bundle'
        varargout{1} = make_synthetic_voltage_bundle_local(varargin{:});
    case 'fit_nested_sg_model'
        varargout{1} = fit_nested_sg_model_local(varargin{:});
    case 'build_jacobian_columns'
        varargout{1} = build_jacobian_columns_local(varargin{:});
    case 'eval_super_gaussian'
        varargout{1} = eval_super_gaussian_local(varargin{:});
    case 'weighted_linear_fit'
        [varargout{1:nargout}] = weighted_linear_fit_local(varargin{:});
    case 'projection_ratio'
        varargout{1} = projection_ratio_local(varargin{:});
    case 'weighted_norm'
        varargout{1} = weighted_norm_local(varargin{:});
    case 'weighted_abs_corr'
        varargout{1} = weighted_abs_corr_local(varargin{:});
    case 'residualize_columns'
        varargout{1} = residualize_columns_local(varargin{:});
    case 'get_sensor_value'
        varargout{1} = get_sensor_value_local(varargin{:});
    otherwise
        error('Unknown synthetic_gap_bias_utils_20250527 action: %s', action);
end
end

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

function fit = fit_nested_sg_model_local(fitBundle, evalBundle, methodName, cfg, prior, initialParams)
if nargin < 6
    initialParams = [];
end

sensorIds = unique(fitBundle.S(:)).';
freeSensorIds = sensorIds(2:end);
nFree = numel(freeSensorIds);
methodName = string(methodName);

initialList = build_initial_list_local(initialParams, nFree, methodName, cfg, prior);
best = struct('obj', inf, 'p', []);

opts = optimset('Display', 'off', 'MaxIter', cfg.fitMaxIter, ...
    'MaxFunEvals', cfg.fitMaxFunEvals, 'TolX', 1e-6, 'TolFun', 1e-8);

for i0 = 1:numel(initialList)
    p0 = initialList{i0};
    objFun = @(p) nested_objective_local(p, fitBundle, methodName, ...
        sensorIds, freeSensorIds, cfg);
    pOpt = fminsearch(objFun, p0, opts);
    obj = objFun(pOpt);
    if obj < best.obj
        best.obj = obj;
        best.p = pOpt;
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

function initialList = build_initial_list_local(initialParams, nFree, methodName, cfg, prior)
initialList = {};
if isempty(initialParams)
    seedAmp = unique([cfg.fitSeedAmpMm(:); prior.main_amp_tip_mm]);
    seedPhi = cfg.fitSeedPhiRad(:);
    for ia = 1:numel(seedAmp)
        for ip = 1:numel(seedPhi)
            initialList{end + 1} = build_initial_parameter_local( ... %#ok<AGROW>
                seedAmp(ia), seedPhi(ip), nFree, methodName);
        end
    end
    return;
end

if iscell(initialParams)
    initialList = initialParams(:).';
else
    for i = 1:size(initialParams, 1)
        initialList{end + 1} = initialParams(i, :); %#ok<AGROW>
    end
end
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

function D = build_jacobian_columns_local(bundle, cfg, prior, baseState)
if nargin < 4 || isempty(baseState)
    baseState = struct();
end

theta = bundle.Theta(:);
x = bundle.X(:);
sensorIds = unique(bundle.S(:)).';
nSensor = numel(sensorIds);

A = get_field_default_local(baseState, 'A_mm', prior.main_amp_tip_mm);
phi = get_field_default_local(baseState, 'phi_rad', cfg.truePhaseRad);
d0 = get_field_default_local(baseState, 'd0_mm', cfg.trueD0Mm);
dgBySensor = get_vector_field_default_local(baseState, 'dg_mm', nSensor);
dmuBySensor = get_vector_field_default_local(baseState, 'dmu_mm_per_mm', nSensor);

sineCol = sin(cfg.targetEO .* theta + phi);
cosCol = cos(cfg.targetEO .* theta + phi);
u = A .* sineCol + d0;

dgSample = zeros(size(x));
dmuSample = zeros(size(x));
for i = 1:nSensor
    mask = bundle.S(:) == sensorIds(i);
    dgSample(mask) = dgBySensor(i);
    dmuSample(mask) = dmuBySensor(i);
end
staticShift = dgSample + dmuSample .* (x - bundle.xc_pt(:));
xIn = x - u - staticShift;
dVdx = super_gaussian_derivative_local(xIn, bundle);

D = struct();
D.sensorIds = sensorIds;
D.weight = max(bundle.W(:), 0);
if ~any(D.weight > 0)
    D.weight = ones(size(x));
end
D.VBase = eval_super_gaussian_local(xIn, bundle);
D.xIn = xIn;
D.sineCol = sineCol;
D.cosCol = cosCol;
D.J_A = dVdx .* (-sineCol);
D.J_phi = dVdx .* (-(A .* cosCol));
D.J_d0 = dVdx .* (-ones(size(x)));
D.J_dg = zeros(numel(x), nSensor);
D.J_dmu = zeros(numel(x), nSensor);

for i = 1:nSensor
    mask = bundle.S(:) == sensorIds(i);
    D.J_dg(mask, i) = dVdx(mask) .* (-1);
    D.J_dmu(mask, i) = dVdx(mask) .* (-(x(mask) - bundle.xc_pt(mask)));
end
end

function val = get_field_default_local(s, name, defaultVal)
if isfield(s, name) && isfinite(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function val = get_vector_field_default_local(s, name, n)
val = zeros(1, n);
if isfield(s, name) && numel(s.(name)) >= n
    raw = s.(name);
    val = raw(1:n);
    val = val(:).';
end
end

function [beta, detail] = weighted_linear_fit_local(X, y, weight)
y = y(:);
weight = max(weight(:), 0);
valid = isfinite(y) & isfinite(weight) & all(isfinite(X), 2);
Xv = X(valid, :);
yv = y(valid);
wv = weight(valid);
if isempty(Xv)
    beta = nan(size(X, 2), 1);
    detail = struct('rmse', NaN, 'weightedRmse', NaN, 'rank', 0, 'cond', NaN, ...
        'residual', nan(size(y)));
    return;
end

sw = sqrt(wv);
Xw = Xv .* sw;
yw = yv .* sw;
badCols = false(1, size(Xw, 2));
for i = 1:size(Xw, 2)
    badCols(i) = norm(Xw(:, i)) <= eps || any(~isfinite(Xw(:, i)));
end
Xkeep = Xw(:, ~badCols);
betaKeep = Xkeep \ yw;
beta = zeros(size(X, 2), 1);
beta(~badCols) = betaKeep;

resValid = yv - Xv * beta;
res = nan(size(y));
res(valid) = resValid;
s = svd(Xkeep, 'econ');
if isempty(s)
    r = 0;
    c = NaN;
elseif numel(s) < 2 || min(s) <= eps
    r = sum(s > max(size(Xkeep)) * eps(max(s)));
    c = Inf;
else
    r = sum(s > max(size(Xkeep)) * eps(max(s)));
    c = max(s) ./ min(s);
end

detail = struct();
detail.rmse = sqrt(mean(resValid.^2, 'omitnan'));
detail.weightedRmse = sqrt(sum(wv .* resValid.^2) ./ max(sum(wv), eps));
detail.rank = r;
detail.cond = c;
detail.residual = res;
end

function R = residualize_columns_local(C, basis, weight)
if isempty(C)
    R = C;
    return;
end
R = nan(size(C));
for i = 1:size(C, 2)
    R(:, i) = residualize_one_column_local(C(:, i), basis, weight);
end
end

function r = residualize_one_column_local(col, basis, weight)
col = col(:);
weight = max(weight(:), 0);
valid = isfinite(col) & all(isfinite(basis), 2) & isfinite(weight);
r = zeros(size(col));
if ~any(valid)
    r(:) = NaN;
    return;
end

cw = col(valid) .* sqrt(weight(valid));
Bw = basis(valid, :) .* sqrt(weight(valid));
Bw = remove_bad_columns_local(Bw);
if isempty(Bw)
    r(valid) = col(valid);
    return;
end
[Q, ~] = qr(Bw, 0);
proj = Q * (Q' * cw);
rw = cw - proj;
r(valid) = rw ./ max(sqrt(weight(valid)), eps);
end

function ratio = projection_ratio_local(col, basis, weight)
col = col(:);
weight = max(weight(:), 0);
valid = isfinite(col) & all(isfinite(basis), 2) & isfinite(weight);
col = col(valid);
basis = basis(valid, :);
weight = weight(valid);
if isempty(col) || weighted_norm_local(col, weight) <= eps || isempty(basis)
    ratio = NaN;
    return;
end
sw = sqrt(weight);
cw = col .* sw;
Bw = basis .* sw;
Bw = remove_bad_columns_local(Bw);
if isempty(Bw)
    ratio = NaN;
    return;
end
[Q, ~] = qr(Bw, 0);
proj = Q * (Q' * cw);
ratio = norm(proj) ./ max(norm(cw), eps);
ratio = min(max(ratio, 0), 1);
end

function rho = weighted_abs_corr_local(a, b, weight)
a = a(:);
b = b(:);
weight = max(weight(:), 0);
valid = isfinite(a) & isfinite(b) & isfinite(weight);
a = a(valid);
b = b(valid);
weight = weight(valid);
if isempty(a)
    rho = NaN;
    return;
end
sw = sqrt(weight);
aw = a .* sw;
bw = b .* sw;
rho = abs(dot(aw, bw)) ./ max(norm(aw) .* norm(bw), eps);
end

function nrm = weighted_norm_local(a, weight)
a = a(:);
weight = max(weight(:), 0);
valid = isfinite(a) & isfinite(weight);
nrm = norm(a(valid) .* sqrt(weight(valid)));
end

function B = remove_bad_columns_local(B)
if isempty(B)
    return;
end
keep = false(1, size(B, 2));
for i = 1:size(B, 2)
    keep(i) = all(isfinite(B(:, i))) && norm(B(:, i)) > eps;
end
B = B(:, keep);
end

function v = eval_super_gaussian_local(xIn, bundle)
v = bundle.B_pt(:) .* exp(-abs((xIn(:) - bundle.xc_pt(:)) ./ bundle.w_pt(:)).^bundle.n_pt(:)) + ...
    bundle.base_pt(:);
end

function dVdx = super_gaussian_derivative_local(xIn, bundle)
z = (xIn(:) - bundle.xc_pt(:)) ./ bundle.w_pt(:);
az = abs(z);
sz = sign(z);
n = bundle.n_pt(:);
baseShape = exp(-(az.^n));
localSlope = -n .* (max(az, eps).^(n - 1)) .* sz ./ bundle.w_pt(:);
dVdx = bundle.B_pt(:) .* baseShape .* localSlope;
dVdx(~isfinite(dVdx)) = 0;
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

function x = wrap_pi_local(x)
x = mod(x + pi, 2*pi) - pi;
end
