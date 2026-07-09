function seedTable = solve_gradient_displacement_vp_seed(bundle, eoCandidates, cfg, fullEvaluator)
%SOLVE_GRADIENT_DISPLACEMENT_VP_SEED Shared displacement-domain VP scan.
%
% The seed observation is q_obs = -(V - F0)./Fx - eta_s. VP only provides
% EO candidates and initial values; when fullEvaluator is supplied, ranking
% uses the complete waveform model weighted RMSE.

if nargin < 4
    fullEvaluator = [];
end

eoCandidates = unique(round(eoCandidates(:).'));

V = get_vector_field_local(bundle, {'V', 'v'});
F0 = get_vector_field_local(bundle, {'F0', 'template_v'});
Fx = get_vector_field_local(bundle, {'Fx', 'template_dv_dx'});
Theta = get_vector_field_local(bundle, {'Theta', 'theta'});
W = get_vector_field_local(bundle, {'W', 'fit_weight'});
sensorIds = get_vector_field_local(bundle, {'sensorIds', 'sensor_ids'}, true);
sensorKey = get_vector_field_local(bundle, {'sensorIndex', 'sensor_index', 'sensor_id'});

pointCount = get_scalar_field_local(bundle, {'pointCount', 'point_count'}, numel(V));
nSensor = numel(sensorIds);
etaStatic = bttcore.fixed_sensor_eta_by_position(sensorIds, cfg);
etaVec0 = bttcore.sensor_eta_vector(sensorKey, sensorIds, etaStatic);

rows = repmat(struct( ...
    'EO', NaN, ...
    'A', NaN, ...
    'phi', NaN, ...
    'dx', NaN, ...
    'linearRmseMv', inf, ...
    'linearRmseMm', inf, ...
    'fallbackRmseMv', inf, ...
    'pointCount', pointCount, ...
    'validPointCount', 0, ...
    'gradientRef', NaN, ...
    'fixedSensorEtaMm', zeros(1, nSensor)), numel(eoCandidates), 1);

if isempty(V) || isempty(F0) || isempty(Fx) || isempty(Theta) || isempty(W)
    seedTable = struct2table(rows);
    return;
end

V = V(:);
F0 = F0(:);
Fx = Fx(:);
Theta = Theta(:);
W = W(:);
etaVec0 = etaVec0(:);

n = min([numel(V), numel(F0), numel(Fx), numel(Theta), numel(W), numel(etaVec0)]);
if n < 5
    seedTable = struct2table(rows);
    return;
end
V = V(1:n);
F0 = F0(1:n);
Fx = Fx(1:n);
Theta = Theta(1:n);
W = W(1:n);
etaVec0 = etaVec0(1:n);

gradAbs = abs(Fx);
gradRef = prctile(gradAbs(isfinite(gradAbs)), ...
    get_cfg_scalar_local(cfg, {'vpGradientReferenceQuantile', ...
    'vp_gradient_reference_quantile'}, 95));
if ~isfinite(gradRef) || gradRef <= 0
    gradRef = max(gradAbs(isfinite(gradAbs)), [], 'omitnan');
end
if ~isfinite(gradRef) || gradRef <= 0
    seedTable = struct2table(rows);
    return;
end

gradMinRatio = get_cfg_scalar_local(cfg, {'vpGradientMinRatio', ...
    'vp_gradient_min_ratio'}, 0.10);
valid = isfinite(V) & isfinite(F0) & isfinite(Fx) & isfinite(Theta) & ...
    isfinite(W) & isfinite(etaVec0) & ...
    gradAbs >= max(max(0, gradMinRatio) * gradRef, eps);

minFitPoints = max(5, floor(get_cfg_scalar_local(cfg, ...
    {'minFitPoints', 'min_fit_points'}, 5)));
if nnz(valid) < minFitPoints
    seedTable = struct2table(rows);
    return;
end

V = V(valid);
F0 = F0(valid);
Fx = Fx(valid);
Theta = Theta(valid);
W = W(valid);
etaVec0 = etaVec0(valid);

qObs = -(V - F0) ./ Fx - etaVec0;
weightFloor = get_cfg_scalar_local(cfg, {'weightFloor', 'weight_floor'}, 0.05);
gradPower = get_cfg_scalar_local(cfg, {'vpGradientWeightPower', ...
    'vp_gradient_weight_power'}, 2);
wq = max(W, weightFloor) .* min(abs(Fx) ./ max(gradRef, eps), 1) .^ ...
    max(0, gradPower);
wq(~isfinite(wq)) = weightFloor;
wq = max(wq, weightFloor);

for i = 1:numel(eoCandidates)
    eo = eoCandidates(i);
    s1 = sin(eo .* Theta);
    c1 = cos(eo .* Theta);
    X = [ones(size(qObs)), s1, c1];
    good = all(isfinite(X), 2) & isfinite(qObs) & isfinite(wq);
    if nnz(good) < max(minFitPoints, size(X, 2) + 1)
        continue;
    end

    Xg = X(good, :);
    yg = qObs(good);
    wg = sqrt(wq(good));
    Xw = Xg .* wg;
    yw = yg .* wg;

    try
        beta = Xw \ yw;
    catch
        beta = pinv(Xw) * yw;
    end

    if isempty(beta) || numel(beta) < 3 || any(~isfinite(beta))
        continue;
    end

    dx = beta(1);
    a = beta(2);
    b = beta(3);
    A = hypot(a, b);
    phi = wrap_to_pi_local(atan2(b, a));

    A = clamp_abs_local(A, get_cfg_scalar_local(cfg, ...
        {'amplitudeLimitMm', 'amplitude_limit_mm'}, inf));
    dx = clamp_scalar_local(dx, get_cfg_scalar_local(cfg, ...
        {'dxLimitMm', 'dx_c_limit_mm', 'dx_c_limit'}, inf));

    linRes = yg - Xg * beta;
    linearRmseMm = sqrt(sum(wq(good) .* linRes.^2) / max(sum(wq(good)), eps));
    fallbackRmse = linearRmseMm;

    if ~isempty(fullEvaluator)
        try
            detail = fullEvaluator([A, phi, dx], eo, bundle, cfg);
            if isstruct(detail) && isfield(detail, 'weightedRmseMv') && ...
                    isfinite(detail.weightedRmseMv)
                fallbackRmse = detail.weightedRmseMv;
            end
        catch
            fallbackRmse = linearRmseMm;
        end
    end

    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx = dx;
    rows(i).linearRmseMv = linearRmseMm;
    rows(i).linearRmseMm = linearRmseMm;
    rows(i).fallbackRmseMv = fallbackRmse;
    rows(i).pointCount = pointCount;
    rows(i).validPointCount = nnz(good);
    rows(i).gradientRef = gradRef;
    rows(i).fixedSensorEtaMm = etaStatic(:).';
end

[~, order] = sortrows([[rows.fallbackRmseMv].', ...
    [rows.linearRmseMm].', [rows.EO].'], [1, 2, 3]);
seedTable = struct2table(rows(order));
end

function value = get_vector_field_local(s, names, allowEmpty)
if nargin < 3
    allowEmpty = false;
end
value = [];
for i = 1:numel(names)
    if isfield(s, names{i})
        value = double(s.(names{i})(:));
        return;
    end
end
if ~allowEmpty
    value = [];
end
end

function value = get_scalar_field_local(s, names, defaultValue)
value = defaultValue;
for i = 1:numel(names)
    if isfield(s, names{i}) && isscalar(s.(names{i})) && isfinite(s.(names{i}))
        value = double(s.(names{i}));
        return;
    end
end
end

function value = get_cfg_scalar_local(cfg, names, defaultValue)
value = defaultValue;
if nargin < 1 || isempty(cfg) || ~isstruct(cfg)
    return;
end
for i = 1:numel(names)
    if isfield(cfg, names{i}) && isscalar(cfg.(names{i})) && ...
            isnumeric(cfg.(names{i})) && isfinite(cfg.(names{i}))
        value = double(cfg.(names{i}));
        return;
    end
end
end

function y = clamp_abs_local(x, limit)
if ~isfinite(x)
    x = 0;
end
if isfinite(limit)
    y = min(max(abs(x), 0), abs(limit));
else
    y = abs(x);
end
end

function y = clamp_scalar_local(x, limit)
if ~isfinite(x)
    x = 0;
end
if isfinite(limit)
    y = max(min(x, abs(limit)), -abs(limit));
else
    y = x;
end
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2 * pi) - pi;
end
