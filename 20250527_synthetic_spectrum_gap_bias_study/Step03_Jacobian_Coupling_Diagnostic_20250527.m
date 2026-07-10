%% Step03: local Jacobian coupling between vibration and static tilt columns.

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

rows = table();
for iSet = 1:numel(cfg.referenceSensorSets)
    sensorSet = cfg.referenceSensorSets{iSet};
    sensorTag = cfg.referenceSensorTags(iSet);
    bundle = filter_bundle_by_sensors_local(baseBundle, sensorSet);
    bundle = thin_bundle_local(bundle, cfg.maxFitPoints);

    D = build_jacobian_columns_local(bundle, cfg, prior);
    vibBasis = [D.J_A, D.J_phi, D.J_d0];
    vibBasis = remove_bad_columns_local(vibBasis);
    [rankVib, condVib] = rank_cond_local(vibBasis);

    sensorIds = unique(bundle.S(:)).';
    for iSensor = 1:numel(sensorIds)
        sid = sensorIds(iSensor);
        Jdg = D.J_dg(:, iSensor);
        Jdmu = D.J_dmu(:, iSensor);

        row = table();
        row.sensors = sensorTag;
        row.sensor_id = sid;
        row.point_count = numel(bundle.T);
        row.A_true_mm = prior.main_amp_tip_mm;
        row.rank_vib_basis = rankVib;
        row.cond_vib_basis = condVib;
        row.rho_A_dg = weighted_abs_corr_local(D.J_A, Jdg, D.weight);
        row.rho_A_dmu = weighted_abs_corr_local(D.J_A, Jdmu, D.weight);
        row.proj_dg_on_vib = projection_ratio_local(Jdg, vibBasis, D.weight);
        row.proj_dmu_on_vib = projection_ratio_local(Jdmu, vibBasis, D.weight);
        row.eta_dg_outside_vib = sqrt(max(0, 1 - row.proj_dg_on_vib.^2));
        row.eta_dmu_outside_vib = sqrt(max(0, 1 - row.proj_dmu_on_vib.^2));
        row.norm_J_A = weighted_norm_local(D.J_A, D.weight);
        row.norm_J_dg = weighted_norm_local(Jdg, D.weight);
        row.norm_J_dmu = weighted_norm_local(Jdmu, D.weight);
        rows = [rows; row]; %#ok<AGROW>
    end
end

csvPath = fullfile(cfg.outputDir, 'Step03_jacobian_coupling.csv');
writetable(rows, csvPath);

disp('Step03 Jacobian coupling:');
disp(rows);
fprintf('Step03 completed.\n  %s\n', csvPath);

function D = build_jacobian_columns_local(bundle, cfg, prior)
theta = bundle.Theta(:);
x = bundle.X(:);
sensorIds = unique(bundle.S(:)).';
nSensor = numel(sensorIds);

A = prior.main_amp_tip_mm;
phi = cfg.truePhaseRad;
d0 = cfg.trueD0Mm;

sineCol = sin(cfg.targetEO .* theta + phi);
cosCol = cos(cfg.targetEO .* theta + phi);
u = A .* sineCol + d0;
xIn = x - u;
dVdx = super_gaussian_derivative_local(xIn, bundle);

D = struct();
D.weight = max(bundle.W(:), 0);
if ~any(D.weight > 0)
    D.weight = ones(size(x));
end
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

function [r, c] = rank_cond_local(B)
B = remove_bad_columns_local(B);
if isempty(B)
    r = 0;
    c = NaN;
    return;
end
s = svd(B, 'econ');
r = sum(s > max(size(B)) * eps(max(s)));
if numel(s) < 2 || min(s) <= eps
    c = Inf;
else
    c = max(s) ./ min(s);
end
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
    'xc_pt','base_pt'};
for i = 1:numel(sampleFields)
    f = sampleFields{i};
    if isfield(bundleIn, f) && numel(bundleIn.(f)) == numel(mask)
        bundleOut.(f) = bundleIn.(f)(mask);
    end
end
bundleOut.point_count = nnz(mask);
end
