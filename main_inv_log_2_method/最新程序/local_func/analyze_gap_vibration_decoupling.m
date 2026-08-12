function diag = analyze_gap_vibration_decoupling(highMap, templateLib, result, opts)
%analyze_gap_vibration_decoupling  Local gap-vibration decoupling diagnostic.
%
% This diagnostic checks whether the gap sensitivity direction can be
% represented by the nuisance subspace spanned by spatial alignment and the
% identified harmonic displacement basis. It does not prove uniqueness of the
% vibration parameters.

if nargin < 4 || isempty(opts)
    opts = struct();
end

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);

g = get_result_field(result, {'g_used', 'g'}, NaN);
dx = get_result_field(result, {'dx_used', 'dx'}, NaN);
if ~isfinite(g) || ~isfinite(dx)
    error('Result must contain finite gap and dx estimates.');
end

p = get_result_p(result);
if isempty(p)
    error('Result must contain harmonic parameter vector p.');
end

freqMode = string(get_opt(opts, 'freqMode', "identified"));
switch freqMode
    case "single"
        p = p(1:3);
    case {"identified", "all"}
        % Use all frequency components available in p.
    otherwise
        error('Unknown freqMode: %s', freqMode);
end

u = local_fit_u(p, t);
xq = x - dx - u;

gapStep = get_opt(opts, 'gapStep', []);
if isempty(gapStep)
    gapStep = max(1e-4, 1e-3 * max(abs(g), 1));
end
gMin = min(templateLib.gapTrain(:));
gMax = max(templateLib.gapTrain(:));
gapStep = min(gapStep, 0.45 * max(g - gMin, eps));
gapStep = min(gapStep, 0.45 * max(gMax - g, eps));
gapStep = max(gapStep, 1e-6);

jg = (eval_gap_template(templateLib, g + gapStep, xq) - ...
    eval_gap_template(templateLib, g - gapStep, xq)) ./ (2 * gapStep);
jg = jg(:);
Fx = eval_gap_derivative(templateLib, g, xq);
Fx = Fx(:);

Ju = build_nuisance_basis(Fx, t, p, opts);
[Q, rankJu, svJu] = orth_basis(Ju);
projJg = Q * (Q' * jg);
resJg = jg - projJg;

Jfull = [normalize_col(jg), normalize_columns(Ju)];
svFull = svd(Jfull, 'econ');

diag = struct();
diag.gap_mm = g;
diag.dx_mm = dx;
diag.freq_Hz = p(3:3:end).';
diag.num_frequency_components = numel(p) / 3;
diag.sample_count = numel(t);
diag.nuisance_column_count = size(Ju, 2);
diag.rank_Ju = rankJu;
diag.rank_Jfull_normalized = sum(svFull > max(size(Jfull)) * eps(max(svFull)));
diag.jg_norm = norm(jg);
diag.jg_projected_norm = norm(projJg);
diag.jg_residual_norm = norm(resJg);
diag.rho_gap_in_nuisance = norm(projJg) / max(norm(jg), eps);
diag.gamma_gap_residual = norm(resJg);
diag.gamma_gap_relative = norm(resJg) / max(norm(jg), eps);
diag.sigma_Ju = svJu(:).';
diag.sigma_Jfull_normalized = svFull(:).';
diag.cond_Ju = local_cond(svJu);
diag.cond_Jfull_normalized = local_cond(svFull);
diag.gap_step_mm = gapStep;
diag.rmse_wave = sqrt(mean((V - eval_gap_template(templateLib, g, xq)).^2));
diag.interpretation = interpret_decoupling(diag);
end

function Ju = build_nuisance_basis(Fx, t, p, opts)
includeDx = get_opt(opts, 'includeDx', true);
includeHarmonics = get_opt(opts, 'includeHarmonics', true);
includeFreqColumns = get_opt(opts, 'includeFreqColumns', false);

cols = {};
if includeDx
    cols{end+1} = -Fx; %#ok<AGROW>
end

if includeHarmonics
    for k = 1:3:numel(p)
        A = p(k);
        phi = p(k + 1);
        f = p(k + 2);
        s = sin(2*pi*f*t);
        c = cos(2*pi*f*t);
        cols{end+1} = -Fx .* s; %#ok<AGROW>
        cols{end+1} = -Fx .* c; %#ok<AGROW>
        if includeFreqColumns
            duDf = A .* cos(2*pi*f*t + phi) .* (2*pi*t);
            cols{end+1} = -Fx .* duDf; %#ok<AGROW>
        end
    end
end

Ju = zeros(numel(t), numel(cols));
for i = 1:numel(cols)
    Ju(:, i) = cols{i}(:);
end
Ju = normalize_columns(Ju);
end

function u = local_fit_u(p, t)
u = zeros(size(t));
for k = 1:3:numel(p)
    u = u + p(k) .* sin(2*pi*p(k + 2).*t + p(k + 1));
end
end

function [Q, r, s] = orth_basis(A)
A = normalize_columns(A);
[Q0, R] = qr(A, 0);
d = abs(diag(R));
tol = max(size(A)) * eps(max(d, [], 'omitnan'));
r = sum(d > tol);
if r == 0
    Q = zeros(size(A, 1), 0);
else
    Q = Q0(:, 1:r);
end
s = svd(A, 'econ');
end

function A = normalize_columns(A)
for i = 1:size(A, 2)
    ni = norm(A(:, i));
    if ni > 0
        A(:, i) = A(:, i) ./ ni;
    end
end
end

function a = normalize_col(a)
n = norm(a);
if n > 0
    a = a ./ n;
end
end

function c = local_cond(s)
s = s(:);
s = s(isfinite(s) & s > 0);
if isempty(s)
    c = Inf;
else
    c = max(s) / max(min(s), eps);
end
end

function txt = interpret_decoupling(diag)
if diag.gamma_gap_relative < 0.05
    txt = "gap direction is almost contained in nuisance subspace";
elseif diag.gamma_gap_relative < 0.20
    txt = "gap direction has weak residual component outside nuisance subspace";
else
    txt = "gap direction has a clear residual component outside nuisance subspace";
end
end

function p = get_result_p(result)
p = [];
if isfield(result, 'fit') && isstruct(result.fit) && isfield(result.fit, 'p')
    p = result.fit.p(:).';
elseif isfield(result, 'p')
    p = result.p(:).';
end
end

function val = get_result_field(S, names, defaultVal)
val = defaultVal;
for i = 1:numel(names)
    if isfield(S, names{i}) && ~isempty(S.(names{i}))
        val = S.(names{i});
        return;
    end
end
end

function val = get_opt(S, name, defaultVal)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end
