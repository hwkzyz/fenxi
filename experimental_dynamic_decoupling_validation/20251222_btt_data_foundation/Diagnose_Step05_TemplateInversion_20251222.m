function Summary = Diagnose_Step05_TemplateInversion_20251222()
%DIAGNOSE_STEP05_TEMPLATEINVERSION_20251222
% Read-only diagnostic for the 20251222 Step05 route.
%
% The super-Gaussian route first maps voltage to a static x estimate and
% then fits the residual displacement. This script applies the same idea to
% the non-parametric Step04 template:
%
%   x_static = T_{s,b}^{-1}(V_obs) on the same side of the pulse center
%   y_obs    = x - x_static
%
% It compares EO ranking for two linear residual models:
%
%   common_offset:        y_obs = A*sin(EO*theta + phi) + d0
%   sensor_static_offset: y_obs = A*sin(EO*theta + phi) + d0 + eta_s
%
% No files used by the formal Step04/Step05 chain are modified.

cfg = BTTProjectConfig_20251222();

result_file = strtrim(getenv('STEP05_DIAG_RESULT_FILE'));
if isempty(result_file)
    result_file = fullfile(cfg.step05_direct_template_output_dir, ...
        cfg.dynamic_cases{1}, ...
        'Result_Step05_FoundationMainPulseAdaptiveNoEta_B1_S123_20251222.mat');
end
if ~isfile(result_file)
    error('Missing Step05 result file for diagnostic: %s', result_file);
end

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_template_inversion_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

loaded = load(result_file, 'Result');
Result = loaded.Result;
S = Result.AnalysisSettings;
if ~isfield(S, 'freq_search_hz') || isempty(S.freq_search_hz)
    S.freq_search_hz = [100 1000];
end
if ~isfield(S, 'eo_pad') || isempty(S.eo_pad)
    S.eo_pad = 2;
end
if ~isfield(S, 'analysis_sensors') || isempty(S.analysis_sensors)
    S.analysis_sensors = Result.SensorIDs;
end
if ~isfield(S, 'weight_floor') || isempty(S.weight_floor)
    S.weight_floor = 0.05;
end
    if ~isfield(S, 'template_interp_method') || isempty(S.template_interp_method)
        S.template_interp_method = 'pchip';
    end
    S.sg_static_calibration = load_sg_static_calibration_local();

    all_rows = table();
    summary_rows = repmat(make_empty_summary_row_local(), ...
        numel(Result.WindowResult), 1);

for iw = 1:numel(Result.WindowResult)
    WR = Result.WindowResult(iw);
    B = choose_diagnostic_bundle_local(WR);
    if isempty(B) || ~isstruct(B) || ~isfield(B, 'x') || isempty(B.x)
        summary_rows(iw).window_id = iw;
        summary_rows(iw).status = "missing_bundle";
        continue;
    end

    eo_candidates = build_eo_candidates_local( ...
        B.rot_freq_mean_hz, S.freq_search_hz, S.eo_pad);

    TinvertCommon = scan_template_inversion_eo_local( ...
        B, eo_candidates, S, false);
    TinvertEta = scan_template_inversion_eo_local( ...
        B, eo_candidates, S, true);
    TSgLike = scan_sg_like_eo_local(B, eo_candidates, S);
    TStep04Sg = scan_step04_parametric_sg_eo_local(B, eo_candidates, S);

    TinvertCommon.model = repmat("template_inversion_common_offset", ...
        height(TinvertCommon), 1);
    TinvertEta.model = repmat("template_inversion_sensor_static_offset", ...
        height(TinvertEta), 1);
    TSgLike.model = repmat("sg_static_model_replay", height(TSgLike), 1);
    TStep04Sg.model = repmat("step04_parametric_sg_replay", ...
        height(TStep04Sg), 1);
    T = [TinvertCommon; TinvertEta; TSgLike; TStep04Sg];
    T.window_id = repmat(iw, height(T), 1);
    T = movevars(T, {'window_id', 'model'}, 'Before', 1);
    all_rows = [all_rows; T]; %#ok<AGROW>

    summary_rows(iw) = summarize_window_local(iw, WR, T);
end

CandidateTable = all_rows;
WindowSummary = struct2table(summary_rows);

candidate_csv = fullfile(out_dir, ...
    'Step05_TemplateInversion_Candidates_B1_S123_20251222.csv');
summary_csv = fullfile(out_dir, ...
    'Step05_TemplateInversion_WindowSummary_B1_S123_20251222.csv');
writetable(CandidateTable, candidate_csv);
writetable(WindowSummary, summary_csv);

Summary = struct();
Summary.ResultFile = result_file;
Summary.OutputDir = out_dir;
Summary.CandidateTable = CandidateTable;
Summary.WindowSummary = WindowSummary;

fprintf('\n=== Step05 template-inversion diagnostic (20251222) ===\n');
fprintf('Source Step05 result:\n  %s\n', result_file);
fprintf('Output dir:\n  %s\n\n', out_dir);
disp(WindowSummary);
print_findings_local(WindowSummary);
end


function B = choose_diagnostic_bundle_local(WR)
B = [];
if isfield(WR, 'CoreBundlePreview') && ~isempty(WR.CoreBundlePreview)
    B = WR.CoreBundlePreview;
elseif isfield(WR, 'BundlePreview') && ~isempty(WR.BundlePreview)
    B = WR.BundlePreview;
elseif isfield(WR, 'CandidateBundlePreview') && ~isempty(WR.CandidateBundlePreview)
    B = WR.CandidateBundlePreview;
end
end


function T = scan_template_inversion_eo_local(B, eo_candidates, S, include_sensor_eta)
x = B.x(:);
v = B.v(:);
theta = B.theta(:);
sensor_id = B.sensor_id(:);
w = B.fit_weight(:);
w(~isfinite(w)) = S.weight_floor;
w = max(w, S.weight_floor);

x_static = invert_template_same_side_local(B);
y_obs = x - x_static;

base_valid = isfinite(x) & isfinite(v) & isfinite(theta) & ...
    isfinite(w) & isfinite(y_obs);

sensors = double(S.analysis_sensors(:).');
n_sensor = numel(sensors);

rows = repmat(struct( ...
    'EO', NaN, ...
    'rank', NaN, ...
    'A_mm', NaN, ...
    'phi_rad', NaN, ...
    'd0_mm', NaN, ...
    'eta_CH1_mm', 0, ...
    'eta_CH2_mm', NaN, ...
    'eta_CH3_mm', NaN, ...
    'weighted_y_rmse_mm', inf, ...
    'weighted_voltage_rmse_v', inf, ...
    'plain_voltage_rmse_v', inf, ...
    'point_count', nnz(base_valid), ...
    'valid_query_fraction', NaN), numel(eo_candidates), 1);

for i = 1:numel(eo_candidates)
    EO = eo_candidates(i);
    s = sin(EO .* theta);
    c = cos(EO .* theta);
    X = [s, c, ones(size(s))];

    if include_sensor_eta
        for is = 2:n_sensor
            X = [X, double(sensor_id == sensors(is))]; %#ok<AGROW>
        end
    end

    valid = base_valid & all(isfinite(X), 2);
    if nnz(valid) < max(20, size(X, 2) + 2)
        continue;
    end

    Xg = X(valid, :);
    yg = y_obs(valid);
    wg = sqrt(w(valid));

    beta = (Xg .* wg) \ (yg .* wg);
    if isempty(beta) || any(~isfinite(beta))
        continue;
    end

    y_fit = X * beta;
    [wrmse_v, prmse_v, valid_q_frac] = evaluate_voltage_rmse_local( ...
        B, y_fit, w, S);

    y_res = yg - Xg * beta;
    y_rmse = sqrt(sum(w(valid) .* y_res.^2) ./ max(sum(w(valid)), eps));

    eta = zeros(1, n_sensor);
    if include_sensor_eta && n_sensor > 1 && numel(beta) >= 3 + n_sensor - 1
        eta(2:end) = beta(4:(3 + n_sensor - 1)).';
    end

    rows(i).EO = EO;
    rows(i).A_mm = hypot(beta(1), beta(2));
    rows(i).phi_rad = wrap_to_pi_local(atan2(beta(2), beta(1)));
    rows(i).d0_mm = beta(3);
    rows(i).eta_CH1_mm = eta_value_local(eta, sensors, 1);
    rows(i).eta_CH2_mm = eta_value_local(eta, sensors, 2);
    rows(i).eta_CH3_mm = eta_value_local(eta, sensors, 3);
    rows(i).weighted_y_rmse_mm = y_rmse;
    rows(i).weighted_voltage_rmse_v = wrmse_v;
    rows(i).plain_voltage_rmse_v = prmse_v;
    rows(i).point_count = nnz(valid);
    rows(i).valid_query_fraction = valid_q_frac;
end

T = struct2table(rows);
T = T(isfinite(T.EO) & isfinite(T.weighted_voltage_rmse_v), :);
if isempty(T)
    return;
end
T = sortrows(T, {'weighted_voltage_rmse_v', 'weighted_y_rmse_mm', 'EO'}, ...
    {'ascend', 'ascend', 'ascend'});
T.rank = (1:height(T)).';
end


function T = scan_sg_like_eo_local(B, eo_candidates, S)
Cal = S.sg_static_calibration;
T = scan_parametric_static_calibration_eo_local(B, eo_candidates, S, Cal);
end


function T = scan_step04_parametric_sg_eo_local(B, eo_candidates, S)
Cal = build_step04_parametric_sg_calibration_local(B, S);
T = scan_parametric_static_calibration_eo_local(B, eo_candidates, S, Cal);
end


function T = scan_parametric_static_calibration_eo_local(B, eo_candidates, S, Cal)
if isempty(Cal)
    T = table();
    return;
end

x = B.x(:);
v = B.v(:);
theta = B.theta(:);
sensor_id = B.sensor_id(:);
w = B.fit_weight(:);
w(~isfinite(w)) = S.weight_floor;
w = max(w, S.weight_floor);

x_stat = nan(size(x));
W_sg = nan(size(x));
for i = 1:numel(Cal)
    sid = Cal(i).sensor_id;
    idx = sensor_id == sid;
    if ~any(idx)
        continue;
    end

    x_centered = x(idx) - Cal(i).xc;
    v_above = max(v(idx) - Cal(i).baseline, 1e-6);
    v_clip = min(max(v_above ./ max(Cal(i).B, 1e-6), 1e-6), 0.999999);
    sgn = sign(x_centered);
    sgn(sgn == 0) = 1;
    x_stat(idx) = Cal(i).xc + sgn .* ...
        Cal(i).w .* (-log(v_clip)).^(1 ./ Cal(i).n);
    if isfield(Cal(i), 'x_ideal') && isfield(Cal(i), 'W_ideal') && ...
            ~isempty(Cal(i).x_ideal) && ~isempty(Cal(i).W_ideal)
        W_sg(idx) = interp1(Cal(i).x_ideal(:), Cal(i).W_ideal(:), ...
            x_centered, 'linear', S.weight_floor);
    end
end

y_obs = x - x_stat;
W_sg(~isfinite(W_sg)) = w(~isfinite(W_sg));
W_sg(~isfinite(W_sg)) = S.weight_floor;
W_sg = max(W_sg, S.weight_floor);

valid_base = isfinite(x) & isfinite(v) & isfinite(theta) & ...
    isfinite(y_obs) & isfinite(W_sg);

rows = repmat(struct( ...
    'EO', NaN, ...
    'rank', NaN, ...
    'A_mm', NaN, ...
    'phi_rad', NaN, ...
    'd0_mm', NaN, ...
    'eta_CH1_mm', NaN, ...
    'eta_CH2_mm', NaN, ...
    'eta_CH3_mm', NaN, ...
    'weighted_y_rmse_mm', inf, ...
    'weighted_voltage_rmse_v', inf, ...
    'plain_voltage_rmse_v', inf, ...
    'point_count', nnz(valid_base), ...
    'valid_query_fraction', NaN), numel(eo_candidates), 1);

for i = 1:numel(eo_candidates)
    EO = eo_candidates(i);
    X = [sin(EO .* theta), cos(EO .* theta), ones(size(theta))];
    valid = valid_base & all(isfinite(X), 2);
    if nnz(valid) < 20
        continue;
    end

    wg = sqrt(W_sg(valid));
    beta = (X(valid, :) .* wg) \ (y_obs(valid) .* wg);
    if numel(beta) < 3 || any(~isfinite(beta))
        continue;
    end

    y_fit = X * beta;
    [wrmse_v, prmse_v, valid_fraction] = evaluate_sg_static_voltage_rmse_local( ...
        B, Cal, y_fit, W_sg);
    y_res = y_obs(valid) - X(valid, :) * beta;
    y_rmse = sqrt(sum(W_sg(valid) .* y_res.^2) ./ ...
        max(sum(W_sg(valid)), eps));

    rows(i).EO = EO;
    rows(i).A_mm = hypot(beta(1), beta(2));
    rows(i).phi_rad = wrap_to_pi_local(atan2(beta(2), beta(1)));
    rows(i).d0_mm = beta(3);
    rows(i).weighted_y_rmse_mm = y_rmse;
    rows(i).weighted_voltage_rmse_v = wrmse_v;
    rows(i).plain_voltage_rmse_v = prmse_v;
    rows(i).point_count = nnz(valid);
    rows(i).valid_query_fraction = valid_fraction;
end

T = struct2table(rows);
T = T(isfinite(T.EO) & isfinite(T.weighted_voltage_rmse_v), :);
if isempty(T)
    return;
end
T = sortrows(T, {'weighted_voltage_rmse_v', 'weighted_y_rmse_mm', 'EO'}, ...
    {'ascend', 'ascend', 'ascend'});
T.rank = (1:height(T)).';
end


function Cal = build_step04_parametric_sg_calibration_local(B, S)
Cal = [];
if ~isfield(B, 'TemplateMini') || isempty(B.TemplateMini)
    return;
end

sensors = double(S.analysis_sensors(:).');
template = struct( ...
    'sensor_id', NaN, ...
    'alpha_ref', NaN, ...
    'dx_trust_left', NaN, ...
    'dx_trust_right', NaN, ...
    'L_opt', NaN, ...
    'xc_seed', NaN, ...
    'B', NaN, ...
    'w', NaN, ...
    'n', NaN, ...
    'xc', NaN, ...
    'baseline', NaN, ...
    'x_ideal', [], ...
    'W_ideal', []);

rows = repmat(template, 0, 1);
for sid = sensors
    idx = find([B.TemplateMini.sensor_id] == sid, 1, 'first');
    if isempty(idx)
        continue;
    end
    row = fit_step04_template_sg_local(B.TemplateMini(idx), template);
    if isfinite(row.B) && isfinite(row.w) && isfinite(row.n) && ...
            isfinite(row.xc) && isfinite(row.baseline)
        rows(end + 1, 1) = row; %#ok<AGROW>
    end
end
Cal = rows;
end


function row = fit_step04_template_sg_local(tpl, template)
row = template;
row.sensor_id = double(tpl.sensor_id);

x = tpl.x_grid(:);
v = tpl.v_grid(:);
if isfield(tpl, 'weight_grid') && ~isempty(tpl.weight_grid)
    wgt = tpl.weight_grid(:);
else
    wgt = ones(size(x));
end

valid = isfinite(x) & isfinite(v) & isfinite(wgt);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && ...
        all(isfinite(tpl.x_domain(1:2)))
    valid = valid & x >= tpl.x_domain(1) & x <= tpl.x_domain(2);
    row.dx_trust_left = abs(tpl.x_domain(1));
    row.dx_trust_right = abs(tpl.x_domain(2));
end

x = x(valid);
v = v(valid);
wgt = max(wgt(valid), 0.05);
if numel(x) < 20
    return;
end

if isfield(tpl, 'baseline') && isfinite(tpl.baseline)
    base0 = tpl.baseline;
else
    base0 = prctile(v, 5);
end
[vmax, imax] = max(v);
if ~isfinite(vmax) || ~isfinite(base0)
    return;
end

B0 = max(vmax - base0, 0.1);
w0 = max(0.25 * range(x), 0.5);
n0 = 2.2;
xc0 = x(imax);
if isfield(tpl, 'xc') && isfinite(tpl.xc)
    xc0 = min(max(xc0, -0.4), 0.4);
end

base_lo = min(v) - 0.5 * max(B0, 0.1);
base_hi = min(v) + 0.25 * max(B0, 0.1);
x_span = max(max(abs(x)), 1);

p0 = [log(B0), log(w0), log(n0), xc0, base0];
obj = @(p) step04_sg_template_fit_objective_local( ...
    p, x, v, wgt, base_lo, base_hi, x_span);
opts = optimset('Display', 'off', 'MaxIter', 600, ...
    'MaxFunEvals', 1600, 'TolX', 1e-8, 'TolFun', 1e-10);

try
    p = fminsearch(obj, p0, opts);
catch
    p = p0;
end

[Bfit, wfit, nfit, xcfit, basefit] = unpack_step04_sg_params_local( ...
    p, base_lo, base_hi, x_span);

row.B = Bfit;
row.w = wfit;
row.n = nfit;
row.xc = xcfit;
row.xc_seed = xc0;
row.baseline = basefit;
row.L_opt = max(row.dx_trust_left, row.dx_trust_right);
row.x_ideal = x(:);
row.W_ideal = wgt(:) ./ max(wgt(:));
end


function obj = step04_sg_template_fit_objective_local( ...
    p, x, v, wgt, base_lo, base_hi, x_span)
[B, w, n, xc, base] = unpack_step04_sg_params_local( ...
    p, base_lo, base_hi, x_span);
v_fit = B .* exp(-abs((x(:) - xc) ./ max(w, 1e-6)).^n) + base;
res = v(:) - v_fit(:);
obj = sum(wgt(:) .* res.^2) ./ max(sum(wgt(:)), eps);
if ~isfinite(obj)
    obj = 1e12;
end
end


function [B, w, n, xc, base] = unpack_step04_sg_params_local( ...
    p, base_lo, base_hi, x_span)
p = p(:).';
B = min(max(exp(p(1)), 0.02), 20);
w = min(max(exp(p(2)), 0.05), 20);
n = min(max(exp(p(3)), 0.5), 10);
xc = min(max(p(4), -x_span), x_span);
base = min(max(p(5), base_lo), base_hi);
end


function [wrmse, prmse, valid_fraction] = evaluate_sg_static_voltage_rmse_local(B, Cal, y_fit, W)
xq = B.x(:) - y_fit(:);
Vpred = nan(size(xq));
sensor_id = B.sensor_id(:);

for i = 1:numel(Cal)
    sid = Cal(i).sensor_id;
    idx = sensor_id == sid;
    if ~any(idx)
        continue;
    end
    Vpred(idx) = Cal(i).B .* ...
        exp(-abs((xq(idx) - Cal(i).xc) ./ Cal(i).w).^Cal(i).n) + ...
        Cal(i).baseline;
end

valid = isfinite(Vpred) & isfinite(B.v(:)) & isfinite(W(:));
valid_fraction = nnz(valid) / max(numel(xq), 1);
if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
    return;
end

res = B.v(:) - Vpred(:);
wrmse = sqrt(sum(W(valid) .* res(valid).^2) ./ max(sum(W(valid)), eps));
prmse = sqrt(mean(res(valid).^2));
end


function x_static = invert_template_same_side_local(B)
x_static = nan(size(B.x(:)));
sensors = unique(B.sensor_id(:).');

for sid = sensors
    if ~isfinite(sid)
        continue;
    end
    idx_sensor = B.sensor_id(:) == sid;
    blades = unique(B.blade_id_vec(idx_sensor).');
    for bid = blades
        idx = idx_sensor & B.blade_id_vec(:) == bid;
        if ~any(idx)
            continue;
        end

        tpl = get_bundle_template_local(B, sid, bid);
        if isempty(tpl)
            continue;
        end

        x_obs = B.x(idx);
        v_obs = B.v(idx);
        side = sign(x_obs);
        side(side == 0) = 1;

        x_out = nan(size(x_obs));
        x_out(side < 0) = inverse_one_side_local( ...
            tpl, v_obs(side < 0), -1);
        x_out(side >= 0) = inverse_one_side_local( ...
            tpl, v_obs(side >= 0), 1);
        x_static(idx) = x_out;
    end
end
end


function xq = inverse_one_side_local(tpl, v_obs, side_sign)
xq = nan(size(v_obs));
if isempty(v_obs)
    return;
end

xg = tpl.x_grid(:);
vg = tpl.v_grid(:);
valid = isfinite(xg) & isfinite(vg);
if isfield(tpl, 'x_domain') && numel(tpl.x_domain) >= 2 && ...
        all(isfinite(tpl.x_domain(1:2)))
    valid = valid & xg >= tpl.x_domain(1) & xg <= tpl.x_domain(2);
end

if side_sign < 0
    valid = valid & xg <= 0;
    [x_path, order] = sort(xg(valid), 'ascend');
else
    valid = valid & xg >= 0;
    [x_path, order] = sort(xg(valid), 'descend');
end
vg_side = vg(valid);
v_path = vg_side(order);

good = isfinite(x_path) & isfinite(v_path);
x_path = x_path(good);
v_path = v_path(good);
if numel(x_path) < 3
    return;
end

v_mono = cummax(v_path);
[v_unique, ia] = unique(v_mono, 'stable');
x_unique = x_path(ia);
if numel(v_unique) < 2 || v_unique(end) <= v_unique(1)
    [~, nearest_idx] = min(abs(v_path(:).' - v_obs(:)), [], 2);
    xq(:) = x_path(nearest_idx);
    return;
end

v_clip = min(max(v_obs(:), v_unique(1)), v_unique(end));
xq(:) = interp1(v_unique, x_unique, v_clip, 'linear', 'extrap');
end


function [wrmse, prmse, valid_fraction] = evaluate_voltage_rmse_local(B, y_fit, w, S)
xq = B.x(:) - y_fit(:);
Vpred = nan(size(xq));

sensors = unique(B.sensor_id(:).');
for sid = sensors
    idx_sensor = B.sensor_id(:) == sid;
    blades = unique(B.blade_id_vec(idx_sensor).');
    for bid = blades
        idx = idx_sensor & B.blade_id_vec(:) == bid;
        tpl = get_bundle_template_local(B, sid, bid);
        if isempty(tpl)
            continue;
        end
        Vpred(idx) = interp1(tpl.x_grid(:), tpl.v_grid(:), xq(idx), ...
            S.template_interp_method, NaN);
    end
end

valid = isfinite(Vpred) & isfinite(B.v(:)) & isfinite(w(:));
valid_fraction = nnz(valid) / max(numel(xq), 1);
if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
    return;
end

res = B.v(:) - Vpred(:);
wrmse = sqrt(sum(w(valid) .* res(valid).^2) ./ max(sum(w(valid)), eps));
prmse = sqrt(mean(res(valid).^2));
end


function tpl = get_bundle_template_local(B, sid, bid)
tpl = [];
if ~isfield(B, 'TemplateMini') || isempty(B.TemplateMini)
    return;
end
for i = 1:numel(B.TemplateMini)
    if double(B.TemplateMini(i).sensor_id) == double(sid) && ...
            double(B.TemplateMini(i).blade_id) == double(bid)
        tpl = B.TemplateMini(i);
        return;
    end
end
end


function row = summarize_window_local(iw, WR, T)
row = make_empty_summary_row_local();
row.window_id = iw;
row.status = "ok";

if isfield(WR, 'Result') && isstruct(WR.Result)
    row.direct_fit_eo = get_field_or_default_local(WR.Result, 'EO_id', NaN);
    row.direct_fit_rmse_v = get_field_or_default_local( ...
        WR.Result, 'weighted_voltage_rmse', NaN);
    row.direct_fit_top_list = string(get_top_list_from_candidate_table_local( ...
        get_field_or_default_local(WR.Result, 'CandidateTable', table())));
end

models = ["template_inversion_common_offset", ...
    "template_inversion_sensor_static_offset"];
for model = models
    M = T(strcmpi(T.model, model), :);
    if isempty(M)
        continue;
    end
    top_list = strjoin(compose('%d', M.EO(1:min(5, height(M))).'), "|");
    eo14_rank = NaN;
    eo14_rmse = NaN;
    idx14 = find(M.EO == 14, 1, 'first');
    if ~isempty(idx14)
        eo14_rank = M.rank(idx14);
        eo14_rmse = M.weighted_voltage_rmse_v(idx14);
    end

    if model == "template_inversion_common_offset"
        row.inv_common_eo = M.EO(1);
        row.inv_common_rmse_v = M.weighted_voltage_rmse_v(1);
        row.inv_common_y_rmse_mm = M.weighted_y_rmse_mm(1);
        row.inv_common_eo14_rank = eo14_rank;
        row.inv_common_eo14_rmse_v = eo14_rmse;
        row.inv_common_top_list = string(top_list);
    else
        row.inv_eta_eo = M.EO(1);
        row.inv_eta_rmse_v = M.weighted_voltage_rmse_v(1);
        row.inv_eta_y_rmse_mm = M.weighted_y_rmse_mm(1);
        row.inv_eta_eo14_rank = eo14_rank;
        row.inv_eta_eo14_rmse_v = eo14_rmse;
        row.inv_eta_top_list = string(top_list);
        row.inv_eta_CH2_mm = M.eta_CH2_mm(1);
        row.inv_eta_CH3_mm = M.eta_CH3_mm(1);
    end
end

M = T(strcmpi(T.model, "sg_static_model_replay"), :);
if ~isempty(M)
    top_list = strjoin(compose('%d', M.EO(1:min(5, height(M))).'), "|");
    idx14 = find(M.EO == 14, 1, 'first');
    row.sg_replay_eo = M.EO(1);
    row.sg_replay_rmse_v = M.weighted_voltage_rmse_v(1);
    row.sg_replay_y_rmse_mm = M.weighted_y_rmse_mm(1);
    row.sg_replay_top_list = string(top_list);
    if ~isempty(idx14)
        row.sg_replay_eo14_rank = M.rank(idx14);
        row.sg_replay_eo14_rmse_v = M.weighted_voltage_rmse_v(idx14);
    end
end

M = T(strcmpi(T.model, "step04_parametric_sg_replay"), :);
if ~isempty(M)
    top_list = strjoin(compose('%d', M.EO(1:min(5, height(M))).'), "|");
    idx14 = find(M.EO == 14, 1, 'first');
    row.step04_sg_eo = M.EO(1);
    row.step04_sg_rmse_v = M.weighted_voltage_rmse_v(1);
    row.step04_sg_y_rmse_mm = M.weighted_y_rmse_mm(1);
    row.step04_sg_top_list = string(top_list);
    if ~isempty(idx14)
        row.step04_sg_eo14_rank = M.rank(idx14);
        row.step04_sg_eo14_rmse_v = M.weighted_voltage_rmse_v(idx14);
    end
end
end


function row = make_empty_summary_row_local()
row = struct( ...
    'window_id', NaN, ...
    'status', "", ...
    'direct_fit_eo', NaN, ...
    'direct_fit_rmse_v', NaN, ...
    'direct_fit_top_list', "", ...
    'inv_common_eo', NaN, ...
    'inv_common_rmse_v', NaN, ...
    'inv_common_y_rmse_mm', NaN, ...
    'inv_common_eo14_rank', NaN, ...
    'inv_common_eo14_rmse_v', NaN, ...
    'inv_common_top_list', "", ...
    'inv_eta_eo', NaN, ...
    'inv_eta_rmse_v', NaN, ...
    'inv_eta_y_rmse_mm', NaN, ...
    'inv_eta_eo14_rank', NaN, ...
    'inv_eta_eo14_rmse_v', NaN, ...
    'inv_eta_top_list', "", ...
    'inv_eta_CH2_mm', NaN, ...
    'inv_eta_CH3_mm', NaN, ...
    'sg_replay_eo', NaN, ...
    'sg_replay_rmse_v', NaN, ...
    'sg_replay_y_rmse_mm', NaN, ...
    'sg_replay_eo14_rank', NaN, ...
    'sg_replay_eo14_rmse_v', NaN, ...
    'sg_replay_top_list', "", ...
    'step04_sg_eo', NaN, ...
    'step04_sg_rmse_v', NaN, ...
    'step04_sg_y_rmse_mm', NaN, ...
    'step04_sg_eo14_rank', NaN, ...
    'step04_sg_eo14_rmse_v', NaN, ...
    'step04_sg_top_list', "");
end


function top_list = get_top_list_from_candidate_table_local(C)
top_list = "";
if isempty(C) || ~istable(C) || ~ismember('EO', C.Properties.VariableNames)
    return;
end
if ismember('weighted_voltage_rmse', C.Properties.VariableNames)
    C = sortrows(C, {'weighted_voltage_rmse', 'EO'}, {'ascend', 'ascend'});
end
top_list = strjoin(compose('%d', C.EO(1:min(5, height(C))).'), "|");
end


function eo_candidates = build_eo_candidates_local(rot_freq_mean_hz, freq_search_hz, eo_pad)
freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min = max(1, ceil(freq_lo / rot_freq_mean_hz) - max(0, eo_pad));
eo_max = max(eo_min, floor(freq_hi / rot_freq_mean_hz) + max(0, eo_pad));
eo_candidates = eo_min:eo_max;
eo_candidates = eo_candidates(eo_candidates > 0);
end


function value = eta_value_local(eta, sensors, sid)
value = NaN;
idx = find(sensors == sid, 1, 'first');
if ~isempty(idx) && idx <= numel(eta)
    value = eta(idx);
end
end


function y = wrap_to_pi_local(x)
y = mod(x + pi, 2 * pi) - pi;
end


function value = get_field_or_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function Cal = load_sg_static_calibration_local()
Cal = [];
sg_root = resolve_sg_root_local(fileparts(mfilename('fullpath')));
if isempty(sg_root) || ~isfolder(sg_root)
    return;
end

sensor_ids = [1 2 3];
template = struct( ...
    'sensor_id', NaN, ...
    'alpha_ref', NaN, ...
    'dx_trust_left', NaN, ...
    'dx_trust_right', NaN, ...
    'L_opt', NaN, ...
    'xc_seed', NaN, ...
    'B', NaN, ...
    'w', NaN, ...
    'n', NaN, ...
    'xc', NaN, ...
    'baseline', NaN, ...
    'x_ideal', [], ...
    'W_ideal', []);
Cal = repmat(template, 0, 1);
for sid = sensor_ids
    file = fullfile(sg_root, sprintf('SGCalib_20251222_B1_S%d.mat', sid));
    if ~isfile(file)
        Cal = [];
        return;
    end
    loaded = load(file, 'Calib');
    C = loaded.Calib.Sensor;
    row = template;
    names = fieldnames(row);
    for i = 1:numel(names)
        if isfield(C, names{i})
            row.(names{i}) = C.(names{i});
        end
    end
    Cal(end + 1, 1) = row; %#ok<AGROW>
end
end


function sg_root = resolve_sg_root_local(route_dir)
env_root = strtrim(getenv('SG_SUCCESS_ROOT'));
if ~isempty(env_root) && isfolder(env_root)
    sg_root = env_root;
    return;
end

validation_root = fileparts(route_dir);
repo_root = fileparts(validation_root);
common_root = fileparts(repo_root);
candidate = fullfile(common_root, ...
    '7超高斯模型-权重-瞬态', '程序', '直叶片验证', '实验验证', ...
    '超高斯', '20251222适配-3号传感间隙变化-改');
if isfolder(candidate)
    sg_root = candidate;
else
    sg_root = '';
end
end


function print_findings_local(WindowSummary)
ok = strcmpi(WindowSummary.status, "ok");
if ~any(ok)
    fprintf('No valid diagnostic windows.\n');
    return;
end

fprintf('\nDiagnostic EO sequences:\n');
fprintf('  direct fit:                    %s\n', ...
    strjoin(compose('%d', WindowSummary.direct_fit_eo(ok).'), "|"));
fprintf('  template inversion common:     %s\n', ...
    strjoin(compose('%d', WindowSummary.inv_common_eo(ok).'), "|"));
fprintf('  template inversion sensor eta: %s\n', ...
    strjoin(compose('%d', WindowSummary.inv_eta_eo(ok).'), "|"));
if ismember('sg_replay_eo', WindowSummary.Properties.VariableNames)
    fprintf('  SG static model replay:       %s\n', ...
        strjoin(compose('%d', WindowSummary.sg_replay_eo(ok).'), "|"));
end
if ismember('step04_sg_eo', WindowSummary.Properties.VariableNames)
    fprintf('  Step04 parametric SG replay:  %s\n', ...
        strjoin(compose('%d', WindowSummary.step04_sg_eo(ok).'), "|"));
end

fprintf('\nEO14 median rank:\n');
fprintf('  common offset:        %.2f\n', ...
    median(WindowSummary.inv_common_eo14_rank(ok), 'omitnan'));
fprintf('  sensor static offset: %.2f\n', ...
    median(WindowSummary.inv_eta_eo14_rank(ok), 'omitnan'));
if ismember('sg_replay_eo14_rank', WindowSummary.Properties.VariableNames)
    fprintf('  SG static replay:     %.2f\n', ...
        median(WindowSummary.sg_replay_eo14_rank(ok), 'omitnan'));
end
if ismember('step04_sg_eo14_rank', WindowSummary.Properties.VariableNames)
    fprintf('  Step04 param SG:      %.2f\n', ...
        median(WindowSummary.step04_sg_eo14_rank(ok), 'omitnan'));
end

fprintf('\nSensor-offset model eta medians:\n');
fprintf('  CH2 %.6f mm, CH3 %.6f mm\n', ...
    median(WindowSummary.inv_eta_CH2_mm(ok), 'omitnan'), ...
    median(WindowSummary.inv_eta_CH3_mm(ok), 'omitnan'));
end
