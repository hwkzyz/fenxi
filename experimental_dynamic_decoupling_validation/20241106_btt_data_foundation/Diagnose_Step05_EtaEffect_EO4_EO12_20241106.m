clc; close all;

%DIAGNOSE_STEP05_ETAEFFECT_EO4_EO12_20241106
% Diagnostic only.  This script tests whether the early-window EO4 result
% is caused by fixing eta_s=0 in the formal Step05 NoEta route.
%
% It does not change the official Step05 result.  It reloads the saved
% Step05 result and refits EO4/EO12 on the stored core/final bundle previews
% with two models:
%   NoEta   : [A, phi, dx_c]
%   WithEta : [A, phi, dx_c, eta_2, eta_3, ...]

cfg = BTTProjectConfig_20241106();

S = struct();
S.case_name = cfg.dynamic_cases{1};
S.blade_id = cfg.step05_target_blades(1);
S.analysis_sensors = cfg.step05_analysis_sensors;
S.target_eos = [4 12];
S.bundle_sources = {'CoreBundlePreview', 'BundlePreview'};
S.result_file = fullfile(cfg.step05_output_dir, S.case_name, ...
    sprintf('Result_Step05_FoundationMainPulseAdaptiveNoEta_B%d_S%s_%s.mat', ...
    S.blade_id, sprintf('%d', S.analysis_sensors), cfg.dataset));
S.output_dir = fullfile(cfg.output_root, 'step05_eta_effect_diagnostics', S.case_name);
S.amplitude_limit_mm = 0.50;
S.dx_c_limit_mm = 0.35;
S.sensor_eta_limit_mm = cfg.step05_sensor_eta_limit_mm;
S.sensor_eta_reg_weight_v_per_mm = cfg.step05_sensor_eta_reg_weight_v_per_mm;
S.overshoot_penalty_weight = 100;
S.invalid_query_penalty_scale = 1e6;
S.weight_floor = 0.05;
S.fminsearch_max_iter = 250;
S.fminsearch_max_fun = 700;
S.template_interp_method = 'pchip';
S.make_plots = true;
S.window_ids = 1:8;

if ~isfile(S.result_file)
    error('Step05 result not found. Run Step05 first:\n  %s', S.result_file);
end

if exist(S.output_dir, 'dir') ~= 7
    mkdir(S.output_dir);
end

loaded = load(S.result_file, 'Result');
Result = loaded.Result;

window_env = strtrim(getenv('ETA_DIAG_WINDOW_IDS'));
if ~isempty(window_env)
    if any(strcmpi(window_env, {'all', 'inf'}))
        S.window_ids = 1:numel(Result.WindowResult);
    else
        tokens = regexp(window_env, '[,;\s]+', 'split');
        ids = [];
        for i = 1:numel(tokens)
            v = str2double(tokens{i});
            if isfinite(v) && v >= 1
                ids(end+1) = round(v); %#ok<AGROW>
            end
        end
        if ~isempty(ids)
            S.window_ids = unique(ids, 'stable');
        end
    end
end

S.window_ids = S.window_ids(S.window_ids >= 1 & ...
    S.window_ids <= numel(Result.WindowResult));

rows = repmat(make_empty_eta_diag_row_local(), 0, 1);

for iw = S.window_ids
    Wres = Result.WindowResult(iw);

    for isource = 1:numel(S.bundle_sources)
        source_name = S.bundle_sources{isource};
        if ~isfield(Wres, source_name) || isempty(Wres.(source_name))
            continue;
        end

        bundle = Wres.(source_name);
        if ~isstruct(bundle) || ~isfield(bundle, 'x') || isempty(bundle.x)
            continue;
        end

        for ieo = 1:numel(S.target_eos)
            eo = S.target_eos(ieo);
            seed = get_candidate_seed_local(Wres.Result, eo, S);

            fit_noeta = fit_one_eo_eta_model_local(bundle, eo, seed, S, false);
            fit_witheta = fit_one_eo_eta_model_local(bundle, eo, seed, S, true);

            rows(end+1, 1) = make_eta_diag_row_local( ...
                iw, source_name, eo, 'NoEta', fit_noeta, bundle, Wres); %#ok<SAGROW>
            rows(end+1, 1) = make_eta_diag_row_local( ...
                iw, source_name, eo, 'WithEta', fit_witheta, bundle, Wres); %#ok<SAGROW>
        end
    end
end

T = struct2table(rows);

if isempty(T)
    error('No diagnostic rows were produced.');
end

T = add_pairwise_winner_columns_local(T, S);

csv_file = fullfile(S.output_dir, 'EtaEffect_EO4_EO12_20241106.csv');
mat_file = fullfile(S.output_dir, 'EtaEffect_EO4_EO12_20241106.mat');
writetable(T, csv_file);
save(mat_file, 'T', 'S', 'Result', '-v7.3');

fprintf('\n=== Step05 eta_s effect diagnostic: 20241106 B%d/S%s ===\n', ...
    S.blade_id, sprintf('%d', S.analysis_sensors));
fprintf('Result source:\n  %s\n', S.result_file);
fprintf('Saved diagnostic table:\n  %s\n', csv_file);

print_eta_effect_summary_local(T, S);

if S.make_plots
    plot_eta_effect_summary_local(T, S);
end


function row = make_empty_eta_diag_row_local()
row = struct( ...
    'window_id', NaN, ...
    'bundle_source', "", ...
    'model', "", ...
    'EO', NaN, ...
    'weighted_rmse', NaN, ...
    'plain_rmse', NaN, ...
    'objective_score', NaN, ...
    'A_mm', NaN, ...
    'phi_rad', NaN, ...
    'dx_c_mm', NaN, ...
    'eta_max_abs_mm', NaN, ...
    'eta_rms_mm', NaN, ...
    'eta_vector', "", ...
    'point_count', NaN, ...
    'pass_count', NaN, ...
    'original_window_eo', NaN, ...
    'original_window_rmse', NaN, ...
    'winner_eo_same_model_source', NaN, ...
    'delta_rmse_EO12_minus_EO4', NaN, ...
    'eta_changes_winner', false);
end


function seed = get_candidate_seed_local(result, eo, S)
seed = struct();
seed.A = 0.05;
seed.phi = 0;
seed.dx_c = 0;
seed.eta = zeros(1, numel(S.analysis_sensors));

if ~isstruct(result) || ~isfield(result, 'CandidateTable') || ...
        isempty(result.CandidateTable) || ~istable(result.CandidateTable)
    return;
end

C = result.CandidateTable;
idx = find(C.EO == eo, 1, 'first');
if isempty(idx)
    return;
end

if ismember('A', C.Properties.VariableNames)
    seed.A = C.A(idx);
end
if ismember('phi', C.Properties.VariableNames)
    seed.phi = C.phi(idx);
end
if ismember('dx_c', C.Properties.VariableNames)
    seed.dx_c = C.dx_c(idx);
end
if ismember('sensor_eta', C.Properties.VariableNames)
    eta = C.sensor_eta(idx, :);
    if isnumeric(eta) && numel(eta) == numel(S.analysis_sensors)
        seed.eta = eta(:).';
    end
end

seed.A = min(max(abs(seed.A), 0), S.amplitude_limit_mm);
if ~isfinite(seed.A)
    seed.A = 0.05;
end
seed.phi = wrap_to_pi_local(seed.phi);
if ~isfinite(seed.phi)
    seed.phi = 0;
end
seed.dx_c = clamp_scalar_local(seed.dx_c, S.dx_c_limit_mm);
if ~isfinite(seed.dx_c)
    seed.dx_c = 0;
end
seed.eta = zeros(1, numel(S.analysis_sensors));
end


function fit = fit_one_eo_eta_model_local(bundle, eo, seed, S, fit_eta)
starts = build_start_points_local(seed, S, fit_eta);
opts = optimset( ...
    'Display', 'off', ...
    'MaxIter', S.fminsearch_max_iter, ...
    'MaxFunEvals', S.fminsearch_max_fun, ...
    'TolX', 1e-7, ...
    'TolFun', 1e-9);

best_p = starts(1, :);
best_obj = inf;
best_exitflag = -1;

for istart = 1:size(starts, 1)
    p0 = starts(istart, :);
    obj = @(p) eta_diag_objective_local(p, bundle, eo, S, fit_eta);

    try
        [p_try, obj_try, exitflag_try] = fminsearch(obj, p0, opts);
    catch
        p_try = p0;
        obj_try = obj(p0);
        exitflag_try = -1;
    end

    if isfinite(obj_try) && obj_try < best_obj
        best_obj = obj_try;
        best_p = p_try;
        best_exitflag = exitflag_try;
    end
end

[A, phi, dx_c, eta] = unpack_eta_diag_params_local(best_p, S, fit_eta);
[wrmse, prmse, objective_score] = evaluate_eta_diag_model_local( ...
    bundle, eo, A, phi, dx_c, eta, S);

fit = struct();
fit.EO = eo;
fit.A = A;
fit.phi = phi;
fit.dx_c = dx_c;
fit.eta = eta;
fit.weighted_rmse = wrmse;
fit.plain_rmse = prmse;
fit.objective_score = objective_score;
fit.exitflag = best_exitflag;
fit.fit_eta = fit_eta;
end


function starts = build_start_points_local(seed, S, fit_eta)
base = [seed.A, seed.phi, seed.dx_c];
if fit_eta
    base = [base, zeros(1, numel(S.analysis_sensors) - 1)];
end

starts = base;

amp_list = unique([seed.A, 0.12]);
phi_list = unique([seed.phi, seed.phi + pi/2]);

for ia = 1:numel(amp_list)
    for ip = 1:numel(phi_list)
        p = base;
        p(1) = min(max(abs(amp_list(ia)), 0), S.amplitude_limit_mm);
        p(2) = wrap_to_pi_local(phi_list(ip));
        starts(end+1, :) = p; %#ok<AGROW>
    end
end

starts = unique(round(starts, 10), 'rows', 'stable');
end


function value = eta_diag_objective_local(p, bundle, eo, S, fit_eta)
[A, phi, dx_c, eta] = unpack_eta_diag_params_local(p, S, fit_eta);
[~, ~, objective_score] = evaluate_eta_diag_model_local( ...
    bundle, eo, A, phi, dx_c, eta, S);

if ~isfinite(objective_score)
    objective_score = 1e12;
end

value = objective_score;
end


function [A, phi, dx_c, eta] = unpack_eta_diag_params_local(p, S, fit_eta)
p = p(:).';

A = min(max(abs(p(1)), 0), S.amplitude_limit_mm);
phi = wrap_to_pi_local(p(2));
dx_c = clamp_scalar_local(p(3), S.dx_c_limit_mm);

eta = zeros(1, numel(S.analysis_sensors));
if fit_eta && numel(p) > 3
    eta(2:end) = p(4:min(numel(p), 2 + numel(S.analysis_sensors)));
end
eta = clamp_vector_local(eta, S.sensor_eta_limit_mm);
eta(1) = 0;
eta(~isfinite(eta)) = 0;
end


function [wrmse, prmse, objective_score] = evaluate_eta_diag_model_local( ...
    bundle, eo, A, phi, dx_c, eta, S)
x = bundle.x(:);
v = bundle.v(:);
theta = bundle.theta(:);
w = bundle.fit_weight(:);
sensors = S.analysis_sensors(:).';
eta_vec = sensor_eta_vector_local(bundle.sensor_id(:), sensors, eta);

xq = x - dx_c - eta_vec - A .* sin(eo .* theta + phi);
Vpred = nan(size(xq));
overshoot = zeros(size(xq));

for sid = sensors
    idx_sensor = bundle.sensor_id(:) == sid;
    if ~any(idx_sensor)
        continue;
    end

    blade_ids = unique(bundle.blade_id_vec(idx_sensor));
    for ib = 1:numel(blade_ids)
        bid = blade_ids(ib);
        idx = idx_sensor & bundle.blade_id_vec(:) == bid;
        tpl = get_bundle_template_local(bundle, sid, bid);
        if isempty(tpl)
            continue;
        end

        left = tpl.x_domain(1);
        right = tpl.x_domain(2);
        overshoot(idx) = max(left - xq(idx), 0) + max(xq(idx) - right, 0);
        Vpred(idx) = interp1(tpl.x_grid(:), tpl.v_grid(:), xq(idx), ...
            S.template_interp_method, NaN);
    end
end

valid = isfinite(Vpred) & isfinite(v) & isfinite(w);

if nnz(valid) < 5
    wrmse = inf;
    prmse = inf;
    objective_score = inf;
    return;
end

w(~isfinite(w)) = S.weight_floor;
w = max(w, S.weight_floor);
w = w ./ max(max(w), eps);
res = v - Vpred;

residual_obj = nansum(w(valid) .* res(valid).^2);
overshoot_obj = S.overshoot_penalty_weight * ...
    nansum(w(valid) .* overshoot(valid).^2);
eta_reg_obj = numel(x) * ...
    (S.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta(:).^2)))^2;
missing_obj = S.invalid_query_penalty_scale * nnz(~valid);

objective_score = residual_obj + overshoot_obj + eta_reg_obj + missing_obj;
wrmse = sqrt(objective_score / max(numel(x), 1));
prmse = sqrt(nanmean(res(valid).^2));
end


function tpl = get_bundle_template_local(bundle, sid, bid)
tpl = [];
if ~isfield(bundle, 'TemplateMini') || isempty(bundle.TemplateMini)
    return;
end

TM = bundle.TemplateMini;
for i = 1:numel(TM)
    if isequal(TM(i).sensor_id, sid) && isequal(TM(i).blade_id, bid)
        tpl = TM(i);
        return;
    end
end
end


function eta_vec = sensor_eta_vector_local(sensor_id, sensors, eta)
eta_vec = zeros(size(sensor_id));
for i = 1:numel(sensors)
    eta_vec(sensor_id == sensors(i)) = eta(i);
end
end


function row = make_eta_diag_row_local(iw, source_name, eo, model_name, fit, bundle, Wres)
row = make_empty_eta_diag_row_local();
row.window_id = iw;
row.bundle_source = string(source_name);
row.model = string(model_name);
row.EO = eo;
row.weighted_rmse = fit.weighted_rmse;
row.plain_rmse = fit.plain_rmse;
row.objective_score = fit.objective_score;
row.A_mm = fit.A;
row.phi_rad = fit.phi;
row.dx_c_mm = fit.dx_c;
row.eta_max_abs_mm = max(abs(fit.eta), [], 'omitnan');
row.eta_rms_mm = sqrt(mean(fit.eta(:).^2, 'omitnan'));
row.eta_vector = string(mat2str(fit.eta, 5));
row.point_count = numel(bundle.x);
row.pass_count = numel(unique(bundle.pass_id));
if isfield(Wres, 'Result') && isstruct(Wres.Result)
    row.original_window_eo = getfield_default_local(Wres.Result, 'EO_id', NaN);
    row.original_window_rmse = getfield_default_local( ...
        Wres.Result, 'weighted_voltage_rmse', NaN);
end
end


function T = add_pairwise_winner_columns_local(T, S)
for i = 1:height(T)
    same = T.window_id == T.window_id(i) & ...
        T.bundle_source == T.bundle_source(i) & ...
        T.model == T.model(i) & ...
        ismember(T.EO, S.target_eos);

    P = T(same, :);
    if height(P) < 2
        continue;
    end

    [~, idx_best] = min(P.weighted_rmse);
    T.winner_eo_same_model_source(i) = P.EO(idx_best);

    i4 = find(P.EO == 4, 1);
    i12 = find(P.EO == 12, 1);
    if ~isempty(i4) && ~isempty(i12)
        T.delta_rmse_EO12_minus_EO4(i) = ...
            P.weighted_rmse(i12) - P.weighted_rmse(i4);
    end
end

for iw = unique(T.window_id(:)).'
    for isource = 1:numel(S.bundle_sources)
        source_name = string(S.bundle_sources{isource});
        noeta = T.window_id == iw & T.bundle_source == source_name & ...
            T.model == "NoEta";
        witheta = T.window_id == iw & T.bundle_source == source_name & ...
            T.model == "WithEta";
        if any(noeta) && any(witheta)
            winner_noeta = T.winner_eo_same_model_source(find(noeta, 1));
            winner_witheta = T.winner_eo_same_model_source(find(witheta, 1));
            changed = isfinite(winner_noeta) && isfinite(winner_witheta) && ...
                winner_noeta ~= winner_witheta;
            T.eta_changes_winner(noeta | witheta) = changed;
        end
    end
end
end


function print_eta_effect_summary_local(T, S)
for isource = 1:numel(S.bundle_sources)
    source_name = string(S.bundle_sources{isource});
    fprintf('\n-- Bundle source: %s --\n', source_name);
    for model_name = ["NoEta", "WithEta"]
        TT = T(T.bundle_source == source_name & T.model == model_name, :);
        if isempty(TT)
            continue;
        end
        winners = groupsummary(TT, 'window_id', 'min', 'weighted_rmse');
        win_eo = nan(height(winners), 1);
        for i = 1:height(winners)
            rows = TT(TT.window_id == winners.window_id(i), :);
            [~, ibest] = min(rows.weighted_rmse);
            win_eo(i) = rows.EO(ibest);
        end
        fprintf('%s winners: EO4=%d, EO12=%d, other=%d\n', ...
            model_name, nnz(win_eo == 4), nnz(win_eo == 12), ...
            nnz(~ismember(win_eo, [4 12])));
    end

    for iw = 1:5
        noeta_delta = get_delta_local(T, iw, source_name, "NoEta");
        eta_delta = get_delta_local(T, iw, source_name, "WithEta");
        eta12 = get_eta_max_local(T, iw, source_name, "WithEta", 12);
        fprintf('W%02d delta(EO12-EO4): NoEta=%+.6g, WithEta=%+.6g, eta12_max=%.4g mm\n', ...
            iw, noeta_delta, eta_delta, eta12);
    end
end
end


function delta = get_delta_local(T, iw, source_name, model_name)
idx = T.window_id == iw & T.bundle_source == source_name & T.model == model_name;
if ~any(idx)
    delta = NaN;
else
    delta = T.delta_rmse_EO12_minus_EO4(find(idx, 1));
end
end


function eta_max = get_eta_max_local(T, iw, source_name, model_name, eo)
idx = T.window_id == iw & T.bundle_source == source_name & ...
    T.model == model_name & T.EO == eo;
if ~any(idx)
    eta_max = NaN;
else
    eta_max = T.eta_max_abs_mm(find(idx, 1));
end
end


function plot_eta_effect_summary_local(T, S)
fig = figure('Color', 'w', 'Name', 'Step05 eta effect EO4 vs EO12');
tiledlayout(fig, 2, numel(S.bundle_sources), 'Padding', 'compact', ...
    'TileSpacing', 'compact');

for isource = 1:numel(S.bundle_sources)
    source_name = string(S.bundle_sources{isource});
    nexttile;
    hold on;
    for model_name = ["NoEta", "WithEta"]
        D = unique(T(T.bundle_source == source_name & ...
            T.model == model_name, {'window_id', 'delta_rmse_EO12_minus_EO4'}), ...
            'rows');
        plot(D.window_id, D.delta_rmse_EO12_minus_EO4, '-o', ...
            'DisplayName', model_name);
    end
    yline(0, 'k--', 'DisplayName', 'tie');
    grid on;
    xlabel('window');
    ylabel('RMSE(EO12)-RMSE(EO4)');
    title(sprintf('%s winner sign', source_name), 'Interpreter', 'none');
    legend('Location', 'best');

    nexttile;
    E = T(T.bundle_source == source_name & T.model == "WithEta", :);
    hold on;
    for eo = S.target_eos
        rows = E(E.EO == eo, :);
        plot(rows.window_id, rows.eta_max_abs_mm, '-o', ...
            'DisplayName', sprintf('EO%d', eo));
    end
    grid on;
    xlabel('window');
    ylabel('max |eta_s| (mm)');
    title(sprintf('%s WithEta eta', source_name), 'Interpreter', 'none');
    legend('Location', 'best');
end

png_file = fullfile(S.output_dir, 'EtaEffect_EO4_EO12_20241106.png');
saveas(fig, png_file);
fprintf('Saved diagnostic figure:\n  %s\n', png_file);
end


function y = clamp_scalar_local(x, limit_abs)
if ~isfinite(x)
    y = x;
else
    y = min(max(x, -abs(limit_abs)), abs(limit_abs));
end
end


function y = clamp_vector_local(x, limit_abs)
y = x;
limit_abs = abs(limit_abs(:).');
if numel(limit_abs) == 1
    limit_abs = repmat(limit_abs, size(y));
end
limit_abs = limit_abs(1:min(numel(limit_abs), numel(y)));
if numel(limit_abs) < numel(y)
    limit_abs(end+1:numel(y)) = limit_abs(end);
end
y = min(max(y, -limit_abs), limit_abs);
end


function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end


function value = getfield_default_local(s, name, default_value)
value = default_value;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end
