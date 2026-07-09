function Summary = Diagnose_Step05_ThetaPhasePolicy_20251222()
%DIAGNOSE_STEP05_THETAPHASEPOLICY_20251222
% Read-only check for the remaining SG-vs-foundation Step05 difference.
%
% SG Step05 evaluates the vibration phase by interpolating between adjacent
% OPR events:
%   theta = 2*pi/epr * ((event_id-1) + local_fraction)
%
% The foundation Step05 bundle currently stores a revolution-anchor phase:
%   theta = interp1(OPR(1:epr:end), 2*pi*rev_id, t)
%
% This diagnostic runs the same direct-template voltage fit on saved
% query-safe core points with both theta definitions.

cfg = BTTProjectConfig_20251222();

result_file = strtrim(getenv('STEP05_DIAG_RESULT_FILE'));
if isempty(result_file)
    result_file = fullfile(cfg.step05_direct_template_output_dir, ...
        cfg.dynamic_cases{1}, ...
        'Result_Step05_FoundationMainPulseAdaptiveFixedJointEta_B1_S123_20251222.mat');
end
if ~isfile(result_file)
    error('Missing Step05 result file: %s', result_file);
end

case_name = cfg.dynamic_cases{1};
step02_file = fullfile(cfg.step02_output_dir, case_name, ...
    'DynamicBTTFeature_20251222.mat');
if ~isfile(step02_file)
    error('Missing Step02 file: %s', step02_file);
end

loaded_result = load(result_file, 'Result');
Result = loaded_result.Result;
loaded_step02 = load(step02_file, 'OPRTable');
OPRTable = loaded_step02.OPRTable;

opr_times = resolve_opr_times_local(OPRTable, cfg);
epr = resolve_epr_local(OPRTable, cfg);

S = Result.AnalysisSettings;
S.analysis_sensors = Result.SensorIDs;
S.weight_floor = get_field_or_default_local(S, 'weight_floor', 0.05);
S.template_interp_method = get_field_or_default_local(S, ...
    'template_interp_method', 'pchip');
S.freq_search_hz = get_field_or_default_local(S, 'freq_search_hz', [100 1000]);
S.eo_pad = get_field_or_default_local(S, 'eo_pad', 2);
S.fminsearch_max_iter = 300;
S.fminsearch_max_fun = 900;

out_dir = fullfile(cfg.output_root, 'diagnostics', ...
    'step05_theta_phase_policy_20251222');
if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

all_rows = table();
summary_rows = repmat(make_empty_summary_row_local(), ...
    numel(Result.WindowResult), 1);

for iw = 1:numel(Result.WindowResult)
    WR = Result.WindowResult(iw);
    B = choose_core_bundle_local(WR);
    if isempty(B) || ~isstruct(B) || ~isfield(B, 'x') || isempty(B.x)
        summary_rows(iw).window_id = iw;
        summary_rows(iw).status = "missing_core_bundle";
        continue;
    end

    eo_candidates = build_eo_candidates_local(B.rot_freq_mean_hz, ...
        S.freq_search_hz, S.eo_pad);

    theta_stored = B.theta(:);
    theta_event = map_time_to_adjacent_opr_phase_local( ...
        opr_times, B.t(:), epr);

    Tstored = run_direct_template_scan_local(B, theta_stored, ...
        eo_candidates, S);
    Tstored.theta_policy = repmat("foundation_rev_anchor", height(Tstored), 1);
    Tevent = run_direct_template_scan_local(B, theta_event, ...
        eo_candidates, S);
    Tevent.theta_policy = repmat("adjacent_opr_event_sg_like", height(Tevent), 1);

    T = [Tstored; Tevent];
    T.window_id = repmat(iw, height(T), 1);
    T = movevars(T, {'window_id', 'theta_policy'}, 'Before', 1);
    all_rows = [all_rows; T]; %#ok<AGROW>

    summary_rows(iw) = summarize_window_local(iw, WR, B, ...
        theta_stored, theta_event, T);
end

CandidateTable = all_rows;
WindowSummary = struct2table(summary_rows);

writetable(CandidateTable, fullfile(out_dir, ...
    'Step05_ThetaPolicy_Candidates_B1_S123_20251222.csv'));
writetable(WindowSummary, fullfile(out_dir, ...
    'Step05_ThetaPolicy_WindowSummary_B1_S123_20251222.csv'));

Summary = struct();
Summary.ResultFile = result_file;
Summary.OutputDir = out_dir;
Summary.CandidateTable = CandidateTable;
Summary.WindowSummary = WindowSummary;

fprintf('\n=== Step05 theta-phase policy diagnostic (20251222) ===\n');
fprintf('Source Step05 result:\n  %s\n', result_file);
fprintf('Output dir:\n  %s\n\n', out_dir);
disp(WindowSummary);
print_findings_local(WindowSummary);
end


function B = choose_core_bundle_local(WR)
B = [];
if isfield(WR, 'CoreBundlePreview') && ~isempty(WR.CoreBundlePreview)
    B = WR.CoreBundlePreview;
elseif isfield(WR, 'BundlePreview') && ~isempty(WR.BundlePreview)
    B = WR.BundlePreview;
elseif isfield(WR, 'CandidateBundlePreview') && ~isempty(WR.CandidateBundlePreview)
    B = WR.CandidateBundlePreview;
end
end


function T = run_direct_template_scan_local(B, theta, eo_candidates, S)
B.theta = theta(:);

rows = repmat(struct( ...
    'EO', NaN, ...
    'rank', NaN, ...
    'A_mm', NaN, ...
    'phi_rad', NaN, ...
    'dx_c_mm', NaN, ...
    'weighted_voltage_rmse_v', inf, ...
    'plain_voltage_rmse_v', inf, ...
    'valid_query_fraction', NaN, ...
    'point_count', numel(B.x)), numel(eo_candidates), 1);

for i = 1:numel(eo_candidates)
    EO = eo_candidates(i);
    seed = seed_from_linearized_template_local(B, EO, S);
    p0 = [seed.A, seed.phi, seed.dx_c];
    obj = @(p) direct_template_objective_local(B, EO, p, S);
    opts = optimset('Display', 'off', ...
        'MaxIter', S.fminsearch_max_iter, ...
        'MaxFunEvals', S.fminsearch_max_fun, ...
        'TolX', 1e-7, ...
        'TolFun', 1e-9);
    try
        [p_best, ~] = fminsearch(obj, p0, opts);
    catch
        p_best = p0;
    end

    p_best(1) = min(abs(p_best(1)), 0.5);
    p_best(2) = wrap_to_pi_local(p_best(2));
    p_best(3) = min(max(p_best(3), -0.35), 0.35);

    [wrmse, prmse, valid_fraction] = evaluate_direct_template_local( ...
        B, EO, p_best, S);

    rows(i).EO = EO;
    rows(i).A_mm = p_best(1);
    rows(i).phi_rad = p_best(2);
    rows(i).dx_c_mm = p_best(3);
    rows(i).weighted_voltage_rmse_v = wrmse;
    rows(i).plain_voltage_rmse_v = prmse;
    rows(i).valid_query_fraction = valid_fraction;
end

T = struct2table(rows);
T = T(isfinite(T.EO) & isfinite(T.weighted_voltage_rmse_v), :);
if isempty(T)
    return;
end
T = sortrows(T, {'weighted_voltage_rmse_v', 'EO'}, ...
    {'ascend', 'ascend'});
T.rank = (1:height(T)).';
end


function seed = seed_from_linearized_template_local(B, EO, S)
v = B.v(:);
theta = B.theta(:);
T0 = B.template_v(:);
Tp = B.template_dv_dx(:);
w = B.fit_weight(:);
w(~isfinite(w)) = S.weight_floor;
w = max(w, S.weight_floor);

valid = isfinite(v) & isfinite(theta) & isfinite(T0) & ...
    isfinite(Tp) & isfinite(w) & abs(Tp) > eps;

seed = struct('A', 0.05, 'phi', 0, 'dx_c', 0);
if nnz(valid) < 20
    return;
end

s = sin(EO .* theta(valid));
c = cos(EO .* theta(valid));
X = [-Tp(valid), -Tp(valid).*s, -Tp(valid).*c];
y = v(valid) - T0(valid);
wg = sqrt(w(valid));

beta = (X .* wg) \ (y .* wg);
if numel(beta) < 3 || any(~isfinite(beta))
    return;
end

seed.dx_c = min(max(beta(1), -0.35), 0.35);
seed.A = min(abs(hypot(beta(2), beta(3))), 0.5);
seed.phi = wrap_to_pi_local(atan2(beta(3), beta(2)));
end


function value = direct_template_objective_local(B, EO, p, S)
[wrmse, ~, valid_fraction] = evaluate_direct_template_local(B, EO, p, S);
if ~isfinite(wrmse)
    value = 1e9;
else
    value = wrmse.^2 * max(numel(B.x), 1) + ...
        5 * (1 - valid_fraction) * max(numel(B.x), 1) * wrmse.^2;
end
end


function [wrmse, prmse, valid_fraction] = evaluate_direct_template_local(B, EO, p, S)
A = min(abs(p(1)), 0.5);
phi = wrap_to_pi_local(p(2));
dx_c = min(max(p(3), -0.35), 0.35);

y = A .* sin(EO .* B.theta(:) + phi);
xq = B.x(:) - dx_c - y(:);
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

w = B.fit_weight(:);
w(~isfinite(w)) = S.weight_floor;
w = max(w, S.weight_floor);
valid = isfinite(Vpred) & isfinite(B.v(:)) & isfinite(w) & ...
    isfinite(B.theta(:));
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


function theta_rot = map_time_to_adjacent_opr_phase_local(opr_times, sample_times, epr)
sample_times = sample_times(:);
opr_times = opr_times(:);
theta_rot = nan(size(sample_times));
if numel(opr_times) < 2
    return;
end

idx = discretize(sample_times, [-inf; opr_times; inf]) - 1;
valid = idx >= 1 & idx < numel(opr_times);
if ~any(valid)
    return;
end

frac = (sample_times(valid) - opr_times(idx(valid))) ./ ...
    max(opr_times(idx(valid) + 1) - opr_times(idx(valid)), eps);
theta_rot(valid) = (2*pi/epr) .* ((idx(valid) - 1) + frac);
end


function opr_times = resolve_opr_times_local(OPRTable, cfg)
if ismember('opr_reference_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_reference_time_s;
elseif isfield(cfg, 'opr_reference_time_column') && ...
        ismember(cfg.opr_reference_time_column, OPRTable.Properties.VariableNames)
    opr_times = OPRTable.(cfg.opr_reference_time_column);
elseif ismember('opr_start_time_s', OPRTable.Properties.VariableNames)
    opr_times = OPRTable.opr_start_time_s;
else
    opr_times = OPRTable.opr_center_time_s;
end
opr_times = opr_times(isfinite(opr_times));
end


function epr = resolve_epr_local(OPRTable, cfg)
epr = NaN;
if ismember('opr_events_per_revolution', OPRTable.Properties.VariableNames)
    v = OPRTable.opr_events_per_revolution;
    v = v(isfinite(v));
    if ~isempty(v)
        epr = mode(v);
    end
end
if ~isfinite(epr) && isfield(cfg, 'step04_opr_events_per_revolution') && ...
        isnumeric(cfg.step04_opr_events_per_revolution)
    epr = cfg.step04_opr_events_per_revolution;
end
if ~isfinite(epr)
    epr = cfg.blades_num;
end
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


function row = summarize_window_local(iw, WR, B, theta_stored, theta_event, T)
row = make_empty_summary_row_local();
row.window_id = iw;
row.status = "ok";
if isfield(WR, 'Result') && isstruct(WR.Result)
    row.formal_eo = get_field_or_default_local(WR.Result, 'EO_id', NaN);
    row.formal_rmse_v = get_field_or_default_local(WR.Result, ...
        'weighted_voltage_rmse', NaN);
end

dtheta = wrap_to_pi_local(theta_event(:) - theta_stored(:));
row.theta_diff_rms_rad = sqrt(mean(dtheta(isfinite(dtheta)).^2));
row.theta_diff_p95_rad = prctile(abs(dtheta(isfinite(dtheta))), 95);
row.point_count = numel(B.x);

policies = ["foundation_rev_anchor", "adjacent_opr_event_sg_like"];
for policy = policies
    M = T(strcmpi(T.theta_policy, policy), :);
    if isempty(M)
        continue;
    end
    top_list = strjoin(compose('%d', M.EO(1:min(5, height(M))).'), "|");
    idx14 = find(M.EO == 14, 1, 'first');
    eo14_rank = NaN;
    eo14_rmse = NaN;
    if ~isempty(idx14)
        eo14_rank = M.rank(idx14);
        eo14_rmse = M.weighted_voltage_rmse_v(idx14);
    end
    if policy == "foundation_rev_anchor"
        row.rev_anchor_eo = M.EO(1);
        row.rev_anchor_rmse_v = M.weighted_voltage_rmse_v(1);
        row.rev_anchor_eo14_rank = eo14_rank;
        row.rev_anchor_eo14_rmse_v = eo14_rmse;
        row.rev_anchor_top_list = string(top_list);
    else
        row.event_phase_eo = M.EO(1);
        row.event_phase_rmse_v = M.weighted_voltage_rmse_v(1);
        row.event_phase_eo14_rank = eo14_rank;
        row.event_phase_eo14_rmse_v = eo14_rmse;
        row.event_phase_top_list = string(top_list);
    end
end
end


function row = make_empty_summary_row_local()
row = struct( ...
    'window_id', NaN, ...
    'status', "", ...
    'formal_eo', NaN, ...
    'formal_rmse_v', NaN, ...
    'point_count', NaN, ...
    'theta_diff_rms_rad', NaN, ...
    'theta_diff_p95_rad', NaN, ...
    'rev_anchor_eo', NaN, ...
    'rev_anchor_rmse_v', NaN, ...
    'rev_anchor_eo14_rank', NaN, ...
    'rev_anchor_eo14_rmse_v', NaN, ...
    'rev_anchor_top_list', "", ...
    'event_phase_eo', NaN, ...
    'event_phase_rmse_v', NaN, ...
    'event_phase_eo14_rank', NaN, ...
    'event_phase_eo14_rmse_v', NaN, ...
    'event_phase_top_list', "");
end


function eo_candidates = build_eo_candidates_local(rot_freq_mean_hz, freq_search_hz, eo_pad)
freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min = max(1, ceil(freq_lo / rot_freq_mean_hz) - max(0, eo_pad));
eo_max = max(eo_min, floor(freq_hi / rot_freq_mean_hz) + max(0, eo_pad));
eo_candidates = eo_min:eo_max;
eo_candidates = eo_candidates(eo_candidates > 0);
end


function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end


function value = get_field_or_default_local(S, field_name, default_value)
if isstruct(S) && isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function print_findings_local(WindowSummary)
ok = strcmpi(WindowSummary.status, "ok");
if ~any(ok)
    fprintf('No valid windows.\n');
    return;
end
fprintf('\nEO sequences:\n');
fprintf('  formal saved result:        %s\n', ...
    strjoin(compose('%d', WindowSummary.formal_eo(ok).'), "|"));
fprintf('  foundation rev-anchor scan: %s\n', ...
    strjoin(compose('%d', WindowSummary.rev_anchor_eo(ok).'), "|"));
fprintf('  adjacent OPR event scan:    %s\n', ...
    strjoin(compose('%d', WindowSummary.event_phase_eo(ok).'), "|"));
fprintf('\nMedian |theta_event - theta_rev| p95 = %.4f rad.\n', ...
    median(WindowSummary.theta_diff_p95_rad(ok), 'omitnan'));
fprintf('Median EO14 rank: rev-anchor %.2f, adjacent-event %.2f.\n', ...
    median(WindowSummary.rev_anchor_eo14_rank(ok), 'omitnan'), ...
    median(WindowSummary.event_phase_eo14_rank(ok), 'omitnan'));
end
