function Summary = Calibrate_Step05_JointStaticEta_FromPreview(dataset_dir, varargin)
%CALIBRATE_STEP05_JOINTSTATICETA_FROMPREVIEW
% Estimate a shared residual sensor eta from saved Step05 bundle previews.
%
% Default diagnostic calibration linearizes around each window's NoEta
% solution:
%   xq0_i = x_i - dx0_w - A0_w*sin(EO*theta_i + phi0_w)
%   q_i = -(V_i - T(xq0_i)) / T'(xq0_i)
%   q_i ~= ddx_w + da_w*sin(EO*theta_i) + db_w*cos(EO*theta_i) + eta_s
%
% The reference sensor eta is fixed to zero.  dx_w, a_w, and b_w are
% window-specific, while eta_s is shared by all selected windows.

p = inputParser;
p.addRequired('dataset_dir', @(x) ischar(x) || isstring(x));
p.addParameter('NoEtaResultFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('FreeEtaResultFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('Step04TemplateFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('GradientMinRatio', 0.08, @(x) isnumeric(x) && isscalar(x));
p.addParameter('QAbsLimitMm', 1.50, @(x) isnumeric(x) && isscalar(x));
p.addParameter('HuberIterations', 4, @(x) isnumeric(x) && isscalar(x));
p.addParameter('HuberK', 1.50, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinWindows', 3, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinPointsPerSensorWindow', 60, @(x) isnumeric(x) && isscalar(x));
p.addParameter('LinearizationMode', 'about_noeta_solution', ...
    @(x) ischar(x) || isstring(x));
p.addParameter('CalibrationEO', 'reported', @(x) ischar(x) || isstring(x));
p.addParameter('ConsensusPlanPolicy', 'all_valid_median', ...
    @(x) ischar(x) || isstring(x));
p.addParameter('ConsensusIQRWarnMm', 0.15, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(dataset_dir, varargin{:});

opt = p.Results;
dataset_dir = char(opt.dataset_dir);
noeta_file = char(opt.NoEtaResultFile);
free_eta_file = char(opt.FreeEtaResultFile);
template_file = char(opt.Step04TemplateFile);
out_dir = char(opt.OutputDir);

if isempty(out_dir)
    out_dir = fullfile(dataset_dir, 'output', 'step05_joint_static_eta_preview');
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

if ~isfile(noeta_file)
    error('NoEta Step05 result file not found: %s', noeta_file);
end

loaded = load(noeta_file, 'Result');
if ~isfield(loaded, 'Result')
    error('NoEta result file must contain variable Result: %s', noeta_file);
end
Result = loaded.Result;

sensors = double(Result.SensorIDs(:).');
target_blade = double(Result.TargetBlade);
reported_eo = resolve_reported_eo_local(Result);
if ~isfinite(reported_eo)
    error('Cannot resolve reported EO from NoEta result: %s', noeta_file);
end

WindowResult = Result.WindowResult;
plans = build_eta_plan_set_local(WindowResult, reported_eo, opt.MinWindows);
TemplateLookup = build_template_lookup_local(template_file, sensors, target_blade);

PlanSummary = repmat(empty_plan_summary_row_local(sensors), numel(plans), 1);
PlanFits = repmat(struct(), numel(plans), 1);
for ip = 1:numel(plans)
    try
        fit = fit_joint_eta_plan_local(WindowResult, sensors, reported_eo, ...
            plans(ip), TemplateLookup, opt);
        PlanFits(ip).name = plans(ip).name;
        PlanFits(ip).fit = fit;
        PlanSummary(ip) = make_plan_summary_row_local(plans(ip), fit, sensors);
    catch ME
        fit = struct();
        fit.status = 'failed';
        fit.failure_reason = ME.message;
        fit.eta_mm = nan(1, numel(sensors));
        fit.eta_se_mm = nan(1, numel(sensors));
        fit.used_window_ids = plans(ip).window_ids;
        PlanFits(ip).name = plans(ip).name;
        PlanFits(ip).fit = fit;
        PlanSummary(ip) = make_plan_summary_row_local(plans(ip), fit, sensors);
    end
end

PlanTable = struct2table(PlanSummary);

[step04_eta, Step04CenterTable] = load_step04_xc_eta_local( ...
    template_file, sensors, target_blade);
[free_eta_median, free_eta_iqr, free_eta_count] = load_free_eta_summary_local( ...
    free_eta_file, sensors);

valid_plan_mask = strcmpi(PlanTable.status, "ok");
eta_mat = nan(0, numel(sensors));
if any(valid_plan_mask)
    eta_mat = vertcat(PlanSummary(valid_plan_mask).eta_mm);
end

Consensus = struct();
Consensus.status = "no_valid_plan";
Consensus.diagnostic_status = "no_valid_plan";
Consensus.eta_median_mm = nan(1, numel(sensors));
Consensus.eta_plan_iqr_mm = nan(1, numel(sensors));
Consensus.max_abs_plan_iqr_mm = NaN;
Consensus.reference_sensor_id = sensors(1);
Consensus.eta_source_policy = string(opt.ConsensusPlanPolicy);
Consensus.primary_plan_name = "";
if ~isempty(eta_mat)
    Consensus.status = "ok";
    Consensus.eta_plan_iqr_mm = iqr_local(eta_mat, 1);
    Consensus.max_abs_plan_iqr_mm = max(abs(Consensus.eta_plan_iqr_mm), [], 'omitnan');
    if isfinite(opt.ConsensusIQRWarnMm) && ...
            Consensus.max_abs_plan_iqr_mm > opt.ConsensusIQRWarnMm
        Consensus.diagnostic_status = "unstable_plan_spread";
    else
        Consensus.diagnostic_status = "ok";
    end

    policy = lower(strtrim(char(opt.ConsensusPlanPolicy)));
    switch policy
        case {'all_ok', 'primary_all_ok', 'full_coverage'}
            primary_idx = find(valid_plan_mask & strcmpi(PlanTable.plan_name, "all_ok"), 1);
            if isempty(primary_idx)
                primary_idx = find(valid_plan_mask, 1);
            end
            Consensus.eta_median_mm = PlanSummary(primary_idx).eta_mm;
            Consensus.primary_plan_name = string(PlanSummary(primary_idx).plan_name);
        case {'all_valid_median', 'median', 'robust_median'}
            Consensus.eta_median_mm = median(eta_mat, 1, 'omitnan');
            Consensus.primary_plan_name = "all_valid_median";
        otherwise
            error('Unsupported ConsensusPlanPolicy: %s', policy);
    end
end

Reference = table(sensors(:), step04_eta(:), free_eta_median(:), free_eta_iqr(:), ...
    'VariableNames', {'sensor_id', 'step04_xc_eta_mm', ...
    'free_window_eta_median_mm', 'free_window_eta_iqr_mm'});

Summary = struct();
Summary.DatasetDir = dataset_dir;
Summary.NoEtaResultFile = noeta_file;
Summary.FreeEtaResultFile = free_eta_file;
Summary.Step04TemplateFile = template_file;
Summary.OutputDir = out_dir;
Summary.TargetBlade = target_blade;
Summary.SensorIDs = sensors;
Summary.ReferenceSensorID = sensors(1);
Summary.ReportedEO = reported_eo;
Summary.Model = ['q = ddx_w + da_w*sin(EO_w*theta) + db_w*cos(EO_w*theta) + eta_s, ' ...
    'linearized about each NoEta window solution by default'];
Summary.SignConvention = 'Step05 forward model uses x_query = x - dx_c - eta_s - A*sin(EO*theta+phi).';
Summary.Options = rmfield(opt, {'dataset_dir', 'NoEtaResultFile', ...
    'FreeEtaResultFile', 'Step04TemplateFile', 'OutputDir'});
Summary.PlanTable = PlanTable;
Summary.PlanFits = PlanFits;
Summary.Consensus = Consensus;
Summary.ReferenceEtaTable = Reference;
Summary.Step04CenterTable = Step04CenterTable;
Summary.FreeEtaWindowCount = free_eta_count;

summary_mat = fullfile(out_dir, sprintf( ...
    'JointStaticEtaPreview_B%d_S%s.mat', target_blade, sensor_tag_local(sensors)));
summary_csv = fullfile(out_dir, sprintf( ...
    'JointStaticEtaPreview_Plans_B%d_S%s.csv', target_blade, sensor_tag_local(sensors)));
reference_csv = fullfile(out_dir, sprintf( ...
    'JointStaticEtaPreview_Reference_B%d_S%s.csv', target_blade, sensor_tag_local(sensors)));
save(summary_mat, 'Summary');
writetable(PlanTable, summary_csv);
writetable(Reference, reference_csv);

if opt.Verbose
    fprintf('\n=== Joint static residual eta from Step05 previews ===\n');
    fprintf('Dataset dir: %s\n', dataset_dir);
fprintf('Target blade: B%d; sensors: %s; reference CH%d; EO=%g\n', ...
        target_blade, mat2str(sensors), sensors(1), reported_eo);
    fprintf('Calibration EO policy: %s\n', char(opt.CalibrationEO));
    fprintf('NoEta result: %s\n', noeta_file);
    if isfile(free_eta_file)
        fprintf('FreeEta result: %s\n', free_eta_file);
    end
    if isfile(template_file)
        fprintf('Step04 template: %s\n', template_file);
    end
    disp(PlanTable);
    fprintf('Consensus eta median across valid plans: %s mm\n', ...
        mat2str(Consensus.eta_median_mm, 8));
    fprintf('Consensus policy: %s; primary plan: %s; diagnostic status: %s\n', ...
        char(Consensus.eta_source_policy), char(Consensus.primary_plan_name), ...
        char(Consensus.diagnostic_status));
    fprintf('Consensus eta plan IQR: %s mm\n', ...
        mat2str(Consensus.eta_plan_iqr_mm, 8));
    fprintf('Reference comparison:\n');
    disp(Reference);
    fprintf('Saved: %s\n', summary_mat);
end
end


function reported_eo = resolve_reported_eo_local(Result)
reported_eo = NaN;
if isfield(Result, 'ResonanceSummary') && isstruct(Result.ResonanceSummary)
    S = Result.ResonanceSummary;
    if isfield(S, 'reported_eo') && isfinite_scalar_local(S.reported_eo)
        reported_eo = double(S.reported_eo);
        return;
    end
    if isfield(S, 'joint_best_eo') && isfinite_scalar_local(S.joint_best_eo)
        reported_eo = double(S.joint_best_eo);
        return;
    end
    if isfield(S, 'dominant_eo') && isfinite_scalar_local(S.dominant_eo)
        reported_eo = double(S.dominant_eo);
        return;
    end
end
if isfield(Result, 'Trend') && istable(Result.Trend) && ...
        any(strcmpi(Result.Trend.Properties.VariableNames, 'EO_id'))
    eo = Result.Trend.EO_id(isfinite(Result.Trend.EO_id));
    if ~isempty(eo)
        reported_eo = mode(round(eo));
    end
end
end


function tf = isfinite_scalar_local(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x);
end


function plans = build_eta_plan_set_local(WindowResult, reported_eo, min_windows)
window_ids = arrayfun(@(w) get_scalar_field_local(w, 'window_id', NaN), WindowResult);
if any(~isfinite(window_ids))
    window_ids = 1:numel(WindowResult);
end

ok = false(size(WindowResult));
eo = nan(size(WindowResult));
for i = 1:numel(WindowResult)
    if isfield(WindowResult(i), 'Result') && isstruct(WindowResult(i).Result)
        r = WindowResult(i).Result;
        ok(i) = isfield(r, 'status') && strcmpi(char(r.status), 'ok');
        eo(i) = get_scalar_field_local(r, 'EO_id', NaN);
    end
end

ok_ids = window_ids(ok);
reported_ids = window_ids(ok & round(eo) == round(reported_eo));
odd_ids = ok_ids(mod(1:numel(ok_ids), 2) == 1);
even_ids = ok_ids(mod(1:numel(ok_ids), 2) == 0);

plans = struct('name', {}, 'description', {}, 'window_ids', {});
plans(end+1) = make_plan_local('all_ok', ...
    'All successful NoEta windows, calibrated at reported joint EO.', ok_ids);
plans(end+1) = make_plan_local('reported_eo_windows', ...
    'Only windows whose NoEta final EO already equals the reported joint EO.', reported_ids);
plans(end+1) = make_plan_local('first6_ok', ...
    'First six successful windows, calibrated at reported joint EO.', first_n_local(ok_ids, 6));
plans(end+1) = make_plan_local('last6_ok', ...
    'Last six successful windows, calibrated at reported joint EO.', last_n_local(ok_ids, 6));
plans(end+1) = make_plan_local('first6_reported_eo', ...
    'First six successful windows whose NoEta final EO equals reported joint EO.', first_n_local(reported_ids, 6));
plans(end+1) = make_plan_local('last6_reported_eo', ...
    'Last six successful windows whose NoEta final EO equals reported joint EO.', last_n_local(reported_ids, 6));
plans(end+1) = make_plan_local('odd_ok', ...
    'Odd-index split of successful windows.', odd_ids);
plans(end+1) = make_plan_local('even_ok', ...
    'Even-index split of successful windows.', even_ids);

keep = arrayfun(@(p) numel(p.window_ids) >= min_windows, plans);
plans = plans(keep);
end


function plan = make_plan_local(name, description, window_ids)
plan = struct();
plan.name = char(name);
plan.description = char(description);
plan.window_ids = double(window_ids(:).');
end


function y = first_n_local(x, n)
x = x(:).';
y = x(1:min(n, numel(x)));
end


function y = last_n_local(x, n)
x = x(:).';
y = x(max(1, numel(x) - n + 1):end);
end


function fit = fit_joint_eta_plan_local(WindowResult, sensors, EO, plan, ...
    TemplateLookup, opt)
n_sensors = numel(sensors);

q_all = [];
theta_all = [];
sensor_idx_all = [];
window_local_all = [];
tprime_all = [];
base_weight_all = [];
design_eo_all = [];
window_rows = repmat(struct( ...
    'window_id', NaN, ...
    'point_count', 0, ...
    'sensor_count', 0, ...
    'q_median_mm', NaN, ...
    'q_iqr_mm', NaN, ...
    'base_eo_id', NaN, ...
    'linearization_note', ""), numel(plan.window_ids), 1);

window_ids_all = arrayfun(@(w) get_scalar_field_local(w, 'window_id', NaN), WindowResult);
if any(~isfinite(window_ids_all))
    window_ids_all = 1:numel(WindowResult);
end

used_window_ids = [];
for iw = 1:numel(plan.window_ids)
    window_id = plan.window_ids(iw);
    widx = find(window_ids_all == window_id, 1, 'first');
    if isempty(widx)
        continue;
    end
    C = get_window_bundle_preview_local(WindowResult(widx));
    if isempty(C)
        continue;
    end

    [rows, diag_row] = extract_linearized_rows_local( ...
        C, WindowResult(widx).Result, sensors, EO, TemplateLookup, opt);
    diag_row.window_id = window_id;
    window_rows(iw) = diag_row;
    if isempty(rows.q)
        continue;
    end
    if diag_row.sensor_count < 2
        continue;
    end

    used_window_ids(end+1) = window_id; %#ok<AGROW>
    local_id = numel(used_window_ids);
    q_all = [q_all; rows.q(:)]; %#ok<AGROW>
    theta_all = [theta_all; rows.theta(:)]; %#ok<AGROW>
    sensor_idx_all = [sensor_idx_all; rows.sensor_idx(:)]; %#ok<AGROW>
    window_local_all = [window_local_all; local_id .* ones(numel(rows.q), 1)]; %#ok<AGROW>
    tprime_all = [tprime_all; rows.tprime(:)]; %#ok<AGROW>
    base_weight_all = [base_weight_all; rows.base_weight(:)]; %#ok<AGROW>
    design_eo_all = [design_eo_all; rows.design_eo(:)]; %#ok<AGROW>
end

if numel(used_window_ids) < opt.MinWindows
    error('Only %d usable windows after preview filtering; minimum is %d.', ...
        numel(used_window_ids), opt.MinWindows);
end
if isempty(q_all)
    error('No usable q rows were extracted from selected previews.');
end

n_windows = numel(used_window_ids);
n_rows = numel(q_all);
n_eta = max(n_sensors - 1, 0);
n_params = n_eta + 3 * n_windows;

row_idx = [];
col_idx = [];
val = [];

eta_col = sensor_idx_all - 1;
eta_mask = sensor_idx_all > 1;
row_idx = [row_idx; find(eta_mask)]; %#ok<AGROW>
col_idx = [col_idx; eta_col(eta_mask)]; %#ok<AGROW>
val = [val; ones(nnz(eta_mask), 1)]; %#ok<AGROW>

base_col = n_eta + (window_local_all - 1) * 3;
row = (1:n_rows).';
row_idx = [row_idx; row; row; row]; %#ok<AGROW>
col_idx = [col_idx; base_col + 1; base_col + 2; base_col + 3]; %#ok<AGROW>
val = [val; ones(n_rows, 1); sin(design_eo_all .* theta_all); ...
    cos(design_eo_all .* theta_all)]; %#ok<AGROW>

X = sparse(row_idx, col_idx, val, n_rows, n_params);

base_weight_all(~isfinite(base_weight_all) | base_weight_all <= 0) = 0;
tprime_all(~isfinite(tprime_all)) = 0;
fit_weight = base_weight_all .* (tprime_all .^ 2);
positive = fit_weight > 0 & isfinite(fit_weight);
if ~any(positive)
    fit_weight = ones(n_rows, 1);
else
    med_w = median(fit_weight(positive));
    fit_weight = fit_weight ./ max(med_w, eps);
    fit_weight(~isfinite(fit_weight) | fit_weight <= 0) = min(fit_weight(positive), [], 'omitnan');
end

robust_weight = ones(n_rows, 1);
beta = nan(n_params, 1);
for iter = 1:max(1, round(opt.HuberIterations))
    sw = sqrt(max(fit_weight .* robust_weight, eps));
    beta = (X .* sw) \ (q_all .* sw);
    rq = q_all - X * beta;
    rv = tprime_all .* rq;
    sigma = robust_sigma_local(rv);
    if sigma <= 0 || ~isfinite(sigma)
        break;
    end
    cutoff = opt.HuberK * sigma;
    robust_weight = min(1, cutoff ./ max(abs(rv), eps));
end

rq = q_all - X * beta;
rv = tprime_all .* rq;
final_weight = max(fit_weight .* robust_weight, eps);
base_weight_report = max(base_weight_all, eps);

eta = zeros(1, n_sensors);
eta_se = nan(1, n_sensors);
if n_eta > 0
    eta(2:end) = beta(1:n_eta).';
end

AtA = full(X' * spdiags(final_weight, 0, n_rows, n_rows) * X);
rc = rcond(AtA);
if isfinite(rc) && rc > 1e-12 && n_rows > n_params
    mse = sum(final_weight .* (rq .^ 2)) / max(n_rows - n_params, 1);
    cov_beta = (AtA \ eye(size(AtA))) .* mse;
    if n_eta > 0
        eta_se(2:end) = sqrt(max(diag(cov_beta(1:n_eta, 1:n_eta)), 0)).';
    end
    eta_se(1) = 0;
end

WindowParam = table();
dx = nan(n_windows, 1);
a = nan(n_windows, 1);
b = nan(n_windows, 1);
amp = nan(n_windows, 1);
phi = nan(n_windows, 1);
for iw = 1:n_windows
    j = n_eta + (iw - 1) * 3;
    dx(iw) = beta(j + 1);
    a(iw) = beta(j + 2);
    b(iw) = beta(j + 3);
    amp(iw) = hypot(a(iw), b(iw));
    phi(iw) = atan2(b(iw), a(iw));
end
WindowParam.window_id = used_window_ids(:);
WindowParam.dx_mm = dx;
WindowParam.a_mm = a;
WindowParam.b_mm = b;
WindowParam.A_mm = amp;
WindowParam.phi_rad = phi;

fit = struct();
fit.status = 'ok';
fit.failure_reason = '';
fit.plan_name = plan.name;
fit.description = plan.description;
fit.reported_eo = EO;
fit.linearization_mode = char(opt.LinearizationMode);
fit.used_window_ids = used_window_ids;
fit.window_count = n_windows;
fit.point_count = n_rows;
fit.sensor_ids = sensors;
fit.eta_mm = eta;
fit.eta_se_mm = eta_se;
fit.max_abs_eta_mm = max(abs(eta), [], 'omitnan');
fit.weighted_q_rmse_mm = sqrt(sum(final_weight .* rq .^ 2) / sum(final_weight));
fit.weighted_linear_voltage_rmse_v = sqrt( ...
    sum(base_weight_report .* rv .^ 2) / sum(base_weight_report));
fit.robust_downweight_fraction = mean(robust_weight < 0.999);
fit.normal_matrix_rcond = rc;
fit.WindowParam = WindowParam;
fit.WindowRows = struct2table(window_rows);
end


function C = get_window_bundle_preview_local(W)
C = [];
names = {'CoreBundlePreview', 'BundlePreview', 'ExpandedBundlePreview'};
for i = 1:numel(names)
    if isfield(W, names{i}) && isstruct(W.(names{i})) && ...
            isfield(W.(names{i}), 'x') && ~isempty(W.(names{i}).x)
        C = W.(names{i});
        return;
    end
end
end


function [rows, diag_row] = extract_linearized_rows_local(C, R, sensors, EO, ...
    TemplateLookup, opt)
rows = struct('q', [], 'theta', [], 'sensor_idx', [], ...
    'tprime', [], 'base_weight', [], 'design_eo', []);
diag_row = struct('window_id', NaN, 'point_count', 0, ...
    'sensor_count', 0, 'q_median_mm', NaN, 'q_iqr_mm', NaN, ...
    'base_eo_id', NaN, 'linearization_note', "");

needed = {'v', 'theta', 'sensor_id', 'template_v', 'template_dv_dx'};
for i = 1:numel(needed)
    if ~isfield(C, needed{i})
        return;
    end
end

v = C.v(:);
theta = C.theta(:);
sensor_id = double(C.sensor_id(:));
n = numel(v);
if any([numel(theta), numel(sensor_id)] ~= n)
    return;
end

[tv, tp, lin_ok, lin_note, base_eo, design_eo] = evaluate_linearization_template_local( ...
    C, R, sensors, EO, TemplateLookup, opt);
diag_row.base_eo_id = base_eo;
diag_row.linearization_note = string(lin_note);
if ~lin_ok
    return;
end
if any([numel(tv), numel(tp)] ~= n)
    return;
end

mask = isfinite(v) & isfinite(theta) & isfinite(sensor_id) & ...
    isfinite(tv) & isfinite(tp) & tp ~= 0;
mask = mask & same_length_mask_local(C, 'finite_mask', n);
mask = mask & same_length_mask_local(C, 'core_mask', n);
mask = mask & same_length_mask_local(C, 'inside_domain_mask', n);
mask = mask & same_length_mask_local(C, 'inside_guard_mask', n);
mask = mask & same_length_mask_local(C, 'main_pulse_mask', n);
mask = mask & same_length_mask_local(C, 'dynamic_effective_mask', n);

grad_keep = false(n, 1);
sensor_idx = nan(n, 1);
for is = 1:numel(sensors)
    sid = sensors(is);
    ms = mask & sensor_id == sid;
    sensor_idx(sensor_id == sid) = is;
    if ~any(ms)
        continue;
    end
    abs_tp = abs(tp(ms));
    gref = quantile_simple_local(abs_tp, 0.95);
    if ~isfinite(gref) || gref <= 0
        gref = max(abs_tp, [], 'omitnan');
    end
    grad_keep(ms) = abs(tp(ms)) >= max(opt.GradientMinRatio * gref, eps);
end

mask = mask & grad_keep & isfinite(sensor_idx);

q = -(v - tv) ./ tp;
mask = mask & isfinite(q) & abs(q) <= opt.QAbsLimitMm;

base_weight = ones(n, 1);
if isfield(C, 'fit_weight') && numel(C.fit_weight) == n
    base_weight = double(C.fit_weight(:));
elseif isfield(C, 'template_weight') && numel(C.template_weight) == n
    base_weight = double(C.template_weight(:));
elseif isfield(C, 'legacy_fit_weight') && numel(C.legacy_fit_weight) == n
    base_weight = double(C.legacy_fit_weight(:));
end
base_weight(~isfinite(base_weight) | base_weight <= 0) = 0;
mask = mask & base_weight > 0;

keep_by_sensor = false(n, 1);
for is = 1:numel(sensors)
    ms = mask & sensor_idx == is;
    if nnz(ms) >= opt.MinPointsPerSensorWindow
        keep_by_sensor(ms) = true;
    end
end
mask = mask & keep_by_sensor;

diag_row.point_count = nnz(mask);
diag_row.sensor_count = numel(unique(sensor_idx(mask)));
if any(mask)
    diag_row.q_median_mm = median(q(mask), 'omitnan');
    diag_row.q_iqr_mm = iqr_local(q(mask), 1);
end

rows.q = q(mask);
rows.theta = theta(mask);
rows.sensor_idx = sensor_idx(mask);
rows.tprime = tp(mask);
rows.base_weight = base_weight(mask);
rows.design_eo = design_eo .* ones(nnz(mask), 1);
end


function [tv, tp, ok, note, base_eo, design_eo] = evaluate_linearization_template_local( ...
    C, R, sensors, EO, TemplateLookup, opt)
tv = [];
tp = [];
ok = false;
note = '';
base_eo = NaN;
design_eo = EO;
mode = lower(strtrim(char(opt.LinearizationMode)));
eo_policy = lower(strtrim(char(opt.CalibrationEO)));

switch mode
    case {'zero', 'raw_x', 'about_raw_x'}
        if ~isfield(C, 'template_v') || ~isfield(C, 'template_dv_dx')
            note = 'missing raw-x template fields';
            return;
        end
        tv = C.template_v(:);
        tp = C.template_dv_dx(:);
        ok = true;
        note = 'linearized about raw mapped x';

    case {'about_noeta_solution', 'noeta', 'residual_noeta'}
        if isempty(TemplateLookup) || ~isfield(C, 'x')
            note = 'missing template lookup or x field';
            return;
        end
        base_eo = get_scalar_field_local(R, 'EO_id', NaN);
        if ~isfinite(base_eo)
            note = 'missing NoEta base EO';
            return;
        end
        if any(strcmpi(eo_policy, {'window_base', 'base_window', ...
                'per_window_base', 'mixed_base'}))
            design_eo = base_eo;
        else
            design_eo = EO;
            if round(base_eo) ~= round(EO)
                note = sprintf('skipped: base EO %.6g differs from reported EO %.6g', ...
                    base_eo, EO);
                return;
            end
        end
        A0 = get_scalar_field_local(R, 'A_id', NaN);
        phi0 = get_scalar_field_local(R, 'phi_id_wrapped', NaN);
        dx0 = get_scalar_field_local(R, 'dx_c_id', NaN);
        if any(~isfinite([A0, phi0, dx0]))
            note = 'missing NoEta base A/phi/dx';
            return;
        end
        x = C.x(:);
        theta = C.theta(:);
        sensor_id = double(C.sensor_id(:));
        if any([numel(theta), numel(sensor_id)] ~= numel(x))
            note = 'x/theta/sensor length mismatch';
            return;
        end
        xq0 = x - dx0 - A0 .* sin(design_eo .* theta + phi0);
        [tv, tp] = evaluate_template_lookup_local(TemplateLookup, sensors, ...
            sensor_id, xq0);
        ok = true;
        note = sprintf('linearized about NoEta solution, design EO %.6g', design_eo);

    otherwise
        error('Unsupported LinearizationMode: %s', mode);
end
end


function [tv, tp] = evaluate_template_lookup_local(TemplateLookup, sensors, ...
    sensor_id, xq)
tv = nan(size(xq));
tp = nan(size(xq));
for is = 1:numel(sensors)
    mask = sensor_id == sensors(is);
    if ~any(mask) || isempty(TemplateLookup(is).x)
        continue;
    end
    tv(mask) = interp1(TemplateLookup(is).x, TemplateLookup(is).v, ...
        xq(mask), 'linear', NaN);
    tp(mask) = interp1(TemplateLookup(is).x, TemplateLookup(is).dv, ...
        xq(mask), 'linear', NaN);
end
end


function m = same_length_mask_local(C, field_name, n)
m = true(n, 1);
if isfield(C, field_name) && numel(C.(field_name)) == n
    x = C.(field_name)(:);
    if islogical(x)
        m = x;
    else
        m = isfinite(double(x)) & double(x) ~= 0;
    end
end
end


function x = quantile_simple_local(v, q)
v = sort(v(isfinite(v)));
if isempty(v)
    x = NaN;
    return;
end
q = min(max(q, 0), 1);
pos = 1 + (numel(v) - 1) * q;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    x = v(lo);
else
    x = v(lo) + (v(hi) - v(lo)) * (pos - lo);
end
end


function s = robust_sigma_local(x)
x = x(isfinite(x));
if isempty(x)
    s = NaN;
    return;
end
med = median(x);
s = 1.4826 * median(abs(x - med));
if s <= 0 || ~isfinite(s)
    s = sqrt(mean((x - med) .^ 2, 'omitnan'));
end
end


function row = empty_plan_summary_row_local(sensors)
row = struct();
row.plan_name = "";
row.status = "";
row.failure_reason = "";
row.window_count = NaN;
row.point_count = NaN;
row.used_windows = "";
row.weighted_q_rmse_mm = NaN;
row.weighted_linear_voltage_rmse_v = NaN;
row.robust_downweight_fraction = NaN;
row.normal_matrix_rcond = NaN;
row.eta_mm = nan(1, numel(sensors));
row.eta_se_mm = nan(1, numel(sensors));
for i = 1:numel(sensors)
    row.(sprintf('eta_CH%d_mm', sensors(i))) = NaN;
    row.(sprintf('eta_CH%d_se_mm', sensors(i))) = NaN;
end
end


function row = make_plan_summary_row_local(plan, fit, sensors)
row = empty_plan_summary_row_local(sensors);
row.plan_name = string(plan.name);
row.status = string(get_char_field_local(fit, 'status', 'failed'));
row.failure_reason = string(get_char_field_local(fit, 'failure_reason', ''));
row.used_windows = string(mat2str(get_numeric_field_local(fit, 'used_window_ids', [])));
if strcmpi(row.status, "ok")
    row.window_count = fit.window_count;
    row.point_count = fit.point_count;
    row.weighted_q_rmse_mm = fit.weighted_q_rmse_mm;
    row.weighted_linear_voltage_rmse_v = fit.weighted_linear_voltage_rmse_v;
    row.robust_downweight_fraction = fit.robust_downweight_fraction;
    row.normal_matrix_rcond = fit.normal_matrix_rcond;
    row.eta_mm = fit.eta_mm;
    row.eta_se_mm = fit.eta_se_mm;
    for i = 1:numel(sensors)
        row.(sprintf('eta_CH%d_mm', sensors(i))) = fit.eta_mm(i);
        row.(sprintf('eta_CH%d_se_mm', sensors(i))) = fit.eta_se_mm(i);
    end
end
end


function x = get_scalar_field_local(S, field_name, default_value)
x = default_value;
if isstruct(S) && isfield(S, field_name) && isnumeric(S.(field_name)) && ...
        isscalar(S.(field_name)) && isfinite(S.(field_name))
    x = double(S.(field_name));
end
end


function x = get_numeric_field_local(S, field_name, default_value)
x = default_value;
if isstruct(S) && isfield(S, field_name) && isnumeric(S.(field_name))
    x = S.(field_name);
end
end


function x = get_char_field_local(S, field_name, default_value)
x = default_value;
if isstruct(S) && isfield(S, field_name)
    v = S.(field_name);
    if ischar(v) || isstring(v)
        x = char(v);
    end
end
end


function TemplateLookup = build_template_lookup_local(template_file, sensors, blade_id)
TemplateLookup = repmat(struct('sensor_id', NaN, 'x', [], 'v', [], 'dv', []), ...
    1, numel(sensors));
for i = 1:numel(sensors)
    TemplateLookup(i).sensor_id = sensors(i);
end
if isempty(template_file) || ~isfile(template_file)
    return;
end
loaded = load(template_file, 'Template');
if ~isfield(loaded, 'Template') || ~isfield(loaded.Template, 'SensorBlade')
    return;
end
SB = loaded.Template.SensorBlade(:);
for i = 1:numel(sensors)
    tpl = [];
    for k = 1:numel(SB)
        if isfield(SB(k), 'sensor_id') && isfield(SB(k), 'blade_id') && ...
                double(SB(k).sensor_id) == sensors(i) && ...
                double(SB(k).blade_id) == blade_id
            tpl = SB(k);
            break;
        end
    end
    if isempty(tpl) || ~isfield(tpl, 'x_grid') || ~isfield(tpl, 'v_grid')
        continue;
    end
    x = double(tpl.x_grid(:));
    v = double(tpl.v_grid(:));
    if isfield(tpl, 'dv_dx') && numel(tpl.dv_dx) == numel(x)
        dv = double(tpl.dv_dx(:));
    else
        dv = gradient(v, x);
    end
    good = isfinite(x) & isfinite(v) & isfinite(dv);
    x = x(good);
    v = v(good);
    dv = dv(good);
    [x, order] = sort(x);
    v = v(order);
    dv = dv(order);
    [x, unique_idx] = unique(x, 'stable');
    TemplateLookup(i).x = x;
    TemplateLookup(i).v = v(unique_idx);
    TemplateLookup(i).dv = dv(unique_idx);
end
end


function [eta, CenterTable] = load_step04_xc_eta_local(template_file, sensors, blade_id)
eta = nan(1, numel(sensors));
CenterTable = table(sensors(:), nan(numel(sensors), 1), strings(numel(sensors), 1), ...
    'VariableNames', {'sensor_id', 'template_xc_mm', 'source'});
if isempty(template_file) || ~isfile(template_file)
    return;
end
loaded = load(template_file, 'Template');
if ~isfield(loaded, 'Template') || ~isfield(loaded.Template, 'SensorBlade')
    return;
end
SB = loaded.Template.SensorBlade(:);
xc = nan(1, numel(sensors));
src = strings(1, numel(sensors));
for i = 1:numel(sensors)
    found = [];
    for k = 1:numel(SB)
        if isfield(SB(k), 'sensor_id') && isfield(SB(k), 'blade_id') && ...
                double(SB(k).sensor_id) == sensors(i) && ...
                double(SB(k).blade_id) == blade_id
            found = SB(k);
            break;
        end
    end
    if isempty(found)
        continue;
    end
    if isfield(found, 'xc') && isfinite(found.xc)
        xc(i) = double(found.xc);
    elseif isfield(found, 'xc_mm') && isfinite(found.xc_mm)
        xc(i) = double(found.xc_mm);
    end
    if isfield(found, 'xc_reference_source')
        src(i) = string(found.xc_reference_source);
    end
end
if all(isfinite(xc))
    eta = xc(1) - xc;
    eta(1) = 0;
end
CenterTable.template_xc_mm = xc(:);
CenterTable.source = src(:);
end


function [eta_median, eta_iqr, count] = load_free_eta_summary_local(free_eta_file, sensors)
eta_median = nan(1, numel(sensors));
eta_iqr = nan(1, numel(sensors));
count = 0;
if isempty(free_eta_file) || ~isfile(free_eta_file)
    return;
end
loaded = load(free_eta_file, 'Result');
if ~isfield(loaded, 'Result') || ~isfield(loaded.Result, 'WindowResult')
    return;
end
WR = loaded.Result.WindowResult;
eta_rows = [];
for i = 1:numel(WR)
    if ~isfield(WR(i), 'Result') || ~isstruct(WR(i).Result)
        continue;
    end
    r = WR(i).Result;
    if isfield(r, 'status') && ~strcmpi(char(r.status), 'ok')
        continue;
    end
    if ~isfield(r, 'sensor_eta_id') || isempty(r.sensor_eta_id)
        continue;
    end
    eta = double(r.sensor_eta_id(:).');
    if numel(eta) ~= numel(sensors) || any(~isfinite(eta))
        continue;
    end
    eta_rows = [eta_rows; eta]; %#ok<AGROW>
end
count = size(eta_rows, 1);
if count > 0
    eta_median = median(eta_rows, 1, 'omitnan');
    eta_iqr = iqr_local(eta_rows, 1);
end
end


function y = iqr_local(x, dim)
if nargin < 2
    dim = 1;
end
q75 = quantile_dim_local(x, 0.75, dim);
q25 = quantile_dim_local(x, 0.25, dim);
y = q75 - q25;
end


function qv = quantile_dim_local(x, q, dim)
if dim == 1
    qv = nan(1, size(x, 2));
    for j = 1:size(x, 2)
        qv(j) = quantile_simple_local(x(:, j), q);
    end
else
    qv = nan(size(x, 1), 1);
    for i = 1:size(x, 1)
        qv(i) = quantile_simple_local(x(i, :).', q);
    end
end
end


function tag = sensor_tag_local(sensors)
tag = sprintf('%d', sensors);
end
