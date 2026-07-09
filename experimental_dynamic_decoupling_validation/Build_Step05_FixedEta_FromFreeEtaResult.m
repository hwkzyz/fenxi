function Summary = Build_Step05_FixedEta_FromFreeEtaResult(result_file, output_dir, varargin)
%BUILD_STEP05_FIXEDETA_FROMFREEETARESULT
% Convert a per-window FreeEta Step05 result into a fixed-eta preview file.
%
% The output MAT contains variable Summary and is compatible with formal
% Step05 fixed_joint_static_eta loaders:
%   Summary.Consensus.eta_median_mm
%   Summary.TargetBlade
%   Summary.SensorIDs

p = inputParser;
p.addRequired('result_file', @(x) ischar(x) || isstring(x));
p.addRequired('output_dir', @(x) ischar(x) || isstring(x));
p.addParameter('DatasetDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('Step04TemplateFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('SourceScript', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConsensusPolicy', 'free_eta_all_ok_median', ...
    @(x) ischar(x) || isstring(x));
p.parse(result_file, output_dir, varargin{:});
opt = p.Results;

result_file = char(result_file);
output_dir = char(output_dir);
if ~isfile(result_file)
    error('FreeEta result file not found: %s', result_file);
end
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

loaded = load(result_file, 'Result');
if ~isfield(loaded, 'Result')
    error('FreeEta result file must contain Result: %s', result_file);
end
Result = loaded.Result;

sensors = double(Result.SensorIDs(:).');
target_blade = double(Result.TargetBlade);
sensor_tag = ['S', sprintf('%d', sensors)];

T = Result.Trend;
ok = strcmpi(string(T.status), "ok");
eta = collect_eta_matrix_local(Result, T, sensors);
ok = ok & all(isfinite(eta), 2);

all_plan = summarize_plan_local("all_ok", T, eta, ok, sensors);

eo_ok = round(T.EO_id(ok));
eo_ok = eo_ok(isfinite(eo_ok));
dominant_eo = NaN;
if ~isempty(eo_ok)
    dominant_eo = mode(eo_ok);
end
dominant_mask = ok & round(T.EO_id) == dominant_eo;
dominant_plan = summarize_plan_local("dominant_eo_ok", ...
    T, eta, dominant_mask, sensors);

PlanTable = struct2table([all_plan; dominant_plan]);

consensus_eta = all_plan.eta_mm;
consensus_iqr = all_plan.eta_iqr_mm;
if numel(consensus_eta) ~= numel(sensors)
    consensus_eta = nan(1, numel(sensors));
end
if numel(consensus_iqr) ~= numel(sensors)
    consensus_iqr = nan(1, numel(sensors));
end
consensus_eta(1) = 0;

Reference = table(sensors(:), nan(numel(sensors), 1), ...
    consensus_eta(:), consensus_iqr(:), ...
    'VariableNames', {'sensor_id', 'step04_xc_eta_mm', ...
    'free_window_eta_median_mm', 'free_window_eta_iqr_mm'});

Consensus = struct();
Consensus.status = "ok";
Consensus.diagnostic_status = "ok";
Consensus.reference_sensor_id = sensors(1);
Consensus.eta_source_policy = string(opt.ConsensusPolicy);
Consensus.primary_plan_name = "all_ok";
Consensus.eta_median_mm = consensus_eta;
Consensus.eta_plan_iqr_mm = consensus_iqr;
Consensus.eta_window_iqr_mm = consensus_iqr;
Consensus.max_abs_plan_iqr_mm = max(abs(consensus_iqr), [], 'omitnan');

Summary = struct();
Summary.DatasetDir = char(opt.DatasetDir);
Summary.NoEtaResultFile = "";
Summary.FreeEtaResultFile = string(result_file);
Summary.Step04TemplateFile = string(opt.Step04TemplateFile);
Summary.OutputDir = string(output_dir);
Summary.TargetBlade = target_blade;
Summary.SensorIDs = sensors;
Summary.ReferenceSensorID = sensors(1);
Summary.ReportedEO = dominant_eo;
Summary.Model = ['fixed eta from per-window FreeEta pilot: ' ...
    'V = T(x - dx_c - eta_s - A*sin(EO*theta+phi))'];
Summary.SignConvention = ...
    'Step05 forward model uses x_query = x - dx_c - eta_s - A*sin(EO*theta+phi).';
Summary.Options = struct( ...
    'source_script', string(opt.SourceScript), ...
    'consensus_policy', string(opt.ConsensusPolicy));
Summary.PlanTable = PlanTable;
Summary.Reference = Reference;
Summary.Consensus = Consensus;

mat_file = fullfile(output_dir, sprintf( ...
    'JointStaticEtaPreview_B%d_%s.mat', target_blade, sensor_tag));
plans_csv = fullfile(output_dir, sprintf( ...
    'JointStaticEtaPreview_Plans_B%d_%s.csv', target_blade, sensor_tag));
reference_csv = fullfile(output_dir, sprintf( ...
    'JointStaticEtaPreview_Reference_B%d_%s.csv', target_blade, sensor_tag));

save(mat_file, 'Summary', '-v7.3');
writetable(PlanTable, plans_csv);
writetable(Reference, reference_csv);

Summary.PreviewMatFile = string(mat_file);
Summary.PlansCsvFile = string(plans_csv);
Summary.ReferenceCsvFile = string(reference_csv);

fprintf('Saved fixed eta preview: %s\n', mat_file);
fprintf('Fixed eta B%d %s: %s mm; IQR: %s mm\n', ...
    target_blade, sensor_tag, mat2str(consensus_eta, 8), ...
    mat2str(consensus_iqr, 8));
end


function eta = collect_eta_matrix_local(Result, T, sensors)
eta = nan(height(T), numel(sensors));
for iw = 1:min(numel(Result.WindowResult), height(T))
    if isfield(Result.WindowResult(iw), 'Result') && ...
            isfield(Result.WindowResult(iw).Result, 'sensor_eta_id')
        e = Result.WindowResult(iw).Result.sensor_eta_id(:).';
        n = min(numel(e), numel(sensors));
        eta(iw, 1:n) = e(1:n);
    end
end
end


function plan = summarize_plan_local(plan_name, T, eta, mask, sensors)
mask = mask(:) & all(isfinite(eta), 2);
eta_plan = nan(1, numel(sensors));
eta_iqr = nan(1, numel(sensors));
if any(mask)
    eta_plan = median(eta(mask, :), 1, 'omitnan');
    eta_iqr = iqr(eta(mask, :), 1);
    eta_plan(1) = 0;
end

plan = struct();
plan.plan_name = string(plan_name);
plan.status = "ok";
plan.failure_reason = "";
if ~any(mask)
    plan.status = "failed";
    plan.failure_reason = "no usable FreeEta windows";
end
plan.window_count = nnz(mask);
plan.used_windows = string(mat2str(T.window_id(mask).'));
plan.eo_sequence = string(strjoin(compose('%d', T.EO_id(mask).'), '|'));
plan.weighted_voltage_rmse_median_v = ...
    median(T.weighted_voltage_rmse(mask), 'omitnan');
plan.A_median_mm = median(T.A_id(mask), 'omitnan');
plan.dx_c_median_mm = median(T.dx_c_id(mask), 'omitnan');
plan.eta_mm = eta_plan;
plan.eta_iqr_mm = eta_iqr;
for is = 1:numel(sensors)
    plan.(sprintf('eta_CH%d_mm', sensors(is))) = eta_plan(is);
    plan.(sprintf('eta_CH%d_iqr_mm', sensors(is))) = eta_iqr(is);
end
end
