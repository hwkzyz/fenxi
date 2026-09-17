function Result = Run_R2_UnifiedBlindBackend(nCases)
%RUN_R2_UNIFIEDBLINDBACKEND  R2 数值后端。
% 主流程保持“配置→建模→逐案例仿真→盲辨识→汇总”的顺序，
% 具体数学运算由 local_func 中的稳定算法函数完成。
if nargin < 1, nCases = 6; end

%% 1. 项目路径与公共数据
root = fileparts(mfilename('fullpath'));
mainDir = fileparts(fileparts(root));
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');
ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
modelDef = get_inv_log_2_model_def();
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    NaN, ctx.cfgAna.xGridN, modelDef, trust);

%% 2. 仿真配置（真值只用于生成观测，不传入辨识器）
base = make_inv_log_2_demo_case(ctx, 0.2);
cfg = base.cfgCase;
cfg.RPM_low = min(cfg.RPM_high, 300);
cfg.NumRevs_low = 20;
cfg.NumRevs_high = 8;
cfg.route30ForwardModel = "low_increment";
cfg.route30FrequencyStructureMode = "single_sync";
eoList = [7 10 13 16 19 26];
ampList = [0.20 0.25 0.16 0.22 0.18 0.24];
gLow = 0.5;
gHigh = 0.6;
truthT = repmat(eval_gap_template(lib, gLow, lib.xGrid), 1, numel(cfg.alpha_k));
truthCal = struct('g0', gLow, 'mu', 0, 'tau', 0, 'zeta', 1, ...
    'kappa', 1, 'b', 0, 'sensorIds', 1:numel(cfg.alpha_k));
truthModel = build_path_template_model(lib, truthT, truthCal);
rows = repmat(row0(), 0, 1);

%% 3. 逐案例生成观测并执行盲辨识
for ic=1:min(nCases,numel(eoList))
    eo = eoList(ic);
    amp = ampList(ic);
    freq = eo * cfg.RPM_high / 60;
    low = simulate_low_speed_template(@(z) eval_gap_template(lib, gLow, z), ...
        cfg, lib.domain, Inf, 88000 + ic, 'fixed_std');
    calibration = calibrate_inv_log_2_low_speed(low, lib, cfg);
    cfg.A_true = amp;
    cfg.f_true = freq;
    cfg.phi_true = pi/4;
    high = simulate_highspeed_from_low_increment(truthModel, gHigh-gLow, ...
        cfg, Inf, 'fixed_std', 89000 + ic);
    mapped = map_highspeed_to_space(high, cfg.alpha_k, cfg.R_tip, ...
        lib.domain, 0.02, 0.30);
    bundle = struct('t', mapped.t_v, 'x', mapped.x_v, 'V', mapped.V_a, ...
        'theta', mapped.theta_v, 'sensorId', mapped.S_v, 'turnId', mapped.rev_v, ...
        'Wvp', ones(size(mapped.V_a)), 'Wfull', ones(size(mapped.V_a)));
    response = struct('evalF', @(g,x) eval_gap_template(calibration.templateModel,g,x), ...
        'evalFx', @(g,x) eval_gap_derivative(calibration.templateModel,g,x));
    % Formal blind protocol: the estimator must scan the complete admissible
    % EO domain.  eoList is used only to select truth cases for waveform
    % generation; passing it here would leak the true EO candidates.
    bc=backend_cfg([min(ctx.gapList),max(ctx.gapList)],cfg.RPM_high/60,1:30);
    B = run_unified_blind_dynamic_backend(bundle, response, bc);
    best = B.best;
    rows(end+1) = struct('case_id', ic, 'eo_true', eo, 'eo_hat', best.eo, ...
        'gap_true_mm', gHigh, 'gap_hat_mm', best.gap, ...
        'gap_error_mm', best.gap-gHigh, 'A_true_mm', amp, ...
        'A_hat_mm', best.A, 'delta_f_hat_hz', best.deltaF, ...
        'rmse_V', best.fullWaveRmseV, 'status', B.status); %#ok<AGROW>
end
%% 4. 保存并返回结果
T = struct2table(rows);
out = fullfile(root, 'output', 'unified_blind_backend');
if ~exist(out, 'dir'), mkdir(out); end
writetable(T, fullfile(out, 'r2_unified_blind_backend.csv'));
save(fullfile(out, 'r2_unified_blind_backend.mat'), 'T', '-v7.3');
Result = struct('table', T, 'outputDir', out);
disp(T);
end
function c = backend_cfg(gb, rh, eo)
c = struct('gapBounds', gb, 'gapGrid', gb(1):0.02:gb(2), 'eoGrid', eo, ...
    'topK', 5, 'rotHz', rh, 'minGradient', 1e-6, 'minSamples', 30, ...
    'dxBound', 0.2, 'ampCoeffBound', 1.5, 'deltaFBound', 2, ...
    'maxIter', 250, 'ambiguityMargin', 0.01);
end
function r = row0()
r = struct('case_id', NaN, 'eo_true', NaN, 'eo_hat', NaN, ...
    'gap_true_mm', NaN, 'gap_hat_mm', NaN, 'gap_error_mm', NaN, ...
    'A_true_mm', NaN, 'A_hat_mm', NaN, 'delta_f_hat_hz', NaN, ...
    'rmse_V', NaN, 'status', "");
end
