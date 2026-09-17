function Result = Run_LowSpeedGapAudit(runMode)
%RUN_LOWSPEEDGAPAUDIT Validate low-speed waveform generation and gap recovery.
% This is an isolated calibration audit before connecting the result to the
% high-speed fast profile bridge.
if nargin < 1, runMode = "smoke"; end
runMode = lower(string(runMode));
root = fileparts(mfilename('fullpath'));
mainDir = fileparts(root);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');
ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
modelDef = get_inv_log_2_model_def();
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    NaN, ctx.cfgAna.xGridN, modelDef, trust);

if runMode == "smoke"
    gaps = [0.2 0.8];
    snrs = [Inf 25];
else
    gaps = [0.2 0.5 0.8 1.1 1.4];
    snrs = [Inf 25 20];
end

rows = repmat(row0(), 0, 1);
for ig = 1:numel(gaps)
    gTrue = gaps(ig);
    idx = find(abs(ctx.gapList - gTrue) < 1e-12, 1);
    if isempty(idx), error('Gap %.3f is not in the static library.', gTrue); end
    FxTrue = griddedInterpolant(ctx.xCell{idx}, ctx.yCell{idx}, ...
        'pchip', 'nearest');
    for isnr = 1:numel(snrs)
        snrDb = snrs(isnr);
        seed = 20260806 + 100 * ig + isnr;
        cfg = ctx.cfgAna;
        cfg.RPM_low = min(cfg.RPM_high, 300);
        cfg.NumRevs_low = 8;
        dataLow = simulate_low_speed_template(FxTrue, cfg, lib.domain, snrDb, seed);
        low = aggregate_low_speed_template(dataLow, cfg.alpha_k, cfg.R_tip, ...
            lib.domain, lib.xGrid);
        state = estimate_highspeed_static_gap_raw(low.mapped, lib, cfg);
        r = row0();
        r.g_true_mm = gTrue;
        r.g_hat_mm = state.gHat;
        r.g_error_mm = state.gHat - gTrue;
        r.dx_hat_mm = state.dx0;
        r.snr_db = snrDb;
        r.seed = seed;
        r.low_rmse_V = sqrt(mean((low.templateLow - eval_gap_template(lib, gTrue, lib.xGrid)).^2));
        r.success = abs(r.g_error_mm) <= 0.04;
        rows(end+1, 1) = r; %#ok<AGROW>
    end
end

Audit = struct2table(rows);
out = fullfile(root, 'output', char(runMode));
if ~exist(out, 'dir'), mkdir(out); end
writetable(Audit, fullfile(out, 'low_speed_gap_audit.csv'));
save(fullfile(out, 'low_speed_gap_audit.mat'), 'Audit', '-v7.3');
Result = struct('Audit', Audit, 'outputDir', out);
end

function r = row0()
r = struct('g_true_mm', NaN, 'g_hat_mm', NaN, 'g_error_mm', NaN, ...
    'dx_hat_mm', NaN, 'snr_db', NaN, 'seed', NaN, 'low_rmse_V', NaN, ...
    'success', false);
end
