function demo = make_inv_log_2_demo_case(ctx, gTrue)
%make_inv_log_2_demo_case  Common high-speed demo data for fixed-trust plots.

if nargin < 2 || isempty(gTrue)
    gTrue = 0.2;
end

cfgCase = ctx.cfgAna;
cfgCase.g_holdout = gTrue;
cfgCase.A_true = [0.25, 0.15];
cfgCase.f_true = [500, 1300];
cfgCase.phi_true = [pi/4, -pi/3];
cfgCase.noiseMode = 'noise_ratio';
cfgCase.snrDb = 10;

templateTruth = build_gap_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    NaN, ctx.cfgAna.xGridN);
cfgCase.noiseRatio = calibrate_demo_snr_ratio(cfgCase, ctx.gapList, ...
    ctx.xCell, ctx.yCell, templateTruth.domain, cfgCase.snrDb);

idxTrue = find(abs(ctx.gapList - gTrue) < 1e-12, 1);
if isempty(idxTrue)
    error('Holdout gap %.6g mm not found in data file.', gTrue);
end
FTrue = griddedInterpolant(ctx.xCell{idxTrue}, ctx.yCell{idxTrue}, ...
    'pchip', 'nearest');

rng(20261205 + round(1000 * gTrue), 'twister');
dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    cfgCase.noiseRatio, templateTruth.domain, cfgCase.A_true, ...
    cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);

demo = struct();
demo.cfgCase = cfgCase;
demo.templateTruth = templateTruth;
demo.dataHigh = dataHigh;
demo.idxTrue = idxTrue;
demo.gTrue = gTrue;
end

function noiseRatio = calibrate_demo_snr_ratio(cfgCase, gapList, xCell, yCell, domain, snrDb)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
dataClean = simulate_rotating_waveform_from_template(Fx, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    0, domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, 'noise_ratio');
signal = dataClean.V_clean(:) - min(dataClean.V_clean(:));
signalRms = sqrt(mean(signal .^ 2));
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, ...
    linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatio = signalRms ./ (templateRange .* 10^(snrDb / 20));
end
