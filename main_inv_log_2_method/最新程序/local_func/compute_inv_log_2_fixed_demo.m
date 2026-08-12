function D = compute_inv_log_2_fixed_demo(ctx, gTrue)
%compute_inv_log_2_fixed_demo  Shared calculation for fixed-trust figures.
% Uses the current default main entry rather than the historical
% single-start implementation.

if nargin < 2 || isempty(gTrue)
    gTrue = 0.2;
end

trustInfo = load_fixed_trust_domain(ctx.thisDir);
demo = make_inv_log_2_demo_case(ctx, gTrue);
modelInv = get_inv_log_2_model_def();
templateInv = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ...
    ctx.yCell, gTrue, ctx.cfgAna.xGridN, modelInv, trustInfo);
highMap = map_highspeed_to_space(demo.dataHigh, demo.cfgCase.alpha_k, ...
    demo.cfgCase.R_tip, templateInv.domain, demo.cfgCase.fitActiveLevel);
staticState = estimate_highspeed_static_gap_raw(highMap, templateInv, demo.cfgCase);
gapState = estimate_gap_init_vib_basis_projected(highMap, templateInv, demo.cfgCase, staticState);
result = run_inv_log_2_main_method(highMap, templateInv, demo.cfgCase, staticState);

D = demo;
D.trustInfo = trustInfo;
D.templateInv = templateInv;
D.highMap = highMap;
D.staticState = staticState;
D.gapState = gapState;
D.result = result;
end
