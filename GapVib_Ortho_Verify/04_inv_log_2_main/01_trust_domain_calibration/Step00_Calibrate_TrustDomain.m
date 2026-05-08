%% Step 00: calibrate and save the inv_log_2 trusted fitting domain.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gHoldout', 'var') || isempty(gHoldout)
    gHoldout = 0.2;
end

modelInv = get_inv_log_2_model_def();
trustOptions = make_gap_trust_options(ctx.cfgAna);
templateInv = make_response_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    gHoldout, ctx.cfgAna.xGridN, modelInv, trustOptions);

save_fixed_trust_domain(ctx.thisDir, templateInv.trustInfo, modelInv, ctx.cfgAna);

disp('Saved fixed trust-domain result:');
disp(get_trust_result_path(ctx.thisDir));
disp(struct2table(rmfield(templateInv.trustInfo, 'metrics')));
