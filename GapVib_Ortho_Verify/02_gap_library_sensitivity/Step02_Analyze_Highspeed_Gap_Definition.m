%% Step 02: analyze high-speed truth-gap definition without vibration/noise.

close all;
clc;

caseDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(caseDir);
rootDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '02_gap_library_sensitivity');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

rows = table();
for ih = 1:numel(gapList)
    gTrue = gapList(ih);
    FTrue = griddedInterpolant(xCell{ih}, yCell{ih}, 'pchip', 'nearest');

    cfgCase = cfgAna;
    cfgCase.g_holdout = gTrue;
    cfgCase.A_true = [0, 0];
    cfgCase.f_true = [500, 1300];
    cfgCase.phi_true = [0, 0];
    cfgCase.noiseMode = 'noise_ratio';
    cfgCase.noiseRatio = 0;

    templateIdeal = build_gap_template_library(gapList, xCell, yCell, ...
        NaN, cfgCase.xGridN);
    rows = [rows; run_one_static_gap_case("ideal_in_library", ...
        FTrue, templateIdeal, cfgCase, gTrue)]; %#ok<AGROW>

    templateLeaveout = build_gap_template_library(gapList, xCell, yCell, ...
        gTrue, cfgCase.xGridN);
    rows = [rows; run_one_static_gap_case("leave_one_out_inv_g_linear", ...
        FTrue, templateLeaveout, cfgCase, gTrue)]; %#ok<AGROW>

    sparseMask = ismember(round(gapList, 12), round([0.2; 0.6; 1.0; 1.4], 12));
    sparseMask(ih) = false;
    if nnz(sparseMask) >= 2
        templateSparse = build_gap_template_library(gapList(sparseMask), ...
            xCell(sparseMask), yCell(sparseMask), NaN, cfgCase.xGridN);
        rows = [rows; run_one_static_gap_case("sparse_leave_one_out_inv_g_linear", ...
            FTrue, templateSparse, cfgCase, gTrue)]; %#ok<AGROW>
    end
end

writetable(rows, fullfile(outDir, 'highspeed_gap_definition_error.csv'));
save(fullfile(outDir, 'highspeed_gap_definition_error.mat'), 'rows', '-v7.3');
disp(rows);

function row = run_one_static_gap_case(libraryMode, FTrue, templateLib, cfgCase, gTrue)
rng(20261001 + round(gTrue * 1000), 'twister');
dataHigh = simulate_rotating_waveform_from_template(FTrue, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, ...
    cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
highMap = map_highspeed_to_space(dataHigh, cfgCase.alpha_k, cfgCase.R_tip, ...
    templateLib.domain, cfgCase.fitActiveLevel);

timerStatic = tic;
staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
elapsed = toc(timerStatic);

V = highMap.V_a(:);
x = highMap.x_v(:);
VFit = eval_gap_template(templateLib, staticState.gHat, x - staticState.dx0);
rmse = sqrt(mean((V - VFit).^2, 'omitnan'));
row = table(string(libraryMode), gTrue, numel(templateLib.gapTrain), ...
    staticState.gHat, staticState.gHat - gTrue, abs(staticState.gHat - gTrue), ...
    staticState.dx0, rmse, elapsed, ...
    'VariableNames', {'library_mode','g_true_mm','num_gap_templates', ...
    'g_est_mm','gap_error_signed_mm','gap_error_abs_mm', ...
    'dx0_mm','static_rmse','elapsed_s'});
end
