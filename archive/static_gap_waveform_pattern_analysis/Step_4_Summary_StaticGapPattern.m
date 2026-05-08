%% Step 4: summarize which pattern is most plausible
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'static_gap_waveform_pattern_results');
load(fullfile(outDir, 'stage1_pointwise_laws.mat'), 'modelSummary');
load(fullfile(outDir, 'stage2_normalization_collapse.mat'), 'stats');
load(fullfile(outDir, 'stage3_feature_laws.mat'), 'featureLawSummary');
load(fullfile(outDir, 'stage4_standard_shape_model.mat'), 'modelFitSummary');
load(fullfile(outDir, 'stage5_leaveoneout_reconstruction.mat'), 'looSummary', 'sparseSummary');

summary = struct();
summary.pointwise_best_mean_model = best_model_name([ ...
    modelSummary.mean_R2_linear_g, ...
    modelSummary.mean_R2_linear_inv_g, ...
    modelSummary.mean_R2_quad_g, ...
    modelSummary.mean_R2_logexp_g], ...
    {'linear g', 'linear 1/g', 'quadratic g', 'log-exp g'});
summary.rank1_energy_raw = stats.raw_rank1_energy;
summary.rank1_energy_norm = stats.norm_rank1_energy;
summary.mean_pair_rmse_raw = stats.raw_mean_pair_rmse;
summary.mean_pair_rmse_norm = stats.norm_mean_pair_rmse;
summary.peak_best_model = best_model_name([featureLawSummary.R2_peak_linear_g, featureLawSummary.R2_peak_linear_inv_g], ...
    {'linear g', 'linear 1/g'});
summary.area_best_model = best_model_name([featureLawSummary.R2_area_linear_g, featureLawSummary.R2_area_linear_inv_g], ...
    {'linear g', 'linear 1/g'});
summary.fwhm_best_model = best_model_name([featureLawSummary.R2_fwhm_linear_g, featureLawSummary.R2_fwhm_linear_inv_g], ...
    {'linear g', 'linear 1/g'});
summary.width_signature_best_model = best_model_name([featureLawSummary.mean_R2_width_linear_g, featureLawSummary.mean_R2_width_linear_inv_g], ...
    {'linear g', 'linear 1/g'});
summary.standard_shape_mean_nrmse = modelFitSummary.mean_nrmse;
summary.standard_shape_max_nrmse = modelFitSummary.max_nrmse;
summary.standard_shape_best_param_laws = sprintf('base:%s, amp:%s, xc:%s, fwhm:%s', ...
    best_model_name([modelFitSummary.R2_base_g, modelFitSummary.R2_base_inv_g], {'g', '1/g'}), ...
    best_model_name([modelFitSummary.R2_amp_g, modelFitSummary.R2_amp_inv_g], {'g', '1/g'}), ...
    best_model_name([modelFitSummary.R2_xc_g, modelFitSummary.R2_xc_inv_g], {'g', '1/g'}), ...
    best_model_name([modelFitSummary.R2_fwhm_g, modelFitSummary.R2_fwhm_inv_g], {'g', '1/g'}));
summary.loo_mean_nrmse = looSummary.mean_nrmse_all;
summary.loo_edge_nrmse = looSummary.mean_nrmse_edge;
summary.loo_interior_nrmse = looSummary.mean_nrmse_interior;
summary.sparse3_mean_nrmse = sparseSummary.mean_nrmse_all;
summary.sparse3_edge_nrmse = sparseSummary.mean_nrmse_edge;
summary.sparse3_interior_nrmse = sparseSummary.mean_nrmse_interior;

save(fullfile(outDir, 'stage4_pattern_summary.mat'), 'summary', '-v7.3');

figure('Name', 'Static Gap Step 4 - Summary', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [3, 3, 20, 12]);
axis off;
text(0.02, 0.90, 'Static Gap Waveform Pattern Summary', 'FontSize', 15, 'FontWeight', 'bold');
text(0.02, 0.76, sprintf('Best pointwise mean model: %s', summary.pointwise_best_mean_model), 'FontSize', 11);
text(0.02, 0.66, sprintf('Rank-1 energy raw -> norm: %.4f -> %.4f', ...
    summary.rank1_energy_raw, summary.rank1_energy_norm), 'FontSize', 11);
text(0.02, 0.56, sprintf('Mean pair RMSE raw -> norm: %.4g -> %.4g', ...
    summary.mean_pair_rmse_raw, summary.mean_pair_rmse_norm), 'FontSize', 11);
text(0.02, 0.46, sprintf('Peak feature best law: %s', summary.peak_best_model), 'FontSize', 11);
text(0.02, 0.36, sprintf('Area feature best law: %s', summary.area_best_model), 'FontSize', 11);
text(0.02, 0.26, sprintf('FWHM feature best law: %s', summary.fwhm_best_model), 'FontSize', 11);
text(0.02, 0.16, sprintf('Width-signature best law: %s', summary.width_signature_best_model), 'FontSize', 11);
text(0.02, 0.08, sprintf('Standard-shape NRMSE mean/max: %.4g / %.4g', ...
    summary.standard_shape_mean_nrmse, summary.standard_shape_max_nrmse), 'FontSize', 11);
text(0.52, 0.76, sprintf('LOO NRMSE all/interior/edge: %.4g / %.4g / %.4g', ...
    summary.loo_mean_nrmse, summary.loo_interior_nrmse, summary.loo_edge_nrmse), 'FontSize', 11);
text(0.52, 0.64, sprintf('Sparse-3 NRMSE all/interior/edge: %.4g / %.4g / %.4g', ...
    summary.sparse3_mean_nrmse, summary.sparse3_interior_nrmse, summary.sparse3_edge_nrmse), 'FontSize', 11);
text(0.52, 0.52, sprintf('Best parameter laws: %s', summary.standard_shape_best_param_laws), 'FontSize', 11);

fprintf('[Static Step 4] Summary:\n');
disp(summary);

function name = best_model_name(values, names)
    [~, idx] = max(values);
    name = names{idx};
end
