%% Step06: gap direction versus tilt direction decoupling.
%
% In a first-order model the two static columns are
%
%   J_g  ~= dV/dg
%   J_mu ~= dV/dg * (x - tau)
%
% If tau is not the weighted center of the observed waveform support, J_mu
% contains a gap-offset component. This script computes the Gram-Schmidt
% correction needed to remove that component, then checks the remaining
% coupling with the vibration subspace.

cfg = study_config_20250527();
priorPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_prior.mat');
if ~isfile(priorPath)
    error('Run Step01 first. Missing file: %s', priorPath);
end
if ~isfile(cfg.resultFileS136)
    error('Missing Step5 S136 result: %s', cfg.resultFileS136);
end

P = load(priorPath, 'prior');
prior = P.prior;
S = load(cfg.resultFileS136, 'Result_Struct');
baseBundle = S.Result_Struct.BestWindow.bundle;

rows = table();
for iSet = 1:numel(cfg.referenceSensorSets)
    sensorSet = cfg.referenceSensorSets{iSet};
    sensorTag = cfg.referenceSensorTags(iSet);
    bundle = synthetic_gap_bias_utils_20250527('filter_bundle_by_sensors', ...
        baseBundle, sensorSet);
    bundle = synthetic_gap_bias_utils_20250527('thin_bundle', ...
        bundle, cfg.maxFitPoints);

    D = synthetic_gap_bias_utils_20250527('build_jacobian_columns', ...
        bundle, cfg, prior, []);
    vibBasis = [D.J_A, D.J_phi, D.J_d0];

    for iSensor = 1:numel(D.sensorIds)
        sid = D.sensorIds(iSensor);
        Jg = D.J_dg(:, iSensor);
        Jmu = D.J_dmu(:, iSensor);
        w = D.weight(:);

        gapAlpha = weighted_inner_local(Jg, Jmu, w) ./ ...
            max(weighted_inner_local(Jg, Jg, w), eps);
        JmuGapPerp = Jmu - gapAlpha .* Jg;

        protectedBasis = [vibBasis, Jg];
        JmuProtectedPerp = synthetic_gap_bias_utils_20250527( ...
            'residualize_columns', Jmu, protectedBasis, w);

        row = table();
        row.sensors = sensorTag;
        row.sensor_id = sid;
        row.point_count = numel(bundle.T);
        row.rho_gap_tilt_raw = synthetic_gap_bias_utils_20250527( ...
            'weighted_abs_corr', Jg, Jmu, w);
        row.gap_projection_of_tilt = synthetic_gap_bias_utils_20250527( ...
            'projection_ratio', Jmu, Jg, w);
        row.gap_equivalent_shift_per_unit_dmu_mm = gapAlpha;
        row.rho_gap_tilt_after_gap_orth = synthetic_gap_bias_utils_20250527( ...
            'weighted_abs_corr', Jg, JmuGapPerp, w);
        row.tilt_norm_fraction_after_gap_orth = synthetic_gap_bias_utils_20250527( ...
            'weighted_norm', JmuGapPerp, w) ./ ...
            max(synthetic_gap_bias_utils_20250527('weighted_norm', Jmu, w), eps);
        row.proj_tilt_on_vib_raw = synthetic_gap_bias_utils_20250527( ...
            'projection_ratio', Jmu, vibBasis, w);
        row.proj_tilt_after_gap_orth_on_vib = synthetic_gap_bias_utils_20250527( ...
            'projection_ratio', JmuGapPerp, vibBasis, w);
        row.proj_tilt_on_protected_vib_gap = synthetic_gap_bias_utils_20250527( ...
            'projection_ratio', Jmu, protectedBasis, w);
        row.tilt_norm_fraction_after_vib_gap_orth = synthetic_gap_bias_utils_20250527( ...
            'weighted_norm', JmuProtectedPerp, w) ./ ...
            max(synthetic_gap_bias_utils_20250527('weighted_norm', Jmu, w), eps);
        rows = [rows; row]; %#ok<AGROW>
    end
end

csvPath = fullfile(cfg.outputDir, 'Step06_gap_tilt_direction_decoupling.csv');
writetable(rows, csvPath);

make_summary_figure_local(rows, cfg);

disp('Step06 gap/tilt direction decoupling:');
disp(rows);
fprintf('Step06 completed.\n  %s\n', csvPath);

function val = weighted_inner_local(a, b, weight)
a = a(:);
b = b(:);
weight = max(weight(:), 0);
valid = isfinite(a) & isfinite(b) & isfinite(weight);
val = sum(weight(valid) .* a(valid) .* b(valid));
end

function make_summary_figure_local(rows, cfg)
idx = rows.sensors == "S136";
if ~any(idx)
    return;
end
R = rows(idx, :);
fig = figure('Name', 'Step06 gap tilt direction decoupling', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 14, 7], ...
    'Visible', 'off');

subplot(1, 2, 1);
bar(categorical(string(R.sensor_id)), ...
    [R.rho_gap_tilt_raw, R.rho_gap_tilt_after_gap_orth], ...
    'LineWidth', 0.4);
ylabel('|weighted correlation|');
xlabel('Sensor');
legend({'raw J_\mu', 'gap-orth J_\mu'}, 'Location', 'northoutside');
title('Gap-tilt direction');
ylim([0, 1]);
box off;

subplot(1, 2, 2);
bar(categorical(string(R.sensor_id)), ...
    [R.proj_tilt_on_vib_raw, R.proj_tilt_after_gap_orth_on_vib, ...
    R.proj_tilt_on_protected_vib_gap], ...
    'LineWidth', 0.4);
ylabel('projection ratio');
xlabel('Sensor');
legend({'raw on vib', 'gap-orth on vib', 'raw on vib+gap'}, ...
    'Location', 'northoutside');
title('Remaining vibration coupling');
ylim([0, 1]);
box off;

pngPath = fullfile(cfg.outputDir, 'Step06_gap_tilt_direction_decoupling.png');
exportgraphics(fig, pngPath, 'Resolution', 220);
close(fig);
end
