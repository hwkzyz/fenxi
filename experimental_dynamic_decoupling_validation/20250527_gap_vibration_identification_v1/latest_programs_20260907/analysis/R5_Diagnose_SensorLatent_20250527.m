function D = R5_Diagnose_SensorLatent_20250527(localizationFile, outputFile, familyFile)
%R5_DIAGNOSE_SENSORLATENT_20250527 Diagnose sensor-wise latent disagreement.
cfg = Config_20250527();
if nargin < 1 || isempty(localizationFile)
    localizationFile = fullfile(cfg.paths.results, 'r5_surface_localization', 'loc.mat');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results, 'r5_sensor_latent_diagnostic.csv');
end
if nargin < 3 || isempty(familyFile)
    familyFile = fullfile(cfg.paths.results, 'r5_experimental_surface_family.mat');
end
S = load(localizationFile, 'T');
T = S.T;
F = load(familyFile, 'F');
familyRange = range(F.F.latent_sorted);
assert(all(ismember({'blade_id','sensor_id','latent_coordinate','gap_or_effective_mm', ...
    'anchor_rmse_mv'}, T.Properties.VariableNames)), ...
    'R5:LocalizationSchema', 'Localization table lacks required diagnostic columns.');
rows = {};
for bid = unique(T.blade_id(:)).'
    useB = T.blade_id == bid;
    s = T(useB, :);
    if ismember('sensorwise_latent_coordinate', T.Properties.VariableNames)
        z = double(s.sensorwise_latent_coordinate);
    else
        z = double(s.latent_coordinate);
    end
    z = z(isfinite(z));
    zRange = range(z);
    rows(end+1, :) = {bid, height(s), min(z), max(z), zRange, ...
        zRange / max(familyRange, eps), mean(s.anchor_rmse_mv, 'omitnan'), ...
        min(s.support_fraction), max(s.support_fraction)}; %#ok<AGROW>
end
D = cell2table(rows, 'VariableNames', {'blade_id','sensor_count', ...
    'latent_min','latent_max','latent_range','latent_range_fraction', ...
    'mean_anchor_rmse_mv','support_fraction_min','support_fraction_max'});
if ~exist(fileparts(outputFile), 'dir'), mkdir(fileparts(outputFile)); end
writetable(D, outputFile);
save(strrep(outputFile, '.csv', '.mat'), 'D', 'localizationFile', 'familyFile', '-v7.3');
end
