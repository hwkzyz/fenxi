function C = R5_CompareSharedIndependentLatent_20241106(localizationFile, outputFile)
%R5_COMPARESHAREDINDEPENDENTLATENT_20241106 Diagnostic for S5/S7 only.
cfg = Config_20241106();
if nargin < 1 || isempty(localizationFile)
    localizationFile = fullfile(cfg.paths.results, 'r5_surface_localization', 'loc.mat');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results, 'r5_shared_vs_independent_latent.csv');
end
S = load(localizationFile, 'T', 'details');
T = S.T; details = S.details;
assert(all(ismember(cfg.case.gapSensors, unique(T.sensor_id))), ...
    'R5:GapSensorMissing', 'Localization file does not contain all configured gap sensors.');
keep = ismember(T.sensor_id, cfg.case.gapSensors);
T = T(keep, :); details = details(keep);
rows = cell(numel(details), 1);
for k = 1:numel(details)
    p = details(k).profile_rmse;
    [independentRmse, linearIndex] = min(p(:));
    [iz, ig] = ind2sub(size(p), linearIndex);
    zIndependent = details(k).profile_z(iz);
    gIndependent = details(k).profile_gap(ig);
    zShared = T.latent_coordinate(k);
    [~, izShared] = min(abs(details(k).profile_z - zShared));
    [sharedRmse, igShared] = min(p(izShared, :));
    rows{k} = table(T.blade_id(k), T.sensor_id(k), zIndependent, ...
        gIndependent, independentRmse, zShared, details(k).profile_gap(igShared), ...
        sharedRmse, sharedRmse - independentRmse, ...
        (sharedRmse - independentRmse) / max(independentRmse, eps) * 100, ...
        'VariableNames', {'blade_id','sensor_id','independent_latent', ...
        'independent_gap_mm','independent_rmse_mv','shared_latent', ...
        'shared_gap_mm','shared_rmse_mv','shared_penalty_mv', ...
        'shared_penalty_percent'});
end
C = vertcat(rows{:});
if ~exist(fileparts(outputFile), 'dir'), mkdir(fileparts(outputFile)); end
writetable(C, outputFile);
save(strrep(outputFile, '.csv', '.mat'), 'C', 'localizationFile', '-v7.3');
end
