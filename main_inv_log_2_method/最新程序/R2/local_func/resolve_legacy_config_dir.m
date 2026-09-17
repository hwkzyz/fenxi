function oldDir = resolve_legacy_config_dir(rootDir)
%resolve_legacy_config_dir  Locate the saved legacy stage0 configuration.

candidates = {
    fullfile(rootDir, 'data', 'gapaware_varpro_holdout_results_modular')
    fullfile(rootDir, 'gapaware_varpro_holdout_results_modular')
    fullfile(rootDir, 'archive', 'gapaware_varpro_holdout_results_modular')
    fullfile(rootDir, 'archive', 'blind_template_highspeed_results')
    };

for ii = 1:numel(candidates)
    if exist(fullfile(candidates{ii}, 'stage0_config.mat'), 'file')
        oldDir = candidates{ii};
        return;
    end
end

error('Could not locate stage0_config.mat under the expected legacy data folders.');
end
