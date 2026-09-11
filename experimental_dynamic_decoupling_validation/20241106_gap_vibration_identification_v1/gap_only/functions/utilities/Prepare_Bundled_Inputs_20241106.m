function Prepare_Bundled_Inputs_20241106(cfg)
%PREPARE_BUNDLED_INPUTS_20241106 Seed active inputs without reading old routes.

copy_tree_local(cfg.paths.bundledPreparedFoundation, cfg.paths.preparedFoundation);
copy_tree_local(cfg.paths.calibrationFoundation, cfg.paths.preparedFoundation);
copy_tree_local(cfg.paths.bundledPreparedGap, cfg.paths.gapRuntime);
ensure_dir_local(cfg.paths.gapRuntime);
ensure_dir_local(cfg.paths.foundationResults);
ensure_dir_local(cfg.paths.gapResults);
ensure_dir_local(cfg.paths.comparison);
ensure_dir_local(cfg.paths.figures);
end

function copy_tree_local(sourceDir, targetDir)
if ~isfolder(sourceDir)
    error('Bundled input folder is missing: %s', sourceDir);
end
ensure_dir_local(targetDir);
[ok, message] = copyfile(fullfile(sourceDir, '*'), targetDir, 'f');
if ~ok
    error('Unable to prepare bundled inputs from %s: %s', sourceDir, message);
end
end

function ensure_dir_local(folder)
if ~isfolder(folder)
    mkdir(folder);
end
end
