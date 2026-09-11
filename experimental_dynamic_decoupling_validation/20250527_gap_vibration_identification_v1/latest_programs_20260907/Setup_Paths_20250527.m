function cfg = Setup_Paths_20250527()
%SETUP_PATHS_20250527 Add only paths owned by the V1 package.
cfg = Config_20250527();
assert(isfolder(cfg.paths.common), 'R5:MissingCommonPath', ...
    'Common R5 worker folder is missing: %s', cfg.paths.common);
addpath(cfg.paths.root, cfg.paths.preparation, cfg.paths.foundation, ...
    cfg.paths.calibration, cfg.paths.gapAware, cfg.paths.utilities);
addpath(cfg.paths.common);
end
