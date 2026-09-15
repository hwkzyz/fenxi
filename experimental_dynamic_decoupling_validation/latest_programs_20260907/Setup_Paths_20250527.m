function cfg = Setup_Paths_20250527()
%SETUP_PATHS_20250527 Add only paths owned by the V1 package.
cfg = Config_20250527();
addpath(cfg.paths.root, cfg.paths.preparation, cfg.paths.foundation, ...
    cfg.paths.calibration, cfg.paths.gapAware, cfg.paths.utilities);
end
