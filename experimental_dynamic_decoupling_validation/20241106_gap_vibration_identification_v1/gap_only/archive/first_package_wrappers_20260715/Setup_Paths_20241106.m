function cfg = Setup_Paths_20241106()
%SETUP_PATHS_20241106 Add only this package's runtime folders to MATLAB path.
cfg = Config_20241106();
addpath(cfg.paths.root);
addpath(cfg.paths.foundation);
addpath(cfg.paths.preparation);
addpath(cfg.paths.gapAware);
addpath(cfg.paths.utilities);
end
