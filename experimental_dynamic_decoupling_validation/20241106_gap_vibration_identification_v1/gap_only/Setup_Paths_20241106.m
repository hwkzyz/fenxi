function cfg = Setup_Paths_20241106()
%SETUP_PATHS_20241106 Add only the migrated package paths to MATLAB.
cfg = Config_20241106();
addpath(cfg.paths.root);
addpath(fullfile(cfg.paths.analysis, 'r5_anchor_guided_experimental_20260903'));
if ~exist(cfg.paths.results, 'dir'), mkdir(cfg.paths.results); end
end
