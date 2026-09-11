function cfg = Setup_Paths_20241106()
%SETUP_PATHS_20241106 Add only the migrated package paths to MATLAB.
cfg = Config_20241106();
addpath(cfg.paths.root);
assert(isfolder(cfg.paths.common), 'R5:MissingCommonPath', ...
    'Common R5 worker folder is missing: %s', cfg.paths.common);
addpath(cfg.paths.common);
addpath(fullfile(cfg.paths.analysis, 'r5_anchor_guided_experimental_20260903'));
if ~exist(cfg.paths.results, 'dir'), mkdir(cfg.paths.results); end
end
