function ctx = load_inv_log_2_project_context(thisDir)
%load_inv_log_2_project_context  Common project setup for inv_log_2 scripts.

if nargin < 1 || isempty(thisDir)
    thisDir = fileparts(mfilename('fullpath'));
end

projectDir = fileparts(thisDir);
rootDir = fileparts(projectDir);

addpath(fullfile(thisDir, 'local_func'));
addpath(thisDir, '-begin');

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(projectDir, 'results', '04_inv_log_2_main');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
cfgAna.twoVpTopKFreqCandidates = 5;
cfgAna.twoVpJointCandidateKeep = 2;
cfgAna.twoVpJointMuList = [0, 0.15, 0.75];
cfgAna.twoVpDxRefineHalfWidth = 0.05;
cfgAna.gapTrustEnable = true;
cfgAna.gapTrustMinHalfWidth = 2.5;
cfgAna.gapTrustStepHalfWidth = 0.25;
cfgAna.gapTrustScoreTolerance = 0.08;
cfgAna.gapTrustDerivativeWeight = 0.20;
cfgAna.gapTrustSensitivityWeight = 0.30;
cfgAna.gapTrustNarrowPenaltyWeight = 0.03;
cfgAna.gapTrustMinResponseFraction = 0.12;
cfgAna.gapTrustMinSensitivityFraction = 0.12;
cfgAna.gapTrustMinGridPoints = 40;

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

ctx = struct();
ctx.thisDir = thisDir;
ctx.projectDir = projectDir;
ctx.rootDir = rootDir;
ctx.oldDir = oldDir;
ctx.outDir = outDir;
ctx.cfg = cfg;
ctx.cfgAna = cfgAna;
ctx.gapList = gapList;
ctx.xCell = xCell;
ctx.yCell = yCell;
end
