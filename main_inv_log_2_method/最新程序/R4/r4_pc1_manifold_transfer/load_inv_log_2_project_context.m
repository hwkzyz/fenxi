function ctx = load_inv_log_2_project_context(thisDir)
%load_inv_log_2_project_context  Common project setup for inv_log_2 scripts.

if nargin < 1 || isempty(thisDir)
    thisDir = fileparts(mfilename('fullpath'));
end

thisDir = char(thisDir);
legacyName = 'GapVib_Ortho_Verify';
mainFolderName = 'main_inv_log_2_method';

mainDir = thisDir;
for ii = 1:8
    hasLocalFunc = exist(fullfile(mainDir, 'local_func'), 'dir') == 7;
    hasContextFile = exist(fullfile(mainDir, 'load_inv_log_2_project_context.m'), 'file') == 2;
    if hasLocalFunc && hasContextFile
        break;
    end
    parentDir = fileparts(mainDir);
    if isempty(parentDir) || strcmp(parentDir, mainDir)
        break;
    end
    mainDir = parentDir;
end

if exist(fullfile(mainDir, 'local_func'), 'dir') ~= 7
    error('Could not locate %s root from: %s', mainFolderName, thisDir);
end

parentDir = fileparts(mainDir);
[~, parentName] = fileparts(parentDir);
if strcmp(parentName, mainFolderName)
    % The maintained program is grouped under main_inv_log_2_method/最新程序.
    % Keep project data and result paths anchored at the original project root.
    rootDir = fileparts(parentDir);
elseif strcmpi(parentName, legacyName)
    rootDir = fileparts(parentDir);
else
    rootDir = parentDir;
end

addpath(fullfile(mainDir, 'local_func'));
addpath(mainDir, '-begin');

oldDir = resolve_legacy_config_dir(rootDir);
outDir = fullfile(rootDir, 'results', mainFolderName);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(rootDir, cfg.dataFile);
cfg.legacyDataFile = cfg.dataFile;
mergedDataFile = fullfile(rootDir, 'data', '间隙的影响', '直叶片2mm_统一间隙库0.2_0.1_1.5mm.txt');
if exist(mergedDataFile, 'file') == 2
    cfg.dataFile = mergedDataFile;
end
cfgAna = make_multicase_cfg(cfg);
cfgAna.twoVpTopKFreqCandidates = 5;
cfgAna.twoVpJointCandidateKeep = 2;
cfgAna.twoVpJointMuList = [0, 0.15, 0.75];
cfgAna.twoVpDxRefineHalfWidth = 0.05;
cfgAna.enhancedUseMidpointSeed = false;
cfgAna.enhancedFreqTopK = 3;
cfgAna.enhancedJointSeedKeep = 2;
cfgAna.finalDesignAmplitudeLimitMm = 1.5;
cfgAna.finalAmplitudeUpperBoundMm = 1.5;
cfgAna.finalNonlinearAmpScaleLower = 0.2;
cfgAna.finalNonlinearAmpScaleUpper = 2.5;
cfgAna.finalNonlinearFreqHalfWidth = 80;
cfgAna.finalVpDxHalfWidth = 0.20;
cfgAna.finalProjectedJointDxHalfWidth = 0.20;
cfgAna.finalProjectedJointGapHalfWidth = 0.12;
cfgAna.querySafeAmplitudeBoundMm = cfgAna.finalAmplitudeUpperBoundMm;
cfgAna.querySafeDxBoundMm = cfgAna.finalProjectedJointDxHalfWidth;
cfgAna.querySafeEtaBoundMm = 0;
cfgAna.querySafeDomainMarginMm = 0.02;
cfgAna.finalDirectABAmplitudeScales = [0.75, 0.85, 1, 1.25, 1.4];
cfgAna.finalDirectABPhaseGrid = (0:23)*pi/12;
cfgAna.finalDirectABPhaseOffsetMode = false;
cfgAna.finalDirectABMaxInitialStartsPerPair = 4;
cfgAna.finalDirectABCoarseAmplitudeScales = [0.75, 1, 1.25];
cfgAna.finalDirectABCoarsePhaseGrid = (0:3)*pi/2;
cfgAna.finalDirectABCoarseMaxInitialStartsPerPair = 1;
cfgAna.finalDirectABCoarseMaxIter = 35;
cfgAna.finalDirectABDenseTopK = 2;
cfgAna.finalDirectABInitialScreenStride = 8;
cfgAna.finalDirectABInitialScreenKeep = 128;
cfgAna.finalDirectABPolish = false;
cfgAna.finalDirectABZeroDxSeed = true;
cfgAna.finalFastModeEnable = true;
cfgAna.finalFastGapDiffThreshold = 0.02;
cfgAna.finalFastVpRhoJThreshold = 0.20;
cfgAna.twoVpTopKFreqCandidatesFast = 1;
cfgAna.twoVpJointCandidateKeepFast = 1;
cfgAna.twoVpJointMuListFast = 0;
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
% Determine the high-speed projected gap interval from the full static
% library profile; no hand-selected gap center is required.
cfgAna.gapProjFullDomain = true;
cfgAna.gapProjFullDomainN = 121;
cfgAna.gapProjProfileRelativeTolerance = 0.10;

[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);

ctx = struct();
ctx.thisDir = mainDir;
ctx.projectDir = rootDir;
ctx.rootDir = rootDir;
ctx.oldDir = oldDir;
ctx.outDir = outDir;
ctx.cfg = cfg;
ctx.cfgAna = cfgAna;
ctx.gapList = gapList;
ctx.xCell = xCell;
ctx.yCell = yCell;
end
