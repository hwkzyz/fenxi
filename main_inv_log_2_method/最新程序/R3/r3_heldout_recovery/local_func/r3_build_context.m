function R3 = r3_build_context(mainDir,P)
%R3_BUILD_CONTEXT Load the FE library and build a calibration-only family.

addpath(mainDir,'-begin');
addpath(fullfile(mainDir,'local_func'),'-begin');
% Optional historical diagnostics are archived; add only when present.
legacyLocal = fullfile(mainDir,'identifiability_mechanism_analysis','local_func');
if exist(legacyLocal,'dir') == 7
    addpath(legacyLocal,'-begin');
end
ctx = load_inv_log_2_project_context(mainDir);
trustInfo = load_fixed_trust_domain(mainDir);

assert_same_grid(ctx.gapList(:).',P.gAllMm,'full FE gap library');
isCal = ismembertol(ctx.gapList,P.gCalMm,1e-12,'DataScale',1);
if nnz(isCal)~=numel(P.gCalMm)
    error('r3:CalibrationSplit','Calibration gap selection is incomplete.');
end
if any(ismembertol(P.gOperatingMm,P.gCalMm,1e-12,'DataScale',1))
    error('r3:SplitLeakage','Operating states overlap calibration states.');
end

modelDef = get_inv_log_2_model_def();
libCal = make_response_template_library(ctx.gapList(isCal),ctx.xCell(isCal),...
    ctx.yCell(isCal),NaN,ctx.cfgAna.xGridN,modelDef,struct('enable',false));
libCal = restrict_gap_template_domain(libCal,trustInfo.domain);
libCal.r3_calibration_gap_states_mm = P.gCalMm(:);
libCal.r3_held_out_operating_states_mm = P.gOperatingMm(:);

cfg = ctx.cfgAna;
cfg.RPM_low = min(cfg.RPM_high,300);
cfg.NumRevs_low = 20;
cfg.NumRevs_high = 8;
cfg.route30ForwardModel = "low_increment";
cfg.route30LowTemplateMethod = "adaptive_sg";
cfg.route30SupportAware = false;
cfg.route30FrequencyStructureMode = "single_sync";
cfg.route30SingleFrequencyRangeHz = P.frequencyRangeHz;
cfg.nearSyncDeltaFreqBoundsHz = [-2 2];
cfg.nearSyncDeltaFreqStepHz = .25;
cfg.nearSyncTopK = 100;
cfg.route30GapHalfWidthMm = max(abs(P.gapBoundsMm-P.gReferenceMm));
cfg.route30GapCount = 61;
cfg.route30SingleGapHalfWidthMm = 0.12;
cfg.route30AmplitudeUpperMm = 1.5;
cfg.route30UseAllTrustedSamples = true;
cfg.syncSingleKeepPerGap = 25;
cfg.syncSingleRefineCount = 12;
cfg.syncSingleReplayMaxSamples = 1024;
cfg.syncSingleIteratedReplayCount = 100;
cfg.route30LowTemplateOptions = struct('gridSpacingMm',.02,...
    'minBinCount',5,'minCoverage',.90,'numFolds',5,...
    'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);

R3 = struct('P',P,'ctx',ctx,'trustInfo',trustInfo,'libCal',libCal,...
    'cfg',cfg,'mainDir',mainDir,'isCalibrationState',isCal);
end

function assert_same_grid(actual,expected,label)
if numel(actual)~=numel(expected)||any(abs(actual(:)-expected(:))>1e-12)
    error('r3:GapGridMismatch','Unexpected %s.',label);
end
end
