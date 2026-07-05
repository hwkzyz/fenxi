function P = NewFlow_Config_20251222()
%NEWFLOW_CONFIG_20251222 Central tuning block for the 20251222 NewFlow route.

routeDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(routeDir);
gapPriorDir = fullfile(validationRoot, '20251222_low_speed_gap_prior_decoupling');
legacyDir = fullfile(gapPriorDir, 'legacy');
if exist(gapPriorDir, 'dir') ~= 7
    error('20251222 gap-prior folder not found: %s', gapPriorDir);
end
if exist(legacyDir, 'dir') ~= 7
    error('20251222 legacy helper folder not found: %s', legacyDir);
end
addpath(gapPriorDir);
addpath(legacyDir);

C = CaseConfig();
baseCfg = Get_20251222_BTT_Config();
syncCfg = build_single_sync_experiment_config_20251222();
[resonanceSelection, resonanceCatalog] = ResonanceRegionCatalog_20251222();
dynamicCase = baseCfg.dynamic_cases{1};

P = struct();
P.routeDir = routeDir;
P.validationRoot = validationRoot;
P.gapPriorDir = gapPriorDir;
P.legacyDir = legacyDir;
P.outputDir = fullfile(routeDir, 'output', 'new_flow');
P.dataset = '20251222';

P.case = C;
P.resonance.selection = resonanceSelection;
P.resonance.catalog = resonanceCatalog;

P.data.dynamicCase = dynamicCase;
P.data.dynamicDataDir = fullfile(baseCfg.dataset_root, dynamicCase);
P.data.caseOutputDir = fullfile(baseCfg.output_root, dynamicCase);
P.data.referenceOutputDir = baseCfg.reference_output_dir;
P.data.sensorConfigFile = fullfile(baseCfg.reference_output_dir, 'Sensor_Config_20251222.mat');

P.machine.bladeCount = baseCfg.blades_num;
P.machine.oprChannel = baseCfg.opr_id;
P.machine.oprPulsesPerRev = baseCfg.blades_num;
P.machine.sampleRateHz = baseCfg.pinlv;
P.machine.tipRadiusMM = baseCfg.r_tip_mm;

P.sensors.analysis = C.sensorIds(:).';
P.identification.blades = C.bladeId;
P.identification.freqSearchHz = [100 1000];
P.identification.eoPad = 2;
P.identification.topKEO = 3;
P.identification.amplitudeLimitMM = 0.50;
P.identification.dxCLimitMM = 0.35;
P.identification.sensorEtaLimitMM = 0.20;
P.identification.sensorEtaAdaptiveLimitEnable = true;
P.identification.sensorEtaAdaptiveQuantile = 85;
P.identification.sensorEtaAdaptiveSafetyFactor = 1.20;
P.identification.sensorEtaAdaptiveMinMM = 0.02;
P.identification.sensorEtaAdaptiveMaxMM = 0.20;
P.identification.sensorEtaAdaptiveScaleTable = [];
P.identification.sensorEtaRegWeightVPerMM = 0.05;
P.identification.overshootPenaltyWeight = 100;
P.identification.pulseSelectionMode = 'single';
P.identification.domainSelectionMode = 'hard';
P.identification.domainSoftMarginMM = 0;
P.identification.queryGuardMM = [];
P.identification.queryGuardMode = 'adaptive';
P.identification.queryGuardQuantile = 95;
P.identification.queryGuardSafetyMM = 0.05;
P.identification.queryGuardMinMM = 0.12;
P.identification.queryGuardMaxMM = [];
P.identification.dynamicEffectiveMode = 'gradient';
P.identification.dynamicTemplateGradientMinRatio = 0.08;
P.identification.dynamicTimeGradientMinRatio = 0.15;
P.identification.dynamicPeakQuantile = 85;
P.identification.phaseSafeExpansion = true;
P.identification.phaseSafeMarginMM = 0.03;
P.identification.phaseSafeReferenceMode = 'prev_window';
P.identification.phaseSafeFallbackMode = 'linear_vp';
P.identification.phaseSafeRefreshEvery = 5;
P.identification.phaseSafePrevMaxClampFraction = 1e-6;
P.identification.phaseSafePrevMinFinalGapRatio = 0.005;
P.identification.debugMaxWindows = inf;
P.identification.allEoWarmupWindows = 0;
P.identification.vpGapRatioFallback = -inf;
P.identification.vpLinearGapRatioFallback = -inf;
P.identification.diagnosticEO = [];

P.region.mode = 'catalog';
P.region.regionId = resonanceSelection.regionId;
P.region.shortTag = resonanceSelection.shortTag;
P.region.tag = resonanceSelection.tag;
P.region.regionStartSec = resonanceSelection.regionStartSec;
P.region.regionEndSec = resonanceSelection.regionEndSec;
P.region.analysisStartTimeSec = resonanceSelection.analysisStartTimeSec;
P.region.targetBladePasses = resonanceSelection.targetBladePasses;
P.region.windowBladePasses = resonanceSelection.windowBladePasses;
P.region.slidingStepBladePasses = resonanceSelection.slidingStepBladePasses;
P.region.dominantOrder = resonanceSelection.dominantOrder;
P.region.dominantFreqHz = resonanceSelection.dominantFreqHz;

P.waveform.windowLaps = resonanceSelection.windowBladePasses;
P.waveform.slidingStepLaps = resonanceSelection.slidingStepBladePasses;
P.waveform.targetLaps = resonanceSelection.targetBladePasses;
P.waveform.pulseWindowSec = syncCfg.pulse_window_sec;
P.waveform.pulsePadSec = syncCfg.pulse_pad_sec;
P.waveform.dynamicWindowMode = syncCfg.dynamic_window_mode; % peak_centered_fixed
P.waveform.windowRule = 'peak_centered_fixed: t_start=t_peak-pulseWindowSec; t_end=t_peak+pulseWindowSec';
P.waveform.saveInlineDynamicMap = true;

P.template.suffix = 'GradientXRange030_OPRCenterStd';

sensorTag = ['S', sprintf('%d', P.sensors.analysis)];
timeLabel = sprintf('%s_T%07.3fs', P.region.shortTag, P.region.analysisStartTimeSec);
timeLabel = regexprep(timeLabel, '[^A-Za-z0-9_\-]', '_');
caseTag = sprintf('B%d_%s', C.bladeId, sensorTag);

P.files.sensorConfig = fullfile(P.outputDir, '00_low_speed_sensor_config', ...
    'Sensor_Config_20251222.mat');
P.files.lowSpeedReference = fullfile(P.outputDir, '01_low_speed_reference', ...
    'LowSpeedReference_20251222.mat');
P.files.lowSpeedNumbering = fullfile(P.outputDir, '02_low_speed_numbering', ...
    'LowSpeedNumbering_20251222.mat');
P.files.regionSelection = fullfile(P.outputDir, '03_region_selection', ...
    sprintf('HighSpeedRegion_%s_20251222.mat', timeLabel));
P.files.peakCache = fullfile(P.outputDir, '04_high_speed_peak_cache', ...
    sprintf('HighSpeedPeakCache_%s_20251222.mat', sensorTag));
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensorTag, ...
    sprintf('HighSpeedNumbering_%s_20251222.mat', timeLabel));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05_low_speed_template_library', ...
    sprintf('LowSpeedTemplateLibrary_%s_%s_20251222.mat', caseTag, P.template.suffix));
P.files.inlineDynamicMap = fullfile(P.outputDir, '05_waveform_cache', sensorTag, ...
    sprintf('DynamicMap_%s_%s_InlineRaw_20251222.mat', caseTag, timeLabel));
P.files.identificationResult = fullfile(P.outputDir, '06_identification', ...
    sprintf('IdentificationResult_%s_%s_20251222.mat', caseTag, timeLabel));
P.files.identificationSummary = fullfile(P.outputDir, '06_identification', ...
    sprintf('IdentificationSummary_%s_%s_20251222.csv', caseTag, timeLabel));
P.files.auditSummary = fullfile(P.outputDir, '07_audit', ...
    sprintf('NumberingAudit_%s_%s_20251222.csv', caseTag, timeLabel));

P.name.resultSuffix = sprintf('%s_NewFlow_InlineRaw_Dx035_WithEta020_Reg005_PrevWinPhaseSafe', ...
    P.region.shortTag);
P.name.dynamicSuffix = sprintf('%s_InlineRaw_Main20L_W3S1_OPRCenterStd', P.region.shortTag);
end
