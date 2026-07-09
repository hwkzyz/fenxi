function P = NewFlow_Config_20241106()
%NEWFLOW_CONFIG_20241106 Central tuning block for the new numbering workflow.
% Edit this file first. Each StepXX script should read parameters from here
% and write one clear artifact for the next step.

routeDir = fileparts(mfilename('fullpath'));
legacyDir = fullfile(fileparts(routeDir), '20241106_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacyDir, 'dir') ~= 7
    error('Legacy helper folder not found: %s', legacyDir);
end
pulsePadSec = 2 / 5e6;

P = struct();
P.routeDir = routeDir;
P.outputDir = fullfile(routeDir, 'output', 'new_flow');

% Raw data and pre-extracted timing products.
P.data.datasetRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
P.data.lowSpeedCase = '900';
P.data.highSpeedCase = '3000_3150';
P.data.lowSpeedDir = fullfile(P.data.datasetRoot, P.data.lowSpeedCase);
P.data.highSpeedDir = fullfile(P.data.datasetRoot, P.data.highSpeedCase);
P.data.legacyOutputRoot = fullfile(legacyDir, 'output');
P.data.lowSpeedReferenceDir = fullfile(routeDir, 'output', 'reference', [P.data.lowSpeedCase, '_reference']);
P.data.highSpeedPulseDir = fullfile(P.data.legacyOutputRoot, P.data.highSpeedCase);

% Channels and machine constants.
P.machine.bladeCount = 6;
P.machine.oprChannel = 1;
P.machine.oprPulsesPerRev = 1;
P.machine.sampleRateHz = 5e6;
P.machine.tipRadiusMM = 65.0;

% Main analysis subset. Numbering is done once for this subset and region.
P.sensors.analysis = [2 3 5 7];
P.sensors.anchorPreference = [5 7 2 3];

% Step03: high-speed region selection.
P.region.mode = 'manual_start';       % 'manual_start' | 'explicit_range' | 'region_plan'
P.region.startTimeSec = 75.0;
P.region.timeRangeSec = [74.998 75.369];
P.region.targetLaps = 20;
P.region.regionId = 2;
P.region.startMode = 'peak';          % only for region_plan: 'peak' | 'start'
P.region.planFile = fullfile(routeDir, 'outputs', 'Step05_BTT_WindowPlan_20241106.csv');

% Step04: high-speed numbering.
P.numbering.seedLaps = 1;
P.numbering.maxCandidateLaps = 5;
P.numbering.matchPolyDegree = 7;
P.numbering.minSeedScoreMargin = 0.02;

% Step06: high-speed waveform-map extraction.
P.waveform.windowLaps = 3;
P.waveform.slidingStepLaps = 1;
P.waveform.pulseWindowSec = NaN;
P.waveform.pulsePadSec = pulsePadSec;
P.waveform.dynamicWindowMode = 'legacy_row_bounds';
P.waveform.source = 'jilublade_row_bounds_from_legacy_timing';
P.waveform.windowRule = 't_start=jilublade(row,1)-pulsePadSec; t_end=jilublade(row,2)+pulsePadSec';

% Step05/Step06: low-speed template calibration and identification settings.
P.template.suffix = 'GradientXRange030_OPRCenterStd';
P.template.blades = 1:P.machine.bladeCount;
P.template.defaultCalibrationLaps = 30;
P.template.templateGridN = 1201;
P.template.splineSmoothing = 0.995;
P.template.weightFloor = 0.05;
P.template.trustQuantile = 0.995;
P.template.trustEdgeMarginMM = 0.02;
P.template.segmentExpandFactor = 0.30;
P.template.gradientEnergyQuantile = 0.995;
P.template.gradientEdgeMarginMM = 0.02;
P.template.xrangeMode = 'threshold';   % 'threshold' | 'energy'
P.template.gradientMinRatio = 0.30;
P.template.amplitudeMinRatio = 0;
P.template.minHalfWidthMM = 2.5;
P.template.maxHalfWidthMM = 4.2;
P.template.sensorDomainExpandMMTable = [];
P.template.centerMode = 'sgfit';       % 'centroid' | 'sgfit'
P.template.preferStableWindowPlan = true;
P.template.minSelectedPulseCount = 3;
P.template.minValidTemplateBins = 30;
P.identification.blades = 1:P.machine.bladeCount;
P.identification.dxCLimitMM = 0.35;
P.identification.amplitudeLimitMM = 0.50;
P.identification.freqSearchHz = [0 1000];
P.identification.eoPad = 2;
P.identification.topKEO = 3;
P.identification.sensorEtaLimitMM = 0;
P.identification.sensorEtaRegWeightVPerMM = 0;
P.identification.overshootPenaltyWeight = 50;
P.identification.pulseSelectionMode = 'single'; % old Step03: 'single' main pulse or 'all'
P.identification.domainSelectionMode = 'hard';
P.identification.domainMarginMM = 0.02;
P.identification.domainSoftMarginMM = 0.10;
P.identification.coverageSafetyMarginMM = 0.05;
P.identification.queryGuardMode = 'adaptive';
P.identification.queryGuardMM = P.identification.amplitudeLimitMM + ...
    P.identification.dxCLimitMM + P.identification.coverageSafetyMarginMM;
P.identification.queryGuardQuantile = 95;
P.identification.queryGuardSafetyMM = 0.05;
P.identification.queryGuardMinMM = 0.12;
P.identification.queryGuardMaxMM = P.identification.queryGuardMM;
P.identification.sensorDomainExpandMMTable = [];
P.identification.sensorQueryGuardScaleTable = [];
P.identification.dynamicEffectiveMode = 'gradient';
P.identification.dynamicTemplateGradientMinRatio = 0.08;
P.identification.dynamicTimeGradientMinRatio = 0.15;
P.identification.dynamicPeakQuantile = 85;
P.identification.defaultSensorThreshold = 0.5;
P.identification.weightFloor = 0.05;
P.identification.vibrationRefineEnable = true;
P.identification.vibrationRefineMaxPass = 1;
P.identification.vibrationRefineMinPointCount = 20;
P.identification.vibrationRefineMinPointCountPerSensor = 6;
P.identification.vibrationRefineDomainMarginMM = 0.02;
P.identification.vibrationRefineUseQueryGuard = true;
P.identification.vibrationRefineUseCompensatedGradient = true;

% Visualization. Fast steps draw figures directly; long steps have separate
% Step04V/Step05V visualizers that read saved artifacts.
P.view.enable = true;
P.view.saveFigures = true;
P.view.figureDir = fullfile(P.outputDir, 'figures');

% Artifact paths.
% Script semantics now follow:
%   Step05 = low-speed template calibration
%   Step06 = direct-template identification with inline high-speed extraction
%   Step07 = audit
% The saved folder names below still keep the verified historical numbering
% (05A/05/06/07) so existing outputs remain readable without migration.
P.files.sensorConfig = fullfile(P.data.lowSpeedReferenceDir, 'Sensor_Config_20241106.mat');
P.files.lowSpeedFingerprint = fullfile(P.outputDir, '01_low_speed_reference', 'LowSpeedReferenceFingerprint_20241106.mat');
P.files.lowSpeedNumbering = fullfile(P.outputDir, '02_low_speed_numbering', sensor_tag_local(P.sensors.analysis), ...
    'LowSpeedNumbering_20241106.mat');
P.files.regionSelection = fullfile(P.outputDir, '03_region_selection', ...
    sprintf('HighSpeedRegion_%s_20241106.mat', time_label_local(P.region.startTimeSec)));
P.files.highSpeedNumbering = fullfile(P.outputDir, '04_high_speed_numbering', sensor_tag_local(P.sensors.analysis), ...
    sprintf('HighSpeedNumbering_%s_20241106.mat', time_label_local(P.region.startTimeSec)));
P.files.waveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensor_tag_local(P.sensors.analysis), ...
    sprintf('WaveformLibrary_%s_20241106.mat', time_label_local(P.region.startTimeSec)));
P.files.filteredWaveformLibrary = fullfile(P.outputDir, '05_waveform_library', sensor_tag_local(P.sensors.analysis), ...
    sprintf('FilteredWaveformLibrary_%s_20241106.mat', time_label_local(P.region.startTimeSec)));
P.files.filteredWaveformSummary = fullfile(P.outputDir, '05_waveform_library', sensor_tag_local(P.sensors.analysis), ...
    sprintf('FilteredWaveformSummary_%s_20241106.csv', time_label_local(P.region.startTimeSec)));
P.files.lowSpeedTemplateLibrary = fullfile(P.outputDir, '05A_low_speed_template_library', sensor_tag_local(P.sensors.analysis), ...
    'LowSpeedTemplateLibrary_20241106.mat');
P.files.lowSpeedTemplateSummary = fullfile(P.outputDir, '05A_low_speed_template_library', sensor_tag_local(P.sensors.analysis), ...
    'LowSpeedTemplateLibrary_20241106.csv');
P.files.lowSpeedTemplateAudit = fullfile(P.outputDir, '05A_low_speed_template_library', sensor_tag_local(P.sensors.analysis), ...
    'LowSpeedTemplateAudit_20241106.mat');
P.files.identificationSummary = fullfile(P.outputDir, '06_identification', sensor_tag_local(P.sensors.analysis), ...
    sprintf('IdentificationSummary_%s_20241106.csv', time_label_local(P.region.startTimeSec)));
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensor_tag_local(P.sensors.analysis), ...
    sprintf('IdentificationResult_%s_20241106.mat', time_label_local(P.region.startTimeSec)));
P.files.auditReport = fullfile(P.outputDir, '07_audit', sensor_tag_local(P.sensors.analysis), ...
    sprintf('NumberingAudit_%s_20241106.csv', time_label_local(P.region.startTimeSec)));
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end
