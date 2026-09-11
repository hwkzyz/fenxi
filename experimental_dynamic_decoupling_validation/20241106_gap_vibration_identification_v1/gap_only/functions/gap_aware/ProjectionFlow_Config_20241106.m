function P = ProjectionFlow_Config_20241106()
%PROJECTIONFLOW_CONFIG_20241106 Central config for the all-blade bank route.
%
% The 20241106 route separates low-speed templates from gap calibration:
%   1) read single blade-sensor low-speed templates for all blades and sensors [2 3 5 7],
%   2) calibrate a gap-correction bank only for capacitive sensors [5 7],
%   3) identify any requested blade by selecting its high-speed windows and
%      loading only the required template/gap entries.

P = struct();
P.dataset = '20241106';

P.machine.bladeCount = 6;
P.machine.oprPulsesPerRev = 6;
P.machine.oprNote = ['This dataset stores six blade-passing slots per rotor ' ...
    'revolution; one fixed target blade contributes one pass per revolution.'];

P.calibration.buildBladeIds = 1:P.machine.bladeCount;
P.calibration.lowSpeedTemplateSensors = [2 3 5 7];
P.calibration.gapSensorIds = [5 7];
P.calibration.directOnlySensorIds = [2 3];
P.calibration.templateBankDirName = 'low_speed_template_bank';
P.calibration.templateBankFile = 'LowSpeedTemplateBank_20241106_B1toB6_S2357.mat';
P.calibration.runtimeTemplateBankFile = 'LowSpeedTemplateBank_20241106_B4_S2357.mat';
P.calibration.gapBankDirName = 'gap_calibration_bank';
P.calibration.gapBankFile = 'GapCalibrationBank_20241106_B1toB6_S57.mat';
P.calibration.runtimeGapBankFile = 'GapCalibrationBank_20241106_B4_S57.mat';
P.calibration.gapEntryPrefix = 'GapCalib';
P.calibration.templateBundleSensorTag = 'S2357';
P.calibration.templateSuffixPreference = 'GradientXRange030_OPRCenterStd';
P.calibration.singleTemplateDirName = 'template_library';
P.calibration.singleTemplatePattern = 'Template_20241106_B%d_S%d.mat';

P.identification.targetBlade = 4;
P.identification.analysisSensors = [2 5 7];
P.identification.gapSensors = [5 7];
P.identification.directOnlySensors = [2];
P.identification.analysisStartTimeSec = 75.0;
P.identification.targetBladePasses = 20;
P.identification.windowBladePasses = 3;
P.identification.slidingStepBladePasses = 1;
P.identification.maxRunWindows = inf;
P.identification.pulseWindowSec = 6e-4;
P.identification.freqSearchHz = [300 1000];
P.run.mode = 'identify';       % calibration, identify, compare, or all
P.run.resultSuffix = '';       % leave empty for Step07J default suffix

P.model.deltaGapLimitMm = 0.25;
P.model.deltaMuLimit = 0.060;
P.model.deltaTauLimitMm = 0.16;

% Compatibility adapter: the root package config is the only formal source
% of case, window, frequency and model settings.
C = Config_20241106();
P.identification.targetBlade = C.case.targetBlade;
P.identification.analysisSensors = C.case.analysisSensors;
P.identification.gapSensors = C.case.gapSensors;
P.identification.directOnlySensors = setdiff(C.case.analysisSensors, ...
    C.case.gapSensors, 'stable');
P.identification.analysisStartTimeSec = C.case.analysisStartTimeSec;
P.identification.targetBladePasses = C.window.targetBladePasses;
P.identification.windowBladePasses = C.window.windowBladePasses;
P.identification.slidingStepBladePasses = C.window.slidingStepBladePasses;
P.identification.freqSearchHz = C.frequency.searchHz;
P.model.deltaGapLimitMm = C.model.deltaGapLimitMm;
P.model.deltaMuLimit = C.model.deltaMuLimit;
P.model.deltaTauLimitMm = C.model.deltaTauLimitMm;
P.package = C;
end
