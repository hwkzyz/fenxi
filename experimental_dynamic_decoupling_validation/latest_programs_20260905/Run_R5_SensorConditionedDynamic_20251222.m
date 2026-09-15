function out = Run_R5_SensorConditionedDynamic_20251222(sidecarFile, foundationFile, outputFile, targetBlade)
%RUN_R5_SENSORCONDITIONEDDYNAMIC_20251222 Unified R5 entry point.
% Static state is frozen in the sidecar; the common runner fits only the
% original dynamic variables and keeps the [300,1000] Hz contract.
rootDir = fileparts(mfilename('fullpath'));
commonDir = fullfile(fileparts(fileparts(rootDir)),'common');
addpath(rootDir,'-begin');
addpath(commonDir,'-begin');

cfg = Config_20251222();
if nargin < 4 || isempty(targetBlade)
    targetBlade = cfg.case.targetBlade;
end
targetBlade = double(targetBlade);
cfg.case.targetBlade = targetBlade;

if nargin < 1 || isempty(sidecarFile)
    sidecarFile = fullfile(cfg.paths.results, ...
        sprintf('R5_SensorConditionedSidecar_B%d_S123.mat',targetBlade));
end
if nargin < 2 || isempty(foundationFile)
    if targetBlade == 1
        foundationFile = fullfile(cfg.paths.preparedFoundation,'step05_r01', ...
            '1000_2500_3500', ...
            'Result_Step05_NoEtaVPTopKDirectTemplate_B1_S123_20251222.mat');
    elseif targetBlade == 5
        foundationFile = fullfile(cfg.paths.preparedFoundation,'step05_r04', ...
            '1000_2500_3500', ...
            'Result_Step05_NoEtaVPTopKDirectTemplate_B5_S123_20251222.mat');
    else
        error('R5:UnknownBlade','No default foundation is registered for B%d.',targetBlade);
    end
end
if nargin < 3 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results, ...
        sprintf('r5_sensor_conditioned_dynamic_B%d_result.mat',targetBlade));
end
assert(exist(sidecarFile,'file')==2,'R5:MissingSidecar', ...
    'Missing 20251222 sidecar: %s',sidecarFile);
assert(exist(foundationFile,'file')==2,'R5:MissingFoundation', ...
    'Missing 20251222 foundation: %s',foundationFile);
assert(exist(cfg.files.lowSpeedTemplate,'file')==2,'R5:MissingTemplate', ...
    'Missing 20251222 low-speed template: %s',cfg.files.lowSpeedTemplate);
assert(isequal(double(cfg.frequency.searchHz),[300 1000]), ...
    'R5:FrequencyContract','The formal 20251222 search band must remain [300,1000] Hz.');

out = R5_Run_SensorConditionedWindowIdentification(cfg, ...
    @R5_Build_SensorConditionedDynamicModels_20251222, ...
    sidecarFile,cfg.files.lowSpeedTemplate,foundationFile,outputFile);
end
