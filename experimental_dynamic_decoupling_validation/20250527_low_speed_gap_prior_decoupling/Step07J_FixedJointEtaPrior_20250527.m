function [eta, info] = Step07J_FixedJointEtaPrior_20250527(targetBlade, sensorIds)
%STEP07J_FIXEDJOINTETAPRIOR_20250527 Folder-local fixed eta prior for Step07J.

expectedBlade = 1;
expectedSensors = [1 3 6];
eta = [0, 0.0033808257, -0.095299809];

if nargin < 1 || isempty(targetBlade)
    targetBlade = expectedBlade;
end
if nargin < 2 || isempty(sensorIds)
    sensorIds = expectedSensors;
end

sensorIds = sensorIds(:).';
if double(targetBlade) ~= expectedBlade || ~isequal(sensorIds, expectedSensors)
    error(['20250527 Step07J local eta prior is defined only for B%d S%s. ' ...
        'Requested B%d S%s. Set STEP07J_STATIC_ETA_FILE for another case.'], ...
        expectedBlade, sprintf('%d', expectedSensors), ...
        double(targetBlade), sprintf('%d', sensorIds));
end

sourceFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '20250527_btt_data_foundation', 'output', ...
    'step05_joint_static_eta_preview', ...
    'JointStaticEtaPreview_B1_S136.mat');

info = struct();
info.file = mfilename('fullpath');
info.sourceFile = sourceFile;
info.sourceField = 'Summary.Consensus.eta_median_mm';
info.sourcePolicy = 'foundation fixed joint eta, frozen into gap folder';
info.targetBlade = expectedBlade;
info.sensorIds = expectedSensors;
info.sensorTag = 'S136';
info.etaMm = eta;
info.table = table(expectedSensors(:), eta(:), ...
    repmat(string(info.file), numel(expectedSensors), 1), ...
    repmat(string(sourceFile), numel(expectedSensors), 1), ...
    'VariableNames', {'sensor_id', 'eta_prior_mm', ...
    'local_prior_file', 'source_joint_static_eta_file'});
end
