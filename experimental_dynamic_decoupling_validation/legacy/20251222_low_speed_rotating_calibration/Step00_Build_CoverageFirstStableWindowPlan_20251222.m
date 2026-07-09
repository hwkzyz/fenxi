%% Step00_Build_CoverageFirstStableWindowPlan_20251222
% Build the stable-window plan required by Step01 for the current CaseConfig.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
C = CaseConfig();
oldProgramDir = fullfile(routeDir, 'archive_previous_routes', 'old_programs');
legacyScript = fullfile(oldProgramDir, 'Step01_00_Build_CoverageFirstStableWindowPlan_20251222.m');

if exist(legacyScript, 'file') ~= 2
    error('Archived stable-window planning script not found: %s', legacyScript);
end

fprintf('\n=== Step00: coverage-first stable-window plan ===\n');
fprintf('Dataset: %s\n', C.dataset);
fprintf('Target blade: %d\n', C.bladeId);
fprintf('Sensors: %s\n', mat2str(C.sensorIds));

setenv('STEP01_TARGET_BLADE', num2str(C.bladeId));
setenv('STEP01_ANALYSIS_SENSORS', sprintf('%d ', C.sensorIds));
run(legacyScript);
setenv('STEP01_TARGET_BLADE', '');
setenv('STEP01_ANALYSIS_SENSORS', '');
