%% Step03_Run_Current_S13_S136_20250527
% Current low-speed-template-only identification entry.
%
% Scope:
%   - Only run S13 and S136.
%   - Do not run the archived EO-lock / old-SG-comparison Step03 route.
%   - Use direct low-speed-template waveform optimization as the primary program.

clear; clc;

route_dir = fileparts(mfilename('fullpath'));
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(route_dir);

sensor_sets = {
    '1 3'
    '1 3 6'
};

primary_script = 'Step03_01_Run_DirectLowSpeedWaveform_20250527.m';

fprintf('\n=== Current Step03 batch: S13 and S136 only ===\n');
fprintf('Primary script: %s\n', primary_script);
fprintf('Archived EO-lock Step03 scripts are intentionally not run.\n');

for k = 1:numel(sensor_sets)
    run_step03_script_local(primary_script, sensor_sets{k});
end

setenv('STEP03_ANALYSIS_SENSORS', '');
fprintf('\nCompleted current Step03 batch for S13 and S136.\n');

function run_step03_script_local(primary_script, sensor_tag)
setenv('STEP03_ANALYSIS_SENSORS', sensor_tag);
fprintf('\n--- Running %s with STEP03_ANALYSIS_SENSORS=\"%s\" ---\n', ...
    primary_script, sensor_tag);
run(primary_script);
end
