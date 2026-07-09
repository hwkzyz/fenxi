%% Main route for 20251222 low-speed gap-prior decoupling
% This is the single entry point for the current paper-facing workflow.
% Method comparison, TRAC and strain-spectrum checks are kept in separate
% comparison/validation programs.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

%% User settings
runFromStep = 1;
runToStep = 8;
stopOnError = true;
runFromStep = parse_numeric_env_local('FLOW_RUN_FROM_STEP', runFromStep);
runToStep = parse_numeric_env_local('FLOW_RUN_TO_STEP', runToStep);

steps = {
    1, 'Extract OPR blade timing',                 'Step01_Extract_OPR_Blade_Timing_20251222.m'
    2, 'Calculate full-time BTT displacement',     'Step02_Calc_FullTime_BTT_Displacement_20251222.m'
    3, 'Locate vibration by strain and OPR',       'Step03_Locate_Vibration_By_Strain_OPR_20251222.m'
    4, 'Detect resonance region by BTT-STE',       'Step04_Detect_Resonance_By_BTT_STE_20251222.m'
    5, 'Build low-speed response surface',         'Step05_Build_Response_Surface_20251222.m'
    5.1, 'Learn shared offset/tilt response',      'Step05I_Learn_OffsetTilt_Shared_Response_Surface_20251222.m'
    6, 'Build high-speed dynamic map',             'Step06_Build_HighMap_20251222.m'
    6.1, 'Calibrate low-speed gap library',        'Step06I_Calibrate_OffsetTiltShared_GapLibrary_20251222.m'
    7, 'Identify vibration with Step07J',          'Step07J_NestedStaticWarp_VPFullWave_20251222.m'
    8, 'Cross-validate Step07J static warp',       'Step08J_CrossValidate_StaticWarp_20251222.m'
    };

fprintf('\n=== 20251222 main low-speed gap-prior route ===\n');
fprintf('Folder: %s\n', thisDir);
fprintf('Step range: %.1f to %.1f\n\n', runFromStep, runToStep);

for i = 1:size(steps, 1)
    stepNo = steps{i, 1};
    stepTitle = steps{i, 2};
    stepFile = fullfile(thisDir, steps{i, 3});
    if stepNo < runFromStep || stepNo > runToStep
        continue;
    end
    if ~isfile(stepFile)
        error('Missing step file: %s', stepFile);
    end

    fprintf('\n--- Step %.1f: %s ---\n', stepNo, stepTitle);
    try
        run(stepFile);
    catch ME
        fprintf(2, 'Step %.1f failed: %s\n', stepNo, ME.message);
        if stopOnError
            rethrow(ME);
        end
    end
end

fprintf('\nMain route complete. Use Step09_Compare_Step07J_Methods_20251222.m for method comparison.\n');

function value = parse_numeric_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value)
    error('%s must be numeric.', name);
end
end
