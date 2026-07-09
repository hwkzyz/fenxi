clc; clear; close all;

%RUN_STEP05_VPTOP3_BOUNDEDRMSE_20251222
% Convenience runner for the 20251222 VP top-3 bounded-RMSE direct-template
% route. The implementation remains in the single official Step05 script;
% this file only sets a reproducible method preset and output tag.
%
% Optional environment overrides before running this script:
%   STEP05_TARGET_BLADES      e.g. '1' or '1 2 3'
%   STEP05_TARGET_LAPS        e.g. '20'
%   STEP05_WINDOW_LAPS        e.g. '3'
%   STEP05_DEBUG_MAX_WINDOWS  e.g. '3' for a quick smoke run

setenv('STEP05_METHOD_PRESET', 'vp_top3_bounded_rmse');
if strlength(strtrim(getenv('STEP05_OUTPUT_TAG'))) == 0
    setenv('STEP05_OUTPUT_TAG', 'vp3_bounded');
end
if strlength(strtrim(getenv('STEP05_FORCE_REBUILD'))) == 0
    setenv('STEP05_FORCE_REBUILD', '1');
end

fprintf('\n=== 20251222 Step05 VP top-3 bounded-RMSE route ===\n');
fprintf('Output tag: %s\n', getenv('STEP05_OUTPUT_TAG'));
fprintf('Target blades env: %s\n', getenv('STEP05_TARGET_BLADES'));
fprintf('Target laps env: %s\n', getenv('STEP05_TARGET_LAPS'));
fprintf('Window laps env: %s\n', getenv('STEP05_WINDOW_LAPS'));
fprintf('Debug windows env: %s\n', getenv('STEP05_DEBUG_MAX_WINDOWS'));

Step05_SingleSync_DirectTemplate_Identification_20251222;
