%% Step07J debug launcher
% Edit only the "User edit area" below, then run this file.
% This launcher sets runtime parameters and then calls the unchanged
% Step07J identification script.

clear; clc;

%% User edit area
regionId = 2;              % 1=R01, 2=R02, ..., 7=R07
analysisStartTimeSec = []; % [] uses ResonanceRegionCatalog default; e.g. 20.10
targetLaps = 20;           % total selected blade passes
windowLaps = 3;            % blade passes per identification window
slidingStepLaps = 1;       % sliding step in blade passes
maxWindows = 10;            % use inf for all windows

eoSearchMode = 'region_order'; % 'region_order' or 'frequency_range'
eoHalfWidth = 10;              % 10 for exploration; 2 for near-region check; 0 for center EO only
freqSearchHz = [];             % [] keeps default; e.g. [500 650]

resultSuffix = ''; % []/'' auto-generates from the settings below

%% Runtime setup
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

[regionSelection, ~] = ResonanceRegionCatalog_20251222(regionId);
if isempty(resultSuffix)
    if isempty(analysisStartTimeSec)
        startTag = sprintf('T%07.3f', regionSelection.analysisStartTimeSec);
    else
        startTag = sprintf('T%07.3f', analysisStartTimeSec);
    end
    startTag = strrep(startTag, '.', 'p');
    resultSuffix = sprintf('%s_debug_%s_L%d_W%d_S%d', ...
        regionSelection.shortTag, startTag, targetLaps, windowLaps, slidingStepLaps);
end

clear_step07j_runtime_env_local();
setenv('BLADE_RESONANCE_REGION_ID', num2str(regionId));
set_or_clear_env_local('BLADE_ANALYSIS_START_TIME_SEC', analysisStartTimeSec);
setenv('BLADE_TARGET_LAPS', num2str(targetLaps));
setenv('BLADE_WINDOW_LAPS', num2str(windowLaps));
setenv('BLADE_SLIDING_STEP_LAPS', num2str(slidingStepLaps));
setenv('STEP07J_MAX_WINDOWS', num2str(maxWindows));
setenv('STEP07J_EO_SEARCH_MODE', eoSearchMode);
setenv('STEP07J_EO_SEARCH_HALF_WIDTH', num2str(eoHalfWidth));
set_or_clear_env_local('STEP07J_FREQ_SEARCH_HZ', freqSearchHz);
setenv('STEP07J_RESULT_SUFFIX', resultSuffix);

fprintf('\n=== Step07J debug launcher ===\n');
fprintf('Region: %s, catalog start %.4f s, EO center %d\n', ...
    regionSelection.shortTag, regionSelection.analysisStartTimeSec, ...
    regionSelection.dominantOrder);
if isempty(analysisStartTimeSec)
    fprintf('Analysis start: catalog default\n');
else
    fprintf('Analysis start override: %.4f s\n', analysisStartTimeSec);
end
fprintf('Laps/window/step/maxWindows: %d/%d/%d/%s\n', ...
    targetLaps, windowLaps, slidingStepLaps, num2str(maxWindows));
fprintf('EO search: %s, half-width %d\n', eoSearchMode, eoHalfWidth);
if ~isempty(freqSearchHz)
    fprintf('Frequency search override: [%.3f %.3f] Hz\n', ...
        freqSearchHz(1), freqSearchHz(2));
end
fprintf('Result suffix: %s\n', resultSuffix);

Step07J_NestedStaticWarp_VPFullWave_20251222;

function set_or_clear_env_local(name, value)
if isempty(value)
    setenv(name, '');
elseif isnumeric(value)
    if isscalar(value)
        setenv(name, num2str(value));
    else
        setenv(name, sprintf('%g %g', value(1), value(2)));
    end
else
    setenv(name, char(value));
end
end

function clear_step07j_runtime_env_local()
names = {
    'BLADE_RESONANCE_REGION_ID'
    'BLADE_ANALYSIS_START_TIME_SEC'
    'BLADE_TARGET_LAPS'
    'BLADE_WINDOW_LAPS'
    'BLADE_SLIDING_STEP_LAPS'
    'BLADE_PULSE_WINDOW_SEC'
    'STEP07J_MAX_WINDOWS'
    'STEP07J_EO_SEARCH_MODE'
    'STEP07J_EO_SEARCH_HALF_WIDTH'
    'STEP07J_FREQ_SEARCH_HZ'
    'STEP07J_RESULT_SUFFIX'
    'STEP07J_DYNAMIC_MAP_FILE'
    'STEP07J_HIGHMAP_FILE'
    'STEP07J_DIRECT_RESULT_FILE'
    };
for i = 1:numel(names)
    setenv(names{i}, '');
end
end
