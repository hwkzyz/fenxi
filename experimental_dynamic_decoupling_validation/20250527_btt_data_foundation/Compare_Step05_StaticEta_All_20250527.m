clc; close all;

%COMPARE_STEP05_STATICETA_ALL_20250527
% Controlled static-eta diagnostic wrapper.
%
% Runs:
%   1) fixed Step04-template-center eta
%   2) bounded Step04-template-center eta, default +/-0.02 mm
%
% The formal Step05 script is not modified by this wrapper. Outputs are
% written under output/step05_static_eta_compare.

old_mode = getenv('STEP05_STATIC_ETA_MODE');
old_halfwidth = getenv('STEP05_STATIC_ETA_HALFWIDTH_MM');
old_show = getenv('STEP05_SHOW_PLOTS');
old_save = getenv('STEP05_SAVE_FIGURES');

cleanup_obj = onCleanup(@() restore_env_local( ...
    old_mode, old_halfwidth, old_show, old_save));

setenv('STEP05_SHOW_PLOTS', '0');
setenv('STEP05_SAVE_FIGURES', '0');

if isempty(strtrim(old_halfwidth))
    setenv('STEP05_STATIC_ETA_HALFWIDTH_MM', '0.02');
end

run_modes = {'fixed', 'bounded'};

for i = 1:numel(run_modes)
    setenv('STEP05_STATIC_ETA_MODE', run_modes{i});
    fprintf('\n\n===== Static eta run %d/%d: %s =====\n', ...
        i, numel(run_modes), run_modes{i});
    Compare_Step05_StaticEta_20250527;
end


function restore_env_local(old_mode, old_halfwidth, old_show, old_save)
setenv('STEP05_STATIC_ETA_MODE', old_mode);
setenv('STEP05_STATIC_ETA_HALFWIDTH_MM', old_halfwidth);
setenv('STEP05_SHOW_PLOTS', old_show);
setenv('STEP05_SAVE_FIGURES', old_save);
end
