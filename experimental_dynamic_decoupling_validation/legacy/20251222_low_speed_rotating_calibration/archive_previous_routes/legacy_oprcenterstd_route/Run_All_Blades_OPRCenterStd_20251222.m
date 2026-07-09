%% Run_All_Blades_OPRCenterStd_20251222
% Batch runner for blades 1..6 using the single-point CaseConfig entry.
% It preserves the existing algorithm route and only automates the
% bootstrap/formal execution sequence for each blade.

clear; clc;

routeDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(routeDir);
oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(projectDir);

bladeIds = 1:6;
sensorIds = [1 2 3];
dataset = '20251222';
summaryRows = cell(numel(bladeIds), 8);

for i = 1:numel(bladeIds)
    bladeId = bladeIds(i);
    caseTag = sprintf('B%d_S%s', bladeId, sprintf('%d', sensorIds));
    summaryFile = fullfile(projectDir, 'output', 'identification', ...
        sprintf('Step04_DirectTemplate_OPRCenterStd_Summary_%s_%s.csv', dataset, caseTag));
    fprintf('\n==============================\n');
    fprintf('Running rotating OPRCenterStd route for %s\n', caseTag);
    fprintf('==============================\n');

    setenv('BLADE_CASE_BLADE_ID', num2str(bladeId));
    setenv('BLADE_CASE_SENSOR_IDS', sprintf('%d ', sensorIds));

    if ~isfile(summaryFile)
        run_reference_bootstrap_local(routeDir);
        run(fullfile(routeDir, 'Step00_Build_CoverageFirstStableWindowPlan_20251222.m'));
        run(fullfile(routeDir, 'Run_Main_OPRCenterStd_20251222.m'));
    else
        fprintf('Summary already exists, skipping recomputation:\n  %s\n', summaryFile);
    end

    T = readtable(summaryFile);
    summaryRows{i, 1} = dataset;
    summaryRows{i, 2} = bladeId;
    summaryRows{i, 3} = ['S', sprintf('%d', sensorIds)];
    summaryRows{i, 4} = T.DominantEO(1);
    summaryRows{i, 5} = T.EOConsistency(1);
    summaryRows{i, 6} = T.MeanFrequencyHz(1);
    summaryRows{i, 7} = T.MeanRMSE(1);
    summaryRows{i, 8} = T.MeanAmplitudeMM(1);
end

setenv('BLADE_CASE_BLADE_ID', '');
setenv('BLADE_CASE_SENSOR_IDS', '');

Summary = cell2table(summaryRows, 'VariableNames', ...
    {'Dataset','BladeId','SensorTag','DominantEO','EOConsistency', ...
     'MeanFrequencyHz','MeanRMSE','MeanAmplitudeMM'});
summaryFile = fullfile(projectDir, 'output', 'identification', ...
    'BatchSummary_AllBlades_OPRCenterStd_20251222.csv');
writetable(Summary, summaryFile);
disp(Summary);
fprintf('\nSaved batch summary:\n  %s\n', summaryFile);

function run_reference_bootstrap_local(routeDir)
setenv('STEP01_ALLOW_REFERENCE_ONLY', '1');
projectDir = fileparts(routeDir);
run(fullfile(projectDir, 'Step01_Main_Build_OPRCenterStd_Template_20251222.m'));

setenv('STEP01_ALLOW_REFERENCE_ONLY', '');
setenv('STEP02_TEMPLATE_SUFFIX', 'OPRCenterStdRef');
setenv('STEP02_DYNAMIC_SUFFIX', 'BootstrapFromRef');
run(fullfile(projectDir, 'Step02_Main_Build_OPRCenterStd_DynamicMap_20251222.m'));

setenv('STEP01_ALLOW_REFERENCE_ONLY', '');
setenv('STEP02_TEMPLATE_SUFFIX', '');
setenv('STEP02_DYNAMIC_SUFFIX', '');
end
