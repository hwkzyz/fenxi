function Summary = Step04_Run_AllBlades_DirectTemplate_20241106()
%STEP04_RUN_ALLBLADES_DIRECTTEMPLATE_20241106 Run the direct route for B1-B6.
% This keeps the same Step01/Step02/Step03 workflow as the validated
% 20250527/20251222 direct-template route, but changes the target blade.

routeDir = fileparts(mfilename('fullpath'));
bladeIds = 1:6;
analysisSensors = [5 7];
regionId = 2;
regionStartMode = 'manual';
useRegionPlan = false;
manualStartTimeSec = 75.0;
calibrationLaps = 30;
sensorTag = ['S', sprintf('%d', analysisSensors)];
runLabel = sprintf('T%03.0fs', manualStartTimeSec);
dynamicSuffix = sprintf('Main20L_W3S1_%s_GradientXRange030_OPRCenterStd', runLabel);
resultSuffix = sprintf('Main_DirectTemplate_OPRCenterStd_%s_Dx035_NoEta_PrevWinPhaseSafe', runLabel);

summaryRows = repmat(empty_summary_row_local(), numel(bladeIds), 1);

fprintf('\n=== Step04: all-blade direct-template identification (20241106) ===\n');
fprintf('Route: %s\n', routeDir);
fprintf('Blades: %s, sensors: %s, regionId=%d, startMode=%s, manualStart=%.3f s\n', ...
    mat2str(bladeIds), mat2str(analysisSensors), regionId, regionStartMode, manualStartTimeSec);

for i = 1:numel(bladeIds)
    bladeId = bladeIds(i);
    fprintf('\n================ Blade %d/%d: B%d ================\n', ...
        i, numel(bladeIds), bladeId);

    ensure_stable_plan_local(routeDir, bladeId, analysisSensors, calibrationLaps);
    set_blade_environment_local(bladeId, analysisSensors, regionId, regionStartMode, ...
        useRegionPlan, manualStartTimeSec, dynamicSuffix, resultSuffix, calibrationLaps);

    run_script_local(fullfile(routeDir, 'Step01_Main_Build_OPRCenterStd_Template_20241106.m'));
    run_script_local(fullfile(routeDir, 'Step02_Main_Build_OPRCenterStd_DynamicMap_20241106.m'));
    run_script_local(fullfile(routeDir, 'Step03_Main_Run_OPRCenterStd_DirectTemplate_20241106.m'));

    summaryRows(i) = summarize_blade_result_local(routeDir, bladeId, sensorTag);
    fprintf('B%d summary: dominant EO=%d, median f=%.3f Hz, dx_c median=%.4f mm, A median=%.4f mm\n', ...
        bladeId, summaryRows(i).dominantEO, summaryRows(i).medianFreqHz, ...
        summaryRows(i).dxCMedianMM, summaryRows(i).AMedianMM);
end

Summary = struct2table(summaryRows);
outDir = fullfile(routeDir, 'output', 'identification');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
outFile = fullfile(outDir, sprintf('Step04_AllBlades_DirectSummary_%s_%s_20241106.csv', sensorTag, runLabel));
writetable(Summary, outFile);

fprintf('\n=== All-blade direct summary ===\n');
disp(Summary);
fprintf('Saved all-blade summary: %s\n', outFile);
end


function row = empty_summary_row_local()
row = struct( ...
    'BladeID', NaN, ...
    'SensorTag', '', ...
    'SelectionMode', '', ...
    'RegionID', NaN, ...
    'RegionStartSec', NaN, ...
    'PeakTimeSec', NaN, ...
    'DynamicStartSec', NaN, ...
    'DynamicEndSec', NaN, ...
    'DominantEOPlan', NaN, ...
    'dominantEO', NaN, ...
    'meanFreqHz', NaN, ...
    'medianFreqHz', NaN, ...
    'BestWindowID', NaN, ...
    'BestEO', NaN, ...
    'BestFreqHz', NaN, ...
    'BestAMM', NaN, ...
    'BestDxCMM', NaN, ...
    'BestWeightedRMSEV', NaN, ...
    'AMedianMM', NaN, ...
    'dxCMedianMM', NaN, ...
    'dxCMinMM', NaN, ...
    'dxCMaxMM', NaN, ...
    'dxCAtLimitCount', NaN, ...
    'CoordinateMaxDeltaMM', NaN);
end


function ensure_stable_plan_local(routeDir, bladeId, analysisSensors, calibrationLaps)
sensorTag = ['S', sprintf('%d', analysisSensors)];
planDir = fullfile(routeDir, 'output', 'stable_window_plan');
if exist(planDir, 'dir') ~= 7
    mkdir(planDir);
end
planFile = fullfile(planDir, sprintf( ...
    'Step01_CoverageFirstStableWindowPlan_B%d_%s_20241106.csv', bladeId, sensorTag));

BladeID = repmat(bladeId, numel(analysisSensors), 1);
SensorID = analysisSensors(:);
LapCount = repmat(calibrationLaps, numel(analysisSensors), 1);
StartOffsetLap = repmat(10, numel(analysisSensors), 1);
T = table(BladeID, SensorID, LapCount, StartOffsetLap);
writetable(T, planFile);
end


function set_blade_environment_local(bladeId, analysisSensors, regionId, regionStartMode, ...
    useRegionPlan, manualStartTimeSec, dynamicSuffix, resultSuffix, calibrationLaps)
setenv('STEP20241106_TARGET_BLADE', num2str(bladeId));
setenv('STEP20241106_ANALYSIS_SENSORS', sprintf('%d ', analysisSensors));
setenv('STEP20241106_REGION_ID', num2str(regionId));
setenv('STEP20241106_REGION_START_MODE', regionStartMode);
setenv('STEP20241106_USE_REGION_PLAN', logical_to_env_local(useRegionPlan));
setenv('STEP20241106_ANALYSIS_START_TIME', sprintf('%.9f', manualStartTimeSec));
setenv('STEP20241106_DYNAMIC_SUFFIX', dynamicSuffix);
setenv('STEP20241106_RESULT_SUFFIX', resultSuffix);
setenv('STEP20241106_CALIBRATION_LAPS', num2str(calibrationLaps));
end


function run_script_local(scriptPath)
[scriptDir, scriptName] = fileparts(scriptPath);
cmd = sprintf('matlab -batch "cd(''%s''); %s"', scriptDir, scriptName);
status = system(cmd);
if status ~= 0
    error('Script failed with status %d: %s', status, scriptPath);
end
end


function text = logical_to_env_local(value)
if value
    text = '1';
else
    text = '0';
end
end


function row = summarize_blade_result_local(routeDir, bladeId, sensorTag)
resultDir = fullfile(routeDir, 'output', 'identification');
pattern = sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s*_20241106.mat', bladeId, sensorTag);
files = dir(fullfile(resultDir, pattern));
if isempty(files)
    error('No Step03 result found for B%d %s in %s.', bladeId, sensorTag, resultDir);
end
[~, idxNewest] = max([files.datenum]);
resultFile = fullfile(files(idxNewest).folder, files(idxNewest).name);

loaded = load(resultFile, 'Result');
Result = loaded.Result;
Trend = Result.Trend;
loadedMap = load(Result.DynamicMapFile, 'DynamicMap');
DynamicMap = loadedMap.DynamicMap;

[~, bestIdx] = min(Trend.weighted_voltage_rmse);
best = Trend(bestIdx, :);

row = empty_summary_row_local();
row.BladeID = bladeId;
row.SensorTag = sensorTag;
if isfield(DynamicMap, 'SelectionInfo')
    row.SelectionMode = DynamicMap.SelectionInfo.start_mode;
    row.RegionID = DynamicMap.SelectionInfo.region_id;
    row.RegionStartSec = DynamicMap.SelectionInfo.region_start_sec;
    row.PeakTimeSec = DynamicMap.SelectionInfo.peak_time_sec;
    row.DominantEOPlan = DynamicMap.SelectionInfo.dominant_order;
end
row.DynamicStartSec = DynamicMap.GlobalTimeWindow(1);
row.DynamicEndSec = DynamicMap.GlobalTimeWindow(2);
row.dominantEO = Result.ResonanceSummary.dominant_eo;
row.meanFreqHz = Result.ResonanceSummary.mean_freq_hz;
row.medianFreqHz = Result.ResonanceSummary.median_freq_hz;
row.BestWindowID = best.window_id;
row.BestEO = best.EO_id;
row.BestFreqHz = best.fn_id;
row.BestAMM = best.A_id;
row.BestDxCMM = best.dx_c_id;
row.BestWeightedRMSEV = best.weighted_voltage_rmse;
row.AMedianMM = median(Trend.A_id, 'omitnan');
row.dxCMedianMM = median(Trend.dx_c_id, 'omitnan');
row.dxCMinMM = min(Trend.dx_c_id, [], 'omitnan');
row.dxCMaxMM = max(Trend.dx_c_id, [], 'omitnan');
row.dxCAtLimitCount = nnz(abs(Trend.dx_c_id) >= 0.349);
row.CoordinateMaxDeltaMM = Result.CoordinateCheck.max_abs_delta_mm;
end
