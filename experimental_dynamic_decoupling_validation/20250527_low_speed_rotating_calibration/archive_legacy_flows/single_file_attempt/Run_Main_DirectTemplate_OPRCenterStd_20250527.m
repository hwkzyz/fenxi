%% Run_Main_DirectTemplate_OPRCenterStd_20250527
% Single-file latest workflow for the 20250527 low-speed rotating calibration.
% Tune the parameters in the block below, then run this file directly.
% No static gap library is used.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));

%% Parameters to tune
P.run.buildTemplate = true;
P.run.buildDynamicMap = true;
P.run.identify = true;
P.run.visualize = true;

P.name.referenceTemplateSuffix = 'OPRCenterStdRef';
P.name.templateSuffix = 'GradientXRange030_OPRCenterStd';
P.name.dynamicSuffix = 'Main20L_W3S1_GradientXRange030_OPRCenterStd';
P.name.resultSuffix = 'Main_DirectTemplate_OPRCenterStd';

P.step01.showLegacyPlots = true;
P.step01.forceLegacyRebuild = false;
P.step01.forceOprCenterRebuild = false;
P.step01.centerMode = 'sgfit';
P.step01.xrangeMode = 'threshold';
P.step01.gradientMinRatio = 0.30;
P.step01.amplitudeMinRatio = 0;
P.step01.minHalfWidthMM = 2.5;
P.step01.maxHalfWidthMM = 4.2;
P.step01.referenceCenterMode = 'sgfit';
P.step01.referenceXrangeMode = 'energy';
P.step01.rebuildReferenceTemplate = false;

P.step02.targetLaps = 20;
P.step02.windowLaps = 3;
P.step02.slidingStepLaps = 1;

P.step03.analysisSensors = [1 3 6];
P.step03.topKEO = 3;
P.step03.dynamicEffectiveMode = 'gradient';
P.step03.pulseMode = 'single';
P.step03.domainSelectionMode = 'hard';
P.step03.domainSoftMarginMM = 0;
P.step03.queryGuardMode = 'adaptive';
P.step03.queryGuardQuantile = 95;
P.step03.queryGuardSafetyMM = 0.05;
P.step03.queryGuardMinMM = 0.12;
P.step03.sensorEtaLimitMM = 0.03;
P.step03.sensorEtaRegWeightVPerMM = 0;
P.step03.overshootPenaltyWeight = 100;

%% Fixed file names
stablePlanFile = fullfile(routeDir, 'output', 'stable_window_plan', ...
    'Step01_CoverageFirstStableWindowPlan_B1_S136_20250527.csv');
referenceCenterTemplateFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('Template_LowSpeedRotating_B1_S136_%s_20250527.mat', P.name.referenceTemplateSuffix));
templateFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('Template_LowSpeedRotating_B1_S136_%s_20250527.mat', P.name.templateSuffix));
dynamicFile = fullfile(routeDir, 'output', 'dynamic_maps', ...
    sprintf('DynamicMap_B1_S136_SlidingWindows_%s_20250527.mat', P.name.dynamicSuffix));

%% Run latest workflow
if P.run.buildTemplate
    fprintf('\n=== Step01: OPR-center timing and low-speed template ===\n');
    apply_step01_preprocess_settings(P);
    step01_extract_opr_blade_timing_embedded(routeDir, P);

    if exist(stablePlanFile, 'file') ~= 2
        error(['Stable-window plan not found:\n  %s\n' ...
            'The old planning script is archived under archive_legacy_flows/old_programs.'], ...
            stablePlanFile);
    end

    if P.step01.rebuildReferenceTemplate || exist(referenceCenterTemplateFile, 'file') ~= 2
        fprintf('Building reference-center template first...\n');
        apply_step01_reference_template_settings(P);
        step01_build_rotating_template_embedded(routeDir);
    end

    apply_step01_main_template_settings(P, stablePlanFile, referenceCenterTemplateFile);
    step01_build_rotating_template_embedded(routeDir);
end

if P.run.buildDynamicMap
    fprintf('\n=== Step02: dynamic sliding-window map ===\n');
    if exist(templateFile, 'file') ~= 2
        error('Template file not found. Run Step01 first: %s', templateFile);
    end
    apply_step02_settings(P);
    step02_build_dynamic_map_embedded(routeDir);
end

if P.run.identify
    fprintf('\n=== Step03: direct-template identification ===\n');
    if exist(templateFile, 'file') ~= 2
        error('Template file not found. Run Step01 first: %s', templateFile);
    end
    if exist(dynamicFile, 'file') ~= 2
        error('Dynamic map file not found. Run Step02 first: %s', dynamicFile);
    end
    apply_step03_settings(P, templateFile, dynamicFile);
    step03_identify_direct_template_embedded(routeDir);
end

if P.run.visualize
    fprintf('\n=== Step04: summary and trend figure ===\n');
    step04_visualize_direct_template_embedded(routeDir, P);
end

%% Parameter adapters
function apply_step01_preprocess_settings(P)
setenv('STEP01_FORCE_LEGACY_REBUILD', logical_to_env_local(P.step01.forceLegacyRebuild));
setenv('STEP01_FORCE_OPR_CENTER_REBUILD', logical_to_env_local(P.step01.forceOprCenterRebuild));
end

function apply_step01_reference_template_settings(P)
setenv('STEP01_STABLE_WINDOW_PLAN', '');
setenv('STEP01_TEMPLATE_SUFFIX', P.name.referenceTemplateSuffix);
setenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE', '');
setenv('STEP01_CENTER_MODE', P.step01.referenceCenterMode);
setenv('STEP01_XRANGE_MODE', P.step01.referenceXrangeMode);
end

function apply_step01_main_template_settings(P, stablePlanFile, referenceCenterTemplateFile)
setenv('STEP01_STABLE_WINDOW_PLAN', stablePlanFile);
setenv('STEP01_TEMPLATE_SUFFIX', P.name.templateSuffix);
setenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE', referenceCenterTemplateFile);
setenv('STEP01_CENTER_MODE', P.step01.centerMode);
setenv('STEP01_XRANGE_MODE', P.step01.xrangeMode);
setenv('STEP01_XRANGE_GRADIENT_MIN_RATIO', num2str(P.step01.gradientMinRatio));
setenv('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', num2str(P.step01.amplitudeMinRatio));
setenv('STEP01_XRANGE_MIN_HALF_WIDTH_MM', num2str(P.step01.minHalfWidthMM));
setenv('STEP01_XRANGE_MAX_HALF_WIDTH_MM', num2str(P.step01.maxHalfWidthMM));
end

function apply_step02_settings(P)
setenv('STEP02_TEMPLATE_SUFFIX', P.name.templateSuffix);
setenv('STEP02_DYNAMIC_SUFFIX', P.name.dynamicSuffix);
setenv('STEP02_TARGET_LAPS', num2str(P.step02.targetLaps));
setenv('STEP02_ANALYSIS_WIN_SIZE', num2str(P.step02.windowLaps));
setenv('STEP02_SLIDING_STEP', num2str(P.step02.slidingStepLaps));
end

function apply_step03_settings(P, templateFile, dynamicFile)
setenv('STEP03_TEMPLATE_FILE', templateFile);
setenv('STEP03_DYNAMIC_MAP_FILE', dynamicFile);
setenv('STEP03_RESULT_SUFFIX', P.name.resultSuffix);
setenv('STEP03_ANALYSIS_SENSORS', sprintf('%d ', P.step03.analysisSensors));
setenv('STEP03_MAIN_TOP_K_EO', num2str(P.step03.topKEO));
setenv('STEP03D_DYNAMIC_EFFECTIVE_MODE', P.step03.dynamicEffectiveMode);
setenv('STEP03D_PULSE_MODE', P.step03.pulseMode);
setenv('STEP03D_DOMAIN_SELECTION_MODE', P.step03.domainSelectionMode);
setenv('STEP03D_DOMAIN_SOFT_MARGIN_MM', num2str(P.step03.domainSoftMarginMM));
setenv('STEP03D_QUERY_GUARD_MODE', P.step03.queryGuardMode);
setenv('STEP03D_QUERY_GUARD_QUANTILE', num2str(P.step03.queryGuardQuantile));
setenv('STEP03D_QUERY_GUARD_SAFETY_MM', num2str(P.step03.queryGuardSafetyMM));
setenv('STEP03D_QUERY_GUARD_MIN_MM', num2str(P.step03.queryGuardMinMM));
setenv('STEP03D_SENSOR_ETA_LIMIT_MM', num2str(P.step03.sensorEtaLimitMM));
setenv('STEP03D_SENSOR_ETA_REG_WEIGHT_V_PER_MM', num2str(P.step03.sensorEtaRegWeightVPerMM));
setenv('STEP03D_OVERSHOOT_PENALTY_WEIGHT', num2str(P.step03.overshootPenaltyWeight));
end

function s = logical_to_env_local(tf)
if tf
    s = '1';
else
    s = '0';
end
end

function step04_visualize_direct_template_embedded(routeDir, P)
resultDir = fullfile(routeDir, 'output', 'identification');
figureDir = fullfile(routeDir, 'output', 'figures', 'main_direct_template_oprcenterstd');
if exist(figureDir, 'dir') ~= 7; mkdir(figureDir); end

resultFile = fullfile(resultDir, sprintf( ...
    'Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_%s_20250527.mat', P.name.resultSuffix));
if exist(resultFile, 'file') ~= 2
    error('Direct OPRCenterStd Step03 result not found: %s', resultFile);
end

loaded = load(resultFile, 'Result');
Result = loaded.Result;
T = Result.Trend;

fig = figure('Name', '20250527 direct OPRCenterStd identification trend', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 17, 11]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(T.window_id, T.EO_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5); hold on;
yline(mode(T.EO_id), ':', 'LineWidth', 1.0);
ylabel('EO');
title('Engine order');
style_axes_for_step04_local();

nexttile;
plot(T.window_id, T.fn_id, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('Frequency (Hz)');
title('Identified frequency');
style_axes_for_step04_local();

nexttile;
yyaxis left;
plot(T.window_id, T.weighted_voltage_rmse, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('RMSE (V)');
yyaxis right;
plot(T.window_id, T.A_id, 's-', 'LineWidth', 1.2, 'MarkerSize', 4.5);
ylabel('A (mm)');
xlabel('Window');
title('Fit quality and amplitude');
style_axes_for_step04_local();

pngFile = fullfile(figureDir, 'Step04_DirectTemplate_OPRCenterStd_Trend_20250527.png');
pdfFile = fullfile(figureDir, 'Step04_DirectTemplate_OPRCenterStd_Trend_20250527.pdf');
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

summary = table();
summary.ResultFile = string(resultFile);
summary.DominantEO = mode(T.EO_id);
summary.EOConsistency = mean(T.EO_id == summary.DominantEO, 'omitnan');
summary.MeanFrequencyHz = mean(T.fn_id, 'omitnan');
summary.MedianFrequencyHz = median(T.fn_id, 'omitnan');
summary.MeanRMSE = mean(T.weighted_voltage_rmse, 'omitnan');
summary.MedianRMSE = median(T.weighted_voltage_rmse, 'omitnan');
summary.MeanAmplitudeMM = mean(T.A_id, 'omitnan');
summaryCsv = fullfile(resultDir, 'Step04_DirectTemplate_OPRCenterStd_Summary_20250527.csv');
writetable(summary, summaryCsv);

fprintf('Loaded result: %s\n', resultFile);
fprintf('Saved figure: %s\n', pngFile);
fprintf('Saved summary: %s\n', summaryCsv);
disp(summary);
end

function style_axes_for_step04_local()
grid on;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function step01_extract_opr_blade_timing_embedded(rootDir, P)
%% Step01: Extract blade timing and upgrade OPR timing to pulse centers
% This is the fixed preprocessing entry for the direct low-speed-template
% route. It first ensures the legacy Step2 products exist, then rewrites
% jiluOPR/omega so jiluOPR(:,1) is the multi-threshold OPR pulse center.

routeDir = rootDir;
validationRoot = fileparts(routeDir);
mainLegacyDir = fullfile(validationRoot, '20250527', 'legacy');
if exist(mainLegacyDir, 'dir') ~= 7
    error('Legacy directory not found: %s', mainLegacyDir);
end
addpath(mainLegacyDir);

cfg = Get_20250527_BTT_Config();
caseName = cfg.dynamic_cases{1};
caseDir = fullfile(cfg.dataset_root, caseName);
caseOutputDir = fullfile(cfg.output_root, caseName);
showPlots = P.step01.showLegacyPlots;
forceLegacyRebuild = parse_logical_env_local('STEP01_FORCE_LEGACY_REBUILD', false);
forceOprCenterRebuild = parse_logical_env_local('STEP01_FORCE_OPR_CENTER_REBUILD', false);

requiredFiles = {
    fullfile(caseOutputDir, 'jiluOPR.mat')
    fullfile(caseOutputDir, 'omega.mat')
    fullfile(caseOutputDir, 'jilublade_probe1.mat')
    };

hasOutputs = all(cellfun(@isfile, requiredFiles));
if hasOutputs && ~forceLegacyRebuild
    fprintf('Legacy Step01/Step2 outputs already exist. Reusing:\n  %s\n', caseOutputDir);
else
    fprintf('Running legacy Step1/Step2 preprocessing for %s...\n', caseName);
    Step2_Extract_JiluBlade_20250527(caseName, showPlots);
end

oprFile = fullfile(caseOutputDir, 'jiluOPR.mat');
omegaFile = fullfile(caseOutputDir, 'omega.mat');
needsOprCenter = true;
if isfile(oprFile) && ~forceOprCenterRebuild
    S = load(oprFile);
    if isfield(S, 'metadata') && isfield(S.metadata, 'method') && ...
            strcmpi(S.metadata.method, 'multi_threshold_width_center')
        needsOprCenter = false;
    end
end

if needsOprCenter
    fprintf('Rebuilding OPR timing with multi-threshold pulse centers...\n');
    fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
    if isempty(fileIds)
        error('No OPR files found under %s.', caseDir);
    end
    [centerTimes, startTimes, endTimes] = extract_opr_multithreshold_centers_local(caseDir, fileIds, cfg);
    jiluOPR = [centerTimes(:), startTimes(:), endTimes(:)];
    [omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(centerTimes, cfg.blades_num);

    metadata = struct();
    metadata.method = 'multi_threshold_width_center';
    metadata.description = ['OPR center is the median over 30/40/50/60/70% pulse-height ' ...
        'rise/fall midpoint. Legacy threshold start/end are retained in jiluOPR columns 2 and 3.'];
    metadata.threshold = cfg.opr_threshold;
    metadata.level_ratios = [0.30 0.40 0.50 0.60 0.70];
    metadata.pinlv = cfg.pinlv;
    metadata.caseName = caseName;
    metadata.createdBy = mfilename;

    if exist(caseOutputDir, 'dir') ~= 7
        mkdir(caseOutputDir);
    end
    backup_file_if_needed_local(oprFile, 'jiluOPR_legacy_start_edge_backup.mat');
    backup_file_if_needed_local(omegaFile, 'omega_legacy_start_edge_backup.mat');
    save(oprFile, 'jiluOPR', 'metadata');
    save(omegaFile, 'omega_time_s', 'omega_rad_s', 'omega_rpm', 'metadata');
else
    fprintf('OPR timing already uses multi-threshold pulse centers. Reusing:\n  %s\n', oprFile);
end

fprintf('\nStep01 complete. Key outputs:\n');
fprintf('  %s\n', oprFile);
fprintf('  %s\n', omegaFile);
fprintf('  %s\n', fullfile(caseOutputDir, 'jilublade_probe*.mat'));

%% Visualization: RPM and blade timing overview
load(omegaFile, 'omega_time_s', 'omega_rpm');

fig = figure('Name', '20250527 Step01 OPR-center and blade timing', ...
    'Color', 'w', 'Position', [80, 80, 1250, 760], 'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(omega_time_s, omega_rpm, 'k-', 'LineWidth', 1.2);
grid on; box on;
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('OPR-center-derived speed: %s', caseName), 'Interpreter', 'none');

nexttile;
hold on; grid on; box on;
for sid = cfg.sensor_ids
    probeFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d.mat', sid));
    if ~isfile(probeFile)
        continue;
    end
    S = load(probeFile, 'jilublade');
    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 4
        continue;
    end
    tBlade = S.jilublade(:, 3);
    bladeId = S.jilublade(:, 4);
    plot(tBlade, bladeId + 0.08 * (sid - cfg.sensor_ids(1)), '.', ...
        'MarkerSize', 4, 'DisplayName', sprintf('CH%d', sid));
end
xlabel('Time (s)');
ylabel('Blade ID');
title('Extracted blade-passing sequence from legacy Step2');
legend('Location', 'best');

outDir = fullfile(routeDir, 'output', 'figures');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
saveas(fig, fullfile(outDir, 'Step01_OPRCenter_Blade_Timing_20250527.png'));

%% Local functions
function [centerTimes, startTimes, endTimes] = extract_opr_multithreshold_centers_local(caseDir, fileIds, cfg)
tail = [];
gapPoints = 1e4;
centerTimes = [];
startTimes = [];
endTimes = [];
for iFile = 1:numel(fileIds)
    raw = load_raw_case_channel_local(caseDir, cfg.opr_id, fileIds(iFile));
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gapPoints] = segment_signal_like_step2_local(raw, cfg.opr_threshold, gapPoints);
    for iSeg = 1:numel(segments.start_idx)
        centerTimes(end + 1, 1) = compute_opr_multithreshold_center_local( ...
            raw, segments.start_idx(iSeg), segments.end_idx(iSeg), cfg); %#ok<AGROW>
        startTimes(end + 1, 1) = segments.start_sample(iSeg) / cfg.pinlv; %#ok<AGROW>
        endTimes(end + 1, 1) = segments.end_sample(iSeg) / cfg.pinlv; %#ok<AGROW>
    end
end
end

function tCenter = compute_opr_multithreshold_center_local(raw, a, b, cfg)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(size(raw, 1), b + pad);
t = raw(a0:b0, 1) / cfg.pinlv;
v = raw(a0:b0, 2);
try
    vs = smooth(v, 16);
catch
    vs = smoothdata(v, 'movmean', 16);
end
[peakVal, iPeak] = max(vs);
baseVal = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vs);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] .* max(peakVal - baseVal, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(t(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = crossing_time_local(t(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
tCenter = median(centers, 'omitnan');
if ~isfinite(tCenter)
    tCenter = 0.5 * (raw(a, 1) + raw(b, 1)) / cfg.pinlv;
end
end

function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < level & v(2:end) >= level, 1, 'first');
else
    idx = find(v(1:end-1) >= level & v(2:end) < level, 1, 'last');
end
if isempty(idx)
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (level - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end

function [segments, tail, nextGap] = segment_signal_like_step2_local(raw, threshold, gapPoints)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];
nextGap = gapPoints;
if isempty(raw)
    return;
end
sig = raw(:, 2);
try
    sigSmooth = smooth(sig, 16);
catch
    sigSmooth = smoothdata(sig, 'movmean', 16);
end
chase = find(sigSmooth > threshold);
if isempty(chase)
    return;
end
sampleOrder = raw(chase, 1);
segStarts = 1;
segEnds = [];
for ii = 1:(numel(sampleOrder) - 1)
    c = sampleOrder(ii + 1) - sampleOrder(ii);
    if c > nextGap
        segEnds(end + 1, 1) = ii; %#ok<AGROW>
        segStarts(end + 1, 1) = ii + 1; %#ok<AGROW>
        nextGap = 0.6 * c;
    end
end
segEnds(end + 1, 1) = numel(sampleOrder);
tailPoint = chase(segStarts(end)) - floor(nextGap / 2);
if tailPoint > 0 && tailPoint < size(raw, 1)
    tail = raw(tailPoint:end, :);
end
if numel(segStarts) < 2
    return;
end
compStarts = segStarts(1:end-1);
compEnds = segEnds(1:end-1);
segments.start_idx = chase(compStarts);
segments.end_idx = chase(compEnds);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end

function raw = load_raw_case_channel_local(caseDir, channelId, fileId)
file = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
if ~isfile(file)
    raw = [];
    return;
end
loaded = load(file);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
end

function [omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(oprTimes, pulsesPerRev)
if numel(oprTimes) <= pulsesPerRev
    omega_time_s = [];
    omega_rad_s = [];
    omega_rpm = [];
    return;
end
revPeriod = oprTimes((pulsesPerRev + 1):end) - oprTimes(1:(end - pulsesPerRev));
omega_time_s = oprTimes(1:(end - pulsesPerRev));
omega_rad_s = 2 * pi ./ max(revPeriod, eps);
omega_rpm = 60 ./ max(revPeriod, eps);
end

function backup_file_if_needed_local(file, backupName)
if exist(file, 'file') ~= 2
    return;
end
backupFile = fullfile(fileparts(file), backupName);
if exist(backupFile, 'file') ~= 2
    copyfile(file, backupFile);
end
end

function fileIds = list_case_file_ids_local(caseDir, channelId)
d = dir(fullfile(caseDir, sprintf('4-%d-*.mat', channelId)));
fileIds = [];
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channelId) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        fileIds(end + 1) = str2double(tok{1}); %#ok<AGROW>
    end
end
fileIds = sort(fileIds);
end

function tf = parse_logical_env_local(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    tf = defaultValue;
else
    tf = any(strcmp(raw, {'1','true','yes','on'}));
end
end
end

function step01_build_rotating_template_embedded(rootDir)
%% Step01_Build_Rotating_Template_20250527
% Build low-speed rotating waveform templates for the proposed method.
%
% This script borrows the super-Gaussian calibration idea only:
%   1) isolate the main pulse from low-speed rotating waveforms,
%   2) keep a reliable spatial trust region,
%   3) use gradient-like weights when fitting the waveform template.
%
% It does not read SGCalib_*.mat and does not use super-Gaussian model
% parameters. The exported template is a non-parametric spline waveform.

%% Settings
route_dir = rootDir;
validation_root = fileparts(route_dir);
legacy_dir = fullfile(validation_root, '20250527', 'legacy');
if exist(legacy_dir, 'dir') ~= 7
    error('Legacy 20250527 helper folder not found: %s', legacy_dir);
end
addpath(legacy_dir);

legacy_cfg = Get_20250527_BTT_Config();
low_speed_case = legacy_cfg.low_speed_case;
low_speed_data_dir = fullfile(legacy_cfg.dataset_root, low_speed_case);
sensor_config_file = fullfile(legacy_cfg.reference_output_dir, 'Sensor_Config_20250527.mat');

target_blade = 1;
analysis_sensors = [1, 3, 6];       % build a reusable sensor library; Step03 can use subsets
sensor_tag = ['S', sprintf('%d', analysis_sensors)];

calibration_laps = 30;
template_grid_n = 1201;
spline_smoothing = 0.995;
weight_floor = 0.05;
trust_quantile = 0.995;
trust_edge_margin_mm = 0.02;
segment_expand_factor = 0.30;
gradient_energy_quantile = 0.995;
gradient_edge_margin_mm = 0.02;
xrange_mode = lower(strtrim(getenv('STEP01_XRANGE_MODE')));
if isempty(xrange_mode)
    xrange_mode = 'energy';
end
if ~ismember(xrange_mode, {'energy', 'threshold'})
    error('STEP01_XRANGE_MODE must be "energy" or "threshold".');
end
xrange_gradient_min_ratio = parse_nonnegative_numeric_env_local('STEP01_XRANGE_GRADIENT_MIN_RATIO', 0.30);
xrange_amplitude_min_ratio = parse_nonnegative_numeric_env_local('STEP01_XRANGE_AMPLITUDE_MIN_RATIO', 0);
xrange_min_half_width_mm = parse_nonnegative_numeric_env_local('STEP01_XRANGE_MIN_HALF_WIDTH_MM', 2.5);
xrange_max_half_width_mm = parse_nonnegative_numeric_env_local('STEP01_XRANGE_MAX_HALF_WIDTH_MM', 4.2);
stable_window_plan_file = strtrim(getenv('STEP01_STABLE_WINDOW_PLAN'));
template_suffix = sanitize_template_suffix_local(strtrim(getenv('STEP01_TEMPLATE_SUFFIX')));
reference_center_template_file = strtrim(getenv('STEP01_REFERENCE_CENTER_TEMPLATE_FILE'));
center_mode = lower(strtrim(getenv('STEP01_CENTER_MODE')));
if isempty(center_mode)
    center_mode = 'centroid';
end
if ~ismember(center_mode, {'centroid', 'sgfit'})
    error('STEP01_CENTER_MODE must be "centroid" or "sgfit".');
end
if ~isempty(stable_window_plan_file) && isempty(template_suffix)
    template_suffix = '_StableWindow';
end

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
figure_dir = fullfile(output_dir, 'figures');
if exist(template_dir, 'dir') ~= 7; mkdir(template_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s%s_20250527.mat', ...
    target_blade, sensor_tag, template_suffix));

fprintf('\n=== Step01: low-speed rotating spline template ===\n');
fprintf('Blade: %d\n', target_blade);
fprintf('Sensors: %s\n', mat2str(analysis_sensors));
fprintf('Low-speed case: %s\n', low_speed_case);
fprintf('Data source: raw low-speed waveforms + Sensor_Config, not SGCalib.\n');
fprintf('Template center mode: %s\n', center_mode);
if ~isempty(stable_window_plan_file)
    fprintf('Stable-window plan: %s\n', stable_window_plan_file);
end
if ~isempty(reference_center_template_file)
    fprintf('Reference center template: %s\n', reference_center_template_file);
end

if ~isfolder(low_speed_data_dir)
    error('Low-speed data folder not found: %s', low_speed_data_dir);
end
if ~isfile(sensor_config_file)
    error('Sensor_Config not found. Run legacy Step1 once: %s', sensor_config_file);
end

loaded_cfg = load(sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;
stable_window_plan = load_stable_window_plan_local(stable_window_plan_file, target_blade, analysis_sensors);
reference_center = load_reference_center_local(reference_center_template_file, analysis_sensors);

%% Load low-speed raw streams and OPR reference
fprintf('\nLoading low-speed raw streams...\n');
[raw_stream, opr_times, opr_reference] = load_low_speed_streams_local(low_speed_data_dir, legacy_cfg, analysis_sensors);
F_omega_deg = build_phase_speed_local(opr_times, legacy_cfg.blades_num);

%% Build template from low-speed target-blade pulses
Template = struct();
Template.Route = 'low_speed_raw_waveform_sg_like_gated_spline_template_no_gap';
Template.TargetBlade = target_blade;
Template.SensorIDs = analysis_sensors;
Template.SensorTag = sensor_tag;
Template.LowSpeedCase = low_speed_case;
Template.LowSpeedDataDir = low_speed_data_dir;
Template.SensorConfigFile = sensor_config_file;
Template.StableWindowPlanFile = stable_window_plan_file;
Template.OPRReference = opr_reference;
Template.CreatedBy = mfilename;
Template.Settings = struct( ...
    'calibration_laps', calibration_laps, ...
    'template_grid_n', template_grid_n, ...
    'spline_smoothing', spline_smoothing, ...
    'weight_floor', weight_floor, ...
    'trust_quantile', trust_quantile, ...
    'trust_edge_margin_mm', trust_edge_margin_mm, ...
    'segment_expand_factor', segment_expand_factor, ...
    'gradient_energy_quantile', gradient_energy_quantile, ...
    'gradient_edge_margin_mm', gradient_edge_margin_mm, ...
    'xrange_mode', xrange_mode, ...
    'xrange_gradient_min_ratio', xrange_gradient_min_ratio, ...
    'xrange_amplitude_min_ratio', xrange_amplitude_min_ratio, ...
    'xrange_min_half_width_mm', xrange_min_half_width_mm, ...
    'xrange_max_half_width_mm', xrange_max_half_width_mm, ...
    'center_mode', center_mode, ...
    'reference_center_template_file', reference_center_template_file);
Template.Sensor = repmat(struct( ...
    'sensor_id', NaN, ...
    'baseline', NaN, ...
    'threshold', NaN, ...
    'xc', NaN, ...
    'xc_detected', NaN, ...
    'xc_reference_source', '', ...
    'x_grid', [], ...
    'v_grid', [], ...
    'dv_dx', [], ...
    'bin_weight', [], ...
    'x_domain', [], ...
    'selection_mode', '', ...
    'point_count', NaN, ...
    'wide_point_count', NaN, ...
    'lap_count', NaN, ...
    'selected_pulse_indices', [], ...
    'fit_method', '', ...
    'spline_smoothing', NaN), numel(analysis_sensors), 1);

plot_cache = repmat(struct('sensor_id', NaN, 'x_wide', [], 'v_wide', [], ...
    'x_selected', [], 'v_selected', [], 'weight_selected', [], ...
    'xc_detected_rel', NaN, 'xc_reference_source', ''), numel(analysis_sensors), 1);

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    threshold = legacy_cfg.sensor_thresholds(sid);
    theta_std = Sensor_Config.Standard_Relative_Angles(sid, target_blade);
    theta_std = convert_standard_angle_to_opr_center_local(theta_std, opr_reference);
    if ~isfinite(theta_std)
        error('Missing standard relative angle for CH%d blade %d.', sid, target_blade);
    end

    R = raw_stream(sid);
    segments = extract_sensor_segments_local(R.T, R.V, threshold, legacy_cfg.gap_points);
    if isempty(segments.start_idx)
        error('No low-speed pulses detected for CH%d.', sid);
    end
    baseline = estimate_background_baseline_local(R.V, segments);

    if isKey(Sensor_Config.Target_Indices, sid)
        start_idx = Sensor_Config.Target_Indices(sid);
    else
        error('Sensor_Config.Target_Indices does not contain CH%d.', sid);
    end
    sensor_laps = calibration_laps;
    start_offset_lap = 0;
    if ~isempty(stable_window_plan)
        plan_row = stable_window_plan(stable_window_plan.SensorID == sid, :);
        if height(plan_row) ~= 1
            error('Stable-window plan must contain exactly one row for CH%d blade %d.', sid, target_blade);
        end
        sensor_laps = round(plan_row.LapCount(1));
        start_offset_lap = round(plan_row.StartOffsetLap(1));
    end
    selected_idx = start_idx + start_offset_lap * legacy_cfg.blades_num + ...
        (0:(sensor_laps - 1)) * legacy_cfg.blades_num + (target_blade - 1);
    selected_idx = selected_idx(selected_idx >= 1 & selected_idx <= numel(segments.start_idx));
    if numel(selected_idx) < 3
        error('CH%d has only %d selected low-speed target-blade pulses.', sid, numel(selected_idx));
    end

    [x_wide, v_wide] = restore_low_speed_point_cloud_local( ...
        R, segments, selected_idx, opr_times, F_omega_deg, theta_std, ...
        legacy_cfg.r_tip_mm, segment_expand_factor);

    select_cfg = struct('threshold', threshold, 'baseline', baseline, ...
        'weight_floor', weight_floor, 'trust_quantile', trust_quantile, ...
        'trust_edge_margin_mm', trust_edge_margin_mm, ...
        'gradient_energy_quantile', gradient_energy_quantile, ...
        'gradient_edge_margin_mm', gradient_edge_margin_mm, ...
        'xrange_mode', xrange_mode, ...
        'xrange_gradient_min_ratio', xrange_gradient_min_ratio, ...
        'xrange_amplitude_min_ratio', xrange_amplitude_min_ratio, ...
        'xrange_min_half_width_mm', xrange_min_half_width_mm, ...
        'xrange_max_half_width_mm', xrange_max_half_width_mm, ...
        'center_mode', center_mode);
    P = select_sg_like_template_points_local(x_wide, v_wide, select_cfg);

    xc_for_coordinate = P.xc;
    xc_reference_source = 'detected_from_current_template_points';
    iref = find(reference_center.sensor_id == sid, 1);
    if ~isempty(iref) && isfinite(reference_center.xc_mm(iref))
        xc_for_coordinate = reference_center.xc_mm(iref);
        xc_reference_source = 'reference_center_template';
    end

    x_rel = P.x_selected(:) - xc_for_coordinate;
    v_raw = P.v_selected(:);
    point_weight = P.weight_selected(:);
    finite_mask = isfinite(x_rel) & isfinite(v_raw) & isfinite(point_weight);
    x_rel = x_rel(finite_mask);
    v_raw = v_raw(finite_mask);
    point_weight = point_weight(finite_mask);

    x_left = min(x_rel);
    x_right = max(x_rel);
    x_grid = linspace(x_left, x_right, template_grid_n).';
    edges = linspace(x_left, x_right, template_grid_n + 1).';

    bin_id = discretize(x_rel, edges);
    valid_bin = ~isnan(bin_id);
    v_med = accumarray(bin_id(valid_bin), v_raw(valid_bin), [template_grid_n, 1], @median, NaN);
    w_bin = accumarray(bin_id(valid_bin), point_weight(valid_bin), [template_grid_n, 1], @mean, NaN);
    bin_count = accumarray(bin_id(valid_bin), 1, [template_grid_n, 1], @sum, 0);
    valid_grid = isfinite(v_med) & bin_count >= 3;
    if nnz(valid_grid) < 30
        error('Too few valid template bins for CH%d.', sid);
    end

    x_fit = x_grid(valid_grid);
    v_fit = v_med(valid_grid);
    w_fit = w_bin(valid_grid);
    [x_unique, ~, ic] = unique(round(x_fit, 6));
    v_unique = accumarray(ic, v_fit, [], @median);
    w_unique = accumarray(ic, w_fit, [], @mean);

    if exist('csaps', 'file') == 2
        try
            v_grid = csaps(x_unique, v_unique, spline_smoothing, x_grid, w_unique);
            fit_method = 'weighted csaps smoothing spline';
        catch
            v_grid = csaps(x_unique, v_unique, spline_smoothing, x_grid);
            fit_method = 'csaps smoothing spline';
        end
    else
        v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
        v_grid = smoothdata(v_fill, 'sgolay', 41);
        fit_method = 'median bins + sgolay fallback';
    end
    v_grid = v_grid(:);
    dv_dx = gradient(v_grid, x_grid);

    Template.Sensor(is).sensor_id = sid;
    Template.Sensor(is).baseline = baseline;
    Template.Sensor(is).threshold = threshold;
    Template.Sensor(is).xc = xc_for_coordinate;
    Template.Sensor(is).xc_detected = P.xc;
    Template.Sensor(is).xc_reference_source = xc_reference_source;
    Template.Sensor(is).x_grid = x_grid;
    Template.Sensor(is).v_grid = v_grid;
    Template.Sensor(is).dv_dx = dv_dx;
    Template.Sensor(is).bin_weight = fillmissing(w_bin, 'linear', 'EndValues', 'nearest');
    Template.Sensor(is).x_domain = [x_left, x_right];
    Template.Sensor(is).selection_mode = P.mode;
    Template.Sensor(is).point_count = numel(x_rel);
    Template.Sensor(is).wide_point_count = numel(x_wide);
    Template.Sensor(is).lap_count = numel(selected_idx);
    Template.Sensor(is).selected_pulse_indices = selected_idx(:);
    Template.Sensor(is).fit_method = fit_method;
    Template.Sensor(is).spline_smoothing = spline_smoothing;

    plot_cache(is).sensor_id = sid;
    plot_cache(is).x_wide = x_wide(:) - xc_for_coordinate;
    plot_cache(is).v_wide = v_wide(:);
    plot_cache(is).x_selected = P.x_selected(:) - xc_for_coordinate;
    plot_cache(is).v_selected = P.v_selected(:);
    plot_cache(is).weight_selected = P.weight_selected(:);
    plot_cache(is).xc_detected_rel = P.xc - xc_for_coordinate;
    plot_cache(is).xc_reference_source = xc_reference_source;

    fprintf(['CH%d: baseline %.5f V, threshold %.3f V, pulses %d, ' ...
        'selected points %d/%d, domain [%.3f, %.3f] mm, xc %.4f mm, ' ...
        'detected xc %.4f mm, source %s\n'], ...
        sid, baseline, threshold, numel(selected_idx), numel(x_rel), numel(x_wide), ...
        x_left, x_right, xc_for_coordinate, P.xc, xc_reference_source);
end

save(template_file, 'Template', '-v7.3');
fprintf('Saved template: %s\n', template_file);

%% Visualization
fig = figure('Name', 'Step01 low-speed gated spline templates', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
tiledlayout(fig, numel(analysis_sensors), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor(is);
    P = plot_cache(is);

    [x_wide_plot, v_wide_plot] = downsample_cloud_local(P.x_wide, P.v_wide, 25000);
    [x_sel_plot, v_sel_plot] = downsample_cloud_local(P.x_selected, P.v_selected, 25000);

    nexttile;
    plot(x_wide_plot, v_wide_plot, '.', 'Color', [0.78 0.78 0.78], 'MarkerSize', 3, ...
        'DisplayName', 'wide pulse cloud'); hold on;
    plot(x_sel_plot, v_sel_plot, 'k.', 'MarkerSize', 3, 'DisplayName', 'selected trust cloud');
    plot(Tpl.x_grid, Tpl.v_grid, 'r-', 'LineWidth', 1.8, 'DisplayName', 'spline template');
    xline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'coordinate center');
    if isfinite(P.xc_detected_rel) && abs(P.xc_detected_rel) > 1e-6
        xline(P.xc_detected_rel, 'm:', 'LineWidth', 1.0, 'DisplayName', 'detected center');
    end
    xline(Tpl.x_domain(1), 'b-.', 'LineWidth', 1.1, 'DisplayName', 'trust boundary');
    xline(Tpl.x_domain(2), 'b-.', 'LineWidth', 1.1, 'HandleVisibility', 'off');
    yline(Tpl.baseline, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'DisplayName', 'baseline');
    xlabel('x relative to template center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d template, %d pulses, %s', sid, Tpl.lap_count, P.xc_reference_source), ...
        'Interpreter', 'none');
    box on; grid on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');
    legend('Location', 'best');

    nexttile;
    yyaxis left;
    plot(Tpl.x_grid, Tpl.dv_dx, 'b-', 'LineWidth', 1.5);
    ylabel('dV/dx (V/mm)');
    yyaxis right;
    plot(Tpl.x_grid, Tpl.bin_weight, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
    ylabel('Template weight');
    xlabel('x relative to template center (mm)');
    title(sprintf('CH%d slope and SG-like weight', sid));
    box on; grid on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');
end

exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_LowSpeedSplineTemplates_B%d_%s%s.png', ...
    target_blade, sensor_tag, template_suffix)), 'Resolution', 300);
exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_LowSpeedSplineTemplates_B%d_%s%s.pdf', ...
    target_blade, sensor_tag, template_suffix)), 'ContentType', 'vector');

%% Local helpers
function suffix = sanitize_template_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function v = parse_nonnegative_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
v = str2double(raw);
if ~isfinite(v) || v < 0
    error('%s must be a nonnegative numeric value.', name);
end
end

function reference_center = load_reference_center_local(template_file, sensor_ids)
reference_center = table(sensor_ids(:), nan(numel(sensor_ids), 1), ...
    'VariableNames', {'sensor_id', 'xc_mm'});
if isempty(template_file)
    return;
end
if exist(template_file, 'file') ~= 2
    error('Reference center template file not found: %s', template_file);
end
loaded = load(template_file, 'Template');
if ~isfield(loaded, 'Template') || ~isfield(loaded.Template, 'Sensor')
    error('Reference center template must contain Template.Sensor: %s', template_file);
end
available = [loaded.Template.Sensor.sensor_id];
for i = 1:numel(sensor_ids)
    idx = find(available == sensor_ids(i), 1);
    if isempty(idx) || ~isfield(loaded.Template.Sensor(idx), 'xc')
        error('Reference center template missing xc for CH%d: %s', sensor_ids(i), template_file);
    end
    reference_center.xc_mm(i) = loaded.Template.Sensor(idx).xc;
end
end

function plan = load_stable_window_plan_local(plan_file, target_blade, analysis_sensors)
plan = table();
if isempty(plan_file)
    return;
end
if exist(plan_file, 'file') ~= 2
    error('Stable-window plan file not found: %s', plan_file);
end
[~, ~, ext] = fileparts(plan_file);
if strcmpi(ext, '.mat')
    loaded = load(plan_file);
    if isfield(loaded, 'window_plan')
        plan = loaded.window_plan;
    else
        error('MAT stable-window plan must contain variable window_plan: %s', plan_file);
    end
else
    plan = readtable(plan_file);
end
required_vars = {'BladeID', 'SensorID', 'LapCount', 'StartOffsetLap'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, plan.Properties.VariableNames)
        error('Stable-window plan missing column %s: %s', required_vars{i}, plan_file);
    end
end
plan = plan(plan.BladeID == target_blade & ismember(plan.SensorID, analysis_sensors), :);
if height(plan) ~= numel(analysis_sensors)
    error('Stable-window plan must contain one row per requested sensor for blade %d.', target_blade);
end
end

function [raw_stream, opr_times, opr_reference] = load_low_speed_streams_local(case_dir, cfg, sensor_ids)
file_ids = list_case_file_ids_local(case_dir, cfg.opr_id);
if isempty(file_ids)
    error('No OPR files found in %s.', case_dir);
end

max_sid = max([sensor_ids(:); cfg.opr_id]);
raw_stream(max_sid) = struct('T', [], 'V', []);
opr_T = [];
opr_V = [];
last_end = [];

for i = 1:numel(file_ids)
    file_id = file_ids(i);
    [t_opr, v_opr] = load_raw_channel_local(case_dir, cfg.opr_id, file_id, cfg.pinlv);
    if isempty(t_opr)
        continue;
    end
    if isempty(last_end)
        offset = 0;
    else
        offset = last_end + 1 / cfg.pinlv - t_opr(1);
    end
    t_opr = t_opr(:) + offset;
    last_end = t_opr(end);

    if i == 1 && cfg.initial_trim_points > 0 && numel(t_opr) > cfg.initial_trim_points
        keep = (cfg.initial_trim_points + 1):numel(t_opr);
        t_opr = t_opr(keep);
        v_opr = v_opr(keep);
    end
    opr_T = [opr_T; t_opr(:)]; %#ok<AGROW>
    opr_V = [opr_V; v_opr(:)]; %#ok<AGROW>

    for sensorId = sensor_ids
        if ~isfile(fullfile(case_dir, sprintf('4-%d-%d.mat', sensorId, file_id)))
            continue;
        end
        [t_local, v_local] = load_raw_channel_local(case_dir, sensorId, file_id, cfg.pinlv);
        t_local = t_local(:) + offset;
        if i == 1 && cfg.initial_trim_points > 0 && numel(t_local) > cfg.initial_trim_points
            keep = (cfg.initial_trim_points + 1):numel(t_local);
            t_local = t_local(keep);
            v_local = v_local(keep);
        end
        raw_stream(sensorId).T = [raw_stream(sensorId).T; t_local(:)]; %#ok<AGROW>
        raw_stream(sensorId).V = [raw_stream(sensorId).V; v_local(:)]; %#ok<AGROW>
    end
end

opr_segments = extract_opr_segments_local(opr_T, opr_V, cfg.opr_threshold, cfg.gap_points);
opr_times = opr_segments.arrival_time(:);
opr_reference = build_opr_reference_local(opr_segments, cfg.blades_num, cfg.r_tip_mm);
fprintf(['Loaded %d OPR pulses from low-speed case. OPR center shift: ' ...
    '%.4f deg, %.4f mm from rising edge.\n'], numel(opr_times), ...
    opr_reference.phase_shift_deg, opr_reference.phase_shift_mm);
end

function theta_std_center = convert_standard_angle_to_opr_center_local(theta_std_start, opr_reference)
theta_std_center = theta_std_start;
if isstruct(opr_reference) && isfield(opr_reference, 'phase_shift_deg') && ...
        isfinite(opr_reference.phase_shift_deg)
    theta_std_center = theta_std_start - opr_reference.phase_shift_deg;
end
end

function opr_reference = build_opr_reference_local(opr_segments, pulses_per_rev, r_tip_mm)
opr_reference = struct('mode', 'multi_threshold_center', ...
    'standard_angle_reference', 'opr_pulse_center', ...
    'phase_shift_deg', 0, ...
    'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if isempty(opr_segments.arrival_time) || isempty(opr_segments.start_time) || pulses_per_rev < 1
    return;
end
center_time = opr_segments.arrival_time(:);
start_time = opr_segments.start_time(:);
n = min(numel(center_time) - pulses_per_rev, numel(start_time));
if n < 1
    return;
end
dt_center = center_time(1:n) - start_time(1:n);
dt_rev = center_time((1:n) + pulses_per_rev) - center_time(1:n);
valid = isfinite(dt_center) & isfinite(dt_rev) & dt_rev > eps;
if ~any(valid)
    return;
end
shift_deg = 360 * dt_center(valid) ./ dt_rev(valid);
opr_reference.phase_shift_deg = median(shift_deg, 'omitnan');
opr_reference.phase_shift_mm = opr_reference.phase_shift_deg * (pi / 180) * r_tip_mm;
opr_reference.median_center_minus_start_s = median(dt_center(valid), 'omitnan');
end

function file_ids = list_case_file_ids_local(case_dir, channel_id)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', channel_id)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channel_id) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        file_ids(i) = str2double(tok{1});
    end
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));
end

function [t_sec, v] = load_raw_channel_local(case_dir, channel_id, file_id, fs)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', channel_id, file_id));
if ~isfile(filepath)
    t_sec = [];
    v = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / fs;
v = raw(:, 2);
end

function segments = extract_opr_segments_local(t, v, threshold, gap_points)
idx = find(v(:) > threshold);
segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, false);
end

function segments = extract_sensor_segments_local(t, v, threshold, gap_points)
if numel(v) >= 21
    v_smooth = sgolayfilt(v(:), 3, 21);
elseif numel(v) >= 5
    v_smooth = smoothdata(v(:), 'movmean', max(3, 2 * floor(numel(v) / 4) + 1));
else
    v_smooth = v(:);
end
idx = find(v_smooth > threshold);
segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, true);
end

function segments = split_segments_to_struct_local(t, v, idx, gap_points, threshold, use_centroid)
segments = struct('start_idx', [], 'end_idx', [], 'start_time', [], ...
    'end_time', [], 'arrival_time', [], 'peak_value', []);
if isempty(idx)
    return;
end
jumps = find(diff(idx) > gap_points);
seg_start_pos = [1; jumps(:) + 1];
seg_end_pos = [jumps(:); numel(idx)];
n = numel(seg_start_pos);
segments.start_idx = zeros(n, 1);
segments.end_idx = zeros(n, 1);
segments.start_time = zeros(n, 1);
segments.end_time = zeros(n, 1);
segments.arrival_time = zeros(n, 1);
segments.peak_value = zeros(n, 1);

for k = 1:n
    a = idx(seg_start_pos(k));
    b = idx(seg_end_pos(k));
    segments.start_idx(k) = a;
    segments.end_idx(k) = b;
    segments.start_time(k) = t(a);
    segments.end_time(k) = t(b);
    [pk, ipk] = max(v(a:b));
    segments.peak_value(k) = pk;
    if use_centroid
        tt = t(a:b);
        vv = max(v(a:b) - threshold, 0);
        if sum(vv) > eps
            segments.arrival_time(k) = sum(tt(:) .* vv(:)) / sum(vv);
        else
            segments.arrival_time(k) = t(a + ipk - 1);
        end
    else
        segments.arrival_time(k) = compute_opr_multithreshold_center_local(t, v, a, b);
    end
end
end

function t_center = compute_opr_multithreshold_center_local(t, v, a, b)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(numel(v), b + pad);
tt = t(a0:b0);
vv = v(a0:b0);
try
    vs = smooth(vv, 16);
catch
    vs = smoothdata(vv, 'movmean', 16);
end
[peakVal, iPeak] = max(vs);
baseVal = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vs);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] .* max(peakVal - baseVal, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(tt(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = crossing_time_local(tt(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
t_center = median(centers, 'omitnan');
if ~isfinite(t_center)
    t_center = 0.5 * (t(a) + t(b));
end
end

function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < level & v(2:end) >= level, 1, 'first');
else
    idx = find(v(1:end-1) >= level & v(2:end) < level, 1, 'last');
end
if isempty(idx)
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (level - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end

function baseline = estimate_background_baseline_local(v, segments)
mask = true(size(v));
for k = 1:numel(segments.start_idx)
    a = segments.start_idx(k);
    b = segments.end_idx(k);
    pad = max(10, round(1.5 * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(v), b + pad);
    mask(aa:bb) = false;
end
if nnz(mask) > 100
    baseline = median(v(mask), 'omitnan');
else
    baseline = median(v(v <= prctile(v, 30)), 'omitnan');
end
end

function F_omega_deg = build_phase_speed_local(opr_times, blades_num)
spd_t = opr_times(1:(end - blades_num));
spd_v = 360 ./ max(opr_times((blades_num + 1):end) - opr_times(1:(end - blades_num)), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end

function theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg)
theta_points = zeros(size(t_seg));
for i = 1:numel(t_seg)
    t_grid = linspace(t_ref, t_seg(i), 10);
    theta_points(i) = trapz(t_grid, F_omega_deg(t_grid));
end
end

function [x_wide, v_wide] = restore_low_speed_point_cloud_local( ...
    R, segments, selected_idx, opr_times, F_omega_deg, theta_std, r_tip_mm, expand_factor)
x_wide = [];
v_wide = [];
for ii = 1:numel(selected_idx)
    pidx = selected_idx(ii);
    a = segments.start_idx(pidx);
    b = segments.end_idx(pidx);
    pad = max(8, round(expand_factor * (b - a + 1)));
    aa = max(1, a - pad);
    bb = min(numel(R.T), b + pad);
    t_seg = R.T(aa:bb);
    v_seg = R.V(aa:bb);
    t_arrival = segments.arrival_time(pidx);
    idx_prev = find(opr_times < t_arrival, 1, 'last');
    if isempty(idx_prev)
        continue;
    end
    theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
    theta_diff_deg = mod(theta_points_deg - theta_std + 180, 360) - 180;
    x_points_mm = theta_diff_deg * (pi / 180) * r_tip_mm;
    x_wide = [x_wide; x_points_mm(:)]; %#ok<AGROW>
    v_wide = [v_wide; v_seg(:)]; %#ok<AGROW>
end
end

function P = select_sg_like_template_points_local(x_wide, v_wide, cfg)
finite = isfinite(x_wide) & isfinite(v_wide);
x_wide = x_wide(finite);
v_wide = v_wide(finite);
v_zero = max(v_wide - cfg.baseline, 0);
if numel(x_wide) < 50
    error('Too few low-speed cloud points for template selection.');
end

high_level = prctile(v_zero, 85);
high_mask = v_zero >= high_level & v_zero > 0;
if nnz(high_mask) >= 5 && sum(v_zero(high_mask)) > eps
    xc = sum(x_wide(high_mask) .* v_zero(high_mask)) / sum(v_zero(high_mask));
else
    [~, imax] = max(v_zero);
    xc = x_wide(imax);
end

strict_mask = v_wide >= cfg.threshold;
dx_strict = x_wide(strict_mask) - xc;
if nnz(strict_mask) < 30
    strict_mask = v_zero >= prctile(v_zero, 60);
    dx_strict = x_wide(strict_mask) - xc;
end

[left_L, right_L, xrange_ok, xrange_mode_used] = estimate_xrange_span_local( ...
    x_wide, v_wide, cfg.baseline, xc, cfg);
if ~xrange_ok
    left_dist = -dx_strict(dx_strict < 0);
    right_dist = dx_strict(dx_strict > 0);
    if isempty(left_dist) || isempty(right_dist)
        L = prctile(abs(dx_strict), 95);
        left_L = L;
        right_L = L;
    else
        left_L = prctile(left_dist, 100 * cfg.trust_quantile);
        right_L = prctile(right_dist, 100 * cfg.trust_quantile);
    end
    left_L = max(left_L - cfg.trust_edge_margin_mm, 0.1);
    right_L = max(right_L - cfg.trust_edge_margin_mm, 0.1);
end

if xrange_ok
    trust_mask = x_wide >= xc - left_L & x_wide <= xc + right_L;
else
    trust_mask = strict_mask & x_wide >= xc - left_L & x_wide <= xc + right_L;
end
if nnz(trust_mask) < 30
    trust_mask = strict_mask;
end

x_selected = x_wide(trust_mask);
v_selected = v_wide(trust_mask);
if isfield(cfg, 'center_mode') && strcmpi(cfg.center_mode, 'sgfit')
    xc = refine_center_by_sg_fit_local(x_selected, v_selected, cfg.baseline, xc);
end
weight_selected = build_gradient_weight_from_wide_cloud_local( ...
    x_selected, x_wide, v_wide, cfg.baseline, cfg.weight_floor);

P = struct();
if xrange_ok
    P.mode = ['raw_low_speed_', xrange_mode_used, '_gate_gradient_weight'];
else
    P.mode = 'raw_low_speed_sg_like_strict_gate_trust_window_gradient_weight';
end
P.xc = xc;
P.x_wide = x_wide;
P.v_wide = v_wide;
P.x_selected = x_selected;
P.v_selected = v_selected;
P.weight_selected = weight_selected;
end

function [left_L, right_L, ok, mode_used] = estimate_xrange_span_local(x, v, baseline, xc, cfg)
if strcmpi(cfg.xrange_mode, 'threshold')
    [left_L, right_L, ok] = estimate_gradient_threshold_span_local( ...
        x, v, baseline, xc, cfg.xrange_gradient_min_ratio, cfg.xrange_amplitude_min_ratio, ...
        cfg.xrange_min_half_width_mm, cfg.xrange_max_half_width_mm);
    mode_used = 'gradient_threshold';
else
    [left_L, right_L, ok] = estimate_gradient_energy_span_local( ...
        x, v, baseline, xc, cfg.gradient_energy_quantile, cfg.gradient_edge_margin_mm);
    mode_used = 'gradient_energy';
end
end

function [left_L, right_L, ok] = estimate_gradient_threshold_span_local( ...
    x, v, baseline, xc, gradient_min_ratio, amplitude_min_ratio, min_half_width, max_half_width)
left_L = NaN;
right_L = NaN;
ok = false;
[x_grid, v_smooth] = build_denoised_static_profile_local(x, v, baseline);
if numel(x_grid) < 20
    return;
end
amp_norm = normalize01_local(v_smooth(:));
g_abs = abs(gradient(v_smooth(:), x_grid(:)));
g_threshold = estimate_noise_aware_gradient_threshold_local(g_abs, amp_norm, gradient_min_ratio);
effective = g_abs >= g_threshold;
if amplitude_min_ratio > 0
    effective = effective & amp_norm >= amplitude_min_ratio;
end
[~, ic] = min(abs(x_grid - xc));
if isempty(ic) || ~isfinite(ic)
    return;
end

left_candidates = find(x_grid(:) < xc & effective(:));
right_candidates = find(x_grid(:) > xc & effective(:));
if isempty(left_candidates) || isempty(right_candidates)
    return;
end

left_idx = left_candidates(1);
right_idx = right_candidates(end);
left_L = max(xc - x_grid(left_idx), 0);
right_L = max(x_grid(right_idx) - xc, 0);
if left_L < min_half_width
    left_L = min_half_width;
end
if right_L < min_half_width
    right_L = min_half_width;
end
if isfinite(max_half_width) && max_half_width > 0
    left_L = min(left_L, max_half_width);
    right_L = min(right_L, max_half_width);
end
ok = isfinite(left_L) && isfinite(right_L) && left_L > 0 && right_L > 0;
end

function [left_L, right_L, ok] = estimate_gradient_energy_span_local( ...
    x, v, baseline, xc, energy_quantile, edge_margin_mm)
left_L = NaN;
right_L = NaN;
ok = false;
[x_grid, v_smooth] = build_denoised_static_profile_local(x, v, baseline);
if numel(x_grid) < 20 || max(v_smooth) <= 0
    return;
end
g = abs(gradient(v_smooth, x_grid));
signal_gate = v_smooth >= prctile(v_smooth, 40);
energy = g(:) .* signal_gate(:) + 0.05 * max(g(:)) * normalize01_local(v_smooth(:));
if sum(energy) <= eps
    return;
end

tail = (1 - energy_quantile) / 2;
cum = cumsum(energy) ./ sum(energy);
[cum_unique, ia] = unique(cum, 'stable');
x_unique = x_grid(ia);
if numel(cum_unique) < 2
    return;
end
x_lo = interp1(cum_unique, x_unique, tail, 'linear', 'extrap');
x_hi = interp1(cum_unique, x_unique, 1 - tail, 'linear', 'extrap');
if ~isfinite(x_lo) || ~isfinite(x_hi) || x_hi <= x_lo || xc <= x_lo || xc >= x_hi
    return;
end
left_L = max(xc - x_lo - edge_margin_mm, 0.1);
right_L = max(x_hi - xc - edge_margin_mm, 0.1);
ok = true;
end

function [x_grid, v_smooth] = build_denoised_static_profile_local(x, v, baseline)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v_zero = max(v(finite) - baseline, 0);
x_grid = [];
v_smooth = [];
if numel(x) < 50 || max(v_zero) <= 0
    return;
end
[x_sort, idx] = sort(x);
v_sort = v_zero(idx);
grid_n = min(700, max(120, round(numel(x_sort) / 120)));
edges = linspace(min(x_sort), max(x_sort), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v_sort(valid), [grid_n, 1], @median, NaN);
if nnz(isfinite(v_med)) < 15
    x_grid = [];
    v_smooth = [];
    return;
end
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
end

function g_threshold = estimate_noise_aware_gradient_threshold_local(g_abs, amp_norm, ratio)
g_abs = g_abs(:);
amp_norm = amp_norm(:);
finite = isfinite(g_abs) & isfinite(amp_norm);
g_abs = g_abs(finite);
amp_norm = amp_norm(finite);
if isempty(g_abs)
    g_threshold = Inf;
    return;
end
noise_gate = amp_norm <= 0.10;
signal_gate = amp_norm >= 0.20;
if nnz(noise_gate) >= 10
    g_noise = prctile(g_abs(noise_gate), 95);
else
    g_noise = prctile(g_abs, 10);
end
if nnz(signal_gate) >= 10
    g_signal = prctile(g_abs(signal_gate), 95);
else
    g_signal = prctile(g_abs, 95);
end
if ~isfinite(g_noise)
    g_noise = 0;
end
if ~isfinite(g_signal) || g_signal <= g_noise
    g_signal = max(g_abs);
end
g_threshold = g_noise + ratio * max(g_signal - g_noise, 0);
end

function y = normalize01_local(x)
x = x(:);
xmin = min(x, [], 'omitnan');
xmax = max(x, [], 'omitnan');
if ~isfinite(xmin) || ~isfinite(xmax) || xmax <= xmin
    y = zeros(size(x));
else
    y = (x - xmin) ./ (xmax - xmin);
end
end

function xc = refine_center_by_sg_fit_local(x, v, baseline, xc_initial)
finite = isfinite(x) & isfinite(v);
x = x(finite);
v = v(finite);
if numel(x) < 50
    xc = xc_initial;
    return;
end
[peak_val, idx_peak] = max(v);
B0 = max(peak_val - baseline, 0.1);
w0 = max(std(x), 0.2);
n0 = 3.0;
xc0 = xc_initial;
if ~isfinite(xc0)
    xc0 = x(idx_peak);
end
sg_model = @(p, xx) p(1) .* exp(-abs((xx - p(4)) ./ max(p(2), 1e-6)).^max(p(3), 1e-6)) + baseline;
obj_fun = @(p) sum((v - sg_model(p, x)).^2);
p0 = [B0, w0, n0, xc0];
lb = [0.05, 0.05, 1.2, min(x) - 0.5];
ub = [max(10, 2 * B0 + 0.5), 6.0, 10.0, max(x) + 0.5];
try
    if exist('fmincon', 'file') == 2
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        p_opt = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], opts);
    else
        penalty_obj = @(p) obj_fun(min(max(p, lb), ub));
        opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
        p_opt = fminsearch(penalty_obj, p0, opts);
        p_opt = min(max(p_opt, lb), ub);
    end
    if isfinite(p_opt(4))
        xc = p_opt(4);
    else
        xc = xc_initial;
    end
catch
    xc = xc_initial;
end
end

function W = build_gradient_weight_from_wide_cloud_local(x_query, x_wide, v_wide, baseline, weight_floor)
finite = isfinite(x_wide) & isfinite(v_wide);
x_wide = x_wide(finite);
v_zero = max(v_wide(finite) - baseline, 0);
if numel(x_wide) < 20 || all(v_zero == 0)
    W = ones(size(x_query));
    return;
end

[x_sort, idx] = sort(x_wide);
v_sort = v_zero(idx);
grid_n = min(500, max(80, round(numel(x_sort) / 200)));
edges = linspace(min(x_sort), max(x_sort), grid_n + 1).';
x_grid = 0.5 * (edges(1:end-1) + edges(2:end));
bin_id = discretize(x_sort, edges);
valid = ~isnan(bin_id);
v_med = accumarray(bin_id(valid), v_sort(valid), [grid_n, 1], @median, NaN);
valid_grid = isfinite(v_med);
if nnz(valid_grid) < 10
    W = ones(size(x_query));
    return;
end
v_fill = fillmissing(v_med, 'linear', 'EndValues', 'nearest');
span = min(41, 2 * floor((numel(v_fill) - 1) / 2) + 1);
if span >= 5
    v_smooth = smoothdata(v_fill, 'sgolay', span);
else
    v_smooth = v_fill;
end
dVdx = abs(gradient(v_smooth, x_grid));
w_grid = weight_floor + (1 - weight_floor) * dVdx ./ max(dVdx + eps);
W = interp1(x_grid, w_grid, x_query(:), 'linear', weight_floor);
W(~isfinite(W)) = weight_floor;
W = max(W, weight_floor);
end

function [x_ds, v_ds] = downsample_cloud_local(x, v, max_points)
keep = isfinite(x) & isfinite(v);
x = x(keep);
v = v(keep);
if numel(x) > max_points
    idx = round(linspace(1, numel(x), max_points));
else
    idx = 1:numel(x);
end
x_ds = x(idx);
v_ds = v(idx);
end
end

function step02_build_dynamic_map_embedded(rootDir)
%% Step02_Build_Dynamic_Map_20250527
% Build sliding-window dynamic waveform maps for the proposed method.
%
% The vibration-region settings are kept the same as the old SG Step5 case:
%   cfg.target_blades        = 1
%   cfg.analysis_sensors     = [1 3 6]
%   cfg.analysis_start_time  = 1.5 s
%   cfg.target_laps          = 20
%   cfg.analysis_win_size    = 3 laps
%   cfg.sliding_step         = 1 lap
%
% This script only prepares waveform windows. The identification model is
% changed in Step03.

%% Settings
route_dir = rootDir;
old_result_file = ...
    'E:\0灏忚鏂?绋嬪簭\0鍗氬＋鏈熼棿灏忚鏂?绋嬪簭\7瓒呴珮鏂ā鍨?鏉冮噸-鐬€乗绋嬪簭\鐩村彾鐗囬獙璇乗瀹為獙楠岃瘉\瓒呴珮鏂痋20250527閫傞厤\Result_single_sync_20250527_B1_S136_20250526_2500_3500_t400_Start1p5s.mat';

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1, 3, 6];
cfg.analysis_start_time = 1.5;
cfg.target_laps = parse_positive_integer_env_local('STEP02_TARGET_LAPS', 20);
cfg.analysis_win_size = parse_positive_integer_env_local('STEP02_ANALYSIS_WIN_SIZE', 3);
cfg.sliding_step = parse_positive_integer_env_local('STEP02_SLIDING_STEP', 1);
sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
source_result_sensor_tag = 'S136';
dynamic_suffix = sanitize_dynamic_suffix_local(strtrim(getenv('STEP02_DYNAMIC_SUFFIX')));

output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
figure_dir = fullfile(output_dir, 'figures');
if exist(dynamic_dir, 'dir') ~= 7; mkdir(dynamic_dir); end
if exist(figure_dir, 'dir') ~= 7; mkdir(figure_dir); end

template_suffix = strtrim(getenv('STEP02_TEMPLATE_SUFFIX'));
if isempty(template_suffix)
    template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', ...
        cfg.target_blades, sensor_tag));
else
    template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_%s_%s_20250527.mat', ...
        cfg.target_blades, sensor_tag, template_suffix));
end
dynamic_file = fullfile(dynamic_dir, sprintf('DynamicMap_B%d_%s_SlidingWindows%s_20250527.mat', ...
    cfg.target_blades, sensor_tag, dynamic_suffix));

fprintf('\n=== Step02: sliding-window dynamic waveform maps ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Sensor tag: %s\n', sensor_tag);
fprintf('Start time: %.3f s, target laps: %d\n', cfg.analysis_start_time, cfg.target_laps);
fprintf('Sliding windows: %d laps, step %d lap(s)\n', cfg.analysis_win_size, cfg.sliding_step);
fprintf('Old SG result: %s\n', old_result_file);

if ~isfile(template_file)
    error('Template file not found. Run Step01 first: %s', template_file);
end
if ~isfile(old_result_file)
    error('Old SG result file not found: %s', old_result_file);
end

loaded_template = load(template_file, 'Template');
Template = loaded_template.Template;
loaded_result = load(old_result_file, 'Result_Struct');
OldResult = loaded_result.Result_Struct;
current_cfg = build_single_sync_experiment_config_20250527();

% Use the current workspace config for OPR timing and raw-data paths.  The
% old SG result is retained only as a trend/window reference.
cfg.dynamic_data_dir = current_cfg.dynamic_data_dir;
cfg.case_output_dir = current_cfg.case_output_dir;
cfg.sensor_config_file = current_cfg.sensor_config_file;
cfg.opr_channel = current_cfg.opr_channel;
cfg.opr_pulses_per_rev = current_cfg.opr_pulses_per_rev;
cfg.pinlv = current_cfg.pinlv;
cfg.r_tip_mm = current_cfg.r_tip_mm;
cfg.pulse_window_sec = current_cfg.pulse_window_sec;
cfg.pulse_pad_sec = current_cfg.pulse_pad_sec;
cfg.dynamic_window_mode = current_cfg.dynamic_window_mode;

loaded_sensor_config = load(cfg.sensor_config_file, 'Sensor_Config');
Sensor_Config = loaded_sensor_config.Sensor_Config;
loaded_opr = load(fullfile(cfg.case_output_dir, 'jiluOPR.mat'), 'jiluOPR');
opr_times = loaded_opr.jiluOPR(:, 1);
opr_reference = build_opr_reference_from_jilu_local(loaded_opr.jiluOPR, cfg.opr_pulses_per_rev, cfg.r_tip_mm);
F_omega_deg = build_phase_speed_local(opr_times, cfg.opr_pulses_per_rev);
fprintf('Dynamic OPR center shift: %.4f deg, %.4f mm from rising edge.\n', ...
    opr_reference.phase_shift_deg, opr_reference.phase_shift_mm);

%% Select the same 20 target-blade laps after the same start time
selection = repmat(struct('sensor_id', NaN, 'selected_rows', []), numel(cfg.analysis_sensors), 1);
global_window = [inf, -inf];
for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_file = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;

    blade_mask = jilublade(:, 4) == cfg.target_blades & jilublade(:, 3) >= cfg.analysis_start_time;
    blade_rows = find(blade_mask);
    if numel(blade_rows) < cfg.target_laps
        error('CH%d has only %d target-blade laps after %.3f s.', sid, numel(blade_rows), cfg.analysis_start_time);
    end
    blade_rows = blade_rows(1:cfg.target_laps);

    selection(is).sensor_id = sid;
    selection(is).selected_rows = blade_rows(:);
    global_window(1) = min(global_window(1), min(jilublade(blade_rows, 1)) - cfg.pulse_window_sec);
    global_window(2) = max(global_window(2), max(jilublade(blade_rows, 2)) + cfg.pulse_window_sec);
end
fprintf('Global raw-data window: %.6f-%.6f s\n', global_window(1), global_window(2));

%% Load raw dynamic voltage streams only once
file_ranges = build_dynamic_file_ranges_local(cfg.dynamic_data_dir, cfg.opr_channel, cfg.pinlv);
selected_file_mask = [file_ranges.t_end] >= global_window(1) & [file_ranges.t_start] <= global_window(2);
selected_file_ranges = file_ranges(selected_file_mask);
if isempty(selected_file_ranges)
    error('No raw dynamic files overlap the requested window.');
end

raw_stream(max(cfg.analysis_sensors)) = struct('T', [], 'V', []);
for ir = 1:numel(selected_file_ranges)
    file_id = selected_file_ranges(ir).file_id;
    offset = selected_file_ranges(ir).offset;
    for sid = cfg.analysis_sensors
        [t_local, v_local] = load_raw_case_channel_local(cfg.dynamic_data_dir, sid, file_id, cfg.pinlv);
        if isempty(t_local)
            continue;
        end
        t_global = t_local(:) + offset;
        keep = t_global >= global_window(1) & t_global <= global_window(2);
        raw_stream(sid).T = [raw_stream(sid).T; t_global(keep)];
        raw_stream(sid).V = [raw_stream(sid).V; v_local(keep)];
    end
end

%% Extract each selected lap waveform once
LapData = repmat(struct('sensor_id', NaN, 'Lap', []), numel(cfg.analysis_sensors), 1);
for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    probe_file = fullfile(cfg.case_output_dir, sprintf('jilublade_probe%d.mat', sid));
    loaded_probe = load(probe_file, 'jilublade');
    jilublade = loaded_probe.jilublade;
    theta_std = Sensor_Config.Standard_Relative_Angles(sid, cfg.target_blades);
    theta_std = convert_standard_angle_to_opr_center_local(theta_std, opr_reference);

    Lap = repmat(struct('lap_id', NaN, 't', [], 'x_abs', [], 'V', [], 'theta', []), cfg.target_laps, 1);
    for lap_id = 1:cfg.target_laps
        row_id = selection(is).selected_rows(lap_id);
        t_peak = jilublade(row_id, 3);
        if strcmpi(cfg.dynamic_window_mode, 'legacy_row_bounds')
            t_start = jilublade(row_id, 1) - cfg.pulse_pad_sec;
            t_end = jilublade(row_id, 2) + cfg.pulse_pad_sec;
        else
            t_start = t_peak - cfg.pulse_window_sec;
            t_end = t_peak + cfg.pulse_window_sec;
        end

        mask = raw_stream(sid).T >= t_start & raw_stream(sid).T <= t_end;
        t_seg = raw_stream(sid).T(mask);
        v_seg = raw_stream(sid).V(mask);
        idx_prev = find(opr_times < t_peak, 1, 'last');
        if isempty(idx_prev) || numel(t_seg) < 5
            continue;
        end

        theta_points_deg = map_segment_to_relative_angle_local(opr_times(idx_prev), t_seg, F_omega_deg);
        theta_diff_deg = mod(theta_points_deg - theta_std + 180, 360) - 180;
        x_abs = theta_diff_deg * (pi / 180) * cfg.r_tip_mm;
        theta_rot = map_time_to_rotor_phase_local(opr_times, t_seg, cfg.opr_pulses_per_rev);
        valid = isfinite(theta_rot);

        Lap(lap_id).lap_id = lap_id;
        Lap(lap_id).t = t_seg(valid);
        Lap(lap_id).x_abs = x_abs(valid);
        Lap(lap_id).V = v_seg(valid);
        Lap(lap_id).theta = theta_rot(valid);
    end

    LapData(is).sensor_id = sid;
    LapData(is).Lap = Lap;
end

%% Build all sliding-window maps
num_windows = floor((cfg.target_laps - cfg.analysis_win_size) / cfg.sliding_step) + 1;
Window = repmat(struct( ...
    'window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, ...
    'Sensor', []), num_windows, 1);

for w_idx = 1:num_windows
    lap_start = 1 + (w_idx - 1) * cfg.sliding_step;
    lap_end = lap_start + cfg.analysis_win_size - 1;
    lap_range = lap_start:lap_end;
    Sensor = repmat(struct( ...
        'sensor_id', NaN, 't', [], 'x_abs', [], 'x_rel', [], ...
        'V', [], 'W', [], 'theta', [], 'point_count', NaN), numel(cfg.analysis_sensors), 1);

    all_t_window = [];
    for is = 1:numel(cfg.analysis_sensors)
        sid = cfg.analysis_sensors(is);
        Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
        t = [];
        x_abs = [];
        V = [];
        theta = [];
        for lap_id = lap_range
            D = LapData(is).Lap(lap_id);
            t = [t; D.t(:)]; %#ok<AGROW>
            x_abs = [x_abs; D.x_abs(:)]; %#ok<AGROW>
            V = [V; D.V(:)]; %#ok<AGROW>
            theta = [theta; D.theta(:)]; %#ok<AGROW>
        end
        x_rel = x_abs - Tpl.xc;
        W = build_simple_waveform_weight_local(V);

        Sensor(is).sensor_id = sid;
        Sensor(is).t = t(:);
        Sensor(is).x_abs = x_abs(:);
        Sensor(is).x_rel = x_rel(:);
        Sensor(is).V = V(:);
        Sensor(is).W = W(:);
        Sensor(is).theta = theta(:);
        Sensor(is).point_count = numel(t);
        all_t_window = [all_t_window; t(:)]; %#ok<AGROW>
    end

    Window(w_idx).window_id = w_idx;
    Window(w_idx).lap_range = lap_range;
    Window(w_idx).time_window = [min(all_t_window), max(all_t_window)];
    Window(w_idx).rot_freq_mean_hz = compute_local_rot_freq_local(opr_times, cfg.opr_pulses_per_rev, Window(w_idx).time_window);
    Window(w_idx).rot_rpm_mean = 60 * Window(w_idx).rot_freq_mean_hz;
    Window(w_idx).Sensor = Sensor;
end

DynamicMap = struct();
DynamicMap.Route = 'low_speed_template_sliding_window_dynamic_map_no_gap';
DynamicMap.TargetBlade = cfg.target_blades;
DynamicMap.SensorIDs = cfg.analysis_sensors;
DynamicMap.SensorTag = sensor_tag;
DynamicMap.SourceSettings = cfg;
DynamicMap.SourceResultFile = old_result_file;
DynamicMap.SourceResultSensorTag = source_result_sensor_tag;
DynamicMap.TemplateFileForXRel = template_file;
DynamicMap.TemplateSuffixForXRel = template_suffix;
DynamicMap.XCenterBySensor = build_xcenter_table_local(Template, cfg.analysis_sensors);
DynamicMap.OPRReference = opr_reference;
if strcmp(sensor_tag, source_result_sensor_tag)
    DynamicMap.OldSGTrends = OldResult.Trends;
    DynamicMap.OldSGBestWindowID = OldResult.BestWindow.window_id;
    DynamicMap.OldSGBestLapRange = OldResult.BestWindow.lap_range;
end
DynamicMap.GlobalTimeWindow = global_window;
DynamicMap.SelectedRawFileIDs = [selected_file_ranges.file_id];
DynamicMap.Window = Window;

save(dynamic_file, 'DynamicMap', '-v7.3');
fprintf('Saved sliding-window dynamic map: %s\n', dynamic_file);
fprintf('Built %d windows.\n', num_windows);
if isfield(DynamicMap, 'OldSGBestWindowID')
    fprintf('Old SG best window: %d, laps %s\n', ...
        DynamicMap.OldSGBestWindowID, mat2str(DynamicMap.OldSGBestLapRange));
else
    fprintf('Old SG comparison skipped: source result is %s, current map is %s.\n', ...
        source_result_sensor_tag, sensor_tag);
end

%% Visualization of one sliding window for checking
if isfield(DynamicMap, 'OldSGBestWindowID')
    best_plot_idx = DynamicMap.OldSGBestWindowID;
else
    best_plot_idx = 1;
end
fig = figure('Name', 'Step02 sliding-window dynamic map check', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.08 0.88 0.78]);
tiledlayout(fig, numel(cfg.analysis_sensors), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(cfg.analysis_sensors)
    sid = cfg.analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = DynamicMap.Window(best_plot_idx).Sensor(is);

    nexttile;
    plot(Tpl.x_grid, Tpl.v_grid, 'k-', 'LineWidth', 1.8, 'DisplayName', 'low-speed spline template'); hold on;
    scatter(D.x_rel, D.V, 8, D.t, 'filled', 'DisplayName', 'dynamic samples');
    xlabel('x relative to template center (mm)');
    ylabel('Voltage (V)');
    title(sprintf('CH%d dynamic samples, window %d, laps %s', ...
        sid, best_plot_idx, mat2str(DynamicMap.Window(best_plot_idx).lap_range)));
    cb = colorbar;
    cb.Label.String = 'Time (s)';
    box on; set(gca, 'TickDir', 'in', 'FontName', 'Times New Roman');
    legend('Location', 'best');
end

exportgraphics(fig, fullfile(figure_dir, sprintf('Step02_DynamicMap_B%d_%s_SlidingWindows%s.png', ...
    cfg.target_blades, sensor_tag, dynamic_suffix)), 'Resolution', 300);

%% Local helpers kept at the end so the main workflow above stays readable
function v = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    error('%s must be a positive integer.', name);
end
v = max(1, floor(tmp));
end

function suffix = sanitize_dynamic_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function T = build_xcenter_table_local(Template, sensor_ids)
sensor_id = sensor_ids(:);
xc_mm = nan(numel(sensor_id), 1);
for i = 1:numel(sensor_id)
    idx = find([Template.Sensor.sensor_id] == sensor_id(i), 1, 'first');
    if ~isempty(idx) && isfield(Template.Sensor(idx), 'xc')
        xc_mm(i) = Template.Sensor(idx).xc;
    end
end
T = table(sensor_id, xc_mm);
end

function theta_std_center = convert_standard_angle_to_opr_center_local(theta_std_start, opr_reference)
theta_std_center = theta_std_start;
if isstruct(opr_reference) && isfield(opr_reference, 'phase_shift_deg') && ...
        isfinite(opr_reference.phase_shift_deg)
    theta_std_center = theta_std_start - opr_reference.phase_shift_deg;
end
end

function opr_reference = build_opr_reference_from_jilu_local(jiluOPR, pulses_per_rev, r_tip_mm)
opr_reference = struct('mode', 'multi_threshold_center', ...
    'standard_angle_reference', 'opr_pulse_center', ...
    'phase_shift_deg', 0, ...
    'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if size(jiluOPR, 2) < 2 || pulses_per_rev < 1
    return;
end
center_time = jiluOPR(:, 1);
if size(jiluOPR, 2) >= 3
    start_time = jiluOPR(:, 2);
else
    start_time = jiluOPR(:, 2);
    if median(center_time - start_time, 'omitnan') < 0
        warning(['Two-column jiluOPR appears to store [start,end], not [center,start]. ' ...
            'OPR center phase shift is left at zero.']);
        return;
    end
end
n = min(numel(center_time) - pulses_per_rev, numel(start_time));
if n < 1
    return;
end
dt_center = center_time(1:n) - start_time(1:n);
dt_rev = center_time((1:n) + pulses_per_rev) - center_time(1:n);
valid = isfinite(dt_center) & isfinite(dt_rev) & dt_rev > eps;
if ~any(valid)
    return;
end
shift_deg = 360 * dt_center(valid) ./ dt_rev(valid);
opr_reference.phase_shift_deg = median(shift_deg, 'omitnan');
opr_reference.phase_shift_mm = opr_reference.phase_shift_deg * (pi / 180) * r_tip_mm;
opr_reference.median_center_minus_start_s = median(dt_center(valid), 'omitnan');
end

function file_ranges = build_dynamic_file_ranges_local(case_dir, opr_channel, pinlv)
d = dir(fullfile(case_dir, sprintf('4-%d-*.mat', opr_channel)));
file_ids = nan(numel(d), 1);
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(opr_channel) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        file_ids(i) = str2double(tok{1});
    end
end
file_ids = sort(unique(file_ids(~isnan(file_ids))));

file_ranges = repmat(struct('file_id', NaN, 'offset', NaN, 't_start', NaN, 't_end', NaN), numel(file_ids), 1);
last_end = [];
for i = 1:numel(file_ids)
    file_id = file_ids(i);
    [t_opr, ~] = load_raw_case_channel_local(case_dir, opr_channel, file_id, pinlv);
    if isempty(t_opr)
        continue;
    end
    if isempty(last_end)
        offset = 0;
    else
        offset = last_end + 1 / pinlv - t_opr(1);
    end
    last_end = t_opr(end) + offset;
    file_ranges(i).file_id = file_id;
    file_ranges(i).offset = offset;
    file_ranges(i).t_start = t_opr(1) + offset;
    file_ranges(i).t_end = t_opr(end) + offset;
end
end

function [t_sec, v] = load_raw_case_channel_local(case_dir, sid, file_id, pinlv)
filepath = fullfile(case_dir, sprintf('4-%d-%d.mat', sid, file_id));
if ~isfile(filepath)
    t_sec = [];
    v = [];
    return;
end
loaded = load(filepath);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
t_sec = raw(:, 1) / pinlv;
v = raw(:, 2);
end

function F_omega_deg = build_phase_speed_local(opr_times, blades_num)
spd_t = opr_times(1:end-blades_num);
spd_v = 360 ./ max(opr_times(blades_num+1:end) - opr_times(1:end-blades_num), eps);
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');
end

function theta_points = map_segment_to_relative_angle_local(t_ref, t_seg, F_omega_deg)
dt_first = linspace(t_ref, t_seg(1), 10);
theta_base = trapz(dt_first, F_omega_deg(dt_first));
w_seg = F_omega_deg(t_seg);
theta_rel = cumtrapz(t_seg, w_seg);
theta_points = theta_base + theta_rel;
end

function theta_rot = map_time_to_rotor_phase_local(opr_times, sample_times, num_blades)
theta_rot = nan(size(sample_times));
for i = 1:numel(opr_times)-1
    t0 = opr_times(i);
    t1 = opr_times(i+1);
    mask = sample_times >= t0 & sample_times <= t1;
    if ~any(mask)
        continue;
    end
    frac = (sample_times(mask) - t0) ./ max(t1 - t0, eps);
    theta_rot(mask) = (2*pi/num_blades) * ((i - 1) + frac);
end
end

function W = build_simple_waveform_weight_local(V)
if isempty(V)
    W = [];
    return;
end
V = V(:);
v_floor = prctile(V, 5);
v_peak = prctile(V, 99);
span = max(v_peak - v_floor, eps);
W = (V - v_floor) ./ span;
W = min(max(W, 0.05), 1.0);
end

function rot_freq_hz = compute_local_rot_freq_local(opr_times, pulses_per_rev, time_window)
mask = opr_times >= time_window(1) & opr_times <= time_window(2);
t = opr_times(mask);
if numel(t) > pulses_per_rev
    rot_freq_hz = median(1 ./ max(t(1+pulses_per_rev:end) - t(1:end-pulses_per_rev), eps), 'omitnan');
else
    rot_freq_hz = NaN;
end
end
end

function step03_identify_direct_template_embedded(rootDir)
%% Step03_Main_VPTop3SynchronousWaveform_20250527
% Main low-speed-template-only identification route.
% First-order template VP screens EO candidates; the final result still
% comes from the full synchronous waveform objective.
%
% Every integer EO candidate is scored by the first-order model
% V - T(x) ~= -T'(x) * [dx_c + a*sin(EO*theta) + b*cos(EO*theta)].
% By default, the top-3 EO candidates enter full waveform refinement in
% every window. For synchronous vibration, frequency is constrained by
% f = EO * rot_freq_mean; the final refinement optimizes amplitude/phase
% and the spatial alignment term dx_c.

%% Settings matching the super-Gaussian Step5 route
route_dir = rootDir;

cfg = struct();
cfg.target_blades = 1;
cfg.analysis_sensors = [1, 3];
sensor_override = strtrim(getenv('STEP03_ANALYSIS_SENSORS'));
if ~isempty(sensor_override)
    parsed_sensors = sscanf(sensor_override, '%d').';
    if isempty(parsed_sensors)
        error('STEP03_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysis_sensors = parsed_sensors;
end
cfg.analysis_start_time = 1.5;
cfg.target_laps = 20;
cfg.analysis_win_size = 3;
cfg.sliding_step = 1;
cfg.freq_search_hz = [100, 1000];
cfg.eo_pad = 2;
cfg.amplitude_limit_mm = 0.50;
cfg.dx_c_limit_mm = 0.20;
cfg.debug_max_windows = parse_positive_integer_env_local('STEP03_04_MAX_WINDOWS', inf);
cfg.vp_top_k_eo = parse_positive_integer_env_local('STEP03_MAIN_TOP_K_EO', 3);
cfg.vp_all_eo_warmup_windows = parse_nonnegative_integer_env_local('STEP03_04_ALL_EO_WARMUP_WINDOWS', 0);
cfg.vp_gap_ratio_fallback = parse_numeric_env_local('STEP03_04_GAP_RATIO_FALLBACK', -inf);
cfg.vp_linear_gap_ratio_fallback = parse_numeric_env_local('STEP03_04_LINEAR_GAP_RATIO_FALLBACK', -inf);
cfg.diagnostic_eo = parse_optional_integer_env_local('STEP03_DIAGNOSTIC_EO');
cfg.sensor_eta_limit_mm = parse_nonnegative_numeric_env_local('STEP03D_SENSOR_ETA_LIMIT_MM', 0.08);
cfg.sensor_eta_reg_weight_v_per_mm = parse_nonnegative_numeric_env_local('STEP03D_SENSOR_ETA_REG_WEIGHT_V_PER_MM', 1.00);
cfg.pulse_selection_mode = lower(strtrim(getenv('STEP03D_PULSE_MODE')));
if isempty(cfg.pulse_selection_mode)
    cfg.pulse_selection_mode = 'single';
end
if ~ismember(cfg.pulse_selection_mode, {'single', 'all'})
    error('STEP03D_PULSE_MODE must be "single" or "all".');
end

method = struct();
method.name = 'template_only_main_vp_top3_synchronous_waveform';
method.selection_rule = 'first_order_template_vp_top3_then_synchronous_waveform_rmse';
method.use_reference_eo_constraint = false;
method.use_vp = true;
method.vp_top_k_eo = cfg.vp_top_k_eo;
method.refine_params = {'A', 'phi', 'dx_c', 'sensor_eta'};
method.template_forward = 'interp1_low_speed_template';
method.vp_seed_model = 'V_minus_Tx_equals_minus_Tprime_times_u';
method.final_waveform_objective = 'fixed_eo_direct_template_voltage_residual_without_sensor_affine_projection';
method.sensor_eta_limit_mm = cfg.sensor_eta_limit_mm;
method.sensor_eta_reg_weight_v_per_mm = cfg.sensor_eta_reg_weight_v_per_mm;
method.weight_floor = 0.05;
method.domain_margin_mm = 0.02;
method.coverage_safety_margin_mm = 0.05;
method.query_guard_mm = cfg.amplitude_limit_mm + cfg.dx_c_limit_mm + method.coverage_safety_margin_mm;
method.query_guard_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MM', method.query_guard_mm);
method.query_guard_mode = lower(strtrim(getenv('STEP03D_QUERY_GUARD_MODE')));
if isempty(method.query_guard_mode)
    method.query_guard_mode = 'adaptive';
end
if ~ismember(method.query_guard_mode, {'fixed', 'adaptive'})
    error('STEP03D_QUERY_GUARD_MODE must be "fixed" or "adaptive".');
end
method.query_guard_quantile = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_QUANTILE', 95);
method.query_guard_safety_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_SAFETY_MM', method.coverage_safety_margin_mm);
method.query_guard_min_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MIN_MM', 0.12);
method.query_guard_max_mm = parse_nonnegative_numeric_env_local('STEP03D_QUERY_GUARD_MAX_MM', method.query_guard_mm);
method.domain_selection_mode = lower(strtrim(getenv('STEP03D_DOMAIN_SELECTION_MODE')));
if isempty(method.domain_selection_mode)
    method.domain_selection_mode = 'hard';
end
if ~ismember(method.domain_selection_mode, {'hard', 'soft'})
    error('STEP03D_DOMAIN_SELECTION_MODE must be "hard" or "soft".');
end
method.domain_soft_margin_mm = parse_nonnegative_numeric_env_local('STEP03D_DOMAIN_SOFT_MARGIN_MM', 0);
method.overshoot_penalty_weight = parse_nonnegative_numeric_env_local('STEP03D_OVERSHOOT_PENALTY_WEIGHT', 100);
method.coverage_mode = 'prefer_points_safe_for_all_bounded_u';
method.extrapolation_mode = 'clamp_to_template_edge';
method.default_sensor_threshold = 0.5;
method.pulse_selection_mode = cfg.pulse_selection_mode;
method.dynamic_effective_mode = lower(strtrim(getenv('STEP03D_DYNAMIC_EFFECTIVE_MODE')));
if isempty(method.dynamic_effective_mode)
    method.dynamic_effective_mode = 'gradient';
end
if ~ismember(method.dynamic_effective_mode, {'legacy', 'gradient'})
    error('STEP03D_DYNAMIC_EFFECTIVE_MODE must be "legacy" or "gradient".');
end
method.dynamic_template_gradient_min_ratio = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_TEMPLATE_GRADIENT_MIN_RATIO', 0.08);
method.dynamic_time_gradient_min_ratio = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_TIME_GRADIENT_MIN_RATIO', 0.15);
method.dynamic_peak_quantile = parse_nonnegative_numeric_env_local('STEP03D_DYNAMIC_PEAK_QUANTILE', 85);

sensor_tag = ['S', sprintf('%d', cfg.analysis_sensors)];
output_dir = fullfile(route_dir, 'output');
template_dir = fullfile(output_dir, 'templates');
dynamic_dir = fullfile(output_dir, 'dynamic_maps');
result_dir = fullfile(output_dir, 'identification');
if exist(result_dir, 'dir') ~= 7; mkdir(result_dir); end

template_file_override = strtrim(getenv('STEP03_TEMPLATE_FILE'));
if ~isempty(template_file_override)
    if exist(template_file_override, 'file') ~= 2
        error('STEP03_TEMPLATE_FILE does not exist: %s', template_file_override);
    end
    template_file = template_file_override;
else
    template_file = find_template_file_for_sensors_local(template_dir, cfg.target_blades, cfg.analysis_sensors);
end
dynamic_file_override = strtrim(getenv('STEP03_DYNAMIC_MAP_FILE'));
if ~isempty(dynamic_file_override)
    if exist(dynamic_file_override, 'file') ~= 2
        error('STEP03_DYNAMIC_MAP_FILE does not exist: %s', dynamic_file_override);
    end
    dynamic_file = dynamic_file_override;
else
    dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, cfg.target_blades, cfg.analysis_sensors);
end
result_suffix = sanitize_result_suffix_local(strtrim(getenv('STEP03_RESULT_SUFFIX')));
result_file = fullfile(result_dir, sprintf('Result_Step03_Main_VPTop3SynchronousWaveform_B%d_%s%s_20250527.mat', ...
    cfg.target_blades, sensor_tag, result_suffix));

loaded_template = load(template_file, 'Template');
loaded_dynamic = load(dynamic_file, 'DynamicMap');
Template = filter_template_sensors_local(loaded_template.Template, cfg.analysis_sensors, sensor_tag);
DynamicMap = filter_dynamic_map_sensors_local(loaded_dynamic.DynamicMap, cfg.analysis_sensors, sensor_tag);
coordinate_check = check_template_dynamic_xcenter_local(Template, DynamicMap, cfg.analysis_sensors, 1e-6);
source_cfg = DynamicMap.SourceSettings;
cfg.analysis_start_time = source_cfg.analysis_start_time;
cfg.target_laps = source_cfg.target_laps;
cfg.analysis_win_size = source_cfg.analysis_win_size;
cfg.sliding_step = source_cfg.sliding_step;

fprintf('\n=== Step03 main: VP top-3 synchronous waveform identification ===\n');
fprintf('Target blade: %d\n', cfg.target_blades);
fprintf('Analysis sensors: %s\n', mat2str(cfg.analysis_sensors));
fprintf('Template source: %s\n', template_file);
fprintf('Dynamic map source: %s\n', dynamic_file);
fprintf('Coordinate check: max |Template xc - DynamicMap xc| = %.6g mm (%s).\n', ...
    coordinate_check.max_abs_delta_mm, coordinate_check.status);
legacy_base_xc_gradient_tpl = contains(template_file, 'GradientXRange030') && ...
    ~contains(dynamic_file, 'GradientXRange030');
if ~coordinate_check.is_consistent && legacy_base_xc_gradient_tpl
    fprintf('Coordinate note: legacy 20250527 route uses base-xc DynamicMap with GradientXRange030 template.\n');
elseif ~coordinate_check.is_consistent
    warning('Template/DynamicMap x-center mismatch. This run is not a clean GradientXRange030 coordinate chain.');
end
fprintf('Frequency search: %.1f-%.1f Hz; no reference EO; first-order VP only screens candidates.\n', ...
    cfg.freq_search_hz(1), cfg.freq_search_hz(2));
fprintf('Model: V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi)); synchronous f = EO*rot_freq_mean.\n');
fprintf(['Template coverage guard: %s, fixed %.3f mm, min %.3f mm, max %.3f mm, ' ...
    'q%.1f + %.3f mm.\n'], method.query_guard_mode, method.query_guard_mm, ...
    method.query_guard_min_mm, method.query_guard_max_mm, method.query_guard_quantile, ...
    method.query_guard_safety_mm);
fprintf('Pulse selection mode: %s.\n', method.pulse_selection_mode);
fprintf('Dynamic effective selection: %s.\n', method.dynamic_effective_mode);
fprintf('Domain selection: %s; soft margin %.3f mm; overshoot penalty %.3g.\n', ...
    method.domain_selection_mode, method.domain_soft_margin_mm, method.overshoot_penalty_weight);
fprintf('First-order VP/adaptive: keep top %d EO candidates in every window by default.\n', ...
    method.vp_top_k_eo);
fprintf('Sensor eta limit: %.4f mm (no EO prior; per-sensor x-zero correction).\n', ...
    method.sensor_eta_limit_mm);
fprintf('Sensor eta regularization: %.4f V/mm.\n', cfg.sensor_eta_reg_weight_v_per_mm);
if cfg.vp_all_eo_warmup_windows > 0 || isfinite(cfg.vp_gap_ratio_fallback) || isfinite(cfg.vp_linear_gap_ratio_fallback)
    fprintf('Optional fallback enabled: warmup=%d, weighted gap=%.6g, linear gap=%.6g.\n', ...
        cfg.vp_all_eo_warmup_windows, cfg.vp_gap_ratio_fallback, cfg.vp_linear_gap_ratio_fallback);
else
    fprintf('Optional fallback disabled: no warmup and no all-EO ambiguity fallback.\n');
end
if isfinite(cfg.diagnostic_eo)
    fprintf('Diagnostic only: report EO%d rank after VP screening and full waveform refinement.\n', cfg.diagnostic_eo);
end

num_windows = numel(DynamicMap.Window);
if isfinite(cfg.debug_max_windows)
    num_windows = min(num_windows, cfg.debug_max_windows);
end

trend_rows = repmat(struct( ...
    'window_id', NaN, 'lap_start', NaN, 'lap_end', NaN, 'window_center_time', NaN, ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, 'A_id', NaN, 'EO_id', NaN, ...
    'fn_id', NaN, 'phi_id_wrapped', NaN, 'dx_c_id', NaN, 'd0_id', NaN, ...
    'sensor_eta_max_abs_mm', NaN, ...
    'weighted_voltage_rmse', NaN, 'plain_voltage_rmse', NaN, ...
    'valid_segment_count', NaN, 'point_count', NaN), num_windows, 1);
window_results = repmat(struct(), num_windows, 1);
best_window = struct('weighted_voltage_rmse', inf);

for w_idx = 1:num_windows
    Wmap = DynamicMap.Window(w_idx);
    bundle = build_template_observation_bundle_local(Wmap, Template, cfg.analysis_sensors, method);
    eo_candidates = build_eo_candidates_local(bundle.rot_freq_mean_hz, cfg.freq_search_hz, cfg.eo_pad);
    seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates);
    [final_eo_candidates, selection_info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, w_idx);
    result = refine_template_waveform_fit_local(bundle, seed_table, final_eo_candidates, cfg);
    result.VPSeedTable = struct2table(seed_table);
    result.VPSelectedEO = final_eo_candidates(:).';
    result.VPSelectionInfo = selection_info;

    trend_rows(w_idx).window_id = w_idx;
    trend_rows(w_idx).lap_start = Wmap.lap_range(1);
    trend_rows(w_idx).lap_end = Wmap.lap_range(end);
    trend_rows(w_idx).window_center_time = median(bundle.T, 'omitnan');
    trend_rows(w_idx).rot_freq_mean_hz = bundle.rot_freq_mean_hz;
    trend_rows(w_idx).rot_rpm_mean = Wmap.rot_rpm_mean;
    trend_rows(w_idx).A_id = result.A_id;
    trend_rows(w_idx).EO_id = result.EO_id;
    trend_rows(w_idx).fn_id = result.fn_id;
    trend_rows(w_idx).phi_id_wrapped = result.phi_id_wrapped;
    trend_rows(w_idx).dx_c_id = result.dx_c_id;
    trend_rows(w_idx).d0_id = result.dx_c_id;
    trend_rows(w_idx).sensor_eta_max_abs_mm = max(abs(result.sensor_eta_id), [], 'omitnan');
    trend_rows(w_idx).weighted_voltage_rmse = result.weighted_voltage_rmse;
    trend_rows(w_idx).plain_voltage_rmse = result.plain_voltage_rmse;
    trend_rows(w_idx).valid_segment_count = result.valid_segment_count;
    trend_rows(w_idx).point_count = result.point_count;

    window_results(w_idx).window_id = w_idx;
    window_results(w_idx).lap_range = Wmap.lap_range;
    window_results(w_idx).time_window = Wmap.time_window;
    window_results(w_idx).bundle = bundle;
    window_results(w_idx).seed_table = seed_table;
    window_results(w_idx).CandidateTable = result.CandidateTable;
    window_results(w_idx).Result = result;

    if result.weighted_voltage_rmse < best_window.weighted_voltage_rmse
        best_window = result;
        best_window.window_id = w_idx;
        best_window.lap_range = Wmap.lap_range;
        best_window.seed_table = seed_table;
    end

    diagnostic_text = '';
    if isfinite(cfg.diagnostic_eo)
        diag_seed_rank = find([seed_table.EO] == cfg.diagnostic_eo, 1);
        diag_final_rank = find(result.CandidateTable.EO == cfg.diagnostic_eo, 1);
        diagnostic_text = sprintf(', EO%d seed/final-rank=%s/%s', ...
            cfg.diagnostic_eo, rank_to_string_local(diag_seed_rank), rank_to_string_local(diag_final_rank));
    end
    fprintf('Window %02d/%02d laps %s: EO=%d, f=%.3f Hz, A=%.4f mm, RMSE=%.5f V, clamp=%.2f%%%s\n', ...
        w_idx, num_windows, mat2str(Wmap.lap_range), result.EO_id, result.fn_id, ...
        result.A_id, result.weighted_voltage_rmse, 100 * result.Coverage.clamp_fraction, diagnostic_text);
end

Trend = struct2table(trend_rows);
Result = struct();
Result.Method = method.name;
Result.AnalysisSettings = cfg;
Result.MethodSettings = method;
Result.TemplateFile = template_file;
Result.DynamicMapFile = dynamic_file;
Result.CoordinateCheck = coordinate_check;
Result.TargetBlade = cfg.target_blades;
Result.SensorIDs = cfg.analysis_sensors;
Result.SensorTag = sensor_tag;
Result.Trend = Trend;
Result.WindowResult = window_results;
Result.BestWindow = best_window;
Result.ResonanceSummary = build_resonance_summary_local(trend_rows);

save(result_file, 'Result', '-v7.3');

fprintf('\n=== Step03 main trend ===\n');
disp(Trend(:, {'window_id','lap_start','lap_end','EO_id','fn_id','A_id','weighted_voltage_rmse'}));
fprintf('Dominant EO: %d; mean f = %.6f Hz; median f = %.6f Hz\n', ...
    Result.ResonanceSummary.dominant_eo, Result.ResonanceSummary.mean_freq_hz, ...
    Result.ResonanceSummary.median_freq_hz);
fprintf('Saved result: %s\n', result_file);

%% Local functions
function v = parse_positive_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = default_value;
else
    v = max(1, floor(tmp));
end
end

function v = parse_nonnegative_integer_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp < 0
    v = default_value;
else
    v = floor(tmp);
end
end

function v = parse_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp)
    v = default_value;
else
    v = tmp;
end
end

function v = parse_nonnegative_numeric_env_local(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    v = default_value;
    return;
end
parsed = str2double(raw);
if ~isfinite(parsed) || parsed < 0
    error('%s must be a nonnegative numeric value.', name);
end
v = parsed;
end

function v = parse_optional_integer_env_local(name)
raw = strtrim(getenv(name));
if isempty(raw)
    v = NaN;
    return;
end
tmp = str2double(raw);
if ~isfinite(tmp) || tmp <= 0
    v = NaN;
else
    v = floor(tmp);
end
end

function suffix = sanitize_result_suffix_local(raw)
if isempty(raw)
    suffix = '';
    return;
end
suffix = regexprep(raw, '[^A-Za-z0-9_\\-]', '_');
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
end

function template_file = find_template_file_for_sensors_local(template_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20250527.mat', target_blade, sensor_tag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', target_blade, sensor_tag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', target_blade, sensor_tag)
    };
for ip = 1:numel(patterns)
    candidate = fullfile(template_dir, patterns{ip});
    if exist(candidate, 'file') == 2
        template_file = candidate;
        return;
    end
end
files = dir(fullfile(template_dir, sprintf('Template_LowSpeedRotating_B%d_*.mat', target_blade)));
if isempty(files)
    error('No low-speed template file found under %s.', template_dir);
end
error('Template for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function dynamic_file = find_dynamic_file_for_sensors_local(dynamic_dir, target_blade, analysis_sensors)
sensor_tag = ['S', sprintf('%d', analysis_sensors)];
patterns = {
    sprintf('DynamicMap_B%d_%s_SlidingWindows_20250527.mat', target_blade, sensor_tag)
    sprintf('DynamicMap_B%d_%s_SlidingWindows_GradientXRange030_20250527.mat', target_blade, sensor_tag)
    };
for ip = 1:numel(patterns)
    candidate = fullfile(dynamic_dir, patterns{ip});
    if exist(candidate, 'file') == 2
        dynamic_file = candidate;
        return;
    end
end
files = dir(fullfile(dynamic_dir, sprintf('DynamicMap_B%d_*.mat', target_blade)));
if isempty(files)
    error('No dynamic map file found under %s.', dynamic_dir);
end
error('DynamicMap for sensor tag %s not found. Available first file: %s', sensor_tag, files(1).name);
end

function Template = filter_template_sensors_local(Template, analysis_sensors, sensor_tag)
available = [Template.Sensor.sensor_id];
sensor_idx = zeros(size(analysis_sensors));
for is = 1:numel(analysis_sensors)
    idx = find(available == analysis_sensors(is), 1);
    if isempty(idx)
        error('Template does not contain sensor %d.', analysis_sensors(is));
    end
    sensor_idx(is) = idx;
end
Template.Sensor = Template.Sensor(sensor_idx);
Template.SensorIDs = analysis_sensors;
Template.SensorTag = sensor_tag;
end

function DynamicMap = filter_dynamic_map_sensors_local(DynamicMap, analysis_sensors, sensor_tag)
for iw = 1:numel(DynamicMap.Window)
    available = [DynamicMap.Window(iw).Sensor.sensor_id];
    sensor_idx = zeros(size(analysis_sensors));
    for is = 1:numel(analysis_sensors)
        idx = find(available == analysis_sensors(is), 1);
        if isempty(idx)
            error('DynamicMap window %d does not contain sensor %d.', iw, analysis_sensors(is));
        end
        sensor_idx(is) = idx;
    end
    DynamicMap.Window(iw).Sensor = DynamicMap.Window(iw).Sensor(sensor_idx);
end
DynamicMap.SensorIDs = analysis_sensors;
DynamicMap.SensorTag = sensor_tag;
end

function check = check_template_dynamic_xcenter_local(Template, DynamicMap, analysis_sensors, tolerance_mm)
rows = repmat(struct('sensor_id', NaN, 'template_xc_mm', NaN, ...
    'dynamic_xc_mm', NaN, 'delta_mm', NaN), numel(analysis_sensors), 1);
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    rows(is).sensor_id = sid;
    if ~isempty(Tpl) && isfield(Tpl, 'xc')
        rows(is).template_xc_mm = Tpl.xc;
    end
    rows(is).dynamic_xc_mm = infer_dynamic_xcenter_local(DynamicMap, sid);
    rows(is).delta_mm = rows(is).template_xc_mm - rows(is).dynamic_xc_mm;
end
delta = [rows.delta_mm];
check = struct();
check.table = struct2table(rows);
check.tolerance_mm = tolerance_mm;
check.max_abs_delta_mm = max(abs(delta), [], 'omitnan');
check.is_consistent = isfinite(check.max_abs_delta_mm) && check.max_abs_delta_mm <= tolerance_mm;
if check.is_consistent
    check.status = 'consistent';
else
    check.status = 'mismatch';
end
if isfield(DynamicMap, 'TemplateFileForXRel')
    check.dynamic_template_file_for_xrel = DynamicMap.TemplateFileForXRel;
end
end

function xc = infer_dynamic_xcenter_local(DynamicMap, sid)
xc = NaN;
for iw = 1:numel(DynamicMap.Window)
    S = DynamicMap.Window(iw).Sensor;
    idx = find([S.sensor_id] == sid, 1, 'first');
    if isempty(idx) || isempty(S(idx).x_abs) || isempty(S(idx).x_rel)
        continue;
    end
    v = S(idx).x_abs(:) - S(idx).x_rel(:);
    xc = median(v(isfinite(v)), 'omitnan');
    return;
end
end

function bundle = build_template_observation_bundle_local(Wmap, Template, analysis_sensors, method)
X = [];
T = [];
V = [];
S = [];
W = [];
Theta = [];
Y_obs = [];
F0_all = [];
Fx_all = [];
sensor_index_pt = [];
window_meta = struct('sensor_id', {}, 'valid_points', {}, 'shape_rmse', {});

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    D = Wmap.Sensor([Wmap.Sensor.sensor_id] == sid);

    v_raw = D.V(:);
    x_raw = D.x_rel(:);
    t_raw = D.t(:);
    theta_raw = D.theta(:);
    f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_raw, 'pchip', NaN);
    fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x_raw, 'pchip', NaN);

    if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
        threshold = Tpl.threshold;
    else
        threshold = method.default_sensor_threshold;
    end
    x_domain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
    if strcmpi(method.pulse_selection_mode, 'all')
        [mask_pulse, pulse_segment_count] = isolate_all_pulses_local(v_raw, threshold);
    else
        [mask_pulse, pulse_segment_count] = isolate_main_pulse_local(v_raw, threshold);
    end
    mask_effective = build_dynamic_effective_mask_local(x_raw, v_raw, t_raw, Tpl, threshold, method);
    mask_domain = x_raw >= Tpl.x_domain(1) + method.domain_margin_mm & ...
                  x_raw <= Tpl.x_domain(2) - method.domain_margin_mm;
    mask_finite = isfinite(f0) & isfinite(fx0) & isfinite(v_raw) & isfinite(theta_raw);
    mask_base = mask_pulse & mask_effective & mask_domain & mask_finite;
    query_guard_mm = resolve_query_guard_mm_local(Tpl, x_raw, v_raw, mask_base, method);
    mask_query_safe = x_raw >= x_domain(1) + query_guard_mm & ...
                      x_raw <= x_domain(2) - query_guard_mm;
    mask = mask_base;
    if strcmpi(method.domain_selection_mode, 'hard')
        safe_mask = mask & mask_query_safe;
        if nnz(safe_mask) >= 8
            mask = safe_mask;
        end
    end
    if nnz(mask) < 8
        mask = mask_effective & mask_domain & mask_finite;
        if strcmpi(method.domain_selection_mode, 'hard')
            safe_mask = mask & mask_query_safe;
            if nnz(safe_mask) >= 8
                mask = safe_mask;
            end
        end
    end
    if nnz(mask) < 8
        continue;
    end

    t_sel = t_raw(mask);
    x_sel = x_raw(mask);
    v_sel = v_raw(mask);
    theta_sel = theta_raw(mask);
    f0_sel = f0(mask);
    fx_sel = fx0(mask);
    y_obs_sel = x_sel - invert_template_voltage_local(Tpl, v_sel, x_sel);

    w_edge = build_edge_weight_local(t_sel, v_sel, method.weight_floor);
    if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
        w_template = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x_sel, 'linear', method.weight_floor);
    else
        w_template = ones(size(x_sel));
    end
    w_domain = build_domain_soft_weight_local(x_sel, x_domain, method.domain_soft_margin_mm, method.weight_floor);
    w_query = build_query_guard_soft_weight_local( ...
        x_sel, x_domain, query_guard_mm, method.domain_soft_margin_mm, method.weight_floor);
    w_gradient = build_template_gradient_weight_local(Tpl, x_sel, method.domain_selection_mode, method.weight_floor);
    w_total = max(method.weight_floor, w_edge .* w_template .* w_domain .* w_query .* w_gradient);
    if max(w_total) > 0
        w_total = max(method.weight_floor, w_total ./ max(w_total));
    end

    v_static = interp1(Tpl.x_grid(:), Tpl.v_grid(:), x_sel, 'pchip', NaN);
    shape_rmse = sqrt(mean((v_sel - v_static).^2, 'omitnan'));

    X = [X; x_sel]; %#ok<AGROW>
    T = [T; t_sel]; %#ok<AGROW>
    V = [V; v_sel]; %#ok<AGROW>
    S = [S; repmat(sid, numel(t_sel), 1)]; %#ok<AGROW>
    W = [W; w_total]; %#ok<AGROW>
    Theta = [Theta; theta_sel]; %#ok<AGROW>
    Y_obs = [Y_obs; y_obs_sel]; %#ok<AGROW>
    F0_all = [F0_all; f0_sel]; %#ok<AGROW>
    Fx_all = [Fx_all; fx_sel]; %#ok<AGROW>
    sensor_index_pt = [sensor_index_pt; repmat(is, numel(t_sel), 1)]; %#ok<AGROW>

    meta_idx = numel(window_meta) + 1;
    window_meta(meta_idx).sensor_id = sid;
    window_meta(meta_idx).valid_points = numel(t_sel);
    window_meta(meta_idx).pulse_segment_count = pulse_segment_count;
    window_meta(meta_idx).dynamic_effective_points = nnz(mask_effective);
    window_meta(meta_idx).query_safe_points = nnz(mask_query_safe(mask));
    window_meta(meta_idx).query_guard_mm = query_guard_mm;
    window_meta(meta_idx).domain_selection_mode = method.domain_selection_mode;
    window_meta(meta_idx).domain_soft_margin_mm = method.domain_soft_margin_mm;
    window_meta(meta_idx).template_domain = x_domain;
    window_meta(meta_idx).selected_x_range = [min(x_sel), max(x_sel)];
    window_meta(meta_idx).shape_rmse = shape_rmse;
end

if isempty(T)
    error('No valid template observation points were constructed for window %d.', Wmap.window_id);
end

interp_v = cell(numel(analysis_sensors), 1);
x_domain_by_sensor = NaN(numel(analysis_sensors), 2);
for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    Tpl = Template.Sensor([Template.Sensor.sensor_id] == sid);
    interp_v{is} = griddedInterpolant(Tpl.x_grid(:), Tpl.v_grid(:), 'pchip', 'none');
    x_domain_by_sensor(is, :) = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];
end

bundle = struct();
bundle.X = X;
bundle.T = T;
bundle.T_rel = T - min(T);
bundle.V = V;
bundle.S = S;
bundle.W = W;
bundle.Theta = Theta;
bundle.Y_obs = Y_obs;
bundle.F0 = F0_all;
bundle.Fx = Fx_all;
bundle.sensor_index = sensor_index_pt;
bundle.sensor_ids = analysis_sensors;
bundle.interp_v = interp_v;
bundle.x_domain_by_sensor = x_domain_by_sensor;
bundle.query_guard_mm = max([window_meta.query_guard_mm], [], 'omitnan');
bundle.overshoot_penalty_weight = method.overshoot_penalty_weight;
bundle.sensor_eta_limit_mm = method.sensor_eta_limit_mm;
bundle.sensor_eta_reg_weight_v_per_mm = method.sensor_eta_reg_weight_v_per_mm;
bundle.window_meta = window_meta;
bundle.valid_segment_count = sum([window_meta.pulse_segment_count]);
bundle.point_count = numel(T);
bundle.time_window = [min(T), max(T)];
bundle.rot_freq_mean_hz = Wmap.rot_freq_mean_hz;
bundle.rot_rpm_mean = Wmap.rot_rpm_mean;
end

function x_stat = invert_template_voltage_local(Tpl, v, x_ref)
x_grid = Tpl.x_grid(:);
v_grid = Tpl.v_grid(:);
x_stat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diff_v = v_grid - vv;
    crossing_x = [];
    exact_idx = find(abs(diff_v) <= 1e-10);
    if ~isempty(exact_idx)
        crossing_x = x_grid(exact_idx);
    end
    for k = 1:numel(diff_v)-1
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k+1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k+1) > 0
            continue;
        end
        denom = v_grid(k+1) - v_grid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - v_grid(k)) / denom;
        crossing_x(end+1, 1) = x_grid(k) + alpha * (x_grid(k+1) - x_grid(k)); %#ok<AGROW>
    end
    if isempty(crossing_x)
        [~, idx] = min(abs(diff_v));
        x_stat(i) = x_grid(idx);
    else
        [~, idx] = min(abs(crossing_x - x_ref(i)));
        x_stat(i) = crossing_x(idx);
    end
end
end

function query_guard_mm = resolve_query_guard_mm_local(Tpl, x, v, base_mask, method)
query_guard_mm = method.query_guard_mm;
if ~strcmpi(method.query_guard_mode, 'adaptive') || nnz(base_mask) < 8
    return;
end
x_sel = x(base_mask);
v_sel = v(base_mask);
x_stat = invert_template_voltage_local(Tpl, v_sel, x_sel);
u_app = abs(x_sel(:) - x_stat(:));
u_app = u_app(isfinite(u_app));
if isempty(u_app)
    return;
end
q = min(max(method.query_guard_quantile, 0), 100);
adaptive_guard = prctile(u_app, q) + method.query_guard_safety_mm;
query_guard_mm = min(max(adaptive_guard, method.query_guard_min_mm), method.query_guard_max_mm);
end

function eo_candidates = build_eo_candidates_local(rot_freq_hz, freq_search_hz, eo_pad)
if nargin < 3 || isempty(eo_pad)
    eo_pad = 0;
end
rot_freq_hz = max(rot_freq_hz, eps);
freq_lo = min(freq_search_hz);
freq_hi = max(freq_search_hz);
eo_min_strict = max(1, ceil(freq_lo / rot_freq_hz));
eo_max_strict = max(eo_min_strict, floor(freq_hi / rot_freq_hz));
eo_candidates = eo_min_strict:eo_max_strict;
if isempty(eo_candidates)
    eo_center = max(1, round(mean([freq_lo, freq_hi]) / rot_freq_hz));
    eo_candidates = max(1, eo_center - eo_pad):max(1, eo_center + eo_pad);
end
freq_candidates = eo_candidates .* rot_freq_hz;
eo_candidates = eo_candidates(freq_candidates >= freq_lo & freq_candidates <= freq_hi);
if isempty(eo_candidates)
    error('No EO candidates remain inside [%.3f, %.3f] Hz at rot_freq %.6f Hz.', ...
        freq_lo, freq_hi, rot_freq_hz);
end
end

function seed_table = solve_template_seed_eo_scan_local(bundle, eo_candidates)
eo_candidates = unique(round(eo_candidates(:).'));
rows = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'weighted_voltage_rmse', inf, 'plain_voltage_rmse', inf, ...
    'linear_vp_rmse', inf), numel(eo_candidates), 1);
n_sensor = numel(bundle.sensor_ids);
for i = 1:numel(eo_candidates)
    eo = eo_candidates(i);
    s1 = sin(eo * bundle.Theta);
    c1 = cos(eo * bundle.Theta);
    basis = [-bundle.Fx(:), -bundle.Fx(:) .* s1, -bundle.Fx(:) .* c1];
    w_sqrt = sqrt(bundle.W(:));
    y = bundle.V(:) - bundle.F0(:);
    coeff = (basis .* w_sqrt) \ (y .* w_sqrt);
    dx_c = coeff(1);
    eta = zeros(1, n_sensor);
    A = hypot(coeff(2), coeff(3));
    phi = atan2(coeff(3), coeff(2));
    lin_res = y - basis * coeff;
    [obj, plain_rmse] = template_synchronous_objective_local([A, phi, dx_c, eta], eo, bundle);
    rows(i).EO = eo;
    rows(i).A = A;
    rows(i).phi = phi;
    rows(i).dx_c = dx_c;
    rows(i).sensor_eta = eta;
    rows(i).sensor_eta_max_abs = max(abs(eta), [], 'omitnan');
    rows(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    rows(i).plain_voltage_rmse = plain_rmse;
    rows(i).linear_vp_rmse = sqrt(sum(bundle.W(:) .* lin_res.^2) / max(sum(bundle.W(:)), eps));
end
score = [[rows.weighted_voltage_rmse].', [rows.linear_vp_rmse].', [rows.plain_voltage_rmse].', [rows.EO].'];
[~, order] = sortrows(score, [1, 2, 3]);
seed_table = rows(order);
end

function eo_keep = select_vp_topk_eo_local(seed_table, top_k)
if nargin < 2 || isempty(top_k) || ~isfinite(top_k)
    eo_keep = unique(round([seed_table.EO]), 'stable');
    return;
end
top_k = min(max(1, floor(top_k)), numel(seed_table));
eo_keep = unique(round([seed_table(1:top_k).EO]), 'stable');
end

function [eo_keep, info] = select_vp_adaptive_eo_local(seed_table, eo_candidates, method, cfg, window_id)
all_eo = unique(round(eo_candidates(:).'));
top_k_eo = select_vp_topk_eo_local(seed_table, method.vp_top_k_eo);
eo_keep = intersect(all_eo, top_k_eo, 'stable');

best_rmse = seed_table(1).weighted_voltage_rmse;
if numel(seed_table) >= 2
    second_rmse = seed_table(2).weighted_voltage_rmse;
else
    second_rmse = inf;
end
vp_gap_ratio = (second_rmse - best_rmse) / max(best_rmse, eps);
linear_gap_ratio = inf;
if numel(seed_table) >= 2
    linear_gap_ratio = (seed_table(2).linear_vp_rmse - seed_table(1).linear_vp_rmse) / ...
        max(seed_table(1).linear_vp_rmse, eps);
end

fallback_reason = '';
if window_id <= cfg.vp_all_eo_warmup_windows
    fallback_reason = 'warmup_all_eo';
elseif vp_gap_ratio < cfg.vp_gap_ratio_fallback
    fallback_reason = 'ambiguous_weighted_vp_gap';
elseif linear_gap_ratio < cfg.vp_linear_gap_ratio_fallback
    fallback_reason = 'ambiguous_linear_vp_gap';
end

if ~isempty(fallback_reason)
    eo_keep = all_eo;
end

info = struct();
info.window_id = window_id;
info.all_eo = all_eo;
info.top_k_eo = top_k_eo;
info.selected_eo = eo_keep;
info.best_seed_eo = seed_table(1).EO;
info.best_seed_weighted_rmse = best_rmse;
info.second_seed_weighted_rmse = second_rmse;
info.vp_gap_ratio = vp_gap_ratio;
info.linear_gap_ratio = linear_gap_ratio;
info.fallback_reason = fallback_reason;
end

function result = refine_template_waveform_fit_local(bundle, seed_table, eo_candidates, cfg)
amp_limit_mm = abs(cfg.amplitude_limit_mm);
dx_c_limit_mm = abs(cfg.dx_c_limit_mm);
sensor_eta_limit_mm = abs(cfg.sensor_eta_limit_mm);
eo_candidates = unique(round(eo_candidates(:).'));
n_sensor = numel(bundle.sensor_ids);
candidate = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx_c', NaN, ...
    'sensor_eta', [], 'sensor_eta_max_abs', NaN, ...
    'fn_hz', NaN, 'objective', inf, 'weighted_voltage_rmse', inf, ...
    'plain_voltage_rmse', inf, 'V_pred', [], 'u_est', []), numel(eo_candidates), 1);
best = struct('objective', inf);

for i = 1:numel(eo_candidates)
    eo = eo_candidates(i);
    seed_idx = find([seed_table.EO] == eo, 1, 'first');
    if isempty(seed_idx)
        seed = seed_table(1);
    else
        seed = seed_table(seed_idx);
    end
    eta0 = zeros(1, n_sensor);
    if isfield(seed, 'sensor_eta') && ~isempty(seed.sensor_eta)
        eta0(1:min(n_sensor, numel(seed.sensor_eta))) = seed.sensor_eta(1:min(n_sensor, numel(seed.sensor_eta)));
    end
    seed_params = [seed.A, seed.phi, seed.dx_c, eta0];
    fun = @(p) template_synchronous_objective_local(p, eo, bundle);
    opts = optimset('Display', 'off', 'MaxIter', 400, 'MaxFunEvals', 1200);
    p_opt = fminsearch(@(p) bounded_synchronous_objective_local(p, fun, amp_limit_mm, dx_c_limit_mm, sensor_eta_limit_mm, n_sensor), seed_params, opts);
    p_opt(1) = min(abs(p_opt(1)), amp_limit_mm);
    p_opt(2) = wrap_to_pi_local(p_opt(2));
    p_opt(3) = max(min(p_opt(3), dx_c_limit_mm), -dx_c_limit_mm);
    p_opt = bound_sensor_eta_params_local(p_opt, sensor_eta_limit_mm, n_sensor);
    [obj, plain_rmse, v_pred, u_est, coverage] = template_synchronous_objective_local(p_opt, eo, bundle);

    candidate(i).EO = eo;
    candidate(i).A = p_opt(1);
    candidate(i).phi = p_opt(2);
    candidate(i).dx_c = p_opt(3);
    candidate(i).sensor_eta = p_opt(4:3+n_sensor);
    candidate(i).sensor_eta_max_abs = max(abs(candidate(i).sensor_eta), [], 'omitnan');
    candidate(i).fn_hz = eo * bundle.rot_freq_mean_hz;
    candidate(i).objective = obj;
    candidate(i).weighted_voltage_rmse = sqrt(obj / bundle.point_count);
    candidate(i).plain_voltage_rmse = plain_rmse;
    candidate(i).clamp_fraction = coverage.clamp_fraction;
    candidate(i).max_query_overshoot_mm = coverage.max_query_overshoot_mm;
    candidate(i).overshoot_penalty = coverage.overshoot_penalty;
    candidate(i).V_pred = v_pred;
    candidate(i).u_est = u_est;

    if obj < best.objective
        best = candidate(i);
    end
end

candidate_table = struct2table(rmfield(candidate, {'V_pred','u_est','sensor_eta'}));
candidate_table = sortrows(candidate_table, {'weighted_voltage_rmse','plain_voltage_rmse','EO'}, ...
    {'ascend','ascend','ascend'});

result = struct();
result.A_id = best.A;
result.EO_id = best.EO;
result.phi_id_wrapped = best.phi;
result.dx_c_id = best.dx_c;
result.d0_id = best.dx_c;
result.sensor_eta_id = best.sensor_eta;
result.sensor_eta_max_abs_mm = best.sensor_eta_max_abs;
result.fn_id = best.fn_hz;
result.weighted_voltage_rmse = best.weighted_voltage_rmse;
result.plain_voltage_rmse = best.plain_voltage_rmse;
result.V_pred = best.V_pred;
result.u_est = best.u_est;
result.objective_score = best.objective;
result.valid_segment_count = bundle.valid_segment_count;
result.point_count = bundle.point_count;
result.CandidateTable = candidate_table;
result.Coverage = struct( ...
    'clamp_fraction', best.clamp_fraction, ...
    'max_query_overshoot_mm', best.max_query_overshoot_mm, ...
    'query_guard_mm', bundle.query_guard_mm, ...
    'overshoot_penalty', best.overshoot_penalty, ...
    'overshoot_penalty_weight', bundle.overshoot_penalty_weight);
result.metric_note = 'Template-only VP top-K synchronous waveform fit with shared dx_c plus per-sensor x-zero eta; no EO prior.';
end

function obj = bounded_synchronous_objective_local(p, fun, amp_limit_mm, dx_c_limit_mm, sensor_eta_limit_mm, n_sensor)
p_use = p;
p_use(1) = min(abs(p_use(1)), amp_limit_mm);
p_use(2) = wrap_to_pi_local(p_use(2));
p_use(3) = max(min(p_use(3), dx_c_limit_mm), -dx_c_limit_mm);
p_use = bound_sensor_eta_params_local(p_use, sensor_eta_limit_mm, n_sensor);
obj = fun(p_use);
end

function p = bound_sensor_eta_params_local(p, sensor_eta_limit_mm, n_sensor)
need = 3 + n_sensor;
if numel(p) < need
    p(numel(p)+1:need) = 0;
end
p = p(1:need);
eta = p(4:end);
eta = max(min(eta, sensor_eta_limit_mm), -sensor_eta_limit_mm);
if ~isempty(eta)
    eta = eta - eta(1);
    eta = max(min(eta, sensor_eta_limit_mm), -sensor_eta_limit_mm);
end
p(4:end) = eta;
end

function [obj, plain_rmse, v_pred, u_est, coverage] = template_synchronous_objective_local(params, eo, bundle)
A = params(1);
phi = params(2);
dx_c = params(3);
params = bound_sensor_eta_params_local(params, abs(bundle.sensor_eta_limit_mm), numel(bundle.sensor_ids));
eta = params(4:end).';
u_est = A .* sin(eo .* bundle.Theta + phi);
eta_sample = eta(bundle.sensor_index(:));
x_in = bundle.X - dx_c - eta_sample - u_est;
v_pred = NaN(size(bundle.V));
clamped = false(size(bundle.V));
overshoot = zeros(size(bundle.V));
for is = 1:numel(bundle.sensor_ids)
    mask = bundle.sensor_index == is;
    x_lo = bundle.x_domain_by_sensor(is, 1);
    x_hi = bundle.x_domain_by_sensor(is, 2);
    x_eval_raw = x_in(mask);
    x_eval = min(max(x_eval_raw, x_lo), x_hi);
    clamped(mask) = x_eval_raw < x_lo | x_eval_raw > x_hi;
    overshoot(mask) = max(x_lo - x_eval_raw, 0) + max(x_eval_raw - x_hi, 0);
    v_pred(mask) = bundle.interp_v{is}(x_eval);
end
valid = isfinite(v_pred);
res = bundle.V(valid) - v_pred(valid);
overshoot_penalty = bundle.overshoot_penalty_weight * ...
    sum(bundle.W(valid) .* (overshoot(valid) .^ 2));
eta_reg_penalty = bundle.point_count * (bundle.sensor_eta_reg_weight_v_per_mm * sqrt(mean(eta .^ 2))) ^ 2;
obj = sum(bundle.W(valid) .* (res .^ 2)) + overshoot_penalty + eta_reg_penalty + 1e6 * nnz(~valid);
plain_rmse = sqrt(mean(res .^2));
coverage = struct();
coverage.clamp_fraction = nnz(clamped) / max(numel(clamped), 1);
coverage.max_query_overshoot_mm = max(overshoot, [], 'omitnan');
coverage.overshoot_penalty = overshoot_penalty;
coverage.sensor_eta_reg_penalty = eta_reg_penalty;
end

function summary = build_resonance_summary_local(trend_rows)
eo_vals = [trend_rows.EO_id];
freq_vals = [trend_rows.fn_id];
amp_vals = [trend_rows.A_id];
wrmse_vals = [trend_rows.weighted_voltage_rmse];
summary = struct();
summary.dominant_eo = mode(eo_vals(isfinite(eo_vals)));
summary.mean_freq_hz = mean(freq_vals, 'omitnan');
summary.median_freq_hz = median(freq_vals, 'omitnan');
summary.mean_amp_mm = mean(amp_vals, 'omitnan');
summary.median_amp_mm = median(amp_vals, 'omitnan');
summary.mean_weighted_rmse_v = mean(wrmse_vals, 'omitnan');
summary.median_weighted_rmse_v = median(wrmse_vals, 'omitnan');
end

function [mask, segment_count] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segment_count = numel(starts);
end

function [mask, segment_count] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segment_count = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segment_count = 1;
    return;
end

segment_gap = idx(starts(2:end)) - idx(ends(1:end-1));
large_gap_threshold = max(10, round(0.02 * numel(v)));
pulse_breaks = find(segment_gap > large_gap_threshold);
group_starts = [1; pulse_breaks + 1];
group_ends = [pulse_breaks; numel(starts)];
for ig = 1:numel(group_starts)
    seg_ids = group_starts(ig):group_ends(ig);
    best_seg = seg_ids(1);
    best_peak = -inf;
    for iseg = seg_ids
        seg = idx(starts(iseg):ends(iseg));
        peak_val = max(v(seg));
        if peak_val > best_peak
            best_peak = peak_val;
            best_seg = iseg;
        end
    end
    mask(idx(starts(best_seg)):idx(ends(best_seg))) = true;
end
segment_count = numel(group_starts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamic_effective_mode, 'gradient')
    mask = finite;
    return;
end

g_tpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    g_tpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
g_tpl_max = max(g_tpl(finite), [], 'omitnan');
if ~isfinite(g_tpl_max) || g_tpl_max <= 0
    mask_tpl = finite;
else
    mask_tpl = g_tpl >= method.dynamic_template_gradient_min_ratio * g_tpl_max;
end

g_time = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    g_time(finite) = abs(gradient(v(finite), t(finite)));
end
g_time_max = max(g_time(finite), [], 'omitnan');
if ~isfinite(g_time_max) || g_time_max <= 0
    mask_time = false(size(v));
else
    mask_time = g_time >= method.dynamic_time_gradient_min_ratio * g_time_max;
end

peak_q = min(max(method.dynamic_peak_quantile, 0), 100);
peak_level = prctile(v(finite), peak_q);
mask_peak = v >= max(threshold, peak_level);

mask = finite & (mask_tpl | mask_time | mask_peak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function w_edge = build_edge_weight_local(t, v, floor_w)
if numel(v) < 3 || range(t) <= 0
    w_edge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    w_edge = dv ./ max(dv);
else
    w_edge = ones(size(dv));
end
w_edge = max(floor_w, w_edge);
end

function w_domain = build_domain_soft_weight_local(x, x_domain, margin_mm, floor_w)
if margin_mm <= 0
    w_domain = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
ratio = min(max(dist_to_edge ./ margin_mm, 0), 1);
w_domain = floor_w + (1 - floor_w) .* ratio;
w_domain(~isfinite(w_domain)) = floor_w;
end

function w_query = build_query_guard_soft_weight_local(x, x_domain, query_guard_mm, margin_mm, floor_w)
if query_guard_mm <= 0 || margin_mm <= 0
    w_query = ones(size(x));
    return;
end
dist_to_edge = min(x(:) - x_domain(1), x_domain(2) - x(:));
soft_start = max(query_guard_mm - margin_mm, 0);
ratio = min(max((dist_to_edge - soft_start) ./ max(margin_mm, eps), 0), 1);
w_query = floor_w + (1 - floor_w) .* ratio;
w_query(~isfinite(w_query)) = floor_w;
end

function w_gradient = build_template_gradient_weight_local(Tpl, x, domain_selection_mode, floor_w)
if ~strcmpi(domain_selection_mode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    w_gradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
g_max = max(g, [], 'omitnan');
if ~isfinite(g_max) || g_max <= 0
    w_gradient = ones(size(x));
    return;
end
w_gradient = floor_w + (1 - floor_w) .* g ./ g_max;
w_gradient(~isfinite(w_gradient)) = floor_w;
end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end

function s = rank_to_string_local(rank_value)
if isfinite(rank_value)
    s = sprintf('%d', rank_value);
else
    s = 'NA';
end
end
end

