%% Step01_Main_Build_OPRCenterStd_Template_20250527
% Build OPR-center timing and the latest OPRCenterStd low-speed template.
% Tune parameters here, then run this file directly.

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
C = CaseConfig();

%% Parameters to tune
P.name.referenceTemplateSuffix = 'OPRCenterStdRef';
P.name.templateSuffix = 'GradientXRange030_OPRCenterStd';

P.step01.showLegacyPlots = true;
P.step01.forceLegacyRebuild = false;
P.step01.forceOprCenterRebuild = false;
P.step01.centerMode = 'sgfit';
P.step01.xrangeMode = 'threshold';
P.step01.gradientMinRatio = 0.30;
P.step01.amplitudeMinRatio = 0;
P.step01.minHalfWidthMM = 2.5;
P.step01.maxHalfWidthMM = 4.2;
P.step01.etaLimitMM = 0.20;
P.step01.etaAdaptiveQuantile = 85;
P.step01.etaAdaptiveSafetyFactor = 1.20;
P.step01.etaAdaptiveMinMM = 0.02;
P.step01.etaAdaptiveMaxMM = P.step01.etaLimitMM;
P.step01.referenceCenterMode = 'sgfit';
P.step01.referenceXrangeMode = 'energy';
P.step01.rebuildReferenceTemplate = false;
allowReferenceOnlyRaw = lower(strtrim(getenv('STEP01_ALLOW_REFERENCE_ONLY')));
P.step01.allowReferenceOnly = any(strcmp(allowReferenceOnlyRaw, {'1','true','yes','on'}));
P.view.saveFigures = false;

%% Files
stablePlanFile = fullfile(routeDir, 'output', 'stable_window_plan', ...
    sprintf('Step01_CoverageFirstStableWindowPlan_%s_%s.csv', C.caseTag, C.dataset));
referenceCenterTemplateFile = fullfile(routeDir, 'output', 'templates', ...
    sprintf('Template_LowSpeedRotating_%s_%s_%s.mat', C.caseTag, P.name.referenceTemplateSuffix, C.dataset));

fprintf('\n=== Step01: OPR-center timing and low-speed template ===\n');
apply_step01_preprocess_settings(P);
step01_extract_opr_blade_timing_embedded(routeDir, P);

if P.step01.allowReferenceOnly
    fprintf('Reference-only bootstrap mode enabled.\n');
    apply_step01_reference_template_settings(P);
    step01_build_rotating_template_embedded(routeDir, C);
    fprintf('Reference-only Step01 complete: %s\n', referenceCenterTemplateFile);
    return;
end

if exist(stablePlanFile, 'file') ~= 2
    error(['Stable-window plan not found:\n  %s\n' ...
        'The old planning script is archived under archive_legacy_flows/old_programs.'], ...
        stablePlanFile);
end

if P.step01.rebuildReferenceTemplate || exist(referenceCenterTemplateFile, 'file') ~= 2
    fprintf('Building reference-center template first...\n');
    apply_step01_reference_template_settings(P);
    step01_build_rotating_template_embedded(routeDir, C);
end

apply_step01_main_template_settings(P, stablePlanFile, referenceCenterTemplateFile);
step01_build_rotating_template_embedded(routeDir, C);

function apply_step01_preprocess_settings(P)
setenv('STEP01_FORCE_LEGACY_REBUILD', logical_to_env_local(P.step01.forceLegacyRebuild));
setenv('STEP01_FORCE_OPR_CENTER_REBUILD', logical_to_env_local(P.step01.forceOprCenterRebuild));
setenv('STEP01_SAVE_FIGURES', logical_to_env_local(P.view.saveFigures));
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
setenv('STEP01_ETA_LIMIT_MM', num2str(P.step01.etaLimitMM));
setenv('STEP01_ETA_ADAPTIVE_QUANTILE', num2str(P.step01.etaAdaptiveQuantile));
setenv('STEP01_ETA_ADAPTIVE_SAFETY_FACTOR', num2str(P.step01.etaAdaptiveSafetyFactor));
setenv('STEP01_ETA_ADAPTIVE_MIN_MM', num2str(P.step01.etaAdaptiveMinMM));
setenv('STEP01_ETA_ADAPTIVE_MAX_MM', num2str(P.step01.etaAdaptiveMaxMM));
end

function s = logical_to_env_local(tf)
if tf
    s = '1';
else
    s = '0';
end
end

function step01_extract_opr_blade_timing_embedded(rootDir, P)
%% Step01: Extract blade timing and upgrade OPR timing to pulse centers
% This is the fixed preprocessing entry for the direct low-speed-template
% route. It first ensures the legacy Step2 products exist, then rewrites
% jiluOPR/omega so jiluOPR(:,1) is the multi-threshold OPR pulse center.

routeDir = rootDir;
validationRoot = fileparts(routeDir);
mainLegacyDir = fullfile(validationRoot, '20250527_low_speed_gap_prior_decoupling', 'legacy');
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
if P.view.saveFigures
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    saveas(fig, fullfile(outDir, 'Step01_OPRCenter_Blade_Timing_20250527.png'));
end

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

function step01_build_rotating_template_embedded(rootDir, C)
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
legacy_dir = fullfile(validation_root, '20250527_low_speed_gap_prior_decoupling', 'legacy');
if exist(legacy_dir, 'dir') ~= 7
    error('Legacy 20250527 helper folder not found: %s', legacy_dir);
end
addpath(legacy_dir);

legacy_cfg = Get_20250527_BTT_Config();
low_speed_case = legacy_cfg.low_speed_case;
low_speed_data_dir = fullfile(legacy_cfg.dataset_root, low_speed_case);
sensor_config_file = fullfile(legacy_cfg.reference_output_dir, 'Sensor_Config_20250527.mat');

target_blade = C.bladeId;
analysis_sensors = C.sensorIds(:).';       % build a reusable sensor library; Step03 can use subsets
sensor_tag = C.sensorTag;

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

template_file = fullfile(template_dir, sprintf('Template_LowSpeedRotating_%s%s_%s.mat', ...
    C.caseTag, template_suffix, C.dataset));

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
    'eta_limit_mm', parse_nonnegative_numeric_env_local('STEP01_ETA_LIMIT_MM', 0.20), ...
    'eta_adaptive_quantile', parse_nonnegative_numeric_env_local('STEP01_ETA_ADAPTIVE_QUANTILE', 85), ...
    'eta_adaptive_safety_factor', parse_nonnegative_numeric_env_local('STEP01_ETA_ADAPTIVE_SAFETY_FACTOR', 1.20), ...
    'eta_adaptive_min_mm', parse_nonnegative_numeric_env_local('STEP01_ETA_ADAPTIVE_MIN_MM', 0.02), ...
    'eta_adaptive_max_mm', parse_nonnegative_numeric_env_local('STEP01_ETA_ADAPTIVE_MAX_MM', 0.20), ...
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
    'spline_smoothing', NaN, ...
    'eta_app_mm', [], ...
    'eta_median_mm', NaN, ...
    'eta_iqr_mm', NaN, ...
    'eta_limit_mm', NaN), numel(analysis_sensors), 1);

plot_cache = repmat(struct('sensor_id', NaN, 'x_wide', [], 'v_wide', [], ...
    'x_selected', [], 'v_selected', [], 'weight_selected', [], ...
    'xc_detected_rel', NaN, 'xc_reference_source', ''), numel(analysis_sensors), 1);

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    threshold = legacy_cfg.sensor_thresholds(sid);
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, target_blade, opr_reference);
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
    [eta_app_mm, eta_median_mm, eta_iqr_mm, eta_limit_mm] = estimate_low_speed_eta_scatter_local( ...
        x_rel, v_raw, x_grid, v_grid, Template.Settings);

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
    Template.Sensor(is).eta_app_mm = eta_app_mm(:);
    Template.Sensor(is).eta_median_mm = eta_median_mm;
    Template.Sensor(is).eta_iqr_mm = eta_iqr_mm;
    Template.Sensor(is).eta_limit_mm = eta_limit_mm;

    plot_cache(is).sensor_id = sid;
    plot_cache(is).x_wide = x_wide(:) - xc_for_coordinate;
    plot_cache(is).v_wide = v_wide(:);
    plot_cache(is).x_selected = P.x_selected(:) - xc_for_coordinate;
    plot_cache(is).v_selected = P.v_selected(:);
    plot_cache(is).weight_selected = P.weight_selected(:);
    plot_cache(is).xc_detected_rel = P.xc - xc_for_coordinate;
    plot_cache(is).xc_reference_source = xc_reference_source;

    fprintf(['CH%d: baseline %.5f V, threshold %.3f V, pulses %d, ' ...
        'selected points %d/%d, domain [%.3f, %.3f] mm, xc %.4f mm, etaLim %.4f mm, ' ...
        'detected xc %.4f mm, source %s\n'], ...
        sid, baseline, threshold, numel(selected_idx), numel(x_rel), numel(x_wide), ...
        x_left, x_right, xc_for_coordinate, eta_limit_mm, P.xc, xc_reference_source);
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

if any(strcmpi(strtrim(getenv('STEP01_SAVE_FIGURES')), {'1','true','yes','on'}))
    exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_LowSpeedSplineTemplates_B%d_%s%s.png', ...
        target_blade, sensor_tag, template_suffix)), 'Resolution', 300);
    exportgraphics(fig, fullfile(figure_dir, sprintf('Step01_LowSpeedSplineTemplates_B%d_%s%s.pdf', ...
        target_blade, sensor_tag, template_suffix)), 'ContentType', 'vector');
end

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

function theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference)
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    theta_std = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, blade_id);
    return;
end
theta_std = Sensor_Config.Standard_Relative_Angles(sid, blade_id);
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
theta_std = convert_standard_angle_to_opr_center_local(theta_std, opr_reference);
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

function [eta_app_mm, eta_median_mm, eta_iqr_mm, eta_limit_mm] = estimate_low_speed_eta_scatter_local( ...
        x_sel, v_sel, x_grid, v_grid, settings)
eta_app_mm = [];
eta_median_mm = NaN;
eta_iqr_mm = NaN;
eta_limit_mm = NaN;
if isempty(x_sel) || isempty(v_sel) || isempty(x_grid) || isempty(v_grid)
    return;
end
keep = isfinite(x_sel(:)) & isfinite(v_sel(:));
if nnz(keep) < 5
    return;
end
x_use = x_sel(keep);
v_use = v_sel(keep);
Tpl = struct('x_grid', x_grid(:), 'v_grid', v_grid(:));
x_stat = invert_template_voltage_local(Tpl, v_use, x_use);
eta_app_mm = x_use - x_stat;
eta_app_mm = eta_app_mm(isfinite(eta_app_mm));
if isempty(eta_app_mm)
    return;
end
eta_median_mm = median(eta_app_mm, 'omitnan');
eta_iqr_mm = iqr(eta_app_mm);
if ~isfield(settings, 'eta_limit_mm') || ~isfinite(settings.eta_limit_mm) || settings.eta_limit_mm <= 0
    return;
end
q = 85;
if isfield(settings, 'eta_adaptive_quantile') && isfinite(settings.eta_adaptive_quantile)
    q = min(max(settings.eta_adaptive_quantile, 0), 100);
end
safety = 1.2;
if isfield(settings, 'eta_adaptive_safety_factor') && isfinite(settings.eta_adaptive_safety_factor)
    safety = max(settings.eta_adaptive_safety_factor, 1);
end
min_mm = 0.02;
if isfield(settings, 'eta_adaptive_min_mm') && isfinite(settings.eta_adaptive_min_mm)
    min_mm = max(settings.eta_adaptive_min_mm, 0);
end
max_mm = abs(settings.eta_limit_mm);
if isfield(settings, 'eta_adaptive_max_mm') && isfinite(settings.eta_adaptive_max_mm)
    max_mm = max(abs(settings.eta_limit_mm), abs(settings.eta_adaptive_max_mm));
end
eta_cand = safety * prctile(abs(eta_app_mm(:)), q);
if ~isfinite(eta_cand) || eta_cand <= 0
    eta_cand = safety * abs(eta_iqr_mm);
end
eta_limit_mm = min(max(eta_cand, min_mm), max_mm);
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
    for k = 1:(numel(diff_v) - 1)
        if ~isfinite(diff_v(k)) || ~isfinite(diff_v(k + 1))
            continue;
        end
        if diff_v(k) == 0 || diff_v(k) * diff_v(k + 1) > 0
            continue;
        end
        denom = v_grid(k + 1) - v_grid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - v_grid(k)) / denom;
        crossing_x(end + 1, 1) = x_grid(k) + alpha * (x_grid(k + 1) - x_grid(k)); %#ok<AGROW>
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
end
