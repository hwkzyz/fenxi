%% Step06: Build dynamic highMap for waveform-domain decoupling
% This step reuses the verified legacy waveform-window extraction, but
% exports the data in the voltage-domain structure required by the new
% clearance-vibration decoupling method.
%
% The highMap coordinate x_v is the dynamic pulse-local circumferential
% displacement in millimetres, mapped onto the same coordinate unit as the
% Step05 static response surface.  Step05 obtains its coordinate from the
% average calibration tip speed; Step06 uses the dynamic OPR phase-restored
% x_points and aligns them with the low-speed template center determined
% from OPR and sensor installation geometry. Dynamic pulse peaks are retained
% only for diagnostics and are not used as the fitting-coordinate origin.

if ~exist('STEP06_SKIP_CLEAR', 'var') || ~STEP06_SKIP_CLEAR
    clear; clc; close all;
else
    clc; close all;
end
forceRebuild = strcmpi(strtrim(getenv('STEP06_FORCE_REBUILD')), '1');

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
C0 = CaseConfig();
legacyDir = fullfile(thisDir, 'legacy');
addpath(legacyDir);

outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

responseFile = fullfile(outDir, 'Step05_Response_Surface_20251222.mat');
speedSummaryFile = fullfile(outDir, 'Step05A_single_gap_speed_summary', ...
    'Step05A_Average_Calibration_Speed.csv');
regionFile = fullfile(outDir, 'Step04_BTT_STE_Resonance_Regions_20251222.mat');
if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
if isfile(speedSummaryFile)
    speedSummary = readtable(speedSummaryFile);
    referenceTipSpeedMmS = speedSummary.averageTangentialSpeedMmS(1);
    referenceRpm = speedSummary.averageRpmCentroid(1);
    referenceSpeedSource = speedSummaryFile;
elseif isfield(responseSurface, 'referenceTangentialSpeedMmS') && ...
        isfield(responseSurface, 'referenceRpm')
    referenceTipSpeedMmS = responseSurface.referenceTangentialSpeedMmS;
    referenceRpm = responseSurface.referenceRpm;
    referenceSpeedSource = responseSurface.referenceTangentialSpeedSource;
else
    error('Missing calibration reference speed in Step05A or responseSurface.');
end
if isfile(regionFile)
    R = load(regionFile, 'regionTable');
    regionTable = R.regionTable;
else
    regionTable = table();
end

cfgLegacy = Get_20251222_BTT_Config();
Pflow = C0.flowConfig;
caseName = Pflow.identification.dynamicCaseName;
targetBlade = C0.bladeId;
analysisSensors = C0.sensorIds;
analysisStartMode = Pflow.identification.analysisStartMode;
manualAnalysisStartTime = Pflow.identification.analysisStartTimeSec;
targetLaps = Pflow.identification.targetBladePasses;
analysisWinSize = Pflow.identification.windowBladePasses;
slidingStep = Pflow.identification.slidingStepBladePasses;
referenceOrder = Pflow.identification.referenceOrder;
pulseWindowSec = Pflow.identification.pulseWindowSec;
resonanceRegionId = Pflow.identification.resonanceRegionId;
resonanceRegionTag = Pflow.identification.resonanceRegionTag;
resonanceRegionShortTag = Pflow.identification.resonanceRegionShortTag;
highMapCacheKey = make_highmap_cache_key_local(resonanceRegionShortTag, ...
    manualAnalysisStartTime, targetLaps, analysisWinSize, slidingStep, pulseWindowSec);
cacheMatFile = fullfile(outDir, sprintf('Step06_HighMap_%s_%s_%s.mat', ...
    C0.dataset, C0.caseTag, highMapCacheKey));
dynamicMapCacheFile = fullfile(outDir, sprintf('Step06_DynamicMap_%s_%s_%s.mat', ...
    C0.dataset, C0.caseTag, highMapCacheKey));
legacyMatFile = fullfile(outDir, sprintf('Step06_HighMap_%s_%s_%s.mat', ...
    C0.dataset, C0.caseTag, resonanceRegionShortTag));
cacheCsvFile = fullfile(outDir, sprintf('Step06_HighMap_Window_Summary_%s_%s_%s.csv', ...
    C0.dataset, C0.caseTag, highMapCacheKey));
legacyCsvFile = fullfile(outDir, sprintf('Step06_HighMap_Window_Summary_%s_%s_%s.csv', ...
    C0.dataset, C0.caseTag, resonanceRegionShortTag));
cacheBaselineCsvFile = fullfile(outDir, sprintf('Step06_Dynamic_Continuous_Baseline_%s_%s_%s.csv', ...
    C0.dataset, C0.caseTag, highMapCacheKey));
legacyBaselineCsvFile = fullfile(outDir, sprintf('Step06_Dynamic_Continuous_Baseline_%s_%s_%s.csv', ...
    C0.dataset, C0.caseTag, resonanceRegionShortTag));

if isfile(cacheMatFile) && ~forceRebuild
    load(cacheMatFile, 'highMap');
    if is_highmap_cache_match_local(highMap, C0, resonanceRegionShortTag, ...
            manualAnalysisStartTime, targetLaps, analysisWinSize, slidingStep, pulseWindowSec)
        fprintf('\nStep06 highMap cache hit:\n  %s\n', cacheMatFile);
        if ~isfile(dynamicMapCacheFile)
            if isfield(highMap, 'lowSpeedTemplateFile') && isfile(highMap.lowSpeedTemplateFile)
                templateFile = highMap.lowSpeedTemplateFile;
            else
                rotDirForTemplate = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
                addpath(rotDirForTemplate);
                templateCase = struct();
                templateCase.dataset = '20251222';
                templateCase.bladeId = targetBlade;
                templateCase.sensorIds = analysisSensors;
                templateCase.sensorTag = ['S', sprintf('%d', analysisSensors)];
                templateCase.caseTag = sprintf('B%d_%s', targetBlade, templateCase.sensorTag);
                [templateFile, ~] = resolveTemplateFile_OPRCenterStd_20251222( ...
                    rotDirForTemplate, templateCase, 'GradientXRange030_OPRCenterStd');
            end
            TplLoad = load(templateFile, 'Template');
            LowSpeedTemplate = TplLoad.Template;
            DynamicMap = build_dynamic_map_from_experiment_case_local( ...
                highMap, LowSpeedTemplate, analysisSensors, C0, templateFile);
            save(dynamicMapCacheFile, 'DynamicMap', '-v7.3');
        end
        matFile = cacheMatFile;
        csvFile = cacheCsvFile;
        baselineCsvFile = cacheBaselineCsvFile;
        if ~isfile(legacyMatFile)
            save(legacyMatFile, 'highMap', '-v7.3');
        end
        return;
    end
end
edgeFraction = 0.10;
minPulsePoints = 12;
voltageAutoScaleThreshold = 20;  % assume V if max abs voltage is below this
dynamicBaselineRevCount = 2;
dynamicBaselinePulsePadS = 8e-4;
dynamicBaselinePulseThresholdFraction = 0.12;
dynamicBaselineSegmentTrimFraction = 0.12;
dynamicBaselineMinSegmentS = 3e-4;

rotDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration');
addpath(rotDir);
templateCase = struct();
templateCase.dataset = '20251222';
templateCase.bladeId = targetBlade;
templateCase.sensorIds = analysisSensors;
templateCase.sensorTag = ['S', sprintf('%d', analysisSensors)];
templateCase.caseTag = sprintf('B%d_%s', targetBlade, templateCase.sensorTag);
[templateFile, templateSourceMode] = resolveTemplateFile_OPRCenterStd_20251222(rotDir, templateCase, 'GradientXRange030_OPRCenterStd');
if ~isfile(templateFile)
    error('Missing low-speed template for coordinate alignment: %s', templateFile);
end
TplLoad = load(templateFile, 'Template');
LowSpeedTemplate = TplLoad.Template;

analysisStartTime = manualAnalysisStartTime;
if ~isempty(regionTable)
    if ismember('regionId', regionTable.Properties.VariableNames)
        regionIndex = find(regionTable.regionId == resonanceRegionId, 1, 'first');
    else
        regionIndex = [];
    end
    if isempty(regionIndex)
        regionIndex = 1;
    end
    region = regionTable(regionIndex, :);
else
    regionIndex = resonanceRegionId;
    region = table(Pflow.identification.resonanceRegionStartSec, ...
        Pflow.identification.resonanceRegionEndSec, analysisStartTime, ...
        resonanceRegionId, NaN, referenceOrder, ...
        'VariableNames', {'regionStart', 'regionEnd', 'regionPeakTime', ...
        'regionId', 'regionPeakFreqHz', 'regionPeakOrder'});
end

%% 2. Build compact experiment case from legacy waveform windows
baseCfg = struct();
baseCfg.dataset_root = cfgLegacy.dataset_root;
baseCfg.low_speed_case = cfgLegacy.low_speed_case;
baseCfg.dynamic_case_name = caseName;
baseCfg.output_root = cfgLegacy.output_root;
baseCfg.reference_output_dir = cfgLegacy.reference_output_dir;
baseCfg.sensor_ids = cfgLegacy.sensor_ids;
baseCfg.opr_id = cfgLegacy.opr_id;
baseCfg.blades_num = cfgLegacy.blades_num;
baseCfg.pinlv = cfgLegacy.pinlv;
baseCfg.r_tip_mm = cfgLegacy.r_tip_mm;
baseCfg.opr_threshold = cfgLegacy.opr_threshold;
baseCfg.sensor_threshold_default = cfgLegacy.sensor_threshold_default;

cfg = struct();
cfg.project_dir = thisDir;
cfg.base_cfg = baseCfg;
cfg.reference_case_name = cfgLegacy.low_speed_case;
cfg.dynamic_case_name = caseName;
cfg.static_data_dir = fullfile(cfgLegacy.dataset_root, cfgLegacy.low_speed_case);
cfg.dynamic_data_dir = fullfile(cfgLegacy.dataset_root, caseName);
cfg.reference_output_dir = cfgLegacy.reference_output_dir;
cfg.case_output_dir = fullfile(cfgLegacy.output_root, caseName);
cfg.target_blade = targetBlade;
cfg.analysis_sensors = analysisSensors;
cfg.analysis_start_time = analysisStartTime;
cfg.target_laps = targetLaps;
cfg.num_blades = cfgLegacy.blades_num;
cfg.opr_pulses_per_rev = cfgLegacy.blades_num;
cfg.pinlv = cfgLegacy.pinlv;
cfg.r_tip_mm = cfgLegacy.r_tip_mm;
cfg.opr_channel = cfgLegacy.opr_id;
cfg.opr_threshold = cfgLegacy.opr_threshold;
cfg.gap_threshold = cfgLegacy.gap_points;
cfg.pulse_window_sec = pulseWindowSec;
cfg.pulse_pad_sec = 2 / cfg.pinlv;
cfg.dynamic_window_mode = 'peak_centered_fixed';
cfg.sensor_config_file = fullfile(cfgLegacy.reference_output_dir, 'Sensor_Config_20251222.mat');
cfg.show_figures = false;
cfg.save_experiment_case = false;

experiment_case = prepare_identification_case_single_sync_20251222(cfg, false, false);

%% 3. Estimate one dynamic continuous-waveform baseline per sensor
dynamicBaselineRows = cell(numel(analysisSensors), 1);
dynamicBaselineBySensor = containers.Map('KeyType', 'double', 'ValueType', 'any');
for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    if sid > numel(experiment_case.Raw_Stream) || isempty(experiment_case.Raw_Stream(sid).T)
        continue;
    end
    tRaw = experiment_case.Raw_Stream(sid).T(:);
    vRaw = experiment_case.Raw_Stream(sid).V(:);
    [baselineInfo, baselineMask] = estimate_dynamic_sensor_baseline_local( ...
        tRaw, vRaw, cfg.pinlv, experiment_case.Diagnostics.rot_freq_mean_hz, ...
        dynamicBaselineRevCount, dynamicBaselinePulsePadS, ...
        dynamicBaselinePulseThresholdFraction, dynamicBaselineSegmentTrimFraction, ...
        dynamicBaselineMinSegmentS, voltageAutoScaleThreshold);
    dynamicBaselineBySensor(sid) = baselineInfo;
    dynamicBaselineRows{is} = table(sid, baselineInfo.baselineOriginalUnit, ...
        baselineInfo.baselineMv, baselineInfo.backgroundMedianMv, ...
        baselineInfo.backgroundFraction, baselineInfo.windowStartS, ...
        baselineInfo.windowEndS, baselineInfo.windowPeakCount, ...
        baselineInfo.voltageScale, ...
        'VariableNames', {'sensorId', 'baselineOriginalUnit', 'baselineMv', ...
        'backgroundMedianMv', 'backgroundFraction', 'baselineWindowStartS', ...
        'baselineWindowEndS', 'baselineWindowPeakCount', 'voltageScale'});
    %#ok<NASGU>
end
dynamicBaselineTable = vertcat(dynamicBaselineRows{~cellfun(@isempty, dynamicBaselineRows)});

%% 4. Convert experiment_case to highMap
xGrid = responseSurface.xGrid(:);
xMin = min(xGrid(responseSurface.effectiveWindow));
xMax = max(xGrid(responseSurface.effectiveWindow));
t_v = [];
x_v = [];
x_fit_v = [];
x_timebased_ref_v = [];
x_mm_v = [];
x_peak_mm_v = [];
V_raw_mV = [];
V_a = [];
baseline_mV_v = [];
S_v = [];
rev_v = [];
W_v = [];
theta_v = [];
pulsePeakTime_v = [];
windowRows = {};

for is = 1:numel(analysisSensors)
    sid = analysisSensors(is);
    dataIdx = find([experiment_case.Extracted_Data.sensor_id] == sid, 1, 'first');
    if isempty(dataIdx)
        continue;
    end
    laps = experiment_case.Extracted_Data(dataIdx).Laps;
    for ilap = 1:numel(laps)
        D = laps(ilap);
        if isempty(D.t_points) || numel(D.t_points) < minPulsePoints
            continue;
        end
        t = D.t_points(:);
        v = D.v_points(:);
        xMm = D.x_points(:);
        [~, idxPeak] = max(v);
        if isfield(D, 'peak_time') && isfinite(D.peak_time)
            peakTime = D.peak_time;
        else
            peakTime = t(idxPeak);
        end

        if ~isKey(dynamicBaselineBySensor, sid)
            error('Missing dynamic baseline for CH%d.', sid);
        end
        baselineInfo = dynamicBaselineBySensor(sid);
        baseline = baselineInfo.baselineOriginalUnit;
        voltageScale = baselineInfo.voltageScale;
        vRawMv = voltageScale * v;
        baselineMv = baselineInfo.baselineMv;
        vCorr = voltageScale * (v - baseline);

        edgeCount = max(3, round(edgeFraction * numel(v)));
        firstEdge = v(1:edgeCount);
        lastEdge = v(end-edgeCount+1:end);
        firstEdgeMedianMv = voltageScale * median(firstEdge, 'omitnan');
        lastEdgeMedianMv = voltageScale * median(lastEdge, 'omitnan');

        tau = t - peakTime;
        xTimeBasedRef = tau .* referenceTipSpeedMmS;
        if all(~isfinite(xMm))
            xPeakMm = 0;
        else
            xPeakMm = interp1(t, xMm, peakTime, 'linear', 'extrap');
        end
        Tpl = get_low_speed_template_sensor_local(LowSpeedTemplate, sid);
        xInstallAlignedMm = xMm - Tpl.xc;
        xPeakRelativeToTemplateMm = xPeakMm - Tpl.xc;
        xTimeBasedLocalMm = xTimeBasedRef - xPeakRelativeToTemplateMm;
        xModelMin = max(xMin, min(Tpl.x_grid(:)));
        xModelMax = min(xMax, max(Tpl.x_grid(:)));
        keep = xInstallAlignedMm >= xModelMin & xInstallAlignedMm <= xModelMax & ...
            isfinite(xInstallAlignedMm) & isfinite(vCorr);
        if nnz(keep) < minPulsePoints
            continue;
        end

        vSmooth = movmedian(vCorr, min(11, nnz(keep)));
        dv = abs(gradient(vSmooth, t));
        w = dv ./ max(dv + eps);
        w = max(0.05, w);

        nKeep = nnz(keep);
        t_v = [t_v; t(keep)]; %#ok<AGROW>
        x_v = [x_v; xInstallAlignedMm(keep)]; %#ok<AGROW>
        x_fit_v = [x_fit_v; xInstallAlignedMm(keep)]; %#ok<AGROW>
        x_timebased_ref_v = [x_timebased_ref_v; xTimeBasedLocalMm(keep)]; %#ok<AGROW>
        x_mm_v = [x_mm_v; xMm(keep)]; %#ok<AGROW>
        x_peak_mm_v = [x_peak_mm_v; xPeakRelativeToTemplateMm * ones(nKeep, 1)]; %#ok<AGROW>
        V_raw_mV = [V_raw_mV; vRawMv(keep)]; %#ok<AGROW>
        V_a = [V_a; vCorr(keep)]; %#ok<AGROW>
        baseline_mV_v = [baseline_mV_v; baselineMv * ones(nKeep, 1)]; %#ok<AGROW>
        S_v = [S_v; sid * ones(nKeep, 1)]; %#ok<AGROW>
        rev_v = [rev_v; ilap * ones(nKeep, 1)]; %#ok<AGROW>
        W_v = [W_v; w(keep)]; %#ok<AGROW>
        theta_v = [theta_v; 2*pi*experiment_case.Diagnostics.rot_freq_mean_hz .* (t(keep) - t(1))]; %#ok<AGROW>
        pulsePeakTime_v = [pulsePeakTime_v; peakTime * ones(nKeep, 1)]; %#ok<AGROW>

        dynamicTipSpeedMmS = 2 * pi * cfg.r_tip_mm * experiment_case.Diagnostics.rot_freq_mean_hz;
        windowRows{end+1, 1} = table(sid, ilap, min(t), max(t), peakTime, xPeakMm, nKeep, ...
            min(xInstallAlignedMm(keep)), max(xInstallAlignedMm(keep)), dynamicTipSpeedMmS, ...
            max(vRawMv, [], 'omitnan'), max(vCorr, [], 'omitnan'), ...
            baseline, baselineMv, firstEdgeMedianMv, lastEdgeMedianMv, ...
            abs(firstEdgeMedianMv - lastEdgeMedianMv), ...
            'VariableNames', {'sensorId', 'lapId', 'timeStart', 'timeEnd', ...
            'peakTime', 'xPeakMm', 'validPoints', 'xMinMm', 'xMaxMm', 'dynamicTipSpeedMmS', ...
            'peakRawMv', 'peakMvBaselineCorrected', 'rawBaselineOriginalUnit', ...
            'baselineMv', 'firstEdgeMedianMv', 'lastEdgeMedianMv', ...
            'edgeMedianDifferenceMv'});
    end
end

if isempty(t_v)
    error('No valid highMap samples were constructed.');
end

W_v = W_v ./ max(W_v);
W_v = max(0.05, W_v);

highMap = struct();
highMap.dataset = '20251222';
highMap.caseName = caseName;
highMap.targetBlade = targetBlade;
highMap.analysisSensors = analysisSensors;
highMap.regionIndex = regionIndex;
highMap.region = region;
highMap.resonanceRegionId = resonanceRegionId;
highMap.resonanceRegionTag = resonanceRegionTag;
highMap.resonanceRegionShortTag = resonanceRegionShortTag;
highMap.analysisStartMode = analysisStartMode;
highMap.analysisStartTime = analysisStartTime;
highMap.targetLaps = targetLaps;
highMap.analysisWinSize = analysisWinSize;
highMap.slidingStep = slidingStep;
highMap.pulseWindowSec = pulseWindowSec;
highMap.cacheKey = highMapCacheKey;
highMap.cacheParameterNote = ['Cache key is based on region, analysis start time, ' ...
    'target laps, window laps, sliding step and pulse window.'];
highMap.t_v = t_v;
highMap.x_v = x_v;
highMap.x_fit_v = x_fit_v;
highMap.x_timebased_ref_v = x_timebased_ref_v;
highMap.x_mm_v = x_mm_v;
highMap.x_peak_mm_v = x_peak_mm_v;
highMap.V_raw_mV = V_raw_mV;
highMap.V_a = V_a;
highMap.baseline_mV_v = baseline_mV_v;
highMap.S_v = S_v;
highMap.rev_v = rev_v;
highMap.W_v = W_v;
highMap.theta_v = theta_v;
highMap.pulsePeakTime_v = pulsePeakTime_v;
highMap.pointCount = numel(t_v);
highMap.rotFreqHz = experiment_case.Diagnostics.rot_freq_mean_hz;
highMap.rotRpm = 60 * highMap.rotFreqHz;
highMap.dynamicTipSpeedMmS = 2 * pi * cfg.r_tip_mm * highMap.rotFreqHz;
highMap.calibrationReferenceRpm = referenceRpm;
highMap.calibrationReferenceTipSpeedMmS = referenceTipSpeedMmS;
highMap.calibrationReferenceTipSpeedSource = referenceSpeedSource;
highMap.dynamicBaselineBySensor = dynamicBaselineTable;
highMap.responseSurfaceXMinMm = xMin;
highMap.responseSurfaceXMaxMm = xMax;
highMap.referenceOrder = referenceOrder;
highMap.referenceFreqHz = referenceOrder * highMap.rotFreqHz;
highMap.referenceNote = ['Reference EO is fixed to EO14 to match the low-speed rotating-template ' ...
    'validation route that starts at 50.2 s and uses 20 laps with 3-lap sliding windows.'];
highMap.lowSpeedTemplateFile = templateFile;
highMap.OPRReference = experiment_case.Diagnostics.OPRReference;
highMap.standardAngleReference = experiment_case.Diagnostics.standard_angle_reference;
highMap.standardAnglePhaseShiftDeg = experiment_case.Diagnostics.standard_angle_phase_shift_deg;
highMap.standardAnglePhaseShiftMm = experiment_case.Diagnostics.standard_angle_phase_shift_mm;
highMap.sensorConfigStandardAngleReference = read_sensor_config_angle_reference_local(experiment_case.Sensor_Config);
highMap.xCoordinateNote = ['x_v and x_fit_v are OPR/installation-based coordinates after subtracting the low-speed ' ...
    'template center Tpl.xc. They are not re-centered by each dynamic pulse peak. x_peak_mm_v stores the dynamic ' ...
    'pulse peak coordinate relative to the same low-speed template center for diagnostics only. ' ...
    'x_timebased_ref_v is a peak-local time-speed coordinate retained for audit only.'];
highMap.voltageUnit = 'mV baseline-corrected; V_raw_mV stores raw voltage after the same mV scaling';
highMap.baselineNote = sprintf(['Dynamic baseline uses the mean of trimmed no-pulse background samples from the first %d revolutions ' ...
    'of each sensor continuous waveform. Step05 static baseline uses the analogous initial-revolution continuous-waveform background method.'], ...
    dynamicBaselineRevCount);
highMap.experiment_case = experiment_case;
if isempty(windowRows)
    highMap.windowTable = table();
else
    highMap.windowTable = vertcat(windowRows{:});
end
DynamicMap = build_dynamic_map_from_experiment_case_local( ...
    highMap, LowSpeedTemplate, analysisSensors, C0, templateFile);

%% 5. Save outputs
matFile = cacheMatFile;
csvFile = cacheCsvFile;
baselineCsvFile = cacheBaselineCsvFile;
save(matFile, 'highMap', '-v7.3');
save(dynamicMapCacheFile, 'DynamicMap', '-v7.3');
writetable(highMap.windowTable, csvFile);
writetable(highMap.dynamicBaselineBySensor, baselineCsvFile);
if ~strcmp(matFile, legacyMatFile)
    save(legacyMatFile, 'highMap', '-v7.3');
    writetable(highMap.windowTable, legacyCsvFile);
    writetable(highMap.dynamicBaselineBySensor, legacyBaselineCsvFile);
end

%% 6. Visualization
fig = figure('Name', '20251222 Step06 highMap', ...
    'Color', 'w', 'Position', [80, 80, 1350, 820], 'NumberTitle', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
for sid = analysisSensors
    mask = highMap.S_v == sid;
    scatter(highMap.t_v(mask), highMap.V_raw_mV(mask), 8, 'filled', ...
        'DisplayName', sprintf('CH%d raw', sid));
    plot(highMap.t_v(mask), highMap.baseline_mV_v(mask), '.', ...
        'MarkerSize', 3, 'Color', [0.35, 0.35, 0.35], ...
        'HandleVisibility', 'off');
end
xlabel('Time (s)');
ylabel('Raw voltage (mV)');
title('Selected dynamic raw samples and continuous no-pulse baselines');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
for sid = analysisSensors
    mask = highMap.S_v == sid;
    scatter(highMap.x_v(mask), highMap.V_a(mask), 8, 'filled', ...
        'DisplayName', sprintf('CH%d', sid));
end
xlabel('Circumferential displacement x (mm)');
ylabel('Voltage (mV)');
title('Mapped highMap samples');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
histogram(highMap.x_v, 50);
xlabel('Circumferential displacement x (mm)');
ylabel('Count');
title('Coordinate coverage');

nexttile; hold on; grid on; box on;
colorsGap = turbo(numel(responseSurface.gTrainMm));
for ig = 1:numel(responseSurface.gTrainMm)
    plot(responseSurface.xGrid, responseSurface.Yfit(:, ig), ...
        '-', 'Color', colorsGap(ig, :), 'LineWidth', 0.9, ...
        'DisplayName', sprintf('g=%.2f mm', responseSurface.gTrainMm(ig)));
end
yline(0, '--k', 'Static/dynamic corrected baseline', ...
    'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
scatter(highMap.x_v, highMap.V_a, 6, [0.2, 0.45, 0.8], 'filled', 'DisplayName', 'Dynamic samples');
xlabel('Circumferential displacement x (mm)');
ylabel('Voltage (mV)');
title('Dynamic samples against all static gap templates');
legend('Location', 'eastoutside');

figFile = fullfile(outDir, sprintf('Step06_HighMap_%s_%s_%s.png', ...
    C0.dataset, C0.caseTag, highMapCacheKey));
saveas(fig, figFile);

fprintf('\nStep06 complete.\n');
fprintf('highMap saved to:\n  %s\n', matFile);
fprintf('Window summary saved to:\n  %s\n', csvFile);
fprintf('Dynamic continuous-baseline summary saved to:\n  %s\n', baselineCsvFile);

%% Local functions
function key = make_highmap_cache_key_local(regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec)
key = sprintf('%s_T%09.3f_L%02d_W%02d_S%02d_PW%07.5f', ...
    regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec);
key = strrep(key, '.', 'p');
key = strrep(key, '+', '');
key = strrep(key, '-', 'm');
end

function tf = is_highmap_cache_match_local(highMap, C0, regionShortTag, startTimeSec, targetLaps, windowLaps, slidingStepLaps, pulseWindowSec)
tol = 1e-9;
tf = isstruct(highMap) && ...
    isfield(highMap, 'dataset') && strcmp(char(string(highMap.dataset)), C0.dataset) && ...
    isfield(highMap, 'targetBlade') && highMap.targetBlade == C0.bladeId && ...
    isfield(highMap, 'analysisSensors') && isequal(reshape(highMap.analysisSensors, 1, []), reshape(C0.sensorIds, 1, [])) && ...
    isfield(highMap, 'resonanceRegionShortTag') && strcmp(char(string(highMap.resonanceRegionShortTag)), regionShortTag) && ...
    isfield(highMap, 'analysisStartTime') && abs(highMap.analysisStartTime - startTimeSec) <= tol && ...
    isfield(highMap, 'targetLaps') && highMap.targetLaps == targetLaps && ...
    isfield(highMap, 'analysisWinSize') && highMap.analysisWinSize == windowLaps && ...
    isfield(highMap, 'slidingStep') && highMap.slidingStep == slidingStepLaps && ...
    isfield(highMap, 'pulseWindowSec') && abs(highMap.pulseWindowSec - pulseWindowSec) <= tol;
end

function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRCenterStd_20251222.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_OPRAnchored_20251222.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_BaseFrame_20251222.mat', targetBlade, exactTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20251222.mat', targetBlade)));
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)(?:_GradientXRange030(?:_[A-Za-z0-9]+)?)?_20251222\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    available = arrayfun(@(c) str2double(c), token{1});
    if all(ismember(sensorIds, available))
        templateFile = fullfile(files(i).folder, files(i).name);
        return;
    end
end
templateFile = fullfile(templateDir, patterns{1});
end

function Tpl = get_low_speed_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Low-speed template does not contain CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function DynamicMap = build_dynamic_map_from_experiment_case_local(highMap, Template, analysisSensors, C0, templateFile)
E = highMap.experiment_case;
windowLaps = highMap.analysisWinSize;
slidingStep = highMap.slidingStep;
targetLaps = highMap.targetLaps;
starts = 1:slidingStep:(targetLaps - windowLaps + 1);
Window = repmat(struct( ...
    'window_id', NaN, 'lap_range', [], 'time_window', [], ...
    'rot_freq_mean_hz', NaN, 'rot_rpm_mean', NaN, ...
    'Sensor', []), numel(starts), 1);

for iw = 1:numel(starts)
    lapRange = starts(iw):(starts(iw) + windowLaps - 1);
    Sensor = repmat(struct( ...
        'sensor_id', NaN, 't', [], 'x_abs', [], 'x_rel', [], ...
        'V', [], 'W', [], 'theta', [], 'point_count', NaN), numel(analysisSensors), 1);
    allT = [];
    for is = 1:numel(analysisSensors)
        sid = analysisSensors(is);
        dataIdx = find([E.Extracted_Data.sensor_id] == sid, 1, 'first');
        if isempty(dataIdx)
            error('experiment_case.Extracted_Data does not contain CH%d.', sid);
        end
        Tpl = get_low_speed_template_sensor_local(Template, sid);
        laps = E.Extracted_Data(dataIdx).Laps;
        t = [];
        xAbs = [];
        V = [];
        theta = [];
        for lapId = lapRange
            if lapId > numel(laps)
                continue;
            end
            D = laps(lapId);
            t = [t; D.t_points(:)]; %#ok<AGROW>
            xAbs = [xAbs; D.x_points(:)]; %#ok<AGROW>
            V = [V; D.v_points(:)]; %#ok<AGROW>
            theta = [theta; infer_lap_theta_for_dynamic_map_local(E, D, lapId)]; %#ok<AGROW>
        end
        Sensor(is).sensor_id = sid;
        Sensor(is).t = t(:);
        Sensor(is).x_abs = xAbs(:);
        Sensor(is).x_rel = xAbs(:) - Tpl.xc;
        Sensor(is).V = V(:);
        Sensor(is).W = build_simple_waveform_weight_for_dynamic_map_local(V);
        Sensor(is).theta = theta(:);
        Sensor(is).point_count = numel(t);
        allT = [allT; t(:)]; %#ok<AGROW>
    end
    finiteT = allT(isfinite(allT));
    Window(iw).window_id = iw;
    Window(iw).lap_range = lapRange;
    if isempty(finiteT)
        Window(iw).time_window = [NaN, NaN];
    else
        Window(iw).time_window = [min(finiteT), max(finiteT)];
    end
    Window(iw).rot_freq_mean_hz = highMap.rotFreqHz;
    Window(iw).rot_rpm_mean = highMap.rotRpm;
    Window(iw).Sensor = Sensor;
end

DynamicMap = struct();
DynamicMap.Route = 'gap_prior_step06_legacy_step02_dynamic_map';
DynamicMap.TargetBlade = highMap.targetBlade;
DynamicMap.SensorIDs = analysisSensors;
DynamicMap.SensorTag = ['S', sprintf('%d', analysisSensors)];
DynamicMap.RegionShortTag = highMap.resonanceRegionShortTag;
DynamicMap.SourceSettings = struct( ...
    'analysis_start_time', highMap.analysisStartTime, ...
    'target_laps', highMap.targetLaps, ...
    'analysis_win_size', highMap.analysisWinSize, ...
    'sliding_step', highMap.slidingStep, ...
    'pulse_window_sec', highMap.pulseWindowSec);
DynamicMap.SourceHighMapCacheKey = highMap.cacheKey;
DynamicMap.TemplateFileForXRel = templateFile;
DynamicMap.XCenterBySensor = build_xcenter_table_for_dynamic_map_local(Template, analysisSensors);
DynamicMap.OPRReference = highMap.OPRReference;
DynamicMap.GlobalTimeWindow = infer_dynamic_map_global_time_window_local(Window);
DynamicMap.Window = Window;
DynamicMap.Dataset = C0.dataset;
end

function theta = infer_lap_theta_for_dynamic_map_local(E, D, lapId)
t = D.t_points(:);
if isfield(D, 'theta_points') && numel(D.theta_points) == numel(t)
    theta = D.theta_points(:);
elseif isfield(E, 'Diagnostics') && isfield(E.Diagnostics, 'rot_freq_mean_hz') && ...
        isfinite(E.Diagnostics.rot_freq_mean_hz)
    theta = 2 * pi * E.Diagnostics.rot_freq_mean_hz .* (t - t(1));
else
    theta = 2 * pi * (lapId - 1) + linspace(0, 2*pi, numel(t)).';
end
end

function W = build_simple_waveform_weight_for_dynamic_map_local(V)
V = V(:);
if isempty(V)
    W = V;
    return;
end
W = ones(size(V));
finite = isfinite(V);
if nnz(finite) >= 5
    g = abs(gradient(V));
    gmax = max(g(finite), [], 'omitnan');
    if isfinite(gmax) && gmax > 0
        W = g ./ gmax;
    end
end
W = max(0.05, W);
end

function T = build_xcenter_table_for_dynamic_map_local(Template, analysisSensors)
rows = repmat(struct('sensor_id', NaN, 'xc_mm', NaN), numel(analysisSensors), 1);
for i = 1:numel(analysisSensors)
    sid = analysisSensors(i);
    Tpl = get_low_speed_template_sensor_local(Template, sid);
    rows(i).sensor_id = sid;
    rows(i).xc_mm = Tpl.xc;
end
T = struct2table(rows);
end

function tw = infer_dynamic_map_global_time_window_local(Window)
t0 = NaN(numel(Window), 1);
t1 = NaN(numel(Window), 1);
for i = 1:numel(Window)
    if numel(Window(i).time_window) == 2
        t0(i) = Window(i).time_window(1);
        t1(i) = Window(i).time_window(2);
    end
end
tw = [min(t0, [], 'omitnan'), max(t1, [], 'omitnan')];
end

function angleRef = read_sensor_config_angle_reference_local(Sensor_Config)
angleRef = 'unknown';
if isstruct(Sensor_Config) && isfield(Sensor_Config, 'Standard_Relative_Angles_Reference')
    angleRef = char(string(Sensor_Config.Standard_Relative_Angles_Reference));
end
end

function [baselineInfo, backgroundMask] = estimate_dynamic_sensor_baseline_local(tRaw, vRaw, fs, rotFreqHz, revCount, pulsePadS, pulseThresholdFraction, segmentTrimFraction, minSegmentS, voltageAutoScaleThreshold)
tRaw = tRaw(:);
vRaw = vRaw(:);
windowStart = tRaw(1);
windowEnd = min(tRaw(end), windowStart + revCount / rotFreqHz);
windowMask = tRaw >= windowStart & tRaw <= windowEnd & isfinite(vRaw);
vWin = vRaw(windowMask);
if isempty(vWin)
    error('Empty dynamic baseline window.');
end
fullMedian = median(vWin, 'omitnan');
upper = prctile(vWin, 99.5);
threshold = fullMedian + pulseThresholdFraction * max(upper - fullMedian, eps);
pulseCore = vRaw > threshold & windowMask;
padN = max(1, round(pulsePadS * fs));
pulseVicinity = conv(double(pulseCore), ones(2 * padN + 1, 1), 'same') > 0;
backgroundMask0 = windowMask & ~pulseVicinity & isfinite(vRaw);
backgroundMask = trim_background_segments_local(backgroundMask0, fs, segmentTrimFraction, minSegmentS);
if nnz(backgroundMask) < 0.10 * nnz(backgroundMask0)
    backgroundMask = backgroundMask0;
end
bg = vRaw(backgroundMask);
baselineOriginal = mean(bg, 'omitnan');
baselineMedianOriginal = median(bg, 'omitnan');
voltageScale = 1;
if max(abs(vWin - baselineOriginal), [], 'omitnan') < voltageAutoScaleThreshold
    voltageScale = 1000;
end
[~, locs] = findpeaks(vWin, 'MinPeakHeight', threshold, ...
    'MinPeakDistance', max(1, round(0.5 / (rotFreqHz * 6) * fs)));
baselineInfo = struct();
baselineInfo.baselineOriginalUnit = baselineOriginal;
baselineInfo.baselineMv = voltageScale * baselineOriginal;
baselineInfo.backgroundMedianMv = voltageScale * baselineMedianOriginal;
baselineInfo.backgroundFraction = nnz(backgroundMask) / nnz(windowMask);
baselineInfo.windowStartS = windowStart;
baselineInfo.windowEndS = windowEnd;
baselineInfo.windowPeakCount = numel(locs);
baselineInfo.thresholdOriginalUnit = threshold;
baselineInfo.voltageScale = voltageScale;
end

function trimmedMask = trim_background_segments_local(mask, fs, trimFraction, minSegmentS)
mask = mask(:);
trimmedMask = false(size(mask));
idx = find(mask);
if isempty(idx)
    return;
end
breaks = [1; find(diff(idx) > 1) + 1; numel(idx) + 1];
minSegmentN = max(3, round(minSegmentS * fs));
for iseg = 1:numel(breaks)-1
    segIdx = idx(breaks(iseg):breaks(iseg+1)-1);
    n = numel(segIdx);
    if n < minSegmentN
        continue;
    end
    trimN = floor(trimFraction * n);
    if 2 * trimN >= n - 2
        trimN = max(0, floor((n - 2) / 2));
    end
    keepIdx = segIdx((1 + trimN):(n - trimN));
    trimmedMask(keepIdx) = true;
end
end
