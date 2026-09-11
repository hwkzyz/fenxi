function report = Check_Program_20241106(mode, verbose)
%CHECK_PROGRAM_20241106 Verify code, bundled inputs and path independence.
if nargin < 1
    mode = 'all';
    verbose = true;
elseif islogical(mode)
    verbose = mode;
    mode = 'all';
elseif nargin < 2
    verbose = true;
end
mode = lower(string(mode));
assert(ismember(mode, ["code","inputs","prepared","results","all"]), ...
    'Check mode must be code, inputs, prepared, results, or all.');

cfg = Config_20241106();
addpath(cfg.paths.root, cfg.paths.foundation, cfg.paths.preparation, ...
    cfg.paths.gapAware, cfg.paths.utilities);
requiredFiles = {
    fullfile(cfg.paths.root, 'Config_20241106.m')
    fullfile(cfg.paths.root, 'Main01_Build_BTT_Displacement_20241106.m')
    fullfile(cfg.paths.root, 'Main02_Prepare_Strain_Resonance_Evidence_20241106.m')
    fullfile(cfg.paths.root, 'Main03_Build_Static_Gap_Response_Surface_20241106.m')
    fullfile(cfg.paths.root, 'Main04_Correct_OffsetTilt_Response_Surface_20241106.m')
    fullfile(cfg.paths.root, 'Main05_Build_LowSpeed_TemplateBank_20241106.m')
    fullfile(cfg.paths.root, 'Main06_Calibrate_GapLibrary_20241106.m')
    fullfile(cfg.paths.root, 'Main07_Foundation_FixedGap_Identification_20241106.m')
    fullfile(cfg.paths.root, 'Main08_Build_GapAware_DynamicMap_20241106.m')
    fullfile(cfg.paths.root, 'Main09_GapAware_FullWave_Identification_20241106.m')
    fullfile(cfg.paths.root, 'Main10_Compare_FixedGap_GapAware_20241106.m')
    fullfile(cfg.paths.root, 'Main11_Estimate_Strain_BTT_TimeAlignment_20241106.m')
    fullfile(cfg.paths.root, 'Main12_Validate_Identification_With_Strain_20241106.m')
    fullfile(cfg.paths.root, 'Main13_Visualize_Strain_BTT_Waveforms_20241106.m')
    fullfile(cfg.paths.utilities, 'Apply_Config_Env_20241106.m')
    fullfile(cfg.paths.utilities, 'Prepare_Bundled_Inputs_20241106.m')
    fullfile(cfg.paths.foundation, 'BTTDataConfig_20241106.m')
    fullfile(cfg.paths.foundation, 'BTTProjectConfig_20241106.m')
    fullfile(cfg.paths.gapAware, 'ProjectionFlow_Config_20241106.m')
    fullfile(cfg.paths.gapAware, '+step07jcore', 'solve_gradient_displacement_vp_seed.m')
    };
bundledInputs = {
    fullfile(cfg.paths.bundledPreparedFoundation, 'step01_low_speed_reference', 'Sensor_Config_20241106.mat')
    fullfile(cfg.paths.bundledPreparedFoundation, 'step02_dynamic_btt', cfg.case.dynamicCase, 'Step02_Dynamic_BTT_Extraction_20241106.mat')
    fullfile(cfg.paths.bundledPreparedFoundation, 'step02_dynamic_btt', cfg.case.dynamicCase, 'jiluOPR.mat')
    fullfile(cfg.paths.calibrationFoundation, 'step04_low_speed_template', 'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat')
    fullfile(cfg.paths.bundledPreparedGap, 'Step02_BTT_Displacement_20241106.mat')
    fullfile(cfg.paths.calibrationGap, 'LowSpeedTemplateBank_20241106_B1toB6_S2357.mat')
    fullfile(cfg.paths.calibrationGap, 'GapCalibrationBank_20241106_B1toB6_S57.mat')
    fullfile(cfg.paths.staticGapWaveformLibrary, 'gap_waveform_manifest.csv')
    cfg.paths.bladeOffsetFile
    };
preparedInputs = {
    fullfile(cfg.paths.preparedFoundation, 'step01_low_speed_reference', 'Sensor_Config_20241106.mat')
    fullfile(cfg.paths.preparedFoundation, 'step02_dynamic_btt', cfg.case.dynamicCase, 'Step02_Dynamic_BTT_Extraction_20241106.mat')
    fullfile(cfg.paths.preparedFoundation, 'step02_dynamic_btt', cfg.case.dynamicCase, 'jiluOPR.mat')
    fullfile(cfg.paths.preparedFoundation, 'step04_low_speed_template', 'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat')
    fullfile(cfg.paths.gapRuntime, 'Step02_BTT_Displacement_20241106.mat')
    };

missingCode = requiredFiles(~cellfun(@isfile, requiredFiles));
missingInputs = bundledInputs(~cellfun(@isfile, bundledInputs));
missingPrepared = preparedInputs(~cellfun(@isfile, preparedInputs));
if ismember(mode, ["code","all"])
    assert(isempty(missingCode), 'Missing package program: %s', strjoin(missingCode, newline));
end

resultFiles = struct('foundationTrend', "", 'gapTrend', "", 'comparisonSummary', "");
if ismember(mode, ["results", "all"])
    timeTag = sprintf('T%03dp%03d', floor(cfg.case.analysisStartTimeSec), ...
        round(1000 * mod(cfg.case.analysisStartTimeSec, 1)));
    foundationHits = dir(fullfile(cfg.paths.foundationResults, '**', ...
        'Trend_Step05_FoundationMainPulseAdaptiveZeroEta_B4_S257_20241106.csv'));
    foundationHits = foundationHits(contains(string({foundationHits.folder}), timeTag));
    gapHits = dir(fullfile(cfg.paths.gapResults, sprintf( ...
        '*Trend*%s*W3S1_F300to1000Hz_DFpm2Hz.csv', timeTag)));
    compareHits = dir(fullfile(cfg.paths.comparison, sprintf( ...
        '*Summary*%s*foundation_vs_gapaware_plainrmse.csv', timeTag)));
    assert(~isempty(foundationHits), 'Missing current-time Foundation trend for %s.', timeTag);
    assert(~isempty(gapHits), 'Missing current-time GapAware trend for %s.', timeTag);
    assert(~isempty(compareHits), 'Missing current-time comparison summary for %s.', timeTag);
    resultFiles.foundationTrend = string(fullfile(foundationHits(1).folder, foundationHits(1).name));
    resultFiles.gapTrend = string(fullfile(gapHits(1).folder, gapHits(1).name));
    resultFiles.comparisonSummary = string(fullfile(compareHits(1).folder, compareHits(1).name));
    Tf = readtable(resultFiles.foundationTrend);
    Tg = readtable(resultFiles.gapTrend);
    Tc = readtable(resultFiles.comparisonSummary, 'TextType', 'string');
    assert(height(Tf) == cfg.window.count, ...
        'Foundation trend has %d rows; expected %d.', height(Tf), cfg.window.count);
    assert(height(Tg) == cfg.window.count, ...
        'GapAware trend has %d rows; expected %d.', height(Tg), cfg.window.count);
    fixedRow = find(strcmpi(string(Tc.method), "fixed"), 1, 'first');
    assert(~isempty(fixedRow), 'Comparison summary has no fixed branch.');
    assert(abs(Tc.eoMatchFractionVsFoundation(fixedRow) - 1) < 1e-12, ...
        'GapAware fixed branch does not reproduce Foundation EO window by window.');
    assert(Tc.meanAbsAmplitudeDeltaVsFoundationMm(fixedRow) < 1e-10 && ...
        abs(Tc.meanRmseDeltaVsFoundationMv(fixedRow)) < 1e-8, ...
        'GapAware fixed branch does not numerically reproduce Foundation.');
end
if ismember(mode, ["inputs","all"])
    assert(isempty(missingInputs), 'Missing bundled input: %s', strjoin(missingInputs, newline));
end
if ismember(mode, ["prepared","all"])
    assert(isempty(missingPrepared), 'Missing active prepared input: %s', ...
        strjoin(missingPrepared, newline));
end
assert(all(ismember(cfg.case.gapSensors, [5 7])), 'Only CH5/CH7 may be gap sensors.');
assert(all(ismember(cfg.case.gapSensors, cfg.case.analysisSensors)), ...
    'Gap sensors must be included in analysisSensors.');
assert(strcmpi(cfg.model.etaPolicy, 'zero'), 'This package requires eta=0.');
assert(strcmpi(cfg.objective.finalType, 'plain_rmse'), ...
    'The formal final objective must be plain_rmse.');
assert(cfg.window.count == 18, 'Default 20/3/1 window settings must produce 18 windows.');

forbiddenGap = ['20241106_low_speed_gap_prior_', 'decoupling'];
forbiddenFoundation = ['20241106_btt_data_', 'foundation'];
forbiddenLegacyCalibration = '20241106_low_speed_rotating_calibration';
if ismember(mode, ["code","all"])
    assert(~isfolder(fullfile(cfg.paths.root, 'src')), ...
        'The hidden src folder must not exist in the formal package.');
    assert(isempty(dir(fullfile(cfg.paths.root, 'Run_*.m'))), ...
        'Legacy Run wrappers must not remain in the package root.');
    for i = 1:numel(requiredFiles)
        codeText = fileread(requiredFiles{i});
        assert(~contains(codeText, forbiddenGap), ...
            'Runtime code refers to the old gap folder: %s', requiredFiles{i});
        assert(~contains(codeText, forbiddenFoundation), ...
            'Runtime code refers to the old Foundation folder: %s', requiredFiles{i});
        assert(~contains(codeText, forbiddenLegacyCalibration), ...
            'Runtime code refers to the old rotating-calibration folder: %s', requiredFiles{i});
    end
    foundationText = fileread(fullfile(cfg.paths.root, ...
        'Main07_Foundation_FixedGap_Identification_20241106.m'));
    assert(contains(foundationText, 'continuousFrequencyRefine'), ...
        'Foundation route has no continuous-frequency refinement code.');
    assert(contains(foundationText, 'residual_obj = nansum(res(valid).^2)'), ...
        'Foundation route is not using the formal plain voltage residual objective.');
    assert(contains(foundationText, 'S.vp_top_k_max = 3;') && ...
        contains(foundationText, 'S.vp_gap_ratio_keep = 0;') && ...
        contains(foundationText, 'S.include_neighbor_eo = false;'), ...
        'Foundation route is not frozen to strict Top-3.');
    assert(contains(foundationText, 'Result.MethodContract = struct('), ...
        'Foundation result does not save MethodContract metadata.');
    gapText = fileread(fullfile(cfg.paths.root, ...
        'Main09_GapAware_FullWave_Identification_20241106.m'));
    assert(contains(gapText, 'assert(~phaseInfo.usedPreviousWindow'), ...
        'GapAware does not guard against previous-window phase information.');
    assert(contains(gapText, 'Result.MethodContract = struct('), ...
        'GapAware result does not save MethodContract metadata.');
    templateBankText = fileread(fullfile(cfg.paths.root, ...
        'Main05_Build_LowSpeed_TemplateBank_20241106.m'));
    assert(contains(templateBankText, 'bundled_step04_common_low_speed_window'), ...
        'Main05 is not using the bundled common-window template source.');
    removedOrderToken = ['reference', 'Order'];
    for i = 1:numel(requiredFiles)
        file = requiredFiles{i};
        assert(~contains(fileread(file), removedOrderToken), ...
            'Removed reference-order prior still exists: %s', file);
    end
    resolved = {which('Config_20241106'), which('BTTDataConfig_20241106'), ...
        which('ProjectionFlow_Config_20241106'), which('Apply_Config_Env_20241106'), ...
        which('Main11_Estimate_Strain_BTT_TimeAlignment_20241106'), ...
        which('Main12_Validate_Identification_With_Strain_20241106'), ...
        which('Main13_Visualize_Strain_BTT_Waveforms_20241106')};
    assert(all(cellfun(@(p) startsWith(p, cfg.paths.root, 'IgnoreCase', true), resolved)), ...
        'MATLAB resolves one or more formal functions outside the new package.');
end

% Validate the provenance of the common low-speed template used by Main05.
templateFile = fullfile(cfg.paths.calibrationFoundation, 'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat');
if ismember(mode, ["inputs", "prepared", "all"])
    assert(isfile(templateFile), 'Missing formal Step04 template: %s', templateFile);
    T0 = load(templateFile, 'Template');
    assert(isfield(T0, 'Template') && isfield(T0.Template, 'SensorBlade'), ...
        'Formal Step04 template has no SensorBlade records.');
    q = T0.Template.SensorBlade;
    for sid = cfg.case.analysisSensors(:).'
        idx = find([q.blade_id] == cfg.case.targetBlade & [q.sensor_id] == sid, 1, 'first');
        assert(~isempty(idx), 'Formal template lacks B%d CH%d.', cfg.case.targetBlade, sid);
        assert(q(idx).stable_window_start_lap == cfg.calibration.commonLowSpeedLapRange(1) && ...
            q(idx).stable_window_end_lap == cfg.calibration.commonLowSpeedLapRange(2), ...
            'B%d CH%d does not use the configured common low-speed lap range.', ...
            cfg.case.targetBlade, sid);
        assert(strcmpi(string(q(idx).stable_window_reason), ...
            "common_window_from_package_config"), ...
            'B%d CH%d common-window provenance is not explicit.', cfg.case.targetBlade, sid);
        assert(strcmpi(string(q(idx).quality_status), "good"), ...
            'B%d CH%d formal low-speed template is not good quality.', cfg.case.targetBlade, sid);
    end
end

report = struct('ok', true, 'root', cfg.paths.root, ...
    'codeFileCount', numel(requiredFiles), 'formalMainCount', 13, ...
    'bundledInputCount', numel(bundledInputs), ...
    'preparedInputCount', numel(preparedInputs), ...
    'finalObjective', cfg.objective.finalType, 'resultFiles', resultFiles);
if verbose
    fprintf('Package check passed.\n');
    fprintf('Root: %s\n', report.root);
    fprintf('Mode: %s; programs: %d, bundled inputs: %d, prepared inputs: %d.\n', ...
        mode, report.codeFileCount, report.bundledInputCount, report.preparedInputCount);
    fprintf('Final objective: ordinary unweighted voltage RMSE.\n');
if ismember(mode, ["results", "all"])
    flowCfg = ProjectionFlow_Config_20241106();
    runtimeTemplateBank = fullfile(cfg.paths.calibrationRuntime, ...
        flowCfg.calibration.runtimeTemplateBankFile);
    runtimeGapBank = fullfile(cfg.paths.calibrationRuntime, ...
        flowCfg.calibration.runtimeGapBankFile);
    assert(isfile(runtimeTemplateBank), ...
        'Missing target-specific runtime template bank: %s', runtimeTemplateBank);
    assert(isfile(runtimeGapBank), ...
        'Missing target-specific runtime gap bank: %s', runtimeGapBank);
    Bt = load(runtimeTemplateBank, 'LowSpeedTemplateBank');
    Bg = load(runtimeGapBank, 'GapCalibrationBank');
    assert(isequal(Bt.LowSpeedTemplateBank.bladeIds(:).', cfg.case.targetBlade), ...
        'Runtime template bank scope is not exactly the configured target blade.');
    assert(isequal(Bg.GapCalibrationBank.bladeIds(:).', cfg.case.targetBlade), ...
        'Runtime gap bank scope is not exactly the configured target blade.');
        fprintf('Matched formal results: Foundation, GapAware and comparison are present.\n');
    end
end
end
