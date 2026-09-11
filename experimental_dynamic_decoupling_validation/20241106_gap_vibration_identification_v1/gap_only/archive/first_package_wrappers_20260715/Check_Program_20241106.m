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
assert(ismember(mode, ["code","inputs","prepared","all"]), ...
    'Check mode must be code, inputs, prepared, or all.');

cfg = Setup_Paths_20241106();
requiredFiles = {
    fullfile(cfg.paths.root, 'Config_20241106.m')
    fullfile(cfg.paths.root, 'Setup_Paths_20241106.m')
    fullfile(cfg.paths.root, 'Run_01_Prepare_Data_20241106.m')
    fullfile(cfg.paths.root, 'Run_02_Foundation_NoGap_20241106.m')
    fullfile(cfg.paths.root, 'Run_03_GapAware_20241106.m')
    fullfile(cfg.paths.root, 'Run_04_Compare_Results_20241106.m')
    fullfile(cfg.paths.utilities, 'Apply_Config_Env_20241106.m')
    fullfile(cfg.paths.utilities, 'Prepare_Bundled_Inputs_20241106.m')
    fullfile(cfg.paths.foundation, 'BTTDataConfig_20241106.m')
    fullfile(cfg.paths.foundation, 'BTTProjectConfig_20241106.m')
    fullfile(cfg.paths.foundation, 'Step05_SingleSync_DirectTemplate_Identification_20241106.m')
    fullfile(cfg.paths.preparation, 'Step06_BuildGapAwareDynamicMap_20241106.m')
    fullfile(cfg.paths.gapAware, 'ProjectionFlow_Config_20241106.m')
    fullfile(cfg.paths.gapAware, 'Step07J_NestedStaticWarp_VPFullWave_20241106.m')
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
if ismember(mode, ["code","all"])
    for i = 1:numel(requiredFiles)
        codeText = fileread(requiredFiles{i});
        assert(~contains(codeText, forbiddenGap), ...
            'Runtime code refers to the old gap folder: %s', requiredFiles{i});
        assert(~contains(codeText, forbiddenFoundation), ...
            'Runtime code refers to the old Foundation folder: %s', requiredFiles{i});
    end
    allCode = dir(fullfile(cfg.paths.root, '**', '*.m'));
    removedOrderToken = ['reference', 'Order'];
    for i = 1:numel(allCode)
        file = fullfile(allCode(i).folder, allCode(i).name);
        assert(~contains(fileread(file), removedOrderToken), ...
            'Removed reference-order prior still exists: %s', file);
    end
    resolved = {which('Config_20241106'), which('BTTDataConfig_20241106'), ...
        which('ProjectionFlow_Config_20241106'), which('Apply_Config_Env_20241106')};
    assert(all(cellfun(@(p) startsWith(p, cfg.paths.root, 'IgnoreCase', true), resolved)), ...
        'MATLAB resolves one or more formal functions outside the new package.');
end

report = struct('ok', true, 'root', cfg.paths.root, ...
    'codeFileCount', numel(requiredFiles), ...
    'bundledInputCount', numel(bundledInputs), ...
    'preparedInputCount', numel(preparedInputs), ...
    'finalObjective', cfg.objective.finalType);
if verbose
    fprintf('Package check passed.\n');
    fprintf('Root: %s\n', report.root);
    fprintf('Mode: %s; programs: %d, bundled inputs: %d, prepared inputs: %d.\n', ...
        mode, report.codeFileCount, report.bundledInputCount, report.preparedInputCount);
    fprintf('Final objective: ordinary unweighted voltage RMSE.\n');
end
end
