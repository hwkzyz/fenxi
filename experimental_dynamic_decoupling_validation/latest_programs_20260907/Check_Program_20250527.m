function report = Check_Program_20250527(mode, verbose)
%CHECK_PROGRAM_20250527 Verify V1 code, inputs and frozen method settings.
if nargin < 1, mode = 'all'; end
if nargin < 2, verbose = true; end
mode = lower(string(mode));
assert(ismember(mode, ["code","inputs","all"]), ...
    'Check mode must be code, inputs, or all.');

cfg = Setup_Paths_20250527();
mainNames = arrayfun(@(i) sprintf('Main%02d_*_20250527.m', i), 1:13, ...
    'UniformOutput', false);
mainFiles = cell(13, 1);
for i = 1:13
    hit = dir(fullfile(cfg.paths.root, mainNames{i}));
    assert(numel(hit) == 1, 'Expected exactly one %s program.', mainNames{i});
    mainFiles{i} = fullfile(hit.folder, hit.name);
    assert(hit.bytes > 500, 'Main%02d is only a wrapper or is incomplete.', i);
end

requiredCode = [mainFiles; {
    fullfile(cfg.paths.root, 'Config_20250527.m')
    fullfile(cfg.paths.root, 'Setup_Paths_20250527.m')
    fullfile(cfg.paths.utilities, 'CaseConfig.m')
    fullfile(cfg.paths.utilities, 'ProjectionFlow_Config_20250527.m')
    fullfile(cfg.paths.foundation, 'BTTDataConfig_20250527.m')
    fullfile(cfg.paths.gapAware, '+step07jcore', 'solve_gradient_displacement_vp_seed.m')
    fullfile(cfg.paths.gapAware, '+step07jcore', 'build_support_aware_adaptive_sg_template.m')
    fullfile(cfg.paths.preparation, 'Prepare_AdaptiveSG_Template_20250527.m')
    }];
requiredInputs = {
    fullfile(cfg.paths.calibrationInputs, 'Step05_Response_Surface_20250527.mat')
    fullfile(cfg.paths.calibrationInputs, 'Step05I_OffsetTilt_Shared_Response_Surface_20250527.mat')
    fullfile(cfg.paths.calibrationInputs, 'Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136.mat')
    fullfile(cfg.paths.calibrationInputs, 'OPR6_SlotCalibration_20250527.mat')
    fullfile(cfg.paths.preparedInputs, 'foundation', 'step01_low_speed_reference', 'Sensor_Config_20250527.mat')
    fullfile(cfg.paths.preparedInputs, 'foundation', 'step02_dynamic_btt', ...
        cfg.case.dynamicCase, 'Step02_Dynamic_BTT_Extraction_20250527.mat')
    fullfile(cfg.paths.preparedInputs, 'foundation', 'step02_dynamic_btt', ...
        cfg.case.dynamicCase, 'jiluOPR.mat')
    fullfile(cfg.paths.preparedInputs, 'foundation', 'step04_low_speed_template', ...
        'Template_OPRCenterStd_LowSpeed_AllBlades_S136_20250527.mat')
    cfg.files.lowSpeedTemplateSource
    cfg.files.lowSpeedTemplate
    fullfile(cfg.paths.preparedInputs, 'gap_aware', 'Step05A_Average_Calibration_Speed.csv')
    };

if ismember(mode, ["code","all"])
    missing = requiredCode(~cellfun(@isfile, requiredCode));
    assert(isempty(missing), 'Missing program: %s', strjoin(missing, newline));
    forbidden = {'referenceOrder', '20250527_low_speed_gap_prior_decoupling', ...
        '20250527_btt_data_foundation', '20250527_low_speed_rotating_calibration'};
    for i = 1:numel(requiredCode)
        text = fileread(requiredCode{i});
        for j = 1:numel(forbidden)
            assert(~contains(text, forbidden{j}), ...
                'Forbidden legacy dependency/prior "%s" in %s.', ...
                forbidden{j}, requiredCode{i});
        end
    end
    gapCode = fileread(mainFiles{10});
    assert(contains(gapCode, "cfg.finalObjective = 'plain'"), ...
        'GapAware final objective is not frozen to plain RMSE.');
    assert(contains(gapCode, "cfg.prevWindowCandidateMode = 'off'"), ...
        'GapAware still shares EO candidates across windows.');
    assert(contains(gapCode, 'cfg.phaseSafeExpansion = false;'), ...
        'GapAware phase-safe expansion is not frozen off.');
    assert(contains(gapCode, 'Result.MethodContract = struct('), ...
        'GapAware result does not save MethodContract metadata.');
    assert(contains(gapCode, 'same corrected forward model') && ...
        ~contains(gapCode, "F0(mask) = voltage_to_template_mv_relative_local(templateV(mask), Tpl)"), ...
        'GapAware still inherits a stale Foundation template linearization.');
    assert(contains(gapCode, "cfg.dxReferenceMode = 'seed_median'"), ...
        'GapAware does not use a window-local dx reference.');
    assert(contains(gapCode, 'refine_continuous_frequency_fit_local'), ...
        'Per-candidate continuous-frequency refinement is missing.');
    foundationCode = fileread(mainFiles{9});
    assert(contains(foundationCode, 'direct_template_contfreq_objective_local'), ...
        'Foundation per-candidate continuous-frequency refinement is missing.');
    assert(contains(foundationCode, 'rmse_objective = prmse.^2'), ...
        'Foundation final candidate objective is not ordinary voltage RMSE.');
    assert(contains(foundationCode, 'map_adjacent_slot_times_local'), ...
        'Foundation does not use the calibrated unequal OPR-slot mapping.');
    assert(contains(foundationCode, 'S.store_bundle_preview_points = inf'), ...
        'Foundation does not retain the full formal waveform bundle for Main10.');
    assert(contains(foundationCode, 'S.vp_gap_ratio_keep = 0;') && ...
        contains(foundationCode, 'S.include_neighbor_eo = false;') && ...
        contains(foundationCode, 'S.phase_safe_expansion = false;'), ...
        'Foundation route is not frozen to strict Top-3/core-only.');
    assert(contains(foundationCode, 'Result.MethodContract = struct('), ...
        'Foundation result does not save MethodContract metadata.');
    resonanceCode = fileread(mainFiles{4});
    assert(contains(resonanceCode, ...
        'strain_time_offset_total_sec - bestTauSec'), ...
        'Strain time-domain evidence is not converted to the BTT clock.');
    strainCode = fileread(mainFiles{12});
    assert(contains(strainCode, ...
        'selectedTimeSec = [min(windowStartSec), max(windowEndSec)]'), ...
        'Strain validation is not tied to the formal identification windows.');
end
if ismember(mode, ["inputs","all"])
    missing = requiredInputs(~cellfun(@isfile, requiredInputs));
    assert(isempty(missing), 'Missing bundled input: %s', strjoin(missing, newline));
end

assert(cfg.window.count == 18, 'Default 20/3/1 settings must produce 18 windows.');
assert(isequal(cfg.frequency.searchHz, [300 1000]), ...
    'Formal search range must be [300 1000] Hz.');
assert(strcmpi(cfg.model.etaPolicy, 'zero') && cfg.model.etaMm == 0, ...
    'Formal eta policy must be eta=0.');
assert(strcmpi(cfg.objective.finalType, 'plain_rmse'), ...
    'Formal final objective must be ordinary voltage RMSE.');
assert(strcmp(cfg.methodFreezeId, '20250527_UNIFIED_STRUCTURED_VOLTAGE_20260811'));
assert(strcmp(cfg.frequency.structureMode, 'single_sync'));
assert(strcmp(cfg.frequency.candidateMethod, 'structured_voltage_vp'));
assert(strcmp(cfg.lowSpeed.templateMethod, 'support_aware_grouped_lap_cv_sg_pchip'));
assert(strcmp(cfg.model.forwardModel, 'low_template_plus_static_gap_increment'));
assert(~cfg.independence.usePreviousWindowCandidate && ...
    ~cfg.independence.useCausalDxState, 'Formal windows must be independent.');

resolved = {which('Config_20250527'), which('CaseConfig'), ...
    which('ProjectionFlow_Config_20250527'), which('BTTDataConfig_20250527')};
assert(all(cellfun(@(p) startsWith(p, cfg.paths.root, 'IgnoreCase', true), resolved)), ...
    'MATLAB resolves one or more formal functions outside this package.');

report = struct('ok', true, 'root', cfg.paths.root, ...
    'formalMainCount', 13, 'requiredInputCount', numel(requiredInputs), ...
    'windowCount', cfg.window.count, 'finalObjective', cfg.objective.finalType, ...
    'etaPolicy', cfg.model.etaPolicy, 'methodFreezeId', cfg.methodFreezeId, ...
    'candidateMethod', cfg.frequency.candidateMethod, ...
    'lowSpeedTemplateMethod', cfg.lowSpeed.templateMethod);
if verbose
    fprintf('20250527 V1 package check passed.\n');
    fprintf('Programs: %d; bundled inputs: %d; windows: %d.\n', ...
        report.formalMainCount, report.requiredInputCount, report.windowCount);
    fprintf('Formal method: adaptive-SG/PCHIP, single-sync structured voltage Top-3, complete voltage optimization, plain RMSE.\n');
end
end
