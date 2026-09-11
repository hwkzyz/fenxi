function report = Check_Program_20251222(mode, verbose)
%CHECK_PROGRAM_20251222 Verify the independent V1 code and bundled inputs.
if nargin < 1, mode = 'all'; end
if nargin < 2, verbose = true; end
mode = lower(string(mode));
assert(ismember(mode, ["code","inputs","all"]), ...
    'Check mode must be code, inputs, or all.');

cfg = Setup_Paths_20251222();
mainFiles = cell(13,1);
for i = 1:13
    hit = dir(fullfile(cfg.paths.root, sprintf('Main%02d_*_20251222.m', i)));
    assert(numel(hit) == 1, 'Expected exactly one Main%02d program.', i);
    assert(hit.bytes > 500, 'Main%02d is incomplete.', i);
    mainFiles{i} = fullfile(hit.folder, hit.name);
end
requiredCode = [mainFiles; {
    fullfile(cfg.paths.root, 'Config_20251222.m')
    fullfile(cfg.paths.root, 'Setup_Paths_20251222.m')
    fullfile(cfg.paths.utilities, 'CaseConfig.m')
    fullfile(cfg.paths.utilities, 'ProjectionFlow_Config_20251222.m')
    fullfile(cfg.paths.utilities, 'ResonanceRegionCatalog_20251222.m')
    fullfile(cfg.paths.foundation, 'BTTDataConfig_20251222.m')
    fullfile(cfg.paths.preparation, 'Get_20251222_BTT_Config.m')
    fullfile(cfg.paths.gapAware, '+step07jcore', ...
        'solve_gradient_displacement_vp_seed.m')}];
requiredInputs = {
    cfg.files.oprSlotCalibration
    cfg.files.responseSurface
    cfg.files.sharedResponseSurface
    cfg.files.gapLibrary
    fullfile(cfg.paths.calibrationInputs, ...
        'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat')
    fullfile(cfg.paths.calibrationInputs, ...
        'Step06I_OffsetTiltShared_GapLibrary_20251222_B5_S123.mat')
    fullfile(cfg.paths.preparedFoundation, 'step01_low_speed_reference', ...
        'Sensor_Config_20251222.mat')
    fullfile(cfg.paths.preparedFoundation, 'step02_dynamic_btt', ...
        cfg.case.dynamicCase, 'DynamicBTTFeature_20251222.mat')
    fullfile(cfg.paths.preparedFoundation, 'step02_dynamic_btt', ...
        cfg.case.dynamicCase, 'jiluOPR.mat')
    cfg.files.lowSpeedTemplate};

if ismember(mode, ["code","all"])
    missing = requiredCode(~cellfun(@isfile, requiredCode));
    assert(isempty(missing), 'Missing program: %s', strjoin(missing,newline));
    forbidden = {'20251222_low_speed_gap_prior_decoupling', ...
        '20251222_btt_data_foundation', '20251222_low_speed_rotating_calibration'};
    for i = 1:numel(requiredCode)
        source = fileread(requiredCode{i});
        for j = 1:numel(forbidden)
            assert(~contains(source, forbidden{j}), ...
                'External program-folder dependency "%s" in %s.', ...
                forbidden{j}, requiredCode{i});
        end
    end
    fixedCode = fileread(mainFiles{9});
    gapCode = fileread(mainFiles{10});
    assert(contains(fixedCode, 'residual_sse = nansum(res(valid).^2)'), ...
        'Foundation final voltage residual is not ordinary unweighted SSE.');
    assert(contains(fixedCode, 'S.store_bundle_preview_points = inf'), ...
        'Foundation does not retain the complete bundle for Main10.');
    assert(contains(fixedCode, 'S.vp_top_k_max = 3;') && ...
        contains(fixedCode, 'S.vp_top_rank_keep = 3;') && ...
        contains(fixedCode, 'S.vp_gap_ratio_keep = 0;') && ...
        contains(fixedCode, 'S.include_neighbor_eo = false;'), ...
        'Foundation route is not frozen to strict Top-3.');
    assert(contains(fixedCode, 'Result.MethodContract = struct('), ...
        'Foundation result does not save MethodContract metadata.');
    assert(contains(gapCode, 'cfg.vpTopK = 3') && ...
        contains(gapCode, 'cfg.vpNeighborRadius = 0'), ...
        'GapAware must use its own strict Top-3 without neighboring EO.');
    assert(contains(gapCode, "cfg.dxReferenceMode = 'seed_median'"), ...
        'GapAware windows do not use the independent seed-median dx reference.');
    assert(contains(gapCode, 'residualObj = sum(res.^2)'), ...
        'GapAware final voltage residual is not ordinary unweighted SSE.');
    assert(contains(gapCode, 'Result.MethodContract = struct('), ...
        'GapAware result does not save MethodContract metadata.');
end
if ismember(mode, ["inputs","all"])
    missing = requiredInputs(~cellfun(@isfile, requiredInputs));
    assert(isempty(missing), 'Missing bundled input: %s', strjoin(missing,newline));
end

assert(cfg.window.count == 18, 'Default 20/3/1 must produce 18 windows.');
assert(strcmp(cfg.methodFreezeId, '20251222_V1_FROZEN_20260717') && ...
    strcmp(cfg.methodStatus, 'frozen_formal_route'), ...
    'The formal method freeze identifier has changed.');
assert(isequal(cfg.frequency.searchHz,[300 1000]) && cfg.frequency.vpTopK == 3);
assert(strcmp(cfg.frequency.foundationCandidatePolicy, 'adaptive_full_fit'));
assert(strcmp(cfg.frequency.gapCandidatePolicy, 'gapaware_vp_top3_only'));
assert(strcmpi(cfg.model.etaPolicy,'zero') && cfg.model.etaMm == 0);
assert(strcmpi(cfg.objective.finalType,'plain_rmse'));
assert(strcmpi(cfg.objective.finalResidualWeighting,'none'));
assert(strcmpi(cfg.model.mainGapModel,'gap_only'));
assert(~cfg.independence.usePreviousWindowCandidate && ...
    ~cfg.independence.useCausalDxState, 'Windows must be independent.');
assert(strcmp(cfg.independence.dxReferenceMode,'seed_median'));
assert(strcmp(cfg.bundle.foundationRoleInGapAware, ...
    'common_raw_waveform_bundle_only') && ...
    ~cfg.bundle.inheritFoundationEOAmplitudePhaseDx, ...
    'GapAware must not inherit fitted Foundation parameters.');
assert(abs(sum(cfg.opr.slotAngleDeg)-360) < 1e-8, ...
    'Unequal OPR slot angles must sum to 360 degrees.');

resolved = {which('Config_20251222'), which('CaseConfig'), ...
    which('ProjectionFlow_Config_20251222'), which('BTTDataConfig_20251222')};
assert(all(cellfun(@(p) startsWith(p,cfg.paths.root,'IgnoreCase',true),resolved)), ...
    'MATLAB resolves one or more formal functions outside this package.');
report = struct('ok',true,'root',cfg.paths.root,'formalMainCount',13, ...
    'requiredInputCount',numel(requiredInputs),'windowCount',cfg.window.count, ...
    'finalObjective',cfg.objective.finalType,'etaPolicy',cfg.model.etaPolicy, ...
    'methodFreezeId',cfg.methodFreezeId);
if verbose
    fprintf('20251222 V1 package check passed.\n');
    fprintf('Programs: %d; bundled inputs: %d; windows: %d.\n', ...
        report.formalMainCount, report.requiredInputCount, report.windowCount);
    fprintf(['Method: unequal OPR slots, eta=0, independent 3-lap windows, ' ...
        'strict method-specific Top-3, plain final voltage residual.\n']);
end
end
