function detailChunk = run_conf_gap_vp_joint_case_range(startCaseIdx, endCaseIdx)
%run_conf_gap_vp_joint_case_range
% Execute a cached subset of the targeted scenarios for conf_gap_vp_joint.

rootDir = fileparts(fileparts(mfilename('fullpath')));
projectDir = rootDir;
workspaceDir = fileparts(projectDir);
addpath(fullfile(projectDir, 'func'));

oldDir = resolve_legacy_config_dir(workspaceDir);
outDir = fullfile(projectDir, 'results', '04_three_method_compare', '06_final_function_ab', 'conf_gap_vp_joint');
caseOutDir = fullfile(outDir, 'cases');
if ~exist(caseOutDir, 'dir')
    mkdir(caseOutDir);
end

load(fullfile(oldDir, 'stage0_config.mat'), 'cfg');
cfg.dataFile = resolve_project_data_file(workspaceDir, cfg.dataFile);
cfgAna = make_multicase_cfg(cfg);
[scenarioList, snrList] = make_multicase_scenarios_method_difference();
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
templateLib = build_gap_template_library(gapList, xCell, yCell, cfg.g_holdout, cfgAna.xGridN);

idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
if isempty(idxHoldout)
    error('Holdout gap %.6g mm not found in data file.', cfg.g_holdout);
end
F_high_true = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');

totalCaseCount = numel(scenarioList) * numel(snrList);
startCaseIdx = max(1, startCaseIdx);
endCaseIdx = min(totalCaseCount, endCaseIdx);

detailChunk = table();
caseCount = 0;
for is = 1:numel(scenarioList)
    sc = scenarioList(is);
    for in = 1:numel(snrList)
        caseCount = caseCount + 1;
        if caseCount < startCaseIdx || caseCount > endCaseIdx
            continue;
        end

        snrTag = sprintf('%02ddB', round(snrList(in)));
        repeatTag = sprintf('_rep%02d', sc.repeat_id);
        caseName = sprintf('%02d_%s%s_snr_%s', caseCount, char(sc.name), repeatTag, snrTag);
        caseFile = fullfile(caseOutDir, sprintf('%02d_%s%s_snr_%s.mat', ...
            caseCount, char(sc.name), repeatTag, snrTag));

        if isfile(caseFile)
            fprintf('[conf-gap] Reuse cached case: %s, SNR %.1f dB\n', sc.name, snrList(in));
            SCase = load(caseFile, 'caseSummary');
            detailChunk = [detailChunk; SCase.caseSummary]; %#ok<AGROW>
            continue;
        end

        fprintf('[conf-gap] Run case %d/%d: %s, SNR %.1f dB\n', ...
            caseCount, totalCaseCount, sc.name, snrList(in));

        cfgCase = cfgAna;
        cfgCase.noiseMode = 'noise_ratio';
        cfgCase.snrDb = snrList(in);
        cfgCase.A_true = sc.A_true;
        cfgCase.f_true = sc.f_true;
        cfgCase.phi_true = sc.phi_true;
        cfgCase.noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, gapList, xCell, yCell, templateLib.domain, cfgCase.snrDb);

        rng(sc.seed + 1000 * in, 'twister');
        Data_High = simulate_rotating_waveform_from_template(F_high_true, cfgCase.RPM_high, ...
            cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
            cfgCase.noiseRatio, templateLib.domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, cfgCase.noiseMode);
        highMap = map_highspeed_to_space(Data_High, cfgCase.alpha_k, cfgCase.R_tip, ...
            templateLib.domain, cfgCase.fitActiveLevel);
        staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
        staticState.gHat = staticState.gHat + sc.gHat_bias_mm;
        staticState.gHat = min(max(staticState.gHat, min(templateLib.gapTrain) - cfgCase.rawGapSearchMargin), ...
            max(templateLib.gapTrain) + cfgCase.rawGapSearchMargin);

        truthCase = struct('g_high', cfg.g_holdout, 'A', cfgCase.A_true, ...
            'f', cfgCase.f_true, 'phi', cfgCase.phi_true);

        timerRun = tic;
        R = run_confidence_gap_vp_joint(highMap, templateLib, cfgCase, staticState);
        R.elapsed_s = toc(timerRun);

        caseSummary = make_three_method_summary(R, truthCase);
        caseSummary = enrich_conf_summary(caseSummary, cfgCase, truthCase, sc, caseName, caseCount);
        detailChunk = [detailChunk; caseSummary]; %#ok<AGROW>

        save(caseFile, 'R', 'caseSummary', 'cfgCase', 'truthCase', 'staticState', 'Data_High', '-v7.3');
        writetable(detailChunk, fullfile(outDir, sprintf('conf_gap_vp_joint_chunk_%03d_%03d.csv', startCaseIdx, endCaseIdx)));
    end
end

append_detail_chunk(outDir, detailChunk);
end

function append_detail_chunk(outDir, detailChunk)
detailFile = fullfile(outDir, 'conf_gap_vp_joint_detail.csv');
if isempty(detailChunk)
    return;
end
if isfile(detailFile)
    oldTable = readtable(detailFile);
    allTable = [oldTable; detailChunk];
    if ismember('case_name', allTable.Properties.VariableNames) && ismember('method', allTable.Properties.VariableNames)
        [~, ia] = unique(strcat(string(allTable.case_name), "||", string(allTable.method)), 'stable');
        allTable = allTable(sort(ia), :);
    end
else
    allTable = detailChunk;
end
writetable(allTable, detailFile);
end

function T = enrich_conf_summary(T, cfgCase, truthCase, sc, caseName, caseIndex)
nRow = height(T);
[fTrue, orderTrue] = sort(truthCase.f(:).');
ATrue = truthCase.A(orderTrue);

T.group = repmat(string(sc.group), nRow, 1);
T.scenario = repmat(string(sc.name), nRow, 1);
T.case_name = repmat(string(caseName), nRow, 1);
T.case_index = repmat(caseIndex, nRow, 1);
T.snr_db = repmat(cfgCase.snrDb, nRow, 1);
T.noise_mode = repmat(string(cfgCase.noiseMode), nRow, 1);
T.gHat_bias_mm = repmat(sc.gHat_bias_mm, nRow, 1);
T.f1_true_Hz = repmat(fTrue(1), nRow, 1);
T.f2_true_Hz = repmat(fTrue(2), nRow, 1);
T.A1_true_mm = repmat(ATrue(1), nRow, 1);
T.A2_true_mm = repmat(ATrue(2), nRow, 1);
T.f1_error_Hz = abs(T.f1_Hz - T.f1_true_Hz);
T.f2_error_Hz = abs(T.f2_Hz - T.f2_true_Hz);
T.A1_error_mm = abs(T.A1_mm - T.A1_true_mm);
T.A2_error_mm = abs(T.A2_mm - T.A2_true_mm);
T.f1_success = double(T.f1_error_Hz <= 20);
T.f2_success = double(T.f2_error_Hz <= 20);
T.double_success = double(T.f1_success > 0 & T.f2_success > 0);
T.branch_jump = double(T.f1_error_Hz >= 30 | T.f2_error_Hz >= 30);
end

function noiseRatio = calibrate_full_waveform_snr_ratio(cfgCase, gapList, xCell, yCell, domain, snrDb)
idxHoldout = find(abs(gapList - cfgCase.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
DataClean = simulate_rotating_waveform_from_template(Fx, cfgCase.RPM_high, ...
    cfgCase.NumRevs_high, cfgCase.fs, cfgCase.R_tip, cfgCase.alpha_k, ...
    0, domain, cfgCase.A_true, cfgCase.f_true, cfgCase.phi_true, 'noise_ratio');
signal = DataClean.V_clean(:) - min(DataClean.V_clean(:));
signalRms = sqrt(mean(signal.^2));
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatio = signalRms ./ (templateRange .* 10^(snrDb / 20));
end
