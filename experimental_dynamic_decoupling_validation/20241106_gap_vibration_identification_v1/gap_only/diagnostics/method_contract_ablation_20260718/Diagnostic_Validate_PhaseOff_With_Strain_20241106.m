function Diagnostic_Validate_PhaseOff_With_Strain_20241106(selectedCase)
%% Step12: validate Step07J against strain-gauge spectrum and waveform, 20241106
% Main figure logic:
%   a) AI1-04 strain spectrum around the selected 12EO resonance band;
%   b) S257 Step07J frequency error relative to the strain FFT peak;
%   c) TRAC trend using AI1-04, with all AI1 channels as sensitivity checks;
%   d) representative-window waveform comparison.
% Strain frequency is identified independently: the complete matched strain
% span supplies one reference peak, then every exact BTT window is fitted in
% a narrow band around that strain-derived reference.
%
% The 20241106 case differs from 20250527: Step04 stores the selected
% AI1-04 strain trace directly instead of a strainSource struct. Therefore
% this script uses the Step04-embedded AI1-04 trace as the primary evidence
% and reads sibling AI1 files only for channel-sensitivity panels.
%
% Strain-to-tip displacement scaling used by this final 20241106 route:
%   FE model: Workbench blade with the real curved root, 632.388232542 Hz;
%   gauge point: Y=4.00 mm from the root and Z=2.80 mm from the edge in
%                Workbench coordinates (X=0 and 1.4 mm are the two faces);
%   response pair: longitudinal strain EPELY / tip displacement UX;
%   surface interpolation: natural interpolation on each physical face;
%   R_front=1.1308777609 1/m, R_back=1.1308564487 1/m;
%   final mean R=1.1308671048 1/m.

clc; close all;

%% 1. Settings
% Diagnostic copy only: validate the phase-off result without touching formal outputs.
thisDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(thisDir, fullfile(thisDir, 'functions', 'foundation'), ...
    fullfile(thisDir, 'functions', 'preparation'), ...
    fullfile(thisDir, 'functions', 'gap_aware'), ...
    fullfile(thisDir, 'functions', 'utilities'));
packageCfg = Config_20241106();
C0 = CaseConfig();
validationSourceDir = packageCfg.paths.strainValidation;
outDir = fullfile(thisDir, 'diagnostics', 'method_contract_ablation_20260718', ...
    'results_strain_validation_phaseoff');
resultDir = packageCfg.paths.gapResults;
figDir = fullfile(outDir, 'figures_step12_trac_strain_spectrum');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

caseCatalog = step12_case_catalog_local(packageCfg);
if nargin < 1
    selectedCase = [];
end
caseSpec = select_step12_case_local(caseCatalog, selectedCase);
if isempty(caseSpec)
    fprintf('Step12 case selection cancelled.\n');
    return;
end
C0.caseTag = sprintf('%s_%s', C0.caseTag, caseSpec.timeTag);

methodName = 'gap_only';
methodLabel = 'S257 gap only';
timeShiftSearchSec = 0;
feModelLabel = '632.388 Hz Workbench curved-root blade';
feModelFrequencyHz = 632.3882325418;
feGaugeLocation = 'Workbench Y=4.00 mm, Z=2.80 mm; faces X=0/1.4 mm';
feResponseDefinition = 'abs(EPELY / UX_tip_m), natural interpolation per face';
feTransferRatioFront_1PerM = 1.1308777609;
feTransferRatioBack_1PerM = 1.13085644867;
feTransferRatio_1PerM = mean([feTransferRatioFront_1PerM, ...
    feTransferRatioBack_1PerM]);
feSourceFile = ['D:\博士-国科\试验台数据\大台子\叶片模型\580Hz传递比分析\' ...
    'workbench_direct_export\workbench_direct_nodes_632Hz.csv'];
K_Strain2MM = 1 / feTransferRatio_1PerM / 1000;
peakSearchHz = [600, 650];
localFrequencyHalfBandHz = 3.0;
localFrequencyStepHz = 0.02;
expectedStrainToBttOffsetSec = packageCfg.strain.timeOffsetPriorSec;
globalAlignmentFile = fullfile(validationSourceDir, sprintf( ...
    'Step12_GlobalTimeAlignment_20241106_B%d_S%s.mat', ...
    packageCfg.case.targetBlade, sprintf('%d', packageCfg.case.analysisSensors)));
step07jResultName = ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B4_S257_diag_contract_phaseoff_formalfoundation.mat';

%% 2. Load Step07J and Step04 evidence
step07jFile = fullfile(resultDir, step07jResultName);
step04File = fullfile(packageCfg.paths.strainEvidence, ...
    'Step04_BTT_STE_Resonance_Regions_20241106.mat');
if ~isfile(step07jFile)
    error('Missing Step07J result: %s', step07jFile);
end
if ~isfile(step04File)
    error('Run Step04 first. Missing file: %s', step04File);
end

S7 = load(step07jFile, 'Result');
R7 = S7.Result;
windowResult = R7.WindowResult(:);
validate_current_step07j_result_local(R7, C0, methodName, step07jFile, ...
    caseSpec.analysisStartTimeSec);
fprintf('Selected Step12 case: %s\n', caseSpec.displayName);
fprintf('Using Step07J result: %s\n', step07jFile);

S4 = load(step04File, 'regionTable', 'localFftTable', 'strain', ...
    'strainToBttOffsetSec', 'selectedStrainChannel', 'channelSummary');
if isempty(S4.regionTable)
    error('Step04 selected region table is empty.');
end
if ~isfinite(S4.strainToBttOffsetSec) || ...
        abs(S4.strainToBttOffsetSec - expectedStrainToBttOffsetSec) > 1e-9
    error('Unexpected strainToBttOffsetSec: %.9f s.', S4.strainToBttOffsetSec);
end
globalAlignment = load_global_alignment_local(globalAlignmentFile, ...
    expectedStrainToBttOffsetSec);
effectiveStrainToBttOffsetSec = globalAlignment.strainToBttOffsetSec;
windowBounds = vertcat(windowResult.timeWindow);
selectedTimeSec = [min(windowBounds(:, 1)), max(windowBounds(:, 2))];
regionIndex = find(selectedTimeSec(1) >= S4.regionTable.bttStartSec & ...
    selectedTimeSec(2) <= S4.regionTable.bttEndSec, 1, 'first');
if isempty(regionIndex)
    error('Step07J time span %.6f-%.6f s is outside the Step04 resonance regions.', ...
        selectedTimeSec(1), selectedTimeSec(2));
end
windowOrderFreqHz = arrayfun(@(w) ...
    w.modelFits.(methodName).EO * w.rotRpmMean / 60, windowResult);
strainOrderFreqHz = mean(windowOrderFreqHz, 'omitnan');

primarySource = make_primary_strain_source_local(S4, ...
    effectiveStrainToBttOffsetSec);
[strainTime, strainMm, strainFile] = load_aligned_strain_mm_local(primarySource, K_Strain2MM);
[strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, strainFftPeakAmpMm] = ...
    compute_selected_strain_fft_mm_local(strainTime, strainMm, selectedTimeSec, peakSearchHz);
strainReferenceFit = fit_strain_window_near_reference_local(strainTime, strainMm, ...
    selectedTimeSec, strainFftPeakHz, localFrequencyHalfBandHz, ...
    localFrequencyStepHz);
if strainReferenceFit.hitFrequencyBoundary
    error('Global strain-frequency refinement reached its search boundary.');
end
strainReferenceHz = strainReferenceFit.frequencyHz;
primaryChannel = primarySource.channel;

strainSources = make_strain_source_options_local(primarySource);

%% 3. Compute primary TRAC table
rows = {};
detailData = {};
detailPickWindow = unique([1, round(numel(windowResult) / 2), R7.BestWindowIndex]);
for w = 1:numel(windowResult)
    tWin = windowResult(w).timeWindow;
    if tWin(2) < selectedTimeSec(1) || tWin(1) > selectedTimeSec(2)
        continue;
    end
    fit = windowResult(w).modelFits.(methodName);
    strainFit = fit_strain_window_near_reference_local(strainTime, strainMm, ...
        tWin, strainReferenceHz, localFrequencyHalfBandHz, localFrequencyStepHz);
    [rawTrac, rawPhase, tRaw, yRaw, yReconRaw] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, 0, fit.amplitudeMm, fit.freqHz);
    [bestTrac, bestPhase, bestShift, tBest, yBest, yReconBest] = ...
        search_best_trac_shift_local(strainTime, strainMm, tWin, ...
        timeShiftSearchSec, fit.amplitudeMm, fit.freqHz);

    rows(end + 1, :) = {windowResult(w).windowId, tWin(1), tWin(2), ...
        tWin(1) - effectiveStrainToBttOffsetSec, ...
        tWin(2) - effectiveStrainToBttOffsetSec, ...
        methodLabel, fit.EO, fit.freqHz, fit.amplitudeMm, fit.weightedRmseMv, ...
        primaryChannel, rawPhase, rawTrac, bestPhase, bestTrac, bestShift, ...
        numel(tRaw), rms(yRaw, 'omitnan'), strainFftPeakHz, ...
        strainFftPeakAmpMm, strainReferenceHz, fit.freqHz - strainReferenceHz, ...
        strainFit.frequencyHz, strainFit.amplitudeMm, ...
        fit.freqHz - strainFit.frequencyHz, ...
        fit.amplitudeMm / strainFit.amplitudeMm, strainFit.rSquared, ...
        strainFit.rmseMm, strainFit.hitFrequencyBoundary}; %#ok<AGROW>

    if ismember(windowResult(w).windowId, detailPickWindow)
        s = struct();
        s.windowId = windowResult(w).windowId;
        s.timeRaw = tRaw;
        s.strainRaw = yRaw;
        s.reconRaw = yReconRaw;
        s.rawTrac = rawTrac;
        s.timeBest = tBest;
        s.strainBest = yBest;
        s.reconBest = yReconBest;
        s.bestTrac = bestTrac;
        s.bestShift = bestShift;
        s.freqHz = fit.freqHz;
        s.ampMm = fit.amplitudeMm;
        detailData{end + 1} = s; %#ok<SAGROW>
    end
end

tracTable = cell2table(rows, 'VariableNames', {'window_id', ...
    'btt_time_start_s', 'btt_time_end_s', 'raw_strain_time_start_s', ...
    'raw_strain_time_end_s', 'method', 'EO', 'frequency_hz', ...
    'amplitude_mm', 'rmse_mV', 'strain_channel', 'raw_phase_rad', ...
    'raw_TRAC', 'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm', 'strain_fft_peak_hz', ...
    'strain_fft_peak_amp_mm', 'strain_global_reference_hz', ...
    'freq_error_vs_strain_global_hz', ...
    'strain_local_frequency_hz', 'strain_local_amplitude_mm', ...
    'btt_minus_strain_local_hz', 'btt_to_strain_amp_ratio', ...
    'strain_fit_r_squared', 'strain_fit_rmse_mm', ...
    'strain_fit_hit_frequency_boundary'});
tracTable.fe_model = repmat(string(feModelLabel), height(tracTable), 1);
tracTable.fe_model_frequency_hz = repmat(feModelFrequencyHz, height(tracTable), 1);
tracTable.fe_gauge_location = repmat(string(feGaugeLocation), height(tracTable), 1);
tracTable.fe_response_definition = repmat(string(feResponseDefinition), height(tracTable), 1);
tracTable.fe_transfer_ratio_1_per_m = repmat(feTransferRatio_1PerM, height(tracTable), 1);
tracTable.strain_to_tip_mm_per_microstrain = repmat(K_Strain2MM, height(tracTable), 1);
tracTable.fe_source_file = repmat(string(feSourceFile), height(tracTable), 1);
tracTable.time_alignment_offset_s = repmat( ...
    effectiveStrainToBttOffsetSec, height(tracTable), 1);
tracTable.time_alignment_candidate_offset_s = repmat( ...
    globalAlignment.estimatedOffsetSec, height(tracTable), 1);
tracTable.time_alignment_accepted = repmat( ...
    globalAlignment.accepted, height(tracTable), 1);
tracTable.time_alignment_decision = repmat( ...
    string(globalAlignment.decision), height(tracTable), 1);
tracTable.case_analysis_start_time_s = repmat( ...
    caseSpec.analysisStartTimeSec, height(tracTable), 1);
tracTable.case_time_tag = repmat(string(caseSpec.timeTag), height(tracTable), 1);
tracTable.step07j_result_file = repmat(string(step07jFile), height(tracTable), 1);
if any(tracTable.strain_fit_hit_frequency_boundary)
    warning('%d strain windows reached the local frequency-search boundary.', ...
        nnz(tracTable.strain_fit_hit_frequency_boundary));
end

csvFile = fullfile(outDir, sprintf('Step12_TRAC_StrainSpectrum_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(tracTable, csvFile);

channelTable = compute_channel_sensitivity_table_local(windowResult, methodName, ...
    methodLabel, strainSources, timeShiftSearchSec, selectedTimeSec, K_Strain2MM);
channelCsvFile = fullfile(outDir, sprintf( ...
    'Step12_TRAC_ChannelSensitivity_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(channelTable, channelCsvFile);

auditTable = build_step07j_order_audit_local(outDir, C0, step07jFile, ...
    methodName, strainReferenceHz, strainOrderFreqHz);
auditCsvFile = fullfile(outDir, sprintf( ...
    'Step12_OrderAudit_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(auditTable, auditCsvFile);

%% 4. Publication-style main figure
natureFigFile = plot_step12_nature_validation_local(tracTable, channelTable, ...
    detailData, strainFftFreqHz, strainFftAmpMm, strainFftPeakHz, ...
    strainReferenceHz, strainOrderFreqHz, figDir, C0, primaryChannel, selectedTimeSec);
auditFigFile = plot_step12_order_audit_local(auditTable, figDir, C0, ...
    strainReferenceHz, strainOrderFreqHz);
[channelTimeFigFile, channelTimeTable] = plot_channel_time_domain_compare_local( ...
    windowResult, methodName, strainSources, channelTable, figDir, C0, ...
    K_Strain2MM, R7.BestWindowIndex);
channelTimeCsvFile = fullfile(outDir, sprintf( ...
    'Step12_ChannelTimeDomainCompare_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(channelTimeTable, channelTimeCsvFile);
[rawStrainFigFile, rawStrainTable, rawStrainWindowSec] = ...
    plot_raw_strain_selection_local(S4.strain, selectedTimeSec, ...
    effectiveStrainToBttOffsetSec, figDir, C0);
rawStrainCsvFile = fullfile(outDir, sprintf( ...
    'Step12_RawStrainSelection_%s_%s.csv', C0.dataset, C0.caseTag));
writetable(rawStrainTable, rawStrainCsvFile);

fprintf('Step12 20241106 validation complete.\n');
fprintf('Primary strain channel: %s (%s)\n', primaryChannel, strainFile);
fprintf(['FE transfer: %s | %.9f Hz | R=%.10f 1/m | ' ...
    'K=%.10g mm/microstrain\n'], feModelLabel, feModelFrequencyHz, ...
    feTransferRatio_1PerM, K_Strain2MM);
fprintf('FE gauge/response: %s | %s\n', feGaugeLocation, feResponseDefinition);
fprintf('FE source: %s\n', feSourceFile);
fprintf(['Global time alignment: applied %.9f s, candidate %.9f s, ' ...
    'accepted=%d (%s).\n'], effectiveStrainToBttOffsetSec, ...
    globalAlignment.estimatedOffsetSec, globalAlignment.accepted, ...
    globalAlignment.decision);
fprintf('Selected time %.3f-%.3f s. Strain FFT peak %.6f Hz, 12EO line %.6f Hz.\n', ...
    selectedTimeSec(1), selectedTimeSec(2), strainFftPeakHz, strainOrderFreqHz);
fprintf('Selected raw strain time %.6f-%.6f s (%d original samples).\n', ...
    rawStrainWindowSec(1), rawStrainWindowSec(2), height(rawStrainTable));
fprintf(['Independent strain fit: global reference %.6f Hz; local search ' ...
    '+/-%.2f Hz with %.3f Hz step.\n'], strainReferenceHz, ...
    localFrequencyHalfBandHz, localFrequencyStepHz);
disp(groupsummary(tracTable, 'strain_channel', 'mean', ...
    {'raw_TRAC', 'best_TRAC', 'best_time_shift_s', ...
    'btt_minus_strain_local_hz', 'btt_to_strain_amp_ratio', 'rmse_mV'}));
amplitudeErrorMm = tracTable.amplitude_mm - tracTable.strain_local_amplitude_mm;
fprintf(['Amplitude comparison: mean BTT %.6f mm, mean strain %.6f mm, ' ...
    'mean bias %.6f mm, MAE %.6f mm, MAPE %.3f%%.\n'], ...
    mean(tracTable.amplitude_mm, 'omitnan'), ...
    mean(tracTable.strain_local_amplitude_mm, 'omitnan'), ...
    mean(amplitudeErrorMm, 'omitnan'), ...
    mean(abs(amplitudeErrorMm), 'omitnan'), ...
    100 * mean(abs(amplitudeErrorMm) ./ tracTable.strain_local_amplitude_mm, 'omitnan'));
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    csvFile, channelCsvFile, auditCsvFile, channelTimeCsvFile, ...
    rawStrainCsvFile, natureFigFile, auditFigFile, channelTimeFigFile, ...
    rawStrainFigFile);

end

%% Helpers
function cases = step12_case_catalog_local(cfg)
times = unique([cfg.strain.alignmentCasesSec(:); cfg.case.analysisStartTimeSec], 'stable');
cases = repmat(struct('analysisStartTimeSec', NaN, 'timeTag', '', ...
    'displayName', '', 'resultFileName', ''), numel(times), 1);
for i = 1:numel(times)
    cases(i).analysisStartTimeSec = times(i);
    cases(i).timeTag = strrep(sprintf('T%07.3f', times(i)), '.', 'p');
    cases(i).displayName = sprintf('%.3f s resonance segment (%s)', ...
        times(i), cases(i).timeTag);
    cases(i).resultFileName = formal_gap_result_name_local(cfg, times(i));
end
end

function name = formal_gap_result_name_local(cfg, analysisTimeSec)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', analysisTimeSec), '.', 'p');
suffix = sprintf('_%s_main_gapfixedtilt_W%dS%d_F%dto%dHz_DFpm%gHz', ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)), ...
    cfg.frequency.refineHalfWidthHz);
name = sprintf('Step07J_NestedStaticWarp_VPFullWave_20241106_B%d_%s%s.mat', ...
    cfg.case.targetBlade, sensorTag, suffix);
end

function selected = select_step12_case_local(cases, requested)
selected = struct([]);
if isempty(requested)
    if ~usejava('desktop')
        error(['No Step12 case was specified. Call the function with 75.0 or 86.6, ' ...
            'for example Main12_Validate_Identification_With_Strain_20241106(86.6).']);
    end
    [index, accepted] = listdlg('PromptString', 'Select the Step12 result to display:', ...
        'SelectionMode', 'single', 'ListString', {cases.displayName}, ...
        'Name', 'Step12 case selection', 'ListSize', [360, 120]);
    if ~accepted
        return;
    end
    selected = cases(index);
    return;
end

if isnumeric(requested) && isscalar(requested) && isfinite(requested)
    index = find(abs([cases.analysisStartTimeSec] - requested) < 1e-9, 1);
elseif ischar(requested) || (isstring(requested) && isscalar(requested))
    token = char(string(requested));
    index = find(strcmpi({cases.timeTag}, token) | ...
        strcmpi({cases.displayName}, token), 1);
else
    index = [];
end
if isempty(index)
    available = strjoin(arrayfun(@(c) sprintf('%.1f (%s)', ...
        c.analysisStartTimeSec, c.timeTag), cases, 'UniformOutput', false), ', ');
    error('Unknown Step12 case. Available cases: %s.', available);
end
selected = cases(index);
end

function [strainTime, strainMm, strainFile] = load_aligned_strain_mm_local(strainSource, strainToMm)
strainFile = strainSource.strain_file;
if isfield(strainSource, 'time') && isfield(strainSource, 'value')
    strainTime = strainSource.time(:) + strainSource.strain_time_offset_total_sec;
    strainMm = (strainSource.value(:) - mean(strainSource.value(:), 'omitnan')) * strainToMm;
else
    [strainTime, strainMm] = load_strain_file_mm_local(strainFile, ...
        strainSource.strain_time_offset_total_sec, strainToMm);
    strainFile = strainSource.strain_file;
end

valid = isfinite(strainTime) & isfinite(strainMm);
strainTime = strainTime(valid);
strainMm = strainMm(valid);
[strainTime, orderIdx] = sort(strainTime);
strainMm = strainMm(orderIdx);
end

function [pngFile, selectionTable, rawWindowSec] = ...
    plot_raw_strain_selection_local(strain, bttWindowSec, offsetSec, figDir, C0)
requiredFields = {'time', 'rawValue', 'value'};
if ~all(isfield(strain, requiredFields))
    error('Step04 strain data lacks time, rawValue, or value for the raw-selection audit.');
end

rawTime = strain.time(:);
rawValue = strain.rawValue(:);
analysisValue = strain.value(:);
rawWindowSec = bttWindowSec - offsetSec;
selected = rawTime >= rawWindowSec(1) & rawTime <= rawWindowSec(2) & ...
    isfinite(rawTime) & isfinite(rawValue) & isfinite(analysisValue);
if nnz(selected) < 20
    error('Raw strain selection %.6f-%.6f s has fewer than 20 valid samples.', ...
        rawWindowSec(1), rawWindowSec(2));
end

selectionTable = table(rawTime(selected), rawTime(selected) + offsetSec, ...
    rawValue(selected), analysisValue(selected), ...
    'VariableNames', {'raw_strain_time_s', 'btt_aligned_time_s', ...
    'raw_strain_microstrain', 'analysis_detrended_strain_microstrain'});

baseName = sprintf('Step12_RawStrainSelection_%s_%s', ...
    C0.dataset, C0.caseTag);
pngFile = fullfile(figDir, [baseName, '.png']);
fig = figure('Name', 'Step12 raw strain selection audit', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 12.5], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
finiteFull = isfinite(rawTime) & isfinite(rawValue);
fullIndex = find(finiteFull);
stride = max(1, ceil(numel(fullIndex) / 50000));
fullIndex = fullIndex(1:stride:end);
plot(ax, rawTime(fullIndex), rawValue(fullIndex), '-', ...
    'Color', [0.30, 0.30, 0.30], 'LineWidth', 0.35);
hold(ax, 'on');
robustLimits = prctile(rawValue(finiteFull), [0.01, 99.99]);
if ~all(isfinite(robustLimits)) || robustLimits(2) <= robustLimits(1)
    robustLimits = [min(rawValue(finiteFull)), max(rawValue(finiteFull))];
end
padding = 0.05 * max(diff(robustLimits), eps);
robustLimits = robustLimits + [-padding, padding];
ylim(ax, robustLimits);
selectionPatch = patch(ax, ...
    [rawWindowSec(1), rawWindowSec(2), rawWindowSec(2), rawWindowSec(1)], ...
    [robustLimits(1), robustLimits(1), robustLimits(2), robustLimits(2)], ...
    [0.85, 0.33, 0.10], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
uistack(selectionPatch, 'bottom');
xline(ax, rawWindowSec(1), '--', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
xline(ax, rawWindowSec(2), '--', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
xlabel(ax, 'Original strain acquisition time (s)');
ylabel(ax, 'Raw strain (microstrain)');
title(ax, sprintf('Complete AI1-04 raw strain; selected %.6f-%.6f s', ...
    rawWindowSec(1), rawWindowSec(2)));
apply_nature_axes_local(ax);
panel_label_local(ax, 'a');

ax = nexttile;
plot(ax, selectionTable.raw_strain_time_s, ...
    selectionTable.raw_strain_microstrain, '-', ...
    'Color', [0.10, 0.35, 0.70], 'LineWidth', 0.65);
xline(ax, rawWindowSec(1), '--', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
xline(ax, rawWindowSec(2), '--', 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
xlim(ax, rawWindowSec);
xlabel(ax, 'Original strain acquisition time (s)');
ylabel(ax, 'Raw strain (microstrain)');
title(ax, sprintf(['Selected original waveform; BTT %.6f-%.6f s, ' ...
    'offset %+.6f s'], bttWindowSec(1), bttWindowSec(2), offsetSec));
apply_nature_axes_local(ax);
panel_label_local(ax, 'b');

export_nature_figure_local(fig, fullfile(figDir, baseName));
% close(fig);
end

function primarySource = make_primary_strain_source_local(S4, effectiveOffsetSec)
primarySource = struct();
primarySource.strain_file = S4.strain.file;
primarySource.channel = string(S4.strain.label);
primarySource.strain_time_offset_total_sec = effectiveOffsetSec;
primarySource.time = S4.strain.time(:);
if isfield(S4.strain, 'value')
    primarySource.value = S4.strain.value(:);
else
    primarySource.value = S4.strain.rawValue(:);
end
if strlength(primarySource.channel) == 0
    primarySource.channel = sprintf('AI1-%02d', S4.selectedStrainChannel);
end
primarySource.channel = char(primarySource.channel);
end

function Alignment = load_global_alignment_local(alignmentFile, priorOffsetSec)
if ~isfile(alignmentFile)
    warning('Global alignment file is missing; retain prior offset %.9f s.', ...
        priorOffsetSec);
    Alignment = struct('strainToBttOffsetSec', priorOffsetSec, ...
        'estimatedOffsetSec', priorOffsetSec, 'accepted', false, ...
        'decision', 'alignment_file_missing_retain_prior');
    return;
end
S = load(alignmentFile, 'Alignment');
if ~isfield(S, 'Alignment') || ...
        ~isfield(S.Alignment, 'strainToBttOffsetSec') || ...
        ~isfinite(S.Alignment.strainToBttOffsetSec)
    error('Invalid global alignment file: %s', alignmentFile);
end
Alignment = S.Alignment;
end

function [freqHz, ampMm, peakFreqHz, peakAmpMm] = compute_selected_strain_fft_mm_local( ...
    strainTime, strainMm, selectedTimeSec, peakSearchHz)
mask = strainTime >= selectedTimeSec(1) & strainTime <= selectedTimeSec(2);
tSeg = strainTime(mask);
ySeg = strainMm(mask);
if numel(tSeg) < 16
    error('Selected strain segment %.3f-%.3f s has too few samples.', ...
        selectedTimeSec(1), selectedTimeSec(2));
end
ySeg = detrend(ySeg(:) - mean(ySeg, 'omitnan'));
dt = median(diff(tSeg), 'omitnan');
fs = 1 / dt;
n = numel(ySeg);
win = hann_window_local(n);
yWin = ySeg .* win;
nfft = 2^nextpow2(n);
Y = fft(yWin, nfft);
ampMm = abs(Y(1:nfft/2 + 1)) / max(sum(win) / 2, eps);
freqHz = (0:nfft/2).' * fs / nfft;

peakMask = freqHz >= peakSearchHz(1) & freqHz <= peakSearchHz(2);
freqBand = freqHz(peakMask);
ampBand = ampMm(peakMask);
[peakAmpMm, peakIdx] = max(ampBand);
peakFreqHz = freqBand(peakIdx);
end

function fit = fit_strain_window_near_reference_local(strainTime, strainMm, ...
    timeWindow, referenceFrequencyHz, halfBandHz, frequencyStepHz)
mask = strainTime >= timeWindow(1) & strainTime <= timeWindow(2);
t = strainTime(mask);
y = strainMm(mask);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
if numel(t) < 20
    error('Strain window %.6f-%.6f s has fewer than 20 valid samples.', ...
        timeWindow(1), timeWindow(2));
end

tLocal = t(:) - mean(t, 'omitnan');
y = detrend(y(:));
frequencyGrid = (referenceFrequencyHz - halfBandHz):frequencyStepHz: ...
    (referenceFrequencyHz + halfBandHz);
if isempty(frequencyGrid)
    error('The local strain-frequency grid is empty.');
end

sse = inf(size(frequencyGrid));
coefficients = nan(numel(frequencyGrid), 3);
for iFrequency = 1:numel(frequencyGrid)
    omegaT = 2 * pi * frequencyGrid(iFrequency) .* tLocal;
    design = [cos(omegaT), sin(omegaT), ones(size(tLocal))];
    coef = design \ y;
    residual = y - design * coef;
    sse(iFrequency) = sum(residual .^ 2);
    coefficients(iFrequency, :) = coef(:).';
end

[bestSse, bestIndex] = min(sse);
bestCoef = coefficients(bestIndex, :);
totalSse = sum((y - mean(y, 'omitnan')) .^ 2);
fit = struct();
fit.frequencyHz = frequencyGrid(bestIndex);
fit.amplitudeMm = hypot(bestCoef(1), bestCoef(2));
fit.phaseRad = atan2(-bestCoef(2), bestCoef(1));
fit.offsetMm = bestCoef(3);
fit.rmseMm = sqrt(bestSse / numel(y));
fit.rSquared = 1 - bestSse / max(totalSse, eps);
fit.sampleCount = numel(y);
fit.hitFrequencyBoundary = bestIndex == 1 || bestIndex == numel(frequencyGrid);
fit.searchRangeHz = [frequencyGrid(1), frequencyGrid(end)];
end

function win = hann_window_local(n)
if n <= 1
    win = ones(n, 1);
else
    k = (0:n-1).';
    win = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
end
end

function [bestTrac, bestPhase, bestShift, bestTime, bestY, bestRecon] = ...
    search_best_trac_shift_local(strainTime, strainMm, tWin, shiftList, ampMm, freqHz)
bestTrac = -inf;
bestPhase = NaN;
bestShift = NaN;
bestTime = [];
bestY = [];
bestRecon = [];
for shiftNow = shiftList
    [tracNow, phaseNow, tNow, yNow, reconNow] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, shiftNow, ampMm, freqHz);
    if isfinite(tracNow) && tracNow > bestTrac
        bestTrac = tracNow;
        bestPhase = phaseNow;
        bestShift = shiftNow;
        bestTime = tNow;
        bestY = yNow;
        bestRecon = reconNow;
    end
end
if isinf(bestTrac)
    bestTrac = NaN;
end
end

function [tracValue, phaseRad, tSeg, ySeg, yRecon] = calc_trac_at_shift_local( ...
    strainTime, strainMm, tWin, shiftSec, ampMm, freqHz)
mask = strainTime >= (tWin(1) + shiftSec) & strainTime <= (tWin(2) + shiftSec);
tSeg = strainTime(mask);
ySeg = strainMm(mask);
if numel(tSeg) < 16
    tracValue = NaN;
    phaseRad = NaN;
    yRecon = [];
    return;
end
ySeg = detrend(ySeg(:) - mean(ySeg, 'omitnan'));
tSeg = tSeg(:);

sinBase = sin(2 * pi * freqHz * tSeg);
cosBase = cos(2 * pi * freqHz * tSeg);
a = ySeg' * sinBase;
b = ySeg' * cosBase;
phaseRad = atan2(b, a);
yRecon = ampMm * sin(2 * pi * freqHz * tSeg + phaseRad);

num = (ySeg' * yRecon)^2;
den = (ySeg' * ySeg) * (yRecon' * yRecon);
if den <= 0 || ~isfinite(den)
    tracValue = NaN;
else
    tracValue = num / den;
end
end

function strainSources = make_strain_source_options_local(primarySource)
strainDir = fileparts(primarySource.strain_file);
files = dir(fullfile(strainDir, 'AI1-*.mat'));
strainSources = struct('channel', {}, 'file', {}, 'offsetSec', {});
for i = 1:numel(files)
    s = struct();
    s.file = fullfile(files(i).folder, files(i).name);
    s.channel = channel_label_from_file_local(s.file);
    s.offsetSec = primarySource.strain_time_offset_total_sec;
    strainSources(end + 1) = s; %#ok<AGROW>
end
if isempty(strainSources)
    s = struct();
    s.file = primarySource.strain_file;
    s.channel = primarySource.channel;
    s.offsetSec = primarySource.strain_time_offset_total_sec;
    strainSources = s;
    warning('No sibling AI1-*.mat strain channel files found in %s; using only %s.', ...
        strainDir, primarySource.channel);
end
end

function T = compute_channel_sensitivity_table_local(windowResult, methodName, ...
    methodLabel, strainSources, shiftList, selectedTimeSec, strainToMm)
rows = {};
for si = 1:numel(strainSources)
    [strainTime, strainMm] = load_strain_file_mm_local( ...
        strainSources(si).file, strainSources(si).offsetSec, strainToMm);
    for w = 1:numel(windowResult)
        tWin = windowResult(w).timeWindow;
        if tWin(2) < selectedTimeSec(1) || tWin(1) > selectedTimeSec(2)
            continue;
        end
        fit = windowResult(w).modelFits.(methodName);
        [rawTrac, rawPhase, tRaw, yRaw] = calc_trac_at_shift_local( ...
            strainTime, strainMm, tWin, 0, fit.amplitudeMm, fit.freqHz);
        [bestTrac, bestPhase, bestShift] = search_best_trac_shift_local( ...
            strainTime, strainMm, tWin, shiftList, fit.amplitudeMm, fit.freqHz);
        rows(end + 1, :) = {windowResult(w).windowId, tWin(1), tWin(2), ...
            methodLabel, fit.EO, fit.freqHz, fit.amplitudeMm, fit.weightedRmseMv, ...
            strainSources(si).channel, strainSources(si).file, rawPhase, rawTrac, ...
            bestPhase, bestTrac, bestShift, numel(tRaw), rms(yRaw, 'omitnan')}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'window_id', 'time_start_s', ...
    'time_end_s', 'method', 'EO', 'frequency_hz', 'amplitude_mm', ...
    'rmse_mV', 'strain_channel', 'strain_file', 'raw_phase_rad', ...
    'raw_TRAC', 'best_phase_rad', 'best_TRAC', 'best_time_shift_s', ...
    'strain_sample_count', 'strain_rms_mm'});
end

function [strainTime, strainMm] = load_strain_file_mm_local(strainFile, offsetSec, strainToMm)
S = load(strainFile, 'Datas');
if ~isfield(S, 'Datas') || size(S.Datas, 2) < 2
    error('Invalid strain Datas in file: %s', strainFile);
end
strainTime = S.Datas(:, 1) + offsetSec;
strainMm = (S.Datas(:, 2) - mean(S.Datas(:, 2), 'omitnan')) * strainToMm;
valid = isfinite(strainTime) & isfinite(strainMm);
strainTime = strainTime(valid);
strainMm = strainMm(valid);
[strainTime, orderIdx] = sort(strainTime);
strainMm = strainMm(orderIdx);
end

function pngFile = plot_step12_nature_validation_local(tracTable, channelTable, ...
    detailData, fftFreqHz, fftAmpMm, strainPeakHz, strainReferenceHz, ...
    orderLineHz, figDir, C0, primaryChannel, selectedTimeSec)
baseName = sprintf('Step12_NatureValidation_%s_%s', C0.dataset, C0.caseTag);
pngFile = fullfile(figDir, [baseName, '.png']);

signalColor = [0.000, 0.447, 0.741];
accentColor = [0.835, 0.369, 0.000];
mutedColor = [0.62, 0.62, 0.62];
darkGray = [0.22, 0.22, 0.22];

fig = figure('Name', 'Step12 20241106 strain validation figure', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 13.2], ...
    'NumberTitle', 'off');
set(fig, 'Renderer', 'painters');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
plot(ax, fftFreqHz, fftAmpMm, 'Color', darkGray, 'LineWidth', 0.9);
xline(ax, strainPeakHz, '-', 'Color', signalColor, 'LineWidth', 1.0);
xline(ax, strainReferenceHz, '-.', 'Color', accentColor, 'LineWidth', 0.9);
xline(ax, orderLineHz, '--', 'Color', mutedColor, 'LineWidth', 0.8);
xlim(ax, [600, 650]);
fftMask = fftFreqHz >= 600 & fftFreqHz <= 650;
ylim(ax, [0, max(fftAmpMm(fftMask)) * 1.18]);
text(ax, 0.64, 0.88, sprintf('FFT peak %.2f Hz', strainPeakHz), ...
    'Units', 'normalized', 'Color', signalColor, 'FontName', 'Arial', ...
    'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
text(ax, 0.64, 0.80, sprintf('12EO line %.2f Hz', orderLineHz), ...
    'Units', 'normalized', 'Color', mutedColor, 'FontName', 'Arial', ...
    'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
text(ax, 0.64, 0.72, sprintf('fit reference %.2f Hz', strainReferenceHz), ...
    'Units', 'normalized', 'Color', accentColor, 'FontName', 'Arial', ...
    'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'FFT amplitude (mm)');
title(ax, sprintf('%s strain spectrum, %.1f-%.1f s', ...
    primaryChannel, selectedTimeSec(1), selectedTimeSec(2)));
panel_label_local(ax, 'a');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
freqError = tracTable.btt_minus_strain_local_hz;
plot(ax, tracTable.window_id, freqError, '-o', 'Color', signalColor, ...
    'LineWidth', 1.0, 'MarkerSize', 3.8, 'MarkerFaceColor', 'w');
yline(ax, 0, '-', 'Color', mutedColor, 'LineWidth', 0.7);
xlim(ax, [min(tracTable.window_id) - 0.5, max(tracTable.window_id) + 0.5]);
ylim(ax, [min(freqError) - 0.06, max(freqError) + 0.06]);
xlabel(ax, 'Step07J window');
ylabel(ax, 'f_{BTT}-f_{strain} (Hz)');
title(ax, 'Independent window-frequency error');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
plot(ax, tracTable.window_id, tracTable.raw_TRAC, '-o', 'Color', signalColor, ...
    'LineWidth', 1.1, 'MarkerSize', 3.8, 'MarkerFaceColor', 'w');
otherChannels = setdiff(unique(channelTable.strain_channel, 'stable'), {primaryChannel});
if ~isempty(otherChannels)
    maskOther = strcmp(channelTable.strain_channel, otherChannels{1});
    plot(ax, channelTable.window_id(maskOther), channelTable.raw_TRAC(maskOther), ...
        '-', 'Color', mutedColor, 'LineWidth', 0.9);
    text(ax, max(tracTable.window_id) + 0.45, ...
        channelTable.raw_TRAC(find(maskOther, 1, 'last')), otherChannels{1}, ...
        'Color', mutedColor, 'FontName', 'Arial', 'FontSize', 6.5);
end
yline(ax, 0.90, ':', 'Color', mutedColor, 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(ax, [min(tracTable.window_id) - 0.5, max(tracTable.window_id) + 2.0]);
ylim(ax, [max(0, min([tracTable.raw_TRAC; channelTable.raw_TRAC]) - 0.08), 1.0]);
xlabel(ax, 'Step07J window');
ylabel(ax, 'TRAC');
title(ax, 'Waveform-shape agreement');
text(ax, max(tracTable.window_id) + 0.45, tracTable.raw_TRAC(end), ...
    primaryChannel, 'Color', signalColor, 'FontName', 'Arial', 'FontSize', 6.5);
panel_label_local(ax, 'c');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
if isempty(detailData)
    text(ax, 0.5, 0.5, 'No waveform detail', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
else
    s = detailData{max(1, round(numel(detailData) / 2))};
    t0 = min(s.timeBest);
    tMs = 1000 * (s.timeBest - t0);
    yComponent = project_single_frequency_component_local( ...
        s.timeBest, s.strainBest, s.freqHz);
    plot(ax, tMs, s.strainBest, '-', 'Color', [0.78, 0.78, 0.78], 'LineWidth', 0.6);
    plot(ax, tMs, yComponent, '-', 'Color', darkGray, 'LineWidth', 0.9);
    plot(ax, tMs, s.reconBest, '-', 'Color', signalColor, 'LineWidth', 1.0);
    xlabel(ax, 'Time from window start (ms)');
    ylabel(ax, 'Displacement (mm)');
    title(ax, sprintf('Representative window W%d, TRAC %.3f', ...
        s.windowId, s.bestTrac));
    text(ax, 0.64, 0.92, ['raw ', primaryChannel], 'Units', 'normalized', ...
        'Color', [0.55, 0.55, 0.55], 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.84, 'synchronous component', 'Units', 'normalized', ...
        'Color', darkGray, 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.76, 'BTT reconstruction', 'Units', 'normalized', ...
        'Color', signalColor, 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
end
panel_label_local(ax, 'd');
apply_nature_axes_local(ax);

export_nature_figure_local(fig, fullfile(figDir, baseName));
end

function [pngFile, T] = plot_channel_time_domain_compare_local(windowResult, ...
    methodName, strainSources, channelTable, figDir, C0, strainToMm, representativeIndex)
if nargin < 8 || ~isfinite(representativeIndex) || representativeIndex < 1 || ...
        representativeIndex > numel(windowResult)
    representativeIndex = max(1, round(numel(windowResult) / 2));
end
wr = windowResult(representativeIndex);
fit = wr.modelFits.(methodName);
tWin = wr.timeWindow;

rows = {};
for si = 1:numel(strainSources)
    [strainTime, strainMm] = load_strain_file_mm_local( ...
        strainSources(si).file, strainSources(si).offsetSec, strainToMm);
    rowMask = strcmp(channelTable.strain_channel, strainSources(si).channel) & ...
        channelTable.window_id == wr.windowId;
    if any(rowMask)
        bestShift = channelTable.best_time_shift_s(find(rowMask, 1, 'first'));
    else
        bestShift = 0;
    end
    [tracVal, ~, tSeg, ySeg, yRecon] = calc_trac_at_shift_local( ...
        strainTime, strainMm, tWin, bestShift, fit.amplitudeMm, fit.freqHz);
    yComponent = project_single_frequency_component_local(tSeg, ySeg, fit.freqHz);
    tMs = 1000 * (tSeg - min(tSeg));
    for ii = 1:numel(tMs)
        rows(end + 1, :) = {wr.windowId, strainSources(si).channel, ...
            strainSources(si).file, bestShift, fit.freqHz, fit.amplitudeMm, ...
            tracVal, tMs(ii), ySeg(ii), yComponent(ii), yRecon(ii)}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'window_id', 'strain_channel', ...
    'strain_file', 'best_time_shift_s', 'frequency_hz', 'amplitude_mm', ...
    'TRAC', 'time_ms', 'raw_mm', 'sync_component_mm', 'btt_reconstruction_mm'});

baseName = sprintf('Step12_ChannelTimeDomainCompare_%s_%s_W%02d', ...
    C0.dataset, C0.caseTag, wr.windowId);
pngFile = fullfile(figDir, [baseName, '.png']);

channels = unique(T.strain_channel, 'stable');
palette = [0.000, 0.447, 0.741; 0.835, 0.369, 0.000; ...
    0.25, 0.25, 0.25; 0.45, 0.45, 0.45];
darkGray = [0.22, 0.22, 0.22];
mutedGray = [0.72, 0.72, 0.72];

fig = figure('Name', 'Step12 channel time-domain comparison', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 13.0], ...
    'NumberTitle', 'off');
set(fig, 'Renderer', 'painters');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
offsetStep = max(0.12, 1.25 * max(abs(T.raw_mm), [], 'omitnan'));
for ci = 1:numel(channels)
    mask = strcmp(T.strain_channel, channels{ci});
    yOffset = (numel(channels) - ci) * offsetStep;
    plot(ax, T.time_ms(mask), T.raw_mm(mask) + yOffset, '-', ...
        'Color', palette(ci, :), 'LineWidth', 0.8);
    text(ax, max(T.time_ms(mask)) * 1.01, yOffset, channels{ci}, ...
        'Color', palette(ci, :), 'FontName', 'Arial', 'FontSize', 6.5, ...
        'VerticalAlignment', 'middle');
end
xlabel(ax, 'Time from window start (ms)');
ylabel(ax, 'Raw displacement + offset (mm)');
title(ax, sprintf('Raw channel waveforms, W%d', wr.windowId));
panel_label_local(ax, 'a');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
for ci = 1:numel(channels)
    mask = strcmp(T.strain_channel, channels{ci});
    plot(ax, T.time_ms(mask), T.sync_component_mm(mask), '-', ...
        'Color', palette(ci, :), 'LineWidth', 1.0);
    text(ax, max(T.time_ms(mask)) * 1.01, ...
        T.sync_component_mm(find(mask, 1, 'last')), channels{ci}, ...
        'Color', palette(ci, :), 'FontName', 'Arial', 'FontSize', 6.5);
end
xlabel(ax, 'Time from window start (ms)');
ylabel(ax, 'Synchronous component (mm)');
title(ax, 'Single-frequency components');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

channelScore = groupsummary(T, 'strain_channel', 'mean', 'TRAC');
[~, scoreOrder] = sort(channelScore.mean_TRAC, 'descend');
detailChannels = channelScore.strain_channel(scoreOrder);
detailChannels = detailChannels(1:min(2, numel(detailChannels)));
for ci = 1:numel(detailChannels)
    ax = nexttile;
    hold(ax, 'on');
    channelNow = detailChannels{ci};
    mask = strcmp(T.strain_channel, channelNow);
    plot(ax, T.time_ms(mask), T.raw_mm(mask), '-', 'Color', mutedGray, 'LineWidth', 0.55);
    plot(ax, T.time_ms(mask), T.sync_component_mm(mask), '-', ...
        'Color', darkGray, 'LineWidth', 0.9);
    colorIdx = find(strcmp(channels, channelNow), 1, 'first');
    plot(ax, T.time_ms(mask), T.btt_reconstruction_mm(mask), '-', ...
        'Color', palette(colorIdx, :), 'LineWidth', 1.0);
    tracVal = T.TRAC(find(mask, 1, 'first'));
    title(ax, sprintf('%s vs BTT, TRAC %.3f', channelNow, tracVal));
    xlabel(ax, 'Time from window start (ms)');
    ylabel(ax, 'Displacement (mm)');
    text(ax, 0.64, 0.92, 'raw', 'Units', 'normalized', 'Color', [0.55, 0.55, 0.55], ...
        'FontName', 'Arial', 'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.84, 'sync component', 'Units', 'normalized', 'Color', darkGray, ...
        'FontName', 'Arial', 'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
    text(ax, 0.64, 0.76, 'BTT reconstruction', 'Units', 'normalized', ...
        'Color', palette(colorIdx, :), 'FontName', 'Arial', 'FontSize', 6.5, ...
        'BackgroundColor', 'w', 'Margin', 1);
    panel_label_local(ax, char('b' + ci));
    apply_nature_axes_local(ax);
end

export_nature_figure_local(fig, fullfile(figDir, baseName));
end

function T = build_step07j_order_audit_local(~, ~, primaryFile, methodName, ...
    strainReferenceHz, orderLineHz)
files = {primaryFile};
labels = {'S257 physical B4'};
rows = {};
for fi = 1:numel(files)
    S = load(files{fi}, 'Result');
    W = S.Result.WindowResult(:);
    for wi = 1:numel(W)
        fit = W(wi).modelFits.(methodName);
        rows(end + 1, :) = {labels{fi}, files{fi}, W(wi).windowId, ...
            W(wi).timeWindow(1), W(wi).timeWindow(2), fit.EO, fit.freqHz, ...
            fit.amplitudeMm, fit.weightedRmseMv, fit.freqHz - strainReferenceHz, ...
            fit.freqHz - orderLineHz}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'result_label', 'result_file', ...
    'window_id', 'time_start_s', 'time_end_s', 'EO', 'frequency_hz', ...
    'amplitude_mm', 'rmse_mV', 'freq_error_vs_strain_reference_hz', ...
    'freq_error_vs_order_line_hz'});
end

function pngFile = plot_step12_order_audit_local(auditTable, figDir, C0, ...
    strainPeakHz, orderLineHz)
baseName = sprintf('Step12_OrderAudit_%s_%s', ...
    C0.dataset, C0.caseTag);
pngFile = fullfile(figDir, [baseName, '.png']);

signalColor = [0.000, 0.447, 0.741];
accentColor = [0.835, 0.369, 0.000];
mutedColor = [0.62, 0.62, 0.62];
darkGray = [0.22, 0.22, 0.22];

fig = figure('Name', 'Step12 20241106 order audit', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18.3, 8.8], ...
    'NumberTitle', 'off');
set(fig, 'Renderer', 'painters');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

labels = unique(auditTable.result_label, 'stable');
colors = [signalColor; accentColor; darkGray];

ax = nexttile;
hold(ax, 'on');
for li = 1:numel(labels)
    mask = strcmp(auditTable.result_label, labels{li});
    plot(ax, auditTable.window_id(mask), auditTable.EO(mask), '-o', ...
        'Color', colors(li, :), 'LineWidth', 1.0, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 3.8);
end
yline(ax, 12, '--', '12EO', 'Color', mutedColor, 'LineWidth', 0.8);
xlabel(ax, 'Step07J window');
ylabel(ax, 'Selected EO');
title(ax, 'Order stability');
legend(ax, labels, 'Location', 'best', 'Box', 'off');
panel_label_local(ax, 'a');
apply_nature_axes_local(ax);

ax = nexttile;
hold(ax, 'on');
for li = 1:numel(labels)
    mask = strcmp(auditTable.result_label, labels{li});
    plot(ax, auditTable.window_id(mask), auditTable.frequency_hz(mask), '-o', ...
        'Color', colors(li, :), 'LineWidth', 1.0, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 3.8);
end
yline(ax, strainPeakHz, '-', 'Color', mutedColor, 'LineWidth', 0.8);
yline(ax, orderLineHz, '--', 'Color', mutedColor, 'LineWidth', 0.8);
text(ax, 0.04, 0.48, sprintf('strain %.2f Hz', strainPeakHz), ...
    'Units', 'normalized', 'Color', mutedColor, 'FontName', 'Arial', ...
    'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
text(ax, 0.04, 0.42, sprintf('12EO %.2f Hz', orderLineHz), ...
    'Units', 'normalized', 'Color', mutedColor, 'FontName', 'Arial', ...
    'FontSize', 6.5, 'BackgroundColor', 'w', 'Margin', 1);
xlabel(ax, 'Step07J window');
ylabel(ax, 'Frequency (Hz)');
title(ax, 'Frequency consistency');
panel_label_local(ax, 'b');
apply_nature_axes_local(ax);

export_nature_figure_local(fig, fullfile(figDir, baseName));
end

function apply_nature_axes_local(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 7, 'LineWidth', 0.6, ...
    'TickDir', 'out', 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off');
ax.Title.FontSize = 7.5;
ax.Title.FontWeight = 'normal';
ax.XLabel.FontSize = 7.5;
ax.YLabel.FontSize = 7.5;
end

function panel_label_local(ax, labelText)
text(ax, -0.12, 1.08, labelText, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function export_nature_figure_local(fig, baseFile)
exportgraphics(fig, [baseFile, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [baseFile, '.png'], 'Resolution', 600);
print(fig, [baseFile, '.tiff'], '-dtiff', '-r600');
try
    saveas(fig, [baseFile, '.emf']);
catch ME
    warning('EMF export failed for %s: %s', baseFile, ME.message);
end
end

function yComponent = project_single_frequency_component_local(t, y, freqHz)
t = t(:);
y = y(:);
valid = isfinite(t) & isfinite(y);
t = t(valid);
y = y(valid);
if isempty(t)
    yComponent = [];
    return;
end
y = detrend(y - mean(y, 'omitnan'));
design = [sin(2 * pi * freqHz * t), cos(2 * pi * freqHz * t)];
coef = design \ y;
yComponent = design * coef;
end

function label = channel_label_from_file_local(pathText)
[~, name] = fileparts(pathText);
token = regexp(name, 'AI1-\d+', 'match', 'once');
if isempty(token)
    label = name;
else
    label = token;
end
end

function validate_current_step07j_result_local(R, C0, methodName, resultFile, ...
    expectedStartTimeSec)
if ~isfield(R, 'flowConfig') || ~isfield(R.flowConfig, 'identification')
    error('Step07J result lacks the saved identification configuration: %s', resultFile);
end
I = R.flowConfig.identification;
if I.targetBlade ~= C0.bladeId || ...
        ~isequal(I.analysisSensors(:).', C0.sensorIds(:).')
    error('Step07J result is not physical B%d with sensors %s.', ...
        C0.bladeId, mat2str(C0.sensorIds));
end
if ~isfield(I, 'analysisStartTimeSec') || ...
        abs(I.analysisStartTimeSec - expectedStartTimeSec) > 1e-9
    error('Step07J result does not use the selected %.6f s start time.', ...
        expectedStartTimeSec);
end
if ~isfield(R, 'cfg') || ~isfield(R.cfg, 'localBundleSource') || ...
        ~strcmpi(R.cfg.localBundleSource, 'foundation_step05_bundle')
    error('Step07J result did not use the Foundation Step05 waveform bundle.');
end
if ~isfield(R, 'FoundationStep05TimeCheck') || ...
        ~isfield(R.FoundationStep05TimeCheck, 'action') || ...
        ~strcmpi(R.FoundationStep05TimeCheck.action, 'foundation_step05_bundle')
    error('Step07J Foundation time-source consistency check was not passed.');
end
if ~isfield(R, 'WindowResult') || isempty(R.WindowResult)
    error('Step07J result has no windows.');
end
if ~isfield(R, 'dynamicFile') || ~isfile(R.dynamicFile)
    error('Step07J DynamicMap source is missing.');
end
D = load(R.dynamicFile, 'DynamicMap');
if ~isfield(D, 'DynamicMap') || ~isfield(D.DynamicMap, 'BladeNumbering') || ...
        ~isfield(D.DynamicMap.BladeNumbering, 'scheme') || ...
        ~strcmpi(D.DynamicMap.BladeNumbering.scheme, ...
        'foundation_physical_blade_id')
    error('Step07J DynamicMap does not use Foundation physical blade IDs.');
end
plannedWindowCount = floor((I.targetBladePasses - I.windowBladePasses) / ...
    I.slidingStepBladePasses) + 1;
if numel(R.WindowResult) ~= plannedWindowCount
    error('Step07J result has %d windows; the saved configuration requires %d.', ...
        numel(R.WindowResult), plannedWindowCount);
end
for iWindow = 1:numel(R.WindowResult)
    if ~isfield(R.WindowResult(iWindow), 'modelFits') || ...
            ~isfield(R.WindowResult(iWindow).modelFits, methodName)
        error('Step07J window %d has no %s result.', iWindow, methodName);
    end
end
fprintf(['Source audit passed: physical B%d, sensors %s, %s, ' ...
    'Foundation Step05, %d windows.\n'], C0.bladeId, ...
    mat2str(C0.sensorIds), methodName, numel(R.WindowResult));
end
