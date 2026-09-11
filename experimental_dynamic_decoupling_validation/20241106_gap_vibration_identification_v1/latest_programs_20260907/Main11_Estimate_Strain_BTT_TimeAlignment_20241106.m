function Alignment = Main11_Estimate_Strain_BTT_TimeAlignment_20241106()
% Estimate one strain-to-BTT clock offset from frequency and TARC evidence.
%
% Sign convention:
%   t_BTT = t_strain_raw + strainToBttOffsetSec.
%
% The variable-speed 84-89 s interval first identifies the correct coarse
% clock-offset family by matching the strain ridge to the OPR 12EO line.
% Phase-consistent TARC then combines the 75.0 s and 86.6 s Step07J cases.
% One unknown global strain/displacement phase is allowed, but each window
% is not allowed to choose an independent phase. This makes TARC useful for
% time alignment while retaining amplitude as an independent diagnostic.

close all;

%% 1. Fixed evidence and search settings
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir, fullfile(thisDir, 'functions', 'foundation'), ...
    fullfile(thisDir, 'functions', 'preparation'), ...
    fullfile(thisDir, 'functions', 'gap_aware'), ...
    fullfile(thisDir, 'functions', 'utilities'));
packageCfg = Config_20241106();
outDir = packageCfg.paths.strainValidation;
figureDir = fullfile(outDir, 'figures_step12_trac_strain_spectrum');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

expectedStartTimeSec = packageCfg.strain.alignmentCasesSec(:);
resultFiles = cell(numel(expectedStartTimeSec), 1);
caseLabels = strings(numel(expectedStartTimeSec), 1);
for iCase = 1:numel(expectedStartTimeSec)
    caseLabels(iCase) = string(strrep(sprintf('T%07.3f', ...
        expectedStartTimeSec(iCase)), '.', 'p'));
    resultFiles{iCase} = formal_gap_result_file_local(packageCfg, ...
        expectedStartTimeSec(iCase));
end
step04File = fullfile(packageCfg.paths.strainEvidence, ...
    'Step04_BTT_STE_Resonance_Regions_20241106.mat');
caseTag = sprintf('B%d_S%s', packageCfg.case.targetBlade, ...
    sprintf('%d', packageCfg.case.analysisSensors));

alignmentFile = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_20241106_%s.mat', caseTag));
summaryCsv = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_Summary_20241106_%s.csv', caseTag));
scanCsv = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_Scan_20241106_%s.csv', caseTag));
fineScanCsv = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_FineScan_20241106_%s.csv', caseTag));
windowCsv = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_WindowDiagnostics_20241106_%s.csv', caseTag));
amplitudeCsv = fullfile(outDir, sprintf( ...
    'Step12_GlobalTimeAlignment_AmplitudeDiagnostics_20241106_%s.csv', caseTag));
figureFile = fullfile(figureDir, sprintf( ...
    'Step12_GlobalTimeAlignment_20241106_%s.png', caseTag));

priorOffsetSec = packageCfg.strain.timeOffsetPriorSec;
coarseSearchHalfWidthSec = 1.0;
coarseSearchStepSec = 0.001;
phaseFamilyGateSec = 0.15;
fineSearchHalfWidthSec = 0.005;
fineSearchStepSec = 0.0001;
alignmentBttSpanSec = [84, 89];

ridgeWindowSec = 0.75;
ridgeStepSec = 0.05;
ridgeFrequencyBandHz = [620, 640];
minimumRidgeQuality = 5;
oprSampleRateHz = 5e6;
oprSmoothRevolutions = 9;

transferRatio_1PerM = 1.130867104785;
strainToMm = 1 / transferRatio_1PerM / 1000;

%% 2. Load the two trusted Step07J results and raw strain evidence
results = cell(numel(resultFiles), 1);
for iCase = 1:numel(resultFiles)
    if ~isfile(resultFiles{iCase})
        error('Missing Step07J result: %s', resultFiles{iCase});
    end
    loaded = load(resultFiles{iCase}, 'Result');
    validate_result_local(loaded.Result, expectedStartTimeSec(iCase));
    results{iCase} = loaded.Result;
end

S4 = load(step04File, 'strain', 'strainToBttOffsetSec', ...
    'selectedStrainChannel');
if abs(S4.strainToBttOffsetSec - priorOffsetSec) > 1e-9
    error('Step04 offset %.9f s does not match prior %.9f s.', ...
        S4.strainToBttOffsetSec, priorOffsetSec);
end

rawTime = S4.strain.time(:);
rawMicrostrain = S4.strain.value(:);
valid = isfinite(rawTime) & isfinite(rawMicrostrain);
rawTime = rawTime(valid);
rawMicrostrain = rawMicrostrain(valid);
[rawTime, order] = sort(rawTime);
rawMicrostrain = rawMicrostrain(order);
sampleDt = median(diff(rawTime), 'omitnan');
rawMm = rawMicrostrain .* strainToMm;

[oprTime, order12FrequencyHz, oprFile] = load_opr_order_frequency_local( ...
    S4.strain.file, oprSampleRateHz, oprSmoothRevolutions);

%% 3. Coarse alignment from the variable-speed strain frequency ridge
rawCenterSpanSec = [ ...
    alignmentBttSpanSec(1) - priorOffsetSec - coarseSearchHalfWidthSec, ...
    alignmentBttSpanSec(2) - priorOffsetSec + coarseSearchHalfWidthSec];
ridge = estimate_strain_frequency_ridge_local(rawTime, rawMicrostrain, ...
    rawCenterSpanSec, ridgeWindowSec, ridgeStepSec, ...
    ridgeFrequencyBandHz);

offsetGrid = (priorOffsetSec - coarseSearchHalfWidthSec: ...
    coarseSearchStepSec:priorOffsetSec + coarseSearchHalfWidthSec).';
frequencyMetrics = evaluate_frequency_alignment_local(offsetGrid, ridge, ...
    oprTime, order12FrequencyHz, alignmentBttSpanSec, minimumRidgeQuality);
[bestFrequencyRmseHz, coarseIndex] = min(frequencyMetrics.rmseHz);
coarseOffsetSec = offsetGrid(coarseIndex);

%% 4. Select the phase family using globally phase-consistent TARC
phaseModels = build_phase_models_local(results, caseLabels);
phaseMetrics = evaluate_phase_alignment_local(offsetGrid, phaseModels, ...
    rawTime, rawMicrostrain, sampleDt);

localPeak = islocalmax(phaseMetrics.combinedConcentration);
eligiblePeak = localPeak & ...
    abs(offsetGrid - coarseOffsetSec) <= phaseFamilyGateSec & ...
    isfinite(phaseMetrics.combinedConcentration);
if ~any(eligiblePeak)
    eligiblePeak = abs(offsetGrid - coarseOffsetSec) <= phaseFamilyGateSec & ...
        isfinite(phaseMetrics.combinedConcentration);
end
eligibleIndex = find(eligiblePeak);
[~, bestEligible] = max(phaseMetrics.combinedConcentration(eligibleIndex));
phaseFamilyIndex = eligibleIndex(bestEligible);
phaseFamilyOffsetSec = offsetGrid(phaseFamilyIndex);

fineOffsetGrid = (phaseFamilyOffsetSec - fineSearchHalfWidthSec: ...
    fineSearchStepSec:phaseFamilyOffsetSec + fineSearchHalfWidthSec).';
finePhaseMetrics = evaluate_phase_alignment_local(fineOffsetGrid, ...
    phaseModels, rawTime, rawMicrostrain, sampleDt);
fineFrequencyMetrics = evaluate_frequency_alignment_local(fineOffsetGrid, ...
    ridge, oprTime, order12FrequencyHz, alignmentBttSpanSec, ...
    minimumRidgeQuality);
[bestPhaseConcentration, fineIndex] = max( ...
    finePhaseMetrics.combinedConcentration);
estimatedOffsetSec = fineOffsetGrid(fineIndex);
estimatedFrequencyRmseHz = fineFrequencyMetrics.rmseHz(fineIndex);
estimatedMedianTarc = finePhaseMetrics.medianTarc(fineIndex, :);

peakValues = phaseMetrics.combinedConcentration(localPeak);
peakValues = sort(peakValues(isfinite(peakValues)), 'descend');
if numel(peakValues) >= 2
    phaseFamilyGap = peakValues(1) - peakValues(2);
else
    phaseFamilyGap = inf;
end

[~, priorIndex] = min(abs(offsetGrid - priorOffsetSec));
priorFrequencyRmseHz = frequencyMetrics.rmseHz(priorIndex);
priorPhaseConcentration = phaseMetrics.combinedConcentration(priorIndex);

accepted = coarseIndex > 1 && coarseIndex < numel(offsetGrid) && ...
    fineIndex > 1 && fineIndex < numel(fineOffsetGrid) && ...
    abs(estimatedOffsetSec - coarseOffsetSec) <= phaseFamilyGateSec && ...
    estimatedFrequencyRmseHz <= bestFrequencyRmseHz + 0.03 && ...
    bestPhaseConcentration >= 0.98 && ...
    all(estimatedMedianTarc >= 0.90) && phaseFamilyGap >= 0.002;

if accepted
    appliedOffsetSec = estimatedOffsetSec;
    decision = 'frequency_coarse_global_phase_tarc_refinement';
else
    appliedOffsetSec = priorOffsetSec;
    decision = 'retain_prior_alignment_checks_failed';
end

%% 5. Independent amplitude and per-window phase diagnostics
[~, windowDiagnostics] = phase_metrics_at_offset_local(appliedOffsetSec, ...
    phaseModels, rawTime, rawMicrostrain, sampleDt, true);
amplitudePrior = amplitude_diagnostics_local(phaseModels, priorOffsetSec, ...
    rawTime, rawMm, sampleDt, "prior");
amplitudeApplied = amplitude_diagnostics_local(phaseModels, appliedOffsetSec, ...
    rawTime, rawMm, sampleDt, "applied");
amplitudeDiagnostics = [amplitudePrior; amplitudeApplied];

prior75 = amplitudePrior(amplitudePrior.case_label == "T075p000", :);
applied75 = amplitudeApplied(amplitudeApplied.case_label == "T075p000", :);
applied86 = amplitudeApplied(amplitudeApplied.case_label == "T086p600", :);

%% 6. Save the shared alignment result and audit tables
Alignment = struct();
Alignment.dataset = '20241106';
Alignment.method = 'variable_speed_frequency_plus_global_phase_tarc';
Alignment.bladeId = 4;
Alignment.sensorIds = [2, 5, 7];
Alignment.strainChannel = string(S4.selectedStrainChannel);
Alignment.signConvention = 't_BTT = t_strain_raw + strainToBttOffsetSec';
Alignment.priorOffsetSec = priorOffsetSec;
Alignment.coarseOffsetSec = coarseOffsetSec;
Alignment.phaseFamilyOffsetSec = phaseFamilyOffsetSec;
Alignment.estimatedOffsetSec = estimatedOffsetSec;
Alignment.fineCorrectionFromPriorSec = estimatedOffsetSec - priorOffsetSec;
Alignment.strainToBttOffsetSec = appliedOffsetSec;
Alignment.accepted = accepted;
Alignment.decision = decision;
Alignment.alignmentBttSpanSec = alignmentBttSpanSec;
Alignment.priorFrequencyRmseHz = priorFrequencyRmseHz;
Alignment.coarseFrequencyRmseHz = bestFrequencyRmseHz;
Alignment.estimatedFrequencyRmseHz = estimatedFrequencyRmseHz;
Alignment.priorPhaseConcentration = priorPhaseConcentration;
Alignment.estimatedPhaseConcentration = bestPhaseConcentration;
Alignment.phaseFamilyGap = phaseFamilyGap;
Alignment.estimatedMedianTarc75 = estimatedMedianTarc(1);
Alignment.estimatedMedianTarc86 = estimatedMedianTarc(2);
Alignment.priorMape = prior75.amplitude_mape_percent / 100;
Alignment.estimatedMape = applied75.amplitude_mape_percent / 100;
Alignment.estimatedAmplitudeMape75Percent = applied75.amplitude_mape_percent;
Alignment.estimatedAmplitudeMape86Percent = applied86.amplitude_mape_percent;
Alignment.searchHalfWidthSec = coarseSearchHalfWidthSec;
Alignment.searchStepSec = coarseSearchStepSec;
Alignment.fineSearchStepSec = fineSearchStepSec;
Alignment.transferRatio_1PerM = transferRatio_1PerM;
Alignment.resultFiles = strjoin(string(resultFiles), '; ');
Alignment.step04File = string(step04File);
Alignment.oprFile = string(oprFile);
Alignment.createdAt = datetime('now');

save(alignmentFile, 'Alignment', 'offsetGrid', 'frequencyMetrics', ...
    'phaseMetrics', 'fineOffsetGrid', 'fineFrequencyMetrics', ...
    'finePhaseMetrics', 'ridge', 'windowDiagnostics', ...
    'amplitudeDiagnostics', '-v7');

summary = table(string(Alignment.dataset), string(Alignment.method), ...
    Alignment.priorOffsetSec, Alignment.coarseOffsetSec, ...
    Alignment.phaseFamilyOffsetSec, Alignment.estimatedOffsetSec, ...
    Alignment.strainToBttOffsetSec, Alignment.accepted, ...
    string(Alignment.decision), Alignment.priorFrequencyRmseHz, ...
    Alignment.coarseFrequencyRmseHz, Alignment.estimatedFrequencyRmseHz, ...
    Alignment.priorPhaseConcentration, ...
    Alignment.estimatedPhaseConcentration, Alignment.phaseFamilyGap, ...
    Alignment.estimatedMedianTarc75, Alignment.estimatedMedianTarc86, ...
    Alignment.estimatedAmplitudeMape75Percent, ...
    Alignment.estimatedAmplitudeMape86Percent, ...
    'VariableNames', {'dataset', 'method', 'prior_offset_s', ...
    'coarse_offset_s', 'phase_family_offset_s', 'estimated_offset_s', ...
    'applied_offset_s', 'accepted', 'decision', ...
    'prior_frequency_rmse_hz', 'coarse_frequency_rmse_hz', ...
    'estimated_frequency_rmse_hz', 'prior_phase_concentration', ...
    'estimated_phase_concentration', 'phase_family_gap', ...
    'median_tarc_75', 'median_tarc_86', ...
    'amplitude_mape_75_percent', 'amplitude_mape_86_percent'});
writetable(summary, summaryCsv);

scanTable = make_scan_table_local(offsetGrid, priorOffsetSec, ...
    frequencyMetrics, phaseMetrics);
fineScanTable = make_scan_table_local(fineOffsetGrid, priorOffsetSec, ...
    fineFrequencyMetrics, finePhaseMetrics);
writetable(scanTable, scanCsv);
writetable(fineScanTable, fineScanCsv);
writetable(windowDiagnostics, windowCsv);
writetable(amplitudeDiagnostics, amplitudeCsv);

plot_alignment_local(ridge, oprTime, order12FrequencyHz, ...
    alignmentBttSpanSec, offsetGrid, frequencyMetrics, phaseMetrics, ...
    priorOffsetSec, coarseOffsetSec, appliedOffsetSec, ...
    windowDiagnostics, estimatedMedianTarc, figureFile);

fprintf('Frequency + global-phase TARC alignment complete.\n');
fprintf('  Prior offset       = %.6f s | frequency RMSE %.6f Hz | phase C %.6f\n', ...
    priorOffsetSec, priorFrequencyRmseHz, priorPhaseConcentration);
fprintf('  Frequency coarse   = %.6f s | frequency RMSE %.6f Hz\n', ...
    coarseOffsetSec, bestFrequencyRmseHz);
fprintf('  Phase-family seed  = %.6f s\n', phaseFamilyOffsetSec);
fprintf(['  Estimated offset   = %.6f s | correction %+.1f ms | ' ...
    'frequency RMSE %.6f Hz\n'], estimatedOffsetSec, ...
    1000 * (estimatedOffsetSec - priorOffsetSec), estimatedFrequencyRmseHz);
fprintf('  Global phase C     = %.6f | median TARC [75, 86.6] = [%.6f, %.6f]\n', ...
    bestPhaseConcentration, estimatedMedianTarc(1), estimatedMedianTarc(2));
fprintf('  Applied offset     = %.6f s | accepted=%d | %s\n', ...
    appliedOffsetSec, accepted, decision);
fprintf('  Amplitude MAPE     = %.3f%% at 75 s | %.3f%% at 86.6 s (diagnostic only)\n', ...
    applied75.amplitude_mape_percent, applied86.amplitude_mape_percent);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    alignmentFile, summaryCsv, scanCsv, fineScanCsv, windowCsv, ...
    amplitudeCsv, figureFile);
end

function file = formal_gap_result_file_local(cfg, analysisTimeSec)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', analysisTimeSec), '.', 'p');
suffix = sprintf('_%s_main_gapfixedtilt_W%dS%d_F%dto%dHz_DFpm%gHz', ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)), ...
    cfg.frequency.refineHalfWidthHz);
file = fullfile(cfg.paths.gapResults, sprintf( ...
    'Step07J_NestedStaticWarp_VPFullWave_20241106_B%d_%s%s.mat', ...
    cfg.case.targetBlade, sensorTag, suffix));
end

function validate_result_local(R, expectedStartTimeSec)
if ~isfield(R, 'flowConfig') || ...
        R.flowConfig.identification.targetBlade ~= 4 || ...
        ~isequal(R.flowConfig.identification.analysisSensors(:).', [2, 5, 7])
    error('Unexpected Step07J blade or sensor set.');
end
if ~isfield(R, 'cfg') || ...
        ~strcmpi(R.cfg.localBundleSource, 'foundation_step05_bundle')
    error('Alignment requires the trusted Foundation Step05 result.');
end
windowBounds = vertcat(R.WindowResult.timeWindow);
if abs(min(windowBounds(:, 1)) - expectedStartTimeSec) > 0.1
    error('Step07J result does not match expected %.1f s case.', ...
        expectedStartTimeSec);
end
end

function [oprTime, order12FrequencyHz, oprFile] = ...
    load_opr_order_frequency_local(strainFile, sampleRateHz, smoothRevolutions)
oprFile = fullfile(fileparts(strainFile), '3000_3150', 'jiluOPR.mat');
if ~isfile(oprFile)
    error('Missing OPR file: %s', oprFile);
end
loaded = load(oprFile, 'jiluOPR');
oprCenter = mean(loaded.jiluOPR(:, 1:2), 2) / sampleRateHz;
oprTime = oprCenter(1:end-1);
order12FrequencyHz = 12 ./ diff(oprCenter);
order12FrequencyHz = movmedian(order12FrequencyHz, smoothRevolutions, ...
    'omitnan');
order12FrequencyHz = movmean(order12FrequencyHz, smoothRevolutions, ...
    'omitnan');
end

function ridge = estimate_strain_frequency_ridge_local(rawTime, rawValue, ...
    centerSpanSec, windowSec, stepSec, frequencyBandHz)
sampleDt = median(diff(rawTime), 'omitnan');
sampleRateHz = 1 / sampleDt;
sampleCount = round(windowSec * sampleRateHz);
if sampleCount < 32
    error('Frequency-ridge window is too short.');
end
window = 0.5 - 0.5 * cos(2 * pi * (0:sampleCount-1).' / ...
    (sampleCount - 1));
nfft = 2^nextpow2(16 * sampleCount);
frequencyHz = (0:nfft/2).' * sampleRateHz / nfft;
bandIndex = find(frequencyHz >= frequencyBandHz(1) & ...
    frequencyHz <= frequencyBandHz(2));
centerTimeSec = (centerSpanSec(1):stepSec:centerSpanSec(2)).';
peakFrequencyHz = nan(size(centerTimeSec));
peakQuality = nan(size(centerTimeSec));

for i = 1:numel(centerTimeSec)
    centerIndex = round((centerTimeSec(i) - rawTime(1)) / sampleDt) + 1;
    firstIndex = centerIndex - floor(sampleCount / 2);
    index = firstIndex:(firstIndex + sampleCount - 1);
    if index(1) < 1 || index(end) > numel(rawValue)
        continue;
    end
    y = detrend(rawValue(index));
    spectrum = abs(fft(y .* window, nfft));
    spectrum = spectrum(1:nfft/2 + 1);
    [peakValue, localIndex] = max(spectrum(bandIndex));
    peakIndex = bandIndex(localIndex);
    peakFrequencyHz(i) = quadratic_peak_frequency_local( ...
        frequencyHz, spectrum, peakIndex);
    peakQuality(i) = peakValue / max(median(spectrum(bandIndex)), eps);
end

ridge = struct('rawCenterTimeSec', centerTimeSec, ...
    'peakFrequencyHz', peakFrequencyHz, 'quality', peakQuality, ...
    'windowSec', windowSec, 'stepSec', stepSec, ...
    'frequencyBandHz', frequencyBandHz);
end

function peakFrequencyHz = quadratic_peak_frequency_local(frequencyHz, ...
    spectrum, peakIndex)
peakFrequencyHz = frequencyHz(peakIndex);
if peakIndex <= 1 || peakIndex >= numel(spectrum)
    return;
end
logAmplitude = log(max(spectrum(peakIndex-1:peakIndex+1), eps));
denominator = logAmplitude(1) - 2 * logAmplitude(2) + logAmplitude(3);
if abs(denominator) <= eps
    return;
end
deltaBin = 0.5 * (logAmplitude(1) - logAmplitude(3)) / denominator;
peakFrequencyHz = peakFrequencyHz + ...
    deltaBin * (frequencyHz(2) - frequencyHz(1));
end

function metrics = evaluate_frequency_alignment_local(offsetGrid, ridge, ...
    oprTime, orderFrequencyHz, bttSpanSec, minimumQuality)
n = numel(offsetGrid);
metrics = struct('rmseHz', nan(n, 1), 'correlation', nan(n, 1), ...
    'meanBiasHz', nan(n, 1), 'sampleCount', zeros(n, 1));
for i = 1:n
    bttTime = ridge.rawCenterTimeSec + offsetGrid(i);
    targetFrequencyHz = interp1(oprTime, orderFrequencyHz, bttTime, ...
        'linear', NaN);
    good = bttTime >= bttSpanSec(1) & bttTime <= bttSpanSec(2) & ...
        isfinite(ridge.peakFrequencyHz) & isfinite(targetFrequencyHz) & ...
        ridge.quality >= minimumQuality;
    if nnz(good) < 20
        continue;
    end
    errorHz = ridge.peakFrequencyHz(good) - targetFrequencyHz(good);
    metrics.rmseHz(i) = sqrt(mean(errorHz .^ 2));
    metrics.correlation(i) = signed_correlation_local( ...
        ridge.peakFrequencyHz(good), targetFrequencyHz(good));
    metrics.meanBiasHz(i) = mean(errorHz);
    metrics.sampleCount(i) = nnz(good);
end
end

function phaseModels = build_phase_models_local(results, caseLabels)
phaseModels = cell(numel(results), 1);
for iCase = 1:numel(results)
    windowResult = results{iCase}.WindowResult(:);
    models = repmat(struct('caseLabel', "", 'windowId', NaN, ...
        'timeWindow', [NaN NaN], 'timeCenterSec', NaN, ...
        'frequencyHz', NaN, 'bttAmplitudeMm', NaN, ...
        'absolutePhaseRad', NaN, 'phaseModelR2', NaN), ...
        numel(windowResult), 1);
    for iWindow = 1:numel(windowResult)
        fit = windowResult(iWindow).modelFits.gap_only;
        bttTime = windowResult(iWindow).bundle.T(:);
        bttDisplacementMm = fit.uMm(:);
        omegaT = 2 * pi * fit.freqHz .* bttTime;
        design = [sin(omegaT), cos(omegaT), ones(size(bttTime))];
        coefficient = design \ bttDisplacementMm;
        reconstructed = design * coefficient;
        totalEnergy = sum((bttDisplacementMm - ...
            mean(bttDisplacementMm)) .^ 2);

        models(iWindow).caseLabel = caseLabels(iCase);
        models(iWindow).windowId = windowResult(iWindow).windowId;
        models(iWindow).timeWindow = windowResult(iWindow).timeWindow;
        models(iWindow).timeCenterSec = mean(windowResult(iWindow).timeWindow);
        models(iWindow).frequencyHz = fit.freqHz;
        models(iWindow).bttAmplitudeMm = fit.amplitudeMm;
        models(iWindow).absolutePhaseRad = atan2( ...
            coefficient(2), coefficient(1));
        models(iWindow).phaseModelR2 = 1 - ...
            sum((bttDisplacementMm - reconstructed) .^ 2) / ...
            max(totalEnergy, eps);
        if models(iWindow).phaseModelR2 < 0.99
            error('BTT absolute-phase reconstruction failed for %s W%d.', ...
                caseLabels(iCase), windowResult(iWindow).windowId);
        end
    end
    phaseModels{iCase} = models;
end
end

function metrics = evaluate_phase_alignment_local(offsetGrid, phaseModels, ...
    rawTime, rawValue, sampleDt)
n = numel(offsetGrid);
nCase = numel(phaseModels);
metrics = struct('combinedConcentration', nan(n, 1), ...
    'caseConcentration', nan(n, nCase), ...
    'medianTarc', nan(n, nCase), ...
    'globalPhaseRad', nan(n, 1));
for i = 1:n
    metric = phase_metrics_at_offset_local(offsetGrid(i), phaseModels, ...
        rawTime, rawValue, sampleDt, false);
    metrics.combinedConcentration(i) = metric.combinedConcentration;
    metrics.caseConcentration(i, :) = metric.caseConcentration;
    metrics.medianTarc(i, :) = metric.medianTarc;
    metrics.globalPhaseRad(i) = metric.globalPhaseRad;
end
end

function [metric, diagnostics] = phase_metrics_at_offset_local(offsetSec, ...
    phaseModels, rawTime, rawValue, sampleDt, makeDiagnostics)
nCase = numel(phaseModels);
caseConcentration = nan(1, nCase);
medianTarc = nan(1, nCase);
allResidual = [];
allWeight = [];
caseColumn = strings(0, 1);
windowColumn = zeros(0, 1);
timeColumn = zeros(0, 1);
frequencyColumn = zeros(0, 1);
tarcColumn = zeros(0, 1);
residualColumn = zeros(0, 1);

for iCase = 1:nCase
    models = phaseModels{iCase};
    residual = nan(numel(models), 1);
    tarc = nan(numel(models), 1);
    for iWindow = 1:numel(models)
        rawWindow = models(iWindow).timeWindow - offsetSec;
        index = time_window_indices_local(rawTime, sampleDt, rawWindow);
        if numel(index) < 20
            continue;
        end
        t = rawTime(index);
        y = detrend(rawValue(index));
        omegaT = 2 * pi * models(iWindow).frequencyHz .* t;
        design = [sin(omegaT), cos(omegaT), ones(size(t))];
        coefficient = design \ y;
        yFit = design * coefficient;
        tarc(iWindow) = squared_correlation_local(y, yFit);
        strainRawPhase = atan2(coefficient(2), coefficient(1));
        residual(iWindow) = angle(exp(1i * (strainRawPhase - ...
            2 * pi * models(iWindow).frequencyHz * offsetSec - ...
            models(iWindow).absolutePhaseRad)));
    end
    weight = sqrt(max(tarc, 0));
    good = isfinite(weight) & isfinite(residual) & weight > 0;
    if any(good)
        caseConcentration(iCase) = abs(sum(weight(good) .* ...
            exp(1i * residual(good)))) / sum(weight(good));
        medianTarc(iCase) = median(tarc(good));
        allResidual = [allResidual; residual(good)]; %#ok<AGROW>
        allWeight = [allWeight; weight(good)]; %#ok<AGROW>
    end
    if makeDiagnostics
        caseColumn = [caseColumn; repmat(models(1).caseLabel, ...
            numel(models), 1)]; %#ok<AGROW>
        windowColumn = [windowColumn; [models.windowId].']; %#ok<AGROW>
        timeColumn = [timeColumn; [models.timeCenterSec].']; %#ok<AGROW>
        frequencyColumn = [frequencyColumn; [models.frequencyHz].']; %#ok<AGROW>
        tarcColumn = [tarcColumn; tarc]; %#ok<AGROW>
        residualColumn = [residualColumn; residual]; %#ok<AGROW>
    end
end

if isempty(allWeight)
    combinedConcentration = NaN;
    globalPhaseRad = NaN;
else
    phaseVector = sum(allWeight .* exp(1i * allResidual));
    combinedConcentration = abs(phaseVector) / sum(allWeight);
    globalPhaseRad = angle(phaseVector);
end
metric = struct('combinedConcentration', combinedConcentration, ...
    'caseConcentration', caseConcentration, 'medianTarc', medianTarc, ...
    'globalPhaseRad', globalPhaseRad);

if makeDiagnostics
    centeredResidual = angle(exp(1i * ...
        (residualColumn - globalPhaseRad)));
    diagnostics = table(caseColumn, windowColumn, timeColumn, ...
        frequencyColumn, tarcColumn, residualColumn, centeredResidual, ...
        repmat(offsetSec, numel(caseColumn), 1), ...
        'VariableNames', {'case_label', 'window_id', 'btt_time_center_s', ...
        'frequency_hz', 'TARC', 'phase_residual_rad', ...
        'centered_phase_residual_rad', 'offset_s'});
else
    diagnostics = table();
end
end

function diagnostics = amplitude_diagnostics_local(phaseModels, offsetSec, ...
    rawTime, rawMm, sampleDt, offsetRole)
nCase = numel(phaseModels);
caseLabel = strings(nCase, 1);
meanBttAmplitudeMm = nan(nCase, 1);
meanStrainAmplitudeMm = nan(nCase, 1);
amplitudeMapePercent = nan(nCase, 1);
amplitudeCorrelation = nan(nCase, 1);

for iCase = 1:nCase
    models = phaseModels{iCase};
    bttAmplitudeMm = [models.bttAmplitudeMm].';
    strainAmplitudeMm = nan(numel(models), 1);
    for iWindow = 1:numel(models)
        rawWindow = models(iWindow).timeWindow - offsetSec;
        index = time_window_indices_local(rawTime, sampleDt, rawWindow);
        if numel(index) < 20
            continue;
        end
        t = rawTime(index);
        y = detrend(rawMm(index));
        omegaT = 2 * pi * models(iWindow).frequencyHz .* t;
        coefficient = [cos(omegaT), sin(omegaT), ones(size(t))] \ y;
        strainAmplitudeMm(iWindow) = hypot(coefficient(1), coefficient(2));
    end
    good = isfinite(bttAmplitudeMm) & isfinite(strainAmplitudeMm) & ...
        strainAmplitudeMm > 0;
    caseLabel(iCase) = models(1).caseLabel;
    meanBttAmplitudeMm(iCase) = mean(bttAmplitudeMm(good));
    meanStrainAmplitudeMm(iCase) = mean(strainAmplitudeMm(good));
    amplitudeMapePercent(iCase) = 100 * mean(abs( ...
        bttAmplitudeMm(good) - strainAmplitudeMm(good)) ./ ...
        strainAmplitudeMm(good));
    amplitudeCorrelation(iCase) = signed_correlation_local( ...
        bttAmplitudeMm(good), strainAmplitudeMm(good));
end

diagnostics = table(repmat(string(offsetRole), nCase, 1), caseLabel, ...
    repmat(offsetSec, nCase, 1), meanBttAmplitudeMm, ...
    meanStrainAmplitudeMm, amplitudeMapePercent, amplitudeCorrelation, ...
    'VariableNames', {'offset_role', 'case_label', 'offset_s', ...
    'mean_btt_amplitude_mm', 'mean_strain_amplitude_mm', ...
    'amplitude_mape_percent', 'amplitude_correlation'});
end

function index = time_window_indices_local(time, sampleDt, timeWindow)
startIndex = max(1, ceil((timeWindow(1) - time(1)) / sampleDt) + 1);
endIndex = min(numel(time), floor((timeWindow(2) - time(1)) / sampleDt) + 1);
if endIndex < startIndex
    index = [];
else
    index = (startIndex:endIndex).';
end
end

function value = squared_correlation_local(a, b)
a = a(:) - mean(a, 'omitnan');
b = b(:) - mean(b, 'omitnan');
good = isfinite(a) & isfinite(b);
a = a(good);
b = b(good);
denominator = (a' * a) * (b' * b);
if numel(a) < 3 || denominator <= 0
    value = NaN;
else
    value = (a' * b) ^ 2 / denominator;
end
end

function value = signed_correlation_local(a, b)
a = a(:) - mean(a, 'omitnan');
b = b(:) - mean(b, 'omitnan');
good = isfinite(a) & isfinite(b);
a = a(good);
b = b(good);
denominator = sqrt((a' * a) * (b' * b));
if numel(a) < 3 || denominator <= 0
    value = NaN;
else
    value = (a' * b) / denominator;
end
end

function scanTable = make_scan_table_local(offsetGrid, priorOffsetSec, ...
    frequencyMetrics, phaseMetrics)
scanTable = table(offsetGrid(:), 1000 * (offsetGrid(:) - priorOffsetSec), ...
    frequencyMetrics.rmseHz, frequencyMetrics.correlation, ...
    frequencyMetrics.meanBiasHz, frequencyMetrics.sampleCount, ...
    phaseMetrics.combinedConcentration, ...
    phaseMetrics.caseConcentration(:, 1), ...
    phaseMetrics.caseConcentration(:, 2), ...
    phaseMetrics.medianTarc(:, 1), phaseMetrics.medianTarc(:, 2), ...
    'VariableNames', {'offset_s', 'correction_from_prior_ms', ...
    'frequency_rmse_hz', 'frequency_correlation', ...
    'frequency_mean_bias_hz', 'frequency_sample_count', ...
    'global_phase_concentration', 'phase_concentration_75', ...
    'phase_concentration_86', 'median_tarc_75', 'median_tarc_86'});
end

function plot_alignment_local(ridge, oprTime, orderFrequencyHz, bttSpanSec, ...
    offsetGrid, frequencyMetrics, phaseMetrics, priorOffsetSec, ...
    coarseOffsetSec, appliedOffsetSec, windowDiagnostics, ...
    estimatedMedianTarc, figureFile)
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18.3, 14.5], 'NumberTitle', 'off', ...
    'Name', 'Step12 global strain-BTT alignment');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
bttTime = ridge.rawCenterTimeSec + appliedOffsetSec;
targetFrequencyHz = interp1(oprTime, orderFrequencyHz, bttTime, ...
    'linear', NaN);
good = bttTime >= bttSpanSec(1) & bttTime <= bttSpanSec(2) & ...
    isfinite(ridge.peakFrequencyHz) & isfinite(targetFrequencyHz);
plot(ax, bttTime(good), targetFrequencyHz(good), 'k-', 'LineWidth', 1.2);
hold(ax, 'on');
plot(ax, bttTime(good), ridge.peakFrequencyHz(good), '-', ...
    'Color', [0.10, 0.35, 0.70], 'LineWidth', 1.1);
xlabel(ax, 'BTT time (s)');
ylabel(ax, 'Frequency (Hz)');
title(ax, 'Variable-speed coarse alignment');
legend(ax, {'OPR 12EO', 'Strain ridge'}, 'Location', 'best');
grid(ax, 'on');

ax = nexttile;
correctionMs = 1000 * (offsetGrid - priorOffsetSec);
plot(ax, correctionMs, frequencyMetrics.rmseHz, ...
    'Color', [0.10, 0.35, 0.70], 'LineWidth', 1.1);
hold(ax, 'on');
xline(ax, 1000 * (coarseOffsetSec - priorOffsetSec), '--', ...
    'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0);
xline(ax, 1000 * (appliedOffsetSec - priorOffsetSec), ':', ...
    'Color', [0.15, 0.55, 0.25], 'LineWidth', 1.2);
xlabel(ax, 'Offset correction from -102.6 s (ms)');
ylabel(ax, 'Frequency RMSE (Hz)');
title(ax, 'Frequency objective');
grid(ax, 'on');

ax = nexttile;
plot(ax, correctionMs, phaseMetrics.combinedConcentration, ...
    'Color', [0.45, 0.20, 0.60], 'LineWidth', 1.2);
hold(ax, 'on');
xline(ax, 1000 * (coarseOffsetSec - priorOffsetSec), '--', ...
    'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0);
xline(ax, 1000 * (appliedOffsetSec - priorOffsetSec), ':', ...
    'Color', [0.15, 0.55, 0.25], 'LineWidth', 1.2);
xlabel(ax, 'Offset correction from -102.6 s (ms)');
ylabel(ax, 'Global phase concentration');
title(ax, 'Phase-consistent TARC family selection');
ylim(ax, [0, 1.02]);
grid(ax, 'on');

ax = nexttile;
caseLabels = unique(windowDiagnostics.case_label, 'stable');
colors = [0.10, 0.35, 0.70; 0.85, 0.33, 0.10];
hold(ax, 'on');
for iCase = 1:numel(caseLabels)
    rows = windowDiagnostics.case_label == caseLabels(iCase);
    plot(ax, windowDiagnostics.btt_time_center_s(rows), ...
        180 / pi * windowDiagnostics.centered_phase_residual_rad(rows), ...
        'o-', 'Color', colors(iCase, :), 'MarkerSize', 3.5, ...
        'LineWidth', 0.8, 'DisplayName', caseLabels(iCase));
end
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'BTT time (s)');
ylabel(ax, 'Centered phase residual (deg)');
title(ax, sprintf('Applied offset; median TARC %.3f / %.3f', ...
    estimatedMedianTarc(1), estimatedMedianTarc(2)));
legend(ax, 'Location', 'best');
grid(ax, 'on');

exportgraphics(fig, figureFile, 'Resolution', 300);
close(fig);
end
