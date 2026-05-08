%% Order/revolution sensitivity comparison on the same simulated waveform
% Scenarios:
%   integer order:     500/1300 Hz at 3000 RPM = 10/26 EO
%   non-integer order: 520/1270 Hz at 3000 RPM = 10.4/25.4 EO
% Methods:
%   A:  low-template guided blind high-speed template refinement
%   A4: low-template frequency locking + post-hoc template transfer
%   A5: alias-expanded warped low-template full-waveform VARPRO
%   A2: constrained-template profile VARPRO
%   A3: physics-guided warped-template profile VARPRO
%   B:  gap-library raw scan + joint refinement
%   C:  same-gap direct fixed-template fit
clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'three_method_comparison_results');
if ~isfolder(outDir)
    mkdir(outDir);
end

cfg = make_default_config(rootDir, outDir);
[gapList, xCell, yCell] = load_stacked_gap_curves(cfg.dataFile);
templateLib = build_gap_template_library(gapList, xCell, yCell, cfg.g_holdout, cfg.xGridN);
templateLibFull = build_gap_template_library(gapList, xCell, yCell, NaN, cfg.xGridN);

nMethod = 7;
nScenario = numel(cfg.orderScenario);
nRev = numel(cfg.numRevsList);
targetSnr = cfg.targetSNRdB;
records = repmat(make_empty_record(), nScenario * nRev * cfg.numRepeat * nMethod, 1);
snrCal = repmat(make_empty_snr_cal_record(), nScenario * nRev, 1);
recId = 0;
calId = 0;

for iScenario = 1:nScenario
    for iRev = 1:nRev
        cfgCase = cfg;
        cfgCase.f_true = cfg.orderScenario(iScenario).f_true;
        cfgCase.NumRevs_high = cfg.numRevsList(iRev);
        [noiseRatioBySnr, cal] = calibrate_full_waveform_snr(cfgCase, gapList, xCell, yCell, templateLibFull.domain);
        noiseRatio = noiseRatioBySnr(1);
        calId = calId + 1;
        snrCal(calId).scenario = string(cfg.orderScenario(iScenario).name);
        snrCal(calId).num_revs_high = cfgCase.NumRevs_high;
        snrCal(calId).target_snr_db = targetSnr;
        snrCal(calId).noise_ratio = noiseRatio;
        snrCal(calId).signal_rms = cal.signalRmsAll;
        snrCal(calId).template_range = cal.templateRange;

        for iRep = 1:cfg.numRepeat
        seed = cfg.baseSeed + 100000 * iScenario + 1000 * iRev + iRep;
        rng(seed, 'twister');
        [Data_High, truth] = simulate_common_highspeed_data(cfgCase, gapList, xCell, yCell, templateLibFull.domain, noiseRatio);
        measuredSnr = measure_full_waveform_snr(Data_High);
        highMap = map_highspeed_to_space(Data_High, cfgCase.alpha_k, cfgCase.R_tip, templateLibFull.domain, cfgCase.fitActiveLevel);
        cloud = make_cloud_from_highmap(highMap, templateLibFull.domain, cfgCase.profile.xGridN);

        tA = tic;
        resultA = run_blind_transfer_method(cloud, gapList, xCell, yCell, cfgCase);
        elapsedA = toc(tA);

        tA4 = tic;
        resultA4 = run_blind_transfer_sensorcv_method(cloud, gapList, xCell, yCell, cfgCase);
        elapsedA4 = toc(tA4);

        tA5 = tic;
        resultA5 = run_a5_warp_only_varpro(cloud, gapList, xCell, yCell, cfgCase);
        elapsedA5 = toc(tA5);

        tA2 = tic;
        resultA2 = run_a2_constrained_varpro(cloud, gapList, xCell, yCell, cfgCase);
        elapsedA2 = toc(tA2);

        tA3 = tic;
        resultA3 = run_a3_physics_template_varpro(cloud, gapList, xCell, yCell, cfgCase);
        elapsedA3 = toc(tA3);

        tB = tic;
        staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfgCase);
        rawScanResult = run_varpro_for_gap(highMap, templateLib, staticState.gHat, cfgCase, 'gap_library_raw');
        resultB = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfgCase, cfgCase.jointRefine);
        elapsedB = toc(tB);

        tC = tic;
        resultC = run_varpro_for_gap(highMap, templateLibFull, cfgCase.g_holdout, cfgCase, 'same_gap_direct');
        elapsedC = toc(tC);

        recId = add_method_record(records, recId, cfgCase, truth, resultA, 'A_blind_transfer', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedA, NaN, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultA4, 'A4_low_locked_transfer', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedA4, NaN, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultA5, 'A5_warp_only_alias', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedA5, NaN, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultA2, 'A2_constrained_varpro', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedA2, NaN, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultA3, 'A3_physics_template', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedA3, NaN, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultB, 'B_gap_library_joint', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedB, resultB.g_used, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        recId = add_method_record(records, recId, cfgCase, truth, resultC, 'C_same_gap_direct', ...
            targetSnr, measuredSnr, noiseRatio, iRep, seed, elapsedC, cfgCase.g_holdout, ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high);
        records = recId.records;
        recId = recId.id;

        fprintf('[OrderRev] %s | %2d rev | SNR %.1f dB | rep %d/%d | A %.4g Hz, A4 %.4g Hz, A5 %.4g Hz, A2 %.4g Hz, A3 %.4g Hz, B %.4g Hz, C %.4g Hz\n', ...
            cfg.orderScenario(iScenario).name, cfgCase.NumRevs_high, targetSnr, iRep, cfg.numRepeat, ...
            mean(abs(sort(resultA.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultA4.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultA5.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultA2.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultA3.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultB.f_id(:).') - sort(truth.f(:).'))), ...
            mean(abs(sort(resultC.f_id(:).') - sort(truth.f(:).'))));
        end
    end
end

records = records(1:recId);
snrCal = snrCal(1:calId);
trialTable = struct2table(records);
summaryTable = summarize_comparison(trialTable, cfg.failureFreqErrorHz);
resultPrefix = sprintf('order_rev_%s', cfg.studyTag);
save(fullfile(outDir, [resultPrefix, '_results.mat']), ...
    'cfg', 'trialTable', 'summaryTable', 'snrCal');
writetable(trialTable, fullfile(outDir, [resultPrefix, '_trials.csv']));
writetable(summaryTable, fullfile(outDir, [resultPrefix, '_summary.csv']));
make_comparison_figure(summaryTable, outDir, resultPrefix);

fprintf('\nOrder/revolution sensitivity summary:\n');
disp(summaryTable);
fprintf('Saved comparison results to: %s\n', outDir);

function cfg = make_default_config(rootDir, outDir)
cfg = struct();
cfg.rootDir = rootDir;
cfg.outDir = outDir;
cfg.studyTag = 'order_rev_sensitivity_10dB_A5_warponly';
txtCandidates = dir(fullfile(rootDir, '**', '*2mm*.txt'));
txtNames = lower(string({txtCandidates.name}));
isTarget = ~startsWith(txtNames, "wave");
targetIdx = find(isTarget, 1, 'first');
if isempty(targetIdx)
    error('Could not find the stacked straight-blade 2 mm gap data file under: %s', rootDir);
end
cfg.dataFile = fullfile(txtCandidates(targetIdx).folder, txtCandidates(targetIdx).name);
cfg.g_low = 0.8;
cfg.g_holdout = 1.0;
cfg.R_tip = 62;
cfg.fs = 2e5;
cfg.RPM_high = 3000;
cfg.NumRevs_high = 8;
cfg.alpha_k = deg2rad([55, 116.76, 178.52]);
cfg.A_true = [0.25, 0.15];
cfg.f_true = [500, 1300];
cfg.phi_true = [pi/4, -pi/3];
cfg.targetSNRdB = 10;
cfg.numRepeat = 3;
cfg.numRevsList = [2, 4, 8, 20];
cfg.orderScenario = struct( ...
    'name', {"integer_500_1300", "noninteger_520_1270"}, ...
    'f_true', {[500, 1300], [520, 1270]});
cfg.baseSeed = 20260430;
cfg.failureFreqErrorHz = 1.0;
cfg.freqBoundaryTolHz = 5.0;
cfg.xGridN = 1201;
cfg.fitActiveLevel = 0.08;
cfg.gapQueryN = 121;
cfg.dxStaticGrid = linspace(-0.6, 0.6, 81);
cfg.rawGapSearchMargin = 0.25;
cfg.f1Grid = 400:5:650;
cfg.f2Grid = 1100:5:1400;
cfg.numVarproCandidates = 8;
cfg.jointRefine.gHalfWidth = 0.05;
cfg.jointRefine.dxHalfWidth = 0.05;
cfg.jointRefine.fHalfWidth = 20;
cfg.jointRefine.ampScaleLo = 0.50;
cfg.jointRefine.ampScaleHi = 1.50;
cfg.profile.xGridN = 1201;
cfg.profile.smoothWindow = 13;
cfg.profile.maxOuterIter = 3;
cfg.profile.tolParam = 1e-5;
cfg.profile.tolTemplate = 1e-5;
cfg.profile.derivativeFloorFrac = 0.03;
cfg.profile.f1Grid = cfg.f1Grid;
cfg.profile.f2Grid = cfg.f2Grid;
cfg.profile.minFreqSeparation = 20;
cfg.profile.ampUpperBound = 0.8;
cfg.profile.maxLsqIter = 300;
cfg.profile.lowInitKeep = 6;
cfg.profile.lowInitRefineCount = 14;
cfg.profile.lowInitUseGainBias = true;
cfg.profile.lowInitAmpScales = 1.0;
cfg.profile.lowInitFreqJitters = 0;
cfg.profile.useLowToHighTemplateTransfer = true;
cfg.profile.transferUseGainBias = true;
cfg.profile.transferResidualSmoothWindow = 21;
cfg.profile.useLocalRefineBounds = true;
cfg.profile.localFreqHalfWidth = 20;
cfg.profile.localAmpScaleLo = 0.60;
cfg.profile.localAmpScaleHi = 1.40;
cfg.profile.edgeGuardHz = 8;
cfg.profile.edgePenaltyWeight = 0.45;
cfg.profile.fitWeightFloor = 0.15;
cfg.profile.fitWeightPower = 1.0;
cfg.profile.transferResidualClipFrac = 0.45;
cfg.profile.selectionCvWeight = 0.65;
cfg.profile.cvCandidateCount = Inf;
cfg.profile.sensorCvWeight = 0.90;
cfg.profile.sensorCvCandidateCount = Inf;
cfg.profile.sensorCvMaxAliasCandidates = 18;
cfg.profile.sensorCvPreselectCount = 4;
cfg.profile.sensorCvRefineCount = 3;
cfg.profile.roughnessWeight = 0.03;
cfg.profile.aliasExpandTurnFreqMultiple = [-1, 0, 1];
cfg.a2.numShapeBasis = 10;
cfg.a2.ridgeLambda = 1e-4;
cfg.a2.numInitCandidates = 8;
cfg.a2.numRefineCandidates = 5;
cfg.a2.localFreqHalfWidth = 20;
cfg.a2.maxLsqIter = 250;
cfg.a2.ampUpperBound = 0.8;
cfg.a2.minValidFraction = 0.75;
cfg.a3.numShapeBasis = 8;
cfg.a3.ridgeLambda = 2e-4;
cfg.a3.numRefineCandidates = 5;
cfg.a3.localFreqHalfWidth = 20;
cfg.a3.maxLsqIter = 250;
cfg.a3.ampUpperBound = 0.8;
cfg.a3.dHalfWidth = 0.25;
cfg.a3.scaleLo = 0.90;
cfg.a3.scaleHi = 1.10;
cfg.a3.shiftPenalty = 1e-3;
cfg.a3.scalePenalty = 1e-2;
cfg.a3.minValidFraction = 0.75;
cfg.a5 = cfg.a3;
cfg.a5.numShapeBasis = 0;
cfg.a5.ridgeLambda = 0;
cfg.a5.numRefineCandidates = 10;
cfg.a5.localFreqHalfWidth = 18;
cfg.a5.dHalfWidth = 0.45;
cfg.a5.scaleLo = 0.82;
cfg.a5.scaleHi = 1.18;
cfg.a5.shiftPenalty = 5e-4;
cfg.a5.scalePenalty = 5e-3;
end

function [gapList, xCell, yCell] = load_stacked_gap_curves(filePath)
data = readmatrix(filePath, 'FileType', 'text', 'CommentStyle', '%');
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
x = data(:, 1);
y = data(:, 2);
breakIdx = [find(diff(x) < 0); numel(x)];
startIdx = [1; breakIdx(1:end-1) + 1];
nCurve = numel(breakIdx);
if nCurve == 7
    gapList = (0.2:0.2:1.4)';
else
    gapList = (1:nCurve)';
end
xCell = cell(nCurve, 1);
yCell = cell(nCurve, 1);
for i = 1:nCurve
    idx = startIdx(i):breakIdx(i);
    [xUnique, ia] = unique(x(idx), 'stable');
    yUnique = y(idx);
    xCell{i} = xUnique(:);
    yCell{i} = yUnique(ia);
end
end

function templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN)
if isnan(gHoldout)
    trainIdx = true(size(gapList));
else
    trainIdx = abs(gapList - gHoldout) > 1e-12;
end
trainIds = find(trainIdx);
G_train = gapList(trainIds);
xMin = -inf;
xMax = inf;
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    xMin = max(xMin, min(xCell{ig}));
    xMax = min(xMax, max(xCell{ig}));
end
xPad = 0.05 * (xMax - xMin);
xGrid = linspace(xMin + xPad, xMax - xPad, xGridN)';
S = zeros(numel(G_train), numel(xGrid));
for ii = 1:numel(trainIds)
    ig = trainIds(ii);
    S(ii, :) = interp1(xCell{ig}, yCell{ig}, xGrid, 'pchip');
end
dx = mean(diff(xGrid));
FxMat = zeros(size(S));
for ii = 1:size(S, 1)
    FxMat(ii, :) = gradient(S(ii, :), dx);
end
templateLib = struct('gapTrain', G_train(:), 'xGrid', xGrid(:), ...
    'S', S, 'FxMat', FxMat, 'domain', [min(xGrid), max(xGrid)]);
end

function [noiseRatioBySnr, cal] = calibrate_full_waveform_snr(cfg, gapList, xCell, yCell, domain)
[DataClean, ~] = simulate_common_highspeed_data(cfg, gapList, xCell, yCell, domain, 0);
signal = DataClean.V_cap_clean(:) - DataClean.baseline;
signalRms = sqrt(mean(signal.^2));
idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
probe = interp1(xCell{idxHoldout}, yCell{idxHoldout}, ...
    linspace(domain(1), domain(2), 400)', 'pchip');
templateRange = range(probe);
noiseRatioBySnr = signalRms ./ (templateRange .* 10.^(cfg.targetSNRdB(:) / 20));
cal = struct('signalRmsAll', signalRms, 'templateRange', templateRange);
end

function [Data_High, truth] = simulate_common_highspeed_data(cfg, gapList, xCell, yCell, domain, noiseRatio)
idxHoldout = find(abs(gapList - cfg.g_holdout) < 1e-12, 1);
Fx = griddedInterpolant(xCell{idxHoldout}, yCell{idxHoldout}, 'pchip', 'nearest');
Data_High = simulate_rotating_waveform_from_template(Fx, cfg.RPM_high, ...
    cfg.NumRevs_high, cfg.fs, cfg.R_tip, cfg.alpha_k, noiseRatio, ...
    domain, cfg.A_true, cfg.f_true, cfg.phi_true);
truth = struct('g_high', cfg.g_holdout, 'g_low', cfg.g_low, ...
    'A', cfg.A_true, 'f', cfg.f_true, 'phi', cfg.phi_true);
end

function Data = simulate_rotating_waveform_from_template(Fx, RPM, NumRevs, fs, R_tip, alpha_k, noiseRatio, domain, A, f, phi)
Omega = RPM * 2*pi / 60;
T = 2*pi / Omega;
dt = 1 / fs;
tEnd = (NumRevs + 1) * T;
t = (0:dt:tEnd)';
oprDelay = 0.08 * T;
T_opr = oprDelay + (0:NumRevs)' * T;
xHalf = 0.52 * diff(domain);
probe = linspace(domain(1), domain(2), 400)';
yProbe = Fx(probe);
baseline = min(yProbe);
V_clean = baseline * ones(size(t));
u_t = zeros(size(t));
for im = 1:numel(A)
    u_t = u_t + A(im) * sin(2*pi*f(im)*t + phi(im));
end
for m = 1:NumRevs
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = abs(t - Te) <= xHalf / (Omega * R_tip);
        xNom = Omega * R_tip * (t(idx) - Te);
        xPhys = xNom - u_t(idx);
        inDomain = xPhys >= domain(1) & xPhys <= domain(2);
        idxAll = find(idx);
        idxUse = idxAll(inDomain);
        if isempty(idxUse)
            continue;
        end
        V_clean(idxUse) = max(V_clean(idxUse), Fx(xPhys(inDomain)));
    end
end
noiseSigma = noiseRatio * max(range(yProbe), eps);
V_cap = V_clean + noiseSigma * randn(size(V_clean));
V_OPR = zeros(size(t));
oprWidth = 0.006 * T;
for m = 1:numel(T_opr)
    V_OPR(abs(t - T_opr(m)) <= oprWidth / 2) = 5.0;
end
V_OPR = V_OPR + 0.02 * randn(size(V_OPR));
Data = struct('t', t, 'V_cap', V_cap, 'V_cap_clean', V_clean, ...
    'V_OPR', V_OPR, 'RPM', RPM, 'T_opr_truth', T_opr(:), ...
    'u_t', u_t, 'baseline', baseline);
end

function snrDb = measure_full_waveform_snr(Data_High)
signal = Data_High.V_cap_clean(:) - Data_High.baseline;
noise = Data_High.V_cap(:) - Data_High.V_cap_clean(:);
snrDb = 20 * log10(sqrt(mean(signal.^2)) / sqrt(mean(noise.^2)));
end

function highMap = map_highspeed_to_space(Data_High, alpha_k, R_tip, domain, activeLevel)
T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
xHalf = 0.48 * diff(domain);
base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + activeLevel * span;
tAll = [];
vAll = [];
xAll = [];
revAll = [];
sAll = [];
for m = 1:numel(T_opr)-1
    Omega = 2*pi / (T_opr(m+1) - T_opr(m));
    for k = 1:numel(alpha_k)
        Te = T_opr(m) + alpha_k(k) / Omega;
        idx = find(abs(t - Te) <= xHalf / (Omega * R_tip));
        if isempty(idx), continue; end
        xLocal = Omega * R_tip * (t(idx) - Te);
        keep = xLocal >= domain(1) & xLocal <= domain(2) & V(idx) > activeThreshold;
        idx = idx(keep);
        xLocal = xLocal(keep);
        if isempty(idx), continue; end
        tAll = [tAll; t(idx)]; %#ok<AGROW>
        vAll = [vAll; V(idx)]; %#ok<AGROW>
        xAll = [xAll; xLocal]; %#ok<AGROW>
        revAll = [revAll; m * ones(numel(idx), 1)]; %#ok<AGROW>
        sAll = [sAll; k * ones(numel(idx), 1)]; %#ok<AGROW>
    end
end
highMap = struct('t_v', tAll, 'V_a', vAll, 'x_v', xAll, ...
    'rev_v', revAll, 'S_v', sAll, 'activeThreshold', activeThreshold);
end

function cloud = make_cloud_from_highmap(highMap, domain, xGridN)
cloud = struct();
cloud.t = highMap.t_v(:);
cloud.x = highMap.x_v(:);
cloud.V = highMap.V_a(:);
cloud.rev = highMap.rev_v(:);
cloud.sensor = highMap.S_v(:);
cloud.domain = domain;
cloud.xGrid = linspace(domain(1), domain(2), xGridN)';
end

function result = run_blind_transfer_method(cloud, gapList, xCell, yCell, cfg)
idxLow = find(abs(gapList(:) - cfg.g_low) < 1e-12, 1);
lowTemplate = struct('x', xCell{idxLow}(:), 'y', yCell{idxLow}(:), 'g', cfg.g_low);
xGrid = cloud.xGrid(:);
T0 = estimate_template_from_cloud(cloud.x, cloud.V, xGrid, cfg.profile);
templateAnchor = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
[pInitList, initInfo] = initial_candidates_by_low_template(cloud, lowTemplate, xGrid, cfg.profile);
pInitList = expand_low_template_starts(pInitList, cfg.profile);
candidateList = cell(size(pInitList, 1), 1);
baseScore = inf(size(pInitList, 1), 1);
for ic = 1:size(pInitList, 1)
    candidate = run_profile_outer_loop(cloud, xGrid, T0, pInitList(ic, :), cfg.profile, templateAnchor);
    candidate.cvRmse = NaN;
    candidate.selectionScore = candidate.rmse * frequency_penalty_multiplier(candidate.p, cfg.profile);
    candidateList{ic} = candidate;
    baseScore(ic) = candidate.selectionScore;
end
[~, candidateOrder] = sort(baseScore, 'ascend');
nCv = min(cfg.profile.cvCandidateCount, numel(candidateOrder));
for ii = 1:nCv
    ic = candidateOrder(ii);
    candidate = candidateList{ic};
    candidate.cvRmse = calc_template_cv_rmse(cloud, templateAnchor, xGrid, candidate.p, cfg.profile);
    candidate.selectionScore = blend_selection_rmse(candidate.rmse, candidate.cvRmse, cfg.profile) * ...
        frequency_penalty_multiplier(candidate.p, cfg.profile);
    candidateList{ic} = candidate;
end
selectionScore = cellfun(@(c) c.selectionScore, candidateList);
[~, bestIdx] = min(selectionScore);
result = candidateList{bestIdx};
result.label = 'blind_transfer';
result.g_used = NaN;
result.initInfo = initInfo;
result.f_id = result.f;
result.A_id = result.A;
end

function result = run_blind_transfer_sensorcv_method(cloud, gapList, xCell, yCell, cfg)
idxLow = find(abs(gapList(:) - cfg.g_low) < 1e-12, 1);
lowTemplate = struct('x', xCell{idxLow}(:), 'y', yCell{idxLow}(:), 'g', cfg.g_low);
xGrid = cloud.xGrid(:);
templateAnchor = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
[pInitList, initInfo] = initial_candidates_by_low_template(cloud, lowTemplate, xGrid, cfg.profile);
if isempty(pInitList)
    error('Low-template initialization did not return any candidate.');
end
p = pInitList(1, :);
T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, cfg.profile);
if p(3) > p(6)
    p = [p(4), p(5), p(6), p(1), p(2), p(3)];
end
[rmse, residual, VFit, xi] = calc_profile_rmse(cloud, T, xGrid, p);
result = struct('p', p, 'A', [p(1), p(4)], 'f', [p(3), p(6)], ...
    'phi', [p(2), p(5)], 'T', T(:), 'xGrid', xGrid(:), 'xi', xi(:), ...
    'VFit', VFit(:), 'residual', residual(:), 'rmse', rmse, ...
    'costHist', rmse, 'label', 'blind_transfer_low_locked');
result.g_used = NaN;
result.initInfo = initInfo;
result.f_id = result.f;
result.A_id = result.A;
end

function result = run_a2_constrained_varpro(cloud, gapList, xCell, yCell, cfg)
idxLow = find(abs(gapList(:) - cfg.g_low) < 1e-12, 1);
lowTemplate = struct('x', xCell{idxLow}(:), 'y', yCell{idxLow}(:));
xGrid = cloud.xGrid(:);
T0 = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
model = make_a2_model(xGrid, T0, cfg.a2);
[pCandidates, initInfo] = a2_initial_frequency_candidates(cloud, model, cfg);
lbWide = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid)];
ubWide = [cfg.a2.ampUpperBound, pi, max(cfg.f1Grid), cfg.a2.ampUpperBound, pi, max(cfg.f2Grid)];
best = struct('cost', inf, 'p', [], 'fit', []);
for ic = 1:min(cfg.a2.numRefineCandidates, size(pCandidates, 1))
    p0 = min(max(pCandidates(ic, :), lbWide), ubWide);
    lb = lbWide;
    ub = ubWide;
    lb(3) = max(lbWide(3), p0(3) - cfg.a2.localFreqHalfWidth);
    ub(3) = min(ubWide(3), p0(3) + cfg.a2.localFreqHalfWidth);
    lb(6) = max(lbWide(6), p0(6) - cfg.a2.localFreqHalfWidth);
    ub(6) = min(ubWide(6), p0(6) + cfg.a2.localFreqHalfWidth);
    resFun = @(p) a2_profile_residual(p, cloud, model, cfg.a2);
    pTry = run_lsqnonlin_or_fminsearch(resFun, p0, lb, ub, cfg.a2.maxLsqIter);
    fitTry = a2_fit_for_p(pTry, cloud, model, cfg.a2);
    costTry = sum(fitTry.profileResidual.^2);
    if costTry < best.cost
        best.cost = costTry;
        best.p = pTry;
        best.fit = fitTry;
    end
end
pOpt = best.p;
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
fit = a2_fit_for_p(pOpt, cloud, model, cfg.a2);
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
Bcloud = a2_static_design(model, cloud.x(:));
Dcloud = -model.dF0(cloud.x(:));
diagInfo = a2_dynamic_diagnostics(Bcloud, a2_dynamic_design(Dcloud, cloud.t(:), f_id(1), f_id(2)));
result = struct('label', 'A2-constrained-template-varpro', ...
    'p', pOpt, 'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(fit.dataResidual(fit.valid).^2)), ...
    'VFit', fit.VFit(:), 'residual', fit.dataResidual(:), ...
    'T', fit.Tgrid(:), 'xGrid', model.xGrid(:), 'thetaTemplate', fit.theta(:), ...
    'eta_dynamic', diagInfo.eta, 'cond_dynamic', diagInfo.condG, ...
    'initInfo', initInfo);
end

function model = make_a2_model(xGrid, T0, a2Cfg)
xGrid = xGrid(:);
T0 = T0(:);
dT0 = gradient(T0, mean(diff(xGrid)));
Phi = eval_a2_shape_basis(xGrid, [min(xGrid), max(xGrid)], a2Cfg.numShapeBasis);
model = struct();
model.xGrid = xGrid;
model.domain = [min(xGrid), max(xGrid)];
model.T0Grid = T0;
model.PhiGrid = Phi;
model.F0 = griddedInterpolant(xGrid, T0, 'pchip', 'none');
model.dF0 = griddedInterpolant(xGrid, dT0, 'pchip', 'nearest');
model.Bx = [model.F0(xGrid), ones(numel(xGrid), 1), Phi];
model.dT0x = model.dF0(xGrid);
end

function Phi = eval_a2_shape_basis(x, domain, nBasis)
x = x(:);
s = (x - domain(1)) ./ max(domain(2) - domain(1), eps);
s = min(max(s, 0), 1);
Phi = zeros(numel(x), nBasis);
for k = 1:nBasis
    Phi(:, k) = cos(pi * k * s);
end
end

function [pCandidates, initInfo] = a2_initial_frequency_candidates(cloud, model, cfg)
a2Cfg = cfg.a2;
x = cloud.x(:);
t = cloud.t(:);
y = cloud.V(:);
B = a2_static_design(model, x);
D = -model.dF0(x);
valid = all(isfinite(B), 2) & isfinite(D) & isfinite(y);
B = B(valid, :);
D = D(valid);
t = t(valid);
y = y(valid);
maxKeep = max(a2Cfg.numInitCandidates * 4, a2Cfg.numInitCandidates);
candCost = inf(maxKeep, 1);
candP = zeros(maxKeep, 6);
candEta = nan(maxKeep, 1);
candCond = nan(maxKeep, 1);
lambda = a2Cfg.ridgeLambda;
penalty = diag([0, 0, lambda * ones(1, a2Cfg.numShapeBasis), 1e-10 * ones(1, 4)]);
for f1 = cfg.f1Grid(:)'
    s1 = sin(2*pi*f1*t);
    c1 = cos(2*pi*f1*t);
    for f2 = cfg.f2Grid(:)'
        if abs(f1 - f2) < cfg.profile.minFreqSeparation
            continue;
        end
        G = [D .* s1, D .* c1, D .* sin(2*pi*f2*t), D .* cos(2*pi*f2*t)];
        X = [B, G];
        theta = (X' * X + penalty) \ (X' * y);
        res = y - X * theta;
        cost = sum(res.^2) + lambda * sum(theta(3:2+a2Cfg.numShapeBasis).^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            beta = theta(end-3:end);
            diagInfo = a2_dynamic_diagnostics(B, G);
            candCost(worstIdx) = cost;
            candP(worstIdx, :) = [hypot(beta(1), beta(2)), atan2(beta(2), beta(1)), f1, ...
                hypot(beta(3), beta(4)), atan2(beta(4), beta(3)), f2];
            candEta(worstIdx) = diagInfo.eta;
            candCond(worstIdx) = diagInfo.condG;
        end
    end
end
[candCost, order] = sort(candCost, 'ascend');
candP = candP(order, :);
candEta = candEta(order);
candCond = candCond(order);
keep = isfinite(candCost);
candCost = candCost(keep);
candP = candP(keep, :);
candEta = candEta(keep);
candCond = candCond(keep);
nOut = min(a2Cfg.numInitCandidates, size(candP, 1));
pCandidates = candP(1:nOut, :);
initInfo = struct('cost', candCost(1:nOut), 'eta', candEta(1:nOut), ...
    'condG', candCond(1:nOut), 'numValidLinearizedSamples', nnz(valid));
end

function r = a2_profile_residual(p, cloud, model, a2Cfg)
fit = a2_fit_for_p(p, cloud, model, a2Cfg);
r = fit.profileResidual;
end

function fit = a2_fit_for_p(p, cloud, model, a2Cfg)
xi = cloud.x(:) - displacement_2freq(p, cloud.t(:));
y = cloud.V(:);
B = a2_static_design(model, xi);
valid = all(isfinite(B), 2) & isfinite(y);
theta = zeros(2 + a2Cfg.numShapeBasis, 1);
if nnz(valid) >= a2Cfg.minValidFraction * numel(y)
    X = B(valid, :);
    yy = y(valid);
    penalty = diag([0, 0, a2Cfg.ridgeLambda * ones(1, a2Cfg.numShapeBasis)]);
    theta = (X' * X + penalty) \ (X' * yy);
    VFit = nan(size(y));
    VFit(valid) = X * theta;
    dataResidual = y - VFit;
else
    dataResidual = 5 * max(range(model.T0Grid), eps) * ones(size(y));
    VFit = y - dataResidual;
end
invalid = ~isfinite(dataResidual);
if any(invalid)
    dataResidual(invalid) = 5 * max(range(model.T0Grid), eps);
end
regResidual = sqrt(a2Cfg.ridgeLambda) * theta(3:end);
profileResidual = [dataResidual; regResidual];
Tgrid = model.Bx * theta;
fit = struct('theta', theta, 'VFit', VFit, 'dataResidual', dataResidual, ...
    'profileResidual', profileResidual, 'valid', valid, 'Tgrid', Tgrid);
end

function B = a2_static_design(model, x)
T0 = model.F0(x(:));
Phi = eval_a2_shape_basis(x(:), model.domain, size(model.PhiGrid, 2));
B = [T0, ones(numel(x), 1), Phi];
end

function G = a2_dynamic_design(D, t, f1, f2)
G = [D .* sin(2*pi*f1*t), D .* cos(2*pi*f1*t), ...
    D .* sin(2*pi*f2*t), D .* cos(2*pi*f2*t)];
end

function diagInfo = a2_dynamic_diagnostics(B, G)
valid = all(isfinite(B), 2) & all(isfinite(G), 2);
B = B(valid, :);
G = G(valid, :);
[Q, ~] = qr(B, 0);
Gperp = G - Q * (Q' * G);
diagInfo = struct();
diagInfo.eta = norm(Gperp, 'fro') / max(norm(G, 'fro'), eps);
s = svd(Gperp, 'econ');
if isempty(s) || min(s) <= eps
    diagInfo.condG = inf;
else
    diagInfo.condG = max(s) / min(s);
end
end

function result = run_a3_physics_template_varpro(cloud, gapList, xCell, yCell, cfg)
idxLow = find(abs(gapList(:) - cfg.g_low) < 1e-12, 1);
lowTemplate = struct('x', xCell{idxLow}(:), 'y', yCell{idxLow}(:));
xGrid = cloud.xGrid(:);
T0 = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
a2Model = make_a2_model(xGrid, T0, cfg.a2);
model = make_a3_model(xGrid, T0, cfg.a3);
[pCandidates, initInfo] = a2_initial_frequency_candidates(cloud, a2Model, cfg);
lbWide = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid), ...
    -cfg.a3.dHalfWidth, cfg.a3.scaleLo];
ubWide = [cfg.a3.ampUpperBound, pi, max(cfg.f1Grid), cfg.a3.ampUpperBound, pi, max(cfg.f2Grid), ...
    cfg.a3.dHalfWidth, cfg.a3.scaleHi];
best = struct('cost', inf, 'theta', [], 'fit', []);
for ic = 1:min(cfg.a3.numRefineCandidates, size(pCandidates, 1))
    p0 = pCandidates(ic, :);
    theta0 = [p0, 0, 1];
    lb = lbWide;
    ub = ubWide;
    lb(3) = max(lbWide(3), p0(3) - cfg.a3.localFreqHalfWidth);
    ub(3) = min(ubWide(3), p0(3) + cfg.a3.localFreqHalfWidth);
    lb(6) = max(lbWide(6), p0(6) - cfg.a3.localFreqHalfWidth);
    ub(6) = min(ubWide(6), p0(6) + cfg.a3.localFreqHalfWidth);
    resFun = @(theta) a3_profile_residual(theta, cloud, model, cfg.a3);
    thetaTry = run_lsqnonlin_or_fminsearch(resFun, theta0, lb, ub, cfg.a3.maxLsqIter);
    fitTry = a3_fit_for_theta(thetaTry, cloud, model, cfg.a3);
    costTry = sum(a3_profile_residual(thetaTry, cloud, model, cfg.a3).^2);
    if costTry < best.cost
        best.cost = costTry;
        best.theta = thetaTry;
        best.fit = fitTry;
    end
end
thetaOpt = best.theta;
thetaOpt(2) = wrap_to_pi_local(thetaOpt(2));
thetaOpt(5) = wrap_to_pi_local(thetaOpt(5));
fit = a3_fit_for_theta(thetaOpt, cloud, model, cfg.a3);
pOpt = thetaOpt(1:6);
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
Bcloud = a3_static_design(model, cloud.x(:), thetaOpt(7), thetaOpt(8));
Dcloud = -model.dF0((cloud.x(:) - thetaOpt(7)) ./ thetaOpt(8));
diagInfo = a2_dynamic_diagnostics(Bcloud, a2_dynamic_design(Dcloud, cloud.t(:), f_id(1), f_id(2)));
result = struct('label', 'A3-physics-guided-template-varpro', ...
    'p', pOpt, 'theta', thetaOpt, 'template_shift', thetaOpt(7), 'template_scale', thetaOpt(8), ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(fit.dataResidual(fit.valid).^2)), ...
    'VFit', fit.VFit(:), 'residual', fit.dataResidual(:), ...
    'T', fit.Tgrid(:), 'xGrid', model.xGrid(:), 'thetaTemplate', fit.thetaTemplate(:), ...
    'eta_dynamic', diagInfo.eta, 'cond_dynamic', diagInfo.condG, ...
    'initInfo', initInfo);
end

function result = run_a5_warp_only_varpro(cloud, gapList, xCell, yCell, cfg)
idxLow = find(abs(gapList(:) - cfg.g_low) < 1e-12, 1);
lowTemplate = struct('x', xCell{idxLow}(:), 'y', yCell{idxLow}(:));
xGrid = cloud.xGrid(:);
T0 = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
initCfg = cfg;
initCfg.a2.numShapeBasis = 0;
initCfg.a2.ridgeLambda = 0;
initModel = make_a2_model(xGrid, T0, initCfg.a2);
model = make_a3_model(xGrid, T0, cfg.a5);
[pCandidates, initInfo] = a2_initial_frequency_candidates(cloud, initModel, initCfg);
pCandidates = expand_order_alias_starts(pCandidates, cfg.profile, cfg.RPM_high / 60);
lbWide = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid), ...
    -cfg.a5.dHalfWidth, cfg.a5.scaleLo];
ubWide = [cfg.a5.ampUpperBound, pi, max(cfg.f1Grid), cfg.a5.ampUpperBound, pi, max(cfg.f2Grid), ...
    cfg.a5.dHalfWidth, cfg.a5.scaleHi];
best = struct('cost', inf, 'theta', [], 'fit', []);
for ic = 1:min(cfg.a5.numRefineCandidates, size(pCandidates, 1))
    p0 = min(max(pCandidates(ic, :), lbWide(1:6)), ubWide(1:6));
    theta0 = [p0, 0, 1];
    lb = lbWide;
    ub = ubWide;
    lb(3) = max(lbWide(3), p0(3) - cfg.a5.localFreqHalfWidth);
    ub(3) = min(ubWide(3), p0(3) + cfg.a5.localFreqHalfWidth);
    lb(6) = max(lbWide(6), p0(6) - cfg.a5.localFreqHalfWidth);
    ub(6) = min(ubWide(6), p0(6) + cfg.a5.localFreqHalfWidth);
    resFun = @(theta) a3_profile_residual(theta, cloud, model, cfg.a5);
    thetaTry = run_lsqnonlin_or_fminsearch(resFun, theta0, lb, ub, cfg.a5.maxLsqIter);
    fitTry = a3_fit_for_theta(thetaTry, cloud, model, cfg.a5);
    costTry = sum(a3_profile_residual(thetaTry, cloud, model, cfg.a5).^2);
    if costTry < best.cost
        best.cost = costTry;
        best.theta = thetaTry;
        best.fit = fitTry;
    end
end
thetaOpt = best.theta;
thetaOpt(2) = wrap_to_pi_local(thetaOpt(2));
thetaOpt(5) = wrap_to_pi_local(thetaOpt(5));
fit = a3_fit_for_theta(thetaOpt, cloud, model, cfg.a5);
pOpt = thetaOpt(1:6);
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
Bcloud = a3_static_design(model, cloud.x(:), thetaOpt(7), thetaOpt(8));
Dcloud = -model.dF0((cloud.x(:) - thetaOpt(7)) ./ thetaOpt(8));
diagInfo = a2_dynamic_diagnostics(Bcloud, a2_dynamic_design(Dcloud, cloud.t(:), f_id(1), f_id(2)));
result = struct('label', 'A5-warp-only-alias-varpro', ...
    'p', pOpt, 'theta', thetaOpt, 'template_shift', thetaOpt(7), 'template_scale', thetaOpt(8), ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(fit.dataResidual(fit.valid).^2)), ...
    'VFit', fit.VFit(:), 'residual', fit.dataResidual(:), ...
    'T', fit.Tgrid(:), 'xGrid', model.xGrid(:), 'thetaTemplate', fit.thetaTemplate(:), ...
    'eta_dynamic', diagInfo.eta, 'cond_dynamic', diagInfo.condG, ...
    'initInfo', initInfo);
end

function model = make_a3_model(xGrid, T0, a3Cfg)
xGrid = xGrid(:);
T0 = T0(:);
dT0 = gradient(T0, mean(diff(xGrid)));
model = struct();
model.xGrid = xGrid;
model.domain = [min(xGrid), max(xGrid)];
model.T0Grid = T0;
model.PhiGrid = eval_a2_shape_basis(xGrid, model.domain, a3Cfg.numShapeBasis);
model.F0 = griddedInterpolant(xGrid, T0, 'pchip', 'none');
model.dF0 = griddedInterpolant(xGrid, dT0, 'pchip', 'nearest');
end

function r = a3_profile_residual(theta, cloud, model, a3Cfg)
fit = a3_fit_for_theta(theta, cloud, model, a3Cfg);
r = fit.profileResidual;
end

function fit = a3_fit_for_theta(theta, cloud, model, a3Cfg)
p = theta(1:6);
d = theta(7);
s = theta(8);
xi = cloud.x(:) - displacement_2freq(p, cloud.t(:));
y = cloud.V(:);
B = a3_static_design(model, xi, d, s);
valid = all(isfinite(B), 2) & isfinite(y);
nTheta = 2 + a3Cfg.numShapeBasis;
thetaTemplate = zeros(nTheta, 1);
if nnz(valid) >= a3Cfg.minValidFraction * numel(y)
    X = B(valid, :);
    yy = y(valid);
    penalty = diag([0, 0, a3Cfg.ridgeLambda * ones(1, a3Cfg.numShapeBasis)]);
    thetaTemplate = (X' * X + penalty) \ (X' * yy);
    VFit = nan(size(y));
    VFit(valid) = X * thetaTemplate;
    dataResidual = y - VFit;
else
    dataResidual = 5 * max(range(model.T0Grid), eps) * ones(size(y));
    VFit = y - dataResidual;
end
invalid = ~isfinite(dataResidual);
if any(invalid)
    dataResidual(invalid) = 5 * max(range(model.T0Grid), eps);
end
regResidual = sqrt(a3Cfg.ridgeLambda) * thetaTemplate(3:end);
shiftResidual = sqrt(a3Cfg.shiftPenalty) * d;
scaleResidual = sqrt(a3Cfg.scalePenalty) * (s - 1);
profileResidual = [dataResidual; regResidual; shiftResidual; scaleResidual];
Tgrid = a3_static_design(model, model.xGrid, d, s) * thetaTemplate;
fit = struct('thetaTemplate', thetaTemplate, 'VFit', VFit, 'dataResidual', dataResidual, ...
    'profileResidual', profileResidual, 'valid', valid, 'Tgrid', Tgrid);
end

function B = a3_static_design(model, x, d, s)
x = x(:);
xWarp = (x - d) ./ max(s, eps);
T0Warp = model.F0(xWarp);
Phi = eval_a2_shape_basis(x, model.domain, size(model.PhiGrid, 2));
B = [T0Warp, ones(numel(x), 1), Phi];
end

function [pList, initInfo] = initial_candidates_by_low_template(cloud, lowTemplate, xGrid, profileCfg)
Tlow = interp1(lowTemplate.x, lowTemplate.y, xGrid, 'pchip', 'extrap');
TAtX = interp1(xGrid, Tlow, cloud.x, 'pchip', 'extrap');
dT = gradient(Tlow(:), mean(diff(xGrid)));
dTAtX = interp1(xGrid, dT, cloud.x, 'pchip', 'extrap');
valid = isfinite(TAtX) & isfinite(dTAtX) & isfinite(cloud.V);
derivScale = max(abs(dTAtX(valid)));
valid = valid & abs(dTAtX) > profileCfg.derivativeFloorFrac * max(derivScale, eps);
theta = [TAtX(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
gain0 = theta(1);
bias0 = theta(2);
f1Grid = profileCfg.f1Grid(:)';
f2Grid = profileCfg.f2Grid(:)';
keepCount = profileCfg.lowInitKeep;
maxKeep = max(keepCount * 6, keepCount);
candCost = inf(maxKeep, 1);
candP = zeros(maxKeep, 6);
lambda = 1e-9 * max(nnz(valid), 1);
t = cloud.t(valid);
y = cloud.V(valid) - (gain0 * TAtX(valid) + bias0);
kx = gain0 * dTAtX(valid);
for f1 = f1Grid
    s1 = sin(2*pi*f1*t);
    c1 = cos(2*pi*f1*t);
    for f2 = f2Grid
        if abs(f1 - f2) < profileCfg.minFreqSeparation
            continue;
        end
        s2 = sin(2*pi*f2*t);
        c2 = cos(2*pi*f2*t);
        H = [-kx .* s1, -kx .* c1, -kx .* s2, -kx .* c2];
        beta = (H' * H + lambda * eye(4)) \ (H' * y);
        res = y - H * beta;
        cost = sum(res.^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            candCost(worstIdx) = cost;
            candP(worstIdx, :) = [hypot(beta(1), beta(2)), atan2(beta(2), beta(1)), f1, ...
                hypot(beta(3), beta(4)), atan2(beta(4), beta(3)), f2];
        end
    end
end
[candCost, order] = sort(candCost, 'ascend');
candP = candP(order, :);
keep = isfinite(candCost);
candP = candP(keep, :);
lb = [0, -pi, min(f1Grid), 0, -pi, min(f2Grid)];
ub = [profileCfg.ampUpperBound, pi, max(f1Grid), profileCfg.ampUpperBound, pi, max(f2Grid)];
numRaw = min(profileCfg.lowInitRefineCount, size(candP, 1));
refinedP = zeros(numRaw, 6);
refinedCost = inf(numRaw, 1);
for i = 1:numRaw
    p0 = min(max(candP(i, :), lb), ub);
    pOpt = refine_p_given_low_template(cloud, Tlow, xGrid, p0, lb, ub, profileCfg);
    refinedP(i, :) = pOpt;
    refinedCost(i) = calc_low_template_rmse(cloud, Tlow, xGrid, pOpt, profileCfg);
end
[refinedCost, order] = sort(refinedCost, 'ascend');
refinedP = refinedP(order, :);
refinedScore = refinedCost .* frequency_penalty_multiplier(refinedP, profileCfg);
[~, scoreOrder] = sort(refinedScore, 'ascend');
refinedP = refinedP(scoreOrder, :);
refinedCost = refinedCost(scoreOrder);
nOut = min(keepCount, size(refinedP, 1));
pList = refinedP(1:nOut, :);
initInfo = struct('lowTemplateGain0', gain0, 'lowTemplateBias0', bias0, ...
    'lowTemplateRmse', refinedCost(1:nOut), 'numValidLinearizedSamples', nnz(valid));
end

function pStart = expand_low_template_starts(pInitList, profileCfg)
ampScales = profileCfg.lowInitAmpScales;
freqJitters = profileCfg.lowInitFreqJitters;
numStarts = profileCfg.lowInitKeep;
starts = [];
baseCount = min(2, size(pInitList, 1));
for ib = 1:baseCount
    p = pInitList(ib, :);
    for ia = 1:numel(ampScales)
        for jf1 = 1:numel(freqJitters)
            for jf2 = 1:numel(freqJitters)
                ps = p;
                ps(1) = min(profileCfg.ampUpperBound, max(0, ampScales(ia) * p(1)));
                ps(4) = min(profileCfg.ampUpperBound, max(0, ampScales(ia) * p(4)));
                ps(3) = min(max(p(3) + freqJitters(jf1), min(profileCfg.f1Grid)), max(profileCfg.f1Grid));
                ps(6) = min(max(p(6) + freqJitters(jf2), min(profileCfg.f2Grid)), max(profileCfg.f2Grid));
                starts = [starts; ps]; %#ok<AGROW>
            end
        end
    end
end
starts = [pInitList; starts]; %#ok<AGROW>
key = round([starts(:, 1), starts(:, 2), starts(:, 3), starts(:, 4), starts(:, 5), starts(:, 6)] * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
starts = starts(ia, :);
pStart = starts(1:min(numStarts, size(starts, 1)), :);
end

function pOut = expand_order_alias_starts(pIn, profileCfg, rotFreqHz)
if isempty(pIn)
    pOut = pIn;
    return;
end
aliasMultiples = profileCfg.aliasExpandTurnFreqMultiple(:).';
starts = [];
for i = 1:size(pIn, 1)
    p = pIn(i, :);
    for m1 = aliasMultiples
        for m2 = aliasMultiples
            q = p;
            q(3) = p(3) + m1 * rotFreqHz;
            q(6) = p(6) + m2 * rotFreqHz;
            if q(3) < min(profileCfg.f1Grid) || q(3) > max(profileCfg.f1Grid)
                continue;
            end
            if q(6) < min(profileCfg.f2Grid) || q(6) > max(profileCfg.f2Grid)
                continue;
            end
            if abs(q(6) - q(3)) < profileCfg.minFreqSeparation
                continue;
            end
            if q(3) > q(6)
                q = [q(4), q(5), q(6), q(1), q(2), q(3)];
            end
            starts = [starts; q]; %#ok<AGROW>
        end
    end
end
if isempty(starts)
    pOut = pIn;
    return;
end
key = round(starts * 1000) / 1000;
[~, ia] = unique(key, 'rows', 'stable');
pOut = starts(ia, :);
end

function pOpt = refine_p_given_low_template(cloud, Tlow, xGrid, p0, lb, ub, profileCfg)
p0 = min(max(p0, lb), ub);
resFun = @(p) low_template_residual(p, cloud, Tlow, xGrid, profileCfg);
pOpt = run_lsqnonlin_or_fminsearch(resFun, p0, lb, ub, profileCfg.maxLsqIter);
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function rmse = calc_low_template_rmse(cloud, Tlow, xGrid, p, profileCfg)
r = low_template_residual(p, cloud, Tlow, xGrid, profileCfg);
rmse = sqrt(mean(r.^2));
end

function r = low_template_residual(p, cloud, Tlow, xGrid, profileCfg)
xi = cloud.x - displacement_2freq(p, cloud.t);
Tfit = interp1(xGrid, Tlow, xi, 'pchip', NaN);
valid = isfinite(Tfit) & isfinite(cloud.V);
if nnz(valid) < 0.75 * numel(cloud.V)
    templateRange = max(Tlow) - min(Tlow);
    r = 5 * max(templateRange, eps) * ones(size(cloud.V));
    return;
end
if profileCfg.lowInitUseGainBias
    theta = [Tfit(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
    Vpred = theta(1) * Tfit + theta(2);
else
    Vpred = Tfit;
end
r = cloud.V - Vpred;
invalid = ~isfinite(r);
if any(invalid)
    templateRange = max(Tlow) - min(Tlow);
    r(invalid) = 5 * max(templateRange, eps);
end
end

function candidate = run_profile_outer_loop(cloud, xGrid, T0, p0, profileCfg, templateAnchor)
p = p0(:).';
T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, profileCfg);
if any(~isfinite(T))
    T = T0(:);
end
[lb, ub] = make_profile_refine_bounds(p, profileCfg);
bestLoop = struct('rmse', inf, 'p', p, 'T', T);
costHist = zeros(profileCfg.maxOuterIter, 1);
for iter = 1:profileCfg.maxOuterIter
    p = refine_p_given_template(cloud, T, xGrid, p, lb, ub, profileCfg);
    T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, profileCfg);
    rmse = calc_profile_rmse(cloud, T, xGrid, p);
    costHist(iter) = rmse;
    if rmse < bestLoop.rmse
        bestLoop.rmse = rmse;
        bestLoop.p = p;
        bestLoop.T = T;
    end
end
p = bestLoop.p;
T = bestLoop.T;
if p(3) > p(6)
    p = [p(4), p(5), p(6), p(1), p(2), p(3)];
end
[rmse, residual, VFit, xi] = calc_profile_rmse(cloud, T, xGrid, p);
candidate = struct('p', p, 'A', [p(1), p(4)], 'f', [p(3), p(6)], ...
    'phi', [p(2), p(5)], 'T', T(:), 'xGrid', xGrid(:), 'xi', xi(:), ...
    'VFit', VFit(:), 'residual', residual(:), 'rmse', rmse, ...
    'costHist', costHist(:));
end

function T = transfer_template_from_anchor(cloud, Tanchor, xGrid, p, profileCfg)
xi = cloud.x - displacement_2freq(p, cloud.t);
TanchorAtXi = interp1(xGrid, Tanchor, xi, 'pchip', NaN);
valid = isfinite(TanchorAtXi) & isfinite(cloud.V);
if nnz(valid) < 0.75 * numel(cloud.V)
    T = estimate_template_from_cloud(xi, cloud.V, xGrid, profileCfg);
    return;
end
theta = [TanchorAtXi(valid), ones(nnz(valid), 1)] \ cloud.V(valid);
residual = cloud.V - (theta(1) * TanchorAtXi + theta(2));
residualGrid = estimate_template_from_cloud(xi(valid), residual(valid), xGrid, profileCfg);
residualGrid = smoothdata(residualGrid, 'sgolay', profileCfg.transferResidualSmoothWindow);
clipLimit = profileCfg.transferResidualClipFrac * max(range(Tanchor), eps);
residualGrid = min(max(residualGrid, -clipLimit), clipLimit);
T = theta(1) * Tanchor(:) + theta(2) + residualGrid(:);
T = fill_template_gaps(T, xGrid, isfinite(T));
T = smoothdata(T, 'sgolay', profileCfg.smoothWindow);
end

function pOpt = refine_p_given_template(cloud, T, xGrid, p0, lb, ub, profileCfg)
p0 = min(max(p0, lb), ub);
templateRange = max(T) - min(T);
templateEval = make_profile_template_eval(T, xGrid, profileCfg);
resFun = @(p) fixed_template_residual(p, cloud, templateRange, templateEval);
pOpt = run_lsqnonlin_or_fminsearch(resFun, p0, lb, ub, profileCfg.maxLsqIter);
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function [lb, ub] = make_profile_refine_bounds(p0, profileCfg)
lbWide = [0, -pi, min(profileCfg.f1Grid), 0, -pi, min(profileCfg.f2Grid)];
ubWide = [profileCfg.ampUpperBound, pi, max(profileCfg.f1Grid), ...
    profileCfg.ampUpperBound, pi, max(profileCfg.f2Grid)];
freqHalfWidth = profileCfg.localFreqHalfWidth;
ampScaleLo = profileCfg.localAmpScaleLo;
ampScaleHi = profileCfg.localAmpScaleHi;
ampFloor = 0.005;
lb = lbWide;
ub = ubWide;
lb(1) = max(lbWide(1), min(p0(1) * ampScaleLo, max(p0(1) - ampFloor, 0)));
ub(1) = min(ubWide(1), max(p0(1) * ampScaleHi, p0(1) + ampFloor));
lb(4) = max(lbWide(4), min(p0(4) * ampScaleLo, max(p0(4) - ampFloor, 0)));
ub(4) = min(ubWide(4), max(p0(4) * ampScaleHi, p0(4) + ampFloor));
lb(3) = max(lbWide(3), p0(3) - freqHalfWidth);
ub(3) = min(ubWide(3), p0(3) + freqHalfWidth);
lb(6) = max(lbWide(6), p0(6) - freqHalfWidth);
ub(6) = min(ubWide(6), p0(6) + freqHalfWidth);
end

function staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg)
gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gapQueryGrid = linspace(gapMin, gapMax, cfg.gapQueryN);
dxGrid = cfg.dxStaticGrid(:)';
y = highMap.V_a(:);
x = highMap.x_v(:);
J2D = zeros(numel(gapQueryGrid), numel(dxGrid));
xShiftGrid = x - dxGrid;
for ig = 1:numel(gapQueryGrid)
    g = gapQueryGrid(ig);
    [curveGrid, ~] = evaluate_gap_curve_grid(templateLib, g);
    yPredMat = interp1(templateLib.xGrid, curveGrid(:), xShiftGrid, 'pchip', 'extrap');
    resMat = y - yPredMat;
    J2D(ig, :) = mean(resMat.^2, 1, 'omitnan');
end
[J_static, dxIdx] = min(J2D, [], 2);
[~, bestIdx] = min(J_static);
staticState = struct('gHat', gapQueryGrid(bestIdx), ...
    'dx0', dxGrid(dxIdx(bestIdx)), 'gapQueryGrid', gapQueryGrid(:), ...
    'dxGrid', dxGrid(:), 'J2D', J2D, 'J_static', J_static(:));
end

function result = run_varpro_for_gap(highMap, templateLib, gFit, cfg, label)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
fixedEval = make_fixed_gap_evaluator(templateLib, gFit);
V0 = fixedEval.F(x);
K1 = fixedEval.dF(x);
DV = V - V0;
[pCandidates, coarseCosts] = coarse_varpro_2freq_candidates(t, DV, K1, ...
    cfg.f1Grid, cfg.f2Grid, cfg.numVarproCandidates);
lb = [0.001, -pi, min(cfg.f1Grid), 0.001, -pi, min(cfg.f2Grid)];
ub = [0.800,  pi, max(cfg.f1Grid), 0.800,  pi, max(cfg.f2Grid)];
resFun = @(p) residual_calc_6d_fixed_gap(p, t, V, x, fixedEval);
bestCost = inf;
pOpt = pCandidates(1, :);
for ic = 1:size(pCandidates, 1)
    pTry = run_lsqnonlin_or_fminsearch(resFun, pCandidates(ic, :), lb, ub, 300);
    costTry = sum(resFun(pTry).^2);
    if costTry < bestCost
        bestCost = costTry;
        pOpt = pTry;
    end
end
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
r = resFun(pOpt);
[f_id, order] = sort([pOpt(3), pOpt(6)]);
A_raw = [abs(pOpt(1)), abs(pOpt(4))];
phi_raw = [pOpt(2), pOpt(5)];
result = struct('label', label, 'g_used', gFit, 'p', pOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)), 'coarseCosts', coarseCosts);
end

function result = run_joint_refinement(highMap, templateLib, staticState, rawScanResult, cfg, jointCfg)
t = highMap.t_v(:);
V = highMap.V_a(:);
x = highMap.x_v(:);
p0 = rawScanResult.p(:).';
theta0 = [staticState.gHat, staticState.dx0, p0];
gLo = max(0.05, staticState.gHat - jointCfg.gHalfWidth);
gHi = staticState.gHat + jointCfg.gHalfWidth;
dxLo = staticState.dx0 - jointCfg.dxHalfWidth;
dxHi = staticState.dx0 + jointCfg.dxHalfWidth;
a1 = max(0.001, jointCfg.ampScaleLo * abs(p0(1)));
a2 = max(0.001, jointCfg.ampScaleLo * abs(p0(4)));
b1 = min(0.8, max(a1 + 1e-6, jointCfg.ampScaleHi * abs(p0(1))));
b2 = min(0.8, max(a2 + 1e-6, jointCfg.ampScaleHi * abs(p0(4))));
f1Lo = max(min(cfg.f1Grid), p0(3) - jointCfg.fHalfWidth);
f1Hi = min(max(cfg.f1Grid), p0(3) + jointCfg.fHalfWidth);
f2Lo = max(min(cfg.f2Grid), p0(6) - jointCfg.fHalfWidth);
f2Hi = min(max(cfg.f2Grid), p0(6) + jointCfg.fHalfWidth);
lb = [gLo, dxLo, a1, -pi, f1Lo, a2, -pi, f2Lo];
ub = [gHi, dxHi, b1,  pi, f1Hi, b2,  pi, f2Hi];
templateEval = make_template_evaluator(templateLib);
resFun = @(theta) residual_calc_joint(theta, t, V, x, templateEval);
thetaOpt = run_lsqnonlin_or_fminsearch(resFun, theta0, lb, ub, 400);
thetaOpt(4) = wrap_to_pi_local(thetaOpt(4));
thetaOpt(7) = wrap_to_pi_local(thetaOpt(7));
r = resFun(thetaOpt);
[f_id, order] = sort([thetaOpt(5), thetaOpt(8)]);
A_raw = [abs(thetaOpt(3)), abs(thetaOpt(6))];
phi_raw = [thetaOpt(4), thetaOpt(7)];
result = struct('label', 'gap_library_joint', ...
    'g_used', thetaOpt(1), 'dx_used', thetaOpt(2), ...
    'p', thetaOpt(3:end), 'theta', thetaOpt, ...
    'A_id', A_raw(order), 'f_id', f_id, 'phi_id', phi_raw(order), ...
    'rmse', sqrt(mean(r.^2)));
end

function [pCandidates, bestCosts] = coarse_varpro_2freq_candidates(t, DV, K1, f1Grid, f2Grid, numCandidates)
lambda = 1e-7 * max(numel(t), 1);
maxKeep = max(numCandidates * 4, numCandidates);
candCost = inf(maxKeep, 1);
candCoef = zeros(maxKeep, 4);
candFreq = zeros(maxKeep, 2);
f1Grid = f1Grid(:)';
f2Grid = f2Grid(:)';
S1 = sin(2*pi*t*f1Grid);
C1 = cos(2*pi*t*f1Grid);
S2 = sin(2*pi*t*f2Grid);
C2 = cos(2*pi*t*f2Grid);
KS1 = -K1 .* S1;
KC1 = -K1 .* C1;
KS2 = -K1 .* S2;
KC2 = -K1 .* C2;
for i1 = 1:numel(f1Grid)
    f1 = f1Grid(i1);
    h1s = KS1(:, i1);
    h1c = KC1(:, i1);
    for i2 = 1:numel(f2Grid)
        f2 = f2Grid(i2);
        if abs(f1 - f2) < 20, continue; end
        H = [h1s, h1c, KS2(:, i2), KC2(:, i2)];
        coef = (H' * H + lambda * eye(4)) \ (H' * DV);
        res = DV - H * coef;
        cost = sum(res.^2);
        [worstCost, worstIdx] = max(candCost);
        if cost < worstCost
            candCost(worstIdx) = cost;
            candCoef(worstIdx, :) = coef(:)';
            candFreq(worstIdx, :) = [f1, f2];
        end
    end
end
[candCost, order] = sort(candCost, 'ascend');
candCoef = candCoef(order, :);
candFreq = candFreq(order, :);
keep = isfinite(candCost);
candCost = candCost(keep);
candCoef = candCoef(keep, :);
candFreq = candFreq(keep, :);
uniqueRows = true(size(candCost));
for i = 2:numel(candCost)
    prev = candFreq(1:i-1, :);
    if any(abs(prev(:, 1) - candFreq(i, 1)) < 8 & abs(prev(:, 2) - candFreq(i, 2)) < 8)
        uniqueRows(i) = false;
    end
end
candCost = candCost(uniqueRows);
candCoef = candCoef(uniqueRows, :);
candFreq = candFreq(uniqueRows, :);
nOut = min(numCandidates, numel(candCost));
pCandidates = zeros(nOut, 6);
for i = 1:nOut
    c = candCoef(i, :);
    pCandidates(i, :) = [hypot(c(1), c(2)), atan2(c(2), c(1)), candFreq(i, 1), ...
        hypot(c(3), c(4)), atan2(c(4), c(3)), candFreq(i, 2)];
end
bestCosts = candCost(1:nOut);
end

function r = residual_calc_6d_fixed_gap(p, t, V, x, fixedEval)
u = displacement_2freq(p, t);
VSim = fixedEval.F(x - u);
r = V - VSim;
end

function r = residual_calc_joint(theta, t, V, x, templateEval)
gFit = theta(1);
dx0 = theta(2);
p = theta(3:end);
u = displacement_2freq(p, t);
VSim = templateEval.F(gFit, x - dx0 - u);
r = V - VSim;
end

function [curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g)
zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
SSort = templateLib.S(order, :);
FxSort = templateLib.FxMat(order, :);
zq = 1 / max(g, 1e-6);
curveGrid = interp1(zSort, SSort, zq, 'linear', 'extrap');
dGrid = interp1(zSort, FxSort, zq, 'linear', 'extrap');
end

function fixedEval = make_fixed_gap_evaluator(templateLib, g)
[curveGrid, dGrid] = evaluate_gap_curve_grid(templateLib, g);
xGrid = templateLib.xGrid(:);
F = griddedInterpolant(xGrid, curveGrid(:), 'pchip', 'linear');
dF = griddedInterpolant(xGrid, dGrid(:), 'pchip', 'linear');
fixedEval = struct();
fixedEval.F = @(xq) F(xq);
fixedEval.dF = @(xq) dF(xq);
end

function templateEval = make_template_evaluator(templateLib)
zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
templateEval = struct();
templateEval.xGrid = templateLib.xGrid(:);
templateEval.zSort = zSort(:);
templateEval.SSort = templateLib.S(order, :);
templateEval.F = @(g, xq) evaluate_template_from_cache(templateEval, g, xq);
end

function yPred = evaluate_template_from_cache(templateEval, g, xq)
zq = 1 / max(g, 1e-6);
curveGrid = interp1(templateEval.zSort, templateEval.SSort, zq, 'linear', 'extrap');
F = griddedInterpolant(templateEval.xGrid, curveGrid(:), 'pchip', 'linear');
yPred = F(xq);
end

function T = estimate_template_from_cloud(xi, V, xGrid, profileCfg)
xi = xi(:);
V = V(:);
xGrid = xGrid(:);
edges = [-inf; 0.5 * (xGrid(1:end-1) + xGrid(2:end)); inf];
bin = discretize(xi, edges);
validBin = isfinite(V) & isfinite(bin);
T = accumarray(bin(validBin), V(validBin), [numel(xGrid), 1], @mean, NaN);
valid = isfinite(T);
if nnz(valid) < 5
    error('Too few samples to estimate the shared template.');
end
T = fill_template_gaps(T, xGrid, valid);
T = smoothdata(T, 'sgolay', profileCfg.smoothWindow);
end

function [rmse, residual, VFit, xi] = calc_profile_rmse(cloud, T, xGrid, p)
xi = cloud.x - displacement_2freq(p, cloud.t);
VFit = interp1(xGrid, T, xi, 'pchip', NaN);
residual = cloud.V - VFit;
valid = isfinite(residual);
rmse = sqrt(mean(residual(valid).^2));
end

function cvRmse = calc_template_cv_rmse(cloud, templateAnchor, xGrid, p, profileCfg)
rev = cloud.rev(:);
splitA = mod(rev, 2) == 0;
splitB = ~splitA;
if nnz(splitA) < 0.25 * numel(rev) || nnz(splitB) < 0.25 * numel(rev)
    cvRmse = NaN;
    return;
end
cloudA = subset_cloud(cloud, splitA);
cloudB = subset_cloud(cloud, splitB);
TA = transfer_template_from_anchor(cloudA, templateAnchor, xGrid, p, profileCfg);
TB = transfer_template_from_anchor(cloudB, templateAnchor, xGrid, p, profileCfg);
rmseAB = calc_profile_rmse(cloudB, TA, xGrid, p);
rmseBA = calc_profile_rmse(cloudA, TB, xGrid, p);
cvRmse = mean([rmseAB, rmseBA], 'omitnan');
end

function cvRmse = calc_sensor_template_cv_rmse(cloud, templateAnchor, xGrid, p, profileCfg)
sensors = unique(cloud.sensor(:).');
scores = nan(numel(sensors), 1);
for i = 1:numel(sensors)
    testIdx = cloud.sensor(:) == sensors(i);
    trainIdx = ~testIdx;
    if nnz(testIdx) < 0.10 * numel(testIdx) || nnz(trainIdx) < 0.50 * numel(trainIdx)
        continue;
    end
    cloudTrain = subset_cloud(cloud, trainIdx);
    cloudTest = subset_cloud(cloud, testIdx);
    Ttrain = transfer_template_from_anchor(cloudTrain, templateAnchor, xGrid, p, profileCfg);
    scores(i) = calc_profile_rmse(cloudTest, Ttrain, xGrid, p);
end
cvRmse = mean(scores, 'omitnan');
end

function score = score_alias_candidate_sensorcv(cloud, templateAnchor, xGrid, p, profileCfg)
T = transfer_template_from_anchor(cloud, templateAnchor, xGrid, p, profileCfg);
trainRmse = calc_profile_rmse(cloud, T, xGrid, p);
sensorCvRmse = calc_sensor_template_cv_rmse(cloud, templateAnchor, xGrid, p, profileCfg);
roughness = calc_template_roughness(T, xGrid, templateAnchor);
score = blend_sensorcv_selection_rmse(trainRmse, sensorCvRmse, profileCfg) + ...
    profileCfg.roughnessWeight * roughness * max(range(templateAnchor), eps);
if ~isfinite(score)
    score = inf;
end
end

function cloudOut = subset_cloud(cloud, idx)
cloudOut = cloud;
cloudOut.t = cloud.t(idx);
cloudOut.x = cloud.x(idx);
cloudOut.V = cloud.V(idx);
cloudOut.rev = cloud.rev(idx);
cloudOut.sensor = cloud.sensor(idx);
end

function scoreRmse = blend_selection_rmse(trainRmse, cvRmse, profileCfg)
if ~isfinite(cvRmse)
    scoreRmse = trainRmse;
    return;
end
w = profileCfg.selectionCvWeight;
scoreRmse = (1 - w) * trainRmse + w * cvRmse;
end

function scoreRmse = blend_sensorcv_selection_rmse(trainRmse, sensorCvRmse, profileCfg)
if ~isfinite(sensorCvRmse)
    scoreRmse = trainRmse;
    return;
end
w = profileCfg.sensorCvWeight;
scoreRmse = (1 - w) * trainRmse + w * sensorCvRmse;
end

function roughness = calc_template_roughness(T, xGrid, templateAnchor)
T = T(:);
xGrid = xGrid(:);
dx = mean(diff(xGrid));
if numel(T) < 5 || ~isfinite(dx) || dx <= 0
    roughness = 0;
    return;
end
d1 = gradient(T, dx);
d2 = gradient(d1, dx);
scale = max(range(templateAnchor), eps);
roughness = sqrt(mean(d2(isfinite(d2)).^2, 'omitnan')) * dx^2 / scale;
if ~isfinite(roughness)
    roughness = 0;
end
end

function r = fixed_template_residual(p, cloud, templateRange, templateEval)
xi = cloud.x - displacement_2freq(p, cloud.t);
VFit = templateEval.F(xi);
w = informative_profile_weights(xi, templateEval);
r = sqrt(w) .* (cloud.V - VFit);
invalid = ~isfinite(r);
if any(invalid)
    r(invalid) = 5 * max(templateRange, eps);
end
end

function templateEval = make_profile_template_eval(T, xGrid, profileCfg)
xGrid = xGrid(:);
T = T(:);
dT = gradient(T, mean(diff(xGrid)));
templateEval = struct();
templateEval.F = griddedInterpolant(xGrid, T, 'pchip', 'none');
templateEval.dF = griddedInterpolant(xGrid, dT, 'pchip', 'nearest');
templateEval.weightFloor = profileCfg.fitWeightFloor;
templateEval.weightPower = profileCfg.fitWeightPower;
end

function w = informative_profile_weights(xi, templateEval)
dTAtXi = templateEval.dF(xi);
slopeScale = max(abs(dTAtXi));
if slopeScale <= eps
    w = ones(size(xi));
    return;
end
wShape = (abs(dTAtXi) ./ slopeScale) .^ templateEval.weightPower;
w = templateEval.weightFloor + (1 - templateEval.weightFloor) * wShape;
w(~isfinite(w)) = templateEval.weightFloor;
end

function mult = frequency_penalty_multiplier(p, profileCfg)
if isempty(p)
    mult = [];
    return;
end
p = reshape(p, [], 6);
f1 = p(:, 3);
f2 = p(:, 6);
hit1 = f1 <= min(profileCfg.f1Grid) + profileCfg.edgeGuardHz | ...
    f1 >= max(profileCfg.f1Grid) - profileCfg.edgeGuardHz;
hit2 = f2 <= min(profileCfg.f2Grid) + profileCfg.edgeGuardHz | ...
    f2 >= max(profileCfg.f2Grid) - profileCfg.edgeGuardHz;
mult = 1 + profileCfg.edgePenaltyWeight * (double(hit1) + double(hit2));
end

function T = fill_template_gaps(T, xGrid, valid)
idx = (1:numel(T))';
valid = valid(:) & isfinite(T(:));
if nnz(valid) < 2
    T(~valid) = 0;
    return;
end
first = find(valid, 1, 'first');
last = find(valid, 1, 'last');
left = idx < first;
right = idx > last;
midMissing = ~valid & idx >= first & idx <= last;
T(left) = T(first);
T(right) = T(last);
T(midMissing) = interp1(xGrid(valid), T(valid), xGrid(midMissing), 'pchip');
end

function pOpt = run_lsqnonlin_or_fminsearch(resFun, p0, lb, ub, maxIter)
p0 = min(max(p0, lb), ub);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, ...
        'MaxIterations', maxIter);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 5 * maxIter, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end

function p = local_percentile(x, pct)
x = sort(x(isfinite(x)));
if isempty(x), p = NaN; return; end
q = 1 + (numel(x) - 1) * pct / 100;
lo = floor(q);
hi = ceil(q);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (q - lo) * (x(hi) - x(lo));
end
end

function record = make_empty_record()
record = struct();
record.scenario = "";
record.num_revs_high = NaN;
record.method = "";
record.target_snr_db = NaN;
record.measured_snr_db = NaN;
record.noise_ratio = NaN;
record.repeat_id = NaN;
record.seed = NaN;
record.g_used_mm = NaN;
record.gap_error_mm = NaN;
record.f1_Hz = NaN;
record.f2_Hz = NaN;
record.A1_mm = NaN;
record.A2_mm = NaN;
record.mean_freq_error_Hz = NaN;
record.mean_amp_error_mm = NaN;
record.rmse = NaN;
record.elapsed_s = NaN;
record.is_success = NaN;
record.f1_boundary_hit = NaN;
record.f2_boundary_hit = NaN;
record.freq_boundary_hit = NaN;
record.a2_eta_dynamic = NaN;
record.a2_cond_dynamic = NaN;
end

function record = make_empty_snr_cal_record()
record = struct();
record.scenario = "";
record.num_revs_high = NaN;
record.target_snr_db = NaN;
record.noise_ratio = NaN;
record.signal_rms = NaN;
record.template_range = NaN;
end

function out = add_method_record(records, recId, cfg, truth, result, methodName, targetSnr, measuredSnr, noiseRatio, repeatId, seed, elapsed, gUsed, scenarioName, numRevsHigh)
recId = recId + 1;
[fTrue, trueOrder] = sort(truth.f(:).');
ATrue = truth.A(trueOrder);
[fEst, estOrder] = sort(result.f_id(:).');
AEst = result.A_id(estOrder);
records(recId).scenario = string(scenarioName);
records(recId).num_revs_high = numRevsHigh;
records(recId).method = string(methodName);
records(recId).target_snr_db = targetSnr;
records(recId).measured_snr_db = measuredSnr;
records(recId).noise_ratio = noiseRatio;
records(recId).repeat_id = repeatId;
records(recId).seed = seed;
records(recId).g_used_mm = gUsed;
if isfinite(gUsed)
    records(recId).gap_error_mm = abs(gUsed - cfg.g_holdout);
end
records(recId).f1_Hz = fEst(1);
records(recId).f2_Hz = fEst(2);
records(recId).A1_mm = AEst(1);
records(recId).A2_mm = AEst(2);
records(recId).mean_freq_error_Hz = mean(abs(fEst - fTrue));
records(recId).mean_amp_error_mm = mean(abs(AEst - ATrue));
records(recId).rmse = result.rmse;
records(recId).elapsed_s = elapsed;
records(recId).is_success = records(recId).mean_freq_error_Hz < cfg.failureFreqErrorHz;
records(recId).f1_boundary_hit = is_near_frequency_boundary(fEst(1), cfg.f1Grid, cfg.freqBoundaryTolHz);
records(recId).f2_boundary_hit = is_near_frequency_boundary(fEst(2), cfg.f2Grid, cfg.freqBoundaryTolHz);
records(recId).freq_boundary_hit = records(recId).f1_boundary_hit || records(recId).f2_boundary_hit;
if isfield(result, 'eta_dynamic')
    records(recId).a2_eta_dynamic = result.eta_dynamic;
end
if isfield(result, 'cond_dynamic')
    records(recId).a2_cond_dynamic = result.cond_dynamic;
end
out = struct('records', records, 'id', recId);
end

function tf = is_near_frequency_boundary(f, grid, tolHz)
tf = f <= min(grid) + tolHz || f >= max(grid) - tolHz;
end

function summaryTable = summarize_comparison(trialTable, failureFreqErrorHz)
groups = unique(trialTable(:, {'scenario', 'num_revs_high', 'method', 'target_snr_db'}), 'rows', 'stable');
summary = repmat(struct('scenario', "", 'num_revs_high', NaN, ...
    'method', "", 'target_snr_db', NaN, ...
    'freq_error_mean_Hz', NaN, 'freq_error_median_Hz', NaN, ...
    'freq_error_p90_Hz', NaN, 'freq_error_std_Hz', NaN, 'amp_error_mean_mm', NaN, ...
    'rmse_mean', NaN, 'elapsed_mean_s', NaN, ...
    'gap_error_mean_mm', NaN, 'success_rate', NaN, ...
    'boundary_hit_rate', NaN, 'num_trials', NaN), height(groups), 1);
for i = 1:height(groups)
    idx = trialTable.scenario == groups.scenario(i) & ...
        trialTable.num_revs_high == groups.num_revs_high(i) & ...
        trialTable.method == groups.method(i) & ...
        trialTable.target_snr_db == groups.target_snr_db(i);
    freqErr = trialTable.mean_freq_error_Hz(idx);
    ampErr = trialTable.mean_amp_error_mm(idx);
    gapErr = trialTable.gap_error_mm(idx);
    summary(i).scenario = groups.scenario(i);
    summary(i).num_revs_high = groups.num_revs_high(i);
    summary(i).method = groups.method(i);
    summary(i).target_snr_db = groups.target_snr_db(i);
    summary(i).freq_error_mean_Hz = mean(freqErr);
    summary(i).freq_error_median_Hz = median(freqErr);
    summary(i).freq_error_p90_Hz = local_percentile(freqErr, 90);
    summary(i).freq_error_std_Hz = std(freqErr);
    summary(i).amp_error_mean_mm = mean(ampErr);
    summary(i).rmse_mean = mean(trialTable.rmse(idx));
    summary(i).elapsed_mean_s = mean(trialTable.elapsed_s(idx));
    summary(i).gap_error_mean_mm = mean(gapErr, 'omitnan');
    summary(i).success_rate = mean(freqErr < failureFreqErrorHz);
    summary(i).boundary_hit_rate = mean(trialTable.freq_boundary_hit(idx));
    summary(i).num_trials = nnz(idx);
end
summaryTable = struct2table(summary);
end

function make_comparison_figure(summaryTable, outDir, resultPrefix)
figure('Name', 'Order/Revolution Method Comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 30, 18]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
methods = unique(summaryTable.method, 'stable');
scenarios = unique(summaryTable.scenario, 'stable');
colors = lines(numel(methods));
lineStyles = {'-', '--', ':', '-.'};

nexttile; hold on;
plot_order_rev_metric(summaryTable, methods, scenarios, colors, lineStyles, ...
    'freq_error_mean_Hz', 'Mean frequency error (Hz)', 'Mean frequency error');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'none', 'FontSize', 7);

nexttile; hold on;
plot_order_rev_metric(summaryTable, methods, scenarios, colors, lineStyles, ...
    'freq_error_median_Hz', 'Median frequency error (Hz)', 'Median frequency error');

nexttile; hold on;
plot_order_rev_metric(summaryTable, methods, scenarios, colors, lineStyles, ...
    'success_rate', 'Success rate', 'Success rate');
ylim([-0.05, 1.05]);

nexttile; hold on;
plot_order_rev_metric(summaryTable, methods, scenarios, colors, lineStyles, ...
    'elapsed_mean_s', 'Elapsed time (s)', 'Runtime');
exportgraphics(gcf, fullfile(outDir, [resultPrefix, '_summary.png']), 'Resolution', 300);
end

function plot_order_rev_metric(summaryTable, methods, scenarios, colors, lineStyles, metricName, yLabelText, titleText)
for i = 1:numel(methods)
    for j = 1:numel(scenarios)
        idx = summaryTable.method == methods(i) & summaryTable.scenario == scenarios(j);
        if ~any(idx)
            continue;
        end
        x = summaryTable.num_revs_high(idx);
        y = summaryTable.(metricName)(idx);
        [x, order] = sort(x);
        y = y(order);
        style = lineStyles{1 + mod(j - 1, numel(lineStyles))};
        label = sprintf('%s / %s', scenarios(j), methods(i));
        plot(x, y, ['o', style], 'LineWidth', 1.2, ...
            'Color', colors(i, :), 'DisplayName', label);
    end
end
grid on;
xlabel('Number of high-speed revolutions');
ylabel(yLabelText);
title(titleText);
xticks(sort(unique(summaryTable.num_revs_high)));
end
