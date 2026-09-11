%% Diagnostic_AdaptiveSG_Template_Comparison_20251222
% Isolated current-template versus grouped-lap-CV adaptive-SG comparison.
% Formal Step04, Main10, configuration, checks, and saved results are not modified.
clear; clc;

thisDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(thisDir));
workspaceRoot = fileparts(fileparts(packageRoot));
foundationRoot = fullfile(workspaceRoot, ...
    'experimental_dynamic_decoupling_validation', '20251222_btt_data_foundation');
addpath(fullfile(packageRoot, 'functions', 'gap_aware'));

resultDir = fullfile(packageRoot, 'diagnostics', ...
    'projection_trust_region_20260729', 'results');
defaultResult = fullfile(resultDir, ...
    'Main_GapAware_VPTopK_FullWave_20251222_B1_S123_diag_projection_allW_fullpoints.mat');
savedResult = strtrim(getenv('STEP07J_DIAG_RESULT_FILE'));
if isempty(savedResult)
    savedResult = defaultResult;
elseif ~isfile(savedResult)
    savedResult = fullfile(resultDir, savedResult);
end
assert(isfile(savedResult), 'Saved diagnostic result not found: %s', savedResult);

sourceFile = fullfile(foundationRoot, 'output', 'step04_low_speed_template', ...
    'LowSpeed_Template_SourceData_20251222.mat');
formalTemplateFile = fullfile(foundationRoot, 'output', ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat');
assert(isfile(sourceFile), 'Low-speed source point cloud not found: %s', sourceFile);
assert(isfile(formalTemplateFile), 'Formal low-speed template not found: %s', formalTemplateFile);

S = load(savedResult, 'Result');
windows = S.Result.WindowResult;
targetBlade = windows(1).bundle.Template.TargetBlade;
sensorIds = windows(1).bundle.sensorIds(:).';
maxWindows = min(numel(windows), parse_positive_integer_env_local( ...
    'STEP07J_DIAG_MAX_WINDOWS', numel(windows)));
windowIds = parse_window_ids_env_local('STEP07J_DIAG_WINDOW_IDS', ...
    1:maxWindows, numel(windows));
maxPoints = parse_positive_integer_env_local('STEP07J_DIAG_MAX_POINTS', 1200);

cfg = struct();
cfg.sgCandidateWindowMm = [0.10 0.18 0.34 0.50 0.82 1.22 1.62 2.02 2.42];
cfg.cvFoldCount = 5;
cfg.minBinCount = 5;
cfg.minFitPoints = 5;
cfg.vpGradientMinRatio = 0.10;
cfg.vpGradientReferenceQuantile = 95;
cfg.vpGradientWeightPower = 2;
cfg.weightFloor = 0.05;
cfg.amplitudeLimitMm = 0.50;
cfg.dxLimitMm = 0.35;
cfg.deltaGapLimitMm = 0.25;
cfg.deltaGapPriorScaleMm = 0.25;
cfg.staticRegWeightMv = 0.75;
cfg.continuousFrequencyRefine = true;
cfg.frequencyRefineHalfWidthHz = 2.0;
cfg.finalCandidateCount = 3;

P = load(sourceFile, 'point_cloud');
T = load(formalTemplateFile, 'Template');
[adaptiveSensor, screenSensor, templateDetail, cvTable, selection] = ...
    build_adaptive_template_local(P.point_cloud, T.Template, targetBlade, ...
    sensorIds, windows(1).bundle.Template.Sensor, cfg);

rows = repmat(empty_window_row_local(), 2 * numel(windowIds), 1);
details = repmat(struct('windowId', NaN, 'method', '', 'seedTable', table(), ...
    'topTable', table(), 'fits', {{}}), 2 * numel(windowIds), 1);
ir = 0;
for ik = 1:numel(windowIds)
    iw = windowIds(ik);
    baseBundle = decimate_bundle_local(windows(iw).bundle, maxPoints);
    currentModel = make_model_local(baseBundle);
    screenBundle = baseBundle;
    screenDerivativeModel = make_model_from_sensor_local(baseBundle, screenSensor);
    screenBundle = update_vp_split_linearization_local( ...
        screenBundle, currentModel, screenDerivativeModel);
    allEO = unique(round(windows(iw).VPSeedTable.EO(:).'));
    splitSeedTable = step07jcore.solve_gradient_displacement_vp_seed( ...
        screenBundle, allEO, cfg, []);
    splitSeedTable = splitSeedTable(isfinite(splitSeedTable.EO) & ...
        isfinite(splitSeedTable.A) & isfinite(splitSeedTable.phi) & ...
        isfinite(splitSeedTable.dx), :);
    methods = {'current_fixed_template', 'adaptive_grouped_cv_sg'};
    templates = {baseBundle.Template.Sensor, adaptiveSensor};
    for im = 1:2
        tStart = tic;
        bundle = baseBundle;
        bundle.Template.Sensor = templates{im};
        model = make_model_local(bundle);
        bundle = update_vp_linearization_local(bundle, model);
        eoCandidates = unique(round(windows(iw).VPSeedTable.EO(:).'));
        diagnosticSeedTable = step07jcore.solve_gradient_displacement_vp_seed( ...
            bundle, eoCandidates, cfg, []);
        diagnosticSeedTable = diagnosticSeedTable(isfinite(diagnosticSeedTable.EO) & ...
            isfinite(diagnosticSeedTable.A) & isfinite(diagnosticSeedTable.phi) & ...
            isfinite(diagnosticSeedTable.dx), :);
        frozenTable = windows(iw).VPSeedTable;
        frozenTable = frozenTable(isfinite(frozenTable.EO) & isfinite(frozenTable.A) & ...
            isfinite(frozenTable.phi) & isfinite(frozenTable.dx), :);
        assert(~isempty(frozenTable), 'No saved formal VP candidates for window %d.', iw);
        nTop = min(cfg.finalCandidateCount, height(frozenTable));
        topTable = frozenTable(1:nTop, :);
        fits = run_complete_fits_local(bundle, model, topTable, cfg);
        [best, margin] = select_best_fit_local(fits);
        querySafeFraction = final_query_safe_fraction_local(best, bundle, model);

        ir = ir + 1;
        row = empty_window_row_local();
        row.region = string(region_tag_local(savedResult));
        row.windowId = windows(iw).windowId;
        row.method = string(methods{im});
        row.targetBlade = targetBlade;
        row.selectedSgWindowMm = selection.selectedWindowMm;
        row.selectedSgSpanBins = selection.selectedSpanBins;
        row.screenSgWindowMm = selection.screenWindowMm;
        row.screenSgSpanBins = selection.screenSpanBins;
        row.candidateEO = join(string(topTable.EO(:).'), ',');
        row.diagnosticVpTopEO = join(string(diagnosticSeedTable.EO( ...
            1:min(cfg.finalCandidateCount, height(diagnosticSeedTable))).'), ',');
        row.finalEO = best.EO;
        row.finalCandidateRank = find(topTable.EO == best.EO, 1, 'first');
        row.finalDiagnosticVpRank = find(diagnosticSeedTable.EO == best.EO, 1, 'first');
        row.splitDerivativeVpTopEO = join(string(splitSeedTable.EO( ...
            1:min(cfg.finalCandidateCount, height(splitSeedTable))).'), ',');
        row.finalSplitDerivativeVpRank = find(splitSeedTable.EO == best.EO, 1, 'first');
        row.finalPlainRmseMv = best.plainRmseMv;
        row.candidateMarginMv = margin;
        row.continuousFrequencyHz = best.frequencyHz;
        row.frequencyDeltaFromEOHz = best.frequencyDeltaFromEOHz;
        row.amplitudeMm = best.amplitudeMm;
        row.phaseRad = best.phaseRad;
        row.dxMm = best.dxMm;
        row.deltaGapMm = join(compose('%.9g', best.deltaGapMm(:).'), ',');
        row.hitBoundary = logical(best.hitAmplitudeBoundary || ...
            best.hitDxBoundary || best.hitGapBoundary);
        row.querySafeFraction = querySafeFraction;
        row.pointCount = bundle.pointCount;
        row.runtimeSec = toc(tStart);
        rows(ir) = row;

        details(ir).windowId = row.windowId;
        details(ir).method = methods{im};
        details(ir).seedTable = diagnosticSeedTable;
        details(ir).topTable = topTable;
        details(ir).fits = fits;
        fprintf('%s W%02d %-24s EO%d %.4f mV | SG %.2f mm | %.1f s\n', ...
            row.region, row.windowId, methods{im}, row.finalEO, ...
            row.finalPlainRmseMv, row.selectedSgWindowMm, row.runtimeSec);
    end
end

windowTable = struct2table(rows(1:ir));
pairedTable = build_paired_table_local(windowTable);
regionTag = region_tag_local(savedResult);
outputTag = output_tag_local(windowIds, numel(windows));
stem = ['AdaptiveSG_Template_Comparison_20251222_', regionTag, outputTag];
writetable(windowTable, fullfile(thisDir, [stem, '_Window.csv']));
writetable(pairedTable, fullfile(thisDir, [stem, '_Paired.csv']));
writetable(templateDetail, fullfile(thisDir, [stem, '_Template.csv']));
writetable(cvTable, fullfile(thisDir, [stem, '_CV.csv']));
save(fullfile(thisDir, [stem, '.mat']), 'windowTable', 'pairedTable', ...
    'templateDetail', 'cvTable', 'selection', 'adaptiveSensor', 'screenSensor', 'details', ...
    'cfg', 'savedResult', 'sourceFile', 'formalTemplateFile', '-v7.3');
fprintf('Saved %s outputs in %s\n', stem, thisDir);


function [adaptiveSensor, screenSensor, detailTable, cvTable, selection] = ...
        build_adaptive_template_local(pointCloud, formalTemplate, bladeId, ...
        sensorIds, currentSensor, cfg)
nSensor = numel(sensorIds);
data = repmat(struct('x', [], 'v', [], 'lap', [], 'grid', [], ...
    'current', [], 'baseline', NaN, 'sensorId', NaN), nSensor, 1);
dxBySensor = nan(nSensor, 1);
for is = 1:nSensor
    sid = sensorIds(is);
    ip = find([pointCloud.sensor_id] == sid & ...
        [pointCloud.blade_id] == bladeId, 1);
    it = find([formalTemplate.SensorBlade.sensor_id] == sid & ...
        [formalTemplate.SensorBlade.blade_id] == bladeId, 1);
    ic = find([currentSensor.sensor_id] == sid, 1);
    assert(~isempty(ip) && ~isempty(it) && ~isempty(ic), ...
        'Missing low-speed sensor/blade entry S%d B%d.', sid, bladeId);
    pc = pointCloud(ip); ft = formalTemplate.SensorBlade(it); ct = currentSensor(ic);
    x = pc.x_mm(:) - ft.xc;
    v = pc.v(:);
    lap = round(pc.lap_index(:));
    grid = ct.x_grid(:);
    keep = isfinite(x) & isfinite(v) & isfinite(lap) & ...
        x >= min(grid) & x <= max(grid);
    if isfinite(pc.stable_window_start_lap) && isfinite(pc.stable_window_end_lap)
        keep = keep & lap >= pc.stable_window_start_lap & ...
            lap <= pc.stable_window_end_lap;
    end
    data(is).x = x(keep); data(is).v = v(keep); data(is).lap = lap(keep);
    data(is).grid = grid; data(is).current = ct.v_grid(:);
    data(is).baseline = ct.baseline; data(is).sensorId = sid;
    dxBySensor(is) = median(diff(grid));
end
dx = median(dxBySensor, 'omitnan');
spans = max(5, 2 .* floor((cfg.sgCandidateWindowMm ./ dx) ./ 2) + 1);
spans = unique(spans, 'stable');
widths = spans .* dx;

allLaps = unique(vertcat(data.lap));
foldId = mod((1:numel(allLaps)) - 1, cfg.cvFoldCount) + 1;
lossAdaptive = nan(cfg.cvFoldCount, nSensor, numel(spans));
lossCurrent = nan(cfg.cvFoldCount, nSensor);
foldAdaptive = cell(cfg.cvFoldCount, nSensor, numel(spans));
foldCurrent = cell(cfg.cvFoldCount, nSensor);
for f = 1:cfg.cvFoldCount
    testLaps = allLaps(foldId == f);
    for is = 1:nSensor
        train = ~ismember(data(is).lap, testLaps);
        test = ~train;
        if nnz(train) < 100 || nnz(test) < 20
            continue;
        end
        raw = bin_template_local(data(is).x(train), data(is).v(train), ...
            data(is).grid, cfg.minBinCount);
        currentFit = fixed_smooth_local(raw, 9);
        foldCurrent{f, is} = currentFit;
        lossCurrent(f, is) = prediction_rmse_local(data(is), test, currentFit);
        for ic = 1:numel(spans)
            adaptiveFit = sg_smooth_local(raw, spans(ic));
            foldAdaptive{f, is, ic} = adaptiveFit;
            lossAdaptive(f, is, ic) = ...
                prediction_rmse_local(data(is), test, adaptiveFit);
        end
    end
end
z = reshape(lossAdaptive, [], numel(spans));
cvMean = mean(z, 1, 'omitnan');
cvSe = std(z, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(z), 1));
[~, choice] = min(cvMean);
[~, iMin] = min(cvMean);
twoSeEligible = cvMean <= cvMean(iMin) + 2 .* cvSe(iMin);
screenChoice = find(twoSeEligible, 1, 'last');
selection = struct('selectedSpanBins', spans(choice), ...
    'selectedWindowMm', widths(choice), 'rule', ...
    'minimum grouped-lap voltage prediction RMSE', ...
    'foldCount', cfg.cvFoldCount, 'screenSpanBins', spans(screenChoice), ...
    'screenWindowMm', widths(screenChoice), 'screenRule', ...
    'largest window within two SE of minimum voltage prediction RMSE');
cvTable = table(spans(:), widths(:), cvMean(:), cvSe(:), twoSeEligible(:), ...
    'VariableNames', {'spanBins','windowMm','cvRmseMv','cvSeMv','twoSeEligible'});

adaptiveSensor = currentSensor;
screenSensor = currentSensor;
detail = repmat(struct('sensorId', NaN, 'bladeId', bladeId, ...
    'pointCount', NaN, 'lapCount', NaN, 'selectedSpanBins', spans(choice), ...
    'selectedWindowMm', widths(choice), 'currentCvRmseMv', NaN, ...
    'adaptiveCvRmseMv', NaN, 'cvRmseDeltaMv', NaN, ...
    'currentDerivativeFoldRmseMvPerMm', NaN, ...
    'adaptiveDerivativeFoldRmseMvPerMm', NaN, ...
    'currentDerivativeFoldCorrelation', NaN, ...
    'adaptiveDerivativeFoldCorrelation', NaN, ...
    'supportLeftMm', NaN, 'supportRightMm', NaN, ...
    'querySafeLeftMm', NaN, 'querySafeRightMm', NaN, ...
    'coverageFraction', NaN), nSensor, 1);
for is = 1:nSensor
    raw = bin_template_local(data(is).x, data(is).v, data(is).grid, cfg.minBinCount);
    adaptive = sg_smooth_local(raw, spans(choice));
    screen = sg_smooth_local(raw, spans(screenChoice));
    ia = find([adaptiveSensor.sensor_id] == data(is).sensorId, 1);
    adaptiveSensor(ia).v_grid = adaptive;
    adaptiveSensor(ia).dv_dx = gradient(adaptive, data(is).grid);
    screenSensor(ia).v_grid = screen;
    screenSensor(ia).dv_dx = gradient(screen, data(is).grid);
    [curDrmse, curCorr] = derivative_repeatability_local( ...
        foldCurrent(:, is), data(is).grid);
    adaptiveFolds = reshape(foldAdaptive(:, is, choice), [], 1);
    [adaDrmse, adaCorr] = derivative_repeatability_local( ...
        adaptiveFolds, data(is).grid);
    detail(is).sensorId = data(is).sensorId;
    detail(is).pointCount = numel(data(is).x);
    detail(is).lapCount = numel(unique(data(is).lap));
    detail(is).currentCvRmseMv = mean(lossCurrent(:, is), 'omitnan');
    detail(is).adaptiveCvRmseMv = mean(lossAdaptive(:, is, choice), 'omitnan');
    detail(is).cvRmseDeltaMv = detail(is).adaptiveCvRmseMv - detail(is).currentCvRmseMv;
    detail(is).currentDerivativeFoldRmseMvPerMm = curDrmse;
    detail(is).adaptiveDerivativeFoldRmseMvPerMm = adaDrmse;
    detail(is).currentDerivativeFoldCorrelation = curCorr;
    detail(is).adaptiveDerivativeFoldCorrelation = adaCorr;
    detail(is).supportLeftMm = min(data(is).grid);
    detail(is).supportRightMm = max(data(is).grid);
    ft = formalTemplate.SensorBlade([formalTemplate.SensorBlade.sensor_id] == ...
        data(is).sensorId & [formalTemplate.SensorBlade.blade_id] == bladeId);
    detail(is).querySafeLeftMm = ft.x_query_safe_domain(1);
    detail(is).querySafeRightMm = ft.x_query_safe_domain(2);
    detail(is).coverageFraction = mean(isfinite(raw));
end
detailTable = struct2table(detail);
end


function y = bin_template_local(x, v, grid, minCount)
dx = median(diff(grid));
edges = [grid(1) - dx / 2; grid(1:end-1) + diff(grid) / 2; grid(end) + dx / 2];
bin = discretize(x, edges);
ok = isfinite(bin) & isfinite(v);
y = accumarray(bin(ok), v(ok), [numel(grid), 1], @median, NaN);
count = accumarray(bin(ok), 1, [numel(grid), 1], @sum, 0);
y(count < minCount) = NaN;
valid = isfinite(y);
assert(nnz(valid) >= 5, 'Insufficient valid bins for low-speed template.');
y = fillmissing(y, 'linear', 'SamplePoints', grid);
y = fillmissing(y, 'nearest');
end


function y = fixed_smooth_local(y, span)
span = valid_odd_span_local(span, numel(y));
y = smoothdata(y(:), 'movmedian', span, 'omitnan');
y = smoothdata(y, 'movmean', span, 'omitnan');
end


function y = sg_smooth_local(y, span)
span = valid_odd_span_local(span, numel(y));
if span >= 5
    y = smoothdata(y(:), 'sgolay', span, 'Degree', min(3, span - 1));
else
    y = y(:);
end
end


function span = valid_odd_span_local(span, n)
span = min(round(span), n);
if mod(span, 2) == 0, span = span - 1; end
span = max(span, min(3, 2 * floor((n - 1) / 2) + 1));
end


function rmse = prediction_rmse_local(data, test, fit)
pred = interp1(data.grid, fit, data.x(test), 'pchip', NaN);
vTest = data.v(test);
ok = isfinite(pred) & isfinite(vTest);
if nnz(ok) < 20
    rmse = NaN;
else
    rmse = 1000 .* sqrt(mean((vTest(ok) - pred(ok)).^2));
end
end


function [rmse, corrMedian] = derivative_repeatability_local(folds, grid)
valid = ~cellfun(@isempty, folds);
folds = folds(valid);
if numel(folds) < 2
    rmse = NaN; corrMedian = NaN; return;
end
D = nan(numel(grid), numel(folds));
for i = 1:numel(folds)
    D(:, i) = gradient(folds{i}, grid);
end
ref = mean(D, 2, 'omitnan');
rmse = 1000 .* sqrt(mean((D - ref).^2, 'all', 'omitnan'));
C = corr(D, 'Rows', 'pairwise');
mask = triu(true(size(C)), 1);
c = C(mask & isfinite(C));
corrMedian = median(c, 'omitnan');
end


function bundle = update_vp_linearization_local(bundle, model)
bundle.F0 = nan(size(bundle.V));
bundle.Fx = nan(size(bundle.V));
for is = 1:numel(bundle.sensorIds)
    mask = bundle.sensorIndex == is;
    [v, aux] = model(is).evaluate(0, bundle.X(mask));
    bundle.F0(mask) = v;
    bundle.Fx(mask) = aux.dVoltageDx;
end
end


function bundle = update_vp_split_linearization_local(bundle, forwardModel, derivativeModel)
bundle.F0 = nan(size(bundle.V));
bundle.Fx = nan(size(bundle.V));
for is = 1:numel(bundle.sensorIds)
    mask = bundle.sensorIndex == is;
    [v, ~] = forwardModel(is).evaluate(0, bundle.X(mask));
    [~, derivativeAux] = derivativeModel(is).evaluate(0, bundle.X(mask));
    bundle.F0(mask) = v;
    bundle.Fx(mask) = derivativeAux.dVoltageDx;
end
end


function fits = run_complete_fits_local(bundle, model, candidateTable, cfg)
fits = cell(height(candidateTable), 1);
for i = 1:height(candidateTable)
    candidate = struct('EO', candidateTable.EO(i), 'A', candidateTable.A(i), ...
        'phi', candidateTable.phi(i), 'dx', candidateTable.dx(i));
    fits{i} = step07jcore.optimize_experimental_full_waveform( ...
        bundle, candidate, model, cfg);
end
end


function [best, margin] = select_best_fit_local(fits)
rmse = cellfun(@(x) x.plainRmseMv, fits);
[rmseSorted, order] = sort(rmse, 'ascend');
best = fits{order(1)};
if numel(order) > 1, margin = rmseSorted(2) - rmseSorted(1); else, margin = NaN; end
end


function fraction = final_query_safe_fraction_local(fit, bundle, model)
u = fit.amplitudeMm .* sin((fit.frequencyHz / bundle.rotFreqMeanHz) .* ...
    bundle.Theta(:) + fit.phaseRad);
valid = false(size(bundle.V(:)));
for is = 1:numel(bundle.sensorIds)
    mask = bundle.sensorIndex(:) == is;
    [~, aux] = model(is).evaluate(fit.deltaGapMm(is), ...
        bundle.X(mask) - fit.dxMm - u(mask));
    valid(mask) = aux.querySafe;
end
fraction = mean(valid);
end


function model = make_model_local(bundle)
nSensor = numel(bundle.sensorIds);
model = repmat(struct('xGrid', [], 'evaluate', []), nSensor, 1);
for is = 1:nSensor
    sid = bundle.sensorIds(is);
    it = find([bundle.Template.Sensor.sensor_id] == sid, 1);
    ic = find([bundle.CorrectedGapLibrary.sensor.sensorId] == sid, 1);
    Tpl = bundle.Template.Sensor(it);
    corr = bundle.CorrectedGapLibrary.sensor(ic);
    rs = bundle.responseSurface;
    model(is).xGrid = Tpl.x_grid(:);
    model(is).evaluate = @(dg, x) local_eval_sensor(Tpl, rs, corr, dg, x);
end
end


function model = make_model_from_sensor_local(bundle, sensorTemplate)
copy = bundle;
copy.Template.Sensor = sensorTemplate;
model = make_model_local(copy);
end


function [v, aux] = local_eval_sensor(Tpl, rs, corr, dg, x)
[v, aux0] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x, 0, 0);
h = 1e-4;
[vp, ~] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x + h, 0, 0);
[vm, ~] = step07jcore.evaluate_experimental_sensor_voltage( ...
    Tpl, rs, corr, dg, x - h, 0, 0);
aux = aux0;
aux.dVoltageDx = (vp - vm) ./ (2 .* h);
end


function b = decimate_bundle_local(b, maxPoints)
if numel(b.X) <= maxPoints, return; end
keep = false(numel(b.X), 1);
for is = 1:numel(b.sensorIds)
    ids = find(b.sensorIndex == is);
    nk = min(numel(ids), max(1, round(maxPoints * numel(ids) / numel(b.X))));
    keep(ids(unique(round(linspace(1, numel(ids), nk))))) = true;
end
fields = {'X','T','TRel','V','W','Theta','S','sensorIndex','F0','Fx'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(b, f) && numel(b.(f)) == numel(keep), b.(f) = b.(f)(keep); end
end
b.pointCount = nnz(keep);
end


function paired = build_paired_table_local(T)
current = T(T.method == "current_fixed_template", :);
adaptive = T(T.method == "adaptive_grouped_cv_sg", :);
paired = table(current.region, current.windowId, current.finalEO, adaptive.finalEO, ...
    current.finalPlainRmseMv, adaptive.finalPlainRmseMv, ...
    adaptive.finalPlainRmseMv - current.finalPlainRmseMv, ...
    current.amplitudeMm, adaptive.amplitudeMm, adaptive.amplitudeMm - current.amplitudeMm, ...
    current.dxMm, adaptive.dxMm, adaptive.dxMm - current.dxMm, ...
    current.querySafeFraction, adaptive.querySafeFraction, ...
    current.hitBoundary, adaptive.hitBoundary, ...
    'VariableNames', {'region','windowId','currentEO','adaptiveEO', ...
    'currentPlainRmseMv','adaptivePlainRmseMv','rmseDeltaMv', ...
    'currentAmplitudeMm','adaptiveAmplitudeMm','amplitudeDeltaMm', ...
    'currentDxMm','adaptiveDxMm','dxDeltaMm', ...
    'currentQuerySafeFraction','adaptiveQuerySafeFraction', ...
    'currentHitBoundary','adaptiveHitBoundary'});
end


function row = empty_window_row_local()
row = struct('region', "", 'windowId', NaN, 'method', "", ...
    'targetBlade', NaN, 'selectedSgWindowMm', NaN, ...
    'selectedSgSpanBins', NaN, 'screenSgWindowMm', NaN, ...
    'screenSgSpanBins', NaN, 'candidateEO', "", 'diagnosticVpTopEO', "", ...
    'finalEO', NaN, 'finalCandidateRank', NaN, ...
    'finalDiagnosticVpRank', NaN, 'splitDerivativeVpTopEO', "", ...
    'finalSplitDerivativeVpRank', NaN, 'finalPlainRmseMv', NaN, ...
    'candidateMarginMv', NaN, 'continuousFrequencyHz', NaN, ...
    'frequencyDeltaFromEOHz', NaN, 'amplitudeMm', NaN, 'phaseRad', NaN, ...
    'dxMm', NaN, 'deltaGapMm', "", 'hitBoundary', false, ...
    'querySafeFraction', NaN, 'pointCount', NaN, 'runtimeSec', NaN);
end


function value = parse_positive_integer_env_local(name, defaultValue)
raw = strtrim(getenv(name));
if isempty(raw), value = defaultValue; return; end
value = str2double(raw);
assert(isfinite(value) && value >= 1, '%s must be a positive integer.', name);
value = floor(value);
end


function ids = parse_window_ids_env_local(name, defaultIds, maxId)
raw = strtrim(getenv(name));
if isempty(raw), ids = defaultIds; return; end
ids = unique(round(sscanf(raw, '%f').'), 'stable');
assert(~isempty(ids) && all(ids >= 1) && all(ids <= maxId), ...
    '%s must contain window IDs from 1 to %d.', name, maxId);
end


function tag = output_tag_local(ids, totalCount)
if isequal(ids(:).', 1:totalCount)
    tag = '';
else
    tag = ['_W', char(join(compose('%02d', ids), '_'))];
end
end


function tag = region_tag_local(fileName)
name = lower(string(fileName));
if contains(name, 'r07'), tag = 'R07';
elseif contains(name, 'r04'), tag = 'R04';
else, tag = 'R01';
end
end
