%% Analyze_CH1_HighStatic_VibrationAware_20250527
% Diagnose whether the CH1 high-static residual is mainly caused by
% unmodeled synchronous vibration in the high-speed samples.
%
% The script compares two models on the same highMap samples:
%   1) static: corrected low-speed gap template with constant gap/dx
%   2) vibration-aware: static model plus A*sin(EO*theta + phi) displacement

clear; clc; close all;

routeDir = fileparts(mfilename('fullpath'));
C0 = CaseConfig();
outputDir = fullfile(routeDir, 'outputs', 'analysis_ch1_highstatic_vibration');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Parameters to tune
P.targetBlade = C0.bladeId;
P.analysisSensors = C0.sensorIds;
P.primarySensor = 1;
P.maxPointsPerSensor = 1800;       % use uniform samples from all highMap points
P.gGridCount = 11;
P.dxSearchMm = -0.35:0.02:0.35;
P.step06ILikeMaxPointsPerSensor = 320;
P.step06ILikeGGridCount = 61;
P.step06ILikeDxSearchMm = -0.08:0.01:0.08;
P.ampStartsMm = [0.03 0.12 0.25 0.38];
P.phaseStartsRad = [0 pi/2 pi 3*pi/2];
P.maxAbsDxMm = 0.60;
P.maxAbsAmpMm = 0.80;
P.minAffineGain = 0.4;
P.maxAffineGain = 1.8;
P.plotMaxPoints = 2500;

%% Input files
highMapFile = fullfile(routeDir, 'outputs', sprintf('Step06_HighMap_%s_%s.mat', C0.dataset, C0.caseTag));
if ~isfile(highMapFile)
    highMapFile = fullfile(routeDir, 'outputs', sprintf('Step06_HighMap_%s_B%d_%s.mat', C0.dataset, C0.bladeId, C0.sensorTag));
end
gapLibFile = fullfile(routeDir, 'outputs', sprintf('Step06I_OffsetTiltShared_GapLibrary_%s_%s.mat', C0.dataset, C0.caseTag));
if ~isfile(gapLibFile)
    gapLibFile = fullfile(routeDir, 'outputs', sprintf('Step06I_OffsetTiltShared_GapLibrary_%s_B%d_%s.mat', C0.dataset, C0.bladeId, C0.sensorTag));
end
step07File = fullfile(routeDir, 'outputs', sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_%s_main_gaptilt.mat', C0.dataset, C0.caseTag));
if ~isfile(step07File)
    step07File = fullfile(routeDir, 'outputs', sprintf('Step07J_NestedStaticWarp_VPFullWave_%s_B%d_%s_main_gaptilt.mat', C0.dataset, C0.bladeId, C0.sensorTag));
end

assert(exist(highMapFile, 'file') == 2, 'Missing highMap file: %s', highMapFile);
assert(exist(gapLibFile, 'file') == 2, 'Missing gap library file: %s', gapLibFile);

Sh = load(highMapFile, 'highMap');
Sg = load(gapLibFile, 'CorrectedGapLibrary', 'highCheckTable');
highMap = Sh.highMap;
CorrectedGapLibrary = Sg.CorrectedGapLibrary;
responseSurface = CorrectedGapLibrary.responseSurface;

St = load(CorrectedGapLibrary.templateFile, 'Template');
Template = St.Template;

eoCandidates = read_eo_candidates_local(step07File, highMap);
fprintf('EO candidates from final Step07J/highMap: %s\n', mat2str(eoCandidates));
fprintf('HighMap: %s\n', highMapFile);
fprintf('Gap library: %s\n', gapLibFile);
fprintf('Template: %s\n\n', CorrectedGapLibrary.templateFile);

%% Run diagnosis
resultRows = {};
fitCache = struct();

for is = 1:numel(P.analysisSensors)
    sid = P.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    corr = get_corrected_sensor_local(CorrectedGapLibrary, sid);
    D = select_highmap_points_local(highMap, sid, P.maxPointsPerSensor);
    D06 = select_highmap_points_local(highMap, sid, P.step06ILikeMaxPointsPerSensor);

    if isempty(D.x)
        warning('No usable highMap points for CH%d.', sid);
        continue;
    end

    gBounds = read_gap_bounds_local(CorrectedGapLibrary, responseSurface);
    gGrid = linspace(gBounds(1), gBounds(2), P.gGridCount);
    gGrid06 = linspace(gBounds(1), gBounds(2), P.step06ILikeGGridCount);
    static06Fit = fit_static_model_local(D06, Tsen, responseSurface, corr, ...
        gGrid06, P.step06ILikeDxSearchMm, "static_step06i_like", P);
    row = make_result_row_local(sid, "static_step06i_like", NaN, static06Fit, D06);
    resultRows{end+1, 1} = row; %#ok<SAGROW>

    staticFit = fit_static_model_local(D, Tsen, responseSurface, corr, ...
        gGrid, P.dxSearchMm, "static_wide", P);

    row = make_result_row_local(sid, "static_wide", NaN, staticFit, D);
    resultRows{end+1, 1} = row; %#ok<SAGROW>

    bestVib = [];
    for ieo = 1:numel(eoCandidates)
        eo = eoCandidates(ieo);
        vibFit = fit_vibration_model_local(D, Tsen, responseSurface, corr, gGrid, eo, staticFit, P);
        row = make_result_row_local(sid, "vibration_aware", eo, vibFit, D);
        resultRows{end+1, 1} = row; %#ok<SAGROW>

        if isempty(bestVib) || vibFit.rmseMv < bestVib.rmseMv
            bestVib = vibFit;
        end
    end

    fitCache(is).sensorId = sid; %#ok<SAGROW>
    fitCache(is).data = D;
    fitCache(is).step06ILikeFit = static06Fit;
    fitCache(is).staticFit = staticFit;
    fitCache(is).bestVibFit = bestVib;

    fprintf('CH%d | static %.2f mV | best vibration-aware EO%d %.2f mV | drop %.1f%% | A %.4f mm\n', ...
        sid, staticFit.rmseMv, bestVib.eo, bestVib.rmseMv, ...
        100 * (staticFit.rmseMv - bestVib.rmseMv) / max(staticFit.rmseMv, eps), ...
        bestVib.ampMm);
    fprintf('      Step06I-like static %.2f mV | dx %.3f mm | g %.4f mm | points %d\n', ...
        static06Fit.rmseMv, static06Fit.dxMm, static06Fit.gMm, static06Fit.validPointCount);
end

ResultTable = struct2table(vertcat(resultRows{:}));
ResultTable.rmseDropFromStaticPct = add_drop_column_local(ResultTable);

csvFile = fullfile(outputDir, 'CH1_HighStatic_VibrationAware_Diagnosis_20250527.csv');
writetable(ResultTable, csvFile);

matFile = fullfile(outputDir, 'CH1_HighStatic_VibrationAware_Diagnosis_20250527.mat');
Diagnostic = struct();
Diagnostic.highMapFile = highMapFile;
Diagnostic.gapLibFile = gapLibFile;
Diagnostic.step07File = step07File;
Diagnostic.templateFile = CorrectedGapLibrary.templateFile;
Diagnostic.eoCandidates = eoCandidates;
Diagnostic.params = P;
if isfield(Sg, 'highCheckTable')
    Diagnostic.originalStep06IHighCheckTable = Sg.highCheckTable;
end
save(matFile, 'Diagnostic', 'ResultTable', 'fitCache');

plotFile = fullfile(outputDir, 'CH1_HighStatic_Static_vs_VibrationAware_20250527.png');
plot_primary_sensor_fit_local(fitCache, P.primarySensor, plotFile, P);

fprintf('\nSaved CSV: %s\n', csvFile);
fprintf('Saved MAT: %s\n', matFile);
fprintf('Saved figure: %s\n', plotFile);

if isfield(Sg, 'highCheckTable')
    fprintf('\nOriginal Step06I high-static check:\n');
    disp(Sg.highCheckTable);
end

fprintf('\nDiagnosis result table:\n');
disp(ResultTable);

%% Local functions
function eoCandidates = read_eo_candidates_local(step07File, highMap)
eoCandidates = [];
if exist(step07File, 'file') == 2
    S = load(step07File, 'Summary');
    if isfield(S, 'Summary') && istable(S.Summary) && any(strcmp(S.Summary.Properties.VariableNames, 'dominant_EO'))
        eoCandidates = [eoCandidates; S.Summary.dominant_EO(:)];
    end
end
if isfield(highMap, 'referenceOrder') && isfinite(highMap.referenceOrder)
    eoCandidates = [eoCandidates; highMap.referenceOrder];
end
eoCandidates = unique(round(eoCandidates(isfinite(eoCandidates) & eoCandidates > 0))).';
if isempty(eoCandidates)
    eoCandidates = [11 14];
end
end

function D = select_highmap_points_local(highMap, sid, maxPoints)
if isfield(highMap, 'x_model_v') && numel(highMap.x_model_v) == numel(highMap.t_v)
    xModel = highMap.x_model_v(:);
elseif isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(highMap.t_v)
    xModel = highMap.x_fit_v(:);
else
    xModel = highMap.x_v(:);
end

if isfield(highMap, 'theta_v') && numel(highMap.theta_v) == numel(highMap.t_v)
    theta = highMap.theta_v(:);
else
    theta = 2 * pi * highMap.rotFreqHz .* (highMap.t_v(:) - min(highMap.t_v(:)));
end

mask = highMap.S_v(:) == sid & isfinite(xModel) & isfinite(highMap.V_a(:)) & ...
    isfinite(highMap.W_v(:)) & isfinite(theta);
idx = find(mask);
if numel(idx) > maxPoints
    idx = idx(round(linspace(1, numel(idx), maxPoints)));
end

D = struct();
D.sensorId = sid;
D.idx = idx;
D.x = xModel(idx);
D.v = highMap.V_a(idx);
D.w = max(highMap.W_v(idx), 0.05);
D.t = highMap.t_v(idx);
D.theta = theta(idx);
D.rev = highMap.rev_v(idx);
D.pointCount = numel(idx);
D.totalAvailable = nnz(mask);
end

function gBounds = read_gap_bounds_local(CorrectedGapLibrary, responseSurface)
if isfield(CorrectedGapLibrary, 'cfg') && isfield(CorrectedGapLibrary.cfg, 'gBoundsMm')
    gBounds = CorrectedGapLibrary.cfg.gBoundsMm;
else
    gBounds = [min(responseSurface.gTrainMm), max(responseSurface.gTrainMm)];
end
gBounds = sort(gBounds(:).');
end

function fit = fit_static_model_local(D, Tsen, responseSurface, corr, gGrid, dxSearchMm, modelType, P)
fit = empty_fit_local(modelType);
for g = gGrid(:).'
    for dx = dxSearchMm(:).'
        raw = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, D.x - dx);
        [rmse, pred, affine] = affine_rmse_local(D.v, raw, D.w, P);
        if rmse < fit.rmseMv
            fit = struct('modelType', modelType, 'eo', NaN, 'rmseMv', rmse, ...
                'gMm', g, 'dxMm', dx, 'ampMm', 0, 'phaseRad', NaN, ...
                'predMv', pred, 'affineGain', affine.gain, ...
                'affineOffsetMv', affine.offsetMv, 'validPointCount', nnz(isfinite(pred)));
        end
    end
end
end

function fit = fit_vibration_model_local(D, Tsen, responseSurface, corr, gGrid, eo, staticFit, P)
fit = empty_fit_local("vibration_aware");
options = optimset('Display', 'off', 'MaxIter', 320, 'MaxFunEvals', 1400, ...
    'TolX', 1e-5, 'TolFun', 1e-5);

for g = gGrid(:).'
    for ia = 1:numel(P.ampStartsMm)
        for ip = 1:numel(P.phaseStartsRad)
            p0 = [staticFit.dxMm, P.ampStartsMm(ia), P.phaseStartsRad(ip)];
            obj = @(p) vibration_objective_local(p, D, Tsen, responseSurface, corr, g, eo, P);
            pOpt = fminsearch(obj, p0, options);
            trial = evaluate_vibration_fit_local(pOpt, D, Tsen, responseSurface, corr, g, eo, P);
            if trial.rmseMv < fit.rmseMv
                fit = trial;
            end
        end
    end
end
end

function obj = vibration_objective_local(p, D, Tsen, responseSurface, corr, g, eo, P)
fit = evaluate_vibration_fit_local(p, D, Tsen, responseSurface, corr, g, eo, P);
penalty = 0;
if abs(p(1)) > P.maxAbsDxMm
    penalty = penalty + 1e5 * (abs(p(1)) - P.maxAbsDxMm)^2;
end
if abs(p(2)) > P.maxAbsAmpMm
    penalty = penalty + 1e5 * (abs(p(2)) - P.maxAbsAmpMm)^2;
end
obj = fit.rmseMv + penalty;
end

function fit = evaluate_vibration_fit_local(p, D, Tsen, responseSurface, corr, g, eo, P)
dx = p(1);
amp = p(2);
phase = wrapToPi_local(p(3));
vibDx = amp .* sin(eo .* D.theta + phase);
raw = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, D.x - dx - vibDx);
[rmse, pred, affine] = affine_rmse_local(D.v, raw, D.w, P);
fit = struct('modelType', "vibration_aware", 'eo', eo, 'rmseMv', rmse, ...
    'gMm', g, 'dxMm', dx, 'ampMm', abs(amp), 'phaseRad', phase, ...
    'predMv', pred, 'affineGain', affine.gain, ...
    'affineOffsetMv', affine.offsetMv, 'validPointCount', nnz(isfinite(pred)));
end

function model = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xOpr)
vLow = interp1(Tsen.x_grid(:), (Tsen.v_grid(:) - Tsen.baseline) * 1000, xOpr, 'pchip', NaN);
[Fdyn, ~] = eval_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, ~] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
end

function [model, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLib = k .* (xOpr - tau);
gEff = g0 + mu .* (xOpr - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
overshoot = max([max(xMin - xLib, [], 'omitnan'), max(xLib - xMax, [], 'omitnan'), ...
    max(gMin - gEff, [], 'omitnan'), max(gEff - gMax, [], 'omitnan'), 0]);
xLib = min(max(xLib, xMin), xMax);
gEff = min(max(gEff, gMin), gMax);
model = eval_response_surface_local(responseSurface, gEff, xLib);
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function [rmse, pred, affine] = affine_rmse_local(v, model, w, P)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model) & isfinite(w);
if nnz(valid) < 8
    rmse = inf;
    affine = struct('gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
beta = (H .* sw) \ (v(valid) .* sw);
if beta(1) < P.minAffineGain || beta(1) > P.maxAffineGain
    rmse = inf;
    affine = struct('gain', beta(1), 'offsetMv', beta(2));
    return;
end
pred(valid) = H * beta;
res = v(valid) - pred(valid);
rmse = sqrt(sum(w(valid) .* res.^2) / max(sum(w(valid)), eps));
affine = struct('gain', beta(1), 'offsetMv', beta(2));
end

function row = make_result_row_local(sid, modelType, eo, fit, D)
row = struct();
row.sensorId = sid;
row.modelType = modelType;
row.EO = eo;
row.weightedRmseMv = fit.rmseMv;
row.gMm = fit.gMm;
row.dxMm = fit.dxMm;
row.vibrationAmpMm = fit.ampMm;
row.phaseRad = fit.phaseRad;
row.affineGain = fit.affineGain;
row.affineOffsetMv = fit.affineOffsetMv;
row.validPointCount = fit.validPointCount;
row.usedPointCount = D.pointCount;
row.totalAvailablePointCount = D.totalAvailable;
end

function drops = add_drop_column_local(T)
drops = nan(height(T), 1);
for i = 1:height(T)
    staticMask = T.sensorId == T.sensorId(i) & T.modelType == "static_wide";
    if any(staticMask)
        staticRmse = T.weightedRmseMv(find(staticMask, 1));
        drops(i) = 100 * (staticRmse - T.weightedRmseMv(i)) / max(staticRmse, eps);
    end
end
end

function fit = empty_fit_local(modelType)
fit = struct('modelType', modelType, 'eo', NaN, 'rmseMv', inf, ...
    'gMm', NaN, 'dxMm', NaN, 'ampMm', NaN, 'phaseRad', NaN, ...
    'predMv', [], 'affineGain', NaN, 'affineOffsetMv', NaN, ...
    'validPointCount', 0);
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function corr = get_corrected_sensor_local(CorrectedGapLibrary, sid)
idx = find([CorrectedGapLibrary.sensor.sensorId] == sid, 1, 'first');
if isempty(idx)
    error('CorrectedGapLibrary does not contain CH%d.', sid);
end
corr = CorrectedGapLibrary.sensor(idx);
end

function angle = wrapToPi_local(angle)
angle = mod(angle + pi, 2 * pi) - pi;
end

function plot_primary_sensor_fit_local(fitCache, primarySensor, figFile, P)
idx = find([fitCache.sensorId] == primarySensor, 1, 'first');
if isempty(idx)
    return;
end
D = fitCache(idx).data;
Sfit = fitCache(idx).staticFit;
Vfit = fitCache(idx).bestVibFit;

plotIdx = 1:numel(D.x);
if numel(plotIdx) > P.plotMaxPoints
    plotIdx = round(linspace(1, numel(D.x), P.plotMaxPoints));
end

fig = figure('Name', 'CH1 high-static static vs vibration-aware diagnosis', ...
    'Color', 'w', 'Position', [80, 80, 1350, 780]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot(D.x(plotIdx), D.v(plotIdx), '.', 'Color', [0.20 0.20 0.20], 'MarkerSize', 6, ...
    'DisplayName', 'highMap data');
plot(D.x(plotIdx), Sfit.predMv(plotIdx), '.', 'Color', [0.85 0.20 0.15], 'MarkerSize', 6, ...
    'DisplayName', sprintf('static %.1f mV', Sfit.rmseMv));
plot(D.x(plotIdx), Vfit.predMv(plotIdx), '.', 'Color', [0.05 0.45 0.80], 'MarkerSize', 6, ...
    'DisplayName', sprintf('EO%d vib %.1f mV', Vfit.eo, Vfit.rmseMv));
xlabel('x_{OPR} (mm)'); ylabel('mV');
title(sprintf('CH%d high-static samples and fitted predictions', primarySensor));
legend('Location', 'best');

nexttile; hold on; grid on; box on;
plot(D.t(plotIdx), D.v(plotIdx) - Sfit.predMv(plotIdx), '.', 'Color', [0.85 0.20 0.15], ...
    'MarkerSize', 6, 'DisplayName', 'static residual');
plot(D.t(plotIdx), D.v(plotIdx) - Vfit.predMv(plotIdx), '.', 'Color', [0.05 0.45 0.80], ...
    'MarkerSize', 6, 'DisplayName', 'vibration-aware residual');
yline(0, 'k--', 'HandleVisibility', 'off');
xlabel('time (s)'); ylabel('residual (mV)');
legend('Location', 'best');

nexttile; hold on; grid on; box on;
histogram(D.v - Sfit.predMv, 55, 'FaceColor', [0.85 0.20 0.15], ...
    'FaceAlpha', 0.42, 'EdgeColor', 'none', 'DisplayName', 'static');
histogram(D.v - Vfit.predMv, 55, 'FaceColor', [0.05 0.45 0.80], ...
    'FaceAlpha', 0.42, 'EdgeColor', 'none', 'DisplayName', 'vibration-aware');
xlabel('residual (mV)'); ylabel('count');
title(sprintf('Best vibration-aware fit: EO%d, A=%.4f mm, dx=%.4f mm', ...
    Vfit.eo, Vfit.ampMm, Vfit.dxMm));
legend('Location', 'best');

exportgraphics(fig, figFile, 'Resolution', 220);
end
