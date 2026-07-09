%% Step06F: Calibrate tilt-corrected per-sensor gap library
% The raw gap library F_raw(g,x) is generated in its own geometry. In the
% rotating experiment, blade/sensor inclination means the blade passage is a
% slanted trajectory in the raw library's (x,g) plane:
%
%   x_lib = k_s * (x_OPR - tau_s)
%   g_eff = g0_s + mu_s * (x_OPR - tau_s)
%
% This step fits each sensor's geometry against the measured low-speed
% template:
%
%   T_low,s(x) ~= a_s * F_raw(g0_s + mu_s*(x-tau_s), k_s*(x-tau_s)) + b_s
%
% Then it builds a corrected incremental library:
%
%   F_corr,s(g,x) = T_low,s(x)
%     + a_s * [F_raw(g + mu_s*(x-tau_s), k_s*(x-tau_s))
%              - F_raw(g0_s + mu_s*(x-tau_s), k_s*(x-tau_s))]

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06f_tilt_corrected_gap_library');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 2, 3];
cfg.gBoundsMm = [];
cfg.tauBoundsMm = [-0.80, 0.80];
cfg.xScaleBounds = [0.65, 1.45];
cfg.muBounds = [-0.35, 0.35];
cfg.minFitPoints = 100;
cfg.highCheckMaxPointsPerSensor = 320;
cfg.highCheckGapGridCount = 61;
cfg.highCheckDxCorrectionGridMm = -0.08:0.01:0.08;
cfg.overshootPenaltyMvPerMm = 600;
cfg.fminOptions = optimset('Display', 'off', 'MaxIter', 500, 'MaxFunEvals', 1800, ...
    'TolX', 1e-6, 'TolFun', 1e-6);

sensorOverride = strtrim(getenv('STEP06F_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06F_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

responseFile = fullfile(outDir, 'Step05_Response_Surface_20251222.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20251222_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
templateDir = fullfile(rootDir, '20251222_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(templateFile)
    error('Missing low-speed template file: %s', templateFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
if isempty(cfg.gBoundsMm)
    cfg.gBoundsMm = [min(responseSurface.gTrainMm), max(responseSurface.gTrainMm)];
end
St = load(templateFile, 'Template');
Template = St.Template;
hasHighMap = isfile(highMapFile);
if hasHighMap
    H = load(highMapFile, 'highMap');
    highMap = attach_model_coordinate_local(H.highMap);
else
    highMap = struct();
end

fprintf('\n=== Step06F: tilt-corrected gap-library calibration ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Sensors: %s\n', mat2str(cfg.analysisSensors));

sensorCorr = repmat(struct('sensorId', NaN, 'g0Mm', NaN, 'tauMm', NaN, ...
    'xScale', NaN, 'muGapPerXMm', NaN, 'tiltAngleDeg', NaN, ...
    'voltageGain', NaN, 'voltageOffsetMv', NaN, 'lowFitRmseMv', NaN, ...
    'overshootRmseMm', NaN, 'pointCount', NaN, 'x', [], 'vLowMv', [], ...
    'vFitMv', [], 'xLib', [], 'gEff', [], 'residualMv', []), numel(cfg.analysisSensors), 1);
rows = cell(numel(cfg.analysisSensors), 1);
highRows = {};
highCache = struct([]);

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    fit = fit_tilt_corrected_low_template_local(Tsen, responseSurface, cfg);
    sensorCorr(is) = fit;
    rows{is} = table(sid, fit.g0Mm, fit.tauMm, fit.xScale, fit.muGapPerXMm, ...
        fit.tiltAngleDeg, fit.voltageGain, fit.voltageOffsetMv, fit.lowFitRmseMv, ...
        fit.overshootRmseMm, fit.pointCount, ...
        'VariableNames', {'sensorId','g0Mm','tauMm','xScale','muGapPerXMm', ...
        'tiltAngleDeg','voltageGain','voltageOffsetMv','lowFitRmseMv', ...
        'overshootRmseMm','pointCount'});
    fprintf('  CH%d: g0 %.4f, tau %.4f, k %.4f, mu %.4f (%.2f deg), low RMSE %.2f mV\n', ...
        sid, fit.g0Mm, fit.tauMm, fit.xScale, fit.muGapPerXMm, ...
        fit.tiltAngleDeg, fit.lowFitRmseMv);

    if hasHighMap
        [hfit, hcache] = high_static_check_tilt_corrected_local( ...
            highMap, sid, Tsen, responseSurface, fit, cfg);
        highRows{end+1, 1} = table(sid, hfit.gMm, hfit.dxCorrectionMm, ...
            hfit.weightedRmseMv, hfit.pointCount, ...
            'VariableNames', {'sensorId','highStaticGMm','dxCorrectionMm', ...
            'weightedRmseMv','pointCount'}); %#ok<SAGROW>
        highCache(end+1).sensorId = sid; %#ok<SAGROW>
        highCache(end).fit = hfit;
        highCache(end).cache = hcache;
    end
end

calibrationTable = vertcat(rows{:});
if isempty(highRows)
    highCheckTable = table();
else
    highCheckTable = vertcat(highRows{:});
end

CorrectedGapLibrary = struct();
CorrectedGapLibrary.dataset = '20251222';
CorrectedGapLibrary.method = 'low_speed_template_tilt_corrected_incremental_gap_library';
CorrectedGapLibrary.description = ['Per-sensor gap library calibrated with low-speed template using a slanted ' ...
    '(x,g) trajectory. The low-speed template is the baseline; the raw gap library supplies clearance-change increments.'];
CorrectedGapLibrary.responseFile = responseFile;
CorrectedGapLibrary.templateFile = templateFile;
CorrectedGapLibrary.targetBlade = cfg.targetBlade;
CorrectedGapLibrary.analysisSensors = cfg.analysisSensors(:).';
CorrectedGapLibrary.responseSurface = responseSurface;
CorrectedGapLibrary.sensor = sensorCorr;
CorrectedGapLibrary.formula = 'T_low(x)+a*(F_raw(g+mu*(x-tau),k*(x-tau))-F_raw(g0+mu*(x-tau),k*(x-tau)))';
CorrectedGapLibrary.cfg = cfg;

matFile = fullfile(outDir, sprintf('Step06F_TiltCorrected_GapLibrary_20251222_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step06F_TiltCorrected_GapLibrary_Calibration_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
highCsvFile = fullfile(outDir, sprintf('Step06F_TiltCorrected_GapLibrary_HighStaticCheck_20251222_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figLow = fullfile(figDir, sprintf('Step06F_LowTemplate_TiltCorrected_Fit_20251222_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figPath = fullfile(figDir, sprintf('Step06F_TiltPath_In_RawGapLibrary_20251222_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figHigh = fullfile(figDir, sprintf('Step06F_HighStatic_TiltCorrected_Check_20251222_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'CorrectedGapLibrary', 'calibrationTable', 'highCheckTable', ...
    'responseFile', 'templateFile', 'highMapFile', '-v7.3');
writetable(calibrationTable, csvFile);
if ~isempty(highCheckTable)
    writetable(highCheckTable, highCsvFile);
end
plot_low_fit_local(sensorCorr, figLow);
plot_tilt_path_local(sensorCorr, responseSurface, figPath);
if hasHighMap && ~isempty(highCache)
    plot_high_static_check_local(highCache, figHigh);
end

fprintf('\nStep06F complete.\n');
disp(calibrationTable);
if ~isempty(highCheckTable)
    disp(highCheckTable);
end
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, csvFile, figLow);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20251222.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20251222.mat', targetBlade, exactTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20251222.mat', targetBlade)));
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)(?:_GradientXRange030)?_20251222\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    available = arrayfun(@(c) str2double(c), token{1});
    if all(ismember(sensorIds, available))
        templateFile = fullfile(files(i).folder, files(i).name);
        return;
    end
end
error('No low-speed template file found under %s for sensors %s.', templateDir, mat2str(sensorIds));
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function highMap = attach_model_coordinate_local(highMap)
if isfield(highMap, 'x_fit_v') && numel(highMap.x_fit_v) == numel(highMap.t_v) && any(isfinite(highMap.x_fit_v))
    highMap.x_model_v = highMap.x_fit_v(:);
else
    highMap.x_model_v = highMap.x_v(:);
end
end

function fit = fit_tilt_corrected_low_template_local(Tsen, responseSurface, cfg)
x = Tsen.x_grid(:);
vLow = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
mask = isfinite(x) & isfinite(vLow);
x = x(mask);
vLow = vLow(mask);
if numel(x) < cfg.minFitPoints
    error('CH%d has too few low-template points.', Tsen.sensor_id);
end

gStarts = unique([min(responseSurface.gTrainMm), median(responseSurface.gTrainMm), max(responseSurface.gTrainMm), 1.04]);
tauStarts = [-0.25, 0, 0.25];
kStarts = [0.85, 1.0, 1.15];
muStarts = [-0.12, 0, 0.12];
best = struct('objective', inf);
for g0 = gStarts(:).'
    for tau = tauStarts
        for k = kStarts
            for mu = muStarts
                p0 = [g0, tau, k, mu];
                fun = @(p) tilt_low_template_objective_scalar_local( ...
                    p, x, vLow, responseSurface, cfg);
                pOpt = fminsearch(fun, p0, cfg.fminOptions);
                [obj, detail] = evaluate_tilt_low_template_local(pOpt, x, vLow, responseSurface, cfg);
                if obj < best.objective
                    best = detail;
                    best.objective = obj;
                end
            end
        end
    end
end
fit = struct('sensorId', Tsen.sensor_id, 'g0Mm', best.g0Mm, 'tauMm', best.tauMm, ...
    'xScale', best.xScale, 'muGapPerXMm', best.muGapPerXMm, ...
    'tiltAngleDeg', atan(best.muGapPerXMm) * 180 / pi, ...
    'voltageGain', best.voltageGain, 'voltageOffsetMv', best.voltageOffsetMv, ...
    'lowFitRmseMv', best.lowFitRmseMv, 'overshootRmseMm', best.overshootRmseMm, ...
    'pointCount', best.pointCount, 'x', x, 'vLowMv', vLow, 'vFitMv', best.vFitMv, ...
    'xLib', best.xLib, 'gEff', best.gEff, 'residualMv', vLow - best.vFitMv);
end

function obj = tilt_low_template_objective_scalar_local(p, x, vLow, responseSurface, cfg)
[obj, ~] = evaluate_tilt_low_template_local(p, x, vLow, responseSurface, cfg);
end

function [obj, detail] = evaluate_tilt_low_template_local(pRaw, x, vLow, responseSurface, cfg)
p = pRaw(:).';
[g0, tau, k, mu, penalty] = clamp_tilt_params_local(p, cfg);
[modelRaw, xLib, gEff, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, x);
[rmse, vFit, affine] = affine_rmse_local(vLow, modelRaw, ones(size(vLow)));
overshootRmse = sqrt(mean(overshoot.^2, 'omitnan'));
obj = rmse + cfg.overshootPenaltyMvPerMm * overshootRmse + penalty;
detail = struct('g0Mm', g0, 'tauMm', tau, 'xScale', k, 'muGapPerXMm', mu, ...
    'voltageGain', affine.gain, 'voltageOffsetMv', affine.offsetMv, ...
    'lowFitRmseMv', rmse, 'overshootRmseMm', overshootRmse, ...
    'pointCount', nnz(isfinite(vFit)), 'vFitMv', vFit, ...
    'xLib', xLib, 'gEff', gEff, 'overshoot', overshoot);
end

function [g0, tau, k, mu, penalty] = clamp_tilt_params_local(p, cfg)
gRaw = p(1); tauRaw = p(2); kRaw = p(3); muRaw = p(4);
g0 = min(max(gRaw, cfg.gBoundsMm(1)), cfg.gBoundsMm(2));
tau = min(max(tauRaw, cfg.tauBoundsMm(1)), cfg.tauBoundsMm(2));
k = min(max(kRaw, cfg.xScaleBounds(1)), cfg.xScaleBounds(2));
mu = min(max(muRaw, cfg.muBounds(1)), cfg.muBounds(2));
penalty = 1e5 * ((gRaw - g0)^2 + (tauRaw - tau)^2 + ...
    (kRaw - k)^2 + (muRaw - mu)^2);
end

function [model, xLib, gEff, overshoot] = eval_raw_tilt_path_local(responseSurface, g0, tau, k, mu, xOpr)
xLibRaw = k .* (xOpr - tau);
gEffRaw = g0 + mu .* (xOpr - tau);
xMin = min(responseSurface.xGrid(:));
xMax = max(responseSurface.xGrid(:));
gMin = min(responseSurface.gTrainMm(:));
gMax = max(responseSurface.gTrainMm(:));
xLib = min(max(xLibRaw, xMin), xMax);
gEff = min(max(gEffRaw, gMin), gMax);
overshootX = max(xMin - xLibRaw, 0) + max(xLibRaw - xMax, 0);
overshootG = max(gMin - gEffRaw, 0) + max(gEffRaw - gMax, 0);
overshoot = hypot(overshootX, overshootG);
model = eval_response_surface_local(responseSurface, gEff, xLib);
end

function [fit, cache] = high_static_check_tilt_corrected_local(highMap, sid, Tsen, responseSurface, corr, cfg)
mask = highMap.S_v == sid & isfinite(highMap.x_model_v) & isfinite(highMap.V_a) & isfinite(highMap.W_v);
idx = find(mask);
if numel(idx) > cfg.highCheckMaxPointsPerSensor
    idx = idx(round(linspace(1, numel(idx), cfg.highCheckMaxPointsPerSensor)));
end
x = highMap.x_model_v(idx);
v = highMap.V_a(idx);
w = max(highMap.W_v(idx), 0.05);
gGrid = linspace(cfg.gBoundsMm(1), cfg.gBoundsMm(2), cfg.highCheckGapGridCount);
fit = struct('gMm', NaN, 'dxCorrectionMm', NaN, 'weightedRmseMv', inf, ...
    'pointCount', 0, 'pred', []);
for g = gGrid(:).'
    for dx = cfg.highCheckDxCorrectionGridMm(:).'
        predRaw = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, x - dx);
        [rmse, pred] = affine_rmse_local(v, predRaw, w);
        if rmse < fit.weightedRmseMv
            fit = struct('gMm', g, 'dxCorrectionMm', dx, ...
                'weightedRmseMv', rmse, 'pointCount', nnz(isfinite(pred)), 'pred', pred);
        end
    end
end
cache = struct('x', x, 'v', v, 'w', w, 'pred', fit.pred);
end

function model = eval_corrected_gap_template_local(Tsen, responseSurface, corr, g, xOpr)
vLow = interp1(Tsen.x_grid(:), (Tsen.v_grid(:) - Tsen.baseline) * 1000, xOpr, 'pchip', NaN);
[Fdyn, ~, ~, ~] = eval_raw_tilt_path_local(responseSurface, g, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
[Fbase, ~, ~, ~] = eval_raw_tilt_path_local(responseSurface, corr.g0Mm, corr.tauMm, ...
    corr.xScale, corr.muGapPerXMm, xOpr);
model = vLow + corr.voltageGain .* (Fdyn - Fbase);
end

function [rmse, pred, affine] = affine_rmse_local(v, model, w)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model) & isfinite(w);
if nnz(valid) < 3
    rmse = inf;
    affine = struct('gain', NaN, 'offsetMv', NaN);
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
beta = (H .* sw) \ (v(valid) .* sw);
pred(valid) = H * beta;
res = v(valid) - pred(valid);
rmse = sqrt(sum(w(valid) .* res.^2) / max(sum(w(valid)), eps));
affine = struct('gain', beta(1), 'offsetMv', beta(2));
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function plot_low_fit_local(sensorCorr, figFile)
fig = figure('Name', 'Step06F low-speed tilt-corrected fit', 'Color', 'w', ...
    'Position', [60, 60, 1450, 340 * numel(sensorCorr)]);
tiledlayout(numel(sensorCorr), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensorCorr)
    S = sensorCorr(i);
    nexttile; hold on; grid on; box on;
    plot(S.x, S.vLowMv, 'k.', 'MarkerSize', 5, 'DisplayName', 'T low');
    plot(S.x, S.vFitMv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'tilt-corrected F raw');
    title(sprintf('CH%d | g0 %.4f | tau %.4f | k %.3f | mu %.3f | RMSE %.2f mV', ...
        S.sensorId, S.g0Mm, S.tauMm, S.xScale, S.muGapPerXMm, S.lowFitRmseMv));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(S.x, S.residualMv, '.', 'MarkerSize', 5);
    yline(0, '--');
    title(sprintf('CH%d residual', S.sensorId));
    xlabel('x_{OPR} (mm)'); ylabel('mV');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_tilt_path_local(sensorCorr, responseSurface, figFile)
fig = figure('Name', 'Step06F tilt path in raw gap library', 'Color', 'w', ...
    'Position', [80, 80, 1300, 310 * numel(sensorCorr)]);
tiledlayout(numel(sensorCorr), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensorCorr)
    S = sensorCorr(i);
    nexttile; hold on; grid on; box on;
    plot(S.xLib, S.gEff, 'r-', 'LineWidth', 1.8, 'DisplayName', 'tilted low-speed path');
    yline(S.g0Mm, 'k--', 'DisplayName', 'mean g0');
    xline(0, ':', 'HandleVisibility', 'off');
    ylim([min(responseSurface.gTrainMm), max(responseSurface.gTrainMm)]);
    xlim([min(responseSurface.xGrid), max(responseSurface.xGrid)]);
    title(sprintf('CH%d raw-library path | tilt %.2f deg | mu %.4f', ...
        S.sensorId, S.tiltAngleDeg, S.muGapPerXMm));
    xlabel('x_{lib} (mm)'); ylabel('g_{eff} (mm)'); legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_high_static_check_local(highCache, figFile)
fig = figure('Name', 'Step06F high-speed static check with tilt-corrected library', ...
    'Color', 'w', 'Position', [80, 80, 1450, 330 * numel(highCache)]);
tiledlayout(numel(highCache), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(highCache)
    C = highCache(i).cache;
    [xSort, order] = sort(C.x);
    nexttile; hold on; grid on; box on;
    plot(C.x, C.v, '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 5, ...
        'DisplayName', 'high-speed samples');
    plot(xSort, C.pred(order), '-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('tilt-corrected F, RMSE %.1f mV', highCache(i).fit.weightedRmseMv));
    title(sprintf('CH%d high-speed static check | g %.4f | dx %.4f mm', ...
        highCache(i).sensorId, highCache(i).fit.gMm, highCache(i).fit.dxCorrectionMm));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

