%% Step06D: Build OPR-anchored gap library from low-speed template
% The raw static gap library has its own calibration x-axis, but no OPR
% reference. The low-speed rotating template has an OPR-referenced x-axis.
% This step estimates a per-sensor x-axis registration:
%
%   x_lib = x_OPR - tau_s
%   F_OPR,s(g, x_OPR) = F_raw(g, x_OPR - tau_s)
%
% The low-speed template is used for x-axis registration only, not as a
% high-speed waveform residual. The script saves diagnostics and figures so
% the registration failure mode is visible.

clear; clc;

%% 1. Paths and settings
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
figDir = fullfile(outDir, 'figures_step06d_opr_anchored_gap_library');
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.gapGridCount = 61;
cfg.tauGridMm = -0.80:0.005:0.80;
cfg.minFitPoints = 80;
cfg.highCheckMaxPointsPerSensor = 320;
cfg.highCheckTauCorrectionGridMm = -0.08:0.01:0.08;
cfg.registrationObjective = 'affine_voltage';

sensorOverride = strtrim(getenv('STEP06D_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06D_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
highMapFile = fullfile(outDir, sprintf('Step06_HighMap_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
templateDir = fullfile(rootDir, '20250527_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(templateFile)
    error('Missing low-speed template file: %s', templateFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
St = load(templateFile, 'Template');
Template = St.Template;
hasHighMap = isfile(highMapFile);
if hasHighMap
    H = load(highMapFile, 'highMap');
    highMap = attach_model_coordinate_local(H.highMap);
else
    highMap = struct();
end

fprintf('\n=== Step06D: OPR-anchored gap library registration ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Sensors: %s\n', mat2str(cfg.analysisSensors));

gGrid = linspace(min(responseSurface.gTrainMm), max(responseSurface.gTrainMm), cfg.gapGridCount);
sensorLib = repmat(struct('sensorId', NaN, 'tauMm', NaN, 'g0Mm', NaN, ...
    'objectiveRmseMv', NaN, 'shapeRmse', NaN, 'affineRmseMv', NaN, 'corr', NaN, ...
    'pointCount', NaN, 'xOpr', [], 'vLowMv', [], 'vFitMv', [], ...
    'dLowNorm', [], 'dFitNorm', [], 'shapeMap', []), numel(cfg.analysisSensors), 1);
rows = cell(numel(cfg.analysisSensors), 1);
highRows = {};
highCache = struct([]);

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    fit = fit_opr_registration_shape_local(Tsen, responseSurface, gGrid, cfg.tauGridMm, cfg);
    sensorLib(is).sensorId = sid;
    sensorLib(is).tauMm = fit.tauMm;
    sensorLib(is).g0Mm = fit.g0Mm;
    sensorLib(is).objectiveRmseMv = fit.objectiveRmseMv;
    sensorLib(is).shapeRmse = fit.shapeRmse;
    sensorLib(is).affineRmseMv = fit.affineRmseMv;
    sensorLib(is).corr = fit.corr;
    sensorLib(is).pointCount = numel(fit.x);
    sensorLib(is).xOpr = fit.x;
    sensorLib(is).vLowMv = fit.vLow;
    sensorLib(is).vFitMv = fit.vFit;
    sensorLib(is).dLowNorm = fit.dLowNorm;
    sensorLib(is).dFitNorm = fit.dFitNorm;
    sensorLib(is).shapeMap = fit.shapeMap;
    rows{is} = table(sid, fit.tauMm, fit.g0Mm, fit.objectiveRmseMv, fit.shapeRmse, ...
        fit.affineRmseMv, fit.corr, numel(fit.x), ...
        'VariableNames', {'sensorId','tauMm','g0Mm','objectiveRmseMv','shapeRmse', ...
        'affineRmseMv','shapeCorrelation','pointCount'});

    fprintf('  CH%d: tau %.4f mm, g0 %.4f mm, objective %.2f mV, shape RMSE %.4f, affine RMSE %.2f mV, corr %.4f\n', ...
        sid, fit.tauMm, fit.g0Mm, fit.objectiveRmseMv, fit.shapeRmse, fit.affineRmseMv, fit.corr);

    if hasHighMap
        [hfit, hcache] = high_static_check_local(highMap, sid, responseSurface, ...
            fit.tauMm, gGrid, cfg.highCheckTauCorrectionGridMm, cfg.highCheckMaxPointsPerSensor);
        highRows{end+1, 1} = table(sid, hfit.gMm, hfit.tauCorrectionMm, ...
            fit.tauMm + hfit.tauCorrectionMm, hfit.weightedRmseMv, hfit.pointCount, ...
            'VariableNames', {'sensorId','highStaticGMm','tauCorrectionMm', ...
            'tauTotalMm','weightedRmseMv','pointCount'}); %#ok<SAGROW>
        highCache(end+1).sensorId = sid; %#ok<SAGROW>
        highCache(end).fit = hfit;
        highCache(end).cache = hcache;
    end
end

registrationTable = vertcat(rows{:});
if isempty(highRows)
    highCheckTable = table();
else
    highCheckTable = vertcat(highRows{:});
end

CompleteGapLibrary = struct();
CompleteGapLibrary.dataset = '20250527';
CompleteGapLibrary.method = 'low_speed_template_opr_anchored_gap_library';
CompleteGapLibrary.description = ['Raw gap response surface with per-sensor OPR x-axis registration estimated ' ...
    'from the low-speed rotating template by affine voltage-profile matching. The fitted affine gain/offset are ' ...
    'used only as nuisance alignment terms; the low-speed waveform residual is not embedded in the library.'];
CompleteGapLibrary.responseFile = responseFile;
CompleteGapLibrary.templateFile = templateFile;
CompleteGapLibrary.targetBlade = cfg.targetBlade;
CompleteGapLibrary.analysisSensors = cfg.analysisSensors(:).';
CompleteGapLibrary.responseSurface = responseSurface;
CompleteGapLibrary.sensor = sensorLib;
CompleteGapLibrary.xRelation = 'x_lib = x_OPR - tauMm';
CompleteGapLibrary.cfg = cfg;

matFile = fullfile(outDir, sprintf('Step06D_OPR_Anchored_GapLibrary_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step06D_OPR_Anchored_GapLibrary_Registration_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
highCsvFile = fullfile(outDir, sprintf('Step06D_OPR_Anchored_GapLibrary_HighStaticCheck_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figOverlay = fullfile(figDir, sprintf('Step06D_LowTemplate_GapLibrary_OPR_Overlay_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figHeat = fullfile(figDir, sprintf('Step06D_Registration_Objective_Heatmap_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
figHigh = fullfile(figDir, sprintf('Step06D_HighStatic_OPRAnchored_Check_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));

save(matFile, 'CompleteGapLibrary', 'registrationTable', 'highCheckTable', ...
    'responseFile', 'templateFile', 'highMapFile', '-v7.3');
writetable(registrationTable, csvFile);
if ~isempty(highCheckTable)
    writetable(highCheckTable, highCsvFile);
end
plot_registration_overlay_local(sensorLib, figOverlay);
plot_registration_heatmap_local(sensorLib, gGrid, cfg.tauGridMm, figHeat);
if hasHighMap && ~isempty(highCache)
    plot_high_static_check_local(highCache, figHigh);
end

fprintf('\nStep06D complete.\n');
disp(registrationTable);
if ~isempty(highCheckTable)
    disp(highCheckTable);
end
fprintf('Saved:\n  %s\n  %s\n  %s\n', matFile, csvFile, figOverlay);

%% Local functions
function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
patterns = {
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_20250527.mat', targetBlade, exactTag)
    sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', targetBlade, exactTag)
    };
for i = 1:numel(patterns)
    candidate = fullfile(templateDir, patterns{i});
    if isfile(candidate)
        templateFile = candidate;
        return;
    end
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20250527.mat', targetBlade)));
if isempty(files)
    files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_GradientXRange030_20250527.mat', targetBlade)));
end
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)(?:_GradientXRange030)?_20250527\.mat$', 'tokens', 'once');
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

function fit = fit_opr_registration_shape_local(Tsen, responseSurface, gGrid, tauGrid, cfg)
xAll = Tsen.x_grid(:);
vAll = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
mask = isfinite(xAll) & isfinite(vAll);
mask = mask & xAll >= min(responseSurface.xGrid) + max(abs(tauGrid)) & ...
    xAll <= max(responseSurface.xGrid) - max(abs(tauGrid));
x = xAll(mask);
vLow = vAll(mask);
if numel(x) < cfg.minFitPoints
    error('CH%d has too few template points for OPR registration.', Tsen.sensor_id);
end
dLow = gradient(vLow, x);
dLowNorm = normalize_shape_local(dLow);

shapeMap = nan(numel(gGrid), numel(tauGrid));
best = struct('objectiveRmseMv', inf);
for ig = 1:numel(gGrid)
    g = gGrid(ig);
    for it = 1:numel(tauGrid)
        tau = tauGrid(it);
        vFitRaw = eval_response_surface_local(responseSurface, g, x - tau);
        if nnz(isfinite(vFitRaw)) < cfg.minFitPoints
            continue;
        end
        dFit = gradient(vFitRaw, x);
        dFitNorm = normalize_shape_local(dFit);
        valid = isfinite(dLowNorm) & isfinite(dFitNorm);
        if nnz(valid) < cfg.minFitPoints
            continue;
        end
        shapeRmse = sqrt(mean((dLowNorm(valid) - dFitNorm(valid)).^2, 'omitnan'));
        [affineRmse, vFitAffine] = affine_rmse_local(vLow, vFitRaw, ones(size(vLow)));
        objectiveRmse = affineRmse;
        shapeMap(ig, it) = objectiveRmse;
        if objectiveRmse < best.objectiveRmseMv
            c = corr(dLowNorm(valid), dFitNorm(valid), 'Rows', 'complete');
            best = struct('sensorId', Tsen.sensor_id, 'tauMm', tau, 'g0Mm', g, ...
                'objectiveRmseMv', objectiveRmse, 'shapeRmse', shapeRmse, ...
                'affineRmseMv', affineRmse, 'corr', c, 'x', x, 'vLow', vLow, 'vFit', vFitAffine, ...
                'dLowNorm', dLowNorm, 'dFitNorm', dFitNorm, 'shapeMap', shapeMap);
        end
    end
end
fit = best;
fit.shapeMap = shapeMap;
end

function [fit, cache] = high_static_check_local(highMap, sid, responseSurface, tauMm, gGrid, tauCorrGrid, maxPoints)
mask = highMap.S_v == sid & isfinite(highMap.x_model_v) & isfinite(highMap.V_a) & isfinite(highMap.W_v);
idx = find(mask);
if numel(idx) > maxPoints
    idx = idx(round(linspace(1, numel(idx), maxPoints)));
end
x = highMap.x_model_v(idx);
v = highMap.V_a(idx);
w = max(highMap.W_v(idx), 0.05);
fit = struct('gMm', NaN, 'tauCorrectionMm', NaN, 'weightedRmseMv', inf, ...
    'pointCount', 0, 'pred', []);
for g = gGrid(:).'
    for dtau = tauCorrGrid(:).'
        model = eval_response_surface_local(responseSurface, g, x - tauMm - dtau);
        [rmse, pred] = affine_rmse_local(v, model, w);
        if rmse < fit.weightedRmseMv
            fit = struct('gMm', g, 'tauCorrectionMm', dtau, ...
                'weightedRmseMv', rmse, 'pointCount', nnz(isfinite(pred)), 'pred', pred);
        end
    end
end
cache = struct('x', x, 'v', v, 'w', w, 'pred', fit.pred);
end

function y = normalize_shape_local(x)
valid = isfinite(x);
y = nan(size(x));
if nnz(valid) < 3
    return;
end
mu = mean(x(valid), 'omitnan');
sg = std(x(valid), 'omitnan');
if ~isfinite(sg) || sg < eps
    sg = 1;
end
y(valid) = (x(valid) - mu) ./ sg;
end

function [rmse, pred] = affine_rmse_local(v, model, w)
pred = nan(size(v));
valid = isfinite(v) & isfinite(model) & isfinite(w);
if nnz(valid) < 3
    rmse = inf;
    return;
end
H = [model(valid), ones(nnz(valid), 1)];
sw = sqrt(max(w(valid), eps));
beta = (H .* sw) \ (v(valid) .* sw);
pred(valid) = H * beta;
res = v(valid) - pred(valid);
rmse = sqrt(sum(w(valid) .* res.^2) / max(sum(w(valid)), eps));
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function plot_registration_overlay_local(sensorLib, figFile)
fig = figure('Name', 'Step06D OPR anchored gap-library overlay', 'Color', 'w', ...
    'Position', [60, 60, 1450, 360 * numel(sensorLib)]);
tiledlayout(numel(sensorLib), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensorLib)
    S = sensorLib(i);
    nexttile; hold on; grid on; box on;
    plot(S.xOpr, normalize_shape_local(S.vLowMv), '.', 'MarkerSize', 5, 'DisplayName', 'T low z');
    plot(S.xOpr, normalize_shape_local(S.vFitMv), '-', 'LineWidth', 1.4, 'DisplayName', 'F aligned z');
    title(sprintf('CH%d voltage shape | tau %.4f mm | g0 %.4f mm', S.sensorId, S.tauMm, S.g0Mm));
    xlabel('x_{OPR} (mm)'); ylabel('z-score'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(S.xOpr, S.dLowNorm, '.', 'MarkerSize', 5, 'DisplayName', 'dT/dx z');
    plot(S.xOpr, S.dFitNorm, '-', 'LineWidth', 1.4, 'DisplayName', 'dF/dx z');
    title(sprintf('Derivative registration | RMSE %.4f | corr %.4f', S.shapeRmse, S.corr));
    xlabel('x_{OPR} (mm)'); ylabel('normalized derivative'); legend('Location', 'best');

    nexttile; hold on; grid on; box on;
    plot(S.xOpr, S.vLowMv - S.vFitMv, '.', 'MarkerSize', 5);
    yline(0, '--');
    title(sprintf('Affine residual | RMSE %.2f mV', S.affineRmseMv));
    xlabel('x_{OPR} (mm)'); ylabel('mV');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_registration_heatmap_local(sensorLib, gGrid, tauGrid, figFile)
fig = figure('Name', 'Step06D registration objective heatmap', 'Color', 'w', ...
    'Position', [80, 80, 1450, 330 * numel(sensorLib)]);
tiledlayout(numel(sensorLib), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sensorLib)
    S = sensorLib(i);
    nexttile; hold on; box on;
    imagesc(tauGrid, gGrid, S.shapeMap);
    set(gca, 'YDir', 'normal');
    colorbar;
    plot(S.tauMm, S.g0Mm, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
    xlabel('\tau_s (mm): x_{lib}=x_{OPR}-\tau_s');
    ylabel('g_0 (mm)');
    title(sprintf('CH%d derivative-shape registration objective', S.sensorId));
end
exportgraphics(fig, figFile, 'Resolution', 220);
end

function plot_high_static_check_local(highCache, figFile)
fig = figure('Name', 'Step06D high-speed static check with OPR anchored library', ...
    'Color', 'w', 'Position', [80, 80, 1450, 330 * numel(highCache)]);
tiledlayout(numel(highCache), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(highCache)
    C = highCache(i).cache;
    [xSort, order] = sort(C.x);
    nexttile; hold on; grid on; box on;
    plot(C.x, C.v, '.', 'Color', [0.35 0.35 0.35], 'MarkerSize', 5, ...
        'DisplayName', 'high-speed samples');
    plot(xSort, C.pred(order), '-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('OPR anchored F, RMSE %.1f mV', highCache(i).fit.weightedRmseMv));
    title(sprintf('CH%d high-speed static check | g %.4f | delta tau %.4f mm', ...
        highCache(i).sensorId, highCache(i).fit.gMm, highCache(i).fit.tauCorrectionMm));
    xlabel('x_{OPR} (mm)'); ylabel('mV'); legend('Location', 'best');
end
exportgraphics(fig, figFile, 'Resolution', 220);
end
