% Fit low-speed no-vibration waveforms with the static gap response surface.
%
% This step does not use the low-speed waveform as the dynamic vibration
% template.  It only estimates a static effective-clearance prior:
%
%   V_low,s(x) = a_s F(g0_s, x - x0_s) + b_s
%
% Step07 then identifies the dynamic vibration while constraining each
% sensor clearance around g0_s.

clear; clc; close all;

%% 1. Paths and parameters
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

cfg = struct();
cfg.targetBlade = 1;
cfg.analysisSensors = [1, 3, 6];
cfg.gapPriorHalfWidthMm = 0.35;
cfg.gapPriorMinHalfWidthMm = 0.15;
cfg.x0LimitMm = 0.80;
cfg.minFitPoints = 80;
cfg.options = optimset('Display', 'off', 'MaxIter', 240, 'MaxFunEvals', 700, ...
    'TolX', 1e-6, 'TolFun', 1e-6);

sensorOverride = strtrim(getenv('STEP06B_ANALYSIS_SENSORS'));
if ~isempty(sensorOverride)
    parsedSensors = sscanf(sensorOverride, '%d').';
    if isempty(parsedSensors)
        error('STEP06B_ANALYSIS_SENSORS must contain integer sensor IDs, for example "1 3 6".');
    end
    cfg.analysisSensors = parsedSensors;
end
halfWidthOverride = str2double(strtrim(getenv('STEP06B_GAP_PRIOR_HALF_WIDTH_MM')));
if isfinite(halfWidthOverride) && halfWidthOverride > 0
    cfg.gapPriorHalfWidthMm = halfWidthOverride;
end
sensorTag = ['S', sprintf('%d', cfg.analysisSensors)];

responseFile = fullfile(outDir, 'Step05_Response_Surface_20250527.mat');
templateDir = fullfile(rootDir, '20250527_low_speed_rotating_calibration', 'output', 'templates');
templateFile = find_low_speed_template_file_local(templateDir, cfg.targetBlade, cfg.analysisSensors);

if ~isfile(responseFile)
    error('Run Step05 first. Missing file: %s', responseFile);
end
if ~isfile(templateFile)
    error('Missing low-speed rotating waveform template file: %s', templateFile);
end

Srf = load(responseFile, 'responseSurface');
responseSurface = Srf.responseSurface;
St = load(templateFile, 'Template');
Template = St.Template;

fprintf('\n=== Step06B: low-speed waveform gap-prior fitting ===\n');
fprintf('Response surface: %s\n', responseFile);
fprintf('Low-speed template: %s\n', templateFile);
fprintf('Sensors: %s, target blade B%d\n', mat2str(cfg.analysisSensors), cfg.targetBlade);
fprintf('Gap-prior half width: %.3f mm\n', cfg.gapPriorHalfWidthMm);

sensorPrior = repmat(struct( ...
    'sensorId', NaN, ...
    'gCenterMm', NaN, ...
    'gLowerMm', NaN, ...
    'gUpperMm', NaN, ...
    'x0Mm', NaN, ...
    'affineGain', NaN, ...
    'affineOffsetMv', NaN, ...
    'plainRmseMv', NaN, ...
    'weightedRmseMv', NaN, ...
    'fitPointCount', NaN, ...
    'templateBaselineV', NaN, ...
    'status', ''), numel(cfg.analysisSensors), 1);
plotCache = repmat(struct('sensorId', NaN, 'x', [], 'vObs', [], 'vFit', []), ...
    numel(cfg.analysisSensors), 1);

for is = 1:numel(cfg.analysisSensors)
    sid = cfg.analysisSensors(is);
    Tsen = get_template_sensor_local(Template, sid);
    [fitResult, cache] = fit_one_sensor_gap_prior_local(Tsen, responseSurface, cfg);

    gLower = max(min(responseSurface.gTrainMm), fitResult.gCenterMm - cfg.gapPriorHalfWidthMm);
    gUpper = min(max(responseSurface.gTrainMm), fitResult.gCenterMm + cfg.gapPriorHalfWidthMm);
    if (gUpper - gLower) < 2 * cfg.gapPriorMinHalfWidthMm
        gLower = max(min(responseSurface.gTrainMm), fitResult.gCenterMm - cfg.gapPriorMinHalfWidthMm);
        gUpper = min(max(responseSurface.gTrainMm), fitResult.gCenterMm + cfg.gapPriorMinHalfWidthMm);
    end

    sensorPrior(is).sensorId = sid;
    sensorPrior(is).gCenterMm = fitResult.gCenterMm;
    sensorPrior(is).gLowerMm = gLower;
    sensorPrior(is).gUpperMm = gUpper;
    sensorPrior(is).x0Mm = fitResult.x0Mm;
    sensorPrior(is).affineGain = fitResult.affineGain;
    sensorPrior(is).affineOffsetMv = fitResult.affineOffsetMv;
    sensorPrior(is).plainRmseMv = fitResult.plainRmseMv;
    sensorPrior(is).weightedRmseMv = fitResult.weightedRmseMv;
    sensorPrior(is).fitPointCount = fitResult.fitPointCount;
    sensorPrior(is).templateBaselineV = Tsen.baseline;
    sensorPrior(is).status = fitResult.status;
    plotCache(is) = cache;

    fprintf('  CH%d: g0 %.4f mm, prior [%.4f, %.4f] mm, x0 %.4f mm, RMSE %.3f mV\n', ...
        sid, sensorPrior(is).gCenterMm, sensorPrior(is).gLowerMm, ...
        sensorPrior(is).gUpperMm, sensorPrior(is).x0Mm, sensorPrior(is).weightedRmseMv);
end

GapPrior = struct();
GapPrior.dataset = '20250527';
GapPrior.route = 'low_speed_waveform_gap_prior_for_dynamic_decoupling';
GapPrior.description = ['Low-speed no-vibration waveform is fitted by the proposed static gap response surface ' ...
    'to estimate per-sensor effective-clearance priors. The low-speed waveform is not used as the dynamic template.'];
GapPrior.targetBlade = cfg.targetBlade;
GapPrior.analysisSensors = cfg.analysisSensors(:).';
GapPrior.sensorTag = sensorTag;
GapPrior.responseFile = responseFile;
GapPrior.templateFile = templateFile;
GapPrior.lowSpeedTemplateRoute = Template.Route;
GapPrior.gapPriorHalfWidthMm = cfg.gapPriorHalfWidthMm;
GapPrior.gapPriorMinHalfWidthMm = cfg.gapPriorMinHalfWidthMm;
GapPrior.x0LimitMm = cfg.x0LimitMm;
GapPrior.responseGapTrainMm = responseSurface.gTrainMm(:);
GapPrior.sensorPrior = sensorPrior;

priorTable = struct2table(sensorPrior);
matFile = fullfile(outDir, sprintf('Step06B_LowSpeed_GapPrior_20250527_B%d_%s.mat', ...
    cfg.targetBlade, sensorTag));
csvFile = fullfile(outDir, sprintf('Step06B_LowSpeed_GapPrior_20250527_B%d_%s.csv', ...
    cfg.targetBlade, sensorTag));
figFile = fullfile(outDir, sprintf('Step06B_LowSpeed_GapPrior_Fit_20250527_B%d_%s.png', ...
    cfg.targetBlade, sensorTag));
save(matFile, 'GapPrior', 'priorTable', '-v7.3');
writetable(priorTable, csvFile);

fig = figure('Name', 'Step06B low-speed gap-prior fit', ...
    'Color', 'w', 'Position', [80, 80, 1450, 320 * numel(cfg.analysisSensors)], ...
    'NumberTitle', 'off');
tiledlayout(numel(cfg.analysisSensors), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for is = 1:numel(plotCache)
    nexttile; hold on; grid on; box on;
    plot(plotCache(is).x, plotCache(is).vObs, '.', 'Color', [0.25, 0.45, 0.85], ...
        'MarkerSize', 5, 'DisplayName', 'low-speed waveform');
    plot(plotCache(is).x, plotCache(is).vFit, '-', 'Color', [0.85, 0.25, 0.15], ...
        'LineWidth', 1.6, 'DisplayName', 'response-surface fit');
    title(sprintf('CH%d | g0 %.4f mm | RMSE %.2f mV', ...
        sensorPrior(is).sensorId, sensorPrior(is).gCenterMm, sensorPrior(is).weightedRmseMv));
    xlabel('x (mm)');
    ylabel('Baseline-corrected voltage (mV)');
    legend('Location', 'best');
end
saveas(fig, figFile);

fprintf('\nStep06B complete.\n');
fprintf('Gap prior saved to:\n  %s\n', matFile);
disp(priorTable);

function templateFile = find_low_speed_template_file_local(templateDir, targetBlade, sensorIds)
exactTag = ['S', sprintf('%d', sensorIds)];
exactFile = fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_%s_20250527.mat', ...
    targetBlade, exactTag));
if isfile(exactFile)
    templateFile = exactFile;
    return;
end
files = dir(fullfile(templateDir, sprintf('Template_LowSpeedRotating_B%d_S*_20250527.mat', targetBlade)));
for i = 1:numel(files)
    token = regexp(files(i).name, '_S([0-9]+)_20250527\.mat$', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    available = arrayfun(@(c) str2double(c), token{1});
    if all(ismember(sensorIds, available))
        templateFile = fullfile(files(i).folder, files(i).name);
        return;
    end
end
templateFile = exactFile;
end

function Tsen = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1, 'first');
if isempty(idx)
    error('Low-speed template does not contain CH%d.', sid);
end
Tsen = Template.Sensor(idx);
end

function [fitResult, cache] = fit_one_sensor_gap_prior_local(Tsen, responseSurface, cfg)
x = Tsen.x_grid(:);
vMv = (Tsen.v_grid(:) - Tsen.baseline) * 1000;
mask = isfinite(x) & isfinite(vMv);
mask = mask & x >= min(responseSurface.xGrid) + cfg.x0LimitMm & ...
    x <= max(responseSurface.xGrid) - cfg.x0LimitMm;
x = x(mask);
vMv = vMv(mask);
if isfield(Tsen, 'bin_weight') && numel(Tsen.bin_weight) == numel(mask)
    w = Tsen.bin_weight(:);
    w = w(mask);
else
    w = ones(size(x));
end
w(~isfinite(w) | w <= 0) = 1;
w = w ./ mean(w, 'omitnan');

if numel(x) < cfg.minFitPoints
    error('CH%d has only %d usable low-speed template points.', Tsen.sensor_id, numel(x));
end

gTrain = responseSurface.gTrainMm(:).';
gMin = min(gTrain);
gMax = max(gTrain);
startsG = unique([gTrain, median(gTrain), mean(gTrain)]);
startsX0 = [-0.20, 0, 0.20];
best = struct('objective', inf);
for g0 = startsG
    for x00 = startsX0
        p0 = [g0, x00];
        objFun = @(p) objective_gap_prior_local(p, x, vMv, w, responseSurface, gMin, gMax, cfg.x0LimitMm);
        pOpt = fminsearch(objFun, p0, cfg.options);
        [obj, vFit, affine, diagInfo] = objective_gap_prior_local(pOpt, x, vMv, w, responseSurface, gMin, gMax, cfg.x0LimitMm);
        if obj < best.objective
            best = struct('objective', obj, 'params', pOpt, 'vFit', vFit, ...
                'affine', affine, 'diagInfo', diagInfo);
        end
    end
end

fitResult = struct();
fitResult.gCenterMm = best.diagInfo.gMm;
fitResult.x0Mm = best.diagInfo.x0Mm;
fitResult.affineGain = best.affine.gain;
fitResult.affineOffsetMv = best.affine.offset;
fitResult.plainRmseMv = best.diagInfo.plainRmseMv;
fitResult.weightedRmseMv = best.diagInfo.weightedRmseMv;
fitResult.fitPointCount = numel(x);
fitResult.status = 'ok';

cache = struct('sensorId', Tsen.sensor_id, 'x', x, 'vObs', vMv, 'vFit', best.vFit);
end

function [obj, vPred, affine, diagInfo] = objective_gap_prior_local(p, x, vObs, w, responseSurface, gMin, gMax, x0LimitMm)
gRaw = p(1);
x0Raw = p(2);
g = min(max(gRaw, gMin), gMax);
x0 = min(max(x0Raw, -x0LimitMm), x0LimitMm);
penalty = 1e12 * ((gRaw - g)^2 + (x0Raw - x0)^2);
modelRaw = eval_response_surface_local(responseSurface, g, x - x0);
valid = isfinite(modelRaw) & isfinite(vObs) & isfinite(w);
if nnz(valid) < 10
    obj = inf;
    vPred = nan(size(vObs));
    affine = struct('gain', NaN, 'offset', NaN);
    diagInfo = struct('gMm', g, 'x0Mm', x0, 'plainRmseMv', inf, 'weightedRmseMv', inf);
    return;
end
A = [modelRaw(valid), ones(nnz(valid), 1)];
ws = sqrt(max(w(valid), eps));
coeff = (A .* ws) \ (vObs(valid) .* ws);
vPred = nan(size(vObs));
vPred(valid) = A * coeff;
res = vObs(valid) - vPred(valid);
obj = sum(w(valid) .* res.^2) + penalty;
affine = struct('gain', coeff(1), 'offset', coeff(2));
diagInfo = struct('gMm', g, 'x0Mm', x0, ...
    'plainRmseMv', sqrt(mean(res.^2, 'omitnan')), ...
    'weightedRmseMv', sqrt(sum(w(valid) .* res.^2) / sum(w(valid))));
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 1), xq, 'linear', NaN);
B1 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 2), xq, 'linear', NaN);
B2 = interp1(responseSurface.xGrid, responseSurface.coeff(:, 3), xq, 'linear', NaN);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end
