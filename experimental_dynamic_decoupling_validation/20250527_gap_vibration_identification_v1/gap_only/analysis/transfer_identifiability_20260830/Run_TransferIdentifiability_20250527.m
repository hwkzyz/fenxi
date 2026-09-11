%% Diagnose low-speed transfer identifiability and incremental-operator stability.
% This script is diagnostic only. It does not replace the formal Main07 or
% Main10 outputs and it never modifies the raw or calibrated source data.

clear; clc; close all;
rng(20260830, 'twister');

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
caseDir = fileparts(fileparts(thisDir));
resultDir = fullfile(caseDir, 'results');
outputDir = fullfile(thisDir, 'outputs');
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

inputFile = fullfile(resultDir, ...
    'Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136_nob_harddomain_fixed12pctroi_formal_20260830.mat');
if ~isfile(inputFile)
    error('Formal low-speed calibration file not found: %s', inputFile);
end
S = load(inputFile, 'CorrectedGapLibrary');
C = S.CorrectedGapLibrary;
rs = C.responseSurface;
cfg0 = C.cfg;

cfg = struct();
cfg.lb = [cfg0.gBoundsMm(1), cfg0.tauBoundsMm(1), cfg0.xScaleBounds(1), cfg0.muBounds(1)];
cfg.ub = [cfg0.gBoundsMm(2), cfg0.tauBoundsMm(2), cfg0.xScaleBounds(2), cfg0.muBounds(2)];
cfg.domainToleranceMm = max(cfg0.domainToleranceMm, 1e-8);
cfg.multistartCount = 24;
cfg.cvStartCount = 8;
cfg.profileStartCount = 5;
cfg.profileGridCount = 11;
cfg.profileHalfSpanFraction = 0.20;
cfg.nearOptimalAbsMv = 0.50;
cfg.nearOptimalRel = 0.01;
cfg.operatorDeltaGapMm = 0.01;
cfg.fminOptions = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxIterations', 800, 'MaxFunctionEvaluations', 12000, ...
    'ConstraintTolerance', cfg.domainToleranceMm, ...
    'OptimalityTolerance', 1e-9, 'StepTolerance', 1e-10);

models = struct( ...
    'name', {'P0_g0', 'P1_g0_gain', 'P2_registration', 'P3_full_transfer'}, ...
    'active', {[1], [1], [1 2 3], [1 2 3 4]}, ...
    'profileGain', {false, true, true, true});

nSensor = numel(C.sensor);
nModel = numel(models);
nestedRows = cell(nSensor * nModel, 1);
cvRows = cell(nSensor * nModel, 1);
svdRows = cell(nSensor, 1);
profileRows = {};
bestP3 = cell(nSensor, 1);
rowNested = 0;
rowCv = 0;

fprintf('Transfer identifiability analysis: %s\n', inputFile);
for is = 1:nSensor
    sensor = C.sensor(is);
    x = sensor.x(:);
    v = sensor.vLowMv(:);
    ok = isfinite(x) & isfinite(v);
    x = x(ok);
    v = v(ok);
    pFormal = [sensor.g0Mm, sensor.tauMm, sensor.xScale, sensor.muGapPerXMm];
    fprintf('\nCH%d: %d transfer-ROI points\n', sensor.sensorId, numel(x));

    for im = 1:nModel
        spec = models(im);
        fit = fit_transfer_model_local(x, v, x, rs, cfg, spec, pFormal, cfg.multistartCount);
        rowNested = rowNested + 1;
        nestedRows{rowNested} = fit_row_local(sensor.sensorId, spec.name, fit);
        fprintf('  %-18s RMSE %.4f mV | g0 %.5f tau %.5f k %.5f mu %.5f a %.5f\n', ...
            spec.name, fit.rmseMv, fit.p(1), fit.p(2), fit.p(3), fit.p(4), fit.gain);
        if strcmp(spec.name, 'P3_full_transfer')
            bestP3{is} = fit;
        end

        foldRmse = nan(3, 1);
        edges = round(linspace(1, numel(x) + 1, 4));
        for ifold = 1:3
            test = false(size(x));
            test(edges(ifold):edges(ifold + 1) - 1) = true;
            train = ~test;
            foldFit = fit_transfer_model_local(x(train), v(train), x, rs, cfg, ...
                spec, fit.p, cfg.cvStartCount);
            if foldFit.feasible
                predTest = predict_transfer_local(foldFit.p, foldFit.gain, x(test), rs);
                foldRmse(ifold) = sqrt(mean((v(test) - predTest).^2, 'omitnan'));
            end
        end
        rowCv = rowCv + 1;
        cvRows{rowCv} = table(sensor.sensorId, string(spec.name), ...
            mean(foldRmse, 'omitnan'), std(foldRmse, 'omitnan'), ...
            foldRmse(1), foldRmse(2), foldRmse(3), ...
            'VariableNames', {'sensorId','model','blockedCvMeanRmseMv','blockedCvStdRmseMv', ...
            'fold1RmseMv','fold2RmseMv','fold3RmseMv'});
    end

    fit = bestP3{is};
    [singularValues, effectiveRank, conditionNumber, weakestDirection] = ...
        scaled_jacobian_svd_local(fit, x, rs, cfg);
    svdRows{is} = table(sensor.sensorId, effectiveRank, conditionNumber, ...
        singularValues(1), singularValues(2), singularValues(3), singularValues(4), singularValues(5), ...
        weakestDirection(1), weakestDirection(2), weakestDirection(3), ...
        weakestDirection(4), weakestDirection(5), ...
        'VariableNames', {'sensorId','effectiveRank','scaledConditionNumber', ...
        'sigma1','sigma2','sigma3','sigma4','sigma5', ...
        'weak_g0','weak_tau','weak_k','weak_mu','weak_gain'});

    baseOperator = incremental_operator_local(fit.p, fit.gain, x, rs, cfg.operatorDeltaGapMm);
    span = cfg.ub - cfg.lb;
    for jp = 1:4
        lo = max(cfg.lb(jp), fit.p(jp) - cfg.profileHalfSpanFraction * span(jp));
        hi = min(cfg.ub(jp), fit.p(jp) + cfg.profileHalfSpanFraction * span(jp));
        profileGrid = unique([linspace(lo, hi, cfg.profileGridCount), fit.p(jp)]);
        for ig = 1:numel(profileGrid)
            profileSpec = models(end);
            profileSpec.active = setdiff(1:4, jp, 'stable');
            pSeed = fit.p;
            pSeed(jp) = profileGrid(ig);
            profileSpec.fixedP = pSeed;
            pfit = fit_transfer_model_local(x, v, x, rs, cfg, profileSpec, ...
                pSeed, cfg.profileStartCount);
            if pfit.feasible
                op = incremental_operator_local(pfit.p, pfit.gain, x, rs, cfg.operatorDeltaGapMm);
                opRel = norm(op - baseOperator) / max(norm(baseOperator), eps);
            else
                opRel = NaN;
            end
            profileRows{end+1, 1} = table(sensor.sensorId, jp, string(param_name_local(jp)), ...
                profileGrid(ig), pfit.rmseMv, pfit.p(1), pfit.p(2), pfit.p(3), pfit.p(4), ...
                pfit.gain, pfit.feasible, pfit.gMarginMm, pfit.xMarginMm, opRel, ...
                'VariableNames', {'sensorId','profileParameterIndex','profileParameter', ...
                'fixedValue','rmseMv','g0Mm','tauMm','xScale','muGapPerXMm', ...
                'voltageGain','domainFeasible','gMarginMm','xMarginMm','operatorRelativeL2'}); %#ok<SAGROW>
        end
    end
end

Nested = vertcat(nestedRows{1:rowNested});
BlockedCV = vertcat(cvRows{1:rowCv});
LocalSVD = vertcat(svdRows{:});
Profile = vertcat(profileRows{:});

nearRows = cell(nSensor, 1);
for is = 1:nSensor
    sid = C.sensor(is).sensorId;
    baseRmse = bestP3{is}.rmseMv;
    tol = max(cfg.nearOptimalAbsMv, cfg.nearOptimalRel * baseRmse);
    P = Profile(Profile.sensorId == sid & Profile.domainFeasible & ...
        Profile.rmseMv <= baseRmse + tol, :);
    values = P.operatorRelativeL2(isfinite(P.operatorRelativeL2));
    if isempty(values)
        values = NaN;
    end
    nearRows{is} = table(sid, baseRmse, tol, height(P), ...
        median(values, 'omitnan'), percentile_local(values, 95), max(values, [], 'omitnan'), ...
        'VariableNames', {'sensorId','bestP3RmseMv','nearOptimalToleranceMv', ...
        'nearOptimalProfileCount','operatorMedianRelativeL2','operatorP95RelativeL2', ...
        'operatorMaxRelativeL2'});
end
OperatorStability = vertcat(nearRows{:});

writetable(Nested, fullfile(outputDir, 'NestedPlatformModels_20250527.csv'));
writetable(BlockedCV, fullfile(outputDir, 'BlockedSpatialCV_20250527.csv'));
writetable(LocalSVD, fullfile(outputDir, 'LocalJacobianSVD_20250527.csv'));
writetable(Profile, fullfile(outputDir, 'P3ParameterProfiles_20250527.csv'));
writetable(OperatorStability, fullfile(outputDir, 'IncrementalOperatorStability_20250527.csv'));
save(fullfile(outputDir, 'TransferIdentifiability_20250527.mat'), ...
    'Nested', 'BlockedCV', 'LocalSVD', 'Profile', 'OperatorStability', ...
    'bestP3', 'cfg', 'inputFile');

fig = figure('Color', 'w', 'Position', [80 80 1280 820]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot_grouped_models_local(Nested, 'rmseMv');
ylabel('Fit RMSE (mV)'); title('Low-speed transfer ROI');
nexttile;
plot_grouped_models_local(BlockedCV, 'blockedCvMeanRmseMv');
ylabel('Blocked-CV RMSE (mV)'); title('Spatial blocked validation');
nexttile;
hold on;
for is = 1:height(LocalSVD)
    semilogy(1:5, table2array(LocalSVD(is, {'sigma1','sigma2','sigma3','sigma4','sigma5'})), ...
        '-o', 'LineWidth', 1.3, 'DisplayName', sprintf('CH%d', LocalSVD.sensorId(is)));
end
grid on; box on; xlabel('Scaled Jacobian singular-value index'); ylabel('Singular value');
title('Local transfer identifiability'); legend('Location', 'southwest');
nexttile;
hold on;
colors = lines(nSensor);
for is = 1:nSensor
    sid = C.sensor(is).sensorId;
    P = Profile(Profile.sensorId == sid & Profile.domainFeasible, :);
    scatter(P.rmseMv - bestP3{is}.rmseMv, 100 * P.operatorRelativeL2, 22, ...
        colors(is, :), 'filled', 'DisplayName', sprintf('CH%d', sid));
end
grid on; box on; xlabel('\DeltaRMSE from P3 optimum (mV)');
ylabel('Incremental-operator change (%)'); title(sprintf('\Delta g = %.3f mm', cfg.operatorDeltaGapMm));
legend('Location', 'best');
exportgraphics(fig, fullfile(outputDir, 'TransferIdentifiability_20250527.png'), 'Resolution', 220);
exportgraphics(fig, fullfile(outputDir, 'TransferIdentifiability_20250527.pdf'), 'ContentType', 'vector');
close(fig);

disp(Nested);
disp(BlockedCV);
disp(LocalSVD);
disp(OperatorStability);
fprintf('\nOutputs: %s\n', outputDir);

function row = fit_row_local(sensorId, modelName, fit)
row = table(sensorId, string(modelName), fit.rmseMv, fit.p(1), fit.p(2), fit.p(3), fit.p(4), ...
    fit.gain, fit.feasible, fit.gMarginMm, fit.xMarginMm, fit.activeBoundary, fit.exitFlag, ...
    'VariableNames', {'sensorId','model','rmseMv','g0Mm','tauMm','xScale','muGapPerXMm', ...
    'voltageGain','domainFeasible','gMarginMm','xMarginMm','activeBoundary','exitFlag'});
end

function fit = fit_transfer_model_local(xFit, vFit, xDomain, rs, cfg, spec, pSeed, nStart)
defaults = [median(rs.gTrainMm(:)), 0, 1, 0];
if isfield(spec, 'fixedP')
    defaults = spec.fixedP;
end
active = spec.active;
fixedP = defaults;
lb = cfg.lb(active);
ub = cfg.ub(active);
starts = zeros(max(nStart, 1), numel(active));
starts(1, :) = min(max(pSeed(active), lb), ub);
for i = 2:size(starts, 1)
    starts(i, :) = lb + rand(size(lb)) .* (ub - lb);
end
fit = empty_fit_local();
for i = 1:size(starts, 1)
    objective = @(z) objective_local(expand_p_local(z, active, fixedP), ...
        spec.profileGain, xFit, vFit, rs);
    nonlcon = @(z) domain_constraints_local(expand_p_local(z, active, fixedP), xDomain, rs);
    try
        [z, ~, exitFlag] = fmincon(objective, starts(i, :), [], [], [], [], ...
            lb, ub, nonlcon, cfg.fminOptions);
    catch
        continue;
    end
    p = expand_p_local(z, active, fixedP);
    [rmse, gain] = objective_local(p, spec.profileGain, xFit, vFit, rs);
    [gMargin, xMargin] = domain_margins_local(p, xDomain, rs);
    feasible = exitFlag > 0 && gMargin >= -cfg.domainToleranceMm && xMargin >= -cfg.domainToleranceMm;
    if feasible && rmse < fit.rmseMv
        fit = struct('p', p, 'gain', gain, 'rmseMv', rmse, 'feasible', true, ...
            'gMarginMm', gMargin, 'xMarginMm', xMargin, ...
            'activeBoundary', min(gMargin, xMargin) <= 10 * cfg.domainToleranceMm, ...
            'exitFlag', exitFlag);
    end
end
end

function fit = empty_fit_local()
fit = struct('p', nan(1,4), 'gain', NaN, 'rmseMv', inf, 'feasible', false, ...
    'gMarginMm', NaN, 'xMarginMm', NaN, 'activeBoundary', false, 'exitFlag', NaN);
end

function p = expand_p_local(z, active, fixedP)
p = fixedP;
p(active) = z;
end

function [rmse, gain] = objective_local(p, profileGain, x, v, rs)
raw = raw_transfer_local(p, x, rs);
valid = isfinite(raw) & isfinite(v);
if nnz(valid) < 8
    rmse = inf;
    gain = NaN;
    return;
end
if profileGain
    gain = raw(valid) \ v(valid);
else
    gain = 1;
end
pred = gain .* raw(valid);
rmse = sqrt(mean((v(valid) - pred).^2));
end

function pred = predict_transfer_local(p, gain, x, rs)
pred = gain .* raw_transfer_local(p, x, rs);
end

function raw = raw_transfer_local(p, x, rs)
g = p(1) + p(4) .* (x - p(2));
xLib = p(3) .* (x - p(2));
[xMin, xMax] = response_x_domain_local(rs);
gMin = min(rs.gTrainMm(:));
gMax = max(rs.gTrainMm(:));
% The hard constraints establish physical feasibility. Clipping here only
% removes round-off at an active boundary, matching formal Main07 behavior.
xLib = min(max(xLib, xMin), xMax);
g = min(max(g, gMin), gMax);
B0 = interp1(rs.xGrid(:), rs.coeff(:,1), xLib, 'linear', NaN);
B1 = interp1(rs.xGrid(:), rs.coeff(:,2), xLib, 'linear', NaN);
B2 = interp1(rs.xGrid(:), rs.coeff(:,3), xLib, 'linear', NaN);
raw = B0 + B1 ./ g + B2 .* log(g ./ rs.g0Mm);
end

function [c, ceq] = domain_constraints_local(p, x, rs)
[xMin, xMax] = response_x_domain_local(rs);
gMin = min(rs.gTrainMm(:));
gMax = max(rs.gTrainMm(:));
xLib = p(3) .* (x - p(2));
g = p(1) + p(4) .* (x - p(2));
c = [xMin - xLib; xLib - xMax; gMin - g; g - gMax];
ceq = [];
end

function [gMargin, xMargin] = domain_margins_local(p, x, rs)
[xMin, xMax] = response_x_domain_local(rs);
gMin = min(rs.gTrainMm(:));
gMax = max(rs.gTrainMm(:));
xLib = p(3) .* (x - p(2));
g = p(1) + p(4) .* (x - p(2));
xMargin = min([xLib - xMin; xMax - xLib]);
gMargin = min([g - gMin; gMax - g]);
end

function [xMin, xMax] = response_x_domain_local(rs)
finiteRows = all(isfinite(rs.coeff), 2);
x = rs.xGrid(finiteRows);
xMin = min(x);
xMax = max(x);
end

function op = incremental_operator_local(p, gain, x, rs, dg)
pDyn = p;
pDyn(1) = pDyn(1) + dg;
op = gain .* (raw_transfer_local(pDyn, x, rs) - raw_transfer_local(p, x, rs));
end

function [s, effectiveRank, conditionNumber, weakestDirection] = scaled_jacobian_svd_local(fit, x, rs, cfg)
q = [fit.p, fit.gain];
steps = [1e-5, 1e-5, 1e-5, 1e-6, 1e-5];
J = nan(numel(x), 5);
base = predict_transfer_local(fit.p, fit.gain, x, rs);
for j = 1:5
    qp = q; qm = q;
    qp(j) = qp(j) + steps(j);
    qm(j) = qm(j) - steps(j);
    predP = predict_q_local(qp, x, rs);
    predM = predict_q_local(qm, x, rs);
    feasibleP = q_feasible_local(qp, x, rs, cfg);
    feasibleM = q_feasible_local(qm, x, rs, cfg);
    if feasibleP && feasibleM
        J(:,j) = (predP - predM) ./ (2 * steps(j));
    elseif feasibleP
        J(:,j) = (predP - base) ./ steps(j);
    elseif feasibleM
        J(:,j) = (base - predM) ./ steps(j);
    end
end
scale = [cfg.ub - cfg.lb, max(abs(fit.gain), 1)];
Jscaled = J .* scale;
valid = all(isfinite(Jscaled), 2);
[~, S, V] = svd(Jscaled(valid, :), 0);
s = diag(S).';
s(end+1:5) = 0;
effectiveRank = nnz(s > max(s(1), eps) * 1e-3);
conditionNumber = s(1) / max(s(5), eps);
weakestDirection = V(:, end).';
end

function pred = predict_q_local(q, x, rs)
pred = predict_transfer_local(q(1:4), q(5), x, rs);
end

function tf = q_feasible_local(q, x, rs, cfg)
[gMargin, xMargin] = domain_margins_local(q(1:4), x, rs);
tf = isfinite(q(5)) && gMargin >= -cfg.domainToleranceMm && xMargin >= -cfg.domainToleranceMm;
end

function name = param_name_local(index)
names = {'g0Mm','tauMm','xScale','muGapPerXMm'};
name = names{index};
end

function value = percentile_local(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
    value = NaN;
    return;
end
q = 1 + (numel(x) - 1) * p / 100;
i0 = floor(q); i1 = ceil(q);
if i0 == i1
    value = x(i0);
else
    value = x(i0) + (q - i0) * (x(i1) - x(i0));
end
end

function plot_grouped_models_local(T, valueName)
models = unique(T.model, 'stable');
sensors = unique(T.sensorId, 'stable');
Y = nan(numel(sensors), numel(models));
for is = 1:numel(sensors)
    for im = 1:numel(models)
        row = T.sensorId == sensors(is) & T.model == models(im);
        if any(row)
            Y(is, im) = T.(valueName)(find(row, 1));
        end
    end
end
bar(Y);
grid on; box on;
xticks(1:numel(sensors)); xticklabels(compose('CH%d', sensors));
legend(strrep(cellstr(models), '_', '\_'), 'Location', 'best');
end
