function Result = Run_01_StaticCalibrationTransferAudit()
%RUN_01_STATICCALIBRATIONTRANSFERAUDIT Audit FE tilt calibration transfer.
% The frozen zero-tilt library is used for C1/C2 inversion. Tilted COMSOL
% surfaces provide independent low- and high-clearance truth waveforms.

root = fileparts(mfilename('fullpath'));
mainDir = fileparts(root);
addpath(mainDir, '-begin');
addpath(fullfile(mainDir, 'local_func'), '-begin');
addpath(fullfile(mainDir, 'identifiability_mechanism_analysis', 'local_func'));

ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
auditGridN = 161;
lib = make_fixed_trust_template_library(ctx.gapList, ctx.xCell, ctx.yCell, ...
    NaN, auditGridN, get_inv_log_2_model_def(), trust);

angles = [0.5 1 1.5 2 2.5 3 3.5];
gapFE = (0.5:0.2:1.5).';
gapLowList = gapFE(1:end-1);
x = lib.xGrid(:);
outDir = fullfile(root, 'output', 'static_calibration_transfer');
if ~exist(outDir, 'dir'), mkdir(outDir); end

lowRows = repmat(low_row(), 0, 1);
transferRows = repmat(transfer_row(), 0, 1);
sourceRows = repmat(source_row(), 0, 1);

for ia = 1:numel(angles)
    angle = angles(ia);
    [filePath, sourceKind] = resolve_tilt_file(ctx.rootDir, angle, gapFE);
    [gTilt, xCell, yCell] = load_stacked_gap_curves(filePath, gapFE);
    tilt = make_response_template_library(gTilt, xCell, yCell, NaN, ...
        auditGridN, get_inv_log_2_model_def(), struct('enable', false));
    tilt = restrict_gap_template_domain(tilt, trust.domain);
    sourceRows(end+1) = struct('tilt_deg', angle, ...
        'source_kind', string(sourceKind), 'file_path', string(filePath)); %#ok<AGROW>

    for ig = 1:numel(gapLowList)
        gLow = gapLowList(ig);
        gHigh = gapFE(ig+1);
        deltaGap = gHigh - gLow;
        yLow = eval_gap_template(tilt, gLow, x);
        yHigh = eval_gap_template(tilt, gHigh, x);
        trueIncrement = yHigh - yLow;
        low = struct('xGrid', x, 'templateLow', yLow);

        opts1 = calibration_options(lib, gLow, [true false false false]);
        [p1, fit1] = calibrate_low_speed_path(lib, low, opts1);
        model1 = build_path_template_model(lib, yLow, p1);
        predIncrement1 = predicted_increment(model1, deltaGap, x);
        metrics1 = transfer_metrics(trueIncrement, predIncrement1, model1, deltaGap, x);
        lowRows(end+1) = make_low_row(angle, gLow, "C1_g0_only", p1, fit1, lib, x); %#ok<AGROW>
        transferRows(end+1) = make_transfer_row(angle, gLow, gHigh, ...
            "C1_g0_only", p1, metrics1); %#ok<AGROW>

        opts2 = calibration_options(lib, gLow, [true true false false]);
        [p2, fit2] = calibrate_low_speed_path(lib, low, opts2);
        model2 = build_path_template_model(lib, yLow, p2);
        predIncrement2 = predicted_increment(model2, deltaGap, x);
        metrics2 = transfer_metrics(trueIncrement, predIncrement2, model2, deltaGap, x);
        lowRows(end+1) = make_low_row(angle, gLow, "C2_g0_mu", p2, fit2, lib, x); %#ok<AGROW>
        transferRows(end+1) = make_transfer_row(angle, gLow, gHigh, ...
            "C2_g0_mu", p2, metrics2); %#ok<AGROW>
    end
end

Low = struct2table(lowRows);
Transfer = struct2table(transferRows);
Sources = struct2table(sourceRows);
Summary = summarize_angles(Transfer, Low, angles, diff(trust.domain));

writetable(Low, fullfile(outDir, 'low_speed_calibration.csv'));
writetable(Transfer, fullfile(outDir, 'static_increment_transfer.csv'));
writetable(Summary, fullfile(outDir, 'angle_summary.csv'));
writetable(Sources, fullfile(outDir, 'data_sources.csv'));
save(fullfile(outDir, 'static_calibration_transfer.mat'), ...
    'Low', 'Transfer', 'Summary', 'Sources', 'angles', 'gapFE', 'trust', '-v7.3');
make_figure(Summary, outDir);
write_report(Summary, Sources, outDir);

Result = struct('low_speed', Low, 'transfer', Transfer, ...
    'summary', Summary, 'sources', Sources, 'output_dir', string(outDir));
disp(Summary);
end

function opts = calibration_options(lib, gLow, fitMask)
opts = struct('g0Init', gLow, 'muInit', 0, 'tauInit', 0, 'zetaInit', 1, ...
    'g0Lower', min(lib.gapTrain), 'g0Upper', max(lib.gapTrain), ...
    'muLower', -0.15, 'muUpper', 0.15, ...
    'tauLower', 0, 'tauUpper', 0, 'zetaLower', 1, 'zetaUpper', 1, ...
    'fitMask', fitMask, 'fitGainOffset', false, 'baseline', 0, ...
    'kappa', 1, 'xDomain', lib.domain, 'maxIter', 500);
end

function y = predicted_increment(model, deltaGap, x)
y = eval_path_increment_template(model, deltaGap, x, 1) - ...
    eval_path_increment_template(model, 0, x, 1);
end

function m = transfer_metrics(yTrue, yPred, model, deltaGap, x)
b = yTrue(:) - yPred(:);
scale = max(norm(yTrue), eps);
h = max(1e-4, abs(deltaGap)*1e-3);
jGap = (predicted_increment(model, deltaGap+h, x) - ...
    predicted_increment(model, deltaGap-h, x))/(2*h);
yHigh = eval_path_increment_template(model, deltaGap, x, 1);
jVib = -gradient(yHigh, median(diff(x)));
m = struct('rmse_V', sqrt(mean(b.^2)), ...
    'relative_norm', norm(b)/scale, ...
    'gap_projection', abs_projection(b, jGap), ...
    'vibration_projection', abs_projection(b, jVib), ...
    'true_increment_rms_V', sqrt(mean(yTrue.^2)), ...
    'pred_increment_rms_V', sqrt(mean(yPred.^2)));
end

function v = abs_projection(a, b)
v = abs(a(:)'*b(:))/max(norm(a)*norm(b), eps);
end

function row = make_low_row(angle, gLow, method, p, fit, lib, x)
[conditionNumber, columnCorrelation] = calibration_geometry(lib, p, x, method);
r = fit.y(:) - fit.fitted(:);
tx = gradient(fit.fitted(:), median(diff(x)));
row = struct('tilt_deg', angle, 'g_low_mm', gLow, 'method', method, ...
    'g0_hat_mm', p.g0, 'mu_hat', p.mu, 'tau_fixed_mm', p.tau, ...
    'zeta_fixed', p.zeta, 'low_rmse_V', fit.rmse_V, ...
    'calibration_condition_number', conditionNumber, ...
    'g0_mu_column_correlation', columnCorrelation, ...
    'translation_residual_correlation', abs_projection(r, tx));
end

function [kappa, rho] = calibration_geometry(lib, p, x, method)
if method == "C1_g0_only"
    kappa = 1;
    rho = NaN;
    return;
end
hg = 1e-4;
hm = 1e-5;
jg = (eval_gap_path_response(lib, p.g0+hg, p.mu, 0, 1, x) - ...
    eval_gap_path_response(lib, p.g0-hg, p.mu, 0, 1, x))/(2*hg);
jm = (eval_gap_path_response(lib, p.g0, p.mu+hm, 0, 1, x) - ...
    eval_gap_path_response(lib, p.g0, p.mu-hm, 0, 1, x))/(2*hm);
J = [jg(:) jm(:)];
Jn = J ./ max(vecnorm(J), eps);
kappa = cond(Jn);
rho = abs(Jn(:,1)'*Jn(:,2));
end

function row = make_transfer_row(angle, gLow, gHigh, method, p, m)
row = struct('tilt_deg', angle, 'g_low_mm', gLow, 'g_high_mm', gHigh, ...
    'delta_gap_mm', gHigh-gLow, 'method', method, 'g0_hat_mm', p.g0, ...
    'mu_hat', p.mu, 'transfer_rmse_V', m.rmse_V, ...
    'transfer_relative_norm', m.relative_norm, ...
    'gap_direction_projection', m.gap_projection, ...
    'vibration_direction_projection', m.vibration_projection, ...
    'true_increment_rms_V', m.true_increment_rms_V, ...
    'pred_increment_rms_V', m.pred_increment_rms_V);
end

function S = summarize_angles(T, L, angles, trustLength)
rows = repmat(summary_row(), numel(angles), 1);
for i = 1:numel(angles)
    a = angles(i);
    q1 = T.tilt_deg==a & T.method=="C1_g0_only";
    q2 = T.tilt_deg==a & T.method=="C2_g0_mu";
    ql1 = L.tilt_deg==a & L.method=="C1_g0_only";
    ql2 = L.tilt_deg==a & L.method=="C2_g0_mu";
    b1 = T.transfer_relative_norm(q1);
    b2 = T.transfer_relative_norm(q2);
    mu = median(L.mu_hat(ql2));
    rows(i) = struct('tilt_deg', a, 'mu_hat_median', mu, ...
        'tilt_severity_mm', abs(mu)*trustLength, ...
        'low_rmse_C1_median_V', median(L.low_rmse_V(ql1)), ...
        'low_rmse_C2_median_V', median(L.low_rmse_V(ql2)), ...
        'transfer_relative_C1_median', median(b1), ...
        'transfer_relative_C2_median', median(b2), ...
        'transfer_relative_C1_max', max(b1), ...
        'transfer_relative_C2_max', max(b2), ...
        'C2_median_reduction_fraction', 1-median(b2)/max(median(b1),eps), ...
        'C2_improves_all_gap_pairs', all(b2 < b1), ...
        'calibration_condition_C2_median', median(L.calibration_condition_number(ql2)));
end
S = struct2table(rows);
end

function make_figure(S, outDir)
fig = figure('Color', 'w', 'Position', [100 100 980 390]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(S.tilt_deg, S.low_rmse_C1_median_V, 'o-', 'LineWidth', 1.3); hold on;
plot(S.tilt_deg, S.low_rmse_C2_median_V, 's-', 'LineWidth', 1.3);
xlabel('FE tilt angle (deg)'); ylabel('Median low-speed RMSE (V)');
legend('C1: g_0 only','C2: g_0, \mu','Location','northwest'); grid on; box on;
nexttile;
plot(S.tilt_deg, S.transfer_relative_C1_median, 'o-', 'LineWidth', 1.3); hold on;
plot(S.tilt_deg, S.transfer_relative_C2_median, 's-', 'LineWidth', 1.3);
xlabel('FE tilt angle (deg)'); ylabel('Median relative increment bias');
legend('C1: g_0 only','C2: g_0, \mu','Location','northwest'); grid on; box on;
exportgraphics(fig, fullfile(outDir, 'static_calibration_transfer.png'), 'Resolution', 220);
close(fig);
end

function write_report(S, Sources, outDir)
fid = fopen(fullfile(outDir, 'README.md'), 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# Static tilt calibration-transfer audit\n\n');
fprintf(fid, 'The tilted COMSOL surfaces generate independent low- and high-clearance truth. ');
fprintf(fid, 'C1 and C2 use the frozen zero-tilt library. No DAQ noise or dynamic solver is included.\n\n');
fprintf(fid, '## Data sources\n\n');
for i=1:height(Sources)
    fprintf(fid, '- %.1f deg: `%s` (%s)\n', Sources.tilt_deg(i), ...
        Sources.file_path(i), Sources.source_kind(i));
end
fprintf(fid, '\n## Gate result\n\n');
allImprove = all(S.C2_improves_all_gap_pairs);
medianReduction = median(S.C2_median_reduction_fraction);
fprintf(fid, '- C2 improves every adjacent gap pair at every angle: %d\n', allImprove);
fprintf(fid, '- Median reduction of relative transfer bias across angles: %.6g\n', medianReduction);
fprintf(fid, '- This audit tests baseline calibration versus finite clearance-increment transfer; it does not test EO recovery.\n');
end

function [filePath, sourceKind] = resolve_tilt_file(rootDir, angle, gapFE)
name = sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt', angle);
candidates = {fullfile(rootDir, 'data', '间隙的影响', name), ...
    fullfile(fileparts(mfilename('fullpath')), 'input_fe_tilt', name), ...
    fullfile('F:\Program Files\comsol_model\电容数据', name)};
for i=1:numel(candidates)
    if exist(candidates{i}, 'file') ~= 2, continue; end
    try
        [g,x,y] = load_stacked_gap_curves(candidates{i}, gapFE); %#ok<ASGLU>
        if numel(g)==numel(gapFE) && all(cellfun(@numel,x)>1) && all(cellfun(@numel,y)>1)
            filePath = candidates{i};
            if i==1, sourceKind = "workspace";
            elseif i==2, sourceKind = "audit_local_copy";
            else, sourceKind = "legacy_original"; end
            return;
        end
    catch
    end
end
error('tilt:MissingReadableFEFile', 'No readable tilted FE file found for %.1f deg.', angle);
end

function r = low_row()
r = struct('tilt_deg',NaN,'g_low_mm',NaN,'method',"",'g0_hat_mm',NaN, ...
    'mu_hat',NaN,'tau_fixed_mm',NaN,'zeta_fixed',NaN,'low_rmse_V',NaN, ...
    'calibration_condition_number',NaN,'g0_mu_column_correlation',NaN, ...
    'translation_residual_correlation',NaN);
end

function r = transfer_row()
r = struct('tilt_deg',NaN,'g_low_mm',NaN,'g_high_mm',NaN,'delta_gap_mm',NaN, ...
    'method',"",'g0_hat_mm',NaN,'mu_hat',NaN,'transfer_rmse_V',NaN, ...
    'transfer_relative_norm',NaN,'gap_direction_projection',NaN, ...
    'vibration_direction_projection',NaN,'true_increment_rms_V',NaN, ...
    'pred_increment_rms_V',NaN);
end

function r = summary_row()
r = struct('tilt_deg',NaN,'mu_hat_median',NaN,'tilt_severity_mm',NaN, ...
    'low_rmse_C1_median_V',NaN,'low_rmse_C2_median_V',NaN, ...
    'transfer_relative_C1_median',NaN,'transfer_relative_C2_median',NaN, ...
    'transfer_relative_C1_max',NaN,'transfer_relative_C2_max',NaN, ...
    'C2_median_reduction_fraction',NaN,'C2_improves_all_gap_pairs',false, ...
    'calibration_condition_C2_median',NaN);
end

function r = source_row()
r = struct('tilt_deg',NaN,'source_kind',"",'file_path',"");
end
