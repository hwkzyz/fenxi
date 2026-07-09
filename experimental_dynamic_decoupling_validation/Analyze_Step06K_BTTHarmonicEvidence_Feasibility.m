%% Analyze_Step06K_BTTHarmonicEvidence_Feasibility.m
% Diagnostic only: test whether Step06K can use BTT as harmonic evidence
% rather than as raw pointwise K calibration.

clear; clc;

root = fileparts(mfilename('fullpath'));
cases = {
    '20241106', '3000_3150';
    '20250527', '20250526_2500-3500_t400';
    '20251222', '1000_2500_3500';
    };

rows = repmat(empty_summary_row_local(), size(cases, 1), 1);
detail = table();

for iCase = 1:size(cases, 1)
    dataset = cases{iCase, 1};
    caseName = cases{iCase, 2};
    baseDir = fullfile(root, [dataset '_btt_data_foundation'], ...
        'output', 'step06k_transfer_ratio_robustness', caseName);
    files.points = fullfile(baseDir, ...
        sprintf('Step06K_TransferRatioRobustness_BTTCorrectionPoints_%s.csv', dataset));
    files.region = fullfile(baseDir, ...
        sprintf('Step06K_TransferRatioRobustness_RegionBand_%s.csv', dataset));
    files.timeseries = fullfile(baseDir, ...
        sprintf('Step06K_TransferRatioRobustness_TimeSeriesBand_%s.csv', dataset));
    files.pointfit = fullfile(baseDir, ...
        sprintf('Step06K_TransferRatioRobustness_KConstrainedFit_%s.csv', dataset));
    files.step05 = fullfile(baseDir, ...
        sprintf('Step06K_TransferRatioRobustness_Step05BandCheck_%s.csv', dataset));

    assert_file_local(files.points);
    assert_file_local(files.region);
    assert_file_local(files.timeseries);
    assert_file_local(files.pointfit);
    assert_file_local(files.step05);

    P = readtable(files.points);
    R = readtable(files.region);
    TS = readtable(files.timeseries);
    PF = readtable(files.pointfit);
    S5 = readtable(files.step05);

    D = analyze_case_local(dataset, caseName, P, R, TS, PF, S5);
    rows(iCase) = D.summary;
    detail = [detail; D.detail]; %#ok<AGROW>
end

Summary = struct2table(rows);
outCsv = fullfile(root, 'Step06K_BTTHarmonicEvidence_Feasibility.csv');
detailCsv = fullfile(root, 'Step06K_BTTHarmonicEvidence_FitDetail.csv');
writetable(Summary, outCsv);
writetable(detail, detailCsv);

disp(Summary);
fprintf('\nSaved:\n  %s\n  %s\n', outCsv, detailCsv);


function D = analyze_case_local(dataset, caseName, P, R, TS, PF, S5)
regionIds = unique(P.region_id(:)).';
detailRows = repmat(empty_detail_row_local(), numel(regionIds) * 2, 1);
row = 0;

for rid = regionIds
    Rp = R(R.region_id == rid, :);
    Tp = P(P.region_id == rid, :);
    TSp = TS(TS.region_id == rid, :);
    if ~isempty(Rp) && all(ismember({'core_time_start_s','core_time_end_s'}, ...
            Rp.Properties.VariableNames)) && isfinite(Rp.core_time_start_s(1))
        mCore = TSp.time_s >= Rp.core_time_start_s(1) & ...
            TSp.time_s <= Rp.core_time_end_s(1);
    else
        mCore = true(height(TSp), 1);
    end

    for modelName = ["shared_offset", "sensor_offset"]
        row = row + 1;
        Fit = fit_harmonic_local(Tp, modelName);
        uCore = Fit.c_re .* TSp.strain_basis_real_microstrain(mCore) + ...
            Fit.c_im .* TSp.strain_basis_imag_microstrain(mCore);
        ACore = sqrt(2) * std(uCore, 'omitnan');
        [ciLow, ciHigh] = bootstrap_A_ci_local(Tp, TSp(mCore, :), modelName, 300);

        detailRows(row).dataset = string(dataset);
        detailRows(row).case_name = string(caseName);
        detailRows(row).region_id = rid;
        detailRows(row).model = modelName;
        detailRows(row).dominant_order = first_or_nan_local(Tp.dominant_order);
        detailRows(row).point_count = height(Tp);
        detailRows(row).sensor_count = numel(unique(Tp.sensor_id));
        detailRows(row).phase_coverage_deg = phase_coverage_deg_local(Tp);
        detailRows(row).K_BTT_um_per_microstrain = 1000 * hypot(Fit.c_re, Fit.c_im);
        detailRows(row).A_BTT_mm = ACore;
        detailRows(row).A_BTT_boot_ci_low_mm = ciLow;
        detailRows(row).A_BTT_boot_ci_high_mm = ciHigh;
        detailRows(row).phase_BTT_rad = atan2(Fit.c_im, Fit.c_re);
        detailRows(row).rmse_mm = Fit.rmse;
        detailRows(row).cond_X = Fit.condX;
        detailRows(row).rank_X = Fit.rankX;
        detailRows(row).sensor_offset_span_mm = Fit.sensorOffsetSpan;
        detailRows(row).status = status_from_fit_local(Fit, Tp);
    end
end
detailRows = detailRows(1:row);
Detail = struct2table(detailRows);

modelA = Detail(Detail.model == "shared_offset", :);
modelB = Detail(Detail.model == "sensor_offset", :);
firstR = R(1, :);
firstPF = PF(1, :);

summary = empty_summary_row_local();
summary.dataset = string(dataset);
summary.case_name = string(caseName);
summary.dominant_order = first_or_nan_local(R.dominant_order);
summary.point_count = sum(modelA.point_count, 'omitnan');
summary.sensor_count = max(modelA.sensor_count, [], 'omitnan');
summary.phase_coverage_deg = median(modelA.phase_coverage_deg, 'omitnan');
summary.FE_A_nominal_mm = firstR.u_frequency_amp_nominal_mm(1);
summary.FE_A_low_mm = firstR.u_frequency_amp_band_low_mm(1);
summary.FE_A_high_mm = firstR.u_frequency_amp_band_high_mm(1);
summary.Step05_A_median_mm = median(S5.Step05_A_mm, 'omitnan');
summary.raw_pointwise_K_um_per_microstrain = firstPF.K_um_per_microstrain(1);
summary.raw_pointwise_A_mm = firstPF.corrected_A_peak_equiv_mm(1);
summary.raw_pointwise_status = string(firstPF.constraint_status(1));
summary.ModelA_A_BTT_mm = median(modelA.A_BTT_mm, 'omitnan');
summary.ModelB_A_BTT_mm = median(modelB.A_BTT_mm, 'omitnan');
summary.ModelA_to_FE_ratio = summary.ModelA_A_BTT_mm / summary.FE_A_nominal_mm;
summary.ModelB_to_FE_ratio = summary.ModelB_A_BTT_mm / summary.FE_A_nominal_mm;
summary.ModelA_to_Step05_ratio = summary.ModelA_A_BTT_mm / summary.Step05_A_median_mm;
summary.ModelB_to_Step05_ratio = summary.ModelB_A_BTT_mm / summary.Step05_A_median_mm;
summary.ModelAB_A_diff_percent = 100 * ...
    (summary.ModelB_A_BTT_mm - summary.ModelA_A_BTT_mm) / summary.ModelA_A_BTT_mm;
summary.ModelA_cond_X = median(modelA.cond_X, 'omitnan');
summary.ModelB_cond_X = median(modelB.cond_X, 'omitnan');
summary.ModelB_sensor_offset_span_mm = median(modelB.sensor_offset_span_mm, 'omitnan');
summary.ModelA_boot_ci_width_mm = median( ...
    modelA.A_BTT_boot_ci_high_mm - modelA.A_BTT_boot_ci_low_mm, 'omitnan');
summary.ModelB_boot_ci_width_mm = median( ...
    modelB.A_BTT_boot_ci_high_mm - modelB.A_BTT_boot_ci_low_mm, 'omitnan');
summary.feasibility_status = feasibility_status_local(summary);

D.summary = summary;
D.detail = Detail;
end


function Fit = fit_harmonic_local(T, modelName)
y = T.displacement_mm(:);
hRe = T.strain_basis_real_microstrain(:);
hIm = T.strain_basis_imag_microstrain(:);
sensorId = T.sensor_id(:);
ok = isfinite(y) & isfinite(hRe) & isfinite(hIm) & isfinite(sensorId);
y = y(ok); hRe = hRe(ok); hIm = hIm(ok); sensorId = sensorId(ok);

if modelName == "shared_offset"
    X = [ones(size(y)), hRe, hIm];
elseif modelName == "sensor_offset"
    sensors = unique(sensorId(:)).';
    S = zeros(numel(y), max(0, numel(sensors) - 1));
    for i = 2:numel(sensors)
        S(:, i - 1) = sensorId == sensors(i);
    end
    X = [ones(size(y)), S, hRe, hIm];
else
    error('Unknown model: %s', modelName);
end

beta = X \ y;
yFit = X * beta;
resid = y - yFit;
idxH = size(X, 2) - 1;
Fit.c_re = beta(idxH);
Fit.c_im = beta(idxH + 1);
Fit.rmse = sqrt(mean(resid .^ 2, 'omitnan'));
Fit.rankX = rank(X);
Fit.condX = normalized_cond_local(X);
Fit.sensorOffsetSpan = sensor_offset_span_local(beta, sensorId, modelName);
end


function c = normalized_cond_local(X)
Xn = X;
for j = 2:size(Xn, 2)
    s = std(Xn(:, j), 'omitnan');
    if isfinite(s) && s > 0
        Xn(:, j) = (Xn(:, j) - mean(Xn(:, j), 'omitnan')) ./ s;
    end
end
sval = svd(Xn, 'econ');
sval = sval(sval > max(size(Xn)) * eps(max(sval)));
if numel(sval) < 2
    c = Inf;
else
    c = max(sval) / min(sval);
end
end


function span = sensor_offset_span_local(beta, sensorId, modelName)
sensors = unique(sensorId(:)).';
if modelName ~= "sensor_offset" || numel(sensors) < 2
    span = 0;
    return;
end
offs = zeros(numel(sensors), 1);
offs(1) = beta(1);
for i = 2:numel(sensors)
    offs(i) = beta(1) + beta(i);
end
span = max(offs) - min(offs);
end


function covDeg = phase_coverage_deg_local(T)
ph = mod(atan2(T.strain_basis_imag_microstrain(:), ...
    T.strain_basis_real_microstrain(:)), 2*pi);
ph = sort(ph(isfinite(ph)));
if numel(ph) < 2
    covDeg = 0;
    return;
end
gaps = diff([ph; ph(1) + 2*pi]);
covDeg = (2*pi - max(gaps)) * 180 / pi;
end


function [ciLow, ciHigh] = bootstrap_A_ci_local(T, TSc, modelName, nBoot)
if height(T) < 10 || isempty(TSc)
    ciLow = NaN; ciHigh = NaN;
    return;
end
A = nan(nBoot, 1);
n = height(T);
for b = 1:nBoot
    idx = randi(n, n, 1);
    Fb = fit_harmonic_local(T(idx, :), modelName);
    u = Fb.c_re .* TSc.strain_basis_real_microstrain + ...
        Fb.c_im .* TSc.strain_basis_imag_microstrain;
    A(b) = sqrt(2) * std(u, 'omitnan');
end
ciLow = percentile_local(A, 2.5);
ciHigh = percentile_local(A, 97.5);
end


function q = percentile_local(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
    q = NaN;
    return;
end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    q = x(lo);
else
    q = x(lo) + (x(hi) - x(lo)) * (pos - lo);
end
end


function s = status_from_fit_local(Fit, T)
if Fit.rankX < 3
    s = "rank_deficient";
elseif Fit.condX > 100
    s = "ill_conditioned";
elseif phase_coverage_deg_local(T) < 120
    s = "poor_phase_coverage";
else
    s = "usable_diagnostic";
end
end


function s = feasibility_status_local(row)
if row.ModelA_cond_X > 100 || row.phase_coverage_deg < 120
    s = "not_identifiable";
elseif abs(row.ModelAB_A_diff_percent) > 25 || ...
        row.ModelB_sensor_offset_span_mm > 0.05
    s = "diagnostic_only_sensor_layered";
elseif row.ModelA_boot_ci_width_mm > 0.05
    s = "diagnostic_only_wide_ci";
else
    s = "usable_harmonic_evidence";
end
end


function row = empty_summary_row_local()
row = struct('dataset', "", 'case_name', "", 'dominant_order', NaN, ...
    'point_count', NaN, 'sensor_count', NaN, 'phase_coverage_deg', NaN, ...
    'FE_A_nominal_mm', NaN, 'FE_A_low_mm', NaN, 'FE_A_high_mm', NaN, ...
    'Step05_A_median_mm', NaN, ...
    'raw_pointwise_K_um_per_microstrain', NaN, ...
    'raw_pointwise_A_mm', NaN, 'raw_pointwise_status', "", ...
    'ModelA_A_BTT_mm', NaN, 'ModelB_A_BTT_mm', NaN, ...
    'ModelA_to_FE_ratio', NaN, 'ModelB_to_FE_ratio', NaN, ...
    'ModelA_to_Step05_ratio', NaN, 'ModelB_to_Step05_ratio', NaN, ...
    'ModelAB_A_diff_percent', NaN, ...
    'ModelA_cond_X', NaN, 'ModelB_cond_X', NaN, ...
    'ModelB_sensor_offset_span_mm', NaN, ...
    'ModelA_boot_ci_width_mm', NaN, ...
    'ModelB_boot_ci_width_mm', NaN, ...
    'feasibility_status', "");
end


function row = empty_detail_row_local()
row = struct('dataset', "", 'case_name', "", 'region_id', NaN, ...
    'model', "", 'dominant_order', NaN, 'point_count', NaN, ...
    'sensor_count', NaN, 'phase_coverage_deg', NaN, ...
    'K_BTT_um_per_microstrain', NaN, 'A_BTT_mm', NaN, ...
    'A_BTT_boot_ci_low_mm', NaN, 'A_BTT_boot_ci_high_mm', NaN, ...
    'phase_BTT_rad', NaN, 'rmse_mm', NaN, 'cond_X', NaN, ...
    'rank_X', NaN, 'sensor_offset_span_mm', NaN, 'status', "");
end


function v = first_or_nan_local(x)
if isempty(x)
    v = NaN;
else
    v = x(1);
end
end


function assert_file_local(file)
if ~isfile(file)
    error('Missing file:\n  %s', file);
end
end
