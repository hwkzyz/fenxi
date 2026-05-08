%% Analyze gap-waveform mapping structure
% Systematic investigation of 5 research questions:
%   Q1: At fixed x, is F(x,g) more linear in g or in 1/g?
%   Q2: What coordinate system makes gap-to-gap mapping simplest?
%   Q3: Is the static waveform family low-rank?
%   Q4: Do geometric features follow simple functional forms in g?
%   Q5: Type A (amplitude-only) vs Type B (+width) vs Type C (+nonlinear warp)?
%
% Input: 直叶片2mm_不同间隙.txt — 7 static waveforms, g = 0.2:0.2:1.4 mm

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'gap_waveform_mapping_structure_results');
if ~isfolder(outDir), mkdir(outDir); end

%% Figure defaults
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 8.5, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.15, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

%% ========================================================================
% LOAD AND PREPROCESS
% ========================================================================
[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);  % largest gap as reference
gRef = gapList(refIdx);

% Preprocess: baseline correction, peak centering
curve = struct();
for ig = 1:nGap
    x = xCell{ig}(:);
    y = yCell{ig}(:);
    base = min(y);
    yBc = y - base;
    amp = max(yBc);
    yNorm = yBc / amp;
    [~, iPeak] = max(yBc);
    xPeak = x(iPeak);
    curve(ig).gap = gapList(ig);
    curve(ig).x = x;
    curve(ig).yRaw = y;
    curve(ig).yBc = yBc;
    curve(ig).base = base;
    curve(ig).amp = amp;
    curve(ig).xPeak = xPeak;
    curve(ig).xCenter = x - xPeak;
    curve(ig).yNorm = yNorm;
end

% Interpolate all to common fine raw-x grid for SVD and fixed-x analysis
xAll = cell2mat(xCell);
xCommon = linspace(min(xAll), max(xAll), 2000)';
YrawCommon = nan(nGap, numel(xCommon));
for ig = 1:nGap
    YrawCommon(ig, :) = interp1(curve(ig).x, curve(ig).yRaw, xCommon, 'pchip');
end

% Interpolate to common peak-centered grid
xcAll = cell2mat({curve.xCenter}');
xcCommon = linspace(min(xcAll), max(xcAll), 2000)';
YbcCommon = nan(nGap, numel(xcCommon));
for ig = 1:nGap
    YbcCommon(ig, :) = interp1(curve(ig).xCenter, curve(ig).yBc, xcCommon, 'pchip', 'extrap');
end

% FWHM for each waveform
FWHM = nan(nGap, 1);
for ig = 1:nGap
    halfMax = curve(ig).amp / 2;
    halfMaxAbs = curve(ig).base + halfMax;
    xRise = curve(ig).x(1:find(curve(ig).yBc == max(curve(ig).yBc), 1));
    yRise = curve(ig).yRaw(1:find(curve(ig).yBc == max(curve(ig).yBc), 1));
    xFall = curve(ig).x(find(curve(ig).yBc == max(curve(ig).yBc), 1):end);
    yFall = curve(ig).yRaw(find(curve(ig).yBc == max(curve(ig).yBc), 1):end);
    xL = interp1(yRise, xRise, halfMaxAbs, 'linear');
    xR = interp1(flipud(yFall), flipud(xFall), halfMaxAbs, 'linear');
    FWHM(ig) = xR - xL;
end

% Equal-level registration to reference
levelMap = linspace(0.005, 0.995, 199)';
nU = 1401;
[curve, uGrid, Cmat, Hmat, Ameas] = ...
    build_equal_level_registration(curve, gapList, refIdx, levelMap, nU);
Cref = Cmat(refIdx, :).';
Href = Hmat(refIdx, :).';
Aref = Ameas(refIdx);

%% ========================================================================
% Q1: FIXED-x POINT-VALUE ANALYSIS
%   For several x positions, examine F(x,g) vs g and vs 1/g.
%   Fit candidate models and compare R².
% ========================================================================
fprintf('\n========== Q1: Fixed-x point-value analysis ==========\n');

% Select representative x positions (in peak-centered coordinate)
% covering left wing, near-peak, right wing
xcSelected = [prctile(xcCommon, 10), prctile(xcCommon, 30), ...
              0, prctile(xcCommon, 70), prctile(xcCommon, 90)];
xcSelected = sort(xcSelected);
nXc = numel(xcSelected);

% Extract F(x,g) at these positions
F_at_xc = nan(nGap, nXc);
for ix = 1:nXc
    for ig = 1:nGap
        F_at_xc(ig, ix) = interp1(curve(ig).xCenter, curve(ig).yBc, xcSelected(ix), 'pchip');
    end
end

% Fit candidate models for each x position
modelNames = {'Linear in g: a + b*g', ...
              'Linear in 1/g: a + b/g', ...
              'Shifted inv: a + b/(g+delta)', ...
              'Power law: a * g^b', ...
              'Exponential: a * exp(b*g)'};
nModels = numel(modelNames);

q1Results = table();
for ix = 1:nXc
    yVec = F_at_xc(:, ix);
    gVec = gapList;
    valid = isfinite(yVec);
    yVal = yVec(valid);
    gVal = gVec(valid);

    R2 = nan(1, nModels);
    params = cell(1, nModels);
    yPred = nan(nGap, nModels);

    % Model 1: Linear in g
    X1 = [ones(sum(valid), 1), gVal];
    b1 = X1 \ yVal;
    yPred(valid, 1) = X1 * b1;
    R2(1) = 1 - sum((yVal - yPred(valid, 1)).^2) / sum((yVal - mean(yVal)).^2);
    params{1} = b1;

    % Model 2: Linear in 1/g
    X2 = [ones(sum(valid), 1), 1./gVal];
    b2 = X2 \ yVal;
    yPred(valid, 2) = X2 * b2;
    R2(2) = 1 - sum((yVal - yPred(valid, 2)).^2) / sum((yVal - mean(yVal)).^2);
    params{2} = b2;

    % Model 3: Shifted inverse a + b/(g+delta)
    % Use fminsearch for delta, then linear for a,b
    obj3 = @(delta) fit_shifted_inv(yVal, gVal, delta);
    [deltaOpt, resOpt] = fminsearch(obj3, 0.5);
    X3 = [ones(sum(valid), 1), 1./(gVal + deltaOpt)];
    b3 = X3 \ yVal;
    yPred(valid, 3) = X3 * b3;
    R2(3) = 1 - sum((yVal - yPred(valid, 3)).^2) / sum((yVal - mean(yVal)).^2);
    params{3} = [b3; deltaOpt];

    % Model 4: Power law a * g^b
    % log(y) = log(a) + b*log(g) — only valid if y > 0
    if all(yVal > 0)
        X4 = [ones(sum(valid), 1), log(gVal)];
        b4 = X4 \ log(yVal);
        yPred(valid, 4) = exp(b4(1)) * gVal.^b4(2);
        R2(4) = 1 - sum((yVal - yPred(valid, 4)).^2) / sum((yVal - mean(yVal)).^2);
        params{4} = [exp(b4(1)); b4(2)];
    else
        R2(4) = NaN;
        params{4} = [NaN; NaN];
    end

    % Model 5: Exponential a*exp(b*g)
    if all(yVal > 0)
        X5 = [ones(sum(valid), 1), gVal];
        b5 = X5 \ log(yVal);
        yPred(valid, 5) = exp(b5(1) + b5(2)*gVal);
        R2(5) = 1 - sum((yVal - yPred(valid, 5)).^2) / sum((yVal - mean(yVal)).^2);
        params{5} = [exp(b5(1)); b5(2)];
    else
        R2(5) = NaN;
        params{5} = [NaN; NaN];
    end

    row = table(xcSelected(ix), R2(1), R2(2), R2(3), R2(4), R2(5), ...
        'VariableNames', {'xc_mm', 'R2_linear_g', 'R2_linear_1overg', ...
        'R2_shifted_inv', 'R2_power', 'R2_exponential'});
    q1Results = [q1Results; row]; %#ok<AGROW>
end

fprintf('Fixed-x model fit R² values:\n');
disp(q1Results);
fprintf('Best model at each x: %s\n', ...
    strjoin(modelNames(q1BestIdx(q1Results)), ', '));

writetable(q1Results, fullfile(outDir, 'q1_fixed_x_point_analysis.csv'));

%% ========================================================================
% Q2: COORDINATE CORRESPONDENCE COMPARISON
%   Compare 4 coordinate systems for shape collapse:
%   (a) Raw x (COMSOL coordinate)
%   (b) Peak-centered x - x_peak
%   (c) Width-normalized (x - x_peak) / FWHM(g)
%   (d) Equal-level registered u
% ========================================================================
fprintf('\n========== Q2: Coordinate correspondence comparison ==========\n');

% (a) Raw x — interpolate to common raw grid
YrawNorm = YrawCommon - min(YrawCommon, [], 2);
YrawNorm = YrawNorm ./ max(YrawNorm, [], 2);

% (b) Peak-centered x — already computed as YbcCommon
YbcNorm = YbcCommon ./ max(YbcCommon, [], 2);

% (c) Width-normalized coordinate
nWn = 2000;
wnCommon = linspace(-3, 3, nWn)';  % in FWHM units
YwnCommon = nan(nGap, nWn);
for ig = 1:nGap
    xScaled = curve(ig).xCenter / FWHM(ig);
    YwnCommon(ig, :) = interp1(xScaled, curve(ig).yNorm, wnCommon, 'linear', 'extrap');
YwnCommon(ig, :) = max(0, min(1, YwnCommon(ig, :)));  % clip extrapolation
end

% (d) Equal-level registered — already computed as Hmat
% (normalized to unit peak)

% Compute shape collapse metrics for each coordinate system
coordNames = {'Raw x', 'Peak-centered x', 'Width-normalized', 'Equal-level registered'};
coordTags = {'raw_x', 'peak_centered', 'width_norm', 'equal_level'};
Ysets = {YrawNorm, YbcNorm, YwnCommon, Hmat};
xGrids = {xCommon, xcCommon, wnCommon, uGrid};

q2ShapeRmse = nan(nGap, numel(coordNames));
q2ShapeCorr = nan(nGap, numel(coordNames));

for ic = 1:numel(coordNames)
    Y = Ysets{ic};
    Yref = Y(refIdx, :);
    for ig = 1:nGap
        valid = isfinite(Y(ig, :)) & isfinite(Yref);
        if sum(valid) > 20
            q2ShapeRmse(ig, ic) = sqrt(mean((Y(ig, valid) - Yref(valid)).^2, 'omitnan'));
            ct = corrcoef(Y(ig, valid)', Yref(valid)');
            if numel(ct) >= 4
                q2ShapeCorr(ig, ic) = ct(1, 2);
            end
        end
    end
end
% Replace any remaining NaN/Inf with median of column
for ic = 1:numel(coordNames)
    col = q2ShapeRmse(:, ic);
    col(~isfinite(col)) = median(col(isfinite(col)));
    q2ShapeRmse(:, ic) = col;
    col = q2ShapeCorr(:, ic);
    col(~isfinite(col)) = median(col(isfinite(col)));
    q2ShapeCorr(:, ic) = col;
end

% Summary table
q2Summary = table(coordNames', coordTags', ...
    mean(q2ShapeRmse, 1)', max(q2ShapeRmse, [], 1)', ...
    mean(q2ShapeCorr, 1)', min(q2ShapeCorr, [], 1)', ...
    'VariableNames', {'coord_system', 'tag', 'mean_shape_rmse', ...
    'max_shape_rmse', 'mean_shape_corr', 'min_shape_corr'});
fprintf('Shape collapse quality by coordinate system:\n');
disp(q2Summary);

writetable(q2Summary, fullfile(outDir, 'q2_coordinate_comparison.csv'));

% Gap-by-gap detail
q2Detail = table(gapList, q2ShapeRmse(:,1), q2ShapeRmse(:,2), q2ShapeRmse(:,3), q2ShapeRmse(:,4), ...
    q2ShapeCorr(:,1), q2ShapeCorr(:,2), q2ShapeCorr(:,3), q2ShapeCorr(:,4), ...
    'VariableNames', {'gap_mm', ...
    'rmse_raw_x', 'rmse_peak_centered', 'rmse_width_norm', 'rmse_equal_level', ...
    'corr_raw_x', 'corr_peak_centered', 'corr_width_norm', 'corr_equal_level'});
writetable(q2Detail, fullfile(outDir, 'q2_coordinate_comparison_by_gap.csv'));

%% ========================================================================
% Q3: LOW-RANK ANALYSIS
%   SVD on: (a) raw waveforms on common x-grid
%           (b) registered waveforms C(u,g)
%           (c) normalized shapes H(u,g)
% ========================================================================
fprintf('\n========== Q3: Low-rank analysis ==========\n');

% (a) Raw baseline-corrected on peak-centered grid
[U_raw, S_raw, V_raw] = svd(YbcCommon, 'econ');
sRaw = diag(S_raw);
varRaw = sRaw.^2 / sum(sRaw.^2);

% (b) Registered waveforms C(u,g)
[U_reg, S_reg, V_reg] = svd(Cmat, 'econ');
sReg = diag(S_reg);
varReg = sReg.^2 / sum(sReg.^2);

% (c) Normalized shapes H(u,g)
[U_shp, S_shp, V_shp] = svd(Hmat, 'econ');
sShp = diag(S_shp);
varShp = sShp.^2 / sum(sShp.^2);

% Rank-k reconstruction error
maxRank = min(nGap, 5);
recErrRaw = nan(maxRank, 1);
recErrReg = nan(maxRank, 1);
recErrShp = nan(maxRank, 1);
for k = 1:maxRank
    recErrRaw(k) = norm(YbcCommon - U_raw(:,1:k) * S_raw(1:k,1:k) * V_raw(:,1:k)', 'fro') ...
        / norm(YbcCommon, 'fro');
    recErrReg(k) = norm(Cmat - U_reg(:,1:k) * S_reg(1:k,1:k) * V_reg(:,1:k)', 'fro') ...
        / norm(Cmat, 'fro');
    recErrShp(k) = norm(Hmat - U_shp(:,1:k) * S_shp(1:k,1:k) * V_shp(:,1:k)', 'fro') ...
        / norm(Hmat, 'fro');
end

q3Summary = table((1:maxRank)', varRaw(1:maxRank), cumsum(varRaw(1:maxRank)), ...
    varReg(1:maxRank), cumsum(varReg(1:maxRank)), ...
    varShp(1:maxRank), cumsum(varShp(1:maxRank)), ...
    recErrRaw, recErrReg, recErrShp, ...
    'VariableNames', {'rank', 'var_raw', 'cumvar_raw', ...
    'var_registered', 'cumvar_registered', 'var_shape', 'cumvar_shape', ...
    'rel_err_raw', 'rel_err_registered', 'rel_err_shape'});
fprintf('Low-rank analysis (variance explained):\n');
disp(q3Summary);

writetable(q3Summary, fullfile(outDir, 'q3_low_rank_svd.csv'));

%% ========================================================================
% Q4: GEOMETRIC FEATURE TRENDS
%   P(g), S(g), W_{1/2}(g), w(q,g) — fit functional forms
% ========================================================================
fprintf('\n========== Q4: Geometric feature trends ==========\n');

% Extract features
Peak = nan(nGap, 1);
Area = nan(nGap, 1);
for ig = 1:nGap
    Peak(ig) = curve(ig).amp;
    Area(ig) = trapz(curve(ig).xCenter, curve(ig).yBc);
end

% Multi-level widths at selected levels
qLevels = [0.1, 0.25, 0.5, 0.75, 0.9]';
nQ = numel(qLevels);
Widths = nan(nGap, nQ);
for ig = 1:nGap
    xc = curve(ig).xCenter;
    yNorm = curve(ig).yNorm;
    [~, iPeak] = max(yNorm);
    for iq = 1:nQ
        xL = interp1(yNorm(1:iPeak), xc(1:iPeak), qLevels(iq), 'linear');
        xR = interp1(flipud(yNorm(iPeak:end)), flipud(xc(iPeak:end)), qLevels(iq), 'linear');
        if isfinite(xL) && isfinite(xR)
            Widths(ig, iq) = xR - xL;
        end
    end
end

% Fit candidate models for each feature
featureNames = {'Peak P(g)', 'Area S(g)', 'FWHM(g)'};
featureTags = {'peak', 'area', 'fwhm'};
featureVals = [Peak, Area, FWHM];

% For widths, also create width-level entries
for iq = 1:nQ
    featureNames{end+1} = sprintf('Width w(q=%.2f,g)', qLevels(iq));
    featureTags{end+1} = sprintf('width_q%02d', round(qLevels(iq)*100));
    featureVals(:, end+1) = Widths(:, iq);
end

nFeat = numel(featureNames);
candidateModels = {'a + b*g', 'a + b/g', 'a + b/(g+delta)', 'a * g^b', 'a * exp(b*g)'};
modelTags = {'linear_g', 'linear_1overg', 'shifted_inv', 'power', 'exponential'};
nModel = numel(candidateModels);

q4Results = table();
for ife = 1:nFeat
    yVec = featureVals(:, ife);
    valid = isfinite(yVec);
    yVal = yVec(valid);
    gVal = gapList(valid);
    R2 = nan(1, nModel);

    % Linear in g
    X = [ones(size(gVal)), gVal];
    b = X \ yVal;
    yhat = X * b;
    R2(1) = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);

    % Linear in 1/g
    X = [ones(size(gVal)), 1./gVal];
    b = X \ yVal;
    yhat = X * b;
    R2(2) = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);

    % Shifted inverse
    obj = @(delta) fit_shifted_inv(yVal, gVal, delta);
    [dOpt, ~] = fminsearch(obj, 0.3);
    X = [ones(size(gVal)), 1./(gVal + dOpt)];
    b = X \ yVal;
    yhat = X * b;
    R2(3) = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);

    % Power law
    if all(yVal > 0)
        X = [ones(size(gVal)), log(gVal)];
        b = X \ log(yVal);
        yhat = exp(b(1) + b(2)*log(gVal));
        R2(4) = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
    end

    % Exponential
    if all(yVal > 0)
        X = [ones(size(gVal)), gVal];
        b = X \ log(yVal);
        yhat = exp(b(1) + b(2)*gVal);
        R2(5) = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
    end

    row = cell2table([{featureNames{ife}, featureTags{ife}}, num2cell(R2)], ...
        'VariableNames', [{'feature', 'tag'}, strcat('R2_', modelTags)]);
    q4Results = [q4Results; row]; %#ok<AGROW>
end

fprintf('Geometric feature functional fit R²:\n');
disp(q4Results);

writetable(q4Results, fullfile(outDir, 'q4_geometric_feature_fits.csv'));

%% ========================================================================
% Q5: MODEL COMPARISON — Type A vs B vs C
%   Type A: F(x,g) ≈ B(g) + A(g)*H(x)            [amplitude only]
%   Type B: F(x,g) ≈ B(g) + A(g)*H((x-x_c)/w(g)) [amplitude + width]
%   Type C: F(x,g) ≈ B(g) + A(g)*H(ψ_g(x))       [full nonlinear warp]
%
% We compare by: how well does each coordinate system collapse normalized
% shapes H across gaps? This directly maps to:
%   Type A ↔ peak-centered x (shape RMSE in x_c domain)
%   Type B ↔ width-normalized x/FWHM (shape RMSE in scaled domain)
%   Type C ↔ equal-level registered u (shape RMSE in u domain)
%
% For reconstruction quality: each type uses reference shape + amplitude
% to predict the reference waveform from a non-reference gap's data.
% ========================================================================
fprintf('\n========== Q5: Model-type comparison ==========\n');

modelLabels = {'Type A: A*H(x) - amplitude only', ...
               'Type B: A*H(x/w) - amplitude + width', ...
               'Type C: A*H(u) - full registration'};

% Shape collapse metrics from Q2 (reinterpreted):
% Type A = peak-centered x (column 2 of q2ShapeRmse/q2ShapeCorr)
% Type B = width-normalized (column 3)
% Type C = equal-level registered (column 4)
q5ShapeRmse = q2ShapeRmse(:, [2, 3, 4]);
q5ShapeCorr = q2ShapeCorr(:, [2, 3, 4]);

% Additional: "reference-gap reconstruction" error for each type
% For each non-ref gap g, predict what reference-gap waveform would be.
% Type A: take g's shape in peak-centered x, give it ref amplitude
% Type B: take g's shape in width-scaled x, give it ref amplitude+width
% Type C: take g's shape in u, give it ref amplitude (i.e. CrefRecon)
q5ReconRmse = nan(nGap, 3);

% Pre-build Type A and B reconstruction grids
xcFine = linspace(min(xcAll), max(xcAll), 2000)';

for ig = 1:nGap
    % True reference waveform (in peak-centered x)
    yRefTrue = interp1(curve(refIdx).xCenter, curve(refIdx).yBc, xcFine, 'pchip');
    yRefTrue = yRefTrue(:);
    validMask = isfinite(yRefTrue);

    % --- Type A reconstruction ---
    % pred_A(x) = A_ref * (yBc(g_i, x) / A(g_i))
    yGi = interp1(curve(ig).xCenter, curve(ig).yBc, xcFine, 'pchip');
    yGi = yGi(:);
    AGi = max(yGi);
    yTypeA = (Aref / max(AGi, eps)) * yGi;
    va = validMask & isfinite(yTypeA);
    if sum(va) > 20
        q5ReconRmse(ig, 1) = sqrt(mean((yTypeA(va) - yRefTrue(va)).^2));
    end

    % --- Type B reconstruction ---
    % pred_B(x) = A_ref * H_gi(x * w_ref/w_gi, g_i)
    wRatio = FWHM(refIdx) / max(FWHM(ig), eps);
    xcWarped = xcFine * wRatio;
    yGiWarped = interp1(curve(ig).xCenter, curve(ig).yBc, xcWarped, 'linear', 'extrap');
    yGiWarped = yGiWarped(:);
    AGiW = max(abs(yGiWarped));
    yTypeB = (Aref / max(AGiW, eps)) * yGiWarped;
    vb = validMask & isfinite(yTypeB);
    if sum(vb) > 20
        q5ReconRmse(ig, 2) = sqrt(mean((yTypeB(vb) - yRefTrue(vb)).^2));
    end

    % --- Type C reconstruction ---
    % CrefRecon(ig,:) = A_ref * H(u,g_i) reconstructed in u-domain
    % Compare in u-domain against Cref
    recInU = (Aref / max(Ameas(ig), eps)) * Cmat(ig, :);
    CrefRow = Cref';  % 1 × nU
    vc = isfinite(recInU) & isfinite(CrefRow);
    if sum(vc) > 20
        q5ReconRmse(ig, 3) = sqrt(mean((recInU(vc) - CrefRow(vc)).^2));
    end
end

q5Summary = table(modelLabels', ...
    mean(q5ShapeRmse, 1)', max(q5ShapeRmse, [], 1)', ...
    mean(q5ShapeCorr, 1)', min(q5ShapeCorr, [], 1)', ...
    mean(q5ReconRmse, 1)', max(q5ReconRmse, [], 1)', ...
    'VariableNames', {'model_type', 'mean_shape_rmse', 'max_shape_rmse', ...
    'mean_shape_corr', 'min_shape_corr', ...
    'mean_recon_rmse_pF', 'max_recon_rmse_pF'});
fprintf('Model-type shape collapse and reconstruction quality:\n');
disp(q5Summary);

writetable(q5Summary, fullfile(outDir, 'q5_model_type_comparison.csv'));

q5Detail = table(gapList, ...
    q5ShapeRmse(:,1), q5ShapeRmse(:,2), q5ShapeRmse(:,3), ...
    q5ShapeCorr(:,1), q5ShapeCorr(:,2), q5ShapeCorr(:,3), ...
    q5ReconRmse(:,1), q5ReconRmse(:,2), q5ReconRmse(:,3), ...
    'VariableNames', {'gap_mm', ...
    'shape_rmse_A', 'shape_rmse_B', 'shape_rmse_C', ...
    'shape_corr_A', 'shape_corr_B', 'shape_corr_C', ...
    'recon_rmse_A', 'recon_rmse_B', 'recon_rmse_C'});
writetable(q5Detail, fullfile(outDir, 'q5_model_type_by_gap.csv'));

%% ========================================================================
% FIGURES
% ========================================================================

%% Figure 1: Q1 — Fixed-x point analysis
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl1 = tiledlayout(fig1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

colors = lines(max(nGap, 14));

% (a) F(x,g) vs x for all gaps, with markers at selected x
nexttile; hold on;
for ig = 1:nGap
    plot(curve(ig).xCenter, curve(ig).yBc, '-', 'Color', colors(ig, :), ...
        'DisplayName', sprintf('g=%.1f', gapList(ig)));
end
for ix = 1:nXc
    xline(xcSelected(ix), 'k--', 'LineWidth', 0.6);
    text(xcSelected(ix), max(YbcCommon(:))*0.95, sprintf('x_%d', ix), ...
        'FontSize', 6.5, 'VerticalAlignment', 'top');
end
xlabel('Peak-centered coordinate (mm)');
ylabel('Baseline-corrected C (pF)');
title('(a) Waveforms with analysis positions');
legend('Location', 'eastoutside', 'FontSize', 6.0);

% (b) F(x,g) vs g at selected positions
nexttile; hold on;
markers = {'o', 's', '^', 'd', 'v'};
for ix = 1:nXc
    plot(gapList, F_at_xc(:, ix), ['-' markers{ix}], 'Color', colors(ix*2, :), ...
        'MarkerFaceColor', colors(ix*2, :), 'MarkerSize', 5, ...
        'DisplayName', sprintf('x_c = %.2f', xcSelected(ix)));
end
xlabel('Gap g (mm)');
ylabel('F(x_c, g) (pF)');
title('(b) F(x,g) vs gap at fixed positions');
legend('Location', 'best', 'FontSize', 6.5);

% (c) F vs 1/g at selected positions
nexttile; hold on;
for ix = 1:nXc
    plot(1./gapList, F_at_xc(:, ix), ['-' markers{ix}], 'Color', colors(ix*2, :), ...
        'MarkerFaceColor', colors(ix*2, :), 'MarkerSize', 5, ...
        'DisplayName', sprintf('x_c = %.2f', xcSelected(ix)));
end
xlabel('1/g (mm^{-1})');
ylabel('F(x_c, g) (pF)');
title('(c) F(x,g) vs 1/g at fixed positions');
legend('Location', 'best', 'FontSize', 6.5);

% (d) Best-fit shifted-inverse model for representative x
nexttile; hold on;
repXc = xcSelected(3);  % near-peak
repY = F_at_xc(:, 3);
yVal = repY(isfinite(repY));
gVal = gapList(isfinite(repY));
[deltaOpt, ~] = fminsearch(@(d) fit_shifted_inv(yVal, gVal, d), 0.3);
Xfit = [ones(size(gVal)), 1./(gVal + deltaOpt)];
bFit = Xfit \ yVal;
gFine = linspace(min(gapList)*0.8, max(gapList)*1.2, 200)';
yFineFit = bFit(1) + bFit(2)./(gFine + deltaOpt);

plot(gapList, repY, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, ...
    'DisplayName', 'Data');
plot(gFine, yFineFit, 'r-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('a+b/(g+%.3f), R^2=%.4f', deltaOpt, ...
    1 - sum((yVal - (bFit(1) + bFit(2)./(gVal + deltaOpt))).^2) / sum((yVal - mean(yVal)).^2)));
xlabel('Gap g (mm)');
ylabel('F(x_c, g) (pF)');
title(sprintf('(d) Shifted-inverse fit at x_c = %.2f', repXc));
legend('Location', 'best', 'FontSize', 7);

% (e) R² comparison across models at each x
nexttile;
R2mat = q1Results{:, 2:end};
bar(categorical(string(round(q1Results.xc_mm, 2))), R2mat);
xlabel('Analysis position x_c (mm)');
ylabel('R^2');
title('(e) Model fit quality by position');
legend(modelNames, 'Location', 'eastoutside', 'FontSize', 6.0);
ylim([0.8, 1.005]);

% (f) Best model per position
nexttile; hold on;
[~, bestIdx] = max(R2mat, [], 2);
xPos = 1:nXc;
for ix = 1:nXc
    bar(xPos(ix), 1, 0.6, 'FaceColor', colors(bestIdx(ix)*2, :), ...
        'EdgeColor', 'k', 'LineWidth', 0.5);
end
set(gca, 'XTick', xPos, 'XTickLabel', string(round(q1Results.xc_mm, 2)));
xlabel('Analysis position x_c (mm)');
ylabel('Best model (categorical)');
title('(f) Best model at each x');
legend(modelNames, 'Location', 'eastoutside', 'FontSize', 6.5);

title(tl1, 'Q1: Fixed-x point-value analysis', 'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_fixed_x_analysis.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_fixed_x_analysis.pdf'), 'ContentType', 'vector');

%% Figure 2: Q2 — Coordinate comparison
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl2 = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

coordPlotNames = {'(a) Raw x', '(b) Peak-centered x', '(c) Width-normalized', '(d) Equal-level registered'};
Yplots = {YrawNorm, YbcNorm, YwnCommon, Hmat};
xPlots = {xCommon, xcCommon, wnCommon, uGrid};
xLabels = {'x (mm)', 'x - x_p (mm)', '(x - x_p) / FWHM', 'u (mm)'};

for ic = 1:4
    nexttile; hold on;
    for ig = 1:nGap
        plot(xPlots{ic}, Yplots{ic}(ig, :), '-', 'Color', colors(ig, :));
    end
    xlabel(xLabels{ic});
    ylabel('Normalized waveform');
    title(coordPlotNames{ic});
end

% (e) RMSE comparison
nexttile;
bar(categorical(coordNames), mean(q2ShapeRmse, 1)', 0.6, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean shape RMSE to reference');
title('(e) Shape collapse RMSE');
xtickangle(15);

% (f) Per-gap RMSE curves
nexttile; hold on;
for ic = 1:4
    plot(gapList, q2ShapeRmse(:, ic), '-o', 'Color', colors(ic*2, :), ...
        'MarkerFaceColor', colors(ic*2, :), 'MarkerSize', 4, ...
        'DisplayName', coordNames{ic});
end
xlabel('Gap g (mm)');
ylabel('Shape RMSE to reference');
title('(f) RMSE by gap and coordinate');
legend('Location', 'best', 'FontSize', 6.5);

title(tl2, 'Q2: Coordinate correspondence comparison', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_coordinate_comparison.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_coordinate_comparison.pdf'), 'ContentType', 'vector');

%% Figure 3: Q3 — Low-rank analysis
fig3 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl3 = tiledlayout(fig3, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Singular value spectrum
nexttile; hold on;
semilogy(1:nGap, sRaw/sRaw(1), 'o-', 'Color', colors(1,:), ...
    'MarkerFaceColor', colors(1,:), 'DisplayName', 'Raw peak-centered');
semilogy(1:nGap, sReg/sReg(1), 's-', 'Color', colors(3,:), ...
    'MarkerFaceColor', colors(3,:), 'DisplayName', 'Registered C(u,g)');
semilogy(1:nGap, sShp/sShp(1), '^-', 'Color', colors(5,:), ...
    'MarkerFaceColor', colors(5,:), 'DisplayName', 'Shape H(u,g)');
xlabel('Singular value index k');
ylabel('Normalized singular value \sigma_k / \sigma_1');
title('(a) Singular value spectrum');
legend('Location', 'northeast', 'FontSize', 7);
set(gca, 'XTick', 1:nGap);

% (b) Cumulative variance
nexttile; hold on;
plot(1:nGap, cumsum(varRaw), 'o-', 'Color', colors(1,:), ...
    'MarkerFaceColor', colors(1,:), 'DisplayName', 'Raw peak-centered');
plot(1:nGap, cumsum(varReg), 's-', 'Color', colors(3,:), ...
    'MarkerFaceColor', colors(3,:), 'DisplayName', 'Registered C(u,g)');
plot(1:nGap, cumsum(varShp), '^-', 'Color', colors(5,:), ...
    'MarkerFaceColor', colors(5,:), 'DisplayName', 'Shape H(u,g)');
yline(0.99, 'k--', 'LineWidth', 0.8);
xlabel('Rank k');
ylabel('Cumulative variance fraction');
title('(b) Cumulative explained variance');
legend('Location', 'southeast', 'FontSize', 7);
set(gca, 'XTick', 1:nGap);

% (c) Rank-1,2,3 reconstructions for raw data
nexttile; hold on;
plot(xcCommon, YbcCommon(refIdx, :), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Reference');
rec1 = U_raw(:,1) * S_raw(1,1) * V_raw(:,1)';
plot(xcCommon, rec1(refIdx, :), '--', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Rank-1 recon');
rec2 = U_raw(:,1:2) * S_raw(1:2,1:2) * V_raw(:,1:2)';
plot(xcCommon, rec2(refIdx, :), '--', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Rank-2 recon');
xlabel('x - x_p (mm)');
ylabel('C (pF)');
title('(c) Raw waveform SVD recon');
legend('Location', 'best', 'FontSize', 7);

% (d) Rank-1,2 reconstructions for registered data
nexttile; hold on;
plot(uGrid, Cref, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Reference');
rec1r = U_reg(:,1) * S_reg(1,1) * V_reg(:,1)';
plot(uGrid, rec1r(refIdx, :), '--', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Rank-1 recon');
rec2r = U_reg(:,1:2) * S_reg(1:2,1:2) * V_reg(:,1:2)';
plot(uGrid, rec2r(refIdx, :), '--', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Rank-2 recon');
xlabel('u (mm)');
ylabel('C(u,g) (pF)');
title('(d) Registered waveform SVD recon');
legend('Location', 'best', 'FontSize', 7);

% (e) Gap coefficients in PC space
nexttile; hold on;
Scores = U_reg * S_reg;  % nGap × nGap: score of gap i on PC k
PC1 = Scores(:, 1);
PC2 = Scores(:, 2);
for ig = 1:nGap
    plot(PC1(ig), PC2(ig), 'o', 'Color', colors(ig, :), ...
        'MarkerFaceColor', colors(ig, :), 'MarkerSize', 8);
    text(PC1(ig)+0.01, PC2(ig)+0.01, sprintf('%.1f', gapList(ig)), 'FontSize', 7);
end
xlabel('PC1 score');
ylabel('PC2 score');
title('(e) Gap states in PC1-PC2 space');

% (f) Shape PC space
nexttile; hold on;
ScoresS = U_shp * S_shp;  % nGap × nGap
PC1s = ScoresS(:, 1);
PC2s = ScoresS(:, 2);
for ig = 1:nGap
    plot(PC1s(ig), PC2s(ig), 'o', 'Color', colors(ig, :), ...
        'MarkerFaceColor', colors(ig, :), 'MarkerSize', 8);
    text(PC1s(ig)+0.001, PC2s(ig)+0.001, sprintf('%.1f', gapList(ig)), 'FontSize', 7);
end
xlabel('Shape PC1 score');
ylabel('Shape PC2 score');
title('(f) Gap states in shape PC space');

title(tl3, 'Q3: Low-rank structure of gap-waveform family', 'FontWeight', 'normal');
exportgraphics(fig3, fullfile(outDir, 'fig3_low_rank_analysis.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'fig3_low_rank_analysis.pdf'), 'ContentType', 'vector');

%% Figure 4: Q4 — Geometric feature trends
fig4 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl4 = tiledlayout(fig4, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Peak amplitude P(g) with best fit
nexttile; hold on;
plot(gapList, Peak, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
[yVal, gVal, bestFit, bestName, bestR2] = find_best_geometric_fit(Peak, gapList, candidateModels);
gFine = linspace(min(gapList)*0.7, max(gapList)*1.3, 200)';
yFine = bestFit(gFine);
plot(gFine, yFine, 'r-', 'LineWidth', 1.2);
xlabel('Gap g (mm)');
ylabel('Peak P(g) (pF)');
title(sprintf('(a) Peak: %s, R^2=%.4f', bestName, bestR2));

% (b) Area S(g) with best fit
nexttile; hold on;
plot(gapList, Area, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
[yVal2, gVal2, bestFit2, bestName2, bestR2_2] = find_best_geometric_fit(Area, gapList, candidateModels);
plot(gFine, bestFit2(gFine), 'r-', 'LineWidth', 1.2);
xlabel('Gap g (mm)');
ylabel('Area S(g) (pF·mm)');
title(sprintf('(b) Area: %s, R^2=%.4f', bestName2, bestR2_2));

% (c) FWHM with best fit
nexttile; hold on;
plot(gapList, FWHM, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
[yVal3, gVal3, bestFit3, bestName3, bestR2_3] = find_best_geometric_fit(FWHM, gapList, candidateModels);
plot(gFine, bestFit3(gFine), 'r-', 'LineWidth', 1.2);
xlabel('Gap g (mm)');
ylabel('FWHM (mm)');
title(sprintf('(c) FWHM: %s, R^2=%.4f', bestName3, bestR2_3));

% (d) Multi-level widths vs g
nexttile; hold on;
for iq = 1:nQ
    plot(gapList, Widths(:, iq), '-o', 'Color', colors(iq*2, :), ...
        'MarkerFaceColor', colors(iq*2, :), 'MarkerSize', 5, ...
        'DisplayName', sprintf('q=%.2f', qLevels(iq)));
end
xlabel('Gap g (mm)');
ylabel('Equal-level width w(q,g) (mm)');
title('(d) Multi-level widths vs gap');
legend('Location', 'best', 'FontSize', 7);

% (e) Width ratios (normalized signature) vs g
nexttile; hold on;
WidthNorm = Widths ./ sum(Widths, 2);
for iq = 1:nQ
    plot(gapList, WidthNorm(:, iq), '-o', 'Color', colors(iq*2, :), ...
        'MarkerFaceColor', colors(iq*2, :), 'MarkerSize', 5, ...
        'DisplayName', sprintf('q=%.2f', qLevels(iq)));
end
xlabel('Gap g (mm)');
ylabel('Normalized width w(q,g) / sum(w)');
title('(e) Normalized width signature vs gap');
legend('Location', 'best', 'FontSize', 7);

% (f) R² comparison for all features
nexttile;
R2feat = q4Results{:, 3:end};
bar(categorical(q4Results.feature), R2feat);
xlabel('Geometric feature');
ylabel('R^2');
title('(f) Model fit quality by feature');
legend(candidateModels, 'Location', 'eastoutside', 'FontSize', 6.0);
xtickangle(20);
ylim([0.8, 1.005]);

title(tl4, 'Q4: Geometric feature trends with gap', 'FontWeight', 'normal');
exportgraphics(fig4, fullfile(outDir, 'fig4_geometric_feature_trends.png'), 'Resolution', 300);
exportgraphics(fig4, fullfile(outDir, 'fig4_geometric_feature_trends.pdf'), 'ContentType', 'vector');

%% Figure 5: Q5 — Model type comparison
fig5 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl5 = tiledlayout(fig5, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a-c) Representative gaps: true vs Type A, B, C reconstruction to reference
repGaps = [0.4, 0.8, 1.2];
xcFine = linspace(min(xcAll), max(xcAll), 1000)';
for ir = 1:3
    [~, ig] = min(abs(gapList - repGaps(ir)));
    yRefTrue = interp1(curve(refIdx).xCenter, curve(refIdx).yBc, xcFine, 'pchip');

    % Type A
    yGi = interp1(curve(ig).xCenter, curve(ig).yBc, xcFine, 'pchip');
    yA = (Aref / max(max(abs(yGi)), eps)) * yGi;

    % Type B
    wRatio = FWHM(refIdx) / max(FWHM(ig), eps);
    xcW = xcFine * wRatio;
    yGiW = interp1(curve(ig).xCenter, curve(ig).yBc, xcW, 'linear', 'extrap');
    yB = (Aref / max(max(abs(yGiW)), eps)) * yGiW;

    % Type C
    yCrec = (Aref / max(Ameas(ig), eps)) * Cmat(ig, :);
    yC = interp1(uGrid, yCrec, xcFine, 'pchip', 'extrap');

    nexttile; hold on;
    plot(xcFine, yRefTrue, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True ref');
    plot(xcFine, yA, '--', 'Color', [0.85 0.35 0.10], 'DisplayName', 'Type A');
    plot(xcFine, yB, '--', 'Color', [0.10 0.60 0.40], 'DisplayName', 'Type B');
    plot(xcFine, yC, '--', 'Color', [0.10 0.45 0.75], 'DisplayName', 'Type C');
    xlabel('x - x_p (mm)');
    ylabel('C (pF)');
    title(sprintf('(%c) g = %.1f mm → ref', char('a'+ir-1), gapList(ig)));
    legend('Location', 'best', 'FontSize', 6.5);
end

% (d) Shape RMSE comparison (shape collapse quality)
nexttile;
bar(categorical(modelLabels), mean(q5ShapeRmse, 1)', 0.6, ...
    'FaceColor', [0.55 0.69 0.89], 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Mean shape RMSE');
title('(d) Shape collapse quality');
xtickangle(12);

% (e) Per-gap shape RMSE
nexttile; hold on;
plot(gapList, q5ShapeRmse(:,1), '-o', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'DisplayName', 'Type A');
plot(gapList, q5ShapeRmse(:,2), '-s', 'Color', [0.10 0.60 0.40], ...
    'MarkerFaceColor', [0.10 0.60 0.40], 'DisplayName', 'Type B');
plot(gapList, q5ShapeRmse(:,3), '-^', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Type C');
xlabel('Gap g (mm)');
ylabel('Shape RMSE');
title('(e) Shape RMSE by gap');
legend('Location', 'northwest', 'FontSize', 7);

% (f) Reconstruction RMSE (predicting reference waveform)
nexttile; hold on;
plot(gapList, q5ReconRmse(:,1), '-o', 'Color', [0.85 0.35 0.10], ...
    'MarkerFaceColor', [0.85 0.35 0.10], 'DisplayName', 'Type A');
plot(gapList, q5ReconRmse(:,2), '-s', 'Color', [0.10 0.60 0.40], ...
    'MarkerFaceColor', [0.10 0.60 0.40], 'DisplayName', 'Type B');
plot(gapList, q5ReconRmse(:,3), '-^', 'Color', [0.10 0.45 0.75], ...
    'MarkerFaceColor', [0.10 0.45 0.75], 'DisplayName', 'Type C');
xlabel('Gap g (mm)');
ylabel('Reconstruction RMSE (pF)');
title('(f) Reference-gap reconstruction error');
legend('Location', 'northwest', 'FontSize', 7);

title(tl5, 'Q5: Model-type comparison — A vs B vs C', 'FontWeight', 'normal');
exportgraphics(fig5, fullfile(outDir, 'fig5_model_type_comparison.png'), 'Resolution', 300);
exportgraphics(fig5, fullfile(outDir, 'fig5_model_type_comparison.pdf'), 'ContentType', 'vector');

%% ========================================================================
% CONSOLE SUMMARY
% ========================================================================
fprintf('\n========================================\n');
fprintf('GAP-WAVEFORM MAPPING STRUCTURE ANALYSIS\n');
fprintf('========================================\n');
fprintf('Data: 7 static gap waveforms, g = 0.2:0.2:1.4 mm\n');
fprintf('Reference gap: %.1f mm\n\n', gRef);

fprintf('--- Q1: Fixed-x point analysis ---\n');
fprintf('Best model counts across all x positions:\n');
[~, bestQ1] = max(R2mat, [], 2);
for im = 1:nModels
    fprintf('  %s: %d positions\n', modelNames{im}, sum(bestQ1 == im));
end
fprintf('Mean R² across positions: linear-g=%.4f, 1/g=%.4f, shifted-inv=%.4f\n', ...
    mean(q1Results.R2_linear_g), mean(q1Results.R2_linear_1overg), ...
    mean(q1Results.R2_shifted_inv));

fprintf('\n--- Q2: Coordinate comparison ---\n');
for ic = 1:numel(coordNames)
    fprintf('  %-25s: mean RMSE=%.4f, mean corr=%.4f\n', ...
        char(coordNames{ic}), q2Summary.mean_shape_rmse(ic), q2Summary.mean_shape_corr(ic));
end

fprintf('\n--- Q3: Low-rank ---\n');
fprintf('Variance in first 2 PCs:\n');
fprintf('  Raw peak-centered:   %.1f%%\n', 100*sum(varRaw(1:2)));
fprintf('  Registered C(u,g):   %.1f%%\n', 100*sum(varReg(1:2)));
fprintf('  Shape H(u,g):        %.1f%%\n', 100*sum(varShp(1:2)));
fprintf('Rank-1 relative error: raw=%.4f, registered=%.4f, shape=%.4f\n', ...
    recErrRaw(1), recErrReg(1), recErrShp(1));

fprintf('\n--- Q4: Geometric features ---\n');
for ife = 1:min(6, nFeat)
    [~, bestM] = max(q4Results{ife, 3:end});
    fprintf('  %-25s: best=%s, R²=%.4f\n', ...
        q4Results.feature{ife}, candidateModels{bestM}, q4Results{ife, 2+bestM});
end

fprintf('\n--- Q5: Model-type ---\n');
for im = 1:3
    fprintf('  %-30s: mean RMSE=%.4f, mean corr=%.4f\n', ...
        modelLabels{im}, q5Summary.mean_rmse_pF(im), q5Summary.mean_corr(im));
end

fprintf('\nResults saved to: %s\n', outDir);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function [gapList, xCell, yCell] = load_stacked_curves(filePath)
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
        xCell{i} = x(idx);
        yCell{i} = y(idx);
    end
end

function [base, amp, yNorm] = normalize_waveform(y)
    y = y(:);
    base = min(y);
    y0 = y - base;
    amp = max(y0);
    yNorm = y0 / max(amp, eps);
end

function [xLeft, xRight] = level_crossings_normalized(x, yNorm, levels)
    x = x(:);
    yNorm = yNorm(:);
    levels = levels(:);
    [~, iPeak] = max(yNorm);
    xRise = x(1:iPeak);
    yRise = yNorm(1:iPeak);
    xFall = x(iPeak:end);
    yFall = yNorm(iPeak:end);
    [yRiseU, idxRiseU] = unique(yRise, 'stable');
    xRiseU = xRise(idxRiseU);
    [yFallU, idxFallU] = unique(yFall, 'stable');
    xFallU = xFall(idxFallU);
    xLeft = nan(size(levels));
    xRight = nan(size(levels));
    for k = 1:numel(levels)
        lv = levels(k);
        if lv >= min(yRiseU) && lv <= max(yRiseU) && lv >= min(yFallU) && lv <= max(yFallU)
            xLeft(k) = interp1(yRiseU, xRiseU, lv, 'linear');
            xRight(k) = interp1(flipud(yFallU), flipud(xFallU), lv, 'linear');
        end
    end
end

function [curve, uGrid, Cmat, Hmat, Ameas] = ...
        build_equal_level_registration(curve, gapList, refIdx, levelMap, nU)
    nGap = numel(gapList);

    % Extract level anchors
    for ig = 1:nGap
        [xL, xR] = level_crossings_normalized(curve(ig).xCenter, curve(ig).yNorm, levelMap);
        curve(ig).xLeftMap = xL;
        curve(ig).xRightMap = xR;
    end

    % Reference grid
    refCurve = curve(refIdx);
    valid = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
    uMin = refCurve.xLeftMap(find(valid, 1, 'first'));
    uMax = refCurve.xRightMap(find(valid, 1, 'first'));
    uGrid = linspace(uMin, uMax, nU)';

    % Register all
    Cmat = nan(nGap, nU);
    for ig = 1:nGap
        valid = isfinite(curve(ig).xLeftMap) & isfinite(curve(ig).xRightMap) & ...
                isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);

        xLeftAnchor = [curve(ig).xLeftMap(valid); 0];
        uLeftAnchor = [refCurve.xLeftMap(valid); 0];
        [xLeftAnchor, leftIdx] = unique(xLeftAnchor, 'stable');
        uLeftAnchor = uLeftAnchor(leftIdx);

        xRightAnchor = [0; flipud(curve(ig).xRightMap(valid))];
        uRightAnchor = [0; flipud(refCurve.xRightMap(valid))];
        [xRightAnchor, rightIdx] = unique(xRightAnchor, 'stable');
        uRightAnchor = uRightAnchor(rightIdx);

        xCenter = curve(ig).xCenter(:);
        yBc = curve(ig).yBc(:);
        [~, iPeak] = max(yBc);

        uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenter(1:iPeak), 'pchip', 'extrap');
        uRight = interp1(xRightAnchor, uRightAnchor, xCenter(iPeak:end), 'pchip', 'extrap');
        uOfX = [uLeft; uRight(2:end)];
        yOfU = [yBc(1:iPeak); yBc(iPeak+1:end)];

        [uSorted, sortIdx] = sort(uOfX, 'ascend');
        ySorted = yOfU(sortIdx);
        [uSorted, uniqIdx] = unique(uSorted, 'stable');
        ySorted = ySorted(uniqIdx);
        Cmat(ig, :) = interp1(uSorted, ySorted, uGrid, 'pchip', 'extrap');
    end

    Ameas = max(Cmat, [], 2);
    Hmat = Cmat ./ max(Ameas, eps);
end

function res = fit_shifted_inv(y, g, delta)
    X = [ones(size(g)), 1./(g + delta)];
    b = X \ y;
    yhat = X * b;
    res = sum((y - yhat).^2);
end

function [bestIdx] = q1BestIdx(q1Results)
    R2mat = q1Results{:, 2:end};
    [~, bestIdx] = max(R2mat, [], 2);
end

function [yVal, gVal, bestFitFn, bestName, bestR2] = find_best_geometric_fit(yVec, gVec, modelNames)
    valid = isfinite(yVec);
    yVal = yVec(valid);
    gVal = gVec(valid);
    bestR2 = -inf;
    bestFitFn = @(g) nan(size(g));
    bestName = '';

    % Linear in g
    X = [ones(size(gVal)), gVal];
    b = X \ yVal;
    yhat = X * b;
    R2 = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
    if R2 > bestR2
        bestR2 = R2;
        bestFitFn = @(g) b(1) + b(2)*g;
        bestName = modelNames{1};
    end

    % Linear in 1/g
    X = [ones(size(gVal)), 1./gVal];
    b = X \ yVal;
    yhat = X * b;
    R2 = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
    if R2 > bestR2
        bestR2 = R2;
        bestFitFn = @(g) b(1) + b(2)./g;
        bestName = modelNames{2};
    end

    % Shifted inverse
    obj = @(d) fit_shifted_inv(yVal, gVal, d);
    [dOpt, ~] = fminsearch(obj, 0.3);
    X = [ones(size(gVal)), 1./(gVal + dOpt)];
    b = X \ yVal;
    yhat = X * b;
    R2 = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
    if R2 > bestR2
        bestR2 = R2;
        bestFitFn = @(g) b(1) + b(2)./(g + dOpt);
        bestName = modelNames{3};
    end

    % Power law
    if all(yVal > 0)
        X = [ones(size(gVal)), log(gVal)];
        b = X \ log(yVal);
        yhat = exp(b(1) + b(2)*log(gVal));
        R2 = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
        if R2 > bestR2
            bestR2 = R2;
            bestFitFn = @(g) exp(b(1) + b(2)*log(g));
            bestName = modelNames{4};
        end
    end

    % Exponential
    if all(yVal > 0)
        X = [ones(size(gVal)), gVal];
        b = X \ log(yVal);
        yhat = exp(b(1) + b(2)*gVal);
        R2 = 1 - sum((yVal - yhat).^2) / sum((yVal - mean(yVal)).^2);
        if R2 > bestR2
            bestR2 = R2;
            bestFitFn = @(g) exp(b(1) + b(2)*g);
            bestName = modelNames{5};
        end
    end
end
