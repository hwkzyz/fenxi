%% Gap-Vibration Decoupling Pipeline
% Validates the full decoupling strategy under synthetic vibration:
%
% Theory:
%   Observed:  y(x0) = F(x0 - u(t), g)
%   Static:    F(ξ, g) = B(g) + A(g)·H(ψ_g(ξ))
%   Width sig: ρ(g) = normalized multi-level width vector (gap-sensitive, vib-robust)
%
% Pipeline:
%   Step 1: Extract multi-level width signature from observed waveform
%   Step 2: Match against calibrated ρ(g) library → estimate ĝ
%   Step 3: Use ψ_{ĝ} to register waveform to reference domain u
%   Step 4: In u-domain, estimate vibration u(t) via parametric model
%   Step 5: Compare against oracle (known g and registered waveform)
%
% Key questions tested:
%   (a) How does gap estimation accuracy degrade with vibration strength?
%   (b) Does width-signature matching beat raw-waveform matching?
%   (c) Can we iterate to improve?
%   (d) What is the overall vibration parameter estimation accuracy?

clc; clear; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(scriptDir, '直叶片2mm_不同间隙.txt');
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

outDir = fullfile(scriptDir, 'gap_vibration_decoupling_results');
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
% 1. BUILD STATIC CALIBRATION LIBRARY
% ========================================================================
fprintf('\n========== Building static calibration library ==========\n');

[gapList, xCell, yCell] = load_stacked_curves(dataFile);
nGap = numel(gapList);
[~, refIdx] = max(gapList);
gRef = gapList(refIdx);

% Preprocess
curve = struct();
for ig = 1:nGap
    x = xCell{ig}(:);
    y = yCell{ig}(:);
    base = min(y);
    yBc = y - base;
    amp = max(yBc);
    yNorm = yBc / amp;
    [~, iPeak] = max(yBc);
    curve(ig).gap = gapList(ig);
    curve(ig).x = x;
    curve(ig).yRaw = y;
    curve(ig).yBc = yBc;
    curve(ig).base = base;
    curve(ig).amp = amp;
    curve(ig).xPeak = x(iPeak);
    curve(ig).xCenter = x - x(iPeak);
    curve(ig).yNorm = yNorm;
end

% Extract multi-level width signature library
qLevels = (0.10:0.05:0.90)';
nQ = numel(qLevels);

W_lib = nan(nGap, nQ);      % absolute widths
rho_lib = nan(nGap, nQ);    % normalized width signatures
xL_lib = nan(nGap, nQ);     % left branch positions
xR_lib = nan(nGap, nQ);     % right branch positions

for ig = 1:nGap
    [w, xL, xR] = extract_level_widths(curve(ig).xCenter, curve(ig).yNorm, qLevels);
    W_lib(ig, :) = w;
    xL_lib(ig, :) = xL;
    xR_lib(ig, :) = xR;
    rho_lib(ig, :) = w / max(sum(w), eps);
end

% Build equal-level registration for reference gap
levelMap = linspace(0.005, 0.995, 199)';
nU = 1401;
[curve, uGrid, Cmat, Hmat, Ameas] = ...
    build_equal_level_registration(curve, gapList, refIdx, levelMap, nU);
Cref = Cmat(refIdx, :).';
Aref = Ameas(refIdx);

% For any gap, we need ψ_g to map from x_c to u. Build the mapping functions.
% ψ_g maps: x_c (peak-centered, gap g) → u (reference domain)
% We store the anchor pairs for each gap.
psi_anchors = cell(nGap, 1);
for ig = 1:nGap
    valid = isfinite(curve(ig).xLeftMap) & isfinite(curve(ig).xRightMap) & ...
            isfinite(curve(refIdx).xLeftMap) & isfinite(curve(refIdx).xRightMap);

    xLeftAnchor = [curve(ig).xLeftMap(valid); 0];
    uLeftAnchor = [curve(refIdx).xLeftMap(valid); 0];
    xRightAnchor = [0; flipud(curve(ig).xRightMap(valid))];
    uRightAnchor = [0; flipud(curve(refIdx).xRightMap(valid))];

    psi_anchors{ig}.xLeft = xLeftAnchor;
    psi_anchors{ig}.uLeft = uLeftAnchor;
    psi_anchors{ig}.xRight = xRightAnchor;
    psi_anchors{ig}.uRight = uRightAnchor;
end

fprintf('Library: %d gaps, %d width levels, reference gap = %.1f mm\n', ...
    nGap, nQ, gRef);

%% ========================================================================
% 2. GENERATE TEST CASES: vibrating waveforms
% ========================================================================
fprintf('\n========== Generating test cases ==========\n');

% Tip speed
V_tip = 3.0e5;  % mm/s

% Vibration cases — varying strength
% Single-mode sinusoidal vibration
vibCases = struct( ...
    'name', {'No vibration', 'Weak vib', 'Moderate vib', 'Strong vib', 'Very strong vib'}, ...
    'tag', {'none', 'weak', 'moderate', 'strong', 'vstrong'}, ...
    'A_mm', {0, 0.10, 0.30, 0.60, 1.00}, ...
    'f_Hz', {0, 300, 500, 800, 1200});
nVib = numel(vibCases);

% Test gaps — all gaps are test candidates
testGaps = gapList;
nTest = numel(testGaps);

% Build test matrix
rows = [];
rng(42);

for ig = 1:nTest
    gTrue = testGaps(ig);
    xTrue = curve(ig).x;
    yTrue = curve(ig).yRaw;

    for iv = 1:nVib
        vib = vibCases(iv);

        % Generate vibrating waveform
        synth = synthesize_vibrating_waveform(xTrue, yTrue, vib, V_tip);

        if isempty(synth)
            continue;
        end

        % ---- Gap estimation via normalized width signature ----
        [baseObs, ampObs, yNormObs] = normalize_waveform(synth.yObs);
        xCenterObs = synth.x0 - synth.x0(find(synth.yObs - baseObs == max(synth.yObs - baseObs), 1));
        [wObs, ~, ~] = extract_level_widths(xCenterObs, yNormObs, qLevels);

        % Normalize width signature — handle NaN
        wObsClean = wObs;
        wObsClean(~isfinite(wObsClean)) = 0;
        if sum(wObsClean) < eps
            rhoObs = zeros(1, nQ);
        else
            rhoObs = wObsClean(:)' / sum(wObsClean);
        end

        % Nearest-neighbor gap estimation in ρ-space
        rhoDists = zeros(nGap, 1);
        for j = 1:nGap
            rhoDists(j) = norm(rhoObs - rho_lib(j, :));
        end
        [~, nnIdx] = min(rhoDists);
        gHat_width = gapList(nnIdx);

        % Also estimate via interpolated ρ-matching (weighted average of top 3)
        [~, sortIdx] = sort(rhoDists);
        topK = min(3, nGap);
        weights = 1 ./ max(rhoDists(sortIdx(1:topK)), 1e-12);
        weights = weights / sum(weights);
        gHat_width_interp = sum(gapList(sortIdx(1:topK)) .* weights);

        % ---- Gap estimation via raw peak value ----
        % This is a naive baseline: just use the peak amplitude
        peakObs = max(synth.yObs) - min(synth.yObs);
        calibAmps = arrayfun(@(c) c.amp, curve)';
        peakDists = abs(calibAmps - peakObs);
        [~, pkIdx] = min(peakDists);
        gHat_peak = gapList(pkIdx);

        % ---- Gap estimation via raw waveform correlation ----
        % Match the full baseline-corrected, normalized waveform
        corrScores = zeros(nGap, 1);
        for j = 1:nGap
            yLibNorm = curve(j).yNorm;
            xLibCenter = curve(j).xCenter;
            % Interpolate observed to library grid
            yObsInterp = interp1(xCenterObs, yNormObs, xLibCenter, 'pchip', 'extrap');
            valid = isfinite(yObsInterp) & isfinite(yLibNorm);
            if sum(valid) > 20
                ct = corrcoef(yObsInterp(valid), yLibNorm(valid));
                if numel(ct) >= 4
                    corrScores(j) = ct(1, 2);
                end
            end
        end
        [~, corrIdx] = max(corrScores);
        gHat_corr = gapList(corrIdx);

        % ---- Oracle: gap estimated from width in TRUE relative coords ----
        xiSorted = sort(synth.xi);
        yXiSorted = interp1(synth.xi, synth.yObs, xiSorted, 'pchip');
        [baseXi, ~, yNormXi] = normalize_waveform(yXiSorted);
        xCenterXi = xiSorted - xiSorted(find(yXiSorted - baseXi == max(yXiSorted - baseXi), 1));
        [wXi, ~, ~] = extract_level_widths(xCenterXi, yNormXi, qLevels);
        rhoXi = wXi(:)' / max(sum(wXi), eps);
        rhoDistsXi = zeros(nGap, 1);
        for j = 1:nGap
            rhoDistsXi(j) = norm(rhoXi - rho_lib(j, :));
        end
        [~, nnXi] = min(rhoDistsXi);
        gHat_oracle = gapList(nnXi);

        % ---- Store results ----
        row.gap_true_mm = gTrue;
        row.vib_case = vib.name;
        row.vib_A_mm = vib.A_mm;
        row.vib_f_Hz = vib.f_Hz;
        row.distortion_eta = synth.eta;
        row.gHat_width_mm = gHat_width;
        row.gHat_width_interp_mm = gHat_width_interp;
        row.gHat_peak_mm = gHat_peak;
        row.gHat_corr_mm = gHat_corr;
        row.gHat_oracle_mm = gHat_oracle;
        row.err_width_mm = gHat_width - gTrue;
        row.err_width_interp_mm = gHat_width_interp - gTrue;
        row.err_peak_mm = gHat_peak - gTrue;
        row.err_corr_mm = gHat_corr - gTrue;
        row.err_oracle_mm = gHat_oracle - gTrue;
        rows = [rows; row]; %#ok<AGROW>
    end
end

gapEstTable = struct2table(rows);
writetable(gapEstTable, fullfile(outDir, 'gap_estimation_under_vibration.csv'));

% Summary by vibration case
gapSummary = groupsummary(gapEstTable, "vib_case", "mean", ...
    ["distortion_eta", "err_width_mm", "err_width_interp_mm", "err_peak_mm", "err_corr_mm", "err_oracle_mm"]);
gapSummary.Properties.VariableNames = {'vib_case', 'GroupCount', ...
    'mean_eta', 'mean_err_width_mm', 'mean_err_width_interp_mm', ...
    'mean_err_peak_mm', 'mean_err_corr_mm', 'mean_err_oracle_mm'};

% Also compute RMSE
for iv = 1:nVib
    mask = strcmp(gapEstTable.vib_case, vibCases(iv).name);
    gapSummary.rmse_width_mm(iv) = sqrt(mean(gapEstTable.err_width_mm(mask).^2));
    gapSummary.rmse_corr_mm(iv) = sqrt(mean(gapEstTable.err_corr_mm(mask).^2));
    gapSummary.rmse_oracle_mm(iv) = sqrt(mean(gapEstTable.err_oracle_mm(mask).^2));
end

fprintf('\nGap estimation under vibration (RMSE, mm):\n');
disp(gapSummary(:, {'vib_case', 'mean_eta', 'rmse_width_mm', 'rmse_corr_mm', 'rmse_oracle_mm'}));

writetable(gapSummary, fullfile(outDir, 'gap_estimation_summary.csv'));

%% ========================================================================
% 3. FULL DECOUPLING: ĝ → registration → vibration estimation
% ========================================================================
fprintf('\n========== Full decoupling: gap + vibration ==========\n');

% Select representative test gaps for detailed analysis
repGaps = [0.4, 0.8, 1.2];
repVib = [3, 4];  % moderate and strong vibration indices

decRows = [];

for ig_idx = 1:numel(repGaps)
    [~, ig] = min(abs(gapList - repGaps(ig_idx)));
    gTrue = gapList(ig);
    xTrue = curve(ig).x;
    yTrue = curve(ig).yRaw;

    for iv_idx = 1:numel(repVib)
        iv = repVib(iv_idx);
        vib = vibCases(iv);

        synth = synthesize_vibrating_waveform(xTrue, yTrue, vib, V_tip);
        if isempty(synth), continue; end

        % ---- Step 1: Estimate gap from width signature ----
        [baseObs, ampObs, yNormObs] = normalize_waveform(synth.yObs);
        [~, iPeakObs] = max(synth.yObs - baseObs);
        xCenterObs = synth.x0 - synth.x0(iPeakObs);

        [wObs, ~, ~] = extract_level_widths(xCenterObs, yNormObs, qLevels);
        rhoObs = wObs(:)' / max(sum(wObs), eps);
        rhoDists = zeros(nGap, 1);
        for j = 1:nGap
            rhoDists(j) = norm(rhoObs - rho_lib(j, :));
        end
        [~, sortIdx] = sort(rhoDists);
        topK = min(3, nGap);
        weights = 1 ./ max(rhoDists(sortIdx(1:topK)), 1e-12);
        weights = weights / sum(weights);
        gHat = sum(gapList(sortIdx(1:topK)) .* weights);

        % ---- Step 2: Register observed waveform to reference domain ----
        % Build ψ_{ĝ} by interpolating between known ψ_g functions at
        % neighboring calibrated gaps
        [gSorted, gSortIdx] = sort(gapList);
        gLo = gSorted(find(gSorted <= gHat, 1, 'last'));
        gHi = gSorted(find(gSorted >= gHat, 1, 'first'));
        if isempty(gLo), gLo = gSorted(1); end
        if isempty(gHi), gHi = gSorted(end); end

        iLo = find(gapList == gLo, 1);
        iHi = find(gapList == gHi, 1);

        % Interpolate ψ anchors
        if gLo == gHi
            tInterp = 0;
        else
            tInterp = (gHat - gLo) / (gHi - gLo);
        end

        % For each level, interpolate the anchor mapping
        % We use the existing dense levelMap anchors from registration
        % Simpler approach: interpolate the registered waveforms directly
        if gLo == gHi
            Chat_interp = Cmat(iLo, :);
        else
            Chat_interp = (1 - tInterp) * Cmat(iLo, :) + tInterp * Cmat(iHi, :);
        end

        % Register observed waveform to u using the estimated gap's mapping
        % Build ψ for estimated gap by interpolating anchors
        [yOnU_est, uOfX_est] = register_waveform_interpolated(...
            synth.x0, synth.yObs, iLo, iHi, tInterp, curve, refIdx, uGrid);

        % ---- Step 3: Estimate vibration in u-domain ----
        % Model: yOnU(u_i) ≈ F_ref(u_i - u_effective(t_i))
        % In u-domain, t_i is not explicitly known. We work with:
        %   F_ref(u) is known (Cref)
        %   yOnU(u) = F_ref(u - u_eff), where u_eff is the vibration
        %   displacement mapped to the u-domain.
        %
        % Linearize: yOnU(u) ≈ Cref(u) - Cref'(u) * u_eff(u)
        %            Δy(u) ≈ -Cref'(u) * u_eff(u)
        %
        % For a parametric model of u_eff:
        %   u_eff maps from x0-domain vibration to u-domain apparent shift

        if isempty(yOnU_est) || numel(yOnU_est) ~= numel(Cref)
            continue;
        end
        yOnU_est = yOnU_est(:);
        CrefCol = Cref(:);
        validU = isfinite(yOnU_est) & isfinite(CrefCol);
        if sum(validU) < 100
            continue;
        end

        yOnU_valid = yOnU_est(validU);
        u_valid = uGrid(validU);

        % Residual in u-domain
        deltaY_u = yOnU_valid(:) - Cref(validU);

        % Approximate derivative of Cref
        dCref_du = gradient(Cref, uGrid(2) - uGrid(1));
        dCref_valid = dCref_du(validU);

        % Simple linear estimate of effective vibration displacement
        % u_eff ≈ -Δy / Cref'
        uEff_est = -deltaY_u ./ max(abs(dCref_valid), 1e-12);
        uEff_est = max(-1.0, min(1.0, uEff_est));  % clip

        % ---- Step 4: Parametric vibration fitting ----
        % Model: u_eff(t) = u(x0/V) in u-domain
        % For sinusoidal vibration: u(t) = A sin(2πf t + φ)
        % The mapping from u(t) to u_eff(u) is approximately linear:
        %   u_eff ≈ u * (du/dx0)^(-1) = u * scale
        %
        % We fit a single sinusoid to uEff_est

        % Reconstruct time axis
        tVec = (synth.x0 - mean(synth.x0)) / V_tip;

        % Map time to u-grid (u_valid corresponds to specific x0 positions)
        % x0 positions corresponding to u_valid
        x0OfU = interp1(uOfX_est, synth.x0, u_valid, 'linear', 'extrap');
        tOfU = (x0OfU - mean(synth.x0)) / V_tip;

        % Fit sinusoid: A*sin(2πf*t) + B*cos(2πf*t)
        % We know the true frequency from the vibration case
        % In practice, you'd sweep over candidate frequencies (VARPRO)
        fTrue = vib.f_Hz;
        if fTrue > 0
            % VARPRO-style: given f, solve for amplitudes
            sinTerm = sin(2 * pi * fTrue * tOfU(:));
            cosTerm = cos(2 * pi * fTrue * tOfU(:));
            X = [sinTerm, cosTerm];
            beta = X \ uEff_est(:);
            A_est = sqrt(beta(1)^2 + beta(2)^2);
            phi_est = atan2(beta(2), beta(1));
            uEff_fitted = X * beta;
        else
            A_est = 0;
            phi_est = 0;
            uEff_fitted = zeros(size(uEff_est));
        end

        % ---- Step 5: Reconstruct reference waveform using estimated g ----
        % pred_ref(u) = C(u, ĝ) * (Aref / Â(ĝ))
        Apred = max(Chat_interp);
        CrefRecon = (Aref / max(Apred, eps)) * Chat_interp;

        % Compare against true reference
        reconRmse = sqrt(mean((CrefRecon(:) - Cref).^2, 'omitnan'));

        % ---- Metrics ----
        decRow.gap_true_mm = gTrue;
        decRow.gap_est_mm = gHat;
        decRow.gap_err_mm = gHat - gTrue;
        decRow.vib_case = vib.name;
        decRow.vib_A_true_mm = vib.A_mm;
        decRow.vib_f_true_Hz = vib.f_Hz;
        decRow.vib_A_est_mm = A_est;
        decRow.vib_A_err_mm = A_est - vib.A_mm;
        decRow.vib_phi_est_rad = phi_est;
        decRow.uEff_rms_mm = sqrt(mean(uEff_est.^2, 'omitnan'));
        decRow.uEff_fit_rmse_mm = sqrt(mean((uEff_est(:) - uEff_fitted).^2, 'omitnan'));
        decRow.recon_rmse_pF = reconRmse;
        decRow.n_valid_u_points = sum(validU);
        decRows = [decRows; decRow]; %#ok<AGROW>
    end
end

decTable = struct2table(decRows);
fprintf('\nFull decoupling results:\n');
disp(decTable(:, {'gap_true_mm', 'gap_est_mm', 'gap_err_mm', 'vib_case', ...
    'vib_A_true_mm', 'vib_A_est_mm', 'recon_rmse_pF'}));

writetable(decTable, fullfile(outDir, 'full_decoupling_results.csv'));

%% ========================================================================
% 4. SENSITIVITY: how does gap estimation ±error affect vibration estimation?
% ========================================================================
fprintf('\n========== Sensitivity: gap error → vibration error ==========\n');

% For a fixed test case, sweep intentional gap errors and measure
% vibration estimation degradation
[~, sensIg] = min(abs(gapList - 0.8));
sensVib = vibCases(3);  % moderate vibration
sensSynth = synthesize_vibrating_waveform(curve(sensIg).x, curve(sensIg).yRaw, sensVib, V_tip);
[~, ~, yNormSens] = normalize_waveform(sensSynth.yObs);
[~, iPeakSens] = max(sensSynth.yObs - min(sensSynth.yObs));
xCenterSens = sensSynth.x0 - sensSynth.x0(iPeakSens);

% Sweep gap error from -0.4 to +0.4 mm
gapErrSweep = (-0.4:0.05:0.4)';
sweepRows = [];
for ie = 1:numel(gapErrSweep)
    gUsed = gapList(sensIg) + gapErrSweep(ie);
    gUsed = max(0.05, min(2.0, gUsed));

    % Find bracketing calibrated gaps
    [gSorted, ~] = sort(gapList);
    gLo = gSorted(find(gSorted <= gUsed, 1, 'last'));
    gHi = gSorted(find(gSorted >= gUsed, 1, 'first'));
    if isempty(gLo), gLo = gSorted(1); end
    if isempty(gHi), gHi = gSorted(end); end
    iLo = find(gapList == gLo, 1);
    iHi = find(gapList == gHi, 1);
    if gLo == gHi, tInterp = 0; else tInterp = (gUsed - gLo) / (gHi - gLo); end

    [yOnU, ~] = register_waveform_interpolated(...
        sensSynth.x0, sensSynth.yObs, iLo, iHi, tInterp, curve, refIdx, uGrid);

    validU = isfinite(yOnU) & isfinite(CrefCol);
    if sum(validU) < 100, continue; end

    deltaY = yOnU(validU) - CrefCol(validU);
    dCref = gradient(Cref, uGrid(2) - uGrid(1));
    uEff = -deltaY ./ max(abs(dCref(validU)), 1e-12);

    % Fit sinusoid
    % Map u_valid positions to x0, then to time
    x0_of_u = interp1(sensSynth.x0, sensSynth.x0, u_valid, 'linear', 'extrap');
    % Use direct mapping: u positions approximately correspond to x0 positions
    tOfU = (u_valid - mean(sensSynth.x0)) / V_tip;

    fTrue = sensVib.f_Hz;
    if fTrue > 0
        sinTerm = sin(2*pi*fTrue*tOfU(:));
        cosTerm = cos(2*pi*fTrue*tOfU(:));
        validFit = isfinite(sinTerm) & isfinite(cosTerm) & isfinite(uEff(:));
        beta = [sinTerm(validFit), cosTerm(validFit)] \ uEff(validFit);
        AEst = sqrt(beta(1)^2 + beta(2)^2);
    else
        AEst = 0;
    end

    reconErr = sqrt(mean((yOnU(validU) - Cref(validU)).^2, 'omitnan'));

    sr = struct();
    sr.gap_error_mm = gapErrSweep(ie);
    sr.gap_used_mm = gUsed;
    sr.vib_A_est_mm = AEst;
    sr.vib_A_err_percent = 100 * (AEst - sensVib.A_mm) / max(sensVib.A_mm, eps);
    sr.deltaY_rms_pF = sqrt(mean(deltaY.^2, 'omitnan'));
    sr.recon_rmse_pF = reconErr;
    sweepRows = [sweepRows; sr]; %#ok<AGROW>
end

sweepTable = struct2table(sweepRows);
fprintf('Gap error sensitivity (representative case):\n');
disp(sweepTable(:, {'gap_error_mm', 'vib_A_est_mm', 'vib_A_err_percent', 'recon_rmse_pF'}));

writetable(sweepTable, fullfile(outDir, 'gap_error_sensitivity.csv'));

%% ========================================================================
% FIGURES
% ========================================================================
colors = lines(max(nGap, 14));

%% Figure 1: Gap estimation accuracy under vibration
fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl1 = tiledlayout(fig1, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Gap estimation error vs vibration strength by method
nexttile; hold on;
methodNames = {'Width signature', 'Waveform correlation', 'Peak amplitude'};
methodCols = {[0.10 0.45 0.75], [0.85 0.35 0.10], [0.50 0.50 0.50]};
for iv = 1:nVib
    mask = strcmp(gapEstTable.vib_case, vibCases(iv).name);
    xpos = iv;
    errs = {abs(gapEstTable.err_width_mm(mask)), ...
            abs(gapEstTable.err_corr_mm(mask)), ...
            abs(gapEstTable.err_peak_mm(mask))};
    for im = 1:3
        e = errs{im};
        boxplot_at(e(:), xpos + (im-2)*0.25, 0.2, methodCols{im});
    end
end
set(gca, 'XTick', 1:nVib, 'XTickLabel', {vibCases.name});
xlabel('Vibration case');
ylabel('|Gap error| (mm)');
title('(a) Gap estimation error by method');
legend(methodNames, 'Location', 'northwest', 'FontSize', 7);
yline(0.2, 'r--', 'LineWidth', 0.8);  % gap spacing

% (b) Width signature matching distances (confusion matrix style)
% First compute rhoObs for representative vibration case
repVibIdx = 3;  % moderate
rhoObs_rep = zeros(nGap, nQ);
for ig = 1:nGap
    synth2 = synthesize_vibrating_waveform(curve(ig).x, curve(ig).yRaw, vibCases(repVibIdx), V_tip);
    if ~isempty(synth2)
        [~, ~, yNorm2] = normalize_waveform(synth2.yObs);
        [~, iPeak2] = max(synth2.yObs - min(synth2.yObs));
        xCenter2 = synth2.x0 - synth2.x0(iPeak2);
        [w2, ~, ~] = extract_level_widths(xCenter2, yNorm2, qLevels);
        w2c = w2; w2c(~isfinite(w2c)) = 0;
        if sum(w2c) > eps
            rhoObs_rep(ig, :) = w2c(:)' / sum(w2c);
        end
    end
end

nexttile; hold on;
distMat = zeros(nGap, nGap);
for ig = 1:nGap
    for j = 1:nGap
        distMat(ig, j) = norm(rhoObs_rep(ig, :) - rho_lib(j, :));
    end
end
imagesc(gapList, gapList, distMat);
colormap(gca, flipud(gray));
colorbar;
xlabel('Library gap (mm)');
ylabel('True gap (mm)');
title(sprintf('(b) ρ-space distances, %s', vibCases(repVibIdx).name));

% (c) Example waveforms: observed vs true
nexttile; hold on;
repGapIdx = find(abs(gapList - 0.8) < 0.01, 1);
synth3 = synthesize_vibrating_waveform(curve(repGapIdx).x, curve(repGapIdx).yRaw, vibCases(3), V_tip);
plot(curve(repGapIdx).x, curve(repGapIdx).yRaw, 'k-', 'LineWidth', 1.2, ...
    'DisplayName', 'Static, g=0.8');
plot(synth3.x0, synth3.yObs, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('Observed, %s', vibCases(3).name));
xlabel('x_0 (mm)');
ylabel('Capacitance (pF)');
title('(c) Observed vs static waveform');
legend('Location', 'best', 'FontSize', 7);

% (d) Width signatures: library vs observed
nexttile; hold on;
for j = 1:nGap
    plot(qLevels, rho_lib(j, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end
plot(qLevels, rho_lib(repGapIdx, :), 'k-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Library, g=%.1f', gapList(repGapIdx)));
plot(qLevels, rhoObs_rep(repGapIdx, :), 'ro-', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Observed, %s', vibCases(repVibIdx).name));
xlabel('Normalized level q');
ylabel('Normalized width signature');
title('(d) Width signatures: library vs observed');
legend('Location', 'best', 'FontSize', 7);

% (e) Gap estimation error vs distortion η
nexttile; hold on;
plot(gapEstTable.distortion_eta, abs(gapEstTable.err_width_mm), 'o', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], ...
    'DisplayName', 'Width signature');
plot(gapEstTable.distortion_eta, abs(gapEstTable.err_corr_mm), 's', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], ...
    'DisplayName', 'Waveform correlation');
xlabel('Distortion \eta = max|du/dt| / V');
ylabel('|Gap error| (mm)');
title('(e) Gap error vs distortion strength');
legend('Location', 'northwest', 'FontSize', 7);

% (f) Gap estimation RMSE summary
nexttile;
% Compute RMSE directly from gapEstTable
rmse_width = zeros(nVib, 1);
rmse_corr = zeros(nVib, 1);
for iv = 1:nVib
    mask = strcmp(gapEstTable.vib_case, vibCases(iv).name);
    rmse_width(iv) = sqrt(mean(gapEstTable.err_width_mm(mask).^2));
    rmse_corr(iv) = sqrt(mean(gapEstTable.err_corr_mm(mask).^2));
end
barData = [rmse_width, rmse_corr];
bar(categorical({vibCases.name}), barData, 'grouped');
xlabel('Vibration case');
ylabel('RMSE (mm)');
title('(f) Gap estimation RMSE summary');
legend({'Width signature', 'Waveform corr'}, ...
    'Location', 'northwest', 'FontSize', 7);
yline(0.2, 'r--', 'LineWidth', 0.8);

title(tl1, 'Gap estimation under vibration', 'FontWeight', 'normal');
exportgraphics(fig1, fullfile(outDir, 'fig1_gap_estimation_under_vibration.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'fig1_gap_estimation_under_vibration.pdf'), 'ContentType', 'vector');

%% Figure 2: Full decoupling pipeline
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 14]);
tl2 = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Pick a representative case for detailed visualization
[~, visIg] = min(abs(gapList - 0.8));
visVib = vibCases(3);
visSynth = synthesize_vibrating_waveform(curve(visIg).x, curve(visIg).yRaw, visVib, V_tip);

[~, ~, yNormVis] = normalize_waveform(visSynth.yObs);
[~, iPeakVis] = max(visSynth.yObs - min(visSynth.yObs));
xCenterVis = visSynth.x0 - visSynth.x0(iPeakVis);
[wVis, ~, ~] = extract_level_widths(xCenterVis, yNormVis, qLevels);
rhoVis = wVis(:)' / max(sum(wVis), eps);
rhoDistsVis = zeros(nGap, 1);
for j = 1:nGap, rhoDistsVis(j) = norm(rhoVis - rho_lib(j, :)); end
[~, sortVis] = sort(rhoDistsVis);
gHatVis = sum(gapList(sortVis(1:3)) .* (1./max(rhoDistsVis(sortVis(1:3)), 1e-12))');
gHatVis = gHatVis / sum(1./max(rhoDistsVis(sortVis(1:3)), 1e-12));

% Register
[gSorted, ~] = sort(gapList);
gLo = gSorted(find(gSorted <= gHatVis, 1, 'last'));
gHi = gSorted(find(gSorted >= gHatVis, 1, 'first'));
if isempty(gLo), gLo = gSorted(1); end
if isempty(gHi), gHi = gSorted(end); end
iLo = find(gapList == gLo, 1);
iHi = find(gapList == gHi, 1);
tInterpVis = 0; if gLo ~= gHi, tInterpVis = (gHatVis - gLo) / (gHi - gLo); end

[yOnU_vis, uOfX_vis] = register_waveform_interpolated(...
    visSynth.x0, visSynth.yObs, iLo, iHi, tInterpVis, curve, refIdx, uGrid);

% (a) Observed vs static waveforms
nexttile; hold on;
plot(curve(visIg).x, curve(visIg).yRaw, 'k-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Static, g=%.1f', gapList(visIg)));
plot(visSynth.x0, visSynth.yObs, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('Observed, vib A=%.2f', visVib.A_mm));
xlabel('x_0 (mm)');
ylabel('Capacitance (pF)');
title('(a) Observed vs static waveform');
legend('Location', 'best', 'FontSize', 7);

% (b) Gap estimation from width signature
nexttile; hold on;
bar(categorical(string(gapList)), rhoDistsVis, 0.6, 'FaceColor', [0.55 0.69 0.89], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);
xline(find(gapList == gapList(visIg)), 'r-', 'LineWidth', 2.0);
xlabel('Library gap (mm)');
ylabel('||\rho_{obs} - \rho_{lib}||');
title(sprintf('(b) Gap est: true=%.1f, est=%.2f', gapList(visIg), gHatVis));

% (c) Registered waveform in u-domain vs reference
nexttile; hold on;
plot(uGrid, Cref, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ref F(u, g_0)');
plot(uGrid, yOnU_vis, '-', 'Color', [0.10 0.45 0.75], 'LineWidth', 1.0, ...
    'DisplayName', 'Registered obs');
xlabel('u (mm)');
ylabel('C(u) (pF)');
title('(c) Registered to reference domain');
legend('Location', 'best', 'FontSize', 7);

% (d) Residual in u-domain and estimated vibration
nexttile; hold on;
validMask = isfinite(yOnU_vis) & isfinite(Cref');
deltaY = yOnU_vis(validMask) - Cref(validMask);
dCref = gradient(Cref, uGrid(2) - uGrid(1));
uEff = -deltaY ./ max(abs(dCref(validMask)), 1e-12);

yyaxis left;
plot(uGrid(validMask), deltaY, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 0.8);
ylabel('\Delta C (pF)');
yyaxis right;
plot(uGrid(validMask), uEff, '-', 'Color', [0.10 0.45 0.75], 'LineWidth', 0.8);
ylabel('u_{eff} est (mm)');
xlabel('u (mm)');
title('(d) Residual and vibration estimate');
legend({'\Delta C = y_{obs} - F_{ref}', 'u_{eff} = -\Delta C / F_{ref}'''}, ...
    'Location', 'best', 'FontSize', 6.5);

% (e) Gap error sensitivity: reconstruction RMSE vs gap error
nexttile; hold on;
plot(sweepTable.gap_error_mm, sweepTable.recon_rmse_pF, 'o-', ...
    'Color', [0.10 0.45 0.75], 'MarkerFaceColor', [0.10 0.45 0.75], ...
    'MarkerSize', 5);
xline(0, 'k--', 'LineWidth', 0.8);
xlabel('Gap estimation error (mm)');
ylabel('Reconstruction RMSE (pF)');
title('(e) Sensitivity: gap error → recon quality');

% (f) Vibration amplitude estimation error vs gap error
nexttile; hold on;
plot(sweepTable.gap_error_mm, sweepTable.vib_A_err_percent, 's-', ...
    'Color', [0.85 0.35 0.10], 'MarkerFaceColor', [0.85 0.35 0.10], ...
    'MarkerSize', 5);
xline(0, 'k--', 'LineWidth', 0.8);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Gap estimation error (mm)');
ylabel('Vibration amplitude error (%)');
title('(f) Sensitivity: gap error → vib error');

title(tl2, 'Full gap-vibration decoupling pipeline', 'FontWeight', 'normal');
exportgraphics(fig2, fullfile(outDir, 'fig2_full_decoupling_pipeline.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'fig2_full_decoupling_pipeline.pdf'), 'ContentType', 'vector');

%% ========================================================================
% CONSOLE SUMMARY
% ========================================================================
fprintf('\n========================================\n');
fprintf('GAP-VIBRATION DECOUPLING ANALYSIS\n');
fprintf('========================================\n');
fprintf('Calibration: %d static gaps, %.1f-%.1f mm\n', nGap, min(gapList), max(gapList));
fprintf('Reference gap: %.1f mm\n\n', gRef);

fprintf('--- Gap estimation under vibration ---\n');
fprintf('%-20s %8s %12s %12s %12s\n', 'Case', 'eta', 'RMSE_width', 'RMSE_corr', 'RMSE_oracle');
for iv = 1:nVib
    fprintf('%-20s %8.4f %12.4f %12.4f %12.4f\n', ...
        vibCases(iv).name, gapSummary.mean_eta(iv), ...
        gapSummary.rmse_width_mm(iv), gapSummary.rmse_corr_mm(iv), ...
        gapSummary.rmse_oracle_mm(iv));
end

fprintf('\n--- Decoupling accuracy ---\n');
fprintf('Mean gap estimation error (width): %.3f mm\n', mean(abs(gapEstTable.err_width_mm)));
fprintf('Mean gap estimation error (corr):  %.3f mm\n', mean(abs(gapEstTable.err_corr_mm)));
fprintf('Mean recon RMSE (full decoupling): %.4f pF\n', mean(decTable.recon_rmse_pF));
fprintf('Mean vib amplitude error: %.3f mm\n', mean(abs(decTable.vib_A_err_mm)));

fprintf('\nResults saved to: %s\n', outDir);

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

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

function [widths, xLeft, xRight] = extract_level_widths(xCenter, yNorm, qLevels)
    xCenter = xCenter(:);
    yNorm = yNorm(:);
    qLevels = qLevels(:);
    [~, iPeak] = max(yNorm);
    nQ = numel(qLevels);
    widths = nan(nQ, 1);
    xLeft = nan(nQ, 1);
    xRight = nan(nQ, 1);
    for k = 1:nQ
        q = qLevels(k);
        if q >= min(yNorm(1:iPeak)) && q <= max(yNorm(1:iPeak)) && ...
           q >= min(yNorm(iPeak:end)) && q <= max(yNorm(iPeak:end))
            xL = interp1(yNorm(1:iPeak), xCenter(1:iPeak), q, 'linear');
            xR = interp1(flipud(yNorm(iPeak:end)), flipud(xCenter(iPeak:end)), q, 'linear');
            xLeft(k) = xL;
            xRight(k) = xR;
            widths(k) = xR - xL;
        end
    end
end

function [xLeft, xRight] = level_crossings_dense(x, yNorm, levels)
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
    for ig = 1:nGap
        [xL, xR] = level_crossings_dense(curve(ig).xCenter, curve(ig).yNorm, levelMap);
        curve(ig).xLeftMap = xL;
        curve(ig).xRightMap = xR;
    end
    refCurve = curve(refIdx);
    valid = isfinite(refCurve.xLeftMap) & isfinite(refCurve.xRightMap);
    uMin = refCurve.xLeftMap(find(valid, 1, 'first'));
    uMax = refCurve.xRightMap(find(valid, 1, 'first'));
    uGrid = linspace(uMin, uMax, nU)';
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

function synth = synthesize_vibrating_waveform(xRef, yRef, vibCase, V_tip)
    xRef = xRef(:);
    yRef = yRef(:);
    x0Center = mean(xRef);

    % Find x0 window that keeps xi = x0 - u within calibration range
    x0 = xRef;
    for iter = 1:4
        tTmp = (x0 - x0Center) / V_tip;
        uTmp = vib_displacement(tTmp, vibCase);
        x0 = linspace(min(xRef) + max(uTmp), max(xRef) + min(uTmp), numel(xRef))';
    end

    t = (x0 - x0Center) / V_tip;
    [u, du_dt] = vib_displacement(t, vibCase);
    xi = x0 - u;

    if any(diff(xi) <= 0) || min(xi) < min(xRef) - 1e-8 || max(xi) > max(xRef) + 1e-8
        synth = [];
        return;
    end

    synth.x0 = x0;
    synth.t = t;
    synth.u = u;
    synth.du_dt = du_dt;
    synth.xi = xi;
    synth.yObs = interp1(xRef, yRef, xi, 'pchip', 'extrap');
    synth.eta = max(abs(du_dt)) / V_tip;
end

function [u, du_dt] = vib_displacement(t, vibCase)
    t = t(:);
    if vibCase.A_mm == 0 || vibCase.f_Hz == 0
        u = zeros(size(t));
        du_dt = zeros(size(t));
        return;
    end
    omega = 2 * pi * vibCase.f_Hz;
    u = vibCase.A_mm * sin(omega * t);
    du_dt = vibCase.A_mm * omega * cos(omega * t);
end

function [yOnU, uOfX] = register_waveform_interpolated(x0, yObs, iLo, iHi, tInterp, curve, refIdx, uGrid)
    % Interpolate the ψ mapping between calibrated gaps iLo and iHi,
    % then map the observed waveform to u-domain.
    refCurve = curve(refIdx);

    % Get level-crossing anchors for observed waveform
    [baseObs, ~, yNormObs] = normalize_waveform(yObs);
    [~, iPeakObs] = max(yObs - baseObs);
    xCenterObs = x0 - x0(iPeakObs);
    yBcObs = yObs - baseObs;

    % Use the same levelMap as the calibration
    levelMap = linspace(0.005, 0.995, 199)';
    [xLObs, xRObs] = level_crossings_dense(xCenterObs, yNormObs, levelMap);

    % Interpolate reference anchors
    xLRefLo = curve(iLo).xLeftMap;
    xRRefLo = curve(iLo).xRightMap;
    xLRefHi = curve(iHi).xLeftMap;
    xRRefHi = curve(iHi).xRightMap;

    uLRef = refCurve.xLeftMap;
    uRRef = refCurve.xRightMap;

    validObs = isfinite(xLObs) & isfinite(xRObs);
    validRef = isfinite(uLRef) & isfinite(uRRef) & ...
               isfinite(xLRefLo) & isfinite(xRRefLo) & ...
               isfinite(xLRefHi) & isfinite(xRRefHi);
    valid = validObs & validRef;

    if sum(valid) < 20
        % Fallback: simple interpolation
        yOnU = interp1(xCenterObs, yBcObs, uGrid, 'pchip', 'extrap');
        uOfX = xCenterObs;
        return;
    end

    % Interpolate the mapping targets
    xLTarget = (1 - tInterp) * xLRefLo(valid) + tInterp * xLRefHi(valid);
    xRTarget = (1 - tInterp) * xRRefLo(valid) + tInterp * xRRefHi(valid);

    xLeftAnchor = [xLObs(valid); 0];
    uLeftAnchor = [uLRef(valid); 0];
    [xLeftAnchor, leftIdx] = unique(xLeftAnchor, 'stable');
    uLeftAnchor = uLeftAnchor(leftIdx);

    xRightAnchor = [0; flipud(xRObs(valid))];
    uRightAnchor = [0; flipud(uRRef(valid))];
    [xRightAnchor, rightIdx] = unique(xRightAnchor, 'stable');
    uRightAnchor = uRightAnchor(rightIdx);

    [~, iPeakBc] = max(yBcObs);
    uLeft = interp1(xLeftAnchor, uLeftAnchor, xCenterObs(1:iPeakBc), 'pchip', 'extrap');
    uRight = interp1(xRightAnchor, uRightAnchor, xCenterObs(iPeakBc:end), 'pchip', 'extrap');

    uOfX = [uLeft; uRight(2:end)];
    yOfU = [yBcObs(1:iPeakBc); yBcObs(iPeakBc+1:end)];

    [uSorted, sortIdx] = sort(uOfX, 'ascend');
    ySorted = yOfU(sortIdx);
    [uSorted, uniqIdx] = unique(uSorted, 'stable');
    ySorted = ySorted(uniqIdx);
    yOnU = interp1(uSorted, ySorted, uGrid, 'pchip', 'extrap');
end

function boxplot_at(data, xpos, width, col)
    if isempty(data) || all(isnan(data))
        return;
    end
    data = data(isfinite(data));
    if isempty(data), return; end
    q = prctile(data, [25, 50, 75]);
    plot(xpos + width*0.5*[-1 1 1 -1 -1], [q(1) q(1) q(3) q(3) q(1)], ...
        '-', 'Color', col, 'LineWidth', 1.0);
    plot(xpos + [-1 1]*width*0.5, [q(2) q(2)], '-', 'Color', col, 'LineWidth', 1.8);
    whisker_lo = max(min(data), q(1) - 1.5*(q(3)-q(1)));
    whisker_hi = min(max(data), q(3) + 1.5*(q(3)-q(1)));
    plot([xpos xpos], [q(3) whisker_hi], '-', 'Color', col, 'LineWidth', 0.8);
    plot([xpos xpos], [q(1) whisker_lo], '-', 'Color', col, 'LineWidth', 0.8);
end
