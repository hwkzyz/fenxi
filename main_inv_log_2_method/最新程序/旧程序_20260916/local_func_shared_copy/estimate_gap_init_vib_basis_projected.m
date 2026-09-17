function gapState = estimate_gap_init_vib_basis_projected(highMap, templateLib, cfg, staticState)
%estimate_gap_init_vib_basis_projected  Gap correction with vibration-shaped nuisance bases.
%
% For fixed dx and candidate g, the first-order vibration term has the form
% F'_g(x-dx) * u(t). This estimator suppresses only the part of the residual
% explained by harmonic vibration bases, instead of removing an unrestricted
% translation direction.

if nargin < 4 || isempty(staticState)
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

x = highMap.x_v(:);
t = highMap.t_v(:);
V = highMap.V_a(:);
dxHat = staticState.dx0;

gapMin = max(0.05, min(templateLib.gapTrain) - cfg.rawGapSearchMargin);
gapMax = max(templateLib.gapTrain) + cfg.rawGapSearchMargin;
gAnchor = staticState.gHat;
if isempty(gAnchor) || ~isfinite(gAnchor)
    gAnchor = 0.5 * (gapMin + gapMax);
end

[freqList, freqMeta] = get_frequency_list(cfg, highMap, templateLib, staticState);
fullDomain = get_cfg_field(cfg, 'gapProjFullDomain', false);
coarseHalfWidth = get_cfg_field(cfg, 'gapProjHalfWidth', 0.08);
coarseN = get_cfg_field(cfg, 'gapProjCoarseN', 61);
fineHalfWidth = get_cfg_field(cfg, 'gapProjFineHalfWidth', 0.02);
fineN = get_cfg_field(cfg, 'gapProjFineN', 31);

if fullDomain
    coarseGapGrid = linspace(gapMin, gapMax, ...
        get_cfg_field(cfg, 'gapProjFullDomainN', max(coarseN, 121)));
else
    coarseGapGrid = linspace(max(gapMin, gAnchor - coarseHalfWidth), ...
        min(gapMax, gAnchor + coarseHalfWidth), coarseN);
end
sampleIdx = pick_evenly_spaced_indices(numel(x), get_cfg_field(cfg, 'gapProjTargetSamples', 5000));
JCoarse = eval_vib_basis_cost_grid(x(sampleIdx), t(sampleIdx), V(sampleIdx), ...
    templateLib, dxHat, coarseGapGrid, freqList, cfg);
[~, iCoarseBest] = min(JCoarse);
gCoarse = coarseGapGrid(iCoarseBest);

coarseStep = coarseGapGrid(min(end, 2)) - coarseGapGrid(1);
fineGapGrid = linspace(max(gapMin, gCoarse - max(fineHalfWidth, coarseStep)), ...
    min(gapMax, gCoarse + max(fineHalfWidth, coarseStep)), fineN);
JFine = eval_vib_basis_cost_grid(x, t, V, templateLib, dxHat, fineGapGrid, freqList, cfg);
[bestCost, iBest] = min(JFine);

gapState = struct();
gapState.gHat = fineGapGrid(iBest);
gapState.g_used = fineGapGrid(iBest);
gapState.dx0 = dxHat;
gapState.dx_used = dxHat;
gapState.g_static = staticState.gHat;
gapState.dx_static = staticState.dx0;
gapState.freq_list = freqList(:);
gapState.freq_meta = freqMeta;
gapState.gap_grid = fineGapGrid(:);
gapState.J = JFine;
gapState.coarse_gap_grid = coarseGapGrid(:);
gapState.coarse_J = JCoarse;
gapState.coarse_best_g = gCoarse;
gapState.best_cost = bestCost;
profileTol = get_cfg_field(cfg, 'gapProjProfileRelativeTolerance', 0.10);
profileFloor = get_cfg_field(cfg, 'gapProjProfileAbsoluteTolerance', 0);
profileKeep = JCoarse <= min(JCoarse) * (1 + profileTol) + profileFloor;
if any(profileKeep)
    gapState.gap_profile_lower = min(coarseGapGrid(profileKeep));
    gapState.gap_profile_upper = max(coarseGapGrid(profileKeep));
else
    gapState.gap_profile_lower = gapState.gHat;
    gapState.gap_profile_upper = gapState.gHat;
end
gapState.gap_profile_half_width = max( ...
    gapState.gHat - gapState.gap_profile_lower, ...
    gapState.gap_profile_upper - gapState.gHat);
gapState.gap_profile_full_domain = fullDomain;
gapState.staticState = staticState;
end

function J = eval_vib_basis_cost_grid(x, t, V, templateLib, dxHat, gapGrid, freqList, cfg)
J = zeros(numel(gapGrid), 1);
includeDx = get_cfg_field(cfg, 'gapVibProjIncludeDx', true);
for ig = 1:numel(gapGrid)
    g = gapGrid(ig);
    xq = x - dxHat;
    F = eval_gap_template(templateLib, g, xq);
    Fx = eval_gap_derivative(templateLib, g, xq);
    valid = isfinite(F) & isfinite(Fx) & isfinite(V);
    r = V(valid) - F(valid);
    B = build_vib_basis(Fx(valid), t(valid), freqList, includeDx);
    if isempty(B)
        rFit = r;
    else
        coef = B \ r;
        rFit = r - B * coef;
    end
    J(ig) = mean(rFit.^2);
end
end

function B = build_vib_basis(Fx, t, freqList, includeDx)
cols = {};
if includeDx
    cols{end+1} = Fx; %#ok<AGROW>
end
for ifr = 1:numel(freqList)
    w = 2 * pi * freqList(ifr);
    cols{end+1} = Fx .* sin(w * t); %#ok<AGROW>
    cols{end+1} = Fx .* cos(w * t); %#ok<AGROW>
end
if isempty(cols)
    B = [];
else
    B = [cols{:}];
end
end

function [freqList, freqMeta] = get_frequency_list(cfg, highMap, templateLib, staticState)
freqMeta = struct('mode', "none", 'J1', NaN, 'J2', NaN, 'deltaJ', NaN, 'rhoJ', NaN);
if isfield(cfg, 'gapInitProjectionFreqs') && ~isempty(cfg.gapInitProjectionFreqs)
    freqList = cfg.gapInitProjectionFreqs(:);
    freqMeta.mode = "manual";
elseif isfield(cfg, 'gapInitFreqMode') && strcmpi(string(cfg.gapInitFreqMode), "vp_fixed_gap")
    [freqList, freqMeta] = estimate_projection_freqs_from_vp(cfg, highMap, templateLib, staticState);
else
    freqList = estimate_projection_freqs_from_residual(cfg, highMap, templateLib, staticState);
    freqMeta.mode = "residual_fft";
end
end

function [freqList, freqMeta] = estimate_projection_freqs_from_vp(cfg, highMap, templateLib, staticState)
t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
cfgVp = cfg;
cfgVp.numVarproCandidates = max(get_cfg_field(cfgVp, 'gapInitVpCandidateCount', 4), 4);
fit = fit_vp_main(t, V, x, templateLib, staticState.gHat, cfgVp);
freqMeta = struct('mode', "vp_fixed_gap", 'J1', fit.J1, 'J2', fit.J2, ...
    'deltaJ', fit.deltaJ, 'rhoJ', fit.rhoJ);
candidateCount = min(get_cfg_field(cfg, 'gapInitVpCandidateCount', 4), numel(fit.candidates));
freqList = zeros(0, 1);
for ic = 1:candidateCount
    freqList = [freqList; fit.candidates(ic).f(:)]; %#ok<AGROW>
end
freqList = unique(round(freqList(:), 6), 'stable');
maxFreqCount = get_cfg_field(cfg, 'gapInitAutoFreqCount', 4);
freqList = freqList(1:min(numel(freqList), maxFreqCount));
freqList = sort(freqList(:));
end

function freqList = estimate_projection_freqs_from_residual(cfg, highMap, templateLib, staticState)
t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
xq = x - staticState.dx0;
F = eval_gap_template(templateLib, staticState.gHat, xq);
Fx = eval_gap_derivative(templateLib, staticState.gHat, xq);
valid = isfinite(F) & isfinite(Fx) & isfinite(V) & abs(Fx) > 0.05 * max(abs(Fx));
if nnz(valid) < 100
    freqList = [];
    return;
end

uApprox = -(V(valid) - F(valid)) ./ Fx(valid);
tValid = t(valid);
[tUnique, ia] = unique(tValid, 'stable');
uUnique = uApprox(ia);
if numel(tUnique) < 100
    freqList = [];
    return;
end

dt = median(diff(tUnique));
if ~isfinite(dt) || dt <= 0
    freqList = [];
    return;
end
fsEff = 1 / dt;
uUnique = uUnique - mean(uUnique, 'omitnan');
uUnique = fillmissing(uUnique, 'linear', 'EndValues', 'nearest');
n = numel(uUnique);
nfft = 2^nextpow2(n);
U = abs(fft(uUnique, nfft));
freqAxis = (0:nfft-1)' * fsEff / nfft;
halfMask = freqAxis > 0 & freqAxis <= fsEff / 2;
freqAxis = freqAxis(halfMask);
U = U(halfMask);

fMin = get_cfg_field(cfg, 'gapInitFreqMin', min(cfg.f1Grid(:)));
fMax = get_cfg_field(cfg, 'gapInitFreqMax', max(cfg.f2Grid(:)));
bandMask = freqAxis >= fMin & freqAxis <= fMax;
freqBand = freqAxis(bandMask);
UBand = U(bandMask);
if isempty(freqBand)
    freqList = [];
    return;
end

peakN = get_cfg_field(cfg, 'gapInitAutoFreqCount', 4);
minSep = get_cfg_field(cfg, 'gapInitAutoFreqMinSep', 80);
[~, order] = sort(UBand, 'descend');
freqList = zeros(0, 1);
for ii = 1:numel(order)
    f = freqBand(order(ii));
    if isempty(freqList) || all(abs(freqList - f) >= minSep)
        freqList(end+1, 1) = f; %#ok<AGROW>
    end
    if numel(freqList) >= peakN
        break;
    end
end
freqList = sort(freqList(:));
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
