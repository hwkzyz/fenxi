%% Step 4: fit the dual-frequency vibration model from recovered local shifts
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage2_segments.mat'), 'segments');
load(fullfile(outDir, 'stage3_blind_template.mat'), 'blindResult');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'truth');

tVec = [segments.tCenter]';
dVec = blindResult.d(:);
dMean = mean(dVec);
dCentered = dVec - dMean;

vibFit = fit_dualfreq_nonlinear_joint(tVec, dCentered, cfg.vibrationFit);
[fTrue, orderTrue] = sort(truth.f(:).');
ATrue = truth.A(orderTrue);
[fEst, orderEst] = sort(vibFit.f(:).');
AEst = vibFit.A(orderEst);

summary = struct();
summary.f_true = fTrue;
summary.f_est = fEst;
summary.A_true = ATrue;
summary.A_est = AEst;
summary.mean_freq_error_Hz = mean(abs(fEst - fTrue));
summary.mean_amp_error_mm = mean(abs(AEst - ATrue));
summary.shift_rmse = sqrt(mean((dCentered - vibFit.dHatCentered).^2));
summary.num_segments = numel(dVec);
summary.num_unique_shift_levels = count_unique_levels(dCentered, cfg.vibrationFit.shiftRoundTol);
summary.shaft_frequency_Hz = cfg.RPM_high / 60;
summary.true_to_shaft_ratio = truth.f(:).' / summary.shaft_frequency_Hz;
summary.true_freq_near_integer_multiple = abs(summary.true_to_shaft_ratio - round(summary.true_to_shaft_ratio)) < 1e-9;
summary.identifiability_note = build_identifiability_note(summary);

save(fullfile(outDir, 'stage4_vibration_fit.mat'), ...
    'vibFit', 'summary', '-v7.3');

fprintf('[Step 4] Blind-template vibration fit summary:\n');
disp(summary);

function vibFit = fit_dualfreq_nonlinear_joint(tVec, dVec, fitCfg)
f1Grid = fitCfg.f1Grid(:)';
f2Grid = fitCfg.f2Grid(:)';
dVec = dVec(:);

best = struct('f1', NaN, 'f2', NaN, 'betaOsc', [], 'cost', inf, 'dOsc', []);
lambda = 1e-9 * max(numel(tVec), 1);
for f1 = f1Grid
    s1 = sin(2*pi*f1*tVec);
    c1 = cos(2*pi*f1*tVec);
    for f2 = f2Grid
        if abs(f1 - f2) < fitCfg.minFreqSeparation
            continue;
        end
        s2 = sin(2*pi*f2*tVec);
        c2 = cos(2*pi*f2*tVec);
        H = [s1, c1, s2, c2];
        betaOsc = (H' * H + lambda * eye(4)) \ (H' * dVec);
        dOsc = H * betaOsc;
        cost = sum((dVec - dOsc).^2);
        if cost < best.cost
            best.f1 = f1;
            best.f2 = f2;
            best.betaOsc = betaOsc;
            best.cost = cost;
            best.dOsc = dOsc;
        end
    end
end

if ~isfinite(best.cost)
    error('Dual-frequency fit failed to find a valid coarse-grid solution.');
end

b = best.betaOsc;
A0 = [hypot(b(1), b(2)), atan2(b(2), b(1)), best.f1, ...
    hypot(b(3), b(4)), atan2(b(4), b(3)), best.f2];

lb = [0, -pi, min(f1Grid), 0, -pi, min(f2Grid)];
ub = [fitCfg.ampUpperBound, pi, max(f1Grid), fitCfg.ampUpperBound, pi, max(f2Grid)];
resFun = @(p) dVec - displacement_2freq(p, tVec);
pOpt = refine_bounded_dualfreq(resFun, A0, lb, ub, fitCfg.maxIter);

if pOpt(3) > pOpt(6)
    pOpt = [pOpt(4), pOpt(5), pOpt(6), pOpt(1), pOpt(2), pOpt(3)];
end

dHatCentered = displacement_2freq(pOpt, tVec);

vibFit = struct();
vibFit.f = [pOpt(3), pOpt(6)];
vibFit.A = [pOpt(1), pOpt(4)];
vibFit.phi = [pOpt(2), pOpt(5)];
vibFit.p0 = A0(:).';
vibFit.p = pOpt(:).';
vibFit.betaOsc = best.betaOsc;
vibFit.coarseFreq = [best.f1, best.f2];
vibFit.coarseCost = best.cost;
vibFit.cost = sum((dVec - dHatCentered).^2);
vibFit.dHatCentered = dHatCentered(:);
vibFit.dHat = dHatCentered(:);
end

function pOpt = refine_bounded_dualfreq(resFun, p0, lb, ub, maxIter)
p0 = min(max(p0, lb), ub);
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'MaxIterations', maxIter);
    pOpt = lsqnonlin(resFun, p0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((p0 - lb + 1e-9) ./ max(ub - p0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    zOpt = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', maxIter * 5, 'TolX', 1e-9));
    pOpt = lb + width ./ (1 + exp(-zOpt));
end
pOpt(2) = wrap_to_pi_local(pOpt(2));
pOpt(5) = wrap_to_pi_local(pOpt(5));
end

function u = displacement_2freq(p, t)
u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end

function a = wrap_to_pi_local(a)
a = mod(a + pi, 2*pi) - pi;
end

function nLevel = count_unique_levels(x, tol)
key = round(x(:) / max(tol, eps)) * max(tol, eps);
nLevel = numel(unique(key));
end

function note = build_identifiability_note(summary)
if all(summary.true_freq_near_integer_multiple) && summary.num_unique_shift_levels <= 4
    note = sprintf(['The blind-shift sequence has only %d unique levels, and the truth frequencies are nearly integer ', ...
        'multiples of the shaft frequency %.3f Hz. With one shift sample per passing, this creates repeated phase sampling ', ...
        'across revolutions, so dual-frequency recovery from d_q alone becomes weakly identifiable.'], ...
        summary.num_unique_shift_levels, summary.shaft_frequency_Hz);
elseif summary.num_unique_shift_levels <= 4
    note = sprintf(['The blind-shift sequence has only %d unique levels. Even though the template reconstruction converged, ', ...
        'the vibration fitting stage may still be under-informed because one shift is extracted from each passing segment.'], ...
        summary.num_unique_shift_levels);
else
    note = 'No obvious repeated-level identifiability issue was detected in the recovered blind-shift sequence.';
end
end
