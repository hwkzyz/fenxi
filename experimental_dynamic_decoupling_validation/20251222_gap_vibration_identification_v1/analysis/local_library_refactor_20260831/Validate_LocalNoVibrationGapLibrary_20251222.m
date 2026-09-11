clc;

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
correctedFile = fullfile(caseDir, 'results', 'prepared', 'calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123_nob_formal_20260830.mat');
localFile = fullfile(thisDir, ...
    'LocalNoVibrationGapLibrary_20251222_B1_S123_nob_20260831.mat');
outputFile = fullfile(thisDir, ...
    'LocalNoVibrationGapLibrary_InterpolationAudit_20251222.csv');

Sc = load(correctedFile, 'CorrectedGapLibrary');
Sl = load(localFile, 'LocalNoVibrationGapLibrary');
C = Sc.CorrectedGapLibrary;
Lall = Sl.LocalNoVibrationGapLibrary;
rows = repmat(struct(), numel(Lall.sensor), 1);

for is = 1:numel(Lall.sensor)
    L = Lall.sensor(is);
    corr = C.sensor([C.sensor.sensorId] == L.sensorId);
    xTest = linspace(min(L.xGridMm), max(L.xGridMm), 1201).';
    deltaLo = max(L.deltaGMinMm, -0.03);
    deltaHi = min(L.deltaGMaxMm, 0.03);
    deltaTest = linspace(deltaLo, deltaHi, 25);
    errors = [];
    for dg = deltaTest
        exact = exact_template_local(C.responseSurface, corr, xTest, dg);
        lookup = lookup_template_local(L, xTest, dg);
        valid = isfinite(exact) & isfinite(lookup);
        errors = [errors; lookup(valid) - exact(valid)]; %#ok<AGROW>
    end
    rows(is).sensorId = L.sensorId;
    rows(is).testPointCount = numel(errors);
    rows(is).deltaTestMinMm = deltaLo;
    rows(is).deltaTestMaxMm = deltaHi;
    rows(is).meanAbsErrorMv = mean(abs(errors), 'omitnan');
    rows(is).rmsErrorMv = sqrt(mean(errors.^2, 'omitnan'));
    rows(is).maxAbsErrorMv = max(abs(errors), [], 'omitnan');
    rows(is).zeroGapEquivalenceErrorMv = L.zeroGapEquivalenceErrorMv;
end

auditTable = struct2table(rows);
writetable(auditTable, outputFile);
disp(auditTable);
fprintf('Saved:\n  %s\n', outputFile);

function model = exact_template_local(R, corr, xq, dg)
t0 = interp1(corr.x(:), corr.vLowMv(:), xq, 'pchip', NaN);
xLib = corr.xScale .* (xq - corr.tauMm);
gBase = corr.g0Mm + corr.muGapPerXMm .* (xq - corr.tauMm);
fBase = eval_response_surface_local(R, gBase, xLib);
fDyn = eval_response_surface_local(R, gBase + dg, xLib);
model = t0 + corr.voltageGain .* (fDyn - fBase);
end

function model = lookup_template_local(L, xq, dg)
gGrid = L.deltaGGridMm(:).';
iHi = find(gGrid >= dg, 1, 'first');
iLo = find(gGrid <= dg, 1, 'last');
yLo = interp1(L.xGridMm, L.noVibrationWaveformMv(:, iLo), xq, 'pchip', NaN);
if iLo == iHi
    model = yLo;
    return;
end
yHi = interp1(L.xGridMm, L.noVibrationWaveformMv(:, iHi), xq, 'pchip', NaN);
alpha = (dg - gGrid(iLo)) / (gGrid(iHi) - gGrid(iLo));
model = (1 - alpha) .* yLo + alpha .* yHi;
end

function F = eval_response_surface_local(R, g, xq)
B0 = interp_coeff_local(R.xGrid, R.coeff(:, 1), xq);
B1 = interp_coeff_local(R.xGrid, R.coeff(:, 2), xq);
B2 = interp_coeff_local(R.xGrid, R.coeff(:, 3), xq);
F = B0 + B1 ./ g + B2 .* log(g ./ R.g0Mm);
end

function yq = interp_coeff_local(xGrid, yGrid, xq)
valid = isfinite(xGrid(:)) & isfinite(yGrid(:));
yq = interp1(xGrid(valid), yGrid(valid), xq, 'pchip', NaN);
end
