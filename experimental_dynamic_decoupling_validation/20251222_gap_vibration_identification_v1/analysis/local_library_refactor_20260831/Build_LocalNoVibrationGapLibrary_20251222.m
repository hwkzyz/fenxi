clc;

% Build a platform-specific no-vibration waveform library without changing
% the existing Step06I/Main10 production files.

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
sourceFile = fullfile(caseDir, 'results', 'prepared', 'calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123_nob_formal_20260830.mat');
outputFile = fullfile(thisDir, ...
    'LocalNoVibrationGapLibrary_20251222_B1_S123_nob_20260831.mat');
summaryFile = fullfile(thisDir, ...
    'LocalNoVibrationGapLibrary_Summary_20251222_B1_S123_nob_20260831.csv');

assert(isfile(sourceFile), 'Missing corrected low-speed library: %s', sourceFile);
S = load(sourceFile, 'CorrectedGapLibrary');
C = S.CorrectedGapLibrary;

directGapSupportMm = [0.90, 2.30];
requestedDeltaRangeMm = [-0.25, 0.25];
deltaStepMm = 5e-4;
xStepMm = 5e-3;

LocalNoVibrationGapLibrary = struct();
LocalNoVibrationGapLibrary.dataset = '20251222';
LocalNoVibrationGapLibrary.method = 'measured_low_speed_anchor_plus_frozen_static_gap_increment';
LocalNoVibrationGapLibrary.description = [ ...
    'Per-sensor measured low-speed waveform anchored at deltaG=0; ' ...
    'the frozen static response supplies only the gap-dependent increment.'];
LocalNoVibrationGapLibrary.sourceCorrectedLibraryFile = sourceFile;
LocalNoVibrationGapLibrary.targetBlade = C.targetBlade;
LocalNoVibrationGapLibrary.analysisSensors = C.analysisSensors;
LocalNoVibrationGapLibrary.directGapSupportMm = directGapSupportMm;
LocalNoVibrationGapLibrary.requestedDeltaRangeMm = requestedDeltaRangeMm;
LocalNoVibrationGapLibrary.deltaStepMm = deltaStepMm;
LocalNoVibrationGapLibrary.xStepMm = xStepMm;
LocalNoVibrationGapLibrary.formula = [ ...
    'Tnv(x,dg)=T0Measured(x)+a*[F(gRef(x)+dg,xLib(x))-F(gRef(x),xLib(x))]'];

rows = repmat(struct(), numel(C.sensor), 1);
sensorLibrary = struct([]);

for is = 1:numel(C.sensor)
    corr = C.sensor(is);
    x = corr.x(:);
    t0 = corr.vLowMv(:);
    xLib = corr.xScale .* (x - corr.tauMm);
    gRef = corr.g0Mm + corr.muGapPerXMm .* (x - corr.tauMm);

    finiteMask = isfinite(x) & isfinite(t0) & isfinite(xLib) & isfinite(gRef);
    x = x(finiteMask);
    t0 = t0(finiteMask);
    xLib = xLib(finiteMask);
    gRef = gRef(finiteMask);
    [x, order] = sort(x);
    t0 = t0(order);
    xSource = x;
    t0Source = t0;
    x = linspace(min(xSource), max(xSource), ...
        ceil((max(xSource) - min(xSource)) / xStepMm) + 1).';
    t0 = interp1(xSource, t0Source, x, 'pchip', NaN);
    xLib = corr.xScale .* (x - corr.tauMm);
    gRef = corr.g0Mm + corr.muGapPerXMm .* (x - corr.tauMm);

    deltaMin = max(directGapSupportMm(1) - gRef);
    deltaMax = min(directGapSupportMm(2) - gRef);
    deltaMin = max(deltaMin, requestedDeltaRangeMm(1));
    deltaMax = min(deltaMax, requestedDeltaRangeMm(2));
    assert(deltaMin < 0 && deltaMax > 0, ...
        'Sensor %d has no two-sided delta-gap support around zero.', corr.sensorId);

    deltaGrid = build_delta_grid_local(deltaMin, deltaMax, deltaStepMm);
    fBase = eval_response_surface_local(C.responseSurface, gRef, xLib);
    dGap = nan(numel(x), numel(deltaGrid));
    for ig = 1:numel(deltaGrid)
        gDyn = gRef + deltaGrid(ig);
        assert(all(gDyn >= directGapSupportMm(1) - 1e-12 & ...
            gDyn <= directGapSupportMm(2) + 1e-12), ...
            'Generated delta-gap state left the direct calibration support.');
        fDyn = eval_response_surface_local(C.responseSurface, gDyn, xLib);
        dGap(:, ig) = corr.voltageGain .* (fDyn - fBase);
    end
    validLibraryRow = isfinite(t0) & isfinite(gRef) & isfinite(xLib) & ...
        all(isfinite(dGap), 2);
    assert(nnz(validLibraryRow) >= 50, ...
        'Sensor %d has insufficient finite local-library support.', corr.sensorId);
    x = x(validLibraryRow);
    t0 = t0(validLibraryRow);
    xLib = xLib(validLibraryRow);
    gRef = gRef(validLibraryRow);
    dGap = dGap(validLibraryRow, :);
    tNv = t0 + dGap;

    zeroIndex = find(abs(deltaGrid) <= 10 * eps, 1, 'first');
    assert(~isempty(zeroIndex), 'deltaGGridMm must contain zero exactly.');
    zeroErrorMv = max(abs(tNv(:, zeroIndex) - t0), [], 'omitnan');

    L = struct();
    L.sensorId = corr.sensorId;
    L.xGridMm = x;
    L.T0MeasuredMv = t0;
    L.gRefPathMm = gRef;
    L.xLibPathMm = xLib;
    L.gCenterMm = corr.g0Mm;
    L.gapSlope = corr.muGapPerXMm;
    L.voltageGain = corr.voltageGain;
    L.registrationTauMm = corr.tauMm;
    L.registrationXScale = corr.xScale;
    L.registrationStatus = 'frozen_from_existing_no_b_low_speed_calibration';
    L.deltaGMinMm = deltaMin;
    L.deltaGMaxMm = deltaMax;
    L.deltaGGridMm = deltaGrid;
    L.gapIncrementWaveformMv = dGap;
    L.noVibrationWaveformMv = tNv;
    L.interpolationMethod = 'linear_delta_pchip_x';
    L.lowSpeedRmseMv = corr.lowFitRmseMv;
    L.zeroGapEquivalenceErrorMv = zeroErrorMv;
    L.gLowerMarginMm = min(gRef) - directGapSupportMm(1);
    L.gUpperMarginMm = directGapSupportMm(2) - max(gRef);
    if is == 1
        sensorLibrary = L;
    else
        sensorLibrary(is) = L;
    end

    rows(is).sensorId = corr.sensorId;
    rows(is).gCenterMm = corr.g0Mm;
    rows(is).gapSlope = corr.muGapPerXMm;
    rows(is).voltageGain = corr.voltageGain;
    rows(is).lowSpeedRmseMv = corr.lowFitRmseMv;
    rows(is).gRefMinMm = min(gRef);
    rows(is).gRefMaxMm = max(gRef);
    rows(is).deltaGMinMm = deltaMin;
    rows(is).deltaGMaxMm = deltaMax;
    rows(is).deltaGGridCount = numel(deltaGrid);
    rows(is).zeroGapEquivalenceErrorMv = zeroErrorMv;
end

LocalNoVibrationGapLibrary.sensor = sensorLibrary;
LocalNoVibrationGapLibrary.createdAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
summaryTable = struct2table(rows);
save(outputFile, 'LocalNoVibrationGapLibrary', 'summaryTable', '-v7.3');
writetable(summaryTable, summaryFile);

disp(summaryTable);
fprintf('Saved:\n  %s\n  %s\n', outputFile, summaryFile);

function grid = build_delta_grid_local(lo, hi, step)
negative = 0:-step:lo;
positive = 0:step:hi;
grid = unique([lo, fliplr(negative), positive, hi], 'sorted');
grid(abs(grid) <= 10 * eps) = 0;
end

function F = eval_response_surface_local(responseSurface, g, xq)
B0 = interp_coeff_local(responseSurface.xGrid, responseSurface.coeff(:, 1), xq);
B1 = interp_coeff_local(responseSurface.xGrid, responseSurface.coeff(:, 2), xq);
B2 = interp_coeff_local(responseSurface.xGrid, responseSurface.coeff(:, 3), xq);
F = B0 + B1 ./ g + B2 .* log(g ./ responseSurface.g0Mm);
end

function yq = interp_coeff_local(xGrid, yGrid, xq)
x = xGrid(:);
y = yGrid(:);
valid = isfinite(x) & isfinite(y);
yq = interp1(x(valid), y(valid), xq, 'pchip', NaN);
end
