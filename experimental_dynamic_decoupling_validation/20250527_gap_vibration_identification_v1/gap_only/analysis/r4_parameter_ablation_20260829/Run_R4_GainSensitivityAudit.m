function T = Run_R4_GainSensitivityAudit()
% Quantify the deterministic effect of replacing fitted a_s by a_s = 1.
% This is a library-level audit; it does not rerun high-speed identification.
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
libFile = fullfile(root, 'inputs', 'calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136.mat');
if ~isfile(libFile)
    error('Missing corrected gap library: %s', libFile);
end
S = load(libFile, 'CorrectedGapLibrary');
L = S.CorrectedGapLibrary;
R = L.responseSurface;
g = R.trueGapMm(:);
x = R.xGrid(:);
rows = cell(numel(L.sensor), 1);
for i = 1:numel(L.sensor)
    C = L.sensor(i);
    % The formal increment is a_s*dF; the a=1 increment is dF.
    [~, Fbase] = rawPath(R, C, x, C.g0Mm);
    dFormal = nan(numel(g), numel(x));
    for j = 1:numel(g)
        [~, Fdyn] = rawPath(R, C, x, g(j));
        dFormal(j,:) = C.voltageGain .* (Fdyn - Fbase);
    end
    dUnit = dFormal ./ C.voltageGain;
    valid = isfinite(dFormal) & isfinite(dUnit);
    absFormal = abs(dFormal(valid));
    absError = abs(dUnit(valid) - dFormal(valid));
    nonzero = valid & abs(dFormal) > 1e-9;
    ratio = abs(dUnit(nonzero) ./ dFormal(nonzero));
    rows{i} = {C.sensorId, C.voltageGain, 100*(C.voltageGain-1), ...
        100*(1-1/C.voltageGain), min(absFormal), max(absFormal), ...
        median(absError), max(absError), min(ratio), max(ratio), nnz(valid)};
end
T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'sensorId','aFitted','aMinus1_percent','incrementUnderestimate_percent', ...
     'formalIncrementAbsMin_mV','formalIncrementAbsMax_mV', ...
     'unitGainAbsErrorMedian_mV','unitGainAbsErrorMax_mV', ...
     'unitToFormalRatioMin','unitToFormalRatioMax','validGridCount'});
out = fullfile(folder, 'formal_vs_a1_gain_sensitivity.csv');
writetable(T, out);
save(fullfile(folder, 'formal_vs_a1_gain_sensitivity.mat'), 'T', 'libFile');
disp(T);
fprintf('Wrote %s\n', out);
end

function [gEff, F] = rawPath(R, C, x, gNominal)
gEff = gNominal + C.muGapPerXMm .* (x - C.tauMm);
xLib = C.xScale .* (x - C.tauMm);
B0 = interp1(R.xGrid, R.coeff(:,1), xLib, 'linear', NaN);
B1 = interp1(R.xGrid, R.coeff(:,2), xLib, 'linear', NaN);
B2 = interp1(R.xGrid, R.coeff(:,3), xLib, 'linear', NaN);
F = B0 + B1 ./ gEff + B2 .* log(gEff ./ R.g0Mm);
end
