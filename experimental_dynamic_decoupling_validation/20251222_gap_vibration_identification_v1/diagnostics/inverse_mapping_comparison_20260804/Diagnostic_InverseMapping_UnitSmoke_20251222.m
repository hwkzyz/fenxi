%% Diagnostic_InverseMapping_UnitSmoke_20251222
% Isolated numerical smoke test for the inverse-mapping core. This does not
% represent experimental results and does not call or modify Main10.
clear; clc;
thisDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(thisDir));
addpath(packageRoot);
addpath(fullfile(packageRoot,'functions','gap_aware'));
cfg = struct('inverseVoltageToleranceMv',1e-9, ...
    'inverseXToleranceMm',1e-10, 'inverseMaxBisectionIter',100, ...
    'inverseMinDerivativeMvPerMm',1e-5, 'weightFloor',0.05, ...
    'freqSearchHz',[300 1000]);
x = linspace(-0.3,0.3,240).';
theta = linspace(0,2*pi,240).';
uTrue = 0.018 + 0.035*sin(7*theta+0.4);
sensorIndex = ones(size(x));
% The bell-shaped response intentionally creates two valid inverse branches.
model = struct('sensor', struct('xGrid',x, 'evaluate',@synthetic_forward));
bundle = struct('X',x+uTrue, 'V',synthetic_forward(0, x), ...
    'W',ones(size(x)), 'Theta',theta, 'sensorIndex',sensorIndex, ...
    'sensorIds',1, 'rotFreqMeanHz',100);
[inv, summary] = step07jcore.inverse_map_experimental_displacement( ...
    bundle, 0, 0, model, cfg);
assert(summary.inverse_valid_fraction > 0.95, 'Too many inverse points rejected.');
assert(summary.inverse_voltage_rmse < 1e-6, 'Inverse voltage residual is too large.');
[seedTable, info] = step07jcore.solve_inverse_mapping_seed(inv, bundle, cfg);
assert(info.candidate_count >= 1, 'No EO candidates generated.');
assert(seedTable.EO(1) == 7, 'Synthetic EO candidate was not ranked first.');
fprintf('Inverse smoke PASS: valid=%.3f, voltage RMSE=%.3g mV, best EO=%d.\n', ...
    summary.inverse_valid_fraction, summary.inverse_voltage_rmse, seedTable.EO(1));

function [v, aux] = synthetic_forward(~, xq)
v = 100 - 500*xq.^2;
aux = struct('dVoltageDx', -1000*xq);
end
