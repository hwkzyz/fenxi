function [scenarioList, snrList] = make_multicase_scenarios_snr_method_compare()
%make_multicase_scenarios_snr_method_compare
% Fixed-amplitude/fixed-frequency scenarios for clean SNR-vs-method comparison.

snrList = [20, 15, 10, 5, 0];
numRepeat = 5;
baseSeed = 20260430;

scenarioList = repmat(struct( ...
    'name', "", ...
    'A_true', [], ...
    'f_true', [], ...
    'phi_true', [], ...
    'repeat_id', NaN, ...
    'seed', NaN), numRepeat, 1);

for i = 1:numRepeat
    scenarioList(i).name = "dual_default_fixed";
    scenarioList(i).A_true = [0.25, 0.15];
    scenarioList(i).f_true = [500, 1300];
    scenarioList(i).phi_true = [pi/4, -pi/3];
    scenarioList(i).repeat_id = i;
    scenarioList(i).seed = baseSeed + i;
end
end
