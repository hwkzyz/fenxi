function [scenarioList, snrList] = make_multicase_scenarios_method_difference()
%make_multicase_scenarios_method_difference
% Targeted scenarios that are designed to expose differences between methods.

snrList = [20, 10, 5];
numRepeat = 5;
baseSeed = 20260504;

baseScenarios = [
    build_scenario("freqsep_f1300", "freq_sep", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], 0)
    build_scenario("freqsep_f1200", "freq_sep", [0.25, 0.15], [500, 1200], [pi/4, -pi/3], 0)
    build_scenario("freqsep_f1100", "freq_sep", [0.25, 0.15], [500, 1100], [pi/4, -pi/3], 0)
    build_scenario("freqsep_f1000", "freq_sep", [0.25, 0.15], [500, 1000], [pi/4, -pi/3], 0)
    build_scenario("freqsep_f0900", "freq_sep", [0.25, 0.15], [500,  900], [pi/4, -pi/3], 0)
    build_scenario("freqsep_f1330", "freq_sep", [0.25, 0.15], [500, 1330], [pi/4, -pi/3], 0)
    build_scenario("weakA2_015", "weak_mode2", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], 0)
    build_scenario("weakA2_010", "weak_mode2", [0.25, 0.10], [500, 1300], [pi/4, -pi/3], 0)
    build_scenario("weakA2_006", "weak_mode2", [0.25, 0.06], [500, 1300], [pi/4, -pi/3], 0)
    build_scenario("weakA2_003", "weak_mode2", [0.25, 0.03], [500, 1300], [pi/4, -pi/3], 0)
    build_scenario("gapbias_p000", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3],  0.00)
    build_scenario("gapbias_p002", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3],  0.02)
    build_scenario("gapbias_m002", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], -0.02)
    build_scenario("gapbias_p005", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3],  0.05)
    build_scenario("gapbias_m005", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], -0.05)
    build_scenario("gapbias_p008", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3],  0.08)
    build_scenario("gapbias_m008", "gap_bias", [0.25, 0.15], [500, 1300], [pi/4, -pi/3], -0.08)
];

scenarioList = repmat(struct( ...
    'name', "", ...
    'group', "", ...
    'A_true', [], ...
    'f_true', [], ...
    'phi_true', [], ...
    'gHat_bias_mm', 0, ...
    'repeat_id', NaN, ...
    'seed', NaN), numel(baseScenarios) * numRepeat, 1);

idx = 0;
for is = 1:numel(baseScenarios)
    for ir = 1:numRepeat
        idx = idx + 1;
        scenarioList(idx) = baseScenarios(is);
        scenarioList(idx).repeat_id = ir;
        scenarioList(idx).seed = baseSeed + 100 * is + ir;
    end
end
end

function sc = build_scenario(name, group, A_true, f_true, phi_true, gHatBias)
sc = struct( ...
    'name', string(name), ...
    'group', string(group), ...
    'A_true', A_true, ...
    'f_true', f_true, ...
    'phi_true', phi_true, ...
    'gHat_bias_mm', gHatBias, ...
    'repeat_id', NaN, ...
    'seed', NaN);
end
