function [scenarioList, noiseList] = make_multicase_scenarios()
%make_multicase_scenarios  Scenario grid used by the saved Step 6G results.

scenarioList = [
    make_scenario("dual_default", [0.25, 0.15], [500, 1300], [pi/4, -pi/3])
    make_scenario("dual_low_amp", [0.12, 0.08], [500, 1300], [pi/4, -pi/3])
    make_scenario("dual_high_amp", [0.45, 0.28], [500, 1300], [pi/4, -pi/3])
    make_scenario("dual_asymmetric", [0.40, 0.06], [500, 1300], [pi/6, -pi/5])
    make_scenario("dual_shifted_freq", [0.25, 0.15], [460, 1180], [pi/7, -pi/4])
];

noiseList = [0, 0.002, 0.005, 0.01, 0.02];
end

function sc = make_scenario(name, A_true, f_true, phi_true)
sc = struct();
sc.name = string(name);
sc.A_true = A_true(:).';
sc.f_true = f_true(:).';
sc.phi_true = phi_true(:).';
end
