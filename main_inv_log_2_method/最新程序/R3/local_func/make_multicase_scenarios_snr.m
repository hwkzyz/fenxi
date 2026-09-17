function [scenarioList, snrList] = make_multicase_scenarios_snr()
%make_multicase_scenarios_snr  Scenario grid using direct SNR-dB definition.

[scenarioList, ~] = make_multicase_scenarios();
snrList = [40, 30, 25, 20, 15, 10];
end
