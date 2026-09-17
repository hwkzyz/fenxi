function [scenarioList, snrList] = make_gap_init_stability_scenarios()
%make_gap_init_stability_scenarios
% Dedicated scenarios for gap-initialization stability only.

% First focus on waveform-structure effects rather than strong noise.
snrList = [20];

baseA = [0.25, 0.15];
A2OnlyList = [0, 0.03, 0.06, 0.10, 0.15, 0.20];

ampScaleList = [0, 0.5, 1.0, 2.0, 3.0, 4.0];
freqCases = {
    'freq_500_1300', [500, 1300]
    'freq_500_1330', [500, 1330]
    'freq_500_1200', [500, 1200]
    'freq_460_1180', [460, 1180]
    };
phaseCases = {
    'phase_default', [pi/4, -pi/3]
    'phase_same', [0, 0]
    'phase_opposite', [0, pi]
    'phase_quadrature', [0, pi/2]
    };
phi2FineList = linspace(0, 2*pi, 13);

rows = {};
for ia = 1:numel(ampScaleList)
    lam = ampScaleList(ia);
    Atrue = lam * baseA;
    for ifc = 1:size(freqCases, 1)
        freqName = freqCases{ifc, 1};
        ftrue = freqCases{ifc, 2};
        for ipc = 1:size(phaseCases, 1)
            phaseName = phaseCases{ipc, 1};
            phitrue = phaseCases{ipc, 2};
            sc.name = string(sprintf('gapinit_%s_%s_amp%03d', freqName, phaseName, round(100 * lam)));
            sc.group = "gap_init_stability";
            sc.A_true = Atrue;
            sc.f_true = ftrue;
            sc.phi_true = phitrue;
            sc.repeat_id = 1;
            sc.seed = 20260520 + ia * 100 + ifc * 10 + ipc;
            rows{end+1, 1} = sc; %#ok<AGROW>
        end
    end
end

% Isolate the impact of the second mode on gap initialization.
% Keep the first mode fixed and only vary A2.
for ia2 = 1:numel(A2OnlyList)
    A2 = A2OnlyList(ia2);
    sc.name = string(sprintf('gapinit_mode2only_A2_%03d', round(1000 * A2)));
    sc.group = "gap_init_mode2_only";
    sc.A_true = [0.25, A2];
    sc.f_true = [500, 1300];
    sc.phi_true = [pi/4, -pi/3];
    sc.repeat_id = 1;
    sc.seed = 20260720 + ia2;
    rows{end+1, 1} = sc; %#ok<AGROW>
end

% Fine phase scan around the frequency pair that showed the strongest
% sensitivity in the coarse screening: [500, 1200].
for ip2 = 1:numel(phi2FineList)
    phi2 = phi2FineList(ip2);
    sc.name = string(sprintf('gapinit_freq_500_1200_phasefine_phi2_%03d', round(rad2deg(phi2))));
    sc.group = "gap_init_phasefine_500_1200";
    sc.A_true = [0.25, 0.15];
    sc.f_true = [500, 1200];
    sc.phi_true = [0, phi2];
    sc.repeat_id = 1;
    sc.seed = 20260820 + ip2;
    rows{end+1, 1} = sc; %#ok<AGROW>
end

scenarioList = vertcat(rows{:});
end
