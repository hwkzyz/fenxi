%% Step01: Extract OPR, rotor speed and blade timing for 20250527
% This step reuses the verified legacy Step1/Step2 preprocessing.
% Output files are written by the legacy workflow under legacy/output/.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
legacyDir = fullfile(thisDir, 'legacy');
addpath(legacyDir);

cfg = Get_20250527_BTT_Config();
caseName = '20250526_2500-3500_t400';
showPlots = true;
forceRebuild = false;

caseOutputDir = fullfile(cfg.output_root, caseName);
requiredFiles = {
    fullfile(caseOutputDir, 'jiluOPR.mat')
    fullfile(caseOutputDir, 'omega.mat')
    fullfile(caseOutputDir, 'jilublade_probe1.mat')
    };

hasOutputs = all(cellfun(@isfile, requiredFiles));
if hasOutputs && ~forceRebuild
    fprintf('Step01 outputs already exist. Reusing:\n  %s\n', caseOutputDir);
else
    fprintf('Running legacy Step1/Step2 preprocessing for %s...\n', caseName);
    Step2_Extract_JiluBlade_20250527(caseName, showPlots);
end

fprintf('\nStep01 complete. Key outputs:\n');
fprintf('  %s\n', fullfile(caseOutputDir, 'jiluOPR.mat'));
fprintf('  %s\n', fullfile(caseOutputDir, 'omega.mat'));
fprintf('  %s\n', fullfile(caseOutputDir, 'jilublade_probe*.mat'));

%% Visualization: RPM and blade timing overview
load(fullfile(caseOutputDir, 'omega.mat'), 'omega_time_s', 'omega_rpm');

figure('Name', '20250527 Step01 OPR and blade timing', ...
    'Color', 'w', 'Position', [80, 80, 1250, 760], 'NumberTitle', 'off');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(omega_time_s, omega_rpm, 'k-', 'LineWidth', 1.2);
grid on; box on;
xlabel('Time (s)');
ylabel('RPM');
title(sprintf('OPR-derived speed: %s', caseName), 'Interpreter', 'none');

nexttile;
hold on; grid on; box on;
for sid = cfg.sensor_ids
    probeFile = fullfile(caseOutputDir, sprintf('jilublade_probe%d.mat', sid));
    if ~isfile(probeFile)
        continue;
    end
    S = load(probeFile, 'jilublade');
    if ~isfield(S, 'jilublade') || size(S.jilublade, 2) < 4
        continue;
    end
    tBlade = S.jilublade(:, 3);
    bladeId = S.jilublade(:, 4);
    plot(tBlade, bladeId + 0.08 * (sid - cfg.sensor_ids(1)), '.', ...
        'MarkerSize', 4, 'DisplayName', sprintf('CH%d', sid));
end
xlabel('Time (s)');
ylabel('Blade ID');
title('Extracted blade-passing sequence from legacy Step2');
legend('Location', 'best');

outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(gcf, fullfile(outDir, 'Step01_OPR_Blade_Timing_20250527.png'));
