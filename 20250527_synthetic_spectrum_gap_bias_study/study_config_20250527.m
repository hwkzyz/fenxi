function cfg = study_config_20250527()
%study_config_20250527  Central config for the synthetic spectrum/bias study.

cfg = struct();
cfg.studyDir = fileparts(mfilename('fullpath'));
cfg.outputDir = fullfile(cfg.studyDir, 'output');
if exist(cfg.outputDir, 'dir') ~= 7
    mkdir(cfg.outputDir);
end

cfg.superGaussianDir = ['E:\0小论文+程序\0博士期间小论文+程序\', ...
    '7超高斯模型-权重-瞬态\程序\直叶片验证\实验验证\超高斯\20250527适配'];

cfg.resultFileS136 = fullfile(cfg.superGaussianDir, ...
    'Result_single_sync_20250527_B1_S136_20250526_2500_3500_t400_Start1p5s.mat');
cfg.step6StrainCsv = fullfile(cfg.superGaussianDir, 'output', ...
    'step6_visualization', 'Step6_B1_StrainReferenceValidation.csv');

cfg.targetBlade = 1;
cfg.targetEO = 14;
cfg.referenceSensorSets = {[1 3], [1 3 6]};
cfg.referenceSensorTags = ["S13", "S136"];

% Synthetic strain evidence.
cfg.syntheticStrainFsHz = 100000;
cfg.syntheticStrainContextSec = 0.08;
cfg.syntheticStrainNoiseStdMicrostrain = 18;
cfg.syntheticStrainSeed = 20260709;
cfg.strainRawFreqBandHz = [520, 640];
cfg.strainRawSpectrumXlimHz = [0, 1000];
cfg.strainRawStftWindowSec = 0.18;
cfg.strainRawStftOverlapRatio = 0.90;

% Step6 currently uses R = 1, so 1 microstrain maps to 0.001 mm. The FE
% ratio can be restored later without changing the synthetic voltage study.
cfg.strainToTipMmPerMicrostrain = 0.001;

% Synthetic voltage cases. Values are equivalent displacement shifts in the
% Step5 super-Gaussian coordinate, not a claim of physical gap truth.
cfg.voltageNoiseStdFractionOfRange = 0.003;
cfg.syntheticVoltageSeed = 20260710;
cfg.truePhaseRad = 0.4;
cfg.trueD0Mm = 0;
cfg.ch6ShiftCases = struct( ...
    'clean',          struct('dg6_mm', 0.00, 'dmu6_mm_per_mm', 0.000, 'noise_scale', 0.0), ...
    'noise_only',     struct('dg6_mm', 0.00, 'dmu6_mm_per_mm', 0.000, 'noise_scale', 1.0), ...
    'ch6_static_dg',  struct('dg6_mm', 0.05, 'dmu6_mm_per_mm', 0.000, 'noise_scale', 1.0), ...
    'ch6_tilt_dmu',   struct('dg6_mm', 0.00, 'dmu6_mm_per_mm', 0.035, 'noise_scale', 1.0), ...
    'ch6_dg_dmu',     struct('dg6_mm', 0.03, 'dmu6_mm_per_mm', 0.030, 'noise_scale', 1.0));

% Fit settings for the diagnostic nested models.
cfg.maxFitPoints = 6000;
cfg.fitMethods = ["fixed", "gap_only", "gap_tilt"];
cfg.fitAmpLimitMm = 0.70;
cfg.fitD0LimitMm = 0.30;
cfg.fitDgLimitMm = 0.20;
cfg.fitDmuLimitMmPerMm = 0.08;
cfg.fitDgRegWeight = 0;
cfg.fitDmuRegWeight = 0;
cfg.fitMaxIter = 260;
cfg.fitMaxFunEvals = 900;
cfg.fitSeedAmpMm = [0.25, 0.37, 0.43];
cfg.fitSeedPhiRad = [-pi/2, 0, pi/2, pi];
end
