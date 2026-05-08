%% Run_All.m — Blind Vibration Parameter Identification (Complete Pipeline)
%
%  CORE IDEA:  Blind affine template reconstruction from multi-revolution
%              high-speed BTT waveform data.
%
%  REQUIRES:   ONLY high-speed multi-revolution data.
%              NO template library, NO low-speed reference, NO gap calibration.
%
%  Forward model:  y_k(x) = T(alpha_k * x - d_k)
%      alpha_k = 1 - dot_d_k / v0   (vibration velocity → spatial scaling)
%      d_k     = vibration displacement at passage centre
%
%  Algorithm:
%    Step 1 — Generate simulated multi-revolution BTT data at hidden gap
%    Step 2 — Extract individual blade-passage waveforms → common x-grid
%    Step 3 — Blind affine registration + averaging → template T(x) + {alpha_k, d_k}
%    Step 4 — Spectral analysis + harmonic fitting of {d_k, dot_d_k}
%    Step 5 — Validation against ground truth
%
%  VALIDATED:  Template shape recovery:  corr > 0.99
%              Velocity (alpha) estimation:  corr > 0.98  (low-freq case)
%              Displacement (d) estimation:   corr ≈ 0.47  (needs improvement)
%
%  NEXT STEPS: Integrate harmonic parametric model into registration loop
%              to improve d_k precision via velocity-aided estimation.

clc; clear; close all;
fprintf('=============================================================\n');
fprintf('  Blind Vibration Parameter Identification\n');
fprintf('  NO template library — ONLY multi-revolution data\n');
fprintf('=============================================================\n\n');

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

%% Config
fprintf('>>> Config ...\n');
Config;

%% Step 1: Generate data
fprintf('\n>>> Step 1: Generate simulation data ...\n');
Step1_GenerateData;

%% Step 2: Extract passages
fprintf('\n>>> Step 2: Extract blade passages ...\n');
Step2_ExtractPassages;

%% Step 3: Blind affine reconstruction
fprintf('\n>>> Step 3: Blind template + affine parameter ID ...\n');
Step3_BlindReconstruct;

%% Step 4: Vibration analysis
fprintf('\n>>> Step 4: Spectral analysis + harmonic fitting ...\n');
Step4_AnalyzeVibration;

%% Step 5: Validation
fprintf('\n>>> Step 5: Validation ...\n');
Step5_Validate;

%% Final summary
fprintf('\n=============================================================\n');
fprintf('  Pipeline complete.\n');
fprintf('  Results saved to: %s\n', fullfile(scriptDir, 'results'));
fprintf('=============================================================\n');
