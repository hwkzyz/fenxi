%% Config.m — Blind Vibration Identification Configuration
%  This method uses ONLY high-speed multi-revolution waveform data.
%  No template library, no low-speed reference, no gap calibration.
%
%  Forward model per passage:
%    y_k(x) = T(alpha_k * x - d_k)
%    alpha_k = 1 - dot_d_k / v0   (vibration velocity → spatial scaling)
%    d_k     = vibration displacement at passage centre
%
%  Gauge conditions:  mean(alpha_k)=1,  mean(d_k)=0  →  unique solution

clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir  = fileparts(scriptDir);

outDir = fullfile(scriptDir, 'results');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- Data source ------------------------------------------------
cfg.dataFile = fullfile(rootDir, '间隙的影响', '直叶片2mm_不同间隙.txt');

%% ---- Operating point --------------------------------------------
cfg.RPM_high    = 3000;         % high-speed rotation (rpm)
cfg.R_tip       = 62;           % blade tip radius (mm)
cfg.fs          = 2e5;          % sampling frequency (Hz)
cfg.NumRevs     = 100;          % number of revolutions (2 s of data)
cfg.alpha_k     = deg2rad([55, 116.76, 178.52]);  % sensor angular positions

%% ---- True vibration (hidden from the algorithm) -----------------
%  Low-frequency test case: f << 75 Hz Nyquist to validate core algorithm.
%  Once validated, extension to high freq via parametric methods follows.
cfg.A_true   = [0.35, 0.20];    % amplitude (mm)
cfg.f_true   = [12, 31];        % frequency (Hz) — well below sensor Nyquist
cfg.phi_true = [pi/4, -pi/3];   % phase (rad)

%% ---- Hidden gap (the gap used to generate high-speed data) ------
cfg.g_hidden = 1.0;             % true high-speed gap, NOT known to algorithm

%% ---- Noise level ------------------------------------------------
cfg.noiseRatio = 0.002;         % relative to waveform span

%% ---- Passage extraction -----------------------------------------
cfg.activeLevel  = 0.08;        % fraction of span above baseline → active region
cfg.xGridN       = 501;         % common spatial grid points per passage

%% ---- Blind reconstruction ---------------------------------------
cfg.blind.maxIter       = 30;       % max affine-registration / averaging cycles
cfg.blind.tolTemplate   = 1e-5;     % relative template change for convergence
cfg.blind.tolParam      = 1e-4;     % mean |delta alpha| for convergence
cfg.blind.alphaSearch   = [0.85, 1.15];  % alpha search range for registration
cfg.blind.dSearch       = [-1.0, 1.0];   % d search range (mm)
cfg.blind.alphaGridN    = 61;       % coarse grid points for alpha
cfg.blind.dGridN        = 81;       % coarse grid points for d
cfg.blind.smoothLambda  = 0.01;     % template smoothing (ridge penalty on 2nd deriv)

%% ---- Vibration analysis -----------------------------------------
cfg.vib.maxModes     = 5;           % max number of modes to consider
cfg.vib.freqSearchHz = [5, 60];     % frequency search range for peak picking

%% ---- Font / plot defaults ---------------------------------------
fontName = 'Times New Roman';
set(groot, 'DefaultAxesFontName', fontName, ...
    'DefaultTextFontName', fontName, ...
    'DefaultLegendFontName', fontName, ...
    'DefaultAxesFontSize', 9, ...
    'DefaultTextFontSize', 9, ...
    'DefaultAxesLineWidth', 0.8, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultAxesTickDir', 'in', ...
    'DefaultAxesBox', 'on');

save(fullfile(outDir, 'cfg.mat'), 'cfg');
fprintf('[Config] Saved to %s\n', fullfile(outDir, 'cfg.mat'));
fprintf('[Config] RPM = %d,  Revs = %d,  hidden gap = %.1f mm\n', ...
    cfg.RPM_high, cfg.NumRevs, cfg.g_hidden);
fprintf('[Config] Sensors: %s deg\n', mat2str(rad2deg(cfg.alpha_k), 3));
fprintf('[Config] True vibration: f=[%.0f, %.0f] Hz, A=[%.2f, %.2f] mm\n', ...
    cfg.f_true(1), cfg.f_true(2), cfg.A_true(1), cfg.A_true(2));
