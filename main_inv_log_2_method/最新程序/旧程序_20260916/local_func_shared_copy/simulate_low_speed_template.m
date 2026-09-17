function dataLow = simulate_low_speed_template(Fx, cfg, domain, noiseLevel, seed, noiseMode)
%SIMULATE_LOW_SPEED_TEMPLATE  Generate a noisy, no-vibration low-speed record.

if nargin < 5 || isempty(seed)
    seed = 1;
end
if nargin < 4 || isempty(noiseLevel)
    noiseLevel = Inf;
end
if nargin < 6 || isempty(noiseMode), noiseMode = 'snr_db'; end
if ~isfield(cfg, 'RPM_low') || isempty(cfg.RPM_low)
    cfg.RPM_low = min(cfg.RPM_high, 300);
end
if ~isfield(cfg, 'NumRevs_low') || isempty(cfg.NumRevs_low)
    cfg.NumRevs_low = max(6, cfg.NumRevs_high);
end

rng(seed, 'twister');
requestedLevel = noiseLevel;
if isinf(noiseLevel)
    noiseLevel = 0;noiseMode = 'noise_ratio';
end

dataLow = simulate_rotating_waveform_from_template(Fx, cfg.RPM_low, ...
    cfg.NumRevs_low, cfg.fs, cfg.R_tip, cfg.alpha_k, noiseLevel, domain, ...
    [], [], [], noiseMode);
dataLow.record_type = "low_speed_no_vibration";
if strcmpi(noiseMode,'snr_db'),dataLow.snr_db_requested=requestedLevel;else,dataLow.snr_db_requested=NaN;end
dataLow.random_seed = seed;
end
