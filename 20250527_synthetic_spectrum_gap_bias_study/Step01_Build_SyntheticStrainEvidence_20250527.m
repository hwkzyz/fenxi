%% Step01: build synthetic strain evidence matching the Step6 spectrum.

cfg = study_config_20250527();
if ~isfile(cfg.step6StrainCsv)
    error('Missing Step6 strain reference CSV: %s', cfg.step6StrainCsv);
end

Tref = readtable(cfg.step6StrainCsv);
requiredVars = ["WindowCenterTimeBTT_s", "WaveformFreq_Hz", ...
    "StrainFreq_Hz", "StrainRawAmplitude"];
for iv = 1:numel(requiredVars)
    if ~ismember(requiredVars(iv), string(Tref.Properties.VariableNames))
        error('Step6 CSV does not contain required variable: %s', requiredVars(iv));
    end
end

valid = isfinite(Tref.WindowCenterTimeBTT_s) & ...
    isfinite(Tref.StrainRawAmplitude) & ...
    isfinite(Tref.StrainFreq_Hz);
if ~any(valid)
    error('No valid strain-reference rows found.');
end

tCenter = Tref.WindowCenterTimeBTT_s(valid);
ampMicro = Tref.StrainRawAmplitude(valid);
freqHz = Tref.StrainFreq_Hz(valid);

tStart = min(tCenter) - cfg.syntheticStrainContextSec;
tEnd = max(tCenter) + cfg.syntheticStrainContextSec;
dt = 1 / cfg.syntheticStrainFsHz;
t = (tStart:dt:tEnd).';
if t(end) < tEnd
    t(end + 1, 1) = tEnd;
end

ampEnv = interp1(tCenter, ampMicro, t, 'pchip', 'extrap');
ampEnv = max(ampEnv, 0);
freqInst = interp1(tCenter, freqHz, t, 'pchip', 'extrap');
freqInst = max(freqInst, eps);
phaseMain = 2 * pi * cumsum(freqInst) * dt;
phaseMain = phaseMain - phaseMain(1) + 0.25;

mainSignal = ampEnv .* sin(phaseMain);
f0 = median(freqHz, 'omitnan');
frot = f0 / cfg.targetEO;

% Weak deterministic components imitate the small side peaks visible in the
% experimental raw-strain spectrum without changing the dominant EO14 truth.
sideSignal = 0.035 * median(ampMicro, 'omitnan') .* sin(2*pi*(f0 - frot).*t + 1.1) + ...
    0.030 * median(ampMicro, 'omitnan') .* sin(2*pi*(f0 + frot).*t - 0.7) + ...
    0.018 * median(ampMicro, 'omitnan') .* sin(2*pi*(2*frot).*t + 0.4);

rng(cfg.syntheticStrainSeed, 'twister');
noiseSignal = cfg.syntheticStrainNoiseStdMicrostrain .* randn(size(t));
strainMicro = mainSignal + sideSignal + noiseSignal;

[freqSpectrum, ampSpectrum] = single_sided_spectrum_local(t, strainMicro);
[tStft, fStft, ampStft] = compute_stft_local(t, strainMicro, ...
    cfg.syntheticStrainFsHz, cfg);

[peakAmp, peakIdx] = max(ampSpectrum(freqSpectrum >= 500 & freqSpectrum <= 650));
freqBand = freqSpectrum(freqSpectrum >= 500 & freqSpectrum <= 650);
if isempty(freqBand)
    peakFreq = NaN;
else
    peakFreq = freqBand(peakIdx);
end

prior = struct();
prior.t_s = t;
prior.strain_microstrain = strainMicro;
prior.fs_hz = cfg.syntheticStrainFsHz;
prior.reference_table = Tref;
prior.main_freq_hz = f0;
prior.main_amp_microstrain_median = median(ampMicro, 'omitnan');
prior.main_amp_tip_mm = prior.main_amp_microstrain_median * cfg.strainToTipMmPerMicrostrain;
prior.peak_freq_hz = peakFreq;
prior.peak_amp_microstrain = peakAmp;
prior.target_eo = cfg.targetEO;
prior.nominal_rot_freq_hz = frot;
prior.strain_to_tip_mm_per_microstrain = cfg.strainToTipMmPerMicrostrain;

matPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_prior.mat');
save(matPath, 'prior');

summary = table();
summary.target_eo = cfg.targetEO;
summary.median_strain_freq_hz = f0;
summary.median_strain_amp_microstrain = prior.main_amp_microstrain_median;
summary.tip_amp_mm_from_current_step6_scale = prior.main_amp_tip_mm;
summary.synthetic_peak_freq_hz = peakFreq;
summary.synthetic_peak_amp_microstrain = peakAmp;
summary.nominal_rot_freq_hz = frot;
summary.synthetic_fs_hz = cfg.syntheticStrainFsHz;
summary.sample_count = numel(t);
writetable(summary, fullfile(cfg.outputDir, 'Step01_synthetic_strain_summary.csv'));

fig = figure('Name', 'Step01 synthetic strain evidence', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 17, 5.3], ...
    'Visible', 'off');

ax1 = subplot(1, 3, 1);
plot(ax1, t, strainMicro, 'k-', 'LineWidth', 0.35);
xlabel(ax1, 'Aligned time (s)');
ylabel(ax1, 'Raw strain (\muepsilon)');
title(ax1, 'Synthetic raw strain');
xlim(ax1, [tStart, tEnd]);
box(ax1, 'off');

ax2 = subplot(1, 3, 2);
imagesc(ax2, tStft, fStft, ampStft);
axis(ax2, 'xy');
colormap(ax2, parula);
xlabel(ax2, 'Aligned time (s)');
ylabel(ax2, 'Frequency (Hz)');
title(ax2, 'Synthetic STFT');
ylim(ax2, cfg.strainRawFreqBandHz);
cb = colorbar(ax2);
cb.Label.String = '\muepsilon';

ax3 = subplot(1, 3, 3);
plot(ax3, freqSpectrum, ampSpectrum, 'k-', 'LineWidth', 0.8);
xlabel(ax3, 'Frequency (Hz)');
ylabel(ax3, 'Amplitude (\muepsilon)');
title(ax3, 'Synthetic spectrum');
xlim(ax3, cfg.strainRawSpectrumXlimHz);
box(ax3, 'off');

pngPath = fullfile(cfg.outputDir, 'Step01_synthetic_strain_evidence.png');
exportgraphics(fig, pngPath, 'Resolution', 220);
close(fig);

fprintf('Step01 completed.\n  %s\n  %s\n', matPath, pngPath);

function [freq_hz, amp_value] = single_sided_spectrum_local(t_sec, y_value)
t_sec = t_sec(:);
y_value = y_value(:);
valid = isfinite(t_sec) & isfinite(y_value);
t_sec = t_sec(valid);
y_value = y_value(valid);
if numel(t_sec) < 4
    freq_hz = NaN;
    amp_value = NaN;
    return;
end
fs_hz = 1 / median(diff(t_sec));
y_value = y_value - mean(y_value, 'omitnan');
n = numel(y_value);
win = hann_window_local(n);
coherent_gain = max(mean(win), eps);
nfft = 2 ^ nextpow2(max(n, 4096));
y_fft = fft(y_value .* win, nfft);
amp_value = abs(y_fft(1:nfft / 2 + 1)) ./ n ./ coherent_gain .* 2;
amp_value(1) = amp_value(1) ./ 2;
freq_hz = (0:(nfft / 2)).' .* fs_hz ./ nfft;
end

function [t_stft, f_stft, amp_stft] = compute_stft_local(t_sec, y_value, fs_hz, cfg)
t_sec = t_sec(:);
y_value = y_value(:);
valid = isfinite(t_sec) & isfinite(y_value);
t_sec = t_sec(valid);
y_value = y_value(valid);
[t_sec, order] = sort(t_sec);
y_value = detrend(y_value(order));

window_len = max(16, round(cfg.strainRawStftWindowSec * fs_hz));
window_len = min(window_len, numel(y_value));
overlap_ratio = min(max(cfg.strainRawStftOverlapRatio, 0), 0.98);
hop_len = max(1, round(window_len * (1 - overlap_ratio)));
start_idx = 1:hop_len:(numel(y_value) - window_len + 1);
if isempty(start_idx)
    start_idx = 1;
end

nfft = 2 ^ nextpow2(max(window_len, 4096));
freq_all = (0:(nfft / 2)).' .* fs_hz ./ nfft;
freq_mask = freq_all >= cfg.strainRawFreqBandHz(1) & ...
    freq_all <= cfg.strainRawFreqBandHz(2);
f_stft = freq_all(freq_mask);
amp_stft = nan(numel(f_stft), numel(start_idx));
t_stft = nan(1, numel(start_idx));

win = hann_window_local(window_len);
coherent_gain = max(mean(win), eps);
for iFrame = 1:numel(start_idx)
    idx = start_idx(iFrame):(start_idx(iFrame) + window_len - 1);
    frame = y_value(idx);
    frame = frame - mean(frame, 'omitnan');
    frame_fft = fft(frame .* win, nfft);
    amp_all = abs(frame_fft(1:nfft / 2 + 1)) ./ window_len ./ coherent_gain .* 2;
    amp_all(1) = amp_all(1) ./ 2;
    amp_stft(:, iFrame) = amp_all(freq_mask);
    t_stft(iFrame) = mean(t_sec(idx), 'omitnan');
end
end

function win = hann_window_local(n)
if n <= 1
    win = ones(n, 1);
    return;
end
idx = (0:(n - 1)).';
win = 0.5 - 0.5 .* cos(2 .* pi .* idx ./ (n - 1));
end
