function result = Step3_Align_StrainSpectrum_With_RPM_20250527(case_name, base_alignment_offset_sec, order_candidates, freq_plot_hz, show_plots)
%STEP3_ALIGN_STRAINSPECTRUM_WITH_RPM_20250527
% Extract the synchronous line directly from the strain STFT, then align
% the BTT RPM curve to that line in time.
%
% Main idea:
%   1. Build full-record strain STFT
%   2. Estimate the strain-side rotational frequency f_rot_strain(t) using
%      harmonic-comb matching across the whole spectrum
%   3. Track a continuous f_rot_strain(t) with dynamic programming
%   4. Scan time shift tau so that RPM_BTT(t-tau)/60 matches f_rot_strain(t)
%   5. Choose the visible black line order from STFT energy along k*f_rot
%   6. Plot STFT + extracted black line + shifted BTT order line

cfg = Get_20250527_BTT_Config();
if nargin < 1 || isempty(case_name)
    case_name = cfg.dynamic_cases{1};
end
if nargin < 2 || isempty(base_alignment_offset_sec)
    base_alignment_offset_sec = cfg.step3_default_alignment_offset_sec;
end
if nargin < 3 || isempty(order_candidates)
    order_candidates = cfg.step3_align_order_candidates;
end
if nargin < 4 || isempty(freq_plot_hz)
    freq_plot_hz = cfg.step3_align_freq_plot_hz;
end
if nargin < 5
    show_plots = true;
end

fprintf('>>> [Step3Align][%s] Preparing Step3 products...\n', case_name);
step3_result = Step3_Displacement_Calculation_20250527(case_name, base_alignment_offset_sec, false);

rpm_time_s = step3_result.rpm_time_s(:);
rpm_values = step3_result.rpm_values(:);
if isempty(rpm_time_s) || isempty(rpm_values)
    error('No RPM trace available from Step3 for case %s.', case_name);
end

[t_strain_raw, v_strain, Fs_strain] = load_strain_signal_local(step3_result.strain_file);
t_strain_aligned = t_strain_raw + base_alignment_offset_sec;

fprintf('>>> [Step3Align][%s] Strain file: %s\n', case_name, step3_result.strain_file);
fprintf('>>> [Step3Align][%s] RPM points: %d | order candidates: [%s]\n', ...
    case_name, numel(rpm_values), num2str(order_candidates));

[stft_time_s, stft_freq_hz, stft_amp] = compute_strain_stft_local(v_strain, Fs_strain, ...
    t_strain_aligned(1), cfg.step3_align_stft_window_sec, cfg.step3_align_stft_overlap_ratio);
stft_time_s = stft_time_s(:);
stft_freq_hz = stft_freq_hz(:);

frot_band_hz = [max(5, min(rpm_values / 60) - cfg.step3_align_frot_margin_hz), ...
    max(rpm_values / 60) + cfg.step3_align_frot_margin_hz];
fprintf('>>> [Step3Align][%s] Strain-side rotational-frequency search band: %.2f to %.2f Hz\n', ...
    case_name, frot_band_hz(1), frot_band_hz(2));

[frot_grid_hz, comb_score] = build_harmonic_comb_score_local( ...
    stft_freq_hz, stft_amp, frot_band_hz, cfg.step3_align_frot_grid_step_hz, ...
    freq_plot_hz(2), cfg.step3_align_frot_min_harmonics);
[frot_strain_hz, frot_confidence] = track_synchronous_frot_local( ...
    frot_grid_hz, comb_score, cfg.step3_align_frot_transition_weight);

tau_candidates = build_tau_candidates_local(stft_time_s, rpm_time_s, cfg.step3_align_tau_step_sec);
[tau_rmse_hz, tau_overlap_count] = scan_tau_alignment_local( ...
    stft_time_s, frot_strain_hz, frot_confidence, rpm_time_s, rpm_values / 60, ...
    tau_candidates, cfg.step3_align_min_valid_points, case_name);
[best_tau_sec, best_tau_rmse_hz, best_tau_overlap_mask] = pick_best_tau_local( ...
    stft_time_s, frot_strain_hz, rpm_time_s, rpm_values / 60, tau_candidates, tau_rmse_hz);

[best_order, order_energy_scores] = pick_visible_order_local( ...
    stft_time_s, frot_strain_hz, stft_freq_hz, stft_amp, order_candidates, freq_plot_hz);

black_line_hz = best_order * frot_strain_hz;
shifted_rpm_time_s = rpm_time_s + best_tau_sec;
shifted_btt_black_hz = best_order * (rpm_values / 60);
btt_shifted_span = [min(shifted_rpm_time_s), max(shifted_rpm_time_s)];

result = struct();
result.case_name = case_name;
result.base_alignment_offset_sec = base_alignment_offset_sec;
result.strain_file = step3_result.strain_file;
result.rpm_time_s = rpm_time_s;
result.rpm_values = rpm_values;
result.stft_time_s = stft_time_s;
result.stft_freq_hz = stft_freq_hz;
result.stft_amp = stft_amp;
result.frot_band_hz = frot_band_hz;
result.frot_grid_hz = frot_grid_hz;
result.comb_score = comb_score;
result.frot_strain_hz = frot_strain_hz;
result.frot_confidence = frot_confidence;
result.tau_candidates = tau_candidates;
result.tau_rmse_hz = tau_rmse_hz;
result.tau_overlap_count = tau_overlap_count;
result.best_tau_sec = best_tau_sec;
result.best_tau_rmse_hz = best_tau_rmse_hz;
result.best_tau_overlap_mask = best_tau_overlap_mask;
result.order_candidates = order_candidates;
result.order_energy_scores = order_energy_scores;
result.best_order = best_order;
result.black_line_hz = black_line_hz;
result.shifted_rpm_time_s = shifted_rpm_time_s;
result.shifted_btt_black_hz = shifted_btt_black_hz;
result.btt_shifted_span = btt_shifted_span;

save(fullfile(step3_result.btt_output_dir, 'Step3_Spectrum_RPM_Alignment_20250527.mat'), 'result');

fprintf('>>> [Step3Align][%s] Best tau = %.3f s | f_rot RMSE = %.3f Hz\n', ...
    case_name, best_tau_sec, best_tau_rmse_hz);
fprintf('>>> [Step3Align][%s] Best visible black-line order = EO%d\n', case_name, best_order);
fprintf('>>> [Step3Align][%s] Shifted BTT window: %.3f s to %.3f s\n', ...
    case_name, btt_shifted_span(1), btt_shifted_span(2));

if show_plots
    plot_alignment_results_local(result, t_strain_aligned, v_strain, freq_plot_hz);
end
end


function [t_strain, v_strain, Fs_strain] = load_strain_signal_local(strain_file)
loaded = load(strain_file);
if ~isfield(loaded, 'Datas')
    error('Strain file %s does not contain Datas.', strain_file);
end

Datas = loaded.Datas;
t_strain = Datas(:, 1);
v_strain = Datas(:, 2);
Fs_strain = 1 / median(diff(t_strain));
end


function [t_stft, f_axis, amp] = compute_strain_stft_local(v_strain, Fs_strain, t_start_aligned, win_sec, overlap_ratio)
win_len = max(256, floor(win_sec * Fs_strain));
noverlap = floor(win_len * overlap_ratio);
nfft = max(1024, 2 ^ nextpow2(win_len));
window = hann(win_len, 'periodic');

[S, f_axis, t_local] = spectrogram(v_strain - mean(v_strain), window, noverlap, nfft, Fs_strain);
t_stft = t_local + t_start_aligned;
amp = abs(S);
end


function [frot_grid_hz, comb_score] = build_harmonic_comb_score_local(f_axis, amp, frot_band_hz, frot_step_hz, max_freq_hz, min_harmonics)
frot_grid_hz = (frot_band_hz(1):frot_step_hz:frot_band_hz(2)).';
num_grid = numel(frot_grid_hz);
num_time = size(amp, 2);
comb_score = -inf(num_grid, num_time);

amp_db = 20 * log10(amp + eps);
for it = 1:num_time
    col = amp_db(:, it);
    col = col - min(col);
    max_col = max(col);
    if max_col > eps
        col = col / max_col;
    end

    for ig = 1:num_grid
        f0 = frot_grid_hz(ig);
        harmonics = f0 * (1:floor(max_freq_hz / f0));
        if numel(harmonics) < min_harmonics
            continue;
        end

        vals = interp1(f_axis, col, harmonics, 'linear', 0);
        weights = 1 ./ sqrt(1:numel(harmonics));
        comb_score(ig, it) = sum(weights(:) .* vals(:)) / sum(weights);
    end
end
end


function [frot_strain_hz, frot_confidence] = track_synchronous_frot_local(frot_grid_hz, comb_score, transition_weight)
num_grid = numel(frot_grid_hz);
num_time = size(comb_score, 2);

dp = -inf(num_grid, num_time);
prev = zeros(num_grid, num_time);
dp(:, 1) = comb_score(:, 1);

freq_diff = abs(frot_grid_hz - frot_grid_hz.');
trans_pen = transition_weight * (freq_diff / max(mean(diff(frot_grid_hz)), eps));

for t = 2:num_time
    for i = 1:num_grid
        prev_vals = dp(:, t - 1) - trans_pen(:, i);
        [best_prev_val, best_prev_idx] = max(prev_vals);
        dp(i, t) = comb_score(i, t) + best_prev_val;
        prev(i, t) = best_prev_idx;
    end
end

[~, end_idx] = max(dp(:, end));
path_idx = zeros(num_time, 1);
path_idx(end) = end_idx;
for t = num_time:-1:2
    path_idx(t - 1) = prev(path_idx(t), t);
end

frot_strain_hz = frot_grid_hz(path_idx);
frot_confidence = zeros(num_time, 1);
for t = 1:num_time
    frot_confidence(t) = comb_score(path_idx(t), t);
end
frot_strain_hz = smoothdata(frot_strain_hz, 'movmean', 5);
frot_confidence = smoothdata(frot_confidence, 'movmean', 5);
end


function tau_candidates = build_tau_candidates_local(stft_time_s, rpm_time_s, tau_step_sec)
tau_min = stft_time_s(1) - rpm_time_s(1);
tau_max = stft_time_s(end) - rpm_time_s(end);
if tau_max < tau_min
    error('The BTT RPM segment is longer than the available strain timeline.');
end
tau_candidates = tau_min:tau_step_sec:tau_max;
end


function [tau_rmse_hz, tau_overlap_count] = scan_tau_alignment_local(stft_time_s, frot_strain_hz, frot_confidence, rpm_time_s, btt_frot_hz, tau_candidates, min_valid_points, case_name)
num_taus = numel(tau_candidates);
tau_rmse_hz = nan(num_taus, 1);
tau_overlap_count = zeros(num_taus, 1);

for iTau = 1:num_taus
    tau = tau_candidates(iTau);
    if mod(iTau, max(1, round(num_taus / 10))) == 1 || iTau == num_taus
        fprintf('>>> [Step3Align][%s] Scanning tau %.3f s (%d/%d)...\n', case_name, tau, iTau, num_taus);
    end

    btt_interp = interp1(rpm_time_s + tau, btt_frot_hz, stft_time_s, 'linear', nan);
    mask = isfinite(btt_interp);
    tau_overlap_count(iTau) = nnz(mask);
    if tau_overlap_count(iTau) < min_valid_points
        continue;
    end

    err_hz = frot_strain_hz(mask) - btt_interp(mask);
    w = frot_confidence(mask) - min(frot_confidence(mask));
    if sum(w) < eps
        w = ones(size(err_hz));
    end
    tau_rmse_hz(iTau) = sqrt(sum(w .* err_hz .^ 2) / sum(w));
end
end


function [best_tau_sec, best_tau_rmse_hz, best_overlap_mask] = pick_best_tau_local(stft_time_s, frot_strain_hz, rpm_time_s, btt_frot_hz, tau_candidates, tau_rmse_hz)
[best_tau_rmse_hz, idx] = min(tau_rmse_hz);
if ~isfinite(best_tau_rmse_hz)
    error('No valid tau solution found when aligning BTT RPM to the strain-side synchronous line.');
end

best_tau_sec = tau_candidates(idx);
btt_interp = interp1(rpm_time_s + best_tau_sec, btt_frot_hz, stft_time_s, 'linear', nan);
best_overlap_mask = isfinite(btt_interp) & isfinite(frot_strain_hz);
end


function [best_order, order_energy_scores] = pick_visible_order_local(stft_time_s, frot_strain_hz, f_axis, amp, order_candidates, freq_plot_hz)
num_orders = numel(order_candidates);
order_energy_scores = nan(num_orders, 1);
amp_db = 20 * log10(amp + eps);

for i = 1:num_orders
    k = order_candidates(i);
    line_hz = k * frot_strain_hz;
    mask = isfinite(line_hz) & line_hz >= freq_plot_hz(1) & line_hz <= freq_plot_hz(2);
    if nnz(mask) < 20
        continue;
    end
    vals = interp2(stft_time_s(:).', f_axis(:), amp_db, stft_time_s(mask).', line_hz(mask).', 'linear', nan);
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    order_energy_scores(i) = mean(vals);
end

[~, idx] = max(order_energy_scores);
best_order = order_candidates(idx);
end


function plot_alignment_results_local(result, t_strain_aligned, v_strain, freq_plot_hz)
log_amp = 20 * log10(result.stft_amp + eps);
full_xlim = [t_strain_aligned(1), t_strain_aligned(end)];

fig1 = figure('Name', sprintf('20250527 Spectrum-RPM Alignment - %s', result.case_name), ...
    'Color', 'w', 'Position', [70, 50, 1500, 980], 'NumberTitle', 'off');
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(tl, [2, 1]);
imagesc(result.stft_time_s, result.stft_freq_hz, log_amp);
axis xy;
hold on;
grid on;
box on;
colormap(turbo);
colorbar;
mask_black = isfinite(result.black_line_hz) & result.black_line_hz >= freq_plot_hz(1) & result.black_line_hz <= freq_plot_hz(2);
plot(result.stft_time_s(mask_black), result.black_line_hz(mask_black), 'k-', 'LineWidth', 2.4, 'DisplayName', 'Extracted black line');
mask_btt = result.shifted_btt_black_hz >= freq_plot_hz(1) & result.shifted_btt_black_hz <= freq_plot_hz(2);
plot(result.shifted_rpm_time_s(mask_btt), result.shifted_btt_black_hz(mask_btt), '-', ...
    'Color', [1.0 0.0 1.0], 'LineWidth', 1.8, 'DisplayName', sprintf('Shifted BTT EO%d', result.best_order));
title(sprintf('Strain STFT + extracted black line | EO=%d | tau=%.3f s | f_{rot} RMSE=%.2f Hz', ...
    result.best_order, result.best_tau_sec, result.best_tau_rmse_hz));
xlabel('Aligned time (s)');
ylabel('Frequency (Hz)');
ylim(freq_plot_hz);
xlim(full_xlim);
legend('Location', 'northeast');

nexttile;
plot(result.tau_candidates, result.tau_rmse_hz, 'LineWidth', 1.6, 'Color', [0 0.45 0.74]);
hold on;
grid on;
box on;
plot(result.best_tau_sec, result.best_tau_rmse_hz, 'rp', 'MarkerFaceColor', 'y', 'MarkerSize', 13);
xlabel('Additional time shift applied to BTT RPM (s)');
ylabel('RMSE of f_{rot} match (Hz)');
title('Tau scan using strain-side synchronous line');

nexttile;
hold on;
grid on;
box on;
yyaxis left;
stride = max(1, floor(numel(t_strain_aligned) / 40000));
plot(t_strain_aligned(1:stride:end), v_strain(1:stride:end), '-', ...
    'Color', [0.20 0.60 0.20 0.35], 'LineWidth', 0.5, 'DisplayName', 'Aligned strain');
ylabel('Strain (V)');
yl = ylim;
plot([result.btt_shifted_span(1), result.btt_shifted_span(1)], yl, '--', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.2, 'DisplayName', 'Shifted BTT start');
plot([result.btt_shifted_span(2), result.btt_shifted_span(2)], yl, '-.', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'DisplayName', 'Shifted BTT end');
text(mean(result.btt_shifted_span), yl(2) - 0.08 * (yl(2) - yl(1)), ...
    sprintf('Best shifted BTT window: %.2f s to %.2f s', result.btt_shifted_span(1), result.btt_shifted_span(2)), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'BackgroundColor', 'w');

yyaxis right;
plot(result.shifted_rpm_time_s, result.rpm_values, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Shifted BTT RPM');
ylabel('RPM');
xlabel('Aligned time (s)');
title('Full-record strain + shifted BTT RPM');
xlim(full_xlim);
legend('Location', 'best');

fig2 = figure('Name', sprintf('20250527 Black-Line Extraction Detail - %s', result.case_name), ...
    'Color', 'w', 'Position', [120, 80, 1400, 760], 'NumberTitle', 'off');
tiledlayout(fig2, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
grid on;
box on;
plot(result.stft_time_s, result.frot_strain_hz, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Estimated strain-side f_{rot}');
plot(result.stft_time_s(result.best_tau_overlap_mask), ...
    interp1(result.shifted_rpm_time_s, result.rpm_values / 60, result.stft_time_s(result.best_tau_overlap_mask), 'linear', nan), ...
    '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'DisplayName', 'Shifted BTT f_{rot}');
xlabel('Aligned time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Synchronous base frequency alignment | tau = %.3f s', result.best_tau_sec));
xlim(full_xlim);
legend('Location', 'best');

nexttile;
bar(result.order_candidates, result.order_energy_scores, 'FaceColor', [0.30 0.60 0.85]);
hold on;
grid on;
box on;
plot(result.best_order, result.order_energy_scores(result.order_candidates == result.best_order), 'rp', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 13);
xlabel('Order candidate');
ylabel('Mean STFT energy along k*f_{rot,strain}');
title(sprintf('Visible-order selection for black line | best = EO%d', result.best_order));
end
