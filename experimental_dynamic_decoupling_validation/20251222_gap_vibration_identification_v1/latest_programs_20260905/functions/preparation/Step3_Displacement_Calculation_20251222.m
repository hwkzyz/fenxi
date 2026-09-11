function result = Step3_Displacement_Calculation_20251222(case_name, alignment_offset_sec, show_plots)
%STEP3_DISPLACEMENT_CALCULATION_20251222 Step3-style synchronized display for 20251222.
% This script follows the 20251222 displacement logic more closely:
%   1. Read Step2 jilublade_probe*.mat arrival-time results
%   2. Use Step1 Standard_Relative_Angles + OPR-integrated relative angle
%      to convert each arrival into vibration displacement (mm)
%   3. Save compatibility files jilublade_probe*_vib_final.mat
%   4. Visualize synchronized BTT displacement / strain, STFT, and local FFT

cfg = Get_20251222_BTT_Config();
if nargin < 1 || isempty(case_name)
    case_name = cfg.dynamic_cases{1};
end
if nargin < 2 || isempty(alignment_offset_sec)
    alignment_offset_sec = cfg.step3_default_alignment_offset_sec;
end
if nargin < 3
    show_plots = true;
end

btt_output_dir = fullfile(cfg.output_root, case_name);
if ~isfolder(btt_output_dir)
    error('Step3 output directory not found: %s', btt_output_dir);
end

strain_case_dir = map_strain_case_dir_local(case_name, cfg.strain_root);
if ~isfolder(strain_case_dir)
    error('Strain case directory not found: %s', strain_case_dir);
end

fprintf('>>> [Step3][%s] BTT output dir: %s\n', case_name, btt_output_dir);
fprintf('>>> [Step3][%s] Strain dir: %s\n', case_name, strain_case_dir);
fprintf('>>> [Step3][%s] Alignment offset applied to strain = %.4f s\n', case_name, alignment_offset_sec);
fprintf('>>> [Step3][%s] Input convention: Step2 jilublade_probe*.mat -> Step3 displacement conversion.\n', case_name);
[btt_time_shift_sec, alignment_meta] = load_step3_alignment_shift_local(btt_output_dir, case_name);
if alignment_meta.used_alignment_file
    fprintf('>>> [Step3][%s] Using Step3Align best_tau = %.4f s to reverse-shift strain timeline.\n', ...
        case_name, btt_time_shift_sec);
else
    fprintf('>>> [Step3][%s] No Step3Align result found. Reverse strain shift = %.4f s.\n', ...
        case_name, btt_time_shift_sec);
end

sensor_config_path = fullfile(cfg.reference_output_dir, 'Sensor_Config_20251222.mat');
if ~isfile(sensor_config_path)
    error('Step3 requires Sensor_Config_20251222.mat. Missing file: %s', sensor_config_path);
end
loaded_cfg = load(sensor_config_path, 'Sensor_Config');
Sensor_Config = loaded_cfg.Sensor_Config;

btt_data = cell(1, numel(cfg.step3_default_btt_sensor_ids));
loaded_sensor_ids = [];
load_modes = strings(1, numel(cfg.step3_default_btt_sensor_ids));
for i = 1:numel(cfg.step3_default_btt_sensor_ids)
    sid = cfg.step3_default_btt_sensor_ids(i);
    btt_data{i} = load_btt_display_data_local(btt_output_dir, sid, Sensor_Config, cfg, case_name);
    if ~isempty(btt_data{i})
        loaded_sensor_ids(end + 1) = sid; %#ok<AGROW>
        load_modes(i) = string(btt_data{i}.mode);
    end
end

if isempty(loaded_sensor_ids)
    error('No Step3-compatible BTT files found under %s', btt_output_dir);
end

strain_file = pick_latest_strain_file_local(strain_case_dir, cfg.step3_default_strain_channel);
[t_strain_raw, v_strain_raw, Fs_strain] = load_strain_signal_local(strain_file);
[v_strain, despike_info] = preprocess_strain_signal_local(v_strain_raw, Fs_strain, cfg);
t_strain_aligned = t_strain_raw + alignment_offset_sec - btt_time_shift_sec;
[rpm_time_s, rpm_values] = load_rpm_trace_local(btt_output_dir, cfg);

fprintf('>>> [Step3][%s] Loaded strain file: %s\n', case_name, strain_file);
if despike_info.enabled
    fprintf('>>> [Step3][%s] Strain despike copy: replaced %d spikes | window=%d samples | sigma=%.1f\n', ...
        case_name, despike_info.num_replaced, despike_info.window_samples, despike_info.sigma);
end
fprintf('>>> [Step3][%s] Loaded BTT sensors: [%s]\n', case_name, num2str(loaded_sensor_ids));
for i = 1:numel(btt_data)
    if isempty(btt_data{i})
        continue;
    end
    fprintf('>>> [Step3][%s] CH%d source: %s\n', case_name, btt_data{i}.sensor_id, btt_data{i}.mode);
end

result = struct();
result.case_name = case_name;
result.btt_output_dir = btt_output_dir;
result.strain_file = strain_file;
result.strain_case_dir = strain_case_dir;
result.loaded_sensor_ids = loaded_sensor_ids;
result.alignment_offset_sec = alignment_offset_sec;
result.btt_time_shift_sec = btt_time_shift_sec;
result.strain_reverse_shift_sec = btt_time_shift_sec;
result.strain_time_offset_total_sec = alignment_offset_sec - btt_time_shift_sec;
result.alignment_meta = alignment_meta;
result.Fs_strain = Fs_strain;
result.load_modes = load_modes(load_modes ~= "");
result.sensor_config_path = sensor_config_path;
result.rpm_time_s = rpm_time_s;
result.rpm_values = rpm_values;
result.despike_info = despike_info;

all_btt_times = [];
for i = 1:numel(btt_data)
    if isempty(btt_data{i})
        continue;
    end
    all_btt_times = [all_btt_times; btt_data{i}.time(:)]; %#ok<AGROW>
end
btt_time_span = [min(all_btt_times), max(all_btt_times)];
[sync_detail_xlim, fft_window, activity_info] = compute_activity_windows_local(btt_data, all_btt_times, t_strain_aligned, cfg);
btt_full_xlim = btt_time_span;
result.sync_xlim = btt_full_xlim;
result.sync_detail_xlim = sync_detail_xlim;
result.btt_full_xlim = btt_full_xlim;
result.fft_window = fft_window;
result.activity_info = activity_info;

if ~show_plots
    save(fullfile(btt_output_dir, 'Step3_Result_20251222.mat'), 'result');
    return;
end

plot_sync_overview_local(btt_data, cfg.step3_default_btt_sensor_ids, t_strain_aligned, v_strain, ...
    cfg.blades_num, sync_detail_xlim, alignment_offset_sec, btt_time_shift_sec, case_name, rpm_time_s, rpm_values, ...
    btt_time_span, alignment_meta, activity_info, despike_info);
plot_strain_stft_local(t_strain_raw, t_strain_aligned, v_strain, Fs_strain, sync_detail_xlim, case_name, despike_info);
plot_local_fft_local(t_strain_aligned, v_strain, Fs_strain, fft_window, case_name);

save(fullfile(btt_output_dir, 'Step3_Result_20251222.mat'), 'result');
fprintf('>>> [Step3][%s] Step3 visualization completed.\n', case_name);
end


function data = load_btt_display_data_local(output_dir, sid, Sensor_Config, cfg, case_name)
data = [];
raw_path = fullfile(output_dir, sprintf('jilublade_probe%d.mat', sid));
vib_path = fullfile(output_dir, sprintf('jilublade_probe%d_vib_final.mat', sid));

if ~isfile(raw_path)
    return;
end

loaded = load(raw_path);
fn = fieldnames(loaded);
jilublade = loaded.(fn{1});

if size(jilublade, 2) < 4
    error('File %s does not contain Step2-compatible jilublade columns.', raw_path);
end

fprintf('>>> [Step3][%s] CH%d: converting Step2 arrivals to displacement...\n', case_name, sid);
[jilublade_vib, disp_summary] = compute_displacement_from_step2_local(jilublade, output_dir, sid, Sensor_Config, cfg);
jilublade = jilublade_vib; %#ok<NASGU>
save(vib_path, 'jilublade');

data = struct();
data.sensor_id = sid;
data.mode = 'step3_displacement_from_step2';
data.file_path = vib_path;
data.matrix = jilublade_vib;
data.time_raw = jilublade_vib(:, 3);
data.time = jilublade_vib(:, 3);
data.blade_id = jilublade_vib(:, 4);
data.value = jilublade_vib(:, 6);
data.series_kind = 'displacement_mm';
data.has_value_series = true;
data.summary = disp_summary;
data.time_shift_sec = 0;
fprintf('>>> [Step3][%s] CH%d: displacement points=%d, std=%.4f mm, saved to %s\n', ...
    case_name, sid, disp_summary.valid_count, disp_summary.std_mm, vib_path);
end


function [jilublade_vib, summary] = compute_displacement_from_step2_local(jilublade, output_dir, sid, Sensor_Config, cfg)
opr_path = fullfile(output_dir, 'jiluOPR.mat');
if ~isfile(opr_path)
    error('Step3 displacement requires jiluOPR.mat. Missing file: %s', opr_path);
end

loaded_opr = load(opr_path);
if ~isfield(loaded_opr, 'jiluOPR')
    error('File %s does not contain jiluOPR.', opr_path);
end

jiluOPR = loaded_opr.jiluOPR;
if size(jiluOPR, 2) < 1
    error('jiluOPR format is invalid in %s.', opr_path);
end

opr_times = jiluOPR(:, 1);
opr_times = opr_times(isfinite(opr_times));
if numel(opr_times) <= cfg.blades_num
    error('Not enough OPR times to build speed interpolant for CH%d.', sid);
end

if ~isfield(Sensor_Config, 'Standard_Relative_Angles') || sid > size(Sensor_Config.Standard_Relative_Angles, 1)
    error('Sensor_Config does not contain Standard_Relative_Angles for CH%d.', sid);
end

spd_t = opr_times(1:(end - cfg.blades_num));
spd_v = 360 ./ (opr_times((cfg.blades_num + 1):end) - opr_times(1:(end - cfg.blades_num)));
F_omega_deg = griddedInterpolant(spd_t, spd_v, 'linear', 'nearest');

n = size(jilublade, 1);
jilublade_vib = nan(n, 6);
jilublade_vib(:, 1:min(4, size(jilublade, 2))) = jilublade(:, 1:min(4, size(jilublade, 2)));

valid_mask = false(n, 1);
disp_mm = nan(n, 1);
theta_act_rel = nan(n, 1);

for i = 1:n
    t_meas = jilublade(i, 3);
    blade_id = jilublade(i, 4);
    if ~isfinite(t_meas) || ~isfinite(blade_id) || blade_id < 1 || blade_id > cfg.blades_num
        continue;
    end

    opr_reference = build_opr_reference_from_jilu_local(jiluOPR, cfg.blades_num, cfg.r_tip_mm);
    theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference);
    if ~isfinite(theta_std)
        continue;
    end

    idx_prev_opr = find(opr_times < t_meas, 1, 'last');
    if isempty(idx_prev_opr)
        continue;
    end

    t_ref = opr_times(idx_prev_opr);
    t_grid = linspace(t_ref, t_meas, 10);
    theta_act = trapz(t_grid, F_omega_deg(t_grid));
    diff_deg = wrap_to_180_local(theta_act - theta_std);

    theta_act_rel(i) = theta_act;
    disp_mm(i) = diff_deg * (pi / 180) * cfg.r_tip_mm;
    valid_mask(i) = true;
end

jilublade_vib(:, 5) = theta_act_rel;
jilublade_vib(:, 6) = disp_mm;

summary = struct();
summary.valid_count = nnz(valid_mask);
summary.std_mm = std(disp_mm(valid_mask));
summary.mean_mm = mean(disp_mm(valid_mask));
summary.min_mm = min(disp_mm(valid_mask));
summary.max_mm = max(disp_mm(valid_mask));
end


function theta_std = read_opr_center_standard_angle_local(Sensor_Config, sid, blade_id, opr_reference)
if isfield(Sensor_Config, 'Standard_Relative_Angles_OPRCenter')
    theta_std = Sensor_Config.Standard_Relative_Angles_OPRCenter(sid, blade_id);
    return;
end
theta_std = Sensor_Config.Standard_Relative_Angles(sid, blade_id);
if isfield(Sensor_Config, 'Standard_Relative_Angles_Reference') && ...
        strcmpi(string(Sensor_Config.Standard_Relative_Angles_Reference), "opr_pulse_center")
    return;
end
if isstruct(opr_reference) && isfield(opr_reference, 'phase_shift_deg') && ...
        isfinite(opr_reference.phase_shift_deg)
    theta_std = theta_std - opr_reference.phase_shift_deg;
end
end


function opr_reference = build_opr_reference_from_jilu_local(jiluOPR, pulses_per_rev, r_tip_mm)
opr_reference = struct('phase_shift_deg', 0, 'phase_shift_mm', 0, ...
    'median_center_minus_start_s', 0);
if size(jiluOPR, 2) < 2 || pulses_per_rev < 1
    return;
end
center_time = jiluOPR(:, 1);
start_time = jiluOPR(:, 2);
n = min(numel(center_time) - pulses_per_rev, numel(start_time));
if n < 1
    return;
end
dt_center = center_time(1:n) - start_time(1:n);
dt_rev = center_time((1:n) + pulses_per_rev) - center_time(1:n);
valid = isfinite(dt_center) & isfinite(dt_rev) & dt_rev > eps;
if ~any(valid)
    return;
end
shift_deg = 360 * dt_center(valid) ./ dt_rev(valid);
opr_reference.phase_shift_deg = median(shift_deg, 'omitnan');
opr_reference.phase_shift_mm = opr_reference.phase_shift_deg * (pi / 180) * r_tip_mm;
opr_reference.median_center_minus_start_s = median(dt_center(valid), 'omitnan');
end


function [btt_time_shift_sec, meta] = load_step3_alignment_shift_local(output_dir, case_name)
btt_time_shift_sec = 0;
meta = struct();
meta.used_alignment_file = false;
meta.alignment_file = '';
meta.best_tau_sec = 0;
meta.best_order = nan;

alignment_path = fullfile(output_dir, 'Step3_Spectrum_RPM_Alignment_20251222.mat');
if ~isfile(alignment_path)
    return;
end

loaded = load(alignment_path);
if ~isfield(loaded, 'result')
    return;
end

alignment_result = loaded.result;
if ~isfield(alignment_result, 'best_tau_sec') || ~isfinite(alignment_result.best_tau_sec)
    return;
end

btt_time_shift_sec = alignment_result.best_tau_sec;
meta.used_alignment_file = true;
meta.alignment_file = alignment_path;
meta.best_tau_sec = btt_time_shift_sec;
if isfield(alignment_result, 'best_order') && isfinite(alignment_result.best_order)
    meta.best_order = alignment_result.best_order;
end

fprintf('>>> [Step3][%s] Alignment file: %s\n', case_name, alignment_path);
if isfinite(meta.best_order)
    fprintf('>>> [Step3][%s] Alignment file best visible order = EO%d\n', case_name, meta.best_order);
end
end


function ang = wrap_to_180_local(ang)
ang = mod(ang + 180, 360) - 180;
end


function strain_case_dir = map_strain_case_dir_local(case_name, strain_root)
switch char(case_name)
    case {'1000_2500_3500', '1000rpm无振动'}
        strain_case_dir = strain_root;
    case '20250526_3150'
        strain_case_dir = fullfile(strain_root, '20251222_3150');
    case '20250526_910'
        strain_case_dir = fullfile(strain_root, '20251222_900');
    case '20250526_2500-3500_t200'
        strain_case_dir = fullfile(strain_root, '20251222_2500_3500_time200');
    case '20250526_2500-3500_t400'
        strain_case_dir = fullfile(strain_root, '20251222_2500_3500_time400');
    case '20250526_2500-3500_t800'
        strain_case_dir = fullfile(strain_root, '20251222_2500_3500_time800');
    otherwise
        candidate = fullfile(strain_root, char(case_name));
        if isfolder(candidate)
            strain_case_dir = candidate;
        elseif has_strain_files_in_root_local(strain_root)
            strain_case_dir = strain_root;
        else
            error('No mapped strain directory for case %s', case_name);
        end
end
end


function tf = has_strain_files_in_root_local(strain_root)
tf = false;
if ~isfolder(strain_root)
    return;
end
d = dir(fullfile(strain_root, 'AI1-*.mat'));
tf = ~isempty(d);
end


function strain_file = pick_latest_strain_file_local(strain_dir, channel_tag)
d = dir(fullfile(strain_dir, sprintf('%s*.mat', channel_tag)));
if isempty(d)
    error('No strain files found under %s with tag %s', strain_dir, channel_tag);
end
[~, idx] = max([d.datenum]);
strain_file = fullfile(strain_dir, d(idx).name);
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


function [v_out, info] = preprocess_strain_signal_local(v_in, Fs_strain, cfg)
v_out = v_in(:);
info = struct('enabled', false, 'num_replaced', 0, 'window_samples', 0, 'sigma', NaN);

if ~isfield(cfg, 'step3_despike_enable') || ~cfg.step3_despike_enable || numel(v_out) < 9
    return;
end

win = max(5, round(cfg.step3_despike_window_sec * Fs_strain));
if mod(win, 2) == 0
    win = win + 1;
end
half_win = floor(win / 2);

med = movmedian(v_out, win, 'omitnan');
resid = v_out - med;
mad_local = movmedian(abs(resid), win, 'omitnan') / 0.6745;
mad_local(mad_local < eps) = eps;

spike_mask = abs(resid) > cfg.step3_despike_sigma * mad_local;
spike_mask(1:half_win) = false;
spike_mask((end - half_win + 1):end) = false;

if any(spike_mask)
    idx = (1:numel(v_out)).';
    keep = ~spike_mask;
    if nnz(keep) >= 2
        v_out(spike_mask) = interp1(idx(keep), v_out(keep), idx(spike_mask), 'linear', 'extrap');
    end
end

info.enabled = true;
info.num_replaced = nnz(spike_mask);
info.window_samples = win;
info.sigma = cfg.step3_despike_sigma;
end


function [rpm_time_s, rpm_values] = load_rpm_trace_local(output_dir, cfg)
rpm_time_s = [];
rpm_values = [];

omega_path = fullfile(output_dir, 'omega.mat');
if isfile(omega_path)
    loaded = load(omega_path);
    if isfield(loaded, 'omega_time_s') && isfield(loaded, 'omega_rpm')
        rpm_time_s = loaded.omega_time_s(:);
        rpm_values = loaded.omega_rpm(:);
        return;
    end
end

opr_path = fullfile(output_dir, 'jiluOPR.mat');
if ~isfile(opr_path)
    return;
end

loaded = load(opr_path);
if ~isfield(loaded, 'jiluOPR')
    return;
end

opr_times = loaded.jiluOPR(:, 1);
opr_times = opr_times(isfinite(opr_times));
if numel(opr_times) <= cfg.blades_num
    return;
end

rev_period = opr_times((cfg.blades_num + 1):end) - opr_times(1:(end - cfg.blades_num));
rpm_time_s = opr_times(1:(end - cfg.blades_num));
rpm_values = 60 ./ rev_period;
end


function sync_xlim = compute_sync_xlim_local(all_btt_times, t_strain_aligned)
if isempty(all_btt_times)
    sync_xlim = [t_strain_aligned(1), t_strain_aligned(end)];
    return;
end

t_min = max(t_strain_aligned(1), min(all_btt_times) - 1);
t_max = min(t_strain_aligned(end), max(all_btt_times) + 1);
if t_max <= t_min
    sync_xlim = [t_strain_aligned(1), t_strain_aligned(end)];
else
    sync_xlim = [t_min, t_max];
end
end


function [sync_xlim, fft_window, activity_info] = compute_activity_windows_local(btt_data, all_btt_times, t_strain_aligned, cfg)
fallback_xlim = compute_sync_xlim_local(all_btt_times, t_strain_aligned);
btt_full_window = [min(all_btt_times), max(all_btt_times)];
activity_info = struct();
activity_info.method = 'fallback';
activity_info.peak_time = mean(fallback_xlim);
activity_info.active_window = fallback_xlim;
activity_info.btt_full_window = btt_full_window;
activity_info.grid_time = [];
activity_info.activity_envelope = [];

if isempty(all_btt_times)
    sync_xlim = fallback_xlim;
    fft_window = compute_stable_fft_window_local(btt_data, sync_xlim, btt_full_window, t_strain_aligned, cfg);
    activity_info.fft_window = fft_window;
    return;
end

t_min = min(all_btt_times);
t_max = max(all_btt_times);
dt = cfg.step3_activity_grid_dt_sec;
if ~isfinite(dt) || dt <= 0
    dt = 0.02;
end

t_grid = (t_min:dt:t_max).';
if numel(t_grid) < 5
    sync_xlim = fallback_xlim;
    fft_window = compute_stable_fft_window_local(btt_data, sync_xlim, btt_full_window, t_strain_aligned, cfg);
    activity_info.fft_window = fft_window;
    return;
end

value_mat = nan(numel(t_grid), numel(btt_data));
for i = 1:numel(btt_data)
    if isempty(btt_data{i}) || ~isfield(btt_data{i}, 'value') || isempty(btt_data{i}.value)
        continue;
    end

    ti = btt_data{i}.time(:);
    vi = abs(btt_data{i}.value(:));
    valid = isfinite(ti) & isfinite(vi);
    if nnz(valid) < 3
        continue;
    end

    [ti_unique, ia] = unique(ti(valid));
    vi_unique = vi(valid);
    vi_unique = vi_unique(ia);
    if numel(ti_unique) < 3
        continue;
    end

    value_mat(:, i) = interp1(ti_unique, vi_unique, t_grid, 'linear', nan);
end

activity = mean(value_mat, 2, 'omitnan');
if all(~isfinite(activity))
    sync_xlim = fallback_xlim;
    fft_window = compute_stable_fft_window_local(btt_data, sync_xlim, btt_full_window, t_strain_aligned, cfg);
    activity_info.fft_window = fft_window;
    return;
end

activity(~isfinite(activity)) = 0;
smooth_n = max(3, round(cfg.step3_activity_smooth_sec / dt));
if mod(smooth_n, 2) == 0
    smooth_n = smooth_n + 1;
end
activity_smooth = smoothdata(activity, 'movmean', smooth_n);

[peak_amp, peak_idx] = max(activity_smooth);
peak_time = t_grid(peak_idx);
threshold = cfg.step3_activity_threshold_ratio * peak_amp;

left_idx = peak_idx;
while left_idx > 1 && activity_smooth(left_idx) >= threshold
    left_idx = left_idx - 1;
end
right_idx = peak_idx;
while right_idx < numel(activity_smooth) && activity_smooth(right_idx) >= threshold
    right_idx = right_idx + 1;
end

active_window = [t_grid(max(1, left_idx)), t_grid(min(numel(t_grid), right_idx))];
active_window = active_window + [-cfg.step3_activity_pad_sec, cfg.step3_activity_pad_sec];

min_width = cfg.step3_activity_min_window_sec;
if diff(active_window) < min_width
    half_extra = 0.5 * (min_width - diff(active_window));
    active_window = active_window + [-half_extra, half_extra];
end

active_window(1) = max(fallback_xlim(1), active_window(1));
active_window(2) = min(fallback_xlim(2), active_window(2));
if active_window(2) <= active_window(1)
    active_window = fallback_xlim;
end

sync_xlim = active_window;
fft_window = compute_stable_fft_window_local(btt_data, sync_xlim, btt_full_window, t_strain_aligned, cfg);

activity_info.method = 'btt_activity_envelope';
activity_info.peak_time = peak_time;
activity_info.peak_amplitude = peak_amp;
activity_info.threshold = threshold;
activity_info.active_window = active_window;
activity_info.btt_full_window = fallback_xlim;
activity_info.grid_time = t_grid;
activity_info.activity_envelope = activity_smooth;
activity_info.fft_window = fft_window;
end


function fft_window = compute_stable_fft_window_local(btt_data, active_window, btt_full_window, t_strain_aligned, cfg)
if isfield(cfg, 'step3_manual_fft_window_sec') && numel(cfg.step3_manual_fft_window_sec) == 2
    fft_window = cfg.step3_manual_fft_window_sec(:).';
    fft_window(1) = max([fft_window(1), btt_full_window(1), t_strain_aligned(1)]);
    fft_window(2) = min([fft_window(2), btt_full_window(2), t_strain_aligned(end)]);
    if fft_window(2) > fft_window(1)
        return;
    end
end

window_sec = cfg.step3_activity_fft_window_sec;
half_width = 0.5 * window_sec;
dt = cfg.step3_activity_grid_dt_sec;
if ~isfinite(dt) || dt <= 0
    dt = 0.02;
end

t_grid = (active_window(1):dt:active_window(2)).';
if numel(t_grid) < 5
    fft_window = [active_window(1), min(active_window(1) + window_sec, active_window(2))];
    return;
end

value_mat = nan(numel(t_grid), numel(btt_data));
for i = 1:numel(btt_data)
    if isempty(btt_data{i}) || ~isfield(btt_data{i}, 'value') || isempty(btt_data{i}.value)
        continue;
    end

    ti = btt_data{i}.time(:);
    vi = abs(btt_data{i}.value(:));
    valid = isfinite(ti) & isfinite(vi) & ti >= active_window(1) & ti <= active_window(2);
    if nnz(valid) < 3
        continue;
    end

    [ti_unique, ia] = unique(ti(valid));
    vi_unique = vi(valid);
    vi_unique = vi_unique(ia);
    if numel(ti_unique) < 3
        continue;
    end

    value_mat(:, i) = interp1(ti_unique, vi_unique, t_grid, 'linear', nan);
end

activity = mean(value_mat, 2, 'omitnan');
activity(~isfinite(activity)) = 0;
activity_smooth = smoothdata(activity, 'movmean', max(3, round(cfg.step3_activity_smooth_sec / dt)));
activity_grad = abs([0; diff(activity_smooth)]) / max(dt, eps);

amp_norm = normalize_01_local(activity_smooth);
grad_norm = normalize_01_local(activity_grad);
score = amp_norm - cfg.step3_stability_derivative_weight * grad_norm;

win_pts = max(3, round(window_sec / dt));
score_mean = movmean(score, win_pts, 'omitnan');
activity_mean = movmean(activity_smooth, win_pts, 'omitnan');
objective = score_mean + 0.3 * normalize_01_local(activity_mean);

[~, best_idx] = max(objective);
center_time = t_grid(best_idx);
fft_window = [center_time - half_width, center_time + half_width];
fft_window(1) = max([fft_window(1), active_window(1), t_strain_aligned(1)]);
fft_window(2) = min([fft_window(2), active_window(2), t_strain_aligned(end)]);
if fft_window(2) <= fft_window(1)
    fft_window = [active_window(1), min(active_window(1) + window_sec, active_window(2))];
    fft_window(1) = max(fft_window(1), t_strain_aligned(1));
    fft_window(2) = min(fft_window(2), t_strain_aligned(end));
end
end


function x = normalize_01_local(x)
x = x(:);
if isempty(x) || all(~isfinite(x))
    x = zeros(size(x));
    return;
end
x(~isfinite(x)) = nan;
xmin = min(x, [], 'omitnan');
xmax = max(x, [], 'omitnan');
if ~isfinite(xmin) || ~isfinite(xmax) || xmax <= xmin
    x = zeros(size(x));
else
    x = (x - xmin) ./ (xmax - xmin);
    x(~isfinite(x)) = 0;
end
end


function plot_sync_overview_local(btt_data, sensor_ids, t_strain_aligned, v_strain, blades_num, sync_xlim, alignment_offset_sec, btt_time_shift_sec, case_name, rpm_time_s, rpm_values, btt_time_span, alignment_meta, activity_info, despike_info)
fig = figure('Name', sprintf('20251222 Step3 Sync Overview - %s', case_name), ...
    'Color', 'w', 'Position', [80, 40, 1500, 1040], 'NumberTitle', 'off');
tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

colors = lines(numel(sensor_ids));
full_xlim = [t_strain_aligned(1), t_strain_aligned(end)];

nexttile([1, 2]);
hold on;
grid on;
box on;
yyaxis left;
stride = max(1, floor(numel(t_strain_aligned) / 40000));
plot(t_strain_aligned(1:stride:end), v_strain(1:stride:end), '-', ...
    'Color', [0.20 0.60 0.20 0.35], 'LineWidth', 0.5, 'DisplayName', 'Despiked strain (full record)');
ylabel('Strain (V)');
yl = ylim;
plot([btt_time_span(1), btt_time_span(1)], yl, '--', 'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.2, 'DisplayName', 'BTT start');
plot([btt_time_span(2), btt_time_span(2)], yl, '-.', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'DisplayName', 'BTT end');
text(mean(btt_time_span), yl(2) - 0.08 * (yl(2) - yl(1)), ...
    sprintf('BTT recorded window: %.2f s to %.2f s | Strain reverse shift = %.3f s', ...
    btt_time_span(1), btt_time_span(2), btt_time_shift_sec), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'BackgroundColor', 'w');
if ~isempty(activity_info) && isfield(activity_info, 'active_window') && numel(activity_info.active_window) == 2
    aw = activity_info.active_window;
    plot([aw(1), aw(1)], yl, ':', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.3, 'DisplayName', 'Active-window start');
    plot([aw(2), aw(2)], yl, ':', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.3, 'HandleVisibility', 'off');
end
if ~isempty(activity_info) && isfield(activity_info, 'fft_window') && numel(activity_info.fft_window) == 2
    fw = activity_info.fft_window;
    plot([fw(1), fw(1)], yl, '-', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.2, 'DisplayName', 'FFT-window start');
    plot([fw(2), fw(2)], yl, '-', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

yyaxis right;
if ~isempty(rpm_time_s)
    plot(rpm_time_s, rpm_values, 'k-', 'LineWidth', 1.1, 'DisplayName', 'BTT-derived RPM (raw time)');
end
ylabel('RPM');
xlabel('Aligned time (s)');
if alignment_meta.used_alignment_file && isfinite(alignment_meta.best_order)
    title(sprintf('Full-record reverse-shifted strain + raw BTT RPM | %s | EO%d | strain shift = -%.3f s', ...
        case_name, alignment_meta.best_order, btt_time_shift_sec), 'Interpreter', 'none');
else
    title(sprintf('Full-record reverse-shifted strain + raw BTT RPM | %s | strain shift = -%.3f s', ...
        case_name, btt_time_shift_sec), 'Interpreter', 'none');
end
if ~isempty(despike_info) && isfield(despike_info, 'enabled') && despike_info.enabled
    text(mean(btt_time_span), yl(1) + 0.10 * (yl(2) - yl(1)), ...
        sprintf('Strain despike copy only | replaced %d points', despike_info.num_replaced), ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
end
xlim(full_xlim);
legend('Location', 'best');

for blade_id = 1:blades_num
    nexttile;
    hold on;
    grid on;
    box on;

    yyaxis left;
    plot(t_strain_aligned, v_strain, '-', 'Color', [0.20 0.60 0.20 0.30], ...
        'LineWidth', 0.6, 'DisplayName', 'Despiked strain');
    ylabel('Strain (V)');

    yyaxis right;
    has_btt = false;
    for i = 1:numel(sensor_ids)
        if isempty(btt_data{i})
            continue;
        end

        data = btt_data{i};
        mask = data.blade_id == blade_id;
        if ~any(mask)
            continue;
        end

        has_btt = true;
        if data.has_value_series
            plot(data.time(mask), data.value(mask), '-', 'LineWidth', 0.9, ...
                'Color', colors(i, :), 'DisplayName', sprintf('CH%d displacement', data.sensor_id));
        else
            scatter(data.time(mask), data.value(mask), 14, colors(i, :), 'filled', ...
                'DisplayName', sprintf('CH%d arrivals', data.sensor_id));
        end
    end

    if has_btt
        yline(0, 'k:', 'HandleVisibility', 'off');
    end
    if ~isempty(activity_info) && isfield(activity_info, 'fft_window') && numel(activity_info.fft_window) == 2
        xline(activity_info.fft_window(1), '-', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        xline(activity_info.fft_window(2), '-', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    ylabel('BTT displacement (mm)');
    xlabel('Time (s)');
    title(sprintf('Blade %d sync | strain offset = %.4f s | strain reverse shift = %.4f s', ...
        blade_id, alignment_offset_sec, btt_time_shift_sec));
    xlim(sync_xlim);
    legend('Location', 'best');
end
end


function plot_strain_stft_local(t_strain_raw, t_strain_aligned, v_strain, Fs_strain, sync_xlim, case_name, despike_info)
fprintf('>>> [Step3][%s] Building strain STFT...\n', case_name);

win_time = 0.2;
overlap_ratio = 0.5;
win_len = max(128, floor(win_time * Fs_strain));
hop_len = max(1, floor(win_len * (1 - overlap_ratio)));
num_segments = floor((numel(v_strain) - win_len) / hop_len) + 1;

w_func = hann(win_len);
acf = 1 / mean(w_func);
f_axis = (0:floor(win_len / 2)) * (Fs_strain / win_len);
amplitudes = zeros(num_segments, numel(f_axis));
t_segments_aligned = zeros(num_segments, 1);

for i = 1:num_segments
    idx_start = (i - 1) * hop_len + 1;
    idx_end = idx_start + win_len - 1;
    seg = v_strain(idx_start:idx_end);
    raw_center_time = t_strain_raw(idx_start + floor(win_len / 2));
    t_segments_aligned(i) = raw_center_time + (t_strain_aligned(1) - t_strain_raw(1));

    seg_detrend = seg - mean(seg);
    seg_windowed = seg_detrend .* w_func;
    Y = fft(seg_windowed);
    P2 = abs(Y / win_len);
    P1 = P2(1:floor(win_len / 2) + 1);
    P1(2:end-1) = 2 * P1(2:end-1);
    amplitudes(i, :) = P1 * acf;
end

figure('Name', sprintf('20251222 Step3 Strain STFT - %s', case_name), ...
    'Color', 'w', 'Position', [160, 120, 1200, 520], 'NumberTitle', 'off');
surf(t_segments_aligned, f_axis, amplitudes', 'EdgeColor', 'none');
shading interp;
axis tight;
% view(0, 90);
colormap(jet);
cb = colorbar;
ylabel(cb, 'Amplitude');
xlabel('Aligned time (s)');
ylabel('Frequency (Hz)');
if ~isempty(despike_info) && isfield(despike_info, 'enabled') && despike_info.enabled
    title(sprintf('Despiked strain STFT (aligned axis): %s', case_name), 'Interpreter', 'none');
else
    title(sprintf('Strain STFT (aligned axis): %s', case_name), 'Interpreter', 'none');
end
ylim([0, 800]);
xlim(sync_xlim);
end


function plot_local_fft_local(t_strain_aligned, v_strain, Fs_strain, fft_window, case_name)
fprintf('>>> [Step3][%s] Building local FFT in [%.3f, %.3f] s...\n', case_name, fft_window(1), fft_window(2));

mask = t_strain_aligned >= fft_window(1) & t_strain_aligned <= fft_window(2);
if nnz(mask) < 32
    warning('Step3 FFT window contains too few points. Skipping local FFT.');
    return;
end

sig = v_strain(mask);
t_seg = t_strain_aligned(mask);
sig = sig - mean(sig);

L = numel(sig);
win = hann(L);
coherent_gain = mean(win);
sig_windowed = sig .* win;
Y = fft(sig_windowed);
P2 = abs(Y / (L * coherent_gain));
P1 = P2(1:floor(L / 2) + 1);
P1(2:end-1) = 2 * P1(2:end-1);
f = Fs_strain * (0:floor(L / 2)) / L;

figure('Name', sprintf('20251222 Step3 Local FFT - %s', case_name), ...
    'Color', 'w', 'Position', [220, 160, 1000, 680], 'NumberTitle', 'off');

subplot(2, 1, 1);
plot(t_seg, v_strain(mask), 'k-', 'LineWidth', 0.8);
grid on;
xlabel('Aligned time (s)');
ylabel('Strain (V)');
title(sprintf('Local time-domain strain | [%.3f, %.3f] s | Hann window for FFT', fft_window(1), fft_window(2)));

subplot(2, 1, 2);
plot(f, P1, 'LineWidth', 1.4, 'Color', [0.85 0.33 0.10]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Local FFT with Hann window');
xlim([0, 800]);
[pk_amp, pk_idx] = max(P1);
hold on;
plot(f(pk_idx), pk_amp, 'ko', 'MarkerFaceColor', 'y');
text(f(pk_idx) + 10, pk_amp, sprintf('Peak %.2f Hz', f(pk_idx)), ...
    'FontSize', 10, 'FontWeight', 'bold');
end

