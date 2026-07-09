function speed_diag = compute_opr_speed_diagnostics_20251222(opr_times, opr_pulses_per_rev, time_window)
%COMPUTE_OPR_SPEED_DIAGNOSTICS_20250527 Build OPR-based speed diagnostics.

if nargin < 2 || isempty(opr_pulses_per_rev)
    opr_pulses_per_rev = 6;
end
if nargin < 3 || isempty(time_window)
    time_window = [min(opr_times), max(opr_times)];
end

opr_times = opr_times(:);
speed_diag = struct();
speed_diag.opr_times = opr_times;
speed_diag.opr_pulses_per_rev = opr_pulses_per_rev;
speed_diag.time_window = time_window(:).';

pulse_mid_time = [];
pulse_rot_freq_hz = [];
pulse_rpm = [];
selected_pulse_mask = [];
if numel(opr_times) >= 2
    pulse_mid_time = 0.5 * (opr_times(1:(end - 1)) + opr_times(2:end));
    pulse_dt = diff(opr_times);
    pulse_rot_freq_hz = 1 ./ max(pulse_dt, eps) ./ opr_pulses_per_rev;
    pulse_rpm = 60 .* pulse_rot_freq_hz;
    selected_pulse_mask = pulse_mid_time >= time_window(1) & pulse_mid_time <= time_window(2);
end

rev_mid_time = [];
rev_rot_freq_hz = [];
rev_rpm = [];
selected_rev_mask = [];
if numel(opr_times) > opr_pulses_per_rev
    rev_t0 = opr_times(1:(end - opr_pulses_per_rev));
    rev_t1 = opr_times((opr_pulses_per_rev + 1):end);
    rev_mid_time = 0.5 * (rev_t0 + rev_t1);
    rev_period = rev_t1 - rev_t0;
    rev_rot_freq_hz = 1 ./ max(rev_period, eps);
    rev_rpm = 60 .* rev_rot_freq_hz;
    selected_rev_mask = rev_mid_time >= time_window(1) & rev_mid_time <= time_window(2);
end

speed_diag.pulse_mid_time = pulse_mid_time;
speed_diag.pulse_rot_freq_hz = pulse_rot_freq_hz;
speed_diag.pulse_rpm = pulse_rpm;
speed_diag.selected_pulse_mask = selected_pulse_mask;
speed_diag.rev_mid_time = rev_mid_time;
speed_diag.rev_rot_freq_hz = rev_rot_freq_hz;
speed_diag.rev_rpm = rev_rpm;
speed_diag.selected_rev_mask = selected_rev_mask;

speed_diag.global_rot_freq_mean_hz = median(rev_rot_freq_hz, 'omitnan');
speed_diag.global_rot_rpm_mean = 60 * speed_diag.global_rot_freq_mean_hz;

selected_rot_freq = rev_rot_freq_hz(selected_rev_mask);
if isempty(selected_rot_freq)
    selected_rot_freq = pulse_rot_freq_hz(selected_pulse_mask);
end
if isempty(selected_rot_freq)
    selected_rot_freq = rev_rot_freq_hz;
end
if isempty(selected_rot_freq)
    selected_rot_freq = pulse_rot_freq_hz;
end

speed_diag.local_rot_freq_mean_hz = median(selected_rot_freq, 'omitnan');
speed_diag.local_rot_rpm_mean = 60 * speed_diag.local_rot_freq_mean_hz;
end
