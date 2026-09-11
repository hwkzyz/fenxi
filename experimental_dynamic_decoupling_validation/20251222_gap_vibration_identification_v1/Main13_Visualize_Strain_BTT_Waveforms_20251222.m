%% Visualize strain-BTT validation (20251222)
% Dedicated visualization program for the identification and strain results.
% Focus:
%   1. vibration parameter trends
%   2. time-domain waveform comparison by sensor

clc
clear
close all

script_dir = fileparts(mfilename('fullpath'));
packageCfg = Setup_Paths_20251222();
caseCfg = CaseConfig();

cfg = struct();
cfg.target_blade = caseCfg.bladeId;
cfg.analysis_sensors = caseCfg.sensorIds;
cfg.dynamic_case_name = '1000_2500_3500';
cfg.analysis_start_time = caseCfg.flowConfig.identification.analysisStartTimeSec;
cfg.method_name = 'gap_only';
cfg.primary_strain_channel = 'AI1-03';
cfg.method_result_file = '';
cfg.method_result_dir = packageCfg.paths.gapResults;
cfg.show_figures = true;
cfg.save_figures = true;
cfg.figure_dir = fullfile(packageCfg.paths.figures, 'strain_btt_waveforms');
cfg.figure_style = 'nature';
cfg.export_pdf = true;
cfg.zoom_half_width_ms = 0.35;
cfg.enable_strain_reference = true;
% Strain-displacement transfer ratio used by the existing 20251222
% strain-spectrum validation route.
% R = abs(strain / u_tip_m), strain is dimensionless and u_tip_m is in m.
% The measured strain channel used below is in microstrain, so:
% u_tip_mm = strain_microstrain * 1e-6 / R * 1000 = strain_microstrain / R / 1000.
cfg.fe_transfer_ratio_label = 'strain-displacement calibration, R=1.266 1/m';
cfg.fe_transfer_ratio_1_per_m = 1.266;
cfg.strain_to_tip_mm = 1 / cfg.fe_transfer_ratio_1_per_m / 1000;
cfg.strain_alignment_file = fullfile(packageCfg.paths.prepared, ...
    'Step04_BTT_STE_Resonance_Regions_20251222.mat');
% This dataset has no separate low-speed strain recording.
cfg.low_speed_strain_dir = '';
cfg.low_speed_strain_file = '';
cfg.strain_reference_window_revs = 3;
cfg.strain_reference_freq_half_band_hz = 8;
cfg.strain_reference_freq_step_hz = 0.1;
cfg.strain_diagnostic_spectrum_half_band_hz = 60;
cfg.strain_diagnostic_time_half_width_ms = 4;
cfg.strain_raw_freq_band_hz = [520, 640];
cfg.strain_raw_spectrum_xlim_hz = [0, 1000];
cfg.strain_operating_freq_xlim_hz = [0, 1000];
cfg.strain_raw_time_context_sec = 0.08;
cfg.strain_raw_stft_window_sec = 0.18;
cfg.strain_raw_stft_overlap_ratio = 0.90;
cfg.strain_component_window_half_width_ms = 5.0;
cfg.strain_component_time_samples = 900;
cfg.export_strain_reference_csv = true;

if ~exist(cfg.figure_dir, 'dir')
    mkdir(cfg.figure_dir);
end

result_file = resolve_method_result_file_local(cfg);
loaded = load(result_file, 'Result');
MethodResult = loaded.Result;
Result = adapt_method_result_local(MethodResult, cfg);
style = build_step6_style_local(cfg);
StrainRef01 = build_step6_strain_reference_local(Result, cfg, 'AI1-01');
StrainRef03 = build_step6_strain_reference_local(Result, cfg);
StrainRef = StrainRef03;

print_method_summary_local(MethodResult, Result, cfg, result_file);
print_step6_strain_reference_summary_local(StrainRef);
print_step6_strain_reference_summary_local(StrainRef01);
if is_valid_strain_reference_local(StrainRef01) && is_valid_strain_reference_local(StrainRef03)
    export_step6_dual_strain_reference_csv_local(StrainRef01, StrainRef03, cfg);
    corr01 = corr(StrainRef01.BTTAmpMM, StrainRef01.AmpTipMM, 'Rows', 'complete');
    corr03 = corr(StrainRef03.BTTAmpMM, StrainRef03.AmpTipMM, 'Rows', 'complete');
    fprintf('  Amplitude-trend Pearson correlation | AI1-01 = %.4f | AI1-03 = %.4f\n', ...
        corr01, corr03);
end

if cfg.show_figures
    fig_param = plot_parameter_trends_local(Result, cfg, style, StrainRef);
    save_figure_step6_local(fig_param, cfg.figure_dir, cfg.target_blade, 'ParameterTrends', cfg);

    if is_valid_strain_reference_local(StrainRef)
        fig_strain_operating = plot_strain_operating_condition_evidence_local(Result, cfg, style, StrainRef);
        save_figure_step6_local(fig_strain_operating, cfg.figure_dir, cfg.target_blade, 'StrainOperatingConditionEvidence', cfg);

        fig_strain_raw = plot_strain_raw_evidence_local(Result, cfg, style, StrainRef);
        save_figure_step6_local(fig_strain_raw, cfg.figure_dir, cfg.target_blade, 'StrainRawEvidence', cfg);

        fig_strain = plot_strain_reference_validation_local(Result, cfg, style, StrainRef);
        save_figure_step6_local(fig_strain, cfg.figure_dir, cfg.target_blade, 'StrainComponentValidation', cfg);

        if is_valid_strain_reference_local(StrainRef01)
            fig_strain_channels = plot_strain_channel_comparison_local( ...
                Result, cfg, style, StrainRef01, StrainRef03);
            save_figure_step6_local(fig_strain_channels, cfg.figure_dir, ...
                cfg.target_blade, 'StrainChannelComparison', cfg);
        end
    end

end


function result_file = resolve_method_result_file_local(cfg)
result_file = cfg.method_result_file;
if ~isfile(result_file)
    sensorTag = ['S', sprintf('%d', cfg.analysis_sensors)];
    files = dir(fullfile(cfg.method_result_dir, sprintf( ...
        'Main_GapAware_VPTopK_FullWave_20251222_B%d_%s*.mat', ...
        cfg.target_blade, sensorTag)));
    if isempty(files)
        error('Run Main10 first; no local gap-aware result was found under %s.', ...
            cfg.method_result_dir);
    end
    [~, idx] = max([files.datenum]);
    result_file = fullfile(files(idx).folder, files(idx).name);
end
end


function Result = adapt_method_result_local(MethodResult, cfg)
window_result = MethodResult.WindowResult(:);
n_window = numel(window_result);
window_id = nan(n_window, 1);
window_center = nan(n_window, 1);
rot_rpm = nan(n_window, 1);
amp = nan(n_window, 1);
freq = nan(n_window, 1);
eo = nan(n_window, 1);
phase = nan(n_window, 1);
d0 = nan(n_window, 1);
weighted_rmse = nan(n_window, 1);
plain_rmse = nan(n_window, 1);
for i_window = 1:n_window
    fit = window_result(i_window).modelFits.(cfg.method_name);
    window_id(i_window) = window_result(i_window).windowId;
    window_center(i_window) = mean(window_result(i_window).timeWindow);
    rot_rpm(i_window) = window_result(i_window).rotRpmMean;
    amp(i_window) = fit.amplitudeMm;
    freq(i_window) = fit.freqHz;
    eo(i_window) = fit.EO;
    phase(i_window) = fit.phaseRad;
    d0(i_window) = fit.dxMm;
    weighted_rmse(i_window) = fit.weightedRmseMv ./ 1000;
    plain_rmse(i_window) = fit.plainRmseMv ./ 1000;
end
Result = struct();
Result.BladeID = cfg.target_blade;
Result.Config = struct('dynamic_case_name', cfg.dynamic_case_name);
Result.Trends = struct('WindowID', window_id, ...
    'WindowCenterTime', window_center, 'RotRPM', rot_rpm, ...
    'Amp', amp, 'Freq', freq, 'EO', eo, 'Phase', phase, 'd0', d0, ...
    'WeightedRMSE', weighted_rmse, 'PlainRMSE', plain_rmse);
end


function print_method_summary_local(MethodResult, Result, cfg, result_file)
best_index = MethodResult.BestWindowIndex;
fprintf(['Strain-BTT validation using the local decoupling method\n' ...
    '  Result file = %s\n  Model = %s | windows = %d | best window = %d\n' ...
    '  Dominant EO = %d | mean frequency = %.4f Hz | mean amplitude = %.6f mm\n'], ...
    result_file, cfg.method_name, numel(Result.Trends.WindowID), best_index, ...
    mode(Result.Trends.EO), mean(Result.Trends.Freq, 'omitnan'), ...
    mean(Result.Trends.Amp, 'omitnan'));
end


function StrainRef = build_step6_strain_reference_local(Result, cfg, requested_channel_tag)
StrainRef = struct();
StrainRef.available = false;
StrainRef.message = '';
if nargin < 3
    requested_channel_tag = '';
end

if ~isfield(cfg, 'enable_strain_reference') || ~cfg.enable_strain_reference
    StrainRef.message = 'Strain reference disabled.';
    return;
end

alignment_file = resolve_step6_strain_alignment_file_local(Result, cfg);
if ~isfile(alignment_file)
    StrainRef.message = sprintf('Strain alignment file not found: %s', alignment_file);
    warning('%s', StrainRef.message);
    return;
end

loaded_alignment = load(alignment_file);
if isfield(loaded_alignment, 'result')
    alignment = loaded_alignment.result;
elseif isfield(loaded_alignment, 'strainSource')
    source = loaded_alignment.strainSource;
    alignment = struct();
    alignment.strain_file = source.strain_file;
    alignment.best_tau_sec = -source.strain_time_offset_total_sec;
else
    StrainRef.message = sprintf('Alignment file has no result or strainSource: %s', alignment_file);
    warning('%s', StrainRef.message);
    return;
end

if ~isempty(requested_channel_tag)
    if ~isfield(alignment, 'strain_file') || isempty(alignment.strain_file)
        StrainRef.message = 'Alignment result does not provide a strain file.';
        warning('%s', StrainRef.message);
        return;
    end
    alignment.strain_file = resolve_strain_channel_peer_file_local( ...
        alignment.strain_file, requested_channel_tag);
end

if ~isfield(alignment, 'strain_file') || ~isfile(alignment.strain_file)
    StrainRef.message = 'Aligned strain file is missing or unavailable.';
    warning('%s', StrainRef.message);
    return;
end
if ~isfield(alignment, 'best_tau_sec') || ~isfinite(alignment.best_tau_sec)
    StrainRef.message = 'Alignment result does not provide a finite best_tau_sec.';
    warning('%s', StrainRef.message);
    return;
end

loaded_strain = load(alignment.strain_file);
if ~isfield(loaded_strain, 'Datas') || size(loaded_strain.Datas, 2) < 2
    StrainRef.message = sprintf('Strain file does not contain Datas(:,1:2): %s', alignment.strain_file);
    warning('%s', StrainRef.message);
    return;
end

t_strain = loaded_strain.Datas(:, 1);
v_strain = loaded_strain.Datas(:, 2);
if isfield(loaded_strain, 'SampleFrequency')
    Fs_strain = str2double(loaded_strain.SampleFrequency);
else
    Fs_strain = 1 / median(diff(t_strain));
end

K_Strain2MM = cfg.strain_to_tip_mm;
xv = Result.Trends.WindowID(:);
btt_center_time = Result.Trends.WindowCenterTime(:);
strain_center_time = btt_center_time + alignment.best_tau_sec;
btt_amp_mm = Result.Trends.Amp(:);
btt_freq_hz = Result.Trends.Freq(:);
rot_rpm = Result.Trends.RotRPM(:);
half_window_sec = 0.5 * cfg.strain_reference_window_revs .* 60 ./ max(rot_rpm, eps);

num_windows = numel(xv);
strain_amp_raw = nan(num_windows, 1);
strain_amp_mm = nan(num_windows, 1);
strain_freq_hz = nan(num_windows, 1);
samples_used = zeros(num_windows, 1);

for iWin = 1:num_windows
    if ~isfinite(strain_center_time(iWin)) || ~isfinite(btt_freq_hz(iWin))
        continue;
    end
    mask = t_strain >= strain_center_time(iWin) - half_window_sec(iWin) & ...
        t_strain <= strain_center_time(iWin) + half_window_sec(iWin);
    samples_used(iWin) = nnz(mask);
    if samples_used(iWin) < 20
        continue;
    end

    t_local = t_strain(mask) - strain_center_time(iWin);
    v_local = detrend(v_strain(mask));
    freq_grid = build_step6_strain_freq_grid_local(btt_freq_hz(iWin), cfg);
    [strain_amp_raw(iWin), strain_freq_hz(iWin)] = estimate_single_tone_amplitude_local( ...
        t_local, v_local, freq_grid);
    strain_amp_mm(iWin) = strain_amp_raw(iWin) * K_Strain2MM;
end

StrainRef.available = any(isfinite(strain_amp_mm));
StrainRef.AlignmentFile = alignment_file;
StrainRef.StrainFile = alignment.strain_file;
[~, strain_name] = fileparts(alignment.strain_file);
channel_match = regexp(strain_name, '^AI\d+-\d+', 'match', 'once');
StrainRef.ChannelTag = string(channel_match);
StrainRef.K_Strain2MM = K_Strain2MM;
StrainRef.FETransferRatioLabel = string(cfg.fe_transfer_ratio_label);
StrainRef.FETransferRatio_1PerM = cfg.fe_transfer_ratio_1_per_m;
StrainRef.StrainRawUnit = 'microstrain';
StrainRef.AlignmentTauSec = alignment.best_tau_sec;
StrainRef.WindowID = xv;
StrainRef.WindowCenterTimeBTT = btt_center_time;
StrainRef.WindowCenterTimeStrain = strain_center_time;
StrainRef.BTTAmpMM = btt_amp_mm;
StrainRef.BTTFreqHz = btt_freq_hz;
StrainRef.AmpRaw = strain_amp_raw;
StrainRef.AmpTipMM = strain_amp_mm;
StrainRef.FreqHz = strain_freq_hz;
StrainRef.AmpDifferenceMM = btt_amp_mm - strain_amp_mm;
StrainRef.AmpRatioBTTtoStrain = btt_amp_mm ./ strain_amp_mm;
StrainRef.FreqDifferenceHz = btt_freq_hz - strain_freq_hz;
StrainRef.SamplesUsed = samples_used;
StrainRef.FsStrain = Fs_strain;

if StrainRef.available && isfield(cfg, 'export_strain_reference_csv') && ...
        cfg.export_strain_reference_csv && isempty(requested_channel_tag)
    export_step6_strain_reference_csv_local(StrainRef, cfg);
end
end


function strain_file = resolve_strain_channel_peer_file_local(reference_file, channel_tag)
[folder, name, ext] = fileparts(reference_file);
peer_name = regexprep(name, '^AI\d+-\d+', channel_tag);
strain_file = fullfile(folder, [peer_name, ext]);
if ~isfile(strain_file)
    strain_file = '';
end
end


function alignment_file = resolve_step6_strain_alignment_file_local(Result, cfg)
if isfield(cfg, 'strain_alignment_file') && ~isempty(cfg.strain_alignment_file)
    alignment_file = cfg.strain_alignment_file;
    return;
end

case_name = cfg.dynamic_case_name;
if isfield(Result, 'Config') && isfield(Result.Config, 'dynamic_case_name') && ...
        ~isempty(Result.Config.dynamic_case_name)
    case_name = Result.Config.dynamic_case_name;
end
alignment_file = fullfile(cfg.project_dir, 'output', case_name, 'Step3_Spectrum_RPM_Alignment_20251222.mat');
end


function freq_grid = build_step6_strain_freq_grid_local(center_freq_hz, cfg)
half_band = cfg.strain_reference_freq_half_band_hz;
step_hz = cfg.strain_reference_freq_step_hz;
freq_grid = (center_freq_hz - half_band):step_hz:(center_freq_hz + half_band);
freq_grid = freq_grid(isfinite(freq_grid) & freq_grid > 0);
if isempty(freq_grid)
    freq_grid = center_freq_hz;
end
end


function [best_amp, best_freq] = estimate_single_tone_amplitude_local(t_local, v_local, freq_grid)
t_local = t_local(:);
v_local = v_local(:);
best_amp = NaN;
best_freq = NaN;
if isempty(t_local) || isempty(v_local) || isempty(freq_grid)
    return;
end

amp_grid = nan(numel(freq_grid), 1);
for iFreq = 1:numel(freq_grid)
    omega_t = 2 * pi * freq_grid(iFreq) * t_local;
    design = [cos(omega_t), sin(omega_t), ones(size(t_local))];
    coef = design \ v_local;
    amp_grid(iFreq) = hypot(coef(1), coef(2));
end

[best_amp, idx_best] = max(amp_grid);
best_freq = freq_grid(idx_best);
end


function export_step6_strain_reference_csv_local(StrainRef, cfg)
if ~exist(cfg.figure_dir, 'dir')
    mkdir(cfg.figure_dir);
end
    csv_file = fullfile(cfg.figure_dir, sprintf('StrainBTT_B%d_StrainReferenceValidation.csv', cfg.target_blade));
T = table(StrainRef.WindowID, ...
    StrainRef.WindowCenterTimeBTT, ...
    StrainRef.WindowCenterTimeStrain, ...
    repmat(StrainRef.FETransferRatio_1PerM, numel(StrainRef.WindowID), 1), ...
    repmat(StrainRef.K_Strain2MM, numel(StrainRef.WindowID), 1), ...
    StrainRef.BTTAmpMM, ...
    StrainRef.AmpTipMM, ...
    StrainRef.AmpDifferenceMM, ...
    StrainRef.AmpRatioBTTtoStrain, ...
    StrainRef.BTTFreqHz, ...
    StrainRef.FreqHz, ...
    StrainRef.FreqDifferenceHz, ...
    StrainRef.AmpRaw, ...
    StrainRef.SamplesUsed, ...
    'VariableNames', {'WindowID', ...
    'WindowCenterTimeBTT_s', ...
    'WindowCenterTimeStrain_s', ...
    'FETransferRatio_1_per_m', ...
    'StrainToTip_mm_per_microstrain', ...
    'WaveformAmp_mm', ...
    'StrainDerivedAmp_mm', ...
    'AmpDifference_mm', ...
    'AmpRatio_Waveform_to_Strain', ...
    'WaveformFreq_Hz', ...
    'StrainFreq_Hz', ...
    'FreqDifference_Hz', ...
    'StrainRawAmplitude', ...
    'SamplesUsed'});
writetable(T, csv_file);
StrainRef.CSVFile = csv_file; %#ok<NASGU>
end


function export_step6_dual_strain_reference_csv_local(Strain01, Strain03, cfg)
if ~exist(cfg.figure_dir, 'dir')
    mkdir(cfg.figure_dir);
end
    csv_file = fullfile(cfg.figure_dir, sprintf( ...
        'StrainBTT_B%d_StrainChannelComparison.csv', cfg.target_blade));
T = table(Strain01.WindowID, Strain01.BTTAmpMM, Strain01.BTTFreqHz, ...
    Strain01.AmpTipMM, Strain01.FreqHz, ...
    Strain03.AmpTipMM, Strain03.FreqHz, ...
    Strain01.AmpDifferenceMM, Strain03.AmpDifferenceMM, ...
    Strain01.AmpRatioBTTtoStrain, Strain03.AmpRatioBTTtoStrain, ...
    'VariableNames', {'WindowID', 'WaveformAmp_mm', 'WaveformFreq_Hz', ...
    'AI1_01_StrainDerivedAmp_mm', 'AI1_01_StrainFreq_Hz', ...
    'AI1_03_StrainDerivedAmp_mm', 'AI1_03_StrainFreq_Hz', ...
    'WaveformMinusAI1_01_mm', 'WaveformMinusAI1_03_mm', ...
    'WaveformToAI1_01_Ratio', 'WaveformToAI1_03_Ratio'});
writetable(T, csv_file);
end


function print_step6_strain_reference_summary_local(StrainRef)
if ~is_valid_strain_reference_local(StrainRef)
    if isfield(StrainRef, 'message') && ~isempty(StrainRef.message)
        fprintf('  Strain-derived displacement reference: unavailable (%s)\n', StrainRef.message);
    else
        fprintf('  Strain-derived displacement reference: unavailable\n');
    end
    return;
end

channel_tag = "unknown";
if isfield(StrainRef, 'ChannelTag') && strlength(StrainRef.ChannelTag) > 0
    channel_tag = StrainRef.ChannelTag;
end
fprintf(['  Strain-derived displacement reference (%s) | R (%s) = %.9g 1/m | K_Strain2MM = %.9g mm/microstrain\n' ...
    '    strain amp = %.6f-%.6f mm | waveform amp = %.6f-%.6f mm | mean ratio = %.3f\n' ...
    '    strain freq = %.4f-%.4f Hz | waveform-strain mean freq diff = %.4f Hz\n'], ...
    channel_tag, ...
    StrainRef.FETransferRatioLabel, ...
    StrainRef.FETransferRatio_1PerM, ...
    StrainRef.K_Strain2MM, ...
    min(StrainRef.AmpTipMM, [], 'omitnan'), max(StrainRef.AmpTipMM, [], 'omitnan'), ...
    min(StrainRef.BTTAmpMM, [], 'omitnan'), max(StrainRef.BTTAmpMM, [], 'omitnan'), ...
    mean(StrainRef.AmpRatioBTTtoStrain, 'omitnan'), ...
    min(StrainRef.FreqHz, [], 'omitnan'), max(StrainRef.FreqHz, [], 'omitnan'), ...
    mean(StrainRef.FreqDifferenceHz, 'omitnan'));
end


function tf = is_valid_strain_reference_local(StrainRef)
tf = isstruct(StrainRef) && isfield(StrainRef, 'available') && StrainRef.available && ...
    isfield(StrainRef, 'AmpTipMM') && any(isfinite(StrainRef.AmpTipMM));
end


function [t_strain_btt, strain_microstrain, Fs_strain] = load_aligned_strain_trace_local(StrainRef)
loaded_strain = load(StrainRef.StrainFile);
if ~isfield(loaded_strain, 'Datas') || size(loaded_strain.Datas, 2) < 2
    error('Strain file does not contain Datas(:,1:2): %s', StrainRef.StrainFile);
end

t_strain = loaded_strain.Datas(:, 1);
strain_microstrain = loaded_strain.Datas(:, 2);
valid = isfinite(t_strain) & isfinite(strain_microstrain);
t_strain = t_strain(valid);
strain_microstrain = strain_microstrain(valid);

if isfield(loaded_strain, 'SampleFrequency')
    Fs_strain = str2double(loaded_strain.SampleFrequency);
else
    Fs_strain = 1 / median(diff(t_strain));
end
if ~isfinite(Fs_strain) || Fs_strain <= 0
    Fs_strain = 1 / median(diff(t_strain));
end

% Convert the strain time axis to the BTT time base using the Step3 lag.
t_strain_btt = t_strain - StrainRef.AlignmentTauSec;
end


function strain_file = resolve_low_speed_strain_file_local(cfg, StrainRef)
if isfield(cfg, 'low_speed_strain_file') && ~isempty(cfg.low_speed_strain_file)
    strain_file = cfg.low_speed_strain_file;
    if isfile(strain_file)
        return;
    end
end

strain_file = '';
if ~isfield(cfg, 'low_speed_strain_dir') || ~isfolder(cfg.low_speed_strain_dir)
    return;
end

mat_files = dir(fullfile(cfg.low_speed_strain_dir, '*.mat'));
if isempty(mat_files)
    return;
end

preferred_channel = '';
if isfield(StrainRef, 'StrainFile') && ~isempty(StrainRef.StrainFile)
    [~, high_name] = fileparts(StrainRef.StrainFile);
    channel_match = regexp(high_name, '^AI\d+-\d+', 'match', 'once');
    if ~isempty(channel_match)
        preferred_channel = channel_match;
    end
end

file_names = {mat_files.name};
selected_idx = [];
if ~isempty(preferred_channel)
    selected_idx = find(startsWith(file_names, preferred_channel, 'IgnoreCase', true), 1, 'first');
end
if isempty(selected_idx)
    selected_idx = find(startsWith(file_names, 'AI1-03', 'IgnoreCase', true), 1, 'first');
end
if isempty(selected_idx)
    selected_idx = 1;
end
strain_file = fullfile(mat_files(selected_idx).folder, mat_files(selected_idx).name);
end


function [t_strain, strain_microstrain, Fs_strain] = load_strain_trace_file_local(strain_file)
loaded_strain = load(strain_file);
if ~isfield(loaded_strain, 'Datas') || size(loaded_strain.Datas, 2) < 2
    error('Strain file does not contain Datas(:,1:2): %s', strain_file);
end

t_strain = loaded_strain.Datas(:, 1);
strain_microstrain = loaded_strain.Datas(:, 2);
valid = isfinite(t_strain) & isfinite(strain_microstrain);
t_strain = t_strain(valid);
strain_microstrain = strain_microstrain(valid);

if isfield(loaded_strain, 'SampleFrequency')
    Fs_strain = str2double(loaded_strain.SampleFrequency);
else
    Fs_strain = 1 / median(diff(t_strain));
end
if ~isfinite(Fs_strain) || Fs_strain <= 0
    Fs_strain = 1 / median(diff(t_strain));
end
end


function [t_evidence_start, t_evidence_end] = get_step6_evidence_time_span_local(Result, cfg)
center_time = Result.Trends.WindowCenterTime(:);
rot_rpm = Result.Trends.RotRPM(:);
half_window_sec = 0.5 .* cfg.strain_reference_window_revs .* 60 ./ max(rot_rpm, eps);
valid = isfinite(center_time) & isfinite(half_window_sec);
if ~any(valid)
    t_evidence_start = min(center_time, [], 'omitnan');
    t_evidence_end = max(center_time, [], 'omitnan');
    return;
end
t_evidence_start = min(center_time(valid) - half_window_sec(valid));
t_evidence_end = max(center_time(valid) + half_window_sec(valid));
end


function [t_stft, f_stft, amp_stft] = compute_step6_strain_stft_local(t_sec, y_value, Fs_hz, cfg)
t_sec = t_sec(:);
y_value = y_value(:);
valid = isfinite(t_sec) & isfinite(y_value);
t_sec = t_sec(valid);
y_value = y_value(valid);

if numel(t_sec) < 16 || ~isfinite(Fs_hz) || Fs_hz <= 0
    t_stft = NaN;
    f_stft = NaN;
    amp_stft = NaN;
    return;
end

[t_sec, sort_idx] = sort(t_sec);
y_value = y_value(sort_idx);
y_value = detrend(y_value);

window_len = max(16, round(cfg.strain_raw_stft_window_sec .* Fs_hz));
window_len = min(window_len, numel(y_value));
overlap_ratio = min(max(cfg.strain_raw_stft_overlap_ratio, 0), 0.98);
hop_len = max(1, round(window_len .* (1 - overlap_ratio)));
start_idx = 1:hop_len:(numel(y_value) - window_len + 1);
if isempty(start_idx)
    start_idx = 1;
end

nfft = 2 ^ nextpow2(max(window_len, 4096));
freq_all = (0:(nfft / 2))' .* Fs_hz ./ nfft;
freq_mask = freq_all >= cfg.strain_raw_freq_band_hz(1) & ...
    freq_all <= cfg.strain_raw_freq_band_hz(2);
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


function fig = plot_parameter_trends_local(Result, cfg, style, StrainRef)
fig = figure('Name', sprintf('Parameter Trends - Blade %d', cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.parameter_width_cm, style.page.parameter_height_cm], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
xv = Result.Trends.WindowID;

ax = nexttile; hold on;
plot(ax, xv, Result.Trends.Amp, '-o', ...
    'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'Waveform identification');
if is_valid_strain_reference_local(StrainRef)
    plot(ax, xv, StrainRef.AmpTipMM, '-s', ...
        'Color', style.colors.strain_ref, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.colors.strain_ref, ...
        'LineWidth', style.line_aux, ...
        'MarkerSize', style.marker_small, ...
        'DisplayName', 'Strain-derived reference');
    lgd = legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
    set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
end
ylabel('Amplitude (mm)');
title('Amplitude');
panel_label_local(ax, 'a', style);

ax = nexttile; hold on;
stairs(ax, xv, Result.Trends.EO, '-', 'Color', style.colors.accent2, 'LineWidth', style.line_main);
plot(ax, xv, Result.Trends.EO, 'o', ...
    'Color', style.colors.accent2, ...
    'MarkerFaceColor', style.colors.accent2, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.8, ...
    'MarkerSize', style.marker_small);
ylabel('Engine order');
title('Engine Order');
panel_label_local(ax, 'b', style);

ax = nexttile; hold on;
plot_metric_series_local(ax, xv, Result.Trends.Freq, style.colors.accent3, 'o', style);
if is_valid_strain_reference_local(StrainRef)
    plot(ax, xv, StrainRef.FreqHz, '-s', ...
        'Color', style.colors.strain_ref, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', style.colors.strain_ref, ...
        'LineWidth', style.line_aux, ...
        'MarkerSize', style.marker_small);
end
ylabel('Frequency (Hz)');
title('Frequency');
panel_label_local(ax, 'c', style);

ax = nexttile; hold on;
plot_metric_series_local(ax, xv, Result.Trends.Phase, style.colors.accent4, 'o', style);
xlabel('Window index');
ylabel('Phase \phi (rad)');
title('Phase');
panel_label_local(ax, 'd', style);

ax = nexttile; hold on;
plot_metric_series_local(ax, xv, Result.Trends.d0, style.colors.accent5, 'o', style);
xlabel('Window index');
ylabel('Static offset d_0 (mm)');
title('Static Offset');
panel_label_local(ax, 'e', style);

ax = nexttile; hold on;
plot(ax, xv, Result.Trends.WeightedRMSE, '-o', ...
    'Color', style.colors.rmse, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.rmse, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'Weighted RMSE');
plot(ax, xv, Result.Trends.PlainRMSE, '-s', ...
    'Color', style.colors.neutral, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.neutral, ...
    'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, ...
    'DisplayName', 'Plain RMSE');
xlabel('Window index');
ylabel('RMSE (V)');
title('Waveform Error');
lgd = legend(ax, 'Location', 'northeast', 'Box', 'off');
set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
panel_label_local(ax, 'f', style);

axes_all = findall(fig, 'Type', 'axes');
for ia = 1:numel(axes_all)
    apply_nature_axes_style_local(axes_all(ia), style);
end
end


function fig = plot_strain_operating_condition_evidence_local(Result, cfg, style, StrainRef)
low_strain_file = resolve_low_speed_strain_file_local(cfg, StrainRef);
if isempty(low_strain_file) || ~isfile(low_strain_file)
    fig = figure('Name', sprintf('Strain Operating Evidence - Blade %d', cfg.target_blade), ...
        'Color', 'w', ...
        'Units', 'centimeters', ...
        'Position', [1.5, 1.5, style.page.operatingevidence_width_cm, style.page.operatingevidence_height_cm], ...
        'NumberTitle', 'off');
    axis off;
    text(0.5, 0.5, 'Low-speed strain file is unavailable.', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', style.font_name, ...
        'FontSize', style.font_size_text);
    return;
end

[t_low, y_low, ~] = load_strain_trace_file_local(low_strain_file);
[t_high, y_high, Fs_high] = load_aligned_strain_trace_local(StrainRef);
[t_evidence_start, t_evidence_end] = get_step6_evidence_time_span_local(Result, cfg);
y_low_detrended = detrend(y_low);
y_high_detrended = detrend(y_high);

fig = figure('Name', sprintf('Strain Operating Evidence - Blade %d', cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.operatingevidence_width_cm, style.page.operatingevidence_height_cm], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile; hold(ax, 'on');
t_low_rel = t_low - t_low(1);
stride = max(1, floor(numel(t_low_rel) / 50000));
plot(ax, t_low_rel(1:stride:end), y_low_detrended(1:stride:end), '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', 0.45);
xlabel('Time (s)');
ylabel('Strain (\muepsilon)');
title('Detrended low-speed strain');
xlim(ax, [t_low_rel(1), t_low_rel(end)]);
panel_label_local(ax, 'a', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
[freq_low, amp_low] = single_sided_spectrum_local(t_low, y_low_detrended);
plot(ax, freq_low, amp_low, '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', style.line_main);
xlabel('Frequency (Hz)');
ylabel('Amplitude (\muepsilon)');
title('Low-speed spectrum');
xlim(ax, cfg.strain_operating_freq_xlim_hz);
panel_label_local(ax, 'b', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
stride = max(1, floor(numel(t_high) / 80000));
plot(ax, t_high(1:stride:end), y_high_detrended(1:stride:end), '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', 0.35);
xlabel('Aligned time (s)');
ylabel('Strain (\muepsilon)');
title('Detrended target-speed strain');
xlim(ax, [t_high(1), t_high(end)]);
mark_time_region_local(ax, t_evidence_start, t_evidence_end, style.colors.strain_ref);
text(ax, 0.58, 0.90, 'Selected region', ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_colorbar, ...
    'Color', style.colors.strain_ref);
panel_label_local(ax, 'c', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
cfg_stft = cfg;
cfg_stft.strain_raw_freq_band_hz = cfg.strain_operating_freq_xlim_hz;
[t_stft, f_stft, amp_stft] = compute_step6_strain_stft_local(t_high, y_high, Fs_high, cfg_stft);
imagesc(ax, t_stft, f_stft, amp_stft);
axis(ax, 'xy');
colormap(ax, parula);
set_stft_color_limits_local(ax, amp_stft);
cb = colorbar(ax);
cb.Label.String = '\muepsilon';
cb.Label.Interpreter = 'tex';
cb.FontName = style.font_name;
cb.FontSize = style.font_size_colorbar;
cb.Label.FontName = style.font_name;
cb.Label.FontSize = style.font_size_colorbar;
xline(ax, t_evidence_start, '--', 'Color', style.colors.strain_ref, 'LineWidth', 0.9);
xline(ax, t_evidence_end, '--', 'Color', style.colors.strain_ref, 'LineWidth', 0.9);
xlabel('Aligned time (s)');
ylabel('Frequency (Hz)');
title('Target-speed STFT');
xlim(ax, [t_high(1), t_high(end)]);
ylim(ax, cfg.strain_operating_freq_xlim_hz);
panel_label_local(ax, 'd', style);
apply_nature_axes_style_local(ax, style);

fprintf('  Strain operating evidence | low-speed file: %s\n', low_strain_file);
fprintf('    Low-speed detrended strain RMS = %.3f microstrain | target-speed detrended strain RMS = %.3f microstrain\n', ...
    sqrt(mean(y_low_detrended .^ 2, 'omitnan')), sqrt(mean(y_high_detrended .^ 2, 'omitnan')));
end


function mark_time_region_local(ax, t_start, t_end, color_value)
if ~isfinite(t_start) || ~isfinite(t_end) || t_end <= t_start
    return;
end
yl = ylim(ax);
patch(ax, [t_start, t_end, t_end, t_start], [yl(1), yl(1), yl(2), yl(2)], color_value, ...
    'FaceAlpha', 0.12, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
xline(ax, t_start, '--', 'Color', color_value, 'LineWidth', 0.9, 'HandleVisibility', 'off');
xline(ax, t_end, '--', 'Color', color_value, 'LineWidth', 0.9, 'HandleVisibility', 'off');
children = ax.Children;
if numel(children) > 1
    ax.Children = [children(2:end); children(1)];
end
ylim(ax, yl);
end


function set_stft_color_limits_local(ax, amp_stft)
amp_color_values = amp_stft(isfinite(amp_stft));
if isempty(amp_color_values)
    return;
end
amp_color_values = sort(amp_color_values(:));
amp_color_max = amp_color_values(max(1, round(0.995 * numel(amp_color_values))));
clim(ax, [0, max(amp_color_max, eps)]);
end


function fig = plot_strain_raw_evidence_local(Result, cfg, style, StrainRef)
[t_strain_btt, strain_microstrain, Fs_strain] = load_aligned_strain_trace_local(StrainRef);
[t_evidence_start, t_evidence_end] = get_step6_evidence_time_span_local(Result, cfg);
t_plot_start = t_evidence_start - cfg.strain_raw_time_context_sec;
t_plot_end = t_evidence_end + cfg.strain_raw_time_context_sec;
mask_plot = t_strain_btt >= t_plot_start & t_strain_btt <= t_plot_end;
mask_evidence = t_strain_btt >= t_evidence_start & t_strain_btt <= t_evidence_end;
if nnz(mask_plot) < 10
    mask_plot = mask_evidence;
end

fig = figure('Name', sprintf('Raw Strain Evidence - Blade %d', cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.rawevidence_width_cm, style.page.rawevidence_height_cm], ...
    'NumberTitle', 'off');
pos_raw = [0.065, 0.24, 0.250, 0.62];
pos_stft = [0.390, 0.24, 0.205, 0.62];
pos_cbar = [0.608, 0.24, 0.010, 0.62];
pos_spectrum = [0.755, 0.24, 0.215, 0.62];

ax = axes('Parent', fig, 'Position', pos_raw); hold(ax, 'on');
t_plot = t_strain_btt(mask_plot);
y_plot = detrend(strain_microstrain(mask_plot));
stride = max(1, floor(numel(t_plot) / 50000));
plot(ax, t_plot(1:stride:end), y_plot(1:stride:end), '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', 0.40);
xline(ax, t_evidence_start, '--', 'Color', style.colors.neutral, 'LineWidth', 0.9);
xline(ax, t_evidence_end, '--', 'Color', style.colors.neutral, 'LineWidth', 0.9);
xlabel('Aligned time (s)');
ylabel('Strain (\muepsilon)');
title('Detrended Strain');
xlim(ax, [t_plot_start, t_plot_end]);
panel_label_local(ax, 'a', style);
apply_nature_axes_style_local(ax, style);

ax = axes('Parent', fig, 'Position', pos_stft); hold(ax, 'on');
[t_stft, f_stft, amp_stft] = compute_step6_strain_stft_local( ...
    t_strain_btt(mask_plot), strain_microstrain(mask_plot), Fs_strain, cfg);
imagesc(ax, t_stft, f_stft, amp_stft);
axis(ax, 'xy');
colormap(ax, parula);
amp_color_values = amp_stft(isfinite(amp_stft));
if ~isempty(amp_color_values)
    amp_color_values = sort(amp_color_values(:));
    amp_color_max = amp_color_values(max(1, round(0.995 * numel(amp_color_values))));
    clim(ax, [0, max(amp_color_max, eps)]);
end
cb = colorbar(ax);
cb.Label.String = '';
cb.FontName = style.font_name;
cb.FontSize = style.font_size_colorbar;
cb.Position = pos_cbar;
ax.Position = pos_stft;
annotation(fig, 'textbox', [pos_cbar(1) - 0.010, pos_cbar(2) + pos_cbar(4) + 0.018, 0.035, 0.035], ...
    'String', '\mu\epsilon', ...
    'EdgeColor', 'none', ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_colorbar, ...
    'Color', style.colors.neutral_dark, ...
    'FitBoxToText', 'off');
xlabel('Aligned time (s)');
ylabel('Frequency (Hz)');
title('Detrended Strain STFT');
xlim(ax, [t_plot_start, t_plot_end]);
ylim(ax, cfg.strain_raw_freq_band_hz);
panel_label_local(ax, 'b', style);
apply_nature_axes_style_local(ax, style);

ax = axes('Parent', fig, 'Position', pos_spectrum); hold(ax, 'on');
[freq_hz, amp_microstrain] = single_sided_spectrum_local( ...
    t_strain_btt(mask_evidence), strain_microstrain(mask_evidence));
plot(ax, freq_hz, amp_microstrain, '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', style.line_main);
xlabel('Frequency (Hz)');
ylabel('Amplitude (\muepsilon)');
title('Detrended Strain Spectrum');
xlim(ax, cfg.strain_raw_spectrum_xlim_hz);
panel_label_local(ax, 'c', style);
apply_nature_axes_style_local(ax, style);
end


function fig = plot_strain_reference_validation_local(Result, cfg, style, StrainRef)
fig = figure('Name', sprintf('Strain Component Validation - Blade %d', cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.componentvalidation_width_cm, style.page.componentvalidation_height_cm], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xv = Result.Trends.WindowID(:);
component_curve = build_representative_component_curve_local(Result, StrainRef, cfg);

ax = nexttile; hold on;
plot(ax, xv, Result.Trends.Amp(:), '-o', ...
    'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'Waveform identification');
plot(ax, xv, StrainRef.AmpTipMM(:), '-s', ...
    'Color', style.colors.strain_ref, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.strain_ref, ...
    'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, ...
    'DisplayName', 'Strain-derived reference');
xlabel('Window index');
ylabel('Amplitude (mm)');
title('Amplitude');
panel_label_local(ax, 'a', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold on;
plot(ax, xv, Result.Trends.Freq(:), '-o', ...
    'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'Waveform identification');
plot(ax, xv, StrainRef.FreqHz(:), '-s', ...
    'Color', style.colors.strain_ref, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.strain_ref, ...
    'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, ...
    'DisplayName', 'Strain-derived reference');
xlabel('Window index');
ylabel('Frequency (Hz)');
title('Frequency');
panel_label_local(ax, 'b', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile([1 2]); hold(ax, 'on');
if component_curve.available
    h_btt = plot(ax, component_curve.time_ms, component_curve.btt_component_mm, '-', ...
        'Color', style.colors.signal, ...
        'LineWidth', 1.05, ...
        'DisplayName', 'BTT');
    h_strain = plot(ax, component_curve.time_ms, component_curve.strain_component_mm, '--', ...
        'Color', style.colors.strain_ref, ...
        'LineWidth', 1.20, ...
        'DisplayName', 'strain-derived');
    xlabel('Time in selected window (ms)');
    ylabel('Dynamic displacement (mm)');
    title(component_curve.title_text);
    xlim(ax, component_curve.xlim_ms);
    ylim(ax, component_curve.ylim_mm);
    lgd = legend(ax, [h_btt, h_strain], ...
        {'BTT', 'strain-derived'}, ...
        'Location', 'southoutside', ...
        'Orientation', 'horizontal', ...
        'Box', 'off');
    set(lgd, ...
        'FontName', style.font_name, ...
        'FontSize', style.font_size_colorbar, ...
        'ItemTokenSize', [12, 4]);
else
    axis(ax, 'off');
    text(ax, 0.5, 0.5, component_curve.message, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', style.font_name, ...
        'FontSize', style.font_size_text, ...
        'Color', style.colors.neutral_dark);
end
panel_label_wide_local(ax, 'c', style);
apply_nature_axes_style_local(ax, style);
end


function fig = plot_strain_channel_comparison_local(Result, cfg, style, Strain01, Strain03)
xv = Result.Trends.WindowID(:);
btt_amp = Result.Trends.Amp(:);
btt_freq = Result.Trends.Freq(:);
color01 = style.colors.strain_ref;
color03 = style.colors.accent4;

fig = figure('Name', sprintf('Strain Channel Comparison - Blade %d', cfg.target_blade), ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.componentvalidation_width_cm, ...
    style.page.componentvalidation_height_cm], 'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile; hold(ax, 'on');
plot(ax, xv, btt_amp, '-o', 'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, 'DisplayName', 'BTT');
plot(ax, xv, Strain01.AmpTipMM, '-s', 'Color', color01, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'AI1-01');
plot(ax, xv, Strain03.AmpTipMM, '-^', 'Color', color03, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'AI1-03');
xlabel(ax, 'Window index'); ylabel(ax, 'Amplitude (mm)'); title(ax, 'Amplitude');
legend(ax, 'Location', 'best', 'Box', 'off');
panel_label_local(ax, 'a', style); apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, xv, btt_freq, '-o', 'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, 'DisplayName', 'BTT');
plot(ax, xv, Strain01.FreqHz, '-s', 'Color', color01, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'AI1-01');
plot(ax, xv, Strain03.FreqHz, '-^', 'Color', color03, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'AI1-03');
xlabel(ax, 'Window index'); ylabel(ax, 'Frequency (Hz)'); title(ax, 'Frequency');
panel_label_local(ax, 'b', style); apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, xv, btt_amp ./ Strain01.AmpTipMM, '-s', 'Color', color01, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'BTT/AI1-01');
plot(ax, xv, btt_amp ./ Strain03.AmpTipMM, '-^', 'Color', color03, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'BTT/AI1-03');
yline(ax, 1, ':', 'Color', style.colors.neutral, 'HandleVisibility', 'off');
xlabel(ax, 'Window index'); ylabel(ax, 'Amplitude ratio'); title(ax, 'BTT-to-strain ratio');
legend(ax, 'Location', 'best', 'Box', 'off');
panel_label_local(ax, 'c', style); apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, xv, btt_amp - Strain01.AmpTipMM, '-s', 'Color', color01, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'BTT-AI1-01');
plot(ax, xv, btt_amp - Strain03.AmpTipMM, '-^', 'Color', color03, ...
    'MarkerFaceColor', 'w', 'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, 'DisplayName', 'BTT-AI1-03');
yline(ax, 0, ':', 'Color', style.colors.neutral, 'HandleVisibility', 'off');
xlabel(ax, 'Window index'); ylabel(ax, 'Amplitude difference (mm)');
title(ax, 'BTT minus strain-derived');
ax.YAxis.Exponent = 0;
ytickformat(ax, '%.3f');
panel_label_local(ax, 'd', style); apply_nature_axes_style_local(ax, style);
end


function component_curve = build_representative_component_curve_local(Result, StrainRef, cfg)
component_curve = struct();
component_curve.available = false;
component_curve.message = 'Medium-error component curve is unavailable.';

btt_amp = Result.Trends.Amp(:);
btt_freq = Result.Trends.Freq(:);
strain_amp = StrainRef.AmpTipMM(:);
strain_freq = StrainRef.FreqHz(:);
window_id = Result.Trends.WindowID(:);
valid = isfinite(btt_amp) & isfinite(btt_freq) & ...
    isfinite(strain_amp) & isfinite(strain_freq) & ...
    isfinite(window_id);

if ~any(valid)
    component_curve.message = 'No valid BTT/strain component pair.';
    return;
end

valid_idx = find(valid);
amp_rel_error = abs(btt_amp(valid_idx) - strain_amp(valid_idx)) ./ ...
    max(abs(strain_amp(valid_idx)), eps);
freq_abs_error = abs(btt_freq(valid_idx) - strain_freq(valid_idx));
amp_score = normalize_error_series_local(amp_rel_error);
freq_score = normalize_error_series_local(freq_abs_error);
validation_score = amp_score + freq_score;
[~, score_order] = sort(validation_score, 'ascend');
selected_idx = valid_idx(score_order(ceil(numel(score_order) / 2)));

num_samples = cfg.strain_component_time_samples;
if ~isfinite(num_samples) || num_samples < 200
    num_samples = 600;
end
num_samples = round(num_samples);
half_width_ms = cfg.strain_component_window_half_width_ms;
if ~isfinite(half_width_ms) || half_width_ms <= 0
    half_width_ms = 5.0;
end
time_ms = linspace(-half_width_ms, half_width_ms, num_samples).';
time_s = time_ms .* 1e-3;

btt_component_mm = btt_amp(selected_idx) .* cos(2 .* pi .* btt_freq(selected_idx) .* time_s);
strain_component_mm = strain_amp(selected_idx) .* cos(2 .* pi .* strain_freq(selected_idx) .* time_s);

all_y = [btt_component_mm(:); strain_component_mm(:)];
max_abs_y = max(abs(all_y), [], 'omitnan');
if ~isfinite(max_abs_y) || max_abs_y <= 0
    max_abs_y = 1;
end

component_curve.available = true;
component_curve.message = '';
component_curve.time_ms = time_ms;
component_curve.btt_component_mm = btt_component_mm;
component_curve.strain_component_mm = strain_component_mm;
component_curve.window_id = window_id(selected_idx);
component_curve.title_text = sprintf('Median-error Window %d EO14 Component', window_id(selected_idx));
component_curve.xlim_ms = [-half_width_ms, half_width_ms];
component_curve.ylim_mm = 1.08 .* [-max_abs_y, max_abs_y];
component_curve.amp_rel_error = amp_rel_error(score_order(ceil(numel(score_order) / 2)));
component_curve.freq_abs_error = freq_abs_error(score_order(ceil(numel(score_order) / 2)));
end


function score = normalize_error_series_local(error_value)
error_value = error_value(:);
finite_value = error_value(isfinite(error_value));
if isempty(finite_value)
    score = zeros(size(error_value));
    return;
end
span_value = max(finite_value) - min(finite_value);
if span_value <= eps
    score = zeros(size(error_value));
else
    score = (error_value - min(finite_value)) ./ span_value;
end
score(~isfinite(score)) = 0;
end


function fig = plot_strain_btt_diagnostic_figure_local(Result, cfg, style, StrainRef)
diag = build_step6_strain_btt_diagnostic_local(Result, cfg, StrainRef);

fig = figure('Name', sprintf('Strain-BTT Diagnostics - Blade %d', cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.diagnostic_width_cm, style.page.diagnostic_height_cm], ...
    'NumberTitle', 'off');
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xv = Result.Trends.WindowID(:);

if ~diag.available
    ax = nexttile([2 2]);
    axis(ax, 'off');
    text(ax, 0.5, 0.5, diag.message, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', style.font_name, ...
        'FontSize', style.font_size_text, ...
        'Color', style.colors.neutral_dark);
    return;
end

ax = nexttile; hold(ax, 'on');
h_btt_amp = plot(ax, xv, Result.Trends.Amp(:), '-o', ...
    'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'BTT');
h_strain_amp = plot(ax, xv, StrainRef.AmpTipMM(:), '-s', ...
    'Color', style.colors.strain_ref, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.strain_ref, ...
    'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, ...
    'DisplayName', 'Strain-derived');
xline(ax, diag.best_window_id, ':', 'Color', style.colors.neutral, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel('Window index');
ylabel('Amplitude (mm)');
title('Amplitude tracking');
lgd = legend(ax, [h_btt_amp, h_strain_amp], {'BTT', 'Strain-derived'}, ...
    'Location', 'northwest', 'Orientation', 'horizontal', 'Box', 'off');
set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
panel_label_local(ax, 'a', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, xv, Result.Trends.Freq(:), '-o', ...
    'Color', style.colors.signal, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main, ...
    'DisplayName', 'BTT');
plot(ax, xv, StrainRef.FreqHz(:), '-s', ...
    'Color', style.colors.strain_ref, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', style.colors.strain_ref, ...
    'LineWidth', style.line_aux, ...
    'MarkerSize', style.marker_small, ...
    'DisplayName', 'Strain');
xline(ax, diag.best_window_id, ':', 'Color', style.colors.neutral, ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel('Window index');
ylabel('Frequency (Hz)');
title('Frequency tracking');
panel_label_local(ax, 'b', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, diag.spectrum_freq_hz, diag.spectrum_amp_mm, '-', ...
    'Color', style.colors.neutral_dark, ...
    'LineWidth', style.line_aux, ...
    'DisplayName', 'Strain spectrum');
xline(ax, diag.btt_freq_hz, '--', 'Color', style.colors.signal, 'LineWidth', 1.2, ...
    'DisplayName', 'BTT f');
xline(ax, diag.strain_freq_hz, '--', 'Color', style.colors.strain_ref, 'LineWidth', 1.2, ...
    'DisplayName', 'Strain f');
xlim(ax, [diag.btt_freq_hz - cfg.strain_diagnostic_spectrum_half_band_hz, ...
    diag.btt_freq_hz + cfg.strain_diagnostic_spectrum_half_band_hz]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (mm)');
title('Strain-derived spectrum');
lgd = legend(ax, 'Location', 'northeast', 'Box', 'off');
set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
panel_label_local(ax, 'c', style);
apply_nature_axes_style_local(ax, style);

ax = nexttile; hold(ax, 'on');
plot(ax, diag.time_strain_ms, diag.u_strain_raw_mm, '-', ...
    'Color', [0.78 0.78 0.78], ...
    'LineWidth', 0.7, ...
    'HandleVisibility', 'off');
h_strain_time = plot(ax, diag.time_strain_lagged_ms, diag.u_strain_fit_lagged_mm, '-', ...
    'Color', style.colors.strain_ref, ...
    'LineWidth', style.line_main, ...
    'DisplayName', 'Strain fit');
h_btt_time = plot(ax, diag.time_btt_fit_ms, diag.u_btt_fit_mm, '--', ...
    'Color', style.colors.signal, ...
    'LineWidth', style.line_main, ...
    'DisplayName', 'BTT fit');
xlabel('Time in best window (ms)');
ylabel('Dynamic displacement (mm)');
title('Time-domain check');
xlim(ax, [-cfg.strain_diagnostic_time_half_width_ms, cfg.strain_diagnostic_time_half_width_ms]);
text(ax, 0.03, 0.94, sprintf('r_0 = %.3f\nr = %.3f\nlag = %.3f ms', ...
    diag.corr_fixed, diag.corr_lagged, diag.best_lag_ms), ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_text, ...
    'Color', style.colors.neutral_dark, ...
    'BackgroundColor', [1 1 1 0.76], ...
    'Margin', 4);
lgd = legend(ax, [h_btt_time, h_strain_time], {'BTT fit', 'Strain fit'}, ...
    'Location', 'southeast', 'Orientation', 'horizontal', 'Box', 'off');
set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
panel_label_local(ax, 'd', style);
apply_nature_axes_style_local(ax, style);

fprintf(['  Strain-BTT diagnostic | best window = %d | amp ratio = %.3f | ' ...
    'freq diff = %.4f Hz | fixed corr = %.3f | lagged corr = %.3f | lag = %.4f ms\n'], ...
    diag.best_window_id, diag.amp_ratio_btt_to_strain, ...
    diag.freq_difference_hz, diag.corr_fixed, diag.corr_lagged, diag.best_lag_ms);
end


function diag = build_step6_strain_btt_diagnostic_local(Result, ~, StrainRef)
diag = struct();
diag.available = false;
diag.message = '';

if ~isfield(Result, 'BestWindow') || ~isfield(Result.BestWindow, 'bundle')
    diag.message = 'Best-window BTT data are unavailable.';
    return;
end
if ~isfield(StrainRef, 'StrainFile') || ~isfile(StrainRef.StrainFile)
    diag.message = 'Strain file is unavailable.';
    return;
end

B = Result.BestWindow;
bundle = B.bundle;
t_btt = bundle.T(:);
if isfield(B, 'u_est') && numel(B.u_est) == numel(t_btt)
    u_btt = B.u_est(:);
else
    t0 = mean(t_btt, 'omitnan');
    u_btt = B.d0_id + B.A_id .* cos(2 * pi * B.fn_id .* (t_btt - t0) + B.phi_id_wrapped);
end
valid_btt = isfinite(t_btt) & isfinite(u_btt);
t_btt = t_btt(valid_btt);
u_btt = u_btt(valid_btt);
if numel(t_btt) < 10
    diag.message = 'Not enough BTT points for diagnostic comparison.';
    return;
end

loaded_strain = load(StrainRef.StrainFile);
if ~isfield(loaded_strain, 'Datas') || size(loaded_strain.Datas, 2) < 2
    diag.message = 'Strain file does not contain Datas(:,1:2).';
    return;
end
t_strain = loaded_strain.Datas(:, 1);
v_strain = loaded_strain.Datas(:, 2);
t_strain_btt = t_strain - StrainRef.AlignmentTauSec;

t_lo = min(t_btt);
t_hi = max(t_btt);
mask_strain = t_strain_btt >= t_lo & t_strain_btt <= t_hi;
if nnz(mask_strain) < 20
    diag.message = 'Not enough aligned strain samples in the best BTT window.';
    return;
end

t_strain_seg = t_strain_btt(mask_strain);
strain_microstrain = detrend(v_strain(mask_strain));
u_strain_raw = strain_microstrain .* StrainRef.K_Strain2MM;
[~, ~, ~, u_strain_fit] = fit_single_tone_series_local( ...
    t_strain_seg - mean(t_strain_seg, 'omitnan'), u_strain_raw, B.fn_id);

[t_btt_sort, idx_sort] = sort(t_btt);
u_btt_sort = u_btt(idx_sort);
[t_btt_unique, idx_unique] = unique(t_btt_sort, 'stable');
u_btt_unique = u_btt_sort(idx_unique);
u_btt_dynamic = u_btt_unique - mean(u_btt_unique, 'omitnan');
t_center = mean([t_lo, t_hi]);
t_btt_fit = linspace(t_lo, t_hi, 2000).';
[~, ~, ~, ~, coef_btt] = fit_single_tone_series_local( ...
    t_btt_unique - t_center, u_btt_dynamic, B.fn_id);
u_btt_fit = evaluate_single_tone_series_local(t_btt_fit - t_center, B.fn_id, coef_btt);
u_btt_fit = u_btt_fit - mean(u_btt_fit, 'omitnan');
u_btt_on_strain = interp1(t_btt_fit, u_btt_fit, t_strain_seg, 'linear', 'extrap');
polarity = choose_polarity_by_correlation_local(u_btt_on_strain, u_strain_fit);
u_strain_raw = polarity .* (u_strain_raw - mean(u_strain_raw, 'omitnan'));
u_strain_fit = polarity .* (u_strain_fit - mean(u_strain_fit, 'omitnan'));
u_btt_on_strain = u_btt_on_strain - mean(u_btt_on_strain, 'omitnan');
[best_lag_sec, corr_lagged, rmse_lagged] = find_best_lag_correlation_local( ...
    t_btt_fit, u_btt_fit, t_strain_seg, u_strain_fit, B.fn_id);

[freq_hz, amp_mm] = single_sided_spectrum_local(t_strain_seg, u_strain_raw);

best_idx = find(StrainRef.WindowID == B.window_id, 1, 'first');
if isempty(best_idx)
    best_idx = max(1, min(numel(StrainRef.WindowID), B.window_id));
end
strain_freq_hz = StrainRef.FreqHz(best_idx);
if ~isfinite(strain_freq_hz)
    strain_freq_hz = B.fn_id;
end

diag.available = true;
diag.best_window_id = B.window_id;
diag.btt_freq_hz = B.fn_id;
diag.strain_freq_hz = strain_freq_hz;
diag.freq_difference_hz = B.fn_id - strain_freq_hz;
diag.amp_ratio_btt_to_strain = StrainRef.BTTAmpMM(best_idx) ./ StrainRef.AmpTipMM(best_idx);
diag.corr_fixed = safe_corr_local(u_btt_on_strain, u_strain_fit);
diag.rmse_fixed_mm = sqrt(mean((u_btt_on_strain - u_strain_fit) .^ 2, 'omitnan'));
diag.best_lag_ms = best_lag_sec .* 1000;
diag.corr_lagged = corr_lagged;
diag.rmse_lagged_mm = rmse_lagged;
diag.spectrum_freq_hz = freq_hz(:);
diag.spectrum_amp_mm = amp_mm(:);
diag.time_btt_sample_ms = (t_btt_unique - t_center) .* 1000;
diag.u_btt_sample_mm = u_btt_dynamic;
diag.time_btt_fit_ms = (t_btt_fit - t_center) .* 1000;
diag.u_btt_fit_mm = u_btt_fit;
diag.time_strain_ms = (t_strain_seg - t_center) .* 1000;
diag.u_strain_raw_mm = u_strain_raw;
diag.u_strain_fit_mm = u_strain_fit;
diag.time_strain_lagged_ms = (t_strain_seg + best_lag_sec - t_center) .* 1000;
diag.u_strain_fit_lagged_mm = u_strain_fit;
end


function [amp_value, phase_rad, offset_value, y_fit, coef] = fit_single_tone_series_local(t_local, y_value, freq_hz)
t_local = t_local(:);
y_value = y_value(:);
omega_t = 2 * pi * freq_hz .* t_local;
design = [cos(omega_t), sin(omega_t), ones(size(t_local))];
coef = design \ y_value;
amp_value = hypot(coef(1), coef(2));
phase_rad = atan2(coef(2), coef(1));
offset_value = coef(3);
y_fit = design * coef;
end


function y_fit = evaluate_single_tone_series_local(t_local, freq_hz, coef)
t_local = t_local(:);
omega_t = 2 * pi * freq_hz .* t_local;
design = [cos(omega_t), sin(omega_t), ones(size(t_local))];
y_fit = design * coef;
end


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
coherent_gain = mean(win);
nfft = 2 ^ nextpow2(max(n, 4096));
y_fft = fft(y_value .* win, nfft);
amp_value = abs(y_fft(1:nfft / 2 + 1)) ./ n ./ max(coherent_gain, eps) .* 2;
amp_value(1) = amp_value(1) ./ 2;
freq_hz = (0:(nfft / 2))' .* fs_hz ./ nfft;
end


function win = hann_window_local(n)
if n <= 1
    win = ones(n, 1);
    return;
end
idx = (0:(n - 1))';
win = 0.5 - 0.5 .* cos(2 .* pi .* idx ./ (n - 1));
end


function polarity = choose_polarity_by_correlation_local(reference_signal, test_signal)
corr_pos = safe_corr_local(reference_signal, test_signal);
corr_neg = safe_corr_local(reference_signal, -test_signal);
if isfinite(corr_neg) && (~isfinite(corr_pos) || corr_neg > corr_pos)
    polarity = -1;
else
    polarity = 1;
end
end


function [best_lag_sec, best_corr, best_rmse] = find_best_lag_correlation_local( ...
    t_reference, y_reference, t_test, y_test, freq_hz)
half_period_sec = 0.5 ./ max(freq_hz, eps);
lag_grid_sec = linspace(-half_period_sec, half_period_sec, 401);
best_lag_sec = 0;
best_corr = -Inf;
best_rmse = Inf;
for iLag = 1:numel(lag_grid_sec)
    lag_sec = lag_grid_sec(iLag);
    y_shifted = interp1(t_test + lag_sec, y_test, t_reference, 'linear', NaN);
    r_value = safe_corr_local(y_reference, y_shifted);
    if ~isfinite(r_value)
        continue;
    end
    valid = isfinite(y_reference) & isfinite(y_shifted);
    rmse_value = sqrt(mean((y_reference(valid) - y_shifted(valid)) .^ 2, 'omitnan'));
    if r_value > best_corr
        best_corr = r_value;
        best_lag_sec = lag_sec;
        best_rmse = rmse_value;
    end
end
if isinf(best_corr)
    best_corr = NaN;
    best_rmse = NaN;
end
end


function r_value = safe_corr_local(x_value, y_value)
x_value = x_value(:);
y_value = y_value(:);
valid = isfinite(x_value) & isfinite(y_value);
x_value = x_value(valid);
y_value = y_value(valid);
if numel(x_value) < 3 || std(x_value) <= eps || std(y_value) <= eps
    r_value = NaN;
    return;
end
c = corrcoef(x_value, y_value);
r_value = c(1, 2);
end


function fig = plot_window_time_compare_figure_local(Result, cfg, mode_name, style)
if strcmp(mode_name, 'BestWindow')
    bundle = Result.BestWindow.bundle;
    v_pred = Result.BestWindow.V_pred;
    figure_name_tag = 'Best Window';
else
    rep_window_id = Result.RepresentativeWindow.window_id;
    bundle = reconstruct_bundle_from_window_local(Result, rep_window_id);
    v_pred = reconstruct_prediction_from_window_local(Result, rep_window_id);
    figure_name_tag = 'Representative Window';
end

sensor_ids = unique(bundle.S(:)).';
fig = figure('Name', sprintf('%s Time Compare - Blade %d', figure_name_tag, cfg.target_blade), ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [1.5, 1.5, style.page.timecompare_width_cm, style.page.timecompare_height_cm], ...
    'NumberTitle', 'off');
tiledlayout(fig, 1, numel(sensor_ids), 'TileSpacing', 'compact', 'Padding', 'compact');
legend_handles = gobjects(1, 2);
legend_axis = [];
for i = 1:numel(sensor_ids)
    sid = sensor_ids(i);
    ax = nexttile;
    hold(ax, 'on');
    mask = bundle.S == sid;
    [t_sid, idx_sid] = sort(bundle.T(mask));
    v_meas = bundle.V(mask);
    v_pred_sid = v_pred(mask);
    [t_ms, v_meas_sorted, v_pred_sorted] = build_zoomed_waveform_view_local( ...
        t_sid, v_meas(idx_sid), v_pred_sid(idx_sid), cfg.zoom_half_width_ms);
    residual = v_meas_sorted - v_pred_sorted;
    sensor_rmse = sqrt(mean(residual .^ 2, 'omitnan'));

    h_meas = plot(ax, t_ms, v_meas_sorted, '-', ...
        'Color', style.colors.neutral_dark, ...
        'LineWidth', style.line_aux, ...
        'DisplayName', 'Measured waveform');
    h_pred = plot(ax, t_ms, v_pred_sorted, '-', ...
        'Color', style.colors.signal, ...
        'LineWidth', style.line_main, ...
        'DisplayName', 'Model prediction');
    if i == 1
        legend_handles = [h_meas, h_pred];
        legend_axis = ax;
    end

    yl = ylim(ax);
    text(ax, 0.02, 0.95, sprintf('CH%d\nN = %d\nRMSE = %.4f V', sid, nnz(mask), sensor_rmse), ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'FontName', style.font_name, ...
        'FontSize', style.font_size_text, ...
        'Color', style.colors.neutral_dark, ...
        'BackgroundColor', [1 1 1 0.72], ...
        'Margin', 4);
    ylim(ax, yl);
    xlabel('Time relative to peak (ms)');
    if i == 1
        ylabel('Voltage (V)');
    end
    title(sprintf('Sensor %d', sid), 'FontName', style.font_name, ...
        'FontSize', style.font_size_title - 1, 'FontWeight', 'bold');
    panel_label_local(ax, char('a' + i - 1), style);
    apply_nature_axes_style_local(ax, style);
end

lgd = legend(legend_axis, legend_handles, {'Measured waveform', 'Model prediction'}, ...
    'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');
set(lgd, 'FontName', style.font_name, 'FontSize', style.font_size_legend);
end


function bundle = reconstruct_bundle_from_window_local(Result, window_id)
if Result.BestWindow.window_id == window_id
    bundle = Result.BestWindow.bundle;
    return;
end
error(['This visualization currently reuses the bundle stored in BestWindow.\n' ...
    'RepresentativeWindow differs from BestWindow in this result file,\n' ...
    'but no standalone bundle was stored for it. Regenerate the identification result with bundle persistence if needed.']);
end


function v_pred = reconstruct_prediction_from_window_local(Result, window_id)
if Result.BestWindow.window_id == window_id
    v_pred = Result.BestWindow.V_pred;
    return;
end
error(['This visualization currently reuses the waveform prediction stored in BestWindow.\n' ...
    'RepresentativeWindow differs from BestWindow in this result file,\n' ...
    'but no standalone prediction was stored for it. Regenerate the identification result with prediction persistence if needed.']);
end


function save_figure_step6_local(fig, figure_dir, blade_id, tag, cfg)
if isempty(figure_dir) || ~ishandle(fig)
    return;
end
png_file = fullfile(figure_dir, sprintf('StrainBTT_B%d_%s.png', blade_id, tag));
exportgraphics(fig, png_file, 'Resolution', 600);
if isfield(cfg, 'export_pdf') && cfg.export_pdf
    pdf_file = fullfile(figure_dir, sprintf('StrainBTT_B%d_%s.pdf', blade_id, tag));
    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
end
end


function sensor_tag = format_sensor_tag_local(sensor_ids)
sensor_ids = sensor_ids(:).';
parts = compose('%d', sensor_ids);
sensor_tag = ['S', strjoin(cellstr(parts), '')];
end


function time_tag = format_time_tag_local(value_sec)
txt = num2str(value_sec, '%.6g');
txt = strrep(txt, '.', 'p');
txt = strrep(txt, '-', 'm');
time_tag = txt;
end


function style = build_step6_style_local(cfg)
style = struct();
style.mode = 'nature';
if isfield(cfg, 'figure_style') && ~isempty(cfg.figure_style)
    style.mode = lower(string(cfg.figure_style));
end
style.font_name = choose_font_local({'Times New Roman', 'Arial', 'Helvetica'});
style.font_size_axis = 9.0;
style.font_size_text = 9.0;
style.font_size_legend = 9.0;
style.font_size_title = 9.0;
style.font_size_panel = 9.0;
style.font_size_colorbar = 7.5;
style.font_size_rawlegend = 5.8;
style.line_main = 1.35;
style.line_aux = 1.0;
style.marker_main = 4.8;
style.marker_small = 4.0;
style.axis_line_width = 0.8;
style.grid_alpha = 0.12;
style.colors = struct( ...
    'signal', [0.00 0.38 0.55], ...
    'accent2', [0.58 0.38 0.00], ...
    'accent3', [0.15 0.53 0.33], ...
    'accent4', [0.51 0.32 0.58], ...
    'accent5', [0.08 0.47 0.64], ...
    'rmse', [0.72 0.18 0.20], ...
    'strain_ref', [0.80 0.22 0.08], ...
    'neutral', [0.45 0.45 0.45], ...
    'neutral_dark', [0.12 0.12 0.12], ...
    'grid', [0.87 0.87 0.87], ...
    'best', [0.76 0.20 0.18], ...
    'representative', [0.25 0.25 0.25]);
style.page = struct( ...
    'parameter_width_cm', 17.0, ...
    'parameter_height_cm', 11.6, ...
    'operatingevidence_width_cm', 17.0, ...
    'operatingevidence_height_cm', 9.0, ...
    'rawevidence_width_cm', 17.0, ...
    'rawevidence_height_cm', 5.4, ...
    'componentvalidation_width_cm', 17.0, ...
    'componentvalidation_height_cm', 10.8, ...
    'diagnostic_width_cm', 17.0, ...
    'diagnostic_height_cm', 11.2, ...
    'timecompare_width_cm', 17.0, ...
    'timecompare_height_cm', 6.8);
end


function font_name = choose_font_local(candidates)
available = listfonts;
font_name = 'Times New Roman';
for i = 1:numel(candidates)
    if any(strcmpi(available, candidates{i}))
        font_name = candidates{i};
        return;
    end
end
end


function plot_metric_series_local(ax, xv, yv, color_value, marker_symbol, style)
plot(ax, xv, yv, ['-' marker_symbol], ...
    'Color', color_value, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color_value, ...
    'LineWidth', style.line_main, ...
    'MarkerSize', style.marker_main);
end


function apply_nature_axes_style_local(ax, style)
set(ax, ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_axis, ...
    'LineWidth', style.axis_line_width, ...
    'Box', 'on', ...
    'TickDir', 'in', ...
    'TickLength', [0.018 0.018], ...
    'XColor', style.colors.neutral_dark, ...
    'YColor', style.colors.neutral_dark, ...
    'Layer', 'top');
grid(ax, 'off');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.Title.FontName = style.font_name;
ax.Title.FontSize = style.font_size_title;
ax.Title.FontWeight = 'normal';
ax.XLabel.FontName = style.font_name;
ax.XLabel.FontSize = style.font_size_axis;
ax.YLabel.FontName = style.font_name;
ax.YLabel.FontSize = style.font_size_axis;
end


function panel_label_local(ax, label_char, style)
text(ax, -0.16, 1.02, sprintf('(%s)', label_char), ...
    'Units', 'normalized', ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_panel, ...
    'FontWeight', 'normal', ...
    'Color', style.colors.neutral_dark, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end


function panel_label_wide_local(ax, label_char, style)
text(ax, 0.00, 1.02, sprintf('(%s)', label_char), ...
    'Units', 'normalized', ...
    'FontName', style.font_name, ...
    'FontSize', style.font_size_panel, ...
    'FontWeight', 'normal', ...
    'Color', style.colors.neutral_dark, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end


function [t_zoom_ms, v_meas_zoom, v_pred_zoom] = build_zoomed_waveform_view_local(t_sec, v_meas, v_pred, half_width_ms)
t_sec = t_sec(:);
v_meas = v_meas(:);
v_pred = v_pred(:);

peak_time = choose_pulse_peak_time_local(t_sec, v_meas);
half_width_sec = half_width_ms * 1e-3;
mask_zoom = t_sec >= (peak_time - half_width_sec) & t_sec <= (peak_time + half_width_sec);
if nnz(mask_zoom) < 10
    [~, idx_peak] = max(v_meas);
    idx_lo = max(1, idx_peak - 1500);
    idx_hi = min(numel(t_sec), idx_peak + 1500);
    mask_zoom = false(size(t_sec));
    mask_zoom(idx_lo:idx_hi) = true;
end

t_zoom = t_sec(mask_zoom);
v_meas_zoom = v_meas(mask_zoom);
v_pred_zoom = v_pred(mask_zoom);
t_zoom_ms = (t_zoom - peak_time) * 1e3;
end


function peak_time = choose_pulse_peak_time_local(t_sec, v_meas)
[peak_val, idx_global] = max(v_meas);
baseline_est = median(v_meas, 'omitnan');
threshold = baseline_est + 0.6 * max(peak_val - baseline_est, eps);
above = v_meas >= threshold;
edges = diff([false; above; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;

if isempty(starts)
    peak_time = t_sec(idx_global);
    return;
end

mid_t = 0.5 * (t_sec(1) + t_sec(end));
peak_times = zeros(numel(starts), 1);
peak_scores = zeros(numel(starts), 1);
for k = 1:numel(starts)
    idx_range = starts(k):stops(k);
    [local_peak, idx_rel] = max(v_meas(idx_range));
    peak_idx = idx_range(idx_rel);
    peak_times(k) = t_sec(peak_idx);
    peak_scores(k) = local_peak;
end

[~, best_group] = min(abs(peak_times - mid_t) - 0.02 * peak_scores / max(peak_scores));
peak_time = peak_times(best_group);
end
