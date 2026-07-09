%% Step06SG_VisualizeStrainGaugeAgainstIdentification_20241106
% Show strain-gauge data directly and compare it with Step06 identification.
% This script is display-first: figures open in MATLAB and are not saved.

clear; close all; clc;

thisFileDir = fileparts(mfilename('fullpath'));
routeDir = fileparts(thisFileDir);
addpath(routeDir);

P = NewFlow_Config_20241106();

% Local parameter block. Edit here directly.
SSG = struct();
SSG.dataRoot = 'D:\鍗氬＋-鍥界\璇曢獙鍙版暟鎹甛鏂拌瘯楠孿20241106_2';
SSG.highDir = fullfile(SSG.dataRoot, '3000_3150');
SSG.strainFiles = { ...
    'AI1-01_20241106170148.mat'
    'AI1-02_20241106170148.mat'
    'AI1-03_20241106170148.mat'
    'AI1-04_20241106170148.mat'};
SSG.strainLabels = {'AI1-01','AI1-02','AI1-03','AI1-04'};
SSG.selectedChannels = [4];          % e.g. [1 2 3 4]
SSG.detailChannel = 4;               % used for STFT + Step06 overlay
SSG.strainSampleRateHz = 10000;
SSG.strainScaleToMicroStrain = 1;    % set here if Datas(:,2) is not already in microstrain
SSG.strainUnitText = 'Microstrain';
SSG.strainUnitSymbol = '\muepsilon';
SSG.strainToBttOffsetSec = -102.6;
SSG.detrendMethod = 'linear';        % 'none' | 'linear' | 'movingmedian'
SSG.stftWindowSec = 2.0;
SSG.stftStepSec = 0.25;
SSG.stftFreqBandHz = [0 1000];
SSG.sgPeakSearchHalfBandHz = 20;
SSG.analysisSensors = [2 5 7];
SSG.startTimeSec = 75.0;
SSG.outputLabel = 'B4_only';
SSG.bladeId = 4;
SSG.showStep06Compare = true;
SSG.showOprSpeedOverview = true;
SSG.overviewTimeMarginSec = 1.0;
SSG.showStrain3D = true;
SSG.spec2DTimeRangeSec = [72 76];   % 2D time-frequency view range
SSG.spec1DTimeRangeSec = [74.75 75.25]; % 1D spectrum accumulation range
SSG.saveFigures = false;

P = apply_step07sg_local_options_local(P, SSG);

StrainGaugeData = repmat(empty_strain_struct_local(), numel(SSG.selectedChannels), 1); %#ok<NASGU>
for i = 1:numel(SSG.selectedChannels)
    ch = SSG.selectedChannels(i);
    StrainGaugeData(i) = load_strain_local( ...
        fullfile(SSG.dataRoot, SSG.strainFiles{ch}), ...
        SSG.strainLabels{ch}, ...
        SSG.strainSampleRateHz, ...
        SSG.detrendMethod, ...
        SSG.strainToBttOffsetSec, ...
        SSG.strainScaleToMicroStrain, ...
        SSG.strainUnitText, ...
        SSG.strainUnitSymbol);
end

IdentificationResult = [];
Step06StrainCompareTable = table(); %#ok<NASGU>
compareBladeId = SSG.bladeId;
if SSG.showStep06Compare && exist(P.files.identificationResult, 'file') == 2
    loaded = load(P.files.identificationResult, 'IdentificationResult');
    IdentificationResult = loaded.IdentificationResult;
    compareBladeId = resolve_compare_blade_local(IdentificationResult, SSG.bladeId);
    Step06StrainCompareTable = build_step06_strain_compare_table_local( ...
        IdentificationResult, StrainGaugeData, SSG, compareBladeId);
end

OprProfile = load_opr_speed_profile_local(P);
plot_opr_speed_and_identification_overview_local( ...
    OprProfile, IdentificationResult, Step06StrainCompareTable, SSG, compareBladeId);
plot_strain_time_series_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId);
plot_strain_stft_3d_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId);
plot_strain_stft_with_step06_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId);
plot_strain_local_spectrogram_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId);
plot_strain_local_frequency_spectrum_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId);
plot_strain_frequency_compare_local(Step06StrainCompareTable, SSG, compareBladeId);

fprintf('\n=== Step06SG: strain gauge vs identification ===\n');
fprintf('Selected channels: %s\n', mat2str(SSG.selectedChannels));
fprintf('Detail channel: %s\n', SSG.strainLabels{SSG.detailChannel});
fprintf('strainToBttOffsetSec = %+0.3f s\n', SSG.strainToBttOffsetSec);
if isempty(IdentificationResult)
    fprintf('Step06 result was not loaded. Only strain-gauge figures are shown.\n');
else
    fprintf('Step06 source: %s\n', P.files.identificationResult);
    fprintf('Compare blade: B%d\n', compareBladeId);
    disp(Step06StrainCompareTable);
end
fprintf('Figures are open in MATLAB. No image file was saved.\n');

function P = apply_step07sg_local_options_local(P, SSG)
P.sensors.analysis = SSG.analysisSensors(:).';
P.region.startTimeSec = SSG.startTimeSec;
P.view.saveFigures = logical(SSG.saveFigures);

sensorTag = sensor_tag_local(P.sensors.analysis);
timeLabel = time_label_local(P.region.startTimeSec);
suffix = output_suffix_local(SSG.outputLabel);
P.files.identificationResult = fullfile(P.outputDir, '06_identification', sensorTag, ...
    sprintf('IdentificationResult_%s_20241106%s.mat', timeLabel, suffix));
end

function S = empty_strain_struct_local()
S = struct('channel_id', NaN, 'label', "", 'file', "", ...
    'time_strain_sec', [], 'time_btt_sec', [], ...
    'raw_value', [], 'value', [], 'detrend_method', "", ...
    'unit_text', "", 'unit_symbol', "");
end

function strain = load_strain_local(filePath, label, Fs, detrendMethod, offsetSec, scaleToMicroStrain, unitText, unitSymbol)
if exist(filePath, 'file') ~= 2
    error('Missing strain file: %s', filePath);
end
S = load(filePath, 'Datas');
strain = empty_strain_struct_local();
strain.label = string(label);
strain.file = string(filePath);
strain.time_strain_sec = S.Datas(:, 1);
strain.time_btt_sec = S.Datas(:, 1) + offsetSec;
strain.raw_value = scaleToMicroStrain * S.Datas(:, 2);
strain.value = preprocess_strain_local(strain.raw_value, Fs, detrendMethod);
strain.detrend_method = string(detrendMethod);
strain.unit_text = string(unitText);
strain.unit_symbol = string(unitSymbol);
end

function y = preprocess_strain_local(x, Fs, detrendMethod)
x = x(:);
finiteMask = isfinite(x);
if any(finiteMask)
    x(~finiteMask) = median(x(finiteMask), 'omitnan');
else
    x(:) = 0;
end
switch lower(string(detrendMethod))
    case "none"
        y = x - median(x, 'omitnan');
    case "linear"
        y = detrend(x, 'linear');
    case "movingmedian"
        win = max(5, round(2.0 * Fs));
        y = x - movmedian(x, win, 'omitnan');
    otherwise
        error('Unknown detrendMethod: %s', detrendMethod);
end
y = y - median(y, 'omitnan');
end

function bladeId = resolve_compare_blade_local(IdentificationResult, requestedBladeId)
bladeId = requestedBladeId;
if isempty(IdentificationResult) || ~isfield(IdentificationResult, 'Trend') || isempty(IdentificationResult.Trend)
    return;
end
available = unique(IdentificationResult.Trend.BladeID(:).');
if isempty(available)
    return;
end
if ~ismember(requestedBladeId, available)
    bladeId = available(1);
end
end

function opr = load_opr_speed_profile_local(P)
opr = struct('time_sec', [], 'rot_hz', [], 'rot_rpm', [], 'source_file', "");
oprFile = fullfile(P.data.highSpeedPulseDir, 'jiluOPR.mat');
if exist(oprFile, 'file') ~= 2
    return;
end
loaded = load(oprFile, 'jiluOPR');
if ~isfield(loaded, 'jiluOPR') || isempty(loaded.jiluOPR)
    return;
end
oprTimes = loaded.jiluOPR(:, 1);
ppr = max(1, round(P.machine.oprPulsesPerRev));
if numel(oprTimes) <= ppr
    return;
end
dt = oprTimes((1 + ppr):end) - oprTimes(1:(end - ppr));
keep = isfinite(dt) & dt > 0;
if ~any(keep)
    return;
end
opr.time_sec = 0.5 * (oprTimes((1 + ppr):end) + oprTimes(1:(end - ppr)));
opr.time_sec = opr.time_sec(keep);
opr.rot_hz = 1 ./ dt(keep);
opr.rot_rpm = 60 * opr.rot_hz;
opr.source_file = string(oprFile);
end

function plot_opr_speed_and_identification_overview_local(opr, IdentificationResult, Tcmp, SSG, compareBladeId)
if ~SSG.showOprSpeedOverview
    return;
end

W = build_identification_window_overview_local(IdentificationResult, compareBladeId);
if isempty(opr.time_sec) && isempty(W)
    return;
end

fig = figure('Name', 'Step06SG OPR speed and identification overview', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 21, 13], 'NumberTitle', 'off');
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

xRange = resolve_overview_time_range_local(opr, W, SSG);

nexttile;
if ~isempty(opr.time_sec)
    plot(opr.time_sec, opr.rot_rpm, 'k-', 'LineWidth', 1.0, 'DisplayName', 'OPR-derived speed'); hold on;
else
    hold on;
end
if ~isempty(W)
    add_window_patches_local(gca, W);
    yRef = max(get(gca, 'YLim'));
    for k = 1:height(W)
        text(W.CenterTimeSec(k), yRef, sprintf('W%d', W.WindowID(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontName', 'Times New Roman', 'FontSize', 7.5, 'Color', [0.25 0.25 0.25]);
    end
end
xlim(xRange);
ylabel('Rotor speed (rpm)');
title(sprintf('B%d overview: OPR-derived speed and all sliding windows', compareBladeId), ...
    'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
style_axes_local();

nexttile;
hold on;
if ~isempty(W)
    yyaxis left;
    plot(W.CenterTimeSec, W.FrequencyHz, 'o-r', 'LineWidth', 1.2, ...
        'MarkerFaceColor', [0.9 0.15 0.15], 'DisplayName', 'Step06 frequency');
    if ~isempty(Tcmp) && height(Tcmp) > 0
        plot(Tcmp.WindowCenterTimeSec, Tcmp.SGPeakFreqNearStep06Hz, 's-b', ...
            'LineWidth', 1.0, 'MarkerFaceColor', [0.2 0.4 0.9], ...
            'DisplayName', 'strain local peak near Step06');
    end
    ylabel('Frequency (Hz)');
    yyaxis right;
    stem(W.CenterTimeSec, W.AmplitudeMM, 'Color', [0.15 0.55 0.25], ...
        'LineWidth', 1.0, 'Marker', 'none', 'DisplayName', 'Step06 amplitude');
    ylabel('Amplitude (mm)');
    for k = 1:height(W)
        yyaxis left;
        text(W.CenterTimeSec(k), W.FrequencyHz(k), sprintf(' EO%d', W.EO(k)), ...
            'Color', [0.2 0.2 0.2], 'FontName', 'Times New Roman', 'FontSize', 7.5, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end
xlim(xRange);
xlabel('BTT-aligned time (s)');
title('All sliding-window identification results in time order', 'FontWeight', 'normal');
grid on;
box on;
legend('Location', 'best', 'Box', 'off');
style_axes_local();
end

function W = build_identification_window_overview_local(IdentificationResult, compareBladeId)
W = table();
if isempty(IdentificationResult) || ~isfield(IdentificationResult, 'Trend') || isempty(IdentificationResult.Trend)
    return;
end
Trend = IdentificationResult.Trend;
Trend = Trend(Trend.BladeID == compareBladeId, :);
if isempty(Trend)
    return;
end

startSec = NaN(height(Trend), 1);
endSec = NaN(height(Trend), 1);
if isfield(IdentificationResult, 'WindowResult') && ~isempty(IdentificationResult.WindowResult)
    WR = IdentificationResult.WindowResult;
    WR = WR([WR.blade_id] == compareBladeId);
    for k = 1:height(Trend)
        idx = find([WR.window_id] == Trend.WindowID(k), 1, 'first');
        if isempty(idx) || ~isfield(WR(idx), 'Bundle') || isempty(WR(idx).Bundle)
            continue;
        end
        tNow = WR(idx).Bundle.T;
        startSec(k) = min(tNow, [], 'omitnan');
        endSec(k) = max(tNow, [], 'omitnan');
    end
end

W = table(Trend.WindowID, Trend.WindowCenterTimeSec, startSec, endSec, ...
    Trend.EO, Trend.FrequencyHz, Trend.AmplitudeMM, Trend.RotRPM, Trend.WeightedVoltageRMSE, ...
    'VariableNames', {'WindowID','CenterTimeSec','StartTimeSec','EndTimeSec', ...
    'EO','FrequencyHz','AmplitudeMM','RotRPM','WeightedVoltageRMSE'});
end

function xRange = resolve_overview_time_range_local(opr, W, SSG)
xAll = [];
if ~isempty(opr.time_sec)
    xAll = [xAll; opr.time_sec(:)];
end
if ~isempty(W)
    xAll = [xAll; W.CenterTimeSec(:)];
    xAll = [xAll; W.StartTimeSec(:)];
    xAll = [xAll; W.EndTimeSec(:)];
end
xAll = xAll(isfinite(xAll));
if isempty(xAll)
    xRange = [0 1];
    return;
end
xRange = [min(xAll) - SSG.overviewTimeMarginSec, max(xAll) + SSG.overviewTimeMarginSec];
end

function add_window_patches_local(ax, W)
if isempty(W)
    return;
end
yl = get(ax, 'YLim');
for k = 1:height(W)
    if ~isfinite(W.StartTimeSec(k)) || ~isfinite(W.EndTimeSec(k))
        continue;
    end
    patch(ax, ...
        [W.StartTimeSec(k) W.EndTimeSec(k) W.EndTimeSec(k) W.StartTimeSec(k)], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        [0.95 0.85 0.25], 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
uistack(findobj(ax, 'Type', 'line'), 'top');
end

function plot_strain_time_series_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId)
fig = figure('Name', 'Step06SG strain time series', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 20, 5 * numel(StrainGaugeData)], ...
    'NumberTitle', 'off');
tiledlayout(fig, numel(StrainGaugeData), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(StrainGaugeData)
    S = StrainGaugeData(i);
    nexttile;
    plot(S.time_btt_sec, S.raw_value, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8, ...
        'DisplayName', 'raw strain'); hold on;
    plot(S.time_btt_sec, S.value, 'b-', 'LineWidth', 0.9, ...
        'DisplayName', sprintf('detrended (%s)', S.detrend_method));
    if ~isempty(IdentificationResult)
        T = IdentificationResult.Trend;
        T = T(T.BladeID == compareBladeId, :);
        for k = 1:height(T)
            xline(T.WindowCenterTimeSec(k), '--', 'Color', [0.85 0.2 0.2], ...
                'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
    ylabel(sprintf('%s (%s)', S.label, S.unit_symbol), 'Interpreter', 'tex');
    title(sprintf('%s in BTT time', S.label), 'Interpreter', 'none', 'FontWeight', 'normal');
    if i == numel(StrainGaugeData)
        xlabel('BTT-aligned time (s)');
    end
    if i == 1
        legend('Location', 'best', 'Box', 'off');
    end
    style_axes_local();
end
end

function plot_strain_stft_with_step06_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId)
detailIdx = find(SSG.selectedChannels == SSG.detailChannel, 1, 'first');
if isempty(detailIdx)
    detailIdx = 1;
end
S = StrainGaugeData(detailIdx);
stft = strain_stft_local(S.time_btt_sec, S.value, SSG.strainSampleRateHz, ...
    SSG.stftWindowSec, SSG.stftStepSec, SSG.stftFreqBandHz);
ampDb = 20 * log10(max(stft.amp, eps));

fig = figure('Name', 'Step06SG STFT with Step06 overlay', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 20, 11], 'NumberTitle', 'off');
imagesc(stft.timeSec, stft.freqHz, ampDb);
axis xy;
colormap(turbo);
colorbar;
hold on;

if ~isempty(IdentificationResult)
    T = IdentificationResult.Trend;
    T = T(T.BladeID == compareBladeId, :);
    plot(T.WindowCenterTimeSec, T.FrequencyHz, 'wo-', 'LineWidth', 1.2, ...
        'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerSize', 6, ...
        'DisplayName', 'Step06 identified frequency');
    for k = 1:height(T)
        text(T.WindowCenterTimeSec(k), T.FrequencyHz(k), sprintf(' W%d', T.WindowID(k)), ...
            'Color', 'w', 'FontName', 'Times New Roman', 'FontSize', 8, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    legend('Location', 'best', 'Box', 'off');
end

xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
title(sprintf('%s detrended STFT with Step06 overlay', S.label), ...
    'Interpreter', 'none', 'FontWeight', 'normal');
style_axes_local();
end

function plot_strain_stft_3d_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId)
if ~SSG.showStrain3D
    return;
end

S = select_detail_strain_local(StrainGaugeData, SSG);
stft = strain_stft_local(S.time_btt_sec, S.value, SSG.strainSampleRateHz, ...
    SSG.stftWindowSec, SSG.stftStepSec, SSG.stftFreqBandHz);

fig = figure('Name', 'Step06SG strain 3D spectrogram', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 21, 13], 'NumberTitle', 'off');
surf(stft.timeSec, stft.freqHz, stft.amp, stft.amp, ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
view(46, 28);
colormap(turbo);
cb = colorbar;
cb.Label.String = sprintf('STFT amplitude (%s)', S.unit_symbol);
cb.Label.Interpreter = 'tex';
hold on;

if ~isempty(IdentificationResult)
    T = IdentificationResult.Trend;
    T = T(T.BladeID == compareBladeId, :);
    if ~isempty(T)
        zLine = max(stft.amp(:), [], 'omitnan');
        if ~isfinite(zLine)
            zLine = 0;
        end
        plot3(T.WindowCenterTimeSec, T.FrequencyHz, zLine * ones(height(T), 1), ...
            'wo-', 'LineWidth', 1.1, 'MarkerFaceColor', [0.9 0.1 0.1], ...
            'DisplayName', 'Step06 windows');
    end
end

xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
zlabel(sprintf('Amplitude (%s)', S.unit_symbol), 'Interpreter', 'tex');
title(sprintf('%s 3D spectrogram', S.label), 'Interpreter', 'none', 'FontWeight', 'normal');
grid on;
box on;
style_axes_local();
end

function plot_strain_local_spectrogram_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId)
S = select_detail_strain_local(StrainGaugeData, SSG);
stft = strain_stft_local(S.time_btt_sec, S.value, SSG.strainSampleRateHz, ...
    SSG.stftWindowSec, SSG.stftStepSec, SSG.stftFreqBandHz);
timeRange = sort(SSG.spec2DTimeRangeSec(:).');
timeMask = stft.timeSec >= timeRange(1) & stft.timeSec <= timeRange(2);
if ~any(timeMask)
    warning('No STFT column falls inside [%.3f, %.3f] s.', timeRange(1), timeRange(2));
    return;
end

fig = figure('Name', 'Step06SG selected-range spectrogram', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 17, 10.5], 'NumberTitle', 'off');
imagesc(stft.timeSec(timeMask), stft.freqHz, stft.amp(:, timeMask));
axis xy;
colormap(turbo);
cb = colorbar;
cb.Label.String = sprintf('STFT amplitude (%s)', S.unit_symbol);
cb.Label.Interpreter = 'tex';
hold on;

if ~isempty(IdentificationResult)
    T = IdentificationResult.Trend;
    T = T(T.BladeID == compareBladeId, :);
    keep = T.WindowCenterTimeSec >= timeRange(1) & T.WindowCenterTimeSec <= timeRange(2);
    T = T(keep, :);
    if ~isempty(T)
        plot(T.WindowCenterTimeSec, T.FrequencyHz, 'wo-', 'LineWidth', 1.2, ...
            'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerSize', 5.5, ...
            'DisplayName', 'Step06 identified frequency');
        for k = 1:height(T)
            text(T.WindowCenterTimeSec(k), T.FrequencyHz(k), sprintf(' W%d', T.WindowID(k)), ...
                'Color', 'w', 'FontName', 'Times New Roman', 'FontSize', 8, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end
        legend('Location', 'best', 'Box', 'off');
    end
end

xlabel('BTT-aligned time (s)');
ylabel('Frequency (Hz)');
title(sprintf('%s spectrogram in [%.3f, %.3f] s', S.label, timeRange(1), timeRange(2)), ...
    'Interpreter', 'none', 'FontWeight', 'normal');
style_axes_local();
end

function plot_strain_local_frequency_spectrum_local(StrainGaugeData, IdentificationResult, SSG, compareBladeId)
S = select_detail_strain_local(StrainGaugeData, SSG);
timeRange = sort(SSG.spec1DTimeRangeSec(:).');
timeMask = S.time_btt_sec >= timeRange(1) & S.time_btt_sec <= timeRange(2);
if nnz(timeMask) < 16
    warning('Too few strain samples inside [%.3f, %.3f] s for local spectrum.', ...
        timeRange(1), timeRange(2));
    return;
end

[freqHz, amp] = single_sided_spectrum_local(S.value(timeMask), SSG.strainSampleRateHz);
freqMask = freqHz >= SSG.stftFreqBandHz(1) & freqHz <= SSG.stftFreqBandHz(2);
freqHz = freqHz(freqMask);
amp = amp(freqMask);

fig = figure('Name', 'Step06SG selected-range frequency spectrum', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 17, 9], 'NumberTitle', 'off');
plot(freqHz, amp, 'k-', 'LineWidth', 1.2, 'DisplayName', 'strain local spectrum'); hold on;

if ~isempty(IdentificationResult)
    T = IdentificationResult.Trend;
    T = T(T.BladeID == compareBladeId, :);
    keep = T.WindowCenterTimeSec >= timeRange(1) & T.WindowCenterTimeSec <= timeRange(2);
    T = T(keep, :);
    if ~isempty(T)
        yTop = max(amp, [], 'omitnan');
        if ~isfinite(yTop)
            yTop = 0;
        end
        for k = 1:height(T)
            xline(T.FrequencyHz(k), '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
            text(T.FrequencyHz(k), 0.92 * yTop, sprintf('W%d EO%d', T.WindowID(k), T.EO(k)), ...
                'Rotation', 90, 'Color', [0.75 0.1 0.1], 'FontName', 'Times New Roman', 'FontSize', 7.5, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
    end
end

xlabel('Frequency (Hz)');
ylabel(sprintf('Amplitude (%s)', S.unit_symbol), 'Interpreter', 'tex');
title(sprintf('%s local spectrum in [%.3f, %.3f] s', S.label, timeRange(1), timeRange(2)), ...
    'Interpreter', 'none', 'FontWeight', 'normal');
box on;
grid on;
style_axes_local();
end

function plot_strain_frequency_compare_local(Tcmp, SSG, compareBladeId)
if isempty(Tcmp) || height(Tcmp) == 0
    return;
end

fig = figure('Name', 'Step06SG frequency comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1, 1, 18, 10], 'NumberTitle', 'off');
plot(Tcmp.WindowCenterTimeSec, Tcmp.Step06FrequencyHz, 'o-r', ...
    'LineWidth', 1.2, 'MarkerFaceColor', [0.9 0.1 0.1], ...
    'DisplayName', 'Step06 frequency'); hold on;
plot(Tcmp.WindowCenterTimeSec, Tcmp.SGPeakFreqNearStep06Hz, 's-b', ...
    'LineWidth', 1.1, 'MarkerFaceColor', [0.2 0.4 0.9], ...
    'DisplayName', sprintf('%s SG local peak', SSG.strainLabels{SSG.detailChannel}));
grid on;
box on;
xlabel('Window center time (s)');
ylabel('Frequency (Hz)');
title(sprintf('B%d: Step06 frequency vs strain local peak', compareBladeId), ...
    'FontWeight', 'normal');
legend('Location', 'best', 'Box', 'off');
style_axes_local();
end

function Tcmp = build_step06_strain_compare_table_local(IdentificationResult, StrainGaugeData, SSG, compareBladeId)
detailIdx = find(SSG.selectedChannels == SSG.detailChannel, 1, 'first');
if isempty(detailIdx)
    detailIdx = 1;
end
S = StrainGaugeData(detailIdx);
stft = strain_stft_local(S.time_btt_sec, S.value, SSG.strainSampleRateHz, ...
    SSG.stftWindowSec, SSG.stftStepSec, SSG.stftFreqBandHz);

Trend = IdentificationResult.Trend;
Trend = Trend(Trend.BladeID == compareBladeId, :);
if isempty(Trend)
    Tcmp = table();
    return;
end

rows = repmat(struct( ...
    'BladeID', NaN, ...
    'WindowID', NaN, ...
    'WindowCenterTimeSec', NaN, ...
    'Step06FrequencyHz', NaN, ...
    'Step06AmplitudeMM', NaN, ...
    'SGPeakFreqNearStep06Hz', NaN, ...
    'SGPeakAmp', NaN, ...
    'FreqDifferenceHz', NaN, ...
    'StrainLabel', ""), height(Trend), 1);

for k = 1:height(Trend)
    [~, tIdx] = min(abs(stft.timeSec - Trend.WindowCenterTimeSec(k)));
    freqMask = stft.freqHz >= Trend.FrequencyHz(k) - SSG.sgPeakSearchHalfBandHz & ...
               stft.freqHz <= Trend.FrequencyHz(k) + SSG.sgPeakSearchHalfBandHz;
    if any(freqMask)
        ampCol = stft.amp(freqMask, tIdx);
        freqCol = stft.freqHz(freqMask);
        [ampPeak, iPeak] = max(ampCol, [], 'omitnan');
        fPeak = freqCol(iPeak);
    else
        [ampPeak, iPeak] = max(stft.amp(:, tIdx), [], 'omitnan');
        fPeak = stft.freqHz(iPeak);
    end

    rows(k).BladeID = Trend.BladeID(k);
    rows(k).WindowID = Trend.WindowID(k);
    rows(k).WindowCenterTimeSec = Trend.WindowCenterTimeSec(k);
    rows(k).Step06FrequencyHz = Trend.FrequencyHz(k);
    rows(k).Step06AmplitudeMM = Trend.AmplitudeMM(k);
    rows(k).SGPeakFreqNearStep06Hz = fPeak;
    rows(k).SGPeakAmp = ampPeak;
    rows(k).FreqDifferenceHz = fPeak - Trend.FrequencyHz(k);
    rows(k).StrainLabel = string(S.label);
end

Tcmp = struct2table(rows);
end

function S = select_detail_strain_local(StrainGaugeData, SSG)
detailIdx = find(SSG.selectedChannels == SSG.detailChannel, 1, 'first');
if isempty(detailIdx)
    detailIdx = 1;
end
S = StrainGaugeData(detailIdx);
end

function stft = strain_stft_local(t, x, Fs, windowSec, stepSec, freqBandHz)
t = t(:);
x = x(:);
winN = max(8, round(windowSec * Fs));
stepN = max(1, round(stepSec * Fs));
nfft = 2 ^ nextpow2(winN);
freq = (0:(nfft/2)).' * Fs / nfft;
freqMask = freq >= freqBandHz(1) & freq <= freqBandHz(2);
freqKeep = freq(freqMask);

starts = 1:stepN:(numel(x) - winN + 1);
amp = NaN(numel(freqKeep), numel(starts));
timeSec = NaN(1, numel(starts));
w = hann_local(winN);
gain = mean(w);

for i = 1:numel(starts)
    idx = starts(i):(starts(i) + winN - 1);
    xw = x(idx);
    xw = xw - mean(xw, 'omitnan');
    y = fft(xw .* w, nfft);
    a2 = abs(y / (winN * gain));
    a1 = a2(1:(nfft/2 + 1));
    if numel(a1) > 2
        a1(2:end-1) = 2 * a1(2:end-1);
    end
    amp(:, i) = a1(freqMask);
    timeSec(i) = mean(t(idx([1, end])));
end

stft.timeSec = timeSec;
stft.freqHz = freqKeep;
stft.amp = amp;
end

function [freqHz, amp] = single_sided_spectrum_local(x, Fs)
x = x(:);
x = x - mean(x, 'omitnan');
n = numel(x);
nfft = 2 ^ nextpow2(max(n, 8));
w = hann_local(n);
gain = mean(w);
y = fft(x .* w, nfft);
a2 = abs(y / (n * gain));
a1 = a2(1:(nfft/2 + 1));
if numel(a1) > 2
    a1(2:end-1) = 2 * a1(2:end-1);
end
freqHz = (0:(nfft/2)).' * Fs / nfft;
amp = a1;
end

function w = hann_local(n)
if n <= 1
    w = ones(n, 1);
    return;
end
k = (0:n-1).';
w = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
end

function style_axes_local()
set(gca, 'FontName', 'Times New Roman', 'FontSize', 8.5, 'TickDir', 'in');
end

function tag = sensor_tag_local(sensorIds)
tag = ['S', sprintf('%d', sensorIds)];
end

function label = time_label_local(tSec)
if ~isfinite(tSec)
    label = 'TUnknown';
    return;
end
tMs = round(1000 * tSec);
if mod(tMs, 1000) == 0
    label = sprintf('T%03ds', tMs / 1000);
else
    label = sprintf('T%03dp%03ds', floor(tMs / 1000), mod(tMs, 1000));
end
end

function suffix = output_suffix_local(label)
suffix = '';
if isempty(label)
    return;
end
label = regexprep(char(label), '[^\w\d-]', '_');
label = regexprep(label, '_+', '_');
label = strtrim(label);
if ~isempty(label)
    suffix = ['_', label];
end
end

