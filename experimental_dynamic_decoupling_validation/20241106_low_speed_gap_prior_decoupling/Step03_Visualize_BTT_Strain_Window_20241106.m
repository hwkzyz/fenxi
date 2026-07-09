%% Step03: visualize BTT displacement and strain spectrum for 20241106
% This step helps choose a stable vibration window before any identification.
%
% Important:
% The default mode now uses the same physical time window for BTT and strain.
% Strain data are first shifted onto the BTT time axis using
%   BTT_time = strain_time + strainToBttOffsetSec
% and then the corresponding raw strain-file window is back-computed.
%
% No image or MAT file is saved.

clear; clc; close all;

%% 1. Simple settings
dataRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
lowDir = fullfile(dataRoot, '900');
highDir = fullfile(dataRoot, '3000_3150');

sampleRateHz = 5e6;
strainSampleRateHz = 10000;
rTipMm = 65.0;
bladeCount = 6;

% Start with the cleaner probes. Add 4 or 6 when debugging those channels.
sensorIds = [2 3 5 7];
targetBladeSlot = 1;

% Main physical window on the BTT axis.
bttWindowSec = [70 80];

% Clock alignment convention used later in Step04/Step05 as well:
%   BTT_time = strain_time + strainToBttOffsetSec
strainToBttOffsetSec = -102.6;

% Default: use the strain segment that corresponds to the same physical
% BTT window. Keep the manual mode for debug if needed.
strainWindowMode = 'aligned_to_btt';   % 'aligned_to_btt' | 'manual_strain_file_time'
manualStrainWindowSec = [150 210];

strainFiles = { ...
    fullfile(dataRoot, 'AI1-01_20241106170148.mat')
    fullfile(dataRoot, 'AI1-02_20241106170148.mat')
    fullfile(dataRoot, 'AI1-03_20241106170148.mat')
    fullfile(dataRoot, 'AI1-04_20241106170148.mat')
    };
strainLabels = {'AI1-01','AI1-02','AI1-03','AI1-04'};

fftBandHz = [0 800];
plotMaxPoints = 6000;
makeStrain3DFigure = true;
segmentDurationSec = 2.0;
segmentStepSec = 0.25;
strain3DTimeMode = 'all_aligned_time';   % 'all_aligned_time' | 'selected_btt_window'

strainWindowSec = resolve_strain_window_local( ...
    bttWindowSec, strainToBttOffsetSec, strainWindowMode, manualStrainWindowSec);

%% 2. Build BTT displacement using low-speed phase reference
low = load_case_timing_local(lowDir, sensorIds, sampleRateHz);
high = load_case_timing_local(highDir, sensorIds, sampleRateHz);
phaseRef = learn_low_speed_phase_reference_local(low, sensorIds, bladeCount);
highDisp = build_displacement_case_local(high, phaseRef, sensorIds, bladeCount, rTipMm);

%% 3. Load strain channels
strain = load_strain_set_local(strainFiles, strainLabels);

%% 4. Print diagnostics
fprintf('\n=== Step03 20241106 BTT/strain window visualization ===\n');
fprintf('BTT window: %.3f to %.3f s, target blade slot: %d, sensors: %s\n', ...
    bttWindowSec(1), bttWindowSec(2), targetBladeSlot, mat2str(sensorIds));
fprintf('Strain-to-BTT offset: %.3f s, using mode: %s\n', ...
    strainToBttOffsetSec, strainWindowMode);
fprintf('Corresponding strain-file window: %.3f to %.3f s\n', ...
    strainWindowSec(1), strainWindowSec(2));
fprintf('The strain plots/FFT are aligned to the same physical BTT window.\n\n');

print_btt_fft_summary_local(highDisp, sensorIds, targetBladeSlot, bttWindowSec, fftBandHz);
print_strain_fft_summary_local(strain, bttWindowSec, strainToBttOffsetSec, ...
    strainSampleRateHz, fftBandHz);

%% 5. Figures
plot_btt_overview_local(high, highDisp, sensorIds, targetBladeSlot, bttWindowSec, plotMaxPoints);
plot_strain_overview_local(strain, bttWindowSec, strainToBttOffsetSec, plotMaxPoints);
plot_fft_summary_local(highDisp, sensorIds, targetBladeSlot, bttWindowSec, ...
    strain, strainWindowSec, strainToBttOffsetSec, strainSampleRateHz, fftBandHz);
if makeStrain3DFigure
    for i = 1:numel(strain)
        plot_strain_3d_fft_local(strain(i), bttWindowSec, strainToBttOffsetSec, ...
            strainSampleRateHz, fftBandHz, segmentDurationSec, segmentStepSec, ...
            strain3DTimeMode);
    end
end

fprintf('\nDone. Figures are open in MATLAB. No image or MAT file was saved.\n');

%% Local functions: load BTT timing
function C = load_case_timing_local(caseDir, sensorIds, sampleRateHz)
C.caseDir = caseDir;
C.opr = load_opr_local(caseDir, sampleRateHz);
C.probe = load_probe_set_local(caseDir, sensorIds);
C.speedTime = C.opr.tCenter(1:end-1);
C.rpm = 60 ./ diff(C.opr.tCenter);
end

function opr = load_opr_local(caseDir, sampleRateHz)
S = load(fullfile(caseDir, 'jiluOPR.mat'), 'jiluOPR');
raw = S.jiluOPR;
opr.tStart = raw(:, 1) / sampleRateHz;
opr.tEnd = raw(:, 2) / sampleRateHz;
opr.tCenter = mean(raw(:, 1:2), 2) / sampleRateHz;
end

function probes = load_probe_set_local(caseDir, sensorIds)
emptyProbe = struct('sensorId', NaN, 'file', '', 'tCenter', []);
probes = repmat(emptyProbe, numel(sensorIds), 1);
for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    f = find_probe_file_local(caseDir, sid);
    if isempty(f)
        error('Missing probe file for CH%d in %s.', sid, caseDir);
    end
    S = load(f, 'jilublade');
    probes(i).sensorId = sid;
    probes(i).file = f;
    probes(i).tCenter = S.jilublade(:, 3);
end
end

function f = find_probe_file_local(caseDir, sensorId)
f1 = fullfile(caseDir, sprintf('jilublade_probe%d_nihe.mat', sensorId));
f2 = fullfile(caseDir, sprintf('jilublade_probe%d.mat', sensorId));
if exist(f1, 'file') == 2
    f = f1;
elseif exist(f2, 'file') == 2
    f = f2;
else
    f = '';
end
end

%% Local functions: BTT displacement
function phaseRef = learn_low_speed_phase_reference_local(C, sensorIds, bladeCount)
phaseRef = repmat(struct('sensorId', NaN, 'alpha', NaN(1, bladeCount)), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    phaseRef(is).sensorId = sensorIds(is);
    slotAlpha = collect_slot_phases_local(C.opr.tCenter, C.probe(is).tCenter, bladeCount);
    phaseRef(is).alpha = median(slotAlpha, 1, 'omitnan');
end
end

function Dcase = build_displacement_case_local(C, phaseRef, sensorIds, bladeCount, rTipMm)
Dcase = repmat(struct('sensorId', NaN, 'slot', []), numel(sensorIds), 1);
for is = 1:numel(sensorIds)
    Dcase(is).sensorId = sensorIds(is);
    Dcase(is).slot = repmat(struct('bladeSlot', NaN, 'time', [], 'xMm', []), 1, bladeCount);
    slotData = collect_slot_phase_times_local(C.opr.tCenter, C.probe(is).tCenter, bladeCount);
    for ib = 1:bladeCount
        alpha = slotData(ib).alpha;
        phaseError = wrap_to_half_turn_local(alpha - phaseRef(is).alpha(ib));
        Dcase(is).slot(ib).bladeSlot = ib;
        Dcase(is).slot(ib).time = slotData(ib).time;
        Dcase(is).slot(ib).xMm = rTipMm * 2 * pi * phaseError;
    end
end
end

function slotAlpha = collect_slot_phases_local(oprTime, probeTime, bladeCount)
slotData = collect_slot_phase_times_local(oprTime, probeTime, bladeCount);
minCount = min(arrayfun(@(s) numel(s.alpha), slotData));
slotAlpha = NaN(minCount, bladeCount);
for ib = 1:bladeCount
    slotAlpha(:, ib) = slotData(ib).alpha(1:minCount);
end
end

function slotData = collect_slot_phase_times_local(oprTime, probeTime, bladeCount)
slotData = repmat(struct('time', [], 'alpha', []), 1, bladeCount);
probeTime = probeTime(:);
for iRev = 1:(numel(oprTime) - 1)
    t0 = oprTime(iRev);
    t1 = oprTime(iRev + 1);
    tRev = sort(probeTime(probeTime > t0 & probeTime < t1));
    if numel(tRev) < bladeCount
        continue;
    end
    alphaRev = (tRev - t0) ./ (t1 - t0);
    [~, order] = sort(abs(alphaRev - 0.5), 'ascend');
    keep = sort(order(1:bladeCount));
    tRev = tRev(keep);
    alphaRev = alphaRev(keep);
    for ib = 1:bladeCount
        slotData(ib).time(end + 1, 1) = tRev(ib); %#ok<AGROW>
        slotData(ib).alpha(end + 1, 1) = alphaRev(ib); %#ok<AGROW>
    end
end
end

function y = wrap_to_half_turn_local(x)
y = mod(x + 0.5, 1.0) - 0.5;
end

%% Local functions: strain
function strain = load_strain_set_local(strainFiles, strainLabels)
strain = repmat(struct('label', '', 'file', '', 'time', [], 'value', []), numel(strainFiles), 1);
for i = 1:numel(strainFiles)
    if exist(strainFiles{i}, 'file') ~= 2
        error('Missing strain file: %s', strainFiles{i});
    end
    S = load(strainFiles{i}, 'Datas');
    strain(i).label = strainLabels{i};
    strain(i).file = strainFiles{i};
    strain(i).time = S.Datas(:, 1);
    strain(i).value = S.Datas(:, 2);
end
end

%% Local functions: diagnostics and plots
function print_btt_fft_summary_local(Dcase, sensorIds, bladeSlot, timeWindow, fftBandHz)
fprintf('BTT dominant frequencies in selected BTT window:\n');
for is = 1:numel(sensorIds)
    S = Dcase(is).slot(bladeSlot);
    [f, a] = local_nonuniform_fft_local(S.time, S.xMm, timeWindow);
    [fPeak, aPeak] = peak_in_band_local(f, a, fftBandHz);
    fprintf('  Probe %d: peak %.2f Hz, amplitude %.4g\n', sensorIds(is), fPeak, aPeak);
end
end

function print_strain_fft_summary_local(strain, bttWindow, offsetSec, fs, fftBandHz)
fprintf('\nStrain dominant frequencies in selected strain window:\n');
for i = 1:numel(strain)
    [f, a] = local_uniform_fft_local(strain(i).time, strain(i).value, bttWindow, offsetSec, fs);
    [fPeak, aPeak] = peak_in_band_local(f, a, fftBandHz);
    fprintf('  %s: peak %.2f Hz, amplitude %.4g\n', strain(i).label, fPeak, aPeak);
end
end

function plot_btt_overview_local(C, Dcase, sensorIds, bladeSlot, timeWindow, maxPts)
figure('Name', 'Step03 20241106 BTT overview', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(C.speedTime, C.rpm, 'k-', 'LineWidth', 1.0);
xline(timeWindow(1), 'r--'); xline(timeWindow(2), 'r--');
grid on; box on;
xlabel('BTT time (s)');
ylabel('Speed (rpm)');
title('High-speed OPR speed');

nexttile; hold on;
for is = 1:numel(sensorIds)
    S = Dcase(is).slot(bladeSlot);
    idx = decimate_index_local(numel(S.time), maxPts);
    plot(S.time(idx), S.xMm(idx), '.', 'MarkerSize', 4, ...
        'DisplayName', sprintf('Probe %d', sensorIds(is)));
end
xline(timeWindow(1), 'r--'); xline(timeWindow(2), 'r--');
grid on; box on;
xlabel('BTT time (s)');
ylabel('Relative displacement (mm)');
title(sprintf('BTT displacement, blade slot %d', bladeSlot));
legend('Location', 'best');
end

function plot_strain_overview_local(strain, bttWindow, offsetSec, maxPts)
figure('Name', 'Step03 20241106 strain overview', 'Color', 'w');
tiledlayout(numel(strain), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(strain)
    nexttile;
    idx = decimate_index_local(numel(strain(i).time), maxPts);
    tAligned = strain(i).time + offsetSec;
    plot(tAligned(idx), strain(i).value(idx), 'k-', 'LineWidth', 0.8);
    xline(bttWindow(1), 'r--'); xline(bttWindow(2), 'r--');
    grid on; box on;
    ylabel(strain(i).label);
    if i == 1
        title(sprintf('Strain channels aligned to BTT axis (offset %.3f s)', offsetSec));
    end
    if i == numel(strain)
        xlabel('Aligned BTT time (s)');
    end
end
end

function plot_fft_summary_local(Dcase, sensorIds, bladeSlot, bttWindow, strain, strainWindow, offsetSec, fs, fftBandHz)
figure('Name', 'Step03 20241106 FFT summary', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
for is = 1:numel(sensorIds)
    S = Dcase(is).slot(bladeSlot);
    [f, a] = local_nonuniform_fft_local(S.time, S.xMm, bttWindow);
    mask = f >= fftBandHz(1) & f <= fftBandHz(2);
    plot(f(mask), a(mask), 'LineWidth', 1.0, 'DisplayName', sprintf('Probe %d', sensorIds(is)));
end
grid on; box on;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(sprintf('BTT FFT, %.2f-%.2f s', bttWindow(1), bttWindow(2)));
legend('Location', 'best');

nexttile; hold on;
for i = 1:numel(strain)
    [f, a] = local_uniform_fft_local(strain(i).time, strain(i).value, bttWindow, offsetSec, fs);
    mask = f >= fftBandHz(1) & f <= fftBandHz(2);
    plot(f(mask), a(mask), 'LineWidth', 1.0, 'DisplayName', strain(i).label);
end
grid on; box on;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(sprintf('Strain FFT aligned to BTT %.2f-%.2f s (raw %.2f-%.2f s)', ...
    bttWindow(1), bttWindow(2), strainWindow(1), strainWindow(2)));
legend('Location', 'best');
end

function [f, amp] = local_nonuniform_fft_local(t, x, timeWindow)
mask = t >= timeWindow(1) & t <= timeWindow(2) & isfinite(x);
t = t(mask);
x = x(mask);
if numel(t) < 8
    f = NaN; amp = NaN; return;
end
dt = median(diff(t), 'omitnan');
fs = 1 / dt;
tq = (t(1):dt:t(end)).';
xq = interp1(t, x, tq, 'linear', 'extrap');
[f, amp] = simple_fft_local(xq, fs);
end

function [f, amp] = local_uniform_fft_local(t, x, bttWindow, offsetSec, fs)
tAligned = t + offsetSec;
mask = tAligned >= bttWindow(1) & tAligned <= bttWindow(2) & isfinite(x);
x = x(mask);
if numel(x) < 8
    f = NaN; amp = NaN; return;
end
[f, amp] = simple_fft_local(x, fs);
end

function strainWindowSec = resolve_strain_window_local(bttWindowSec, offsetSec, modeName, manualWindowSec)
switch lower(strtrim(modeName))
    case 'aligned_to_btt'
        strainWindowSec = bttWindowSec - offsetSec;
    case 'manual_strain_file_time'
        strainWindowSec = manualWindowSec;
    otherwise
        error('Unsupported strainWindowMode: %s', modeName);
end
if numel(strainWindowSec) ~= 2 || any(~isfinite(strainWindowSec))
    error('Resolved strain window must be a finite 1x2 vector.');
end
strainWindowSec = reshape(strainWindowSec, 1, 2);
if strainWindowSec(2) <= strainWindowSec(1)
    error('Resolved strain window must satisfy t2 > t1.');
end
end

function [f, amp] = simple_fft_local(x, fs)
x = x(:);
x = x - mean(x, 'omitnan');
n = numel(x);
y = fft(x);
amp2 = abs(y / n);
nHalf = floor(n / 2) + 1;
amp = amp2(1:nHalf);
if nHalf > 2
    amp(2:end-1) = 2 * amp(2:end-1);
end
if ~isempty(amp)
    amp(1) = 0;
end
f = (0:nHalf-1).' * fs / n;
end

function [fPeak, aPeak] = peak_in_band_local(f, a, band)
mask = isfinite(f) & isfinite(a) & f >= band(1) & f <= band(2);
if ~any(mask)
    fPeak = NaN; aPeak = NaN; return;
end
ff = f(mask);
aa = a(mask);
if ~isempty(aa)
    aa(1) = 0;
end
[aPeak, idx] = max(aa);
fPeak = ff(idx);
end

function idx = decimate_index_local(n, maxPts)
if n <= maxPts
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxPts))).';
end
end

function plot_strain_3d_fft_local(strainOne, bttWindowSec, offsetSec, fs, fftBandHz, segmentDurationSec, segmentStepSec, timeMode)
tAligned = strainOne.time(:) + offsetSec;
x = strainOne.value(:);
switch lower(strtrim(timeMode))
    case 'all_aligned_time'
        keep = isfinite(tAligned) & isfinite(x);
        titleSuffix = sprintf('all aligned time | highlight BTT %.2f-%.2f s', ...
            bttWindowSec(1), bttWindowSec(2));
    case 'selected_btt_window'
        keep = tAligned >= bttWindowSec(1) & tAligned <= bttWindowSec(2) & isfinite(x);
        titleSuffix = sprintf('selected BTT %.2f-%.2f s', bttWindowSec(1), bttWindowSec(2));
    otherwise
        error('Unsupported strain3DTimeMode: %s', timeMode);
end
if nnz(keep) < 8
    warning('Skipping 3D FFT for %s because the aligned window has too few samples.', strainOne.label);
    return;
end

tSel = tAligned(keep);
xSel = x(keep);
rawSelTime = strainOne.time(keep);

segLen = max(8, round(segmentDurationSec * fs));
segStep = max(1, round(segmentStepSec * fs));
n = numel(xSel);
startIdx = 1:segStep:(n - segLen + 1);
if isempty(startIdx)
    startIdx = 1;
    segLen = n;
end

freq = [];
ampMat = [];
segTime = zeros(numel(startIdx), 1);
rawTime = zeros(numel(startIdx), 1);
for k = 1:numel(startIdx)
    idx = startIdx(k):(startIdx(k) + segLen - 1);
    idx(idx > n) = [];
    xSeg = xSel(idx);
    [fSeg, aSeg] = simple_fft_local(xSeg, fs);
    if isempty(freq)
        freq = fSeg(:);
        ampMat = NaN(numel(freq), numel(startIdx));
    end
    ampMat(:, k) = aSeg(:);
    segTime(k) = mean(tSel(idx), 'omitnan');
    rawTime(k) = mean(rawSelTime(idx), 'omitnan');
end

freqMask = freq >= fftBandHz(1) & freq <= fftBandHz(2);
freqPlot = freq(freqMask);
ampPlot = ampMat(freqMask, :);
[tMesh, fMesh] = meshgrid(segTime, freqPlot);

figure('Name', sprintf('Step03 %s 3D segmented FFT', strainOne.label), ...
    'Color', 'w', 'NumberTitle', 'off');
surf(tMesh, fMesh, ampPlot, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(35, 30);
grid on; box on;
colormap(jet);
colorbar;
xline(bttWindowSec(1), 'w--', 'LineWidth', 1.0);
xline(bttWindowSec(2), 'w--', 'LineWidth', 1.0);
xlabel('Aligned BTT time (s)');
ylabel('Frequency (Hz)');
zlabel('FFT amplitude');
title(sprintf('%s 3D segmented FFT | %s | raw %.2f-%.2f s', ...
    strainOne.label, titleSuffix, ...
    min(rawTime, [], 'omitnan'), max(rawTime, [], 'omitnan')));
end
