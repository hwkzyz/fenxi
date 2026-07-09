%% Step02: calculate simple BTT displacement for 20241106
% This script uses the 900 rpm data as the no-vibration timing reference.
% It does not reuse the old 20250527/20251222 OPR channel geometry.
%
% For each probe:
%   1. Convert each blade-passing time to a normalized OPR phase in one rev.
%   2. Use the low-speed 900 rpm run to learn six reference phases.
%   3. Convert high-speed phase deviations to circumferential displacement:
%
%        x_mm = r_tip_mm * 2*pi * wrap(alpha_high - alpha_ref)
%
% Figures are opened directly. No image or MAT file is saved.

clear; clc; close all;

%% 1. Simple settings
dataRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
lowDir = fullfile(dataRoot, '900');
highDir = fullfile(dataRoot, '3000_3150');

sensorIds = [2 3 4 5 6 7];
bladeCount = 6;
sampleRateHz = 5e6;
rTipMm = 65.0;
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% The blade slot is relative to this dataset's OPR pulse, not an old-date
% absolute blade label. Use 1 first; later steps can compare all slots.
targetBladeSlot = 1;

%% 2. Load low-speed and high-speed timing
low = load_case_timing_local(lowDir, sensorIds, sampleRateHz);
high = load_case_timing_local(highDir, sensorIds, sampleRateHz);

%% 3. Learn low-speed reference phases
phaseRef = learn_low_speed_phase_reference_local(low, sensorIds, bladeCount);

%% 4. Build high-speed displacement curves
highDisp = build_displacement_case_local(high, phaseRef, sensorIds, bladeCount, rTipMm);
lowCheck = build_displacement_case_local(low, phaseRef, sensorIds, bladeCount, rTipMm);

%% 5. Print diagnostics
fprintf('\n=== Step02 20241106 simple BTT displacement ===\n');
fprintf('Reference run: %s\n', lowDir);
fprintf('Dynamic run  : %s\n', highDir);
fprintf('Sensors: %s, blade slots: 1-%d, r_tip = %.2f mm\n', ...
    mat2str(sensorIds), bladeCount, rTipMm);
fprintf('Displacement sign follows x = R*2*pi*(phase - low-speed reference phase).\n');
fprintf('Blade slots are relative to the 20241106 OPR pulse.\n\n');

print_phase_reference_local(phaseRef, sensorIds, bladeCount);
print_displacement_summary_local('Low-speed self-check', lowCheck, sensorIds, bladeCount);
print_displacement_summary_local('High-speed result', highDisp, sensorIds, bladeCount);

%% 6. Figures
plot_phase_reference_local(phaseRef, sensorIds, bladeCount);
plot_speed_and_target_blade_local(high, highDisp, sensorIds, targetBladeSlot);
plot_sensor_slot_overview_local(highDisp, sensorIds, bladeCount);
plot_low_speed_self_check_local(lowCheck, sensorIds, targetBladeSlot);

%% 7. Save reusable BTT displacement result for later steps
matFile = fullfile(outDir, 'Step02_BTT_Displacement_20241106.mat');
save(matFile, 'highDisp', 'lowCheck', 'phaseRef', 'sensorIds', ...
    'bladeCount', 'sampleRateHz', 'rTipMm', 'lowDir', 'highDir');
fprintf('\nSaved reusable Step02 BTT result:\n  %s\n', matFile);

fprintf('\nDone. Figures are open in MATLAB. No image or MAT file was saved.\n');

%% Local functions
function C = load_case_timing_local(caseDir, sensorIds, sampleRateHz)
C.caseDir = caseDir;
C.opr = load_opr_local(caseDir, sampleRateHz);
C.probe = load_probe_set_local(caseDir, sensorIds);
C.speedTime = C.opr.tCenter(1:end-1);
C.rpm = 60 ./ diff(C.opr.tCenter);
end

function opr = load_opr_local(caseDir, sampleRateHz)
oprFile = fullfile(caseDir, 'jiluOPR.mat');
if exist(oprFile, 'file') ~= 2
    error('Missing OPR file: %s', oprFile);
end

S = load(oprFile, 'jiluOPR');
raw = S.jiluOPR;
opr.file = oprFile;
opr.raw = raw;
opr.tStart = raw(:, 1) / sampleRateHz;
opr.tEnd = raw(:, 2) / sampleRateHz;
opr.tCenter = mean(raw(:, 1:2), 2) / sampleRateHz;
end

function probes = load_probe_set_local(caseDir, sensorIds)
emptyProbe = struct('sensorId', NaN, 'file', '', 'data', [], 'tCenter', []);
probes = repmat(emptyProbe, numel(sensorIds), 1);

for i = 1:numel(sensorIds)
    sid = sensorIds(i);
    probeFile = find_probe_file_local(caseDir, sid);
    if isempty(probeFile)
        warning('Probe %d file not found in %s.', sid, caseDir);
        probes(i).sensorId = sid;
        continue;
    end

    S = load(probeFile, 'jilublade');
    probes(i).sensorId = sid;
    probes(i).file = probeFile;
    probes(i).data = S.jilublade;
    probes(i).tCenter = S.jilublade(:, 3);
end
end

function probeFile = find_probe_file_local(caseDir, sensorId)
preferred = fullfile(caseDir, sprintf('jilublade_probe%d_nihe.mat', sensorId));
fallback = fullfile(caseDir, sprintf('jilublade_probe%d.mat', sensorId));

if exist(preferred, 'file') == 2
    probeFile = preferred;
elseif exist(fallback, 'file') == 2
    probeFile = fallback;
else
    probeFile = '';
end
end

function phaseRef = learn_low_speed_phase_reference_local(C, sensorIds, bladeCount)
phaseRef = repmat(struct('sensorId', NaN, 'alpha', NaN(1, bladeCount), ...
    'phaseDeg', NaN(1, bladeCount), 'usableRevCount', 0), numel(sensorIds), 1);

for is = 1:numel(sensorIds)
    P = C.probe(is);
    phaseRef(is).sensorId = sensorIds(is);
    if isempty(P.tCenter)
        continue;
    end

    slotAlpha = collect_slot_phases_local(C.opr.tCenter, P.tCenter, bladeCount);
    phaseRef(is).usableRevCount = size(slotAlpha, 1);
    phaseRef(is).alpha = median(slotAlpha, 1, 'omitnan');
    phaseRef(is).phaseDeg = 360 * phaseRef(is).alpha;
end
end

function Dcase = build_displacement_case_local(C, phaseRef, sensorIds, bladeCount, rTipMm)
Dcase = repmat(struct('sensorId', NaN, 'slot', []), numel(sensorIds), 1);

for is = 1:numel(sensorIds)
    P = C.probe(is);
    Dcase(is).sensorId = sensorIds(is);
    Dcase(is).slot = repmat(struct('bladeSlot', NaN, 'time', [], 'xMm', [], ...
        'phaseError', []), 1, bladeCount);
    for ib = 1:bladeCount
        Dcase(is).slot(ib).bladeSlot = ib;
    end
    if isempty(P.tCenter)
        continue;
    end

    slotData = collect_slot_phase_times_local(C.opr.tCenter, P.tCenter, bladeCount);
    alphaRef = phaseRef(is).alpha;

    for ib = 1:bladeCount
        alpha = slotData(ib).alpha;
        t = slotData(ib).time;
        phaseError = wrap_to_half_turn_local(alpha - alphaRef(ib));
        Dcase(is).slot(ib).time = t;
        Dcase(is).slot(ib).phaseError = phaseError;
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
    inRev = probeTime > t0 & probeTime < t1;
    tRev = sort(probeTime(inRev));

    if numel(tRev) < bladeCount
        continue;
    end

    % Keep the six most central events if an interval has extra edge points.
    alphaRev = (tRev - t0) ./ (t1 - t0);
    [~, order] = sort(abs(alphaRev - 0.5), 'ascend');
    keep = sort(order(1:bladeCount));
    tRev = tRev(keep);
    alphaRev = alphaRev(keep);

    for ib = 1:bladeCount
        slotData(ib).time(end + 1, 1) = tRev(ib);
        slotData(ib).alpha(end + 1, 1) = alphaRev(ib);
    end
end
end

function y = wrap_to_half_turn_local(x)
y = mod(x + 0.5, 1.0) - 0.5;
end

function print_phase_reference_local(phaseRef, sensorIds, bladeCount)
fprintf('Low-speed reference phase, degrees inside one OPR period:\n');
header = sprintf('Probe  usableRev');
for ib = 1:bladeCount
    header = sprintf('%s   B%d_deg', header, ib);
end
fprintf('%s\n', header);
for is = 1:numel(sensorIds)
    row = sprintf('%5d  %9d', sensorIds(is), phaseRef(is).usableRevCount);
    for ib = 1:bladeCount
        row = sprintf('%s  %7.2f', row, phaseRef(is).phaseDeg(ib));
    end
    fprintf('%s\n', row);
end
end

function print_displacement_summary_local(name, Dcase, sensorIds, bladeCount)
fprintf('\n%s displacement std, mm:\n', name);
header = sprintf('Probe');
for ib = 1:bladeCount
    header = sprintf('%s      B%d', header, ib);
end
fprintf('%s\n', header);

for is = 1:numel(sensorIds)
    row = sprintf('%5d', sensorIds(is));
    for ib = 1:bladeCount
        x = Dcase(is).slot(ib).xMm;
        row = sprintf('%s  %7.4f', row, std(x, 'omitnan'));
    end
    fprintf('%s\n', row);
end
end

function plot_phase_reference_local(phaseRef, sensorIds, bladeCount)
figure('Name', 'Step02 20241106 low-speed reference phase', 'Color', 'w');
phaseMat = NaN(numel(sensorIds), bladeCount);
for is = 1:numel(sensorIds)
    phaseMat(is, :) = phaseRef(is).phaseDeg;
end
plot(1:bladeCount, phaseMat.', 'o-', 'LineWidth', 1.0);
grid on; box on;
xlabel('Blade slot relative to 20241106 OPR');
ylabel('Reference phase in one OPR period (deg)');
legend(compose('Probe %d', sensorIds), 'Location', 'best');
title('Low-speed 900 rpm timing reference');
end

function plot_speed_and_target_blade_local(C, Dcase, sensorIds, targetBladeSlot)
figure('Name', 'Step02 20241106 high-speed target blade displacement', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(C.speedTime, C.rpm, 'k-', 'LineWidth', 1.0);
grid on; box on;
xlabel('Time (s)');
ylabel('Speed (rpm)');
title('High-speed OPR speed');

nexttile; hold on;
for is = 1:numel(sensorIds)
    S = Dcase(is).slot(targetBladeSlot);
    plot(S.time, S.xMm, '.', 'MarkerSize', 4, ...
        'DisplayName', sprintf('Probe %d', sensorIds(is)));
end
grid on; box on;
xlabel('Time (s)');
ylabel('Relative displacement (mm)');
title(sprintf('High-speed displacement, blade slot %d', targetBladeSlot));
legend('Location', 'best');
end

function plot_sensor_slot_overview_local(Dcase, sensorIds, bladeCount)
figure('Name', 'Step02 20241106 high-speed all blade slots', 'Color', 'w');
tiledlayout(numel(sensorIds), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for is = 1:numel(sensorIds)
    nexttile; hold on;
    for ib = 1:bladeCount
        S = Dcase(is).slot(ib);
        plot(S.time, S.xMm, '.', 'MarkerSize', 2, ...
            'DisplayName', sprintf('B%d', ib));
    end
    grid on; box on;
    ylabel(sprintf('P%d mm', sensorIds(is)));
    if is == 1
        title('High-speed displacement by probe and blade slot');
    end
    if is == numel(sensorIds)
        xlabel('Time (s)');
    end
end
legend('Location', 'bestoutside');
end

function plot_low_speed_self_check_local(Dcase, sensorIds, targetBladeSlot)
figure('Name', 'Step02 20241106 low-speed self-check', 'Color', 'w');
hold on;
for is = 1:numel(sensorIds)
    S = Dcase(is).slot(targetBladeSlot);
    plot(S.time, S.xMm, '.', 'MarkerSize', 4, ...
        'DisplayName', sprintf('Probe %d', sensorIds(is)));
end
grid on; box on;
xlabel('Time (s)');
ylabel('Relative displacement (mm)');
title(sprintf('Low-speed self-check, blade slot %d should stay near 0', targetBladeSlot));
legend('Location', 'best');
end
