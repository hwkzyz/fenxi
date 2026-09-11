%% Step01: visualize OPR and blade timing for 20241106
% This step uses the already extracted jiluOPR/jilublade files.
% It only prints diagnostics and opens figures. It does not save files.

clear; clc; close all;

%% 1. Simple settings
dataRoot = 'D:\博士-国科\试验台数据\新试验\20241106_2';
lowDir = fullfile(dataRoot, '900');
highDir = fullfile(dataRoot, '3000_3150');

sensorIds = [2 3 4 5 6 7];
bladeCount = 6;
sampleRateHz = 5e6;

% Important:
% The 20241106 OPR sensor/channel is not assumed to be the same as the
% 20250527/20251222 routes. This Step01 intentionally reads the already
% extracted jiluOPR.mat files. If raw PCI extraction is rebuilt later, define
% a separate 20241106 OPR channel and threshold instead of reusing old cfg.opr_id.
oprSource = 'pre-extracted jiluOPR.mat';

cases = struct( ...
    'name', {'Low speed 900 rpm', 'High speed 3000-3150 rpm'}, ...
    'dir', {lowDir, highDir});

%% 2. Load data
for ic = 1:numel(cases)
    cases(ic).opr = load_opr_local(cases(ic).dir, sampleRateHz);
    cases(ic).probe = load_probe_set_local(cases(ic).dir, sensorIds);
end

%% 3. Print quick diagnostics
fprintf('\n=== Step01 20241106 OPR and blade timing diagnostics ===\n');
fprintf('Data root: %s\n', dataRoot);
fprintf('Sensors: %s, blade count: %d, sample rate: %.0f Hz\n', ...
    mat2str(sensorIds), bladeCount, sampleRateHz);
fprintf('OPR source: %s. Raw PCI OPR channel is not reused from older routes.\n', oprSource);

for ic = 1:numel(cases)
    print_case_summary_local(cases(ic), bladeCount);
end

%% 4. Figures
plot_opr_speed_local(cases);
plot_probe_time_coverage_local(cases, sensorIds);
plot_probe_interval_quality_local(cases, sensorIds);

fprintf('\nDone. Figures are open in MATLAB. No image or MAT file was saved.\n');

%% Local functions
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
opr.rpmTime = opr.tCenter(1:end-1);
opr.rpm = 60 ./ diff(opr.tCenter);
end

function probes = load_probe_set_local(caseDir, sensorIds)
emptyProbe = struct('sensorId', NaN, 'file', '', 'data', [], ...
    'tStart', [], 'tEnd', [], 'tCenter', []);
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
    data = S.jilublade;
    probes(i).sensorId = sid;
    probes(i).file = probeFile;
    probes(i).data = data;
    probes(i).tStart = data(:, 1);
    probes(i).tEnd = data(:, 2);
    probes(i).tCenter = data(:, 3);
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

function print_case_summary_local(caseData, bladeCount)
opr = caseData.opr;
rpmMed = median(opr.rpm, 'omitnan');
rpmMin = min(opr.rpm, [], 'omitnan');
rpmMax = max(opr.rpm, [], 'omitnan');

fprintf('\n--- %s ---\n', caseData.name);
fprintf('OPR count = %d, time = %.3f to %.3f s, rpm median/min/max = %.1f / %.1f / %.1f\n', ...
    numel(opr.tCenter), opr.tCenter(1), opr.tCenter(end), rpmMed, rpmMin, rpmMax);

for i = 1:numel(caseData.probe)
    P = caseData.probe(i);
    if isempty(P.tCenter)
        fprintf('Probe %d: missing\n', P.sensorId);
        continue;
    end
    expectedProbeCount = numel(opr.tCenter) * bladeCount;
    countRatio = numel(P.tCenter) / expectedProbeCount;
    fprintf('Probe %d: count = %d, time = %.3f to %.3f s, count/(OPR*blade) = %.3f\n', ...
        P.sensorId, numel(P.tCenter), P.tCenter(1), P.tCenter(end), countRatio);
end
end

function plot_opr_speed_local(cases)
figure('Name', 'Step01 20241106 OPR speed', 'Color', 'w');
tiledlayout(numel(cases), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for ic = 1:numel(cases)
    nexttile;
    plot(cases(ic).opr.rpmTime, cases(ic).opr.rpm, 'k-', 'LineWidth', 1.0);
    grid on; box on;
    xlabel('Time (s)');
    ylabel('Speed (rpm)');
    title(cases(ic).name);
end
end

function plot_probe_time_coverage_local(cases, sensorIds)
figure('Name', 'Step01 20241106 probe timing coverage', 'Color', 'w');
tiledlayout(numel(cases), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for ic = 1:numel(cases)
    nexttile; hold on;
    for i = 1:numel(sensorIds)
        P = cases(ic).probe(i);
        if isempty(P.tCenter)
            continue;
        end
        idx = decimate_index_local(numel(P.tCenter), 5000);
        y = P.sensorId * ones(numel(idx), 1);
        plot(P.tCenter(idx), y, '.', 'MarkerSize', 3);
    end
    grid on; box on;
    xlabel('Blade passing time (s)');
    ylabel('Probe ID');
    yticks(sensorIds);
    title([cases(ic).name, ': detected blade-passing points']);
end
end

function plot_probe_interval_quality_local(cases, sensorIds)
figure('Name', 'Step01 20241106 probe interval quality', 'Color', 'w');
tiledlayout(numel(cases), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for ic = 1:numel(cases)
    nexttile; hold on;
    for i = 1:numel(sensorIds)
        P = cases(ic).probe(i);
        if numel(P.tCenter) < 2
            continue;
        end
        dt = diff(P.tCenter);
        idx = decimate_index_local(numel(dt), 5000);
        plot(P.tCenter(idx), dt(idx) * 1000, '.', 'MarkerSize', 3, ...
            'DisplayName', sprintf('Probe %d', P.sensorId));
    end
    grid on; box on;
    xlabel('Time (s)');
    ylabel('Adjacent blade interval (ms)');
    title([cases(ic).name, ': adjacent blade-passing interval']);
    legend('Location', 'best');
end
end

function idx = decimate_index_local(n, maxPoints)
if n <= maxPoints
    idx = (1:n).';
else
    idx = unique(round(linspace(1, n, maxPoints))).';
end
end
