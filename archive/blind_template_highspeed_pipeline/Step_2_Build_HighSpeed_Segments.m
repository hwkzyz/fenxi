%% Step 2: cut local passing segments from the high-speed waveform
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'Data_High', 'truth');

[segments, xRef, activeThreshold] = extract_highspeed_segments(Data_High, truth, cfg);
segTable = make_segment_table(segments);

save(fullfile(outDir, 'stage2_segments.mat'), ...
    'segments', 'segTable', 'xRef', 'activeThreshold', '-v7.3');

figure('Name', 'Blind-Template Step 2 - Local Segments', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 10]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
nShow = min(12, numel(segments));
colors = lines(max(nShow, 3));
for i = 1:nShow
    plot(segments(i).xRef, segments(i).yInterp, '-', 'Color', colors(mod(i-1, size(colors, 1)) + 1, :), ...
        'DisplayName', sprintf('rev %d / sensor %d', segments(i).rev, segments(i).sensor));
end
xlabel('Local coordinate x (mm)');
ylabel('Voltage');
title('Example local passing segments');

nexttile;
bar(categorical(string(segTable.sensor)), segTable.count);
xlabel('Sensor id');
ylabel('Number of valid segments');
title('Valid passing segments per sensor');

fprintf('[Step 2] Valid segments extracted: %d\n', numel(segments));
disp(segTable);

function [segments, xRef, activeThreshold] = extract_highspeed_segments(Data_High, truth, cfg)
T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
OmegaNom = cfg.RPM_high * 2*pi / 60;
xRef = linspace(truth.domain(1), truth.domain(2), cfg.segment.xGridN)';
xHalf = 0.50 * diff(truth.domain);
base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + cfg.segment.activeLevel * span;
segments = struct('rev', {}, 'sensor', {}, 'tCenter', {}, 'xLocal', {}, ...
    'yLocal', {}, 'xRef', {}, 'yInterp', {}, 'mask', {}, 'weight', {}, 'trueShift', {});

for m = 1:numel(T_opr)-1
    Omega = 2*pi / (T_opr(m+1) - T_opr(m));
    if ~isfinite(Omega) || Omega <= 0
        Omega = OmegaNom;
    end
    for s = 1:numel(cfg.alpha_k)
        Te = T_opr(m) + cfg.alpha_k(s) / Omega;
        idx = find(abs(t - Te) <= xHalf / (Omega * cfg.R_tip));
        if isempty(idx), continue; end
        xLocal = Omega * cfg.R_tip * (t(idx) - Te);
        keep = xLocal >= truth.domain(1) & xLocal <= truth.domain(2) & V(idx) > activeThreshold;
        idx = idx(keep);
        xLocal = xLocal(keep);
        if numel(idx) < cfg.segment.minPointCount, continue; end
        [xUnique, ia] = unique(xLocal, 'stable');
        yUnique = V(idx);
        yUnique = yUnique(ia);
        yInterp = interp1(xUnique, yUnique, xRef, 'pchip', NaN);
        mask = isfinite(yInterp);
        if mean(mask) < cfg.segment.minValidFrac
            continue;
        end
        peakWeight = max(range(yInterp(mask)), eps);
        segWeight = max(cfg.segment.weightFloor, peakWeight * sqrt(nnz(mask) / numel(mask)));

        [~, tIdx] = min(abs(t - Te));
        trueShift = Data_High.u_t(tIdx);

        segments(end+1).rev = m; %#ok<SAGROW>
        segments(end).sensor = s;
        segments(end).tCenter = Te;
        segments(end).xLocal = xUnique(:);
        segments(end).yLocal = yUnique(:);
        segments(end).xRef = xRef;
        segments(end).yInterp = yInterp(:);
        segments(end).mask = mask(:);
        segments(end).weight = segWeight;
        segments(end).trueShift = trueShift;
    end
end
end

function segTable = make_segment_table(segments)
if isempty(segments)
    segTable = table([], [], 'VariableNames', {'sensor', 'count'});
    return;
end
sensors = [segments.sensor]';
uSensors = unique(sensors);
counts = zeros(size(uSensors));
for i = 1:numel(uSensors)
    counts(i) = nnz(sensors == uSensors(i));
end
segTable = table(uSensors, counts, 'VariableNames', {'sensor', 'count'});
end

function p = local_percentile(x, pct)
x = sort(x(isfinite(x)));
if isempty(x), p = NaN; return; end
q = 1 + (numel(x) - 1) * pct / 100;
lo = floor(q);
hi = ceil(q);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (q - lo) * (x(hi) - x(lo));
end
end
