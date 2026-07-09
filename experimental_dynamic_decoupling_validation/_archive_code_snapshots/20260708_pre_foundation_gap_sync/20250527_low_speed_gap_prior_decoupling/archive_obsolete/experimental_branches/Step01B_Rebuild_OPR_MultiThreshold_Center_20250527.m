function rebuildInfo = Step01B_Rebuild_OPR_MultiThreshold_Center_20250527()
%STEP01B_REBUILD_OPR_MULTITHRESHOLD_CENTER_20250527 Refresh OPR timing centers.
% jiluOPR(:,1) is the robust OPR center used by downstream coordinate code.
% jiluOPR(:,2:3) keep the legacy threshold start/end times for diagnostics.

thisDir = fileparts(mfilename('fullpath'));
validationRoot = fileparts(thisDir);
mainLegacyDir = fullfile(validationRoot, '20250527', 'legacy');
localLegacyDir = fullfile(thisDir, 'legacy');
if exist(localLegacyDir, 'dir') ~= 7
    error('Local legacy directory not found: %s', localLegacyDir);
end
addpath(localLegacyDir);

cfg = Get_20250527_BTT_Config();
caseName = cfg.dynamic_cases{1};
caseDir = fullfile(cfg.dataset_root, caseName);
fileIds = list_case_file_ids_local(caseDir, cfg.opr_id);
if isempty(fileIds)
    error('No OPR files found under %s.', caseDir);
end

[centerTimes, startTimes, endTimes] = extract_opr_multithreshold_centers_local(caseDir, fileIds, cfg);
jiluOPR = [centerTimes(:), startTimes(:), endTimes(:)];
[omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(centerTimes, cfg.blades_num);
omega_zong = omega_rad_s;

metadata = struct();
metadata.method = 'multi_threshold_width_center';
metadata.description = ['OPR center is median over 30/40/50/60/70% pulse-height ' ...
    'rise/fall midpoint. Legacy threshold start/end are retained in jiluOPR columns 2 and 3.'];
metadata.threshold = cfg.opr_threshold;
metadata.level_ratios = [0.30 0.40 0.50 0.60 0.70];
metadata.pinlv = cfg.pinlv;
metadata.caseName = caseName;
metadata.createdBy = mfilename;

targetLegacyDirs = unique_existing_dirs_local({mainLegacyDir, localLegacyDir});
for i = 1:numel(targetLegacyDirs)
    outDir = fullfile(targetLegacyDirs{i}, 'output', caseName);
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    backup_file_if_needed_local(fullfile(outDir, 'jiluOPR.mat'), 'jiluOPR_legacy_start_edge_backup.mat');
    backup_file_if_needed_local(fullfile(outDir, 'omega.mat'), 'omega_legacy_start_edge_backup.mat');
    save(fullfile(outDir, 'jiluOPR.mat'), 'jiluOPR', 'metadata');
    save(fullfile(outDir, 'omega.mat'), 'omega_zong', 'omega_time_s', 'omega_rad_s', 'omega_rpm', 'metadata');
    fprintf('Saved OPR-center timing to:\n  %s\n', outDir);
end

fprintf('\nOPR multi-threshold center complete.\n');
fprintf('Pulses: %d\n', numel(centerTimes));
fprintf('Median center-start shift: %.3f us\n', median(1e6 * (centerTimes - startTimes), 'omitnan'));
fprintf('jiluOPR columns: [center_time_s, legacy_start_time_s, legacy_end_time_s]\n');

rebuildInfo = struct();
rebuildInfo.caseName = caseName;
rebuildInfo.pulseCount = numel(centerTimes);
rebuildInfo.medianCenterStartShiftUs = median(1e6 * (centerTimes - startTimes), 'omitnan');
rebuildInfo.method = metadata.method;
rebuildInfo.outputDirs = targetLegacyDirs;
end

%% Local functions
function [centerTimes, startTimes, endTimes] = extract_opr_multithreshold_centers_local(caseDir, fileIds, cfg)
tail = [];
gapPoints = 1e4;
centerTimes = [];
startTimes = [];
endTimes = [];
for iFile = 1:numel(fileIds)
    raw = load_raw_case_channel_local(caseDir, cfg.opr_id, fileIds(iFile));
    if isempty(raw)
        continue;
    end
    if iFile == 1 && cfg.initial_trim_points > 0 && size(raw, 1) > cfg.initial_trim_points
        raw(1:cfg.initial_trim_points, :) = [];
    end
    if ~isempty(tail)
        raw = [tail; raw]; %#ok<AGROW>
    end
    [segments, tail, gapPoints] = segment_signal_like_step2_local(raw, cfg.opr_threshold, gapPoints);
    for iSeg = 1:numel(segments.start_idx)
        centerTimes(end + 1, 1) = compute_opr_multithreshold_center_local( ...
            raw, segments.start_idx(iSeg), segments.end_idx(iSeg), cfg); %#ok<AGROW>
        startTimes(end + 1, 1) = segments.start_sample(iSeg) / cfg.pinlv; %#ok<AGROW>
        endTimes(end + 1, 1) = segments.end_sample(iSeg) / cfg.pinlv; %#ok<AGROW>
    end
end
end

function t_center = compute_opr_multithreshold_center_local(raw, a, b, cfg)
pad = max(8, round(0.25 * (b - a + 1)));
a0 = max(1, a - pad);
b0 = min(size(raw, 1), b + pad);
t = raw(a0:b0, 1) / cfg.pinlv;
v = raw(a0:b0, 2);
try
    vs = smooth(v, 16);
catch
    vs = smoothdata(v, 'movmean', 16);
end
[peakVal, iPeak] = max(vs);
baseVal = median(vs(vs <= prctile(vs, 30)), 'omitnan');
if ~isfinite(baseVal)
    baseVal = min(vs);
end
levels = baseVal + [0.30 0.40 0.50 0.60 0.70] .* max(peakVal - baseVal, eps);
centers = nan(numel(levels), 1);
for k = 1:numel(levels)
    tr = crossing_time_local(t(1:iPeak), vs(1:iPeak), levels(k), 'rising');
    tf = crossing_time_local(t(iPeak:end), vs(iPeak:end), levels(k), 'falling');
    if isfinite(tr) && isfinite(tf)
        centers(k) = 0.5 * (tr + tf);
    end
end
t_center = median(centers, 'omitnan');
if ~isfinite(t_center)
    t_center = 0.5 * (raw(a, 1) + raw(b, 1)) / cfg.pinlv;
end
end

function tc = crossing_time_local(t, v, level, direction)
tc = NaN;
t = t(:);
v = v(:);
if numel(t) < 2
    return;
end
if strcmpi(direction, 'rising')
    idx = find(v(1:end-1) < level & v(2:end) >= level, 1, 'first');
else
    idx = find(v(1:end-1) >= level & v(2:end) < level, 1, 'last');
end
if isempty(idx)
    return;
end
dv = v(idx + 1) - v(idx);
if abs(dv) < eps
    tc = t(idx);
else
    alpha = (level - v(idx)) / dv;
    tc = t(idx) + alpha * (t(idx + 1) - t(idx));
end
end

function [segments, tail, nextGap] = segment_signal_like_step2_local(raw, threshold, gapPoints)
segments = struct('start_idx', [], 'end_idx', [], 'start_sample', [], 'end_sample', []);
tail = [];
nextGap = gapPoints;
if isempty(raw)
    return;
end
sig = raw(:, 2);
try
    sigSmooth = smooth(sig, 16);
catch
    sigSmooth = smoothdata(sig, 'movmean', 16);
end
chase = find(sigSmooth > threshold);
if isempty(chase)
    return;
end
sampleOrder = raw(chase, 1);
segStarts = 1;
segEnds = [];
for ii = 1:(numel(sampleOrder) - 1)
    c = sampleOrder(ii + 1) - sampleOrder(ii);
    if c > nextGap
        segEnds(end + 1, 1) = ii; %#ok<AGROW>
        segStarts(end + 1, 1) = ii + 1; %#ok<AGROW>
        nextGap = 0.6 * c;
    end
end
segEnds(end + 1, 1) = numel(sampleOrder);
tailPoint = chase(segStarts(end)) - floor(nextGap / 2);
if tailPoint > 0 && tailPoint < size(raw, 1)
    tail = raw(tailPoint:end, :);
end
if numel(segStarts) < 2
    return;
end
compStarts = segStarts(1:end-1);
compEnds = segEnds(1:end-1);
segments.start_idx = chase(compStarts);
segments.end_idx = chase(compEnds);
segments.start_sample = raw(segments.start_idx, 1);
segments.end_sample = raw(segments.end_idx, 1);
end

function raw = load_raw_case_channel_local(caseDir, channelId, fileId)
file = fullfile(caseDir, sprintf('4-%d-%d.mat', channelId, fileId));
if ~isfile(file)
    raw = [];
    return;
end
loaded = load(file);
fn = fieldnames(loaded);
raw = loaded.(fn{1});
raw(raw(:, 1) == 0, :) = [];
end

function [omega_time_s, omega_rad_s, omega_rpm] = compute_speed_local(oprTimes, pulsesPerRev)
if numel(oprTimes) <= pulsesPerRev
    omega_time_s = [];
    omega_rad_s = [];
    omega_rpm = [];
    return;
end
revPeriod = oprTimes((pulsesPerRev + 1):end) - oprTimes(1:(end - pulsesPerRev));
omega_time_s = oprTimes(1:(end - pulsesPerRev));
omega_rad_s = 2 * pi ./ max(revPeriod, eps);
omega_rpm = 60 ./ max(revPeriod, eps);
end

function dirs = unique_existing_dirs_local(candidates)
dirs = {};
for i = 1:numel(candidates)
    if exist(candidates{i}, 'dir') == 7 && ~any(strcmp(dirs, candidates{i}))
        dirs{end + 1} = candidates{i}; %#ok<AGROW>
    end
end
end

function backup_file_if_needed_local(file, backupName)
if exist(file, 'file') ~= 2
    return;
end
backupFile = fullfile(fileparts(file), backupName);
if exist(backupFile, 'file') ~= 2
    copyfile(file, backupFile);
end
end

function fileIds = list_case_file_ids_local(caseDir, channelId)
d = dir(fullfile(caseDir, sprintf('4-%d-*.mat', channelId)));
fileIds = [];
for i = 1:numel(d)
    tok = regexp(d(i).name, ['4-' num2str(channelId) '-(\d+)\.mat'], 'tokens', 'once');
    if ~isempty(tok)
        fileIds(end + 1) = str2double(tok{1}); %#ok<AGROW>
    end
end
fileIds = sort(fileIds);
end
