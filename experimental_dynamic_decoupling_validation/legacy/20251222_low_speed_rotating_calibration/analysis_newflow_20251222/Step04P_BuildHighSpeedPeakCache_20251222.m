%% Step04P_BuildHighSpeedPeakCache_20251222
% Build a lightweight cache of existing high-speed blade timing rows.

clear; clc;
P = NewFlow_Config_20251222();
ensure_parent_dir_local(P.files.peakCache);

HighSpeedPeakCache = struct();
HighSpeedPeakCache.dataset = P.dataset;
HighSpeedPeakCache.mode = 'existing_jilublade_timing_cache';
HighSpeedPeakCache.case_output_dir = P.data.caseOutputDir;
HighSpeedPeakCache.analysis_sensors = P.sensors.analysis;
HighSpeedPeakCache.sensor = repmat(struct('sensor_id', NaN, 'jilublade_file', '', ...
    'row_count', NaN, 'time_min_sec', NaN, 'time_max_sec', NaN), numel(P.sensors.analysis), 1);

for i = 1:numel(P.sensors.analysis)
    sid = P.sensors.analysis(i);
    probeFile = fullfile(P.data.caseOutputDir, sprintf('jilublade_probe%d.mat', sid));
    require_file_local(probeFile, sprintf('jilublade timing for CH%d', sid));
    S = load(probeFile, 'jilublade');
    HighSpeedPeakCache.sensor(i).sensor_id = sid;
    HighSpeedPeakCache.sensor(i).jilublade_file = probeFile;
    HighSpeedPeakCache.sensor(i).row_count = size(S.jilublade, 1);
    HighSpeedPeakCache.sensor(i).time_min_sec = min(S.jilublade(:, 3), [], 'omitnan');
    HighSpeedPeakCache.sensor(i).time_max_sec = max(S.jilublade(:, 3), [], 'omitnan');
end

save(P.files.peakCache, 'HighSpeedPeakCache');

fprintf('\n=== Step04P: high-speed peak cache wrapper ===\n');
fprintf('Saved: %s\n', P.files.peakCache);

function require_file_local(pathText, label)
if exist(pathText, 'file') ~= 2
    error('Missing %s: %s', label, pathText);
end
end

function ensure_parent_dir_local(pathText)
folder = fileparts(pathText);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
