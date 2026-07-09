%% Step04_NumberHighSpeedSensorsInRegion_20251222
% Build current-region selected_rows_by_physical from existing jilublade IDs.

clear; clc;
P = NewFlow_Config_20251222();
ensure_parent_dir_local(P.files.highSpeedNumbering);
require_file_local(P.files.regionSelection, 'Step03 region selection');
require_file_local(P.files.peakCache, 'Step04P peak cache');

loadedRegion = load(P.files.regionSelection, 'RegionSelection');
RegionSelection = loadedRegion.RegionSelection;

HighSpeedNumbering = struct();
HighSpeedNumbering.dataset = P.dataset;
HighSpeedNumbering.mode = 'existing_jilublade_physical_blade_numbering';
HighSpeedNumbering.region_selection = RegionSelection;
HighSpeedNumbering.analysis_sensors = P.sensors.analysis;
HighSpeedNumbering.blade_count = P.machine.bladeCount;
HighSpeedNumbering.sensor = repmat(struct('sensor_id', NaN, ...
    'selected_rows_by_physical', [], 'selected_peak_times_by_physical', []), ...
    numel(P.sensors.analysis), 1);

for is = 1:numel(P.sensors.analysis)
    sid = P.sensors.analysis(is);
    probeFile = fullfile(P.data.caseOutputDir, sprintf('jilublade_probe%d.mat', sid));
    S = load(probeFile, 'jilublade');
    selectedRows = nan(P.region.targetBladePasses, P.machine.bladeCount);
    selectedTimes = nan(P.region.targetBladePasses, P.machine.bladeCount);
    for bladeId = 1:P.machine.bladeCount
        mask = S.jilublade(:, 4) == bladeId & ...
               S.jilublade(:, 3) >= RegionSelection.analysis_start_time_sec;
        rows = find(mask);
        if numel(rows) < P.region.targetBladePasses
            error('CH%d B%d has only %d passes after %.4f s; need %d.', ...
                sid, bladeId, numel(rows), RegionSelection.analysis_start_time_sec, ...
                P.region.targetBladePasses);
        end
        rows = rows(1:P.region.targetBladePasses);
        selectedRows(:, bladeId) = rows(:);
        selectedTimes(:, bladeId) = S.jilublade(rows, 3);
    end
    HighSpeedNumbering.sensor(is).sensor_id = sid;
    HighSpeedNumbering.sensor(is).selected_rows_by_physical = selectedRows;
    HighSpeedNumbering.sensor(is).selected_peak_times_by_physical = selectedTimes;
end

save(P.files.highSpeedNumbering, 'HighSpeedNumbering');

fprintf('\n=== Step04: high-speed numbering in selected region ===\n');
fprintf('Region: %s, target passes: %d\n', RegionSelection.tag, P.region.targetBladePasses);
fprintf('Saved:  %s\n', P.files.highSpeedNumbering);

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
