%% Step03_SelectHighSpeedRegion_20251222
% Select a resonance-region-driven high-speed interval.

clear; clc;
P = NewFlow_Config_20251222();
ensure_parent_dir_local(P.files.regionSelection);

RegionSelection = struct();
RegionSelection.dataset = P.dataset;
RegionSelection.selection_source = 'ResonanceRegionCatalog_20251222';
RegionSelection.region_id = P.region.regionId;
RegionSelection.short_tag = P.region.shortTag;
RegionSelection.tag = P.region.tag;
RegionSelection.region_start_sec = P.region.regionStartSec;
RegionSelection.region_end_sec = P.region.regionEndSec;
RegionSelection.analysis_start_time_sec = P.region.analysisStartTimeSec;
RegionSelection.start_time_sec = P.region.analysisStartTimeSec;
RegionSelection.target_laps = P.region.targetBladePasses;
RegionSelection.window_laps = P.region.windowBladePasses;
RegionSelection.sliding_step_laps = P.region.slidingStepBladePasses;
RegionSelection.dominant_order = P.region.dominantOrder;
RegionSelection.dominant_freq_hz = P.region.dominantFreqHz;
RegionSelection.analysis_sensors = P.sensors.analysis;
RegionSelection.target_blade = P.case.bladeId;

save(P.files.regionSelection, 'RegionSelection');
writetable(struct2table(RegionSelection), strrep(P.files.regionSelection, '.mat', '.csv'));

fprintf('\n=== Step03: select high-speed resonance region ===\n');
fprintf('Region: %s [%g, %g] s, analysis start %.4f s\n', ...
    RegionSelection.tag, RegionSelection.region_start_sec, ...
    RegionSelection.region_end_sec, RegionSelection.analysis_start_time_sec);
fprintf('Saved:  %s\n', P.files.regionSelection);

function ensure_parent_dir_local(pathText)
folder = fileparts(pathText);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
