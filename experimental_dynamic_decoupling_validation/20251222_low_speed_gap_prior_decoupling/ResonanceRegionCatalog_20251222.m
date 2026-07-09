function [selection, catalog] = ResonanceRegionCatalog_20251222(regionId)
%RESONANCEREGIONCATALOG_20251222 Fixed resonance-region picks for 20251222.
%
% Region bounds come from Step04_BTT_STE_Resonance_Regions_20251222.csv
% after the raw segmented-FFT 100 microstrain threshold update. The
% analysisStartTimeSec values are fixed manual picks used to build
% comparable dynamic-map windows for subsequent identification.
%
% Optional runtime selector:
%   BLADE_RESONANCE_REGION_ID = 1, 2, ...

catalog = make_region(1, 44.7418, 59.4098, 50.2000, 20, 14, 582.22, ...
    'R01_EO14_main_50p2', 'Main EO14 region; keeps the historical 50.2 s start.');
catalog(end+1) = make_region(2, 19.8062, 21.2730, 20.0000, 20, 20, 580.86, ...
    'R02_EO20_20p0', 'Short EO20 region near 20.17 s.');
catalog(end+1) = make_region(3, 13.9390, 14.6724, 13.9400, 20, 20, 629.94, ...
    'R03_EO20_13p94', 'Short EO20 region near 14.31 s.');
catalog(end+1) = make_region(4, 31.5406, 33.0074, 31.7000, 20, 18, 629.94, ...
    'R04_EO18_31p7', 'EO18 region near 31.91 s.');
catalog(end+1) = make_region(5, 143.7508, 145.9510, 144.2000, 20, 20, 579.49, ...
    'R05_EO20_144p2', 'Low-speed EO20 region near 145.58 s.');
catalog(end+1) = make_region(6, 16.8726, 18.3394, 17.0000, 20, 20, 579.49, ...
    'R06_EO20_17p0', 'EO20 region near 17.24 s.');
catalog(end+1) = make_region(7, 38.8746, 40.3414, 39.1000, 20, 15, 579.49, ...
    'R07_EO15_39p1', 'EO15 region near 39.24 s.');

if nargin < 1 || isempty(regionId) || ~isfinite(regionId)
    raw = strtrim(getenv('BLADE_RESONANCE_REGION_ID'));
    if isempty(raw)
        regionId = 1;
    else
        regionId = str2double(raw);
    end
end

idx = find([catalog.regionId] == round(regionId), 1, 'first');
if isempty(idx)
    error('Unknown BLADE_RESONANCE_REGION_ID=%s. Available IDs are %s.', ...
        num2str(regionId), mat2str([catalog.regionId]));
end
selection = catalog(idx);
end

function r = make_region(regionId, regionStartSec, regionEndSec, analysisStartTimeSec, ...
    targetBladePasses, dominantOrder, dominantFreqHz, tag, note)
r = struct();
r.regionId = regionId;
r.tag = tag;
r.shortTag = sprintf('R%02d', regionId);
r.regionStartSec = regionStartSec;
r.regionEndSec = regionEndSec;
r.analysisStartTimeSec = analysisStartTimeSec;
r.targetBladePasses = targetBladePasses;
r.windowBladePasses = 3;
r.slidingStepBladePasses = 1;
r.dominantOrder = dominantOrder;
r.dominantFreqHz = dominantFreqHz;
r.note = note;
end
