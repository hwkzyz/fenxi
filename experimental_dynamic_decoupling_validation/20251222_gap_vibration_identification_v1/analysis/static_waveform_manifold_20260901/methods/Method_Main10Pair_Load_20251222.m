function [A,B] = Method_Main10Pair_Load_20251222(dataDir)
% Adapter for the free-b_s versus b_s=0 Main10 result pair.
% Both tables are vibration-identification outputs; B is not an independently
% measured no-vibration waveform truth.
freeFile=fullfile(dataDir,'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_unified_structured_voltage_r04_20260811.csv');
nobFile=fullfile(dataDir,'Main_GapAware_VPTopK_FullWave_Trend_20251222_B5_S123_nob_formal_B5_20260830.csv');
if ~isfile(freeFile)||~isfile(nobFile), A=table; B=table; return; end
A=readtable(freeFile); B=readtable(nobFile);
end
