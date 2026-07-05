# Improved 20250527 BTT Data Foundation Step01

This revision updates the uploaded Step01 low-speed reference layer while keeping
it method-neutral.

## Main changes

1. `BTTDataConfig_20250527.m`
   - Adds `cfg.opr_timing_method = 'multi_threshold_center'` by default.
   - Adds configurable probe timing via `cfg.probe_arrival_method = 'half_area'`.
   - Adds `cfg.time_index_mode = 'auto'` for local/absolute MAT time-index handling.
   - Makes the Blade-1 convention explicit: largest reference-sensor peak in the first low-speed search revolutions.
   - Adds fingerprint quality gates and angle-statistics settings.

2. `extract_low_speed_btt_features_20250527.m`
   - Supports OPR rising-edge or multi-threshold pulse-center timing.
   - Supports probe half-area, threshold-centroid, peak, and polynomial-centroid arrival time.
   - Records timing method and detected time-index mode in `case_data`.

3. `Step01_Build_LowSpeed_Reference_20250527.m`
   - Keeps the max-peak blade convention, but records it clearly in `Sensor_Config.Blade_ID_Definition`.
   - Adds best/second-best fingerprint correlation and correlation gap.
   - Saves angle mean/median/std/IQR/count/min/max.
   - Writes `Standard_Relative_Angles_Quality_20250527.csv`.

## Suggested run

Copy these three files into your MATLAB route folder, replacing the current Step01 data-foundation files, then run:

```matlab
Step01_Build_LowSpeed_Reference_20250527(false)
```

Use the generated `Sensor_Config_20250527.mat` as the common low-speed reference for later dynamic extraction and observation-bundle construction.
