# 20251222 OPRCenterStd Low-Speed Rotating Calibration

This folder now keeps only the current OPRCenterStd direct-template route in
the top level. Historical GradientXRange030, BaseFrame, scan, comparison, and
diagnostic programs were moved to `archive_previous_routes/`.

## Run

Run the complete current route:

```matlab
Run_Main_OPRCenterStd_20251222
```

Run step by step:

```matlab
Step01_Main_Build_OPRCenterStd_Template_20251222
Step02_Main_Build_OPRCenterStd_DynamicMap_20251222
Step03_Main_Run_OPRCenterStd_DirectTemplate_20251222
Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222
```

## Current File Roles

- `Run_Main_OPRCenterStd_20251222.m`: one-command current route.
- `Step01_Main_Build_OPRCenterStd_Template_20251222.m`: OPR pulse-center timing and low-speed OPRCenterStd template.
- `Step02_Main_Build_OPRCenterStd_DynamicMap_20251222.m`: dynamic sliding-window map in the same OPRCenterStd coordinate.
- `Step03_Main_Run_OPRCenterStd_DirectTemplate_20251222.m`: direct-template synchronous vibration identification.
- `Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20251222.m`: summary table and trend figure.

Each Step file is intentionally self-contained: tune parameters at the top,
then run the file directly. The implementation is embedded below as local
functions so the route is easier to read and adjust.

## Method

Relative to the older 20251222 GradientXRange030 route, the current route uses:

- multi-threshold OPR pulse-center timing;
- low-speed template and dynamic map in the same OPR-center standard-angle frame;
- direct template-based synchronous identification:

```text
V = T_low(x - dx_c - eta_s - A*sin(EO*theta + phi))
f = EO * rot_freq_mean
```

The default Step03 data selection keeps the 20251222 settings that were stable
in prior tests:

```matlab
analysis_start_time = 50.2
analysis_sensors = [1 2 3]
top_k_eo = 3
dynamic_effective_mode = gradient
pulse_mode = all
domain_selection_mode = soft
domain_soft_margin_mm = 0.20
sensor_eta_limit_mm = 0.03
sensor_eta_reg_weight_v_per_mm = 0
overshoot_penalty_weight = 100
```

## Outputs

- Template:
  `output/templates/Template_LowSpeedRotating_B1_S123_GradientXRange030_OPRCenterStd_20251222.mat`
- Dynamic map:
  `output/dynamic_maps/DynamicMap_B1_S123_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20251222.mat`
- Identification result:
  `output/identification/Result_Step03_Main_VPTop3SynchronousWaveform_B1_S123_Main_DirectTemplate_OPRCenterStd_20251222.mat`
- Summary:
  `output/identification/Step04_DirectTemplate_OPRCenterStd_Summary_20251222.csv`
- Trend figure:
  `output/figures/main_direct_template_oprcenterstd/Step04_DirectTemplate_OPRCenterStd_Trend_20251222.png`

## Archive

Historical programs are kept for traceability:

- `archive_previous_routes/gradient_xrange030_route/`: previous named main route and its reusable core scripts.
- `archive_previous_routes/diagnostics_and_comparisons/`: BaseFrame route, all-EO baseline, start-time scan, and soft-vs-hard visualization.
- `archive_previous_routes/old_programs/`: older experimental and diagnostic scripts.

These archived scripts are not needed for the current OPRCenterStd run.
