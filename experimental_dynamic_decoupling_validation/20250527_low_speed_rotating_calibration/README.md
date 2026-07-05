# 20250527 Low-Speed Rotating Calibration

This folder keeps the latest OPRCenterStd direct-template workflow as separated
Step main programs. Each Step file can be opened, tuned, and run directly.

Run in order:

```matlab
Step01_Main_Build_OPRCenterStd_Template_20250527
Step02_Main_Build_OPRCenterStd_DynamicMap_20250527
Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20250527
Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20250527
```

## File Roles

- `Step01_Main_Build_OPRCenterStd_Template_20250527.m`: OPR-center timing and low-speed OPRCenterStd template.
- `Step02_Main_Build_OPRCenterStd_DynamicMap_20250527.m`: dynamic sliding-window waveform map.
- `Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20250527.m`: main direct-template vibration identification with `eta_s`.
- `Step03_Compare_Run_OPRCenterStd_DirectTemplate_NoEta_20250527.m`: no-`eta_s` comparison/baseline route.
- `Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20250527.m`: summary table and trend figure.

The tuning parameters are at the top of each Step file. The implementation used
by that Step is embedded below as local functions, so the Step files no longer
act as wrappers around other Step/core scripts.

## Current Route

The latest route uses:

- multi-threshold OPR pulse center timing;
- standard sensor angle converted into the OPR-center frame;
- non-parametric low-speed rotating templates;
- dynamic sliding windows built in the same OPRCenterStd coordinate;
- direct template-based synchronous identification without a static gap library.

## Key Outputs

- Template:
  `output/templates/Template_LowSpeedRotating_B1_S136_GradientXRange030_OPRCenterStd_20250527.mat`
- Dynamic map:
  `output/dynamic_maps/DynamicMap_B1_S136_SlidingWindows_Main20L_W3S1_GradientXRange030_OPRCenterStd_20250527.mat`
- Main identification result:
  `output/identification/Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_DirectTemplate_OPRCenterStd_Dx035_WithEta020_Reg002_PrevWinPhaseSafe_20250527.mat`
- No-eta comparison result:
  `output/identification/Result_Step03_Main_VPTop3SynchronousWaveform_B1_S136_Main_DirectTemplate_OPRCenterStd_Dx035_NoEta_PrevWinPhaseSafe_20250527.mat`
- Summary:
  `output/identification/Step04_DirectTemplate_OPRCenterStd_Summary_20250527.csv`
- Trend figure:
  `output/figures/main_direct_template_oprcenterstd/Step04_DirectTemplate_OPRCenterStd_Trend_20250527.png`

## Folder Layout

```text
.
|-- Step01_Main_Build_OPRCenterStd_Template_20250527.m
|-- Step02_Main_Build_OPRCenterStd_DynamicMap_20250527.m
|-- Step03_Main_Run_OPRCenterStd_DirectTemplate_WithEta_20250527.m
|-- Step03_Compare_Run_OPRCenterStd_DirectTemplate_NoEta_20250527.m
|-- Step04_Main_Visualize_OPRCenterStd_DirectTemplate_20250527.m
|-- archive_legacy_flows/
`-- output/
```

Historical GradientXRange, BaseFrame, VP/AllEO, old wrappers, previous
single-file attempts, and separated core files are kept only for traceability in
`archive_legacy_flows`.
