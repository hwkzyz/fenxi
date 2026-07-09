# Obsolete and diagnostic archive

These files are kept for traceability. They are not part of the current main
route.

The current runnable route lives in the parent folder. Do not run archived
scripts unless you are reproducing an old diagnostic.

## experimental_branches

These scripts were tried during method development and were replaced by the
current Step05I/Step06I/Step07J route:

```text
Step01B_Rebuild_OPR_MultiThreshold_Center_20250527.m
Step05G_Build_PerBlade_Response_Surface_20250527.m
Step05H_Build_OffsetCorrected_Shared_Response_Surface_20250527.m
Step05J_Learn_BladeGeometry_Modes_20250527.m
Step06F_Calibrate_TiltCorrected_GapLibrary_20250527.m
Step06G_Calibrate_PerBladeWeighted_GapLibrary_20250527.m
Step07C_LowTemplateTiltGap_VPFullWave_20250527.m
Step07D_PerBladeWeightedGap_VPFullWave_20250527.m
Step07I_OffsetTiltSharedGap_VPFullWave_20250527.m
```

## obsolete_routes

These early routes were replaced before the final Step07J route.

```text
Step06B_Fit_LowSpeed_Gap_Prior_20250527.m
```

Early low-speed gap-prior fit. Superseded by Step06I, which calibrates the
per-sensor tilted `(x,g)` path using the Step05I offset-tilt shared library.

```text
Step07_Run_Decoupled_Identification_20250527.m
```

Old hybrid low-speed/gap-library identification route. It is not the current
result path.

```text
Step07A_GapLibrary_StaticTemplate_FullWave_20250527.m
```

Pure gap-library full-wave diagnostic. Useful only to show why the raw gap
library alone is not enough for the experiment.

```text
Step07B_TiltCorrected_FullWave_Identification_20250527.m
```

Intermediate tilted-library full-wave route. It did not reuse the direct
low-speed-template Step03 bundle, so its comparison was not clean enough.

## diagnostics

These scripts are visualization or diagnosis helpers, not the current main
pipeline:

```text
Step05A_Visualize_Raw_Calibration_Centroids_20250527.m
Step05B_Visualize_Static_Continuous_Baseline_20250527.m
Step06A_Visualize_Dynamic_Baseline_20250527.m
Step06C_Diagnose_Static_Template_Correction_20250527.m
Step06D_Build_OPR_Anchored_GapLibrary_20250527.m
Step06E_Plot_LowTemplate_GapLibrary_Overlay_20250527.m
Step08_Compare_Decoupling_Models_20250527.m
Step09_Visualize_Decoupled_Results_20250527.m
```

The current Step07J and Step08J scripts export the comparison CSV files and
figures used by the main result, so the old Step08/Step09 scripts are no
longer needed for the main route.
