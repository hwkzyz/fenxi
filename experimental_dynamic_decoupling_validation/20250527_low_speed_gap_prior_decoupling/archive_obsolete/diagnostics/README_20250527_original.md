# 20250527 Experimental Validation Folder

This folder uses the verified processing scripts from the original
`20250527适配` workflow.

## Folder Layout

```text
20250527/
  legacy/   copied original preprocessing scripts
    output/ copied old Step2/Step3 products
  outputs/  new validation-stage summaries
  Step01_Extract_OPR_Blade_Timing_20250527.m
  Step02_Calc_FullTime_BTT_Displacement_20250527.m
  Step03_Locate_Vibration_By_Strain_OPR_20250527.m
  Step04_Detect_Resonance_By_BTT_STE_20250527.m
```

## Current Steps

Run these main scripts in order:

```matlab
Step01_Extract_OPR_Blade_Timing_20250527
Step02_Calc_FullTime_BTT_Displacement_20250527
Step03_Locate_Vibration_By_Strain_OPR_20250527
Step04_Detect_Resonance_By_BTT_STE_20250527
```

Each file is an independent main program. Avoid hiding new validation logic
inside extra subfunctions; later steps should continue as `Step04_...`,
`Step05_...`, etc.

These steps keep the old preprocessing logic intact:

1. Step01 checks and reuses copied old Step2 products: `jiluOPR.mat`,
   `omega.mat`, `jilublade_probe*.mat`;
2. Step02 calculates full-time BTT displacement once and saves a validation
   copy;
3. Step03 uses strain spectrum and OPR/RPM to locate the vibration window.

Old Step4 static calibration and Step5 super-Gaussian identification are
not copied into this folder. The later calibration/identification stage will
be rewritten for the proposed decoupling method.

Key validation copies are saved to:

```text
outputs/Step01_OPR_Blade_Timing_20250527.png
outputs/Step02_FullTime_BTT_Displacement_20250527.mat
outputs/Step02_FullTime_BTT_Displacement_20250527.png
outputs/Step03_Strain_OPR_Vibration_Window_20250527.mat
outputs/Step03_Strain_OPR_Vibration_Window_20250527.png
outputs/Step04_BTT_STE_Resonance_Regions_20250527.mat
outputs/Step04_BTT_STE_Resonance_Regions_20250527.csv
outputs/Step04_BTT_STE_Resonance_Regions_20250527.png
outputs/Step04_Aligned_Strain_STFT_Validation_20250527.png
```

The copied legacy products are under:

```text
legacy/output/20250526_2500-3500_t400/
```

By default, the Step scripts reuse these products and do not recompute them.
Set `forceRebuild = true` inside a Step script only when the old products
need to be regenerated.
