# 20250527 Synthetic Spectrum and Gap-Bias Study

This folder isolates the current amplitude-bias question from the rest of the
experimental pipelines.

## Question

The 20250527 strain evidence shows a stable EO14 component near 580 Hz. The
waveform identification amplitudes differ across sensor combinations and gap
models, especially when CH6 is included. This study tests whether the bias is
caused by:

- static waveform error from CH6 being interpreted as vibration in the
  super-Gaussian route;
- residual tilt freedom being able to absorb EO14-like waveform changes in the
  gap_tilt route;
- or ordinary spectral/noise effects.

## Principle

Use the experimental Step5 bundle as the sampling skeleton. The synthetic study
does not invent ideal probe timing. It reuses:

- `bundle.X`
- `bundle.T`
- `bundle.S`
- `bundle.Theta`
- `bundle.W`
- `bundle.B_pt`, `bundle.w_pt`, `bundle.n_pt`, `bundle.xc_pt`, `bundle.base_pt`

Only the voltage observation is replaced by a controlled synthetic waveform.

Step08 adds the stricter physical version of the same idea. It reuses the
Step07J physical bundle, low-speed template, corrected gap library, and
response surface, then injects controlled vibration, gap change, initial tilt,
tilt-calibration error, and residual tilt-change cases. This keeps the main
question focused on gap change and calibrated initial tilt, not on assuming
large free high-speed dmu.

## Scripts

Run in order:

```matlab
Step00_Run_All_20250527
```

Or run separately:

```matlab
Step01_Build_SyntheticStrainEvidence_20250527
Step02_BundleVoltageBiasStudy_20250527
Step03_Jacobian_Coupling_Diagnostic_20250527
Step04_StaticDg_IdentifiabilitySweep_20250527
Step05_OrthogonalizedResidualTilt_20250527
Step06_GapTiltDirection_Decoupling_20250527
Step07_ProtectedTwoStageResidualTilt_20250527
Step08_PhysicalGapInitialTilt_SuperGaussianComparison_20250527
```

## Outputs

All outputs are written to:

```text
output/
```

Important files:

- `Step01_synthetic_strain_prior.mat`
- `Step01_synthetic_strain_evidence.png`
- `Step02_voltage_bias_detail.csv`
- `Step02_voltage_bias_summary.csv`
- `Step03_jacobian_coupling.csv`
- `Step04_static_dg_identifiability_sweep.csv`
- `Step05_orthogonalized_residual_tilt.csv`
- `Step06_gap_tilt_direction_decoupling.csv`
- `Step06_gap_tilt_direction_decoupling.png`
- `Step07_protected_two_stage_residual_tilt.csv`
- `Step08_physical_gap_initial_tilt_detail.csv`
- `Step08_physical_gap_initial_tilt_summary.csv`
- `Step08_physical_gap_initial_tilt_qa.csv`
- `Step08_physical_gap_initial_tilt_amplitude_error.png`
- `超高斯模型对比分析设计.md`

## Interpretation

The most useful diagnostics are:

- `amp_error_mm`: recovered amplitude minus true amplitude.
- `sensors`: compare `S13` against `S136`.
- `case_name`: compare clean, noise-only, CH6 static shift, and CH6 tilt shift.
- `rho_A_dmu`: direct correlation between amplitude and tilt columns.
- `proj_dmu_on_vib`: how much of a tilt column lies in the vibration subspace.
- `seed_policy` in Step04: whether the nonlinear gap_tilt fit was started
  independently, from the gap_only zero-dmu point, or from the synthetic truth.
- `dmu_column_meaning` in Step05: raw physical residual dmu versus an
  orthogonal residual column that is protected from stealing first-order
  vibration/static-gap content.
- `gap_equivalent_shift_per_unit_dmu_mm` in Step06: how much of the raw
  tilt column is actually a gap-offset component under the weighted voltage
  metric.
- `A_policy` in Step07: whether the vibration amplitude is locked to the
  main gap-fixed-tilt solution or allowed to move in a diagnostic linear
  refit.
- `Step08_physical_gap_initial_tilt_qa.csv`: correctness checks for the
  physical simulation. The clean no-gap and zero-amplitude cases must recover
  the truth, and the clean CH6 gap cases must recover the injected `dg6`.
- `Step08_physical_gap_initial_tilt_summary.csv`: paper-facing comparison of
  `superGaussian`, exact no-gap template, `gap_fixed_tilt`, and selected
  `gap_tilt_residual` diagnostics under no gap/no tilt, gap-only, gap plus
  calibrated initial tilt, tilt-calibration error, and tilt-change cases.

If CH6 static shift increases S136 amplitude in fixed/super-Gaussian-style
fits, while the tilt column has a large projection onto the vibration subspace,
then the experimental amplitude discrepancy has a structural mechanism rather
than being only a numerical accident.
