# Gap Library And High-Speed Gap Definition Sensitivity

This folder isolates the influence of the static gap library from vibration,
noise, and frequency-branch errors.

## Purpose

The new two-VP method depends on both:

- the static template `F_g(x)`;
- the derivative template `F'_g(x)`.

Therefore, library interpolation error and the definition of the high-speed
truth gap can affect both gap estimation and vibration identification. The
first analysis should not include vibration or noise, otherwise library error
and frequency-branch error become mixed.

## Main Note

- `RESULTS_CN.md`
  Chinese summary of the present closed-loop findings for the equivalent
  response models. It connects the static-library theory in
  `03_static_template_library_theory` with the unchanged high-speed two-VP
  method.

## Scripts

- `Step01_Analyze_Template_Library_Error.m`
  Leave-one-gap-out template reconstruction error. It compares the true
  holdout curve with template curves reconstructed from the remaining gap
  library.

- `Step02_Analyze_Highspeed_Gap_Definition.m`
  No-vibration and no-noise high-speed waveform test. It checks whether the
  static gap estimator can recover the high-speed truth gap under different
  library definitions.

- `Step03_Run_New_Method_Library_Sensitivity.m`
  Small vibration-identification test for the standalone new method under
  different library definitions.

- `Step04_Run_New_Method_Gap_Sweep.m`
  Gap-sweep vibration-identification test for the standalone new method. It
  repeats the library sensitivity test for every measured `g_true` and records
  both template and derivative RMSE, so boundary-gap failures can be separated
  from optimizer or frequency-branch failures.
  The default `runMode = "smoke"` runs the most sensitive `g_true = 0.2 mm`
  case at 10 dB for a fast check. Set `runMode = "quick"` for representative
  boundary and middle gaps, or `runMode = "full"` for all measured gaps and
  `[20, 10, 5]` dB.

- `Step05_Run_Equivalent_Response_Model_ClosedLoop.m`
  Closed-loop verification for the static-template response models proposed in
  `03_static_template_library_theory`. It keeps the current two-VP method
  unchanged and only replaces the continuous static-template reconstruction
  rule, so the effect of `1/g`, calibrated power-law, or exponential response
  coordinates can be tested directly on high-speed identification.
  The script accepts `runMode = "smoke"`, `"quick"`, or `"full"`.
  `smoke` tests the most sensitive `g_true = 0.2 mm` case at 10 dB.
  `quick` adds representative middle and large gaps.
  `full` sweeps all measured gaps at `[20, 10, 5]` dB.

## Outputs

Results are written to:

```text
GapVib_Ortho_Verify/results/02_gap_library_sensitivity
```

Key output files:

- `template_library_error.csv`
- `highspeed_gap_definition_error.csv`
- `new_method_library_sensitivity.csv`
- `new_method_gap_sweep_library_sensitivity.csv`
- `equivalent_response_closed_loop.csv`
- `equivalent_response_closed_loop_summary.csv`

If `Step05` is run with `runMode = "quick"` or `"full"`, the outputs are
written with a suffix to avoid overwriting the smoke results:

- `equivalent_response_closed_loop_quick.csv`
- `equivalent_response_closed_loop_quick_summary.csv`
- `equivalent_response_closed_loop_full.csv`
- `equivalent_response_closed_loop_full_summary.csv`

## Interpretation

Use this order:

1. If `template_library_error.csv` already shows large template or derivative
   error, the template library itself is a limiting factor.
2. If the no-vibration high-speed test still has gap bias under
   `leave_one_out`, the bias is caused by library interpolation or gap truth
   definition, not by vibration identification.
3. If `ideal_in_library` is much better than `leave_one_out`, then adding or
   improving gap templates is more important than modifying the optimizer.
4. If `Step05` shows that a calibrated response coordinate improves the
   high-speed gap and frequency errors while the optimizer is unchanged, the
   main limitation is the static-template reconstruction rule rather than the
   two-VP identification logic.
