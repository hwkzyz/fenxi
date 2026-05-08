# GapVib_Ortho_Verify

This project has been reorganized around the current two-VP method and the
static gap-template library analysis.

## Current Directories

- `01_current_two_vp_main`
  Main implementation and key-case validation for the latest two-VP +
  vibration-basis gap correction + final joint refinement method.

- `02_gap_library_sensitivity`
  Sensitivity analysis for finite static gap libraries, interpolation error,
  high-speed gap definition, and boundary small-gap cases.

- `03_static_template_library_theory`
  Theory notes for reconstructing the continuous static template family
  `F_g(x)` and derivative `F'_g(x)` from finite calibrated static curves.

- `04_inv_log_2_main`
  Standalone entry folder for the current recommended
  `inv_log_2 = [1, 1/g, log(g)]` method, including key-case runs, compact
  closed-loop comparison, and direct-display plotting scripts.

- `func`
  Shared MATLAB functions used by the current method and analysis scripts.

- `results`
  Output folders using the same `01/02/03/04` numbering as the analysis
  modules.

## Recommended Run Order

Validate the current method first:

```matlab
run('01_current_two_vp_main/Step01_Run_Current_TwoVP_Key_Cases.m')
```

Run the standalone `inv_log_2` workflow here:

```matlab
run('04_inv_log_2_main/Step01_Run_InvLog2_Key_Cases.m')
run('04_inv_log_2_main/Step02_Run_InvLog2_ClosedLoop_Quick.m')
run('04_inv_log_2_main/Plot01_Show_InvLog2_Template_Reconstruction.m')
```

Then analyze the static gap-library limitation:

```matlab
run('02_gap_library_sensitivity/Step01_Analyze_Template_Library_Error.m')
run('02_gap_library_sensitivity/Step02_Analyze_Highspeed_Gap_Definition.m')
run('02_gap_library_sensitivity/Step03_Run_New_Method_Library_Sensitivity.m')
```

For a fast boundary small-gap check:

```matlab
run('02_gap_library_sensitivity/Step04_Run_New_Method_Gap_Sweep.m')
```

`Step04` defaults to `runMode = "smoke"`. Set it to `"quick"` or `"full"` in
the script for broader sweeps.

## Archive

Old comparison and development outputs are kept under `results/archive` when
needed for history. The current top-level workflow should use only the
directories listed above.
