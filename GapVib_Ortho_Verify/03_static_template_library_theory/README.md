# Static Template Library Theory

This folder analyzes the theoretical role of the static gap-template library
in full-waveform gap-vibration decoupling.

## Purpose

The two-VP gap-vibration method estimates an effective gap by matching a
continuous static template family to the high-speed waveform. Therefore, its
accuracy depends on whether the finite static library can reconstruct the
static no-vibration waveform at the true high-speed gap.

The key question is:

```text
Given finite calibrated static curves F_{g_i}(x), how should we reconstruct
the continuous template F_g(x) and its derivative F'_g(x) for high-speed gap
estimation and vibration separation?
```

## Main Documents

- `THEORY_CN.md`
  Chinese theory note. It describes the background, difficulty, physical
  distance-response model, gap-coordinate derivation, and validation criteria.

- `THEORY_FIELD_RESPONSE_CN.md`
  Chinese note on the capacitor-sensor field/capacitance/voltage-response
  chain. It explains what is already physically justified and what is still an
  equivalent low-dimensional response model rather than a full field solution.

- `THEORY_COMPLETE_FIELD_CALIBRATED_TEMPLATE_CN.md`
  Complete Chinese theory note for the current best interpretation: field
  mechanism plus calibrated static-template expansion. It positions
  `[1, 1/g, log(g)]` as the main physically supported candidate and
  explains how the reconstructed no-vibration template enters the complete
  gap-vibration identification theory.

- `THEORY_FIELD_BASIS_CN.md`
  Chinese note explaining the field-inspired basis functions such as
  `[1, 1/g, log(g)]`, `[1, 1/g, log(g), 1/g^2]`, and inverse-polynomial
  variants.

- `LITERATURE_SENSOR_RESPONSE_CN.md`
  Chinese literature note on how capacitive, inductive, and eddy-current tip
  clearance sensors are usually modeled. It separates ideal field formulas,
  conditioning-circuit relations, calibration functions, and finite-element
  electromagnetic models.

- `RESULTS_CN.md`
  Chinese summary of the first-round leave-one-gap-out results and the current
  interpretation of the best gap-coordinate candidates, including the current
  closed-loop correction after `02_gap_library_sensitivity/Step05`.

- `THEORY_RESPONSE_GEOMETRY_CN.md`
  Chinese note on the more fundamental response-geometry model
  `F_g(x)=H(g+s(x))`. It explains why this route is theoretically stronger
  than fixed gap coordinates, and why it should first be treated as an
  independent analysis branch rather than an immediate replacement.

- `RESULTS_RESPONSE_GEOMETRY_CN.md`
  Chinese summary of the first-round response-geometry test results.

- `RESULTS_FIELD_CALIBRATED_TEMPLATE_CN.md`
  Chinese result note for the field-calibrated template model. It records the
  leave-one-gap-out, high-order inverse-term, and high-speed closed-loop
  comparisons without mixing program-specific results into the theory note.

## Scripts

- `Step01_Analyze_Gap_Coordinate_Models.m`
  Leave-one-gap-out validation for physically motivated gap-coordinate models.
  It compares `g`, `1/g`, calibrated `g^{-p}`, second-order power-law
  expansion, and exponential response coordinates. The evaluation includes
  both static template error `F_g(x)` and derivative error `F'_g(x)`.

- `Step02_Verify_Equivalent_Response_Assumptions.m`
  Verification layer for the equivalent-response assumption. It selects
  representative models from Step 01, compares them against the current `1/g`
  baseline on matched holdouts, and checks full-library linearity using `R^2`
  for both `F_g(x)` and `F'_g(x)`.

- `Step03_Analyze_Response_Geometry_Model.m`
  Independent test for the response-geometry model `F_g(x)=H(g+s(x))`. It
  fits a nonparametric monotone response and a smooth geometry-offset profile
  from the static template library, then compares this more fundamental model
  against the current `1/g` baseline and the best calibrated power-law model.

- `Step04_Analyze_Field_Inspired_Basis_Models.m`
  Leave-one-gap-out comparison for field-inspired basis functions. It tests
  whether physically interpretable expansions such as `[1, 1/g, log(g)]` can
  match or improve the current best calibrated power-law model.

## Relation To Existing Modules

- `01_current_two_vp_main`
  Uses the static template family during static gap initialization, vibration
  basis gap correction, and final joint refinement.

- `02_gap_library_sensitivity`
  Quantifies the effect of finite-library interpolation and boundary
  extrapolation on gap and vibration identification.

This folder provides the theoretical basis for improving the template
reconstruction step rather than adding empirical penalties to the optimizer.

## Run

Run:

```matlab
run('03_static_template_library_theory/Step01_Analyze_Gap_Coordinate_Models.m')
run('03_static_template_library_theory/Step02_Verify_Equivalent_Response_Assumptions.m')
```

Outputs are written to:

```text
GapVib_Ortho_Verify/results/03_static_template_library_theory
```

Key output files:

- `gap_coordinate_model_leave_one_out.csv`
  Per-holdout reconstruction errors.

- `gap_coordinate_model_summary.csv`
  Model-level ranking by combined template and derivative reconstruction error.

- `gap_coordinate_model_analysis.mat`
  MATLAB workspace data for later plotting and paper-table generation.

- `equivalent_response_validation_summary.csv`
  Summary of selected models against the current `1/g` baseline.

- `equivalent_response_baseline_comparison.csv`
  Per-holdout paired comparison between selected models and `1/g`.

- `equivalent_response_full_library_linearity.csv`
  Full-library `R^2` and fit RMSE for selected equivalent-response models.
- `02_gap_library_sensitivity/equivalent_response_closed_loop*.csv`
  Closed-loop evidence for whether the best static-library model also improves
  high-speed identification when the two-VP optimizer is unchanged.

## Workflow

1. Use `Step01` to identify whether the current `1/g` coordinate is justified
   by the static library.
2. Use `Step02` to verify paired improvements against the current `1/g`
   baseline for both `F_g(x)` and `F'_g(x)`.
3. Use `Step03` to test whether a more fundamental `H(d)+s(x)` model can
   beat the current low-dimensional calibrated response models under the
   existing template-library size.
4. Use `Step04` to compare physically interpretable field-inspired basis
   functions against the calibrated power-law models.
5. Feed representative candidates into
   `02_gap_library_sensitivity/Step05_Run_Equivalent_Response_Model_ClosedLoop.m`
   to decide which model still wins after entering the high-speed
   gap-vibration identification loop.
