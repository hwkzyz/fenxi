# EO29 clearance-confounding attribution

Standalone, frozen-scope audit for the zero-vibration `g_low=0.5 mm -> g_high=0.8 mm` clearance-confounding question. It neither changes the production solver nor writes into the theory directory.

Run `Run_EO29_ClearanceConfoundingAttribution` from MATLAB. Outputs are written only to `output/` in this directory:

- `per_sensor_clearance_residual.csv`: static clearance residual and tangent-projection diagnostic.
- `deprecated_apparent_shift_scores.csv`: retained only to document the invalid apparent-shift diagnostic; it is not causal evidence.
- `operator_ablation_gamma.csv`: energy-consistent Full / no-intra-passage / no-angular-diversity / flattened-sensitivity Gamma table with rank and conditioning checks.
- `frozen_operator_span_audit.csv`: direct per-sensor audit showing why frozen-angle projection is null after global registration residualization.
- `fixed_reference_nonlinear_fits.csv`: fair EO23/29/30 nonlinear fits at 5/50/300 um with fixed `g_ref=0.5 mm`.
- `ANALYSIS_REPORT.md`: numerical conclusion and interpretation limits.
- `eo29_clearance_confounding_attribution.mat`: reproducibility payload.

The apparent shift is explicitly a local translation-equivalent projection, not a claim that clearance change caused true vibration. Its negligible explained energy invalidates it as a geometry-attribution target; the operator ablation retains the complete clearance residual instead.
