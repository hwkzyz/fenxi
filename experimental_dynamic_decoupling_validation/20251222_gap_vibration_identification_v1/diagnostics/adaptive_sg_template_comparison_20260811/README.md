# 20251222 adaptive-SG low-speed template diagnostic

This directory is an isolated comparison. It does not modify formal Step04,
Main10, configuration, check scripts, or frozen result files.

The diagnostic changes only the low-speed template smoother:

1. Current route: fixed span 9, median followed by mean smoothing.
2. Candidate route: a common physical Savitzky-Golay window selected by
   complete-lap grouped cross-validation, followed by the same downstream
   PCHIP evaluation already used by the experimental voltage model.

Both routes use the same raw low-speed point cloud, target blade, sensor
coordinates, support/query-safe contract, static gap-response increment,
high-speed bundles, frozen formal voltage-domain VP Top-3 candidates and
initial values, full voltage optimizer, parameter bounds, regularization,
and final plain voltage RMSE. Freezing the saved formal Top-3 prevents a
simplified diagnostic VP implementation from changing the baseline. Linear
VP rankings are recomputed for both templates and saved only as diagnostics.

Run `Diagnostic_AdaptiveSG_Template_Comparison_20251222.m`. Optional
environment variables:

- `STEP07J_DIAG_RESULT_FILE`: another saved R01/R07/R04 result bundle.
- `STEP07J_DIAG_WINDOW_IDS="1 7 13"`: selected windows.
- `STEP07J_DIAG_MAX_WINDOWS=1`: first-window smoke test.
- `STEP07J_DIAG_MAX_POINTS=1200`: common per-window diagnostic point cap.

The template CSV reports grouped-lap voltage prediction error and derivative
repeatability across folds. Derivative errors are repeatability metrics, not
errors against an unavailable experimental ground truth.

## 2026-08-11 result

The comparison was completed for all 18 windows in R01, R07, and R04 (54
windows total). Grouped-lap CV selected the same 1.62 mm (81-bin) SG window
for blade 1 and blade 5.

| Region | EO changed | Diagnostic VP Top-3 set changed | Mean adaptive-current RMSE | Amplitude SD ratio |
|---|---:|---:|---:|---:|
| R01 | 0/18 | 0/18 | +0.02994 mV | 0.9994 |
| R07 | 0/18 | 0/18 | +0.05744 mV | 1.0013 |
| R04 | 0/18 | 0/18 | +0.05801 mV | 0.9974 |

At low speed, adaptive SG reduced held-out-lap voltage RMSE for every tested
sensor by 0.146--0.286 mV. Cross-fold derivative RMSE also fell from
9.38--15.97 to 2.35--10.78 mV/mm, with derivative correlation remaining
above 0.9998 for the adaptive templates.

At high speed, the final EO never changed and the diagnostic linear-VP Top-3
set never changed. The order of two secondary R07 candidates swapped in one
window. Adaptive SG improved final voltage RMSE in 6/54 windows, but increased
the mean final RMSE by 0.0485 mV overall. No adaptive fit hit a parameter
bound, every fit retained a query-safe fraction of 1.0, and amplitude
dispersion was effectively unchanged.

The amplitude change is negligible. Mean adaptive-minus-current amplitude
was +0.0000676 mm in R01, -0.000105 mm in R07, and -0.000149 mm in R04.
The largest absolute change was 0.000369 mm (R04 window 7), or 0.2922%.
Across all 54 windows, the mean signed change was -0.0000619 mm and the mean
absolute change was 0.000107 mm. Regional amplitude standard-deviation ratios
were 0.9994, 1.0013, and 0.9974, so adaptive SG neither introduces a
meaningful amplitude bias nor reduces amplitude dispersion.

## Split-derivative screening check

A second isolated check retained the current formal forward template and
used a separately smoothed derivative only in the linear EO scan. The
two-standard-error low-speed guard selected 2.42 mm (121 bins). It moved the
winning R01 EO from rank 2 to rank 1 in windows 5 and 10, while retaining the
winner in Top-3 for all R01 and R07 windows. R04 window 18 remained rank 6;
the split derivative did not repair that known linear-VP miss. Across all 54
windows it changed one Top-3 set, involving only secondary R07 candidates.

This is a small ranking improvement, not evidence for changing the formal
candidate route: formal Top-3 already retained the final winner, and the
split derivative adds no new winning-EO recall.

Conclusion: adaptive SG clearly improves low-speed repeatability, but this
version does not improve the downstream experimental voltage fit. It should
not replace the formal template builder on the current evidence. The next
isolated test should examine an independently regularized derivative for EO
screening while keeping the forward voltage template fixed; changing the
candidate funnel at the same time would confound that test.

Aggregated results are in
`AdaptiveSG_Template_RegionSummary_20251222.csv`.
