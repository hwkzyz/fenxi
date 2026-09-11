# 20251222 structured voltage funnel diagnostic

This directory is an isolated comparison. It does not modify `Main10`,
`Config_20251222.m`, `Check_Program_20251222.m`, or formal result files.

The diagnostic changes only the EO candidate funnel:

1. Baseline: the saved voltage-domain VP Top-3 candidates.
2. Funnel: every saved VP EO receives one short joint voltage-model update;
   the leading candidates receive a second short update; the best three then
   enter the same complete voltage-model optimizer used by the baseline.

Both routes use the same saved bundle, calibrated experimental forward model,
parameter bounds, regularization, continuous frequency refinement, and final
plain unweighted voltage RMSE. The low-speed template is intentionally held
fixed in this first comparison so that any difference is attributable to the
candidate funnel rather than to simultaneous template changes.

Run `Diagnostic_SingleSync_VoltageFunnel_20251222.m`. Set
`STEP07J_DIAG_WINDOW_IDS="1 7 13"` for selected windows, or
`STEP07J_DIAG_MAX_WINDOWS=1` for a first-window smoke test. The default runs
all saved windows. Subset runs use tagged output names and do not overwrite
the complete result.

Set `STEP07J_DIAG_RESULT_FILE` to another saved projection diagnostic result
to evaluate another region. The output name includes `R01`, `R07`, or `R04`
derived from the selected filename.

## 2026-08-11 result

The comparison was completed for 18 saved windows in each of R01, R07, and
R04 (54 windows total):

| Region | Top-3 set changed | Final EO changed | Mean funnel-baseline RMSE |
|---|---:|---:|---:|
| R01 | 0/18 | 0/18 | +0.0001435 mV |
| R07 | 4/18 | 0/18 | +0.00000018 mV |
| R04 | 7/18 | 0/18 | -0.00000577 mV |

No funnel fit hit the amplitude, dx, or gap bounds. The largest worsening was
0.00258 mV (R01 window 10), and the largest improvement was 0.000104 mV (R04
window 17). These changes are negligible relative to the 30--45 mV fitted
waveform RMSE.

Conclusion: the saved single-synchronous VP already recalls the winning EO.
The short-joint funnel can reorder secondary candidates but does not improve
the selected EO or final fit. It should not replace the formal single-sync
candidate route on this evidence. The simulation Funnel V2 benefit concerns
dual-frequency pair recall and should be tested only when an experimental
dual-frequency structure is justified.

Aggregated results are in:

- `SingleSync_VoltageFunnel_vs_VP_20251222_R01.csv`
- `SingleSync_VoltageFunnel_vs_VP_20251222_R07.csv`
- `SingleSync_VoltageFunnel_vs_VP_20251222_R04.csv`
- `SingleSync_VoltageFunnel_RegionSummary_20251222.csv`
