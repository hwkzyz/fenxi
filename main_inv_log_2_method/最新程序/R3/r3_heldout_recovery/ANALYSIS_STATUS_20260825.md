# R3 implementation and gate status (2026-08-25)

## Frozen scientific claim

Conditional on one frozen low-speed calibration, clearance-state adaptation
prevents mismatch at held-out operating clearance states from being attributed
to synchronous blade motion. The main endpoint is `P_vib`; `P_FP` and
`P_joint` remain distinct supporting endpoints.

## Completed gates

- Calibration states: `0.2, 0.3, 0.4, 0.6, 0.8, 1.1, 1.3, 1.4, 1.5 mm`.
- Held-out operating states: `0.5, 0.7, 0.9, 1.0, 1.2 mm`.
- Held-out static interpolation NRMSE: `0.80%, 1.18%, 1.54%, 1.63%, 0.79%`;
  maximum `1.63%`, below the predeclared `5%` gate.
- Noiseless calibrated state-matched gate: all 15 nonzero cases passed for
  `A = 0.10, 0.25, 0.37 mm`; every recovered order was `EO = 10`.
- Frozen noisy low-speed calibration: voltage noise standard deviation
  `0.000416872 V`; estimated reference coordinate `0.811123 mm` for the
  `0.8 mm` raw reference state.
- Sparse 15 dB smoke case at `g = 0.5 mm`, `A = 0.10 mm`, three paired seeds:
  fixed succeeded `0/3`, adaptive `3/3`, calibrated state-matched `3/3`.
  The fixed fit selected `EO = 23` with amplitude about `1.06 mm`; adaptive
  recovered `EO = 10`, amplitude about `0.1004 mm`, and gap about `0.505 mm`.

## Deliberately nonformal output

The detection-threshold workflow was exercised with only two null replicates.
It is marked `formal = false`; therefore `Run_05_Main15dB` correctly refuses to
start. Those provisional thresholds cannot be reported as `P_FP` evidence.

## Remaining formal computation

1. Run all 500 predeclared null realizations in
   `Run_04_CalibrateDetectionThreshold`.
2. Run `Run_05_Main15dB` over the six clearance states and seven amplitude
   levels (including the separate null level).
3. Generate Fig. 3 with `Plot_R3_Fig3` and supplementary null/joint evidence
   with `Plot_R3_Supplementary`.

No formal probability-domain claim is supported until those steps finish.
