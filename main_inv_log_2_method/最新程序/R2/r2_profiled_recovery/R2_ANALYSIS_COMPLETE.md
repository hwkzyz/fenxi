# R2 profiled recovery mechanism: completed evidence audit

## Scope

R2 tests whether a fixed low-speed reference forces operating-state mismatch into motion parameters, and whether a calibrated state-dependent family restores the state-matched motion solution. The main mechanism uses a shared record-level gap, fixed correct EO, quadrature motion coefficients, and common interior support.

## Completed gates

1. Nested-objective identity: over 100 legal parameter draws, the maximum absolute difference between the fixed objective and the adaptive objective at the reference gap was `2.25e-10`; maximum relative difference was `9.01e-7`.
2. Canonical noiseless profile: for `g_low=0.50 mm`, `g_true=0.60 mm`, clearance-only fixed inference produced `A=0.08291 mm`, whereas adaptive/state-matched inference produced approximately zero motion. With `A_true=0.25 mm`, fixed inference produced `A=0.43165 mm`, while adaptive/state-matched inference produced `A=0.24999 mm`.
3. Off-grid interpolation check: `g_true=0.575 mm` was absent from the discrete state grid but shared the same static response library and interpolation rule as the inverse family. It is therefore an off-grid consistency check, not an independent-generator test.
4. Leave-one-state-out validation: the measured `0.60 mm` state was removed from the inverse response family, which was rebuilt from all remaining discrete states. Truth used the withheld `0.60 mm` response from the complete library. Clearance-only fixed inference produced `A=0.08359 mm`; leave-one-out adaptive inference selected `g=0.60 mm` and `A=0.00019 mm`. With `A_true=0.25 mm`, fixed inference produced `A=0.43340 mm`; leave-one-out adaptive inference selected `g=0.60 mm` and `A=0.25035 mm`. The residual at the interpolated state-matched solution was approximately `1.28e-4 V`, making the state-family interpolation error explicit.
5. Multi-start audit: 12 independent starts at five representative gaps converged to the same profiled solution. Maximum relative SSE spread was `1.81e-8` at the true gap and `1.63e-10` at the far endpoint; no material alternative basin was found.
6. Curvature audit: empirical profile curvature was `0.15793`; the nuisance-profiled Jacobian curvature was `0.15565`, giving a ratio of `1.0147`. Thus the local Schur-complement geometry explains the nonlinear profile curvature.
7. High-speed observation noise: at 15 dB, all 12 paired seeds selected `g_hat=0.575 mm`. Estimated amplitudes ranged from `0.24922` to `0.25053 mm`; the fixed-to-adaptive objective improvement remained positive for every seed.
8. Per-sensor extension smoke test: independent true gaps `[0.40,0.50,0.60] mm` were estimated as `[0.40417,0.51932,0.59643] mm`; the fixed-frequency pair was recovered as `[700.15,1199.73] Hz`. This is confirmatory SI evidence, not the main shared-gap profile.
9. Dual-sync production pilot: at 15 dB, the production solver recovered EO `[10,26]`, maximum amplitude error `0.000795 mm`, and clearance-increment error `2.99e-5 mm` for both tested search strategies.

## Supported claim

Under the declared mechanism conditions, fixed inference is the adaptive objective constrained to the reference state. Operating-state mismatch is therefore absorbed by motion parameters on the fixed slice. Allowing the calibrated state coordinate produces a distinct profiled minimum near the correct state and restores the state-matched motion solution. The result survives interpolation to a withheld discrete calibration state, solver multi-start audit, and 15 dB high-speed observation noise.

## Boundary

R2 does not establish a full operating-domain success probability, unknown-frequency model selection, weak-component limits, or combined low-speed/high-speed uncertainty. Those belong to R3/R4 or supplementary calibration studies. Passage-level free-gap over-flexibility is an optional SI ablation, not required for the R2 main claim.
