# Static Waveform Manifold Analysis

This independent analysis branch tests the static-library premise before it is
connected to low-speed or high-speed production identification.

Run `run_static_waveform_manifold_analysis.m` from this folder.

Pipeline (the regression gate used for every method revision):

1. `Step01_BuildStaticWaveformFamily.m`: fixes one B2 reference coordinate per
   nominal gap state and applies that same correction to every blade; independent
   non-reference peak alignment is diagnostic only.
2. `Step02_BuildBladeWaveformSurfaces.m`: fits one continuous gap surface per blade,
   `S_b(g,xi)`, with the established three gap basis functions.
3. `Step06_V1_SameBladeTransition_20251222.m`: V1 local adjacent-gap transition.
4. `Step07_V2_LeaveOneGap_20251222.m`: V2 same-blade leave-one-gap interpolation.
5. `Step14_StrictLeaveBladeOut_MethodComparison_20251222.m`: strict V3 method
   comparison. The held blade contributes only its low-state waveform; its other
   gaps and target waveform are hidden.
6. `Step09_V4_ExperimentalClosure_20251222.m`: V4 per-window comparison of the
   existing formal and no-intercept high-speed runs.
7. `Step10_ValidationGate_20251222.m`: writes the machine-readable pass/fail gate.

The scripts read the B2-anchored Step05 output and write only under this
folder's `results` directory. They do not change Main07, Main10, or production
gap libraries.

Current evidence: V1 PASS (90 adjacent transitions, P90 relative RMSE 8.959%)
and V2 PASS (96 held-out states, P90 7.405%). In strict blind V3, the best
deployable method is shifted nearest-surface localization plus a donor-surface
increment: all 90 transfers improve over zero increment, P90 relative RMSE is
8.559%, and P90 transfer/baseline is 0.528. It is therefore close to, but does
not pass, the original 0.50 gate. The oracle base-gap-conditioned result is
0.470 and shows that low-state localization, rather than the increment law, is
the remaining V3 bottleneck.

V4 is split into a method adapter and an invariant per-window comparison
framework. The existing free-`b_s` versus no-`b_s` pair is an engineering
closure regression, not an independent no-vibration truth test. The B2-donor
candidate preserves all 18 EO decisions relative to the formal shared model,
with mean absolute differences of 0.0194 Hz in frequency and 0.238 mV in RMSE.
Independent dynamic-gap or no-vibration truth is still unavailable.
