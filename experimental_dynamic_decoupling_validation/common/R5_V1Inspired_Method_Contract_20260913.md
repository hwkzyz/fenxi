# V1-inspired R5 method contract

This document fixes the method boundary for the next experimental runs. It is
not a second forward model and it does not merge the V1 and R5 evaluators.

## What is retained from V1

The V1 `gap_only` route is retained as the dynamic identification contract:

- the same observation windows and sensor samples;
- the same rotor-angle phase law and EO candidate generation;
- the same continuous-frequency refinement around each EO;
- the same dynamic parameters `A`, `phase`, `dx`, and sensor-wise `dg`;
- the same voltage-only objective and final EO selection rule.

This is the part that already recovered the 20241106 and 20250527 frequency
branches. It must be replayed before any new static operator is evaluated.

## What R5 contributes

R5 contributes only a replacement for the calibrated static gap response
operator. For a fixed V1 dynamic state, its output is the sensor-wise
differential increment

```text
DeltaV_R5(x_V1, g0_V1, dg) = R5(x_V1, g0_V1 + dg) - R5(x_V1, g0_V1)
```

The V1 template, coordinates, gains, phase convention, and reference gap are
not re-applied by R5. In particular, no R5 template, nested Foundation anchor,
or second registration is part of the production method.

## Required staged verification

The implementation is accepted only in this order:

1. V1 self-replay reproduces the stored V1 waveform and frequency results.
2. With `dg = 0`, the R5 static operator returns exactly zero increment.
3. For fixed V1 dynamic parameters, the R5 increment is inspected sensor by
   sensor for sign, support, magnitude, and coordinate validity.
4. Only then is the V1 dynamic optimizer rerun with the R5 increment.

If stage 1 fails, the problem is in data adaptation. If stage 2 fails, the
problem is in the R5 differential operator. If stages 1--3 pass but frequency
changes, the candidate/phase contract was changed and the run is invalid.

## Prohibited shortcuts

- Do not lock EO to the V1 result; it remains a voltage-selected candidate.
- Do not use strain in the objective or in initialization.
- Do not add V1 and R5 gap increments together.
- Do not add a Foundation residual correction to hide a registration mismatch.
- Do not select the method because its amplitude is closer to strain.

Strain, amplitude transfer ratios, profile likelihood, and cross-condition
comparisons are posterior validation outputs only.

## Required run ledger

Every staged run must save the method contract, code revision, case, window,
EO, refined frequency, amplitude, `dx`, `dg`, RMSE, support/boundary flags, and
the corresponding V1 result. A change is not considered an improvement until
its first changed quantity and its physical interpretation are recorded.

## Verification record (2026-09-13)

The canonical voltage-only synthetic recovery was rerun after the contract was
written.  It recovered `EO=14`, `f=581.4000 Hz`, `A=0.2800 mm`,
`phase=0.3700 rad`, and `dx=0.0450 mm`, with `RMSE=1.27e-7 mV`.
This verifies numerical self-consistency of the current R5 backend only; it is
not evidence that the experimental static registration or amplitude scale is
correct.

## External input inventory (2026-09-13)

`R5_Build_ExternalInputManifest` indexed the existing external `gap_only`
folders without copying any data into the source tree. It found 94, 44, and
114 MAT artifacts for 20241106, 20250527, and 20251222 respectively. The
manifests classify response surfaces, templates, gap libraries, Foundation
results, V1 gap results, and sensor sidecars. These inventories are an input
selection aid only; no file is promoted to a formal R5 input until its role and
coordinate contract are reviewed.

For the 20250527 candidate set, top-level contracts were inspected before any
dynamic run: the response surface contains `responseSurface`, the low-speed
template contains `Template`, the V1 result contains `Result/Summary/Trend`,
and the sensor sidecar contains `SensorConditionedLibrary`. The gap-library
file is a calibration bundle, not a dynamic Foundation input. This distinction
prevents passing a calibration bundle as the window result and silently
changing the frequency-search contract.

## Current implementation status (2026-09-13)

### Direct-channel protection

The canonical evaluator now enforces the sensor-role boundary at runtime:
`M(i).isGapSensor == false` forces `dg_i=0`, so a direct V1 channel keeps its
original template response and cannot receive an R5 gap increment. This is a
structural invariant, not an EO or strain constraint. For 20241106, S2 is
therefore unchanged and only S5/S7 can use the R5 gap operator. MATLAB syntax
checking passed; the existing 20241106 replay remains the required regression
before accepting new three-condition results.

### Dimension-contract correction (2026-09-13)

The 20241106 sidecar stores several baseline and registration fields as
duplicated two-element vectors, while the dynamic operator requires one
effective scalar per sensor. The gap adapter now collapses these legacy fields
before response-surface evaluation. The canonical evaluator normalizes
observation vectors to columns and expands scalar increments per observation.
A one-window replay now reaches optimization without the previous oversized
array failure. After correcting the diagnostic field names, window 1 gives
`EO=12`, `f=631.696034 Hz`, `A=0.348004 mm`, and `RMSE=44.731417 mV`
(`status=pass`). This is only a one-window regression; the full 18-window
result remains pending.

The subsequent 18-window run completed without the dimension failure. Windows
1--17 selected EO12 with `f=631.48--632.41 Hz` and `A=0.318--0.348 mm`.
Window 18 selected EO6 (`f=316.787 Hz`, `A=0.179 mm`, `RMSE=44.070 mV`),
while its EO12 alternative had `RMSE=45.284 mV`. This is a near-tie between
an order-12 resonance and a subharmonic candidate, not a reason to lock EO.
The run is diagnostic-only for frequency stability; candidate discrimination
must be improved using the V1 voltage contract, without strain or reference-EO
constraints.

## Three-condition status table (2026-09-13)

| condition | V1 gap_only reference | current V1-shell/R5-gap diagnostic | acceptance |
|---|---|---|---|
| 20241106 B4 | EO12, about 632 Hz, 0.337--0.386 mm | windows 1--17 EO12, 631.48--632.41 Hz, 0.318--0.348 mm; window 18 EO6/EO12 near-tie | diagnostic until window-18 ambiguity is resolved |
| 20250527 B1 | EO14, about 580.2 Hz, 0.397--0.406 mm | previous R5-shell run EO19, about 787.6--787.8 Hz, 0.367--0.374 mm | rejected: frequency contract not reproduced |
| 20251222 B1 | V1 reference pending | one-window EO10, 416.089 Hz, 0.134 mm | rejected: one-window diagnostic only |

The table is deliberately a status record, not a choice of EO or a strain
calibration. Formal acceptance requires the current shell to reproduce the V1
frequency family for all three conditions while changing only the gap-sensor
operator.

The common worker now records `eo_near_tie`, `eo_near_tie_rmse_mv`, and
`eo_ambiguity` whenever refined EO candidates are within 3% RMSE of the
winner. This changes reporting only; it does not lock EO or alter the voltage
objective. It makes the 20241106 window-18 EO6/EO12 ambiguity explicit in every
future run.

The three external `gap_only` directories are present, but their artifacts
are not yet consumable by the common worker through one verified adapter. The
current `latest_programs` configurations resolve inputs relative to the
workspace case directory, while the authoritative historical files are under
the external `gap_only` trees. Consequently no new experimental R5 result is
reported here. The next admissible run is a V1 waveform self-replay followed
by the zero-`dg` differential check; a frequency or amplitude result produced
before those checks is diagnostic only.

The first external builder probe also found a concrete contract gap for
20241106: the V1 window declares analysis sensors `[2 5 7]`, whereas both
available R5 sensor sidecars contain only sensors `[5 7]`. The common worker
correctly rejects this as a sensor-role mismatch. Sensor 2 must therefore be
reconstructed or explicitly excluded by a reviewed case contract; silently
changing `analysisSensors` would invalidate the V1 frequency comparison.

After using the existing low-speed-bank adapter to restore sensor 2, a full
20241106 bridge run was attempted with the external V1 windows and R5 sidecar.
The worker completed, but the result is an invalid diagnostic: windows 1--5
selected EO16 at about 842 Hz, windows 6--9 and 12--18 selected EO17 at about
894 Hz, and amplitudes hit the 0.5 mm bound; RMSE was about 2.3e6 mV. Two
windows were invalid. Restoring the sensor list alone is therefore not
enough: the adapted template/response-surface voltage or coordinate contract
still has a major scale or registration mismatch. These values must not be
compared with V1 or used as experimental conclusions.

The fixed-state audit was run for the restored 20241106 model. With `dg=0`,
all three sensors had zero discrepancy and full support for S2/S7; S5 had
85.1% support because its response/template domains do not fully overlap. A
test increment `dg=[0,0.05,0.05]` produced negative increments of about
10.8--217.1 mV (S5) and 12.0--133.1 mV (S7), with no unit conversion inside
the audit. This confirms the zero-increment invariant, but also identifies
S5 domain overlap as a required gate before dynamic fitting.

## Fixed-state audit tool (2026-09-13)

`R5_AuditFixedStateIncrement` was added as an optimizer-free diagnostic. It
reports, per sensor, the template range, the maximum zero-`dg` discrepancy,
the finite-`dg` increment range, prediction range, and valid support fraction.
This separates unit/registration failures from EO and vibration optimization
and is now required before another experimental bridge run.

A one-window 20251222 B1 diagnostic selected EO10 at `416.089 Hz`, with
`A=0.134 mm` and `RMSE=775.875 mV`. This is not accepted as an experimental
result; it confirms that the unresolved dynamic EO-ranking problem also affects
20251222, so the three conditions must not yet share a production R5 dynamic
operator.

The 20250527 static observation audit passed for S1/S3/S6: zero-gap
equivalence, template replay, coordinate consistency, derivative finiteness,
and full 101-point support all passed. The remaining EO19 branch is therefore
not a basic template, unit, or static-domain failure; the next diagnostic must
compare each EO under the same fixed dynamic state.

The shared evaluator now applies the same unit conversion and normalizes
cell-valued MATLAB unit metadata, preventing equivalent `V`/`mV` contracts in
the 20250527 and 20251222 folders from being interpreted inconsistently.

The first post-correction 20250527 B1 run remains invalid: all 18 windows
selected EO19 with `f=787.609--787.775 Hz`, `A=0.367--0.374 mm`, and
`RMSE=55.832--62.403 mV` (the older sidecar variant similarly selected EO19
near 788 Hz and hit the amplitude bound). This is not a unit-scale failure;
it shows that the 20250527 R5 static operator changes the EO ranking relative
to V1. The result is retained as a diagnostic, and no EO or frequency value
from it is accepted as formal output.

The 20250527 builder was also made unit-aware at both its template preview and
dynamic baseline paths. Re-running all 18 windows produced the same EO19 and
787.6--787.8 Hz branch, confirming that its EO shift is not caused by a hidden
template `mV`/`V` conversion. The remaining discrepancy is therefore in the
20250527 static response-surface registration/scale or its migrated V1
candidate contract, and must be isolated with fixed-state per-EO replay.

## Unit correction (2026-09-13)

The 20241106 V1 template contract was checked directly: `Template.Sensor`
stores `v_grid` in `mV`, while the R5 builder previously multiplied every
template by `1000` as if it were in volts. This produced the observed
million-mV bridge residual. The builder now converts according to the explicit
`voltageUnit` field: `V` is converted once to `mV`, and `mV` is used directly.
S2 remains the unchanged V1 direct/template channel; the change is only the
unit boundary shared by the gap-sensor bridge.

The corrected 20241106 bridge was rerun on all 18 windows. It returned EO12
for every window, `f=631.406--632.379 Hz`, `A=0.335--0.384 mm`, and
`RMSE=35.096--44.155 mV`. This is consistent with the historical V1
`gap_only` frequency branch and its amplitude/RMSE scale. The prior EO16/17,
842/894-Hz, and million-mV result was therefore confirmed as a template-unit
failure, not an experimental EO jump. This run is the first valid 20241106
V1-contract/R5-gap-operator bridge baseline.

The same legacy-vector normalization was applied to the 20250527 gap adapter.
The complete rerun remains `EO19`, `787.609--787.775 Hz`,
`A=0.367--0.374 mm`, and `RMSE=55.832--62.403 mV`; the scalar-shape defect
was real but is not the cause of the 20250527 frequency-family mismatch.

### Dynamic-contract correction (2026-09-13)

The evaluator was rewritten to include the V1 spatial template-shift term for
gap sensors: `Vlow(x_current) + [VnoGap(x_current)-VnoGap(x_base)] +
DeltaV_gap`. Direct channels remain `Vlow(x_current)`. This is the intended
V1 dynamic contract with only the gap operator replaced. The 20241106 rerun
still gives EO12 for windows 1--17; window 18 remains an EO6/EO12 near-tie
(`f=317.060` versus `632.472 Hz`), so the correction is theoretically
necessary but does not by itself remove the subharmonic ambiguity.

The corresponding 20250527 rerun was completed after this correction. All 18
windows still select EO19, with `f=787.609--787.775 Hz`, `A=0.367--0.374 mm`,
and `RMSE=55.832--62.403 mV`. This branch remains diagnostic only. The next
required experiment is fixed-state replay of every V1 candidate EO (especially
EO14, EO11, and EO19) with the identical V1 window and support contract, so
that the ranking change can be attributed to the R5 gap increment itself.

The fixed-state complete-EO diagnostic for 20250527 window 1 confirms this
interpretation. EO14 gives `f=580.22 Hz`, `A=0.500 mm`, `RMSE=107.13 mV`;
EO11 gives `455.74 Hz`, `0.332 mm`, `100.68 mV`; EO19 gives `787.76 Hz`,
`0.197 mm`, `61.80 mV`; and EO10 gives `414.55 Hz`, `0.321 mm`, `74.35 mV`.
EO19 is therefore selected by the R5 voltage operator itself, not because the
correct V1 candidate was omitted. Forcing EO14 would be an invalid frequency
constraint. Before production results, the R5 increment must be audited at
fixed V1 states for sign, scale, x/g registration, and gain.

The increment bridge was then corrected to pass the actual target blade into
the R5 model builder (the old diagnostic accidentally passed the first sensor
ID). On 20250527, fixed V1 gap-only states give zero increment on S1/S6
because their V1 `dg` is zero; S3 is the active gap channel. For S3 the R5/V1
increment gain ratio is `0.709`, correlation `0.986`, RMS difference `22.7 mV`
(V1 RMS `77.3 mV`, R5 RMS `54.8 mV`). The waveforms therefore have nearly the
same spatial shape but R5 under-scales the active increment by about 29%.
This is the first direct numerical evidence explaining why the R5 operator
changes the EO ranking. The next correction must address the response-surface
scale/registration, not frequency or EO constraints.
The bridge audit also identified a concrete static-calibration mismatch. For
20250527 S3, V1 uses `tau=-0.181630 mm`, `xScale=0.989717`, and
`voltageGain=1.098428`, while the R5 sidecar uses `tau=-0.067074 mm`,
`xScale=1.005232`, and `voltageGain=0.877992`; S1/S6 show the same pattern.
Thus the current R5 increment is evaluated through a different coordinate and
gain contract, not merely a different response surface. V1 must own these
static quantities, with only the R5 surface coefficients replaced.

The first implementation of that interface was tested on all 18 windows of
20250527. The builder now accepts the V1 calibration object and overrides
`tau`, `xScale`, `muGapPerXMm`, and `voltageGain` before evaluating the R5
surface. The result remains EO19 (`787.649--787.803 Hz`), with
`A=0.191--0.198 mm` and `RMSE=58.814--65.819 mV`. This is useful evidence:
static calibration mismatch explains the amplitude scale change but does not
alone restore the V1 EO14 frequency family. The R5 surface coefficients or
their coordinate definition still require a fixed-state audit.
