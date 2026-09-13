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

## Fixed-state audit tool (2026-09-13)

`R5_AuditFixedStateIncrement` was added as an optimizer-free diagnostic. It
reports, per sensor, the template range, the maximum zero-`dg` discrepancy,
the finite-`dg` increment range, prediction range, and valid support fraction.
This separates unit/registration failures from EO and vibration optimization
and is now required before another experimental bridge run.
