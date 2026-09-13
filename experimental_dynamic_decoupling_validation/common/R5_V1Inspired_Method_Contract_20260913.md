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
