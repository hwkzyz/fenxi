# R5 Sensor-Conditioned Decoupling Package

Package version: `R5_SensorConditioned_Decoupling_20260916`

This directory is the independent execution package for the two formal cases:

- `20250527`: sensors `[1 3 6]`, target blade B1
- `20251222`: sensors `[1 2 3]`, target blade B1 by default

The package does not add, call, or depend on the historical `latest_programs_*`
directories. All executable MATLAB files used by the packaged dynamic route are
under `core/` and `adapters/`. Required MAT inputs are copied under `inputs/`
and are checked before execution.

## Method

The formal route is:

```text
low-speed waveform
  -> low-speed template and response-surface state
  -> PC1/anchor/R4 absolute static localization
  -> frozen sensor-conditioned sidecar
  -> high-speed low_increment fitting
  -> complete nonlinear voltage replay
  -> Foundation replay audit
  -> independent strain/BTT posterior audit
```

The low-speed stage determines the frozen static observation operator. For each
sensor it supplies the template, coordinate registration, response-surface
state, reference gap `g0`, voltage gain, and declared support domains.

The high-speed stage fits only:

```text
EO, f, A, phi, dx, delta_g_s
```

The formal frequency search band remains `[300, 1000] Hz`. For each gap sensor:

```matlab
gLow  = g0_s;
gHigh = g0_s + delta_g_s;
Vpred = low_speed_template(x_current) + ...
        gain_s * (F_s(x_registered,gHigh) - F_s(x_registered,gLow));
```

The dynamic stage does not reconstruct a spatial gap slope using
`g0 + mu .* (x - tau)`. Any `mu` or historical tilt fields retained in a
sidecar are audit-only and do not enter the formal dynamic forward model.

Strain and BTT are independent posterior evidence only. The strain transfer
ratio has not been independently validated; amplitude is therefore provisional
and is not a hard acceptance criterion.

## Directory Layout

```text
R5_SensorConditioned_Decoupling_20260916/
  adapters/       case configurations, builders, entry adapters
  core/           common runner, canonical forward, audits, MAT staging
  inputs/         packaged calibration, templates, and Foundation MAT files
  results/        outputs, separated by case
  tests/          package-local regression tests
  tools/          dependency and input tools
  legacy/         package-specific legacy notes only
  run_R5_20250527.m
  run_R5_20251222.m
  run_R5_AllCases.m
  verify_R5_Installation.m
```

`inputs/` contains frozen calibration/template artifacts and prepared Foundation
window files. Raw experimental data are not required for the packaged dynamic
rerun. Rebuilding the low-speed templates or sidecar from raw waveforms is a
separate preparation workflow and must use an explicitly supplied input path.

## Running

From any current MATLAB working directory:

```matlab
restoredefaultpath;
rehash toolboxcache;
addpath('E:/.../experimental_dynamic_decoupling_validation/R5_SensorConditioned_Decoupling_20260916');
verify_R5_Installation('case','all');
run_R5_20250527;
run_R5_20251222;
```

The package entry points derive all paths from their own file location. They do
not use `pwd` and do not require the legacy directories on the MATLAB path.

To run both cases:

```matlab
restoredefaultpath;
addpath('E:/.../R5_SensorConditioned_Decoupling_20260916');
run_R5_AllCases;
```

A full run reads the frozen sidecar, low-speed template, and prepared Foundation
window bundle, performs the static contract audit, fits all available windows,
and writes the result under `results/<case>/`.

## Inputs

### 20250527

```text
inputs/20250527/frozen/r5_sensor_conditioned_sidecar.mat
inputs/20250527/frozen/Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S136_20250527.mat
inputs/20250527/foundation/Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat
inputs/20250527/calibration/Step05_Response_Surface_20250527.mat
inputs/20250527/calibration/Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136.mat
```

### 20251222

```text
inputs/20251222/frozen/R5_SensorConditionedSidecar_B1_S123.mat
inputs/20251222/frozen/Template_AdaptiveSG.mat
inputs/20251222/foundation/Result_Step05_NoEtaVPTopKDirectTemplate_B1_S123_20251222.mat
inputs/20251222/calibration/Step05_Response_Surface_20251222.mat
inputs/20251222/calibration/Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat
```

The input check fails explicitly if any required file is absent. No missing
Foundation, template, or gap-library path is guessed.

## Outputs

Each case writes:

```text
results/<case>/r5_sensor_conditioned_dynamic*.mat
```

The result contains the schema, source and staged paths, static audit, method
contract snapshot, window rows, support diagnostics, EO audit, dynamic
parameters, and Foundation replay fields.

Important status meanings:

```text
support_fail          declared support is insufficient
replay_unavailable    Foundation baseline diagnostic unavailable
replay_fail           Foundation baseline does not reproduce the observation
replay_pass           Foundation baseline reproduces the observation
 diagnostic_only      R5 fit retained for diagnosis, not formal acceptance
usable                R5 own support, frequency, EO, and fit checks pass

Foundation replay is diagnostic only. It is not an R5 veto condition and does
not determine `usable` versus `rejected`. The method comparison is the separate
`R5_Compare_Foundation_vs_R5` output.
```

A short Foundation preview cannot substitute for a complete point-aligned
waveform. In particular, a 40-point preview cannot be compared directly with a
full 13,000-point observation bundle.

## Verification

Run the package checks from a clean MATLAB path:

```matlab
restoredefaultpath;
addpath('E:/.../R5_SensorConditioned_Decoupling_20260916');
verify_R5_Installation('case','all');
Test_R5_PackageContracts;
```

The contract tests verify the fixed frequency band, frozen absolute static
state, forbidden dynamic tilt, and the existence of all package-local input
files. A successful MATLAB process exit alone is not sufficient: inspect the
saved MAT result and its `out.schema`, `out.rows`, and replay fields.

## Current Limitations

- The strain-to-displacement transfer ratio is provisional. R5 amplitude `A`
  is a model fit and is not declared correct or incorrect from the current
  strain amplitude conversion.
- 20250527 must remain `replay_unavailable` unless a complete point-aligned
  Foundation replay vector is provided. A preview vector is not enough.
- 20251222 independent strain/BTT posterior input is not fabricated when it is
  absent.
- The package contains frozen prepared inputs for reproducible dynamic reruns;
  raw-data preparation is intentionally not silently triggered by a dynamic
  entry point.

## Legacy Programs

Historical `Main*`, old GapAware, old explicit-profile, and exploratory
programs are not part of this package. They remain archived in their original
case directories under `legacy_programs/` and must not be added to the formal
R5 MATLAB path or mixed with the packaged sidecars/results.
