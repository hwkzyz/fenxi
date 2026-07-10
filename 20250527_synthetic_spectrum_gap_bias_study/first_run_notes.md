# First Run Notes

Generated on 2026-07-09 from the current S136 Step5 bundle:

```text
Result_single_sync_20250527_B1_S136_20250526_2500_3500_t400_Start1p5s.mat
```

## Step01 Strain Prior

The Step6 strain-reference CSV gives:

```text
median strain frequency = 580.1171 Hz
target EO              = 14
nominal rotor freq     = 41.4369 Hz
median strain amplitude= 371.0926 microstrain
current Step6 scale    = 0.3711 mm
synthetic spectrum peak= 579.8340 Hz, 365.3118 microstrain
```

This preserves the experimental conclusion that a plausible EO14 tip
amplitude is about `0.35-0.39 mm`, depending on the strain-to-tip scale.

## Step02 Voltage Bias Study

In clean/noise-only cases, S13 and S136 recover the injected amplitude
`A_true = 0.3711 mm` with errors below about `0.1%`.

When CH6 has a static equivalent shift:

```text
case: ch6_static_dg, S136
fixed    A = 0.3596 mm, error = -3.10%
gap_only A = 0.3699 mm, error = -0.32%
gap_tilt A = 0.3549 mm, error = -4.36%
```

When CH6 has a tilt-like equivalent shift:

```text
case: ch6_tilt_dmu, S136
fixed    A = 0.4422 mm, error = +19.15%
gap_only A = 0.4165 mm, error = +12.24%
gap_tilt A = 0.3824 mm, error = +3.05%
```

When CH6 has both static shift and tilt-like shift:

```text
case: ch6_dg_dmu, S136
fixed    A = 0.4260 mm, error = +14.81%
gap_only A = 0.3980 mm, error = +7.26%
gap_tilt A = 0.3420 mm, error = -7.84%
```

This last case is closest to the experimental concern: adding a flexible
tilt-like correction can flip the amplitude bias downward even when fixed and
gap_only are biased upward.

## Step03 Jacobian Coupling

For S136, CH6 has:

```text
rho_A_dg        = 0.7327
proj_dg_on_vib  = 0.7699
rho_A_dmu       = 0.0523
proj_dmu_on_vib = 0.6600
```

The direct pairwise `A`-`dmu` correlation is not large, but about `66%` of the
CH6 tilt column lies in the vibration subspace `[A, phi, d0]`. This supports
the mechanism that dmu can trade against the vibration fit through a
multi-parameter path, not only through direct A-column correlation.

## Current Interpretation

The synthetic evidence supports this working hypothesis:

```text
CH6 static/tilt mismatch can push fixed or super-Gaussian-style amplitude high.
Allowing tilt correction can reduce the voltage residual but may also pull
the recovered vibration amplitude downward when shift and tilt coexist.
```

The next useful extension is a parameter sweep over `dg6`, `dmu6`, and the
allowed/regularized `dmu` range, rather than relying on one hand-picked CH6
scenario.

## Step04 Static-dg Identifiability Sweep

The pure static CH6 shift case clarifies an important point. With clean data,
the correct `gap_tilt` solution exists:

```text
case: dg_clean, truth_start
gap_tilt A = 0.37109 mm, dg6 = 0.05000 mm, dmu6 = 0
```

Starting `gap_tilt` from the `gap_only` zero-dmu point also returns the same
amplitude to numerical accuracy:

```text
case: dg_clean, gap_only_zero_dmu_start
gap_tilt A = 0.37116-0.37127 mm, dg6 ~= 0.04894-0.04995 mm, dmu6 ~= 0
```

However, the independent nonlinear `gap_tilt` starts often converge to a poor
local basin even in the clean pure-dg case:

```text
case: dg_clean, independent gap_tilt
A ~= 0.3429-0.3568 mm, error ~= -7.6% to -3.9%
dg6 ~= 0.014-0.026 mm, dmu6 may move away from zero
```

So the pure-dg failure is not a physical necessity of the model. It is a
combination of over-parameterized static freedom and a difficult nonlinear
optimization landscape. The correct theoretical expectation is:

```text
pure dg6 and exact model:
gap_tilt with dmu6 = 0 should reduce to gap_only / gap_fixed_tilt.
```

It should equal `fixed` only if `fixed` already contains the true CH6 static
offset. In the current synthetic and Step07J naming, `fixed` has no high-speed
`dg_s`, so pure CH6 dg can bias `fixed`.

## Step05 Orthogonalized Residual Tilt

The linearized diagnostic protects the subspace
`[A, phi, d0, dg_s]` and projects residual `dmu_s` onto its orthogonal
complement. For CH6:

```text
proj_dmu6_on_protected = 0.7955
dmu6_perp_norm_fraction = 0.6060
```

Thus about 80% of the CH6 residual-tilt column lies in the vibration/static-gap
subspace. In the noisy pure-dg case:

```text
linear_gap_only       A = 0.37150 mm, error = +0.11%
linear_gap_tilt_raw   A = 0.36915 mm, error = -0.52%, dmu6 = 0.00364
linear_gap_tilt_orth  A = 0.37150 mm, error = +0.11%, dmu6 = 0.00364
```

The orthogonalized residual dmu can still reduce residual waveform error, but
it cannot change the first-order amplitude/static-gap solution. This is the
right diagnostic behavior if residual dmu is retained only as a post-fit
explanation of remaining waveform structure.

## Updated Working Conclusion

The main paper-facing amplitude should come from the model that estimates
`dg_s` while locking the tilt slope to the low-speed calibrated
`corr.muGapPerXMm`:

```text
main:       gap_fixed_tilt / current gap_only
diagnostic: gap_tilt_residual / current gap_tilt
```

If a residual high-speed `dmu_s` is kept, it should be protected from stealing
vibration by either:

1. orthogonalizing the residual dmu column against `[A, phi, dx/d0, dg_s]`;
2. fitting residual dmu only after the main protected fit;
3. adding a physical prior/significance test based on low-speed uncertainty
   and accepting dmu only when the residual reduction exceeds the noise level.

## Step06 Gap/Tilt Direction Decoupling

For static gap/tilt only, CH6 in S136 gives:

```text
rho_gap_tilt_raw                       = 0.0639
gap_equivalent_shift_per_unit_dmu      = -0.1213 mm
tilt_norm_fraction_after_gap_orth      = 0.9980
proj_tilt_on_vib_raw                  = 0.6600
proj_tilt_on_protected_vib_gap         = 0.7180
```

Thus the raw CH6 tilt column contains only a small direct `dg` component in
the current weighted observation metric. Moving the tilt pivot to the
weighted center can cleanly remove that small static `dg` component, but the
dominant risk remains the projection of the tilt column onto the vibration
subspace.

## Step07 Protected Two-Stage Residual Tilt

A two-stage diagnostic was tested:

```text
Stage 1: nonlinear gap_fixed_tilt / current gap_only, report A
Stage 2: fit residual dmu on the Stage-1 residual, with A locked
```

For pure CH6 static offset:

```text
stage1_gap_fixed_tilt        A = 0.37293 mm, error = +0.50%
current_gap_tilt_independent A = 0.35634 mm, error = -3.97%
stage2_safe_dmu_A_locked     A = 0.37293 mm, dmu6 = 0.00030
```

This is the desired behavior: residual dmu does not pull the reported
amplitude down when the true mismatch is only `dg6`.

For true CH6 tilt mismatch:

```text
ch6_tilt_dmu:
stage1_gap_fixed_tilt        A = 0.41917 mm, error = +12.96%
stage2_safe_dmu_A_locked     A = 0.41917 mm, dmu6 = 0.03573
```

The safe residual dmu finds the injected tilt but cannot repair the part of
the tilt effect that is indistinguishable from `[A, phi, d0, dg]`. This is an
identifiability limit, not just an optimization issue. Recovering the buried
component requires an external prior or extra experiment design, such as
low-speed tilt calibration, multi-window shared static dmu, excluding suspect
sensors, or an independent strain/FE amplitude reference.
