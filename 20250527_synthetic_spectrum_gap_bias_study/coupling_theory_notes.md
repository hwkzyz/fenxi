# Gap/tilt coupling theory note

## 1. Model hierarchy

Use the following naming to avoid ambiguity:

```text
fixed:
  no high-speed dg, no high-speed residual dmu

gap_fixed_tilt / current gap_only:
  estimate dg_s
  lock mu_s to the low-speed calibrated corr.muGapPerXMm

gap_tilt_residual / current gap_tilt:
  estimate dg_s and a residual high-speed dmu_s
  use only as diagnostic unless dmu_s is protected by constraints
```

If CH6 has only a static equivalent offset `dg6`, the expected relationship is:

```text
gap_tilt_residual with dmu6 = 0 == gap_fixed_tilt / gap_only
```

It is not expected to equal `fixed` unless the `fixed` template already
contains the true CH6 static offset.

## 2. First-order coupling

Near a current solution, the voltage residual can be written as

```text
r ~= J_A dA + J_phi dphi + J_dx ddx + J_g dg + J_mu dmu
```

where the residual tilt column is not independent:

```text
J_mu ~= F_g (x_eval - tau)
x_eval = X - dx - A sin(EO theta + phi) - eta
```

Because `x_eval` contains the vibration term, `J_mu` contains a component that
is synchronous with the vibration. Therefore an unconstrained optimizer can
reduce residual by trading:

```text
A, phi, dx, dg_s <-> dmu_s
```

This trade is especially dangerous when only one sensor, such as CH6, carries
the mismatch. The multi-sensor bundle may improve residual while shifting the
reported global vibration amplitude.

## 3. Why pure dg can still fail numerically

For a clean pure-dg case, the physical optimum exists:

```text
A = A_true
dg6 = dg6_true
dmu6 = 0
```

The Step04 sweep confirms this when the solver is started at the truth or at
the zero-dmu nested point. The independent `gap_tilt` optimization can still
fall into a local basin where `dg6` is too small and `dmu6` is nonzero. That is
an optimizer/identifiability symptom, not evidence that pure `dg6` requires
residual tilt.

## 4. Coupling control

The safest main result is:

```text
estimate q = [A, phi, dx/d0, dg_s]
fix dmu_s = 0
```

If residual dmu is needed for diagnostics, protect the main solution:

```text
J_mu_perp = (I - P_Jq) J_mu
```

Then fit residual dmu using `J_mu_perp`, not raw `J_mu`. In a nonlinear
implementation, recompute the projection at the current protected solution and
apply a trust region and physical prior:

```text
min ||W^(1/2) (r(q) - J_mu_perp dmu)||^2 + lambda ||dmu / sigma_mu||^2
```

The reported vibration amplitude remains the protected `q` solution. The
residual `dmu_s` is reported as a diagnostic residual shape parameter, not as
an independent physical tilt unless it passes a noise/significance test.
