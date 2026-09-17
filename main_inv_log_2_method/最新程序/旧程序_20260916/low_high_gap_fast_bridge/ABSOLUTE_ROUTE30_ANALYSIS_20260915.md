# Absolute Route-30 A/B conclusion

## Controlled comparison

The two smoke tests use the same three-sensor, 20-revolution low-speed and
8-revolution high-speed contract, with `g_low=0.50 mm`, `g_high=0.70 mm`,
`EO=11`, `A=0.20 mm`, and `phi=0.4 rad`. The only changed item is the
high-speed forward model and its inverse evaluator.

| route | gap estimate (mm) | EO | amplitude (mm) | phase (rad) | voltage RMSE (V) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `low_increment` | 0.70175 | 11 | 0.15681 | not retained by the legacy smoke output | 1.2012e-3 |
| absolute Route-30 | 0.70001 | 11 | 0.19981 | 0.40093 | 5.7922e-5 |

The absolute result was obtained through the public
`run_inv_log_2_low_high_main` entry, with `cfg.route30ForwardModel="absolute"`
and `cfg.route30FrequencyStructureMode="general"`; it is therefore not a
private diagnostic shortcut.

## Interpretation

The forward surface and the sampled data are valid. In the low-increment
case, replaying the true parameters gives `5.62e-15 V`, and the fixed-EO,
fixed-gap amplitude profile has its minimum at the true amplitude. The error
appears only when the final optimizer jointly releases gap, displacement
offset, amplitude, phase, and frequency. The VP seed is biased low and the
joint local basin compensates that bias with gap/offset changes.

The original effective method was the absolute response surface

```text
V(t) = R(g, x(t) - dx - u(t))
u(t) = a sin(2*pi*f*t) + b cos(2*pi*f*t)
```

It scans `(g,f)` in the voltage domain, obtains `(a,b)` by linear projection,
then refines `[g,dx,a,b,f]` against the complete waveform. Amplitude and phase
are reported only after refinement:

```text
A = hypot(a,b)
phi = atan2(b,a)
```

## Recommended production setting

For the current single/unknown-frequency identification, use the absolute
Route-30 path and keep the general frequency router:

```matlab
cfg.route30ForwardModel = "absolute";
cfg.route30FrequencyStructureMode = "general";
```

Do not pass a low-speed calibrated path model to the high-speed inverse in
this mode, and do not use `simulate_highspeed_from_low_increment` for the
corresponding validation. Existing low-increment outputs are retained as an
ablation/diagnostic control only.

The reproducible control script is
`Run_AbsoluteSingleBlindSmoke.m`; its result is written to
`output/absolute_single_blind_smoke/absolute_single_blind_smoke.csv`.
