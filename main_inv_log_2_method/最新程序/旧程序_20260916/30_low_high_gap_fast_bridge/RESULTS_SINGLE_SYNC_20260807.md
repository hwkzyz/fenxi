# Structured single-synchronous audit (2026-08-07)

The prototype uses the measured revolution timing to construct shaft angle and
fits integer engine-order bases `sin(k*theta)` and `cos(k*theta)`. All 8912
trusted samples are retained. Large-amplitude protection evaluates every
integer order and gap candidate with exact-forward amplitude/phase replay before
the final fixed-order nonlinear optimization.

## Results

| SNR | EO 10 / 500 Hz / 0.10 mm | EO 14 / 700 Hz / 0.15 mm | EO 18 / 900 Hz / 0.25 mm |
|---:|---|---|---|
| 20 dB | correct | correct | correct |
| 10 dB | correct | correct | correct |
| 5 dB | correct | correct | selected EO 12 / 600 Hz |

At 5 dB in the EO-18 case, the selected EO-12 solution has voltage RMSE
`0.036005 V`, while the generating EO-18 model has `0.036008 V`. The shaft-angle
representation error is below `2.3e-14 mm`, so the failure is not caused by
angle conversion or the integer-order basis. The wrong order is slightly better
under this noise realization; this is an objective-resolution boundary.

At 10 and 20 dB, EO 18 is recovered exactly, with estimated amplitude `0.25039`
and `0.24996 mm`, respectively. Thus a 0.25 mm amplitude does not make VP or the
synchronous model fail at every SNR. The 5 dB failure results from the combination
of noise and a nearly tied competing order for the tested gap/window.

## Implication

Integer-order enumeration is useful because it is finite, but it cannot force
the generating order to have the minimum noisy voltage objective. A practical
synchronous solver should use a cheap VP pass for normal cases and trigger exact
full-order replay only when the nonlinear replay margin is small. At a near tie,
the result should include competing orders and a confidence flag rather than
claim a unique order.
