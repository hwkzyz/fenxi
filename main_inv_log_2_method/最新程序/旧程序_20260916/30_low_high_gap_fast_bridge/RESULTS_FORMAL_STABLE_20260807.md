# Formal shared-gap stabilization results (2026-08-07)

## Failure chain

The unified entry originally applied the query-safe mask to the absolute
template. This changed the VP screening objective and reduced the dual truth
recall. Restricting that mask to the finite-support low-increment model
restored two of the three original dual cases.

The remaining failures were candidate-coverage and candidate-transition
problems, not gap-estimation or forward-model failures:

- `900+1400 Hz`: the true full-voltage parameters had lower RMSE than the
  selected wrong valley, but the retained pool was initially 50 Hz away.
- `700+1200 Hz` with `g_low=0.8 mm`, `g_high=0.2 mm`: conditional completion
  inserted the true pair, but its complete-waveform linear replay rank was 4.

## Implemented bridge

```text
single VP + baseline dual VP
-> provisional BIC order
-> dual records only: one-dimensional conditional completion
-> complete-waveform linear reprojection
-> Top-10 short nonlinear voltage refinement
-> one full nonlinear polish
```

The conditional anchor grid is `300:100:1500 Hz`, while the scanned missing
frequency retains the production 5 Hz grid. Therefore the additional work is
linear in the frequency-grid size and does not repeat the full pair grid.

## Current evidence

`Run_FormalMainRegression`:

- single: 3/3 successes (`531`, `733`, `917 Hz`);
- dual: 3/3 successes (`500+1300`, `700+1200`, `900+1400 Hz`);
- maximum dual frequency error: 0.092 Hz;
- maximum single frequency error: 0.065 Hz.

`Run_FormalDualRepeatRegression(1:3)`:

- 9/9 dual successes over three independent high-speed noise seeds;
- maximum frequency error: 0.201 Hz;
- gap error remained within 0.001 mm in the tested cases;
- elapsed high-speed identification time: 11.4--21.4 s per case.

These results establish a stable matched-template, 20 dB shared-gap simulation
route for the tested cases. They do not yet prove stability under model
mismatch, the low-speed increment template, arbitrary SNR, or independent
per-sensor gaps.
